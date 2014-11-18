#include "Host.h"
#include "Simulation.h"

#include <iostream>
#include <sstream>

using namespace std;
using namespace zppsim;

Host::Host(
	Population * popPtr, int64_t id, double birthTime, double deathTime,
	bool writeToDatabase,
	Database & db,
	zppdb::Table<HostRow> & table
):
	id(id), popPtr(popPtr),
	birthTime(birthTime), deathTime(deathTime), nextInfectionId(0),
	deathEvent(new DeathEvent(this)),
	immunity(this, false),
	clinicalImmunity(this, true)
{
//	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	if(writeToDatabase) {
		HostRow row;
		row.hostId = id;
		row.birthTime = birthTime;
		row.deathTime = deathTime;
		db.insert(table, row);
	}
	addEvent(deathEvent.get());
}

double Host::getAge()
{
	return getTime() - birthTime;
}

int64_t Host::getActiveInfectionCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			count++;
		}
	}
	return count;
}

std::vector<GenePtr> Host::getActiveInfectionGenes()
{
	vector<GenePtr> genes;
	genes.reserve(infections.size());
	for(auto & infection : infections) {
		if(infection.isActive()) {
			genes.push_back(infection.getCurrentGene());
		}
	}
	return genes;
}

std::vector<int64_t> Host::getActiveInfectionGeneIds()
{
	vector<int64_t> ids;
	ids.reserve(infections.size());
	for(auto & infection : infections) {
		if(infection.isActive()) {
			ids.push_back(infection.getCurrentGeneId());
		}
	}
	return ids;
}

int64_t Host::getActiveInfectionImmunityCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			if(immunity.isImmune(infection.getCurrentGene())) {
				count++;
			}
		}
	}
	return count;
}

int64_t Host::getActiveInfectionClinicalImmunityCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			if(clinicalImmunity.isImmune(infection.getCurrentGene())) {
				count++;
			}
		}
	}
	return count;
}

void Host::prepareToDie()
{
//	cerr << popPtr->getTime() << ", host going to die: " << toString() << '\n';
	
	// Remove all events
	for(auto & infection : infections) {
		infection.prepareToEnd();
	}
	immunity.prepareToDie();
	clinicalImmunity.prepareToDie();
	
	removeEvent(deathEvent.get());
}

void Host::transmitTo(Host & dstHost)
{
	rng_t * rngPtr = getRngPtr();
	
	if(infections.size() == 0) {
	//	cerr << "No infections to transmit" << endl;
		return;
	}
//	cerr << "Transmitting to " << dstHost.popPtr->id << ", " << dstHost.id << '\n';
	
	// Get some current infections according to transmission probability
	vector<StrainPtr> originalStrains;
	for(auto infItr = infections.begin(); infItr != infections.end(); infItr++) {
		if(infItr->isActive()) {
			bernoulli_distribution flipCoin(infItr->transmissionProbability());
			if(flipCoin(*rngPtr)) {
				originalStrains.push_back(
					popPtr->simPtr->mutateStrain(infItr->strainPtr)
				);
			}
		}
	}
	
	vector<StrainPtr> strainsToTransmit;
	if(originalStrains.size() > 1) {
		// Take each original strain with probability 1.0 - pRecombinant
		double pRecombinant = popPtr->simPtr->parPtr->pRecombinant;
		assert(pRecombinant >= 0.0 && pRecombinant <= 1.0);
		strainsToTransmit = drawMultipleBernoulliRandomSubset(
			*rngPtr, originalStrains, 1.0 - pRecombinant, false
		);
		
		// Complete a set of size originalStrains.size() using recombinants
		int64_t nRecombinants = originalStrains.size() - strainsToTransmit.size();
		for(int64_t i = 0; i < nRecombinants; i++) {
			uniform_int_distribution<int64_t> indDist(0, originalStrains.size() - 1);
			int64_t ind1 = indDist(*rngPtr);
			int64_t ind2 = indDist(*rngPtr);
			strainsToTransmit.push_back(
				popPtr->simPtr->recombineStrains(
					originalStrains[ind1],
					originalStrains[ind2]
				)
			);
		}
	}
	else {
		strainsToTransmit = originalStrains;
	}
	
	// Transmit all strains
	for(auto & strainPtr : strainsToTransmit) {
		dstHost.receiveInfection(strainPtr);
	}
}

void Host::receiveInfection(StrainPtr & strainPtr)
{
	assert(strainPtr->size() > 0);
	
	popPtr->countTransmission();
	
	rng_t * rngPtr = getRngPtr();
	double time = popPtr->getTime();
	
	double tLiverStage = popPtr->simPtr->parPtr->tLiverStage;
	int64_t initialGeneIndex = tLiverStage == 0 ? 0 : WAITING_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, initialGeneIndex, time);
	
	
	// Get an iterator to the infection so that events can
	// point back to their infection
	list<Infection>::reverse_iterator reverseItr = infections.rbegin();
	reverseItr++;
	list<Infection>::iterator infectionItr = reverseItr.base();
	
	// If starting in liver stage, create a fixed-time transition event
	// (liver stage -> first gene not yet active)
	if(initialGeneIndex == WAITING_STAGE) {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(infectionItr, time + tLiverStage)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	// Otherwise create a rate-based transition event
	// (first gene not yet active -> first gene active)
	else {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(
				infectionItr,
				infectionItr->transitionRate(),
				time,
				*rngPtr
			)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	
	// Create a clearance event
	// (rate will depend on state as determined in clearanceRate()
	// and may be zero)
	infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
		new ClearanceEvent(
			infectionItr, infectionItr->clearanceRate(), time, *rngPtr
		)
	);
	addEvent(infectionItr->clearanceEvent.get());
	
//	cerr << time << ": " << infectionItr->toString() << " begun" << '\n';
}

void Host::updateInfectionRates()
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		if(itr->geneIndex != WAITING_STAGE) {
			itr->updateClearanceRate();
			itr->updateTransitionRate();
		}
	}
}

void Host::clearInfection(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	bool shouldUpdateAllRates = infectionItr->transitionAffectsAllInfections();
	
//	double time = popPtr->getTime();
//	cerr << time << ": " << infectionItr->toString() << " clearing" << '\n';
	
	// Gain immunity to active gene
	if(infectionItr->active) {
		GenePtr genePtr = infectionItr->getCurrentGene();
		immunity.gainImmunity(genePtr);
		if(getSimulationParametersPtr()->trackClinicalImmunity) {
			clinicalImmunity.gainImmunity(genePtr);
		}
	}
	
	// Remove infection
	infectionItr->prepareToEnd();
	infections.erase(infectionItr);
	
	if(shouldUpdateAllRates) {
		updateInfectionRates();
	}
}

double Host::getTime()
{
	return popPtr->getTime();
}

rng_t * Host::getRngPtr()
{
	return popPtr->rngPtr;
}

SimParameters * Host::getSimulationParametersPtr()
{
	return popPtr->simPtr->parPtr;
}

PopulationParameters * Host::getPopulationParametersPtr()
{
	return popPtr->parPtr;
}

void Host::addEvent(Event * event)
{
	popPtr->addEvent(event);
}

void Host::removeEvent(Event * event)
{
	popPtr->removeEvent(event);
}

void Host::setEventRate(RateEvent * event, double rate)
{
	popPtr->setEventRate(event, rate);
}

void Host::writeInfections(Database & db, Table<InfectionRow> & table)
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		itr->write(db, table);
	}
}

std::string Host::toString()
{
	stringstream ss;
	ss << popPtr->toString() << ".h" << id;
	return ss.str();
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->prepareToDie();
	hostPtr->popPtr->removeHost(hostPtr);
}
