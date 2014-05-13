#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;
using namespace zppsim;

Host::Host(Population * popPtr, size_t id, double birthTime, double deathTime, DBTable * table) :
	id(id), popPtr(popPtr),
	birthTime(birthTime), deathTime(deathTime), nextInfectionId(0),
	deathEvent(new DeathEvent(this)),
	immunity(this, false),
	clinicalImmunity(this, true)
{
//	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	if(table != nullptr) {
		DBRow row;
		row.set("hostId", int64_t(id));
		row.set("birthTime", birthTime);
		row.set("deathTime", deathTime);
		table->insert(row);
	}
	
	addEvent(deathEvent.get());
}

double Host::getAge()
{
	return getTime() - birthTime;
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
		size_t nRecombinants = originalStrains.size() - strainsToTransmit.size();
		for(size_t i = 0; i < nRecombinants; i++) {
			uniform_int_distribution<size_t> indDist(0, originalStrains.size() - 1);
			size_t ind1 = indDist(*rngPtr);
			size_t ind2 = indDist(*rngPtr);
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
	size_t initialGeneIndex = tLiverStage == 0 ? 0 : WAITING_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, initialGeneIndex);
	
	
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
	
	double time = popPtr->getTime();
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

void Host::writeInfections(DBTable * table)
{
	if(table != nullptr) {
		for(auto itr = infections.begin(); itr != infections.end(); itr++) {
			itr->write(table);
		}
	}
}

std::string Host::toString()
{
	return popPtr->toString() + ".h" + strprintf("%u", id);
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

/*** INFECTION PROCESS EVENTS ***/

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	RateEvent(time), infectionItr(infectionItr)
{
}

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	RateEvent(rate, time, rng), infectionItr(infectionItr)
{
}

TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	InfectionProcessEvent(infectionItr, time)
{
}


TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

ClearanceEvent::ClearanceEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

void TransitionEvent::performEvent(zppsim::EventQueue &queue)
{
	// If this is the final deactivation, then it's equivalent to clearing
	if(infectionItr->active
		&& infectionItr->geneIndex == infectionItr->strainPtr->size() - 1
	) {
		infectionItr->hostPtr->clearInfection(infectionItr);
	}
	// Otherwise, actually perform a transition
	else {
		infectionItr->performTransition();
	}
}

void ClearanceEvent::performEvent(zppsim::EventQueue &queue)
{
	infectionItr->hostPtr->clearInfection(infectionItr);
}
