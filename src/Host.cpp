#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;
using namespace zppsim;

Host::Host(Population * popPtr, size_t id, double deathTime) :
	popPtr(popPtr), id(id), deathTime(deathTime), nextInfectionId(0),
	deathEvent(new DeathEvent(this))
{
	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	popPtr->addEvent(deathEvent.get());
}

void Host::die()
{
	cerr << popPtr->getTime() << ", death event: " << popPtr->id << ", " << id << '\n';
	
	id = popPtr->nextHostId++;
	deathTime = popPtr->simPtr->getTime() + popPtr->simPtr->drawHostLifetime();
	
	cerr << "new id: " << id << '\n';
	cerr << "new death time: " << deathTime << '\n';
	
	// Remove infections and immune history
	for(auto & infection : infections) {
		if(infection.transitionEvent != nullptr) {
			popPtr->removeEvent(infection.transitionEvent.get());
		}
		if(infection.clearanceEvent != nullptr) {
			popPtr->removeEvent(infection.clearanceEvent.get());
		}
	}
	infections.clear();
	immunity.clear();
	
	popPtr->setEventTime(deathEvent.get(), deathTime);
}

void Host::transmitTo(Host & dstHost, rng_t & rng, double pRecombination)
{
	if(infections.size() == 0) {
		cerr << "No infections to transmit" << endl;
		return;
	}
	
	// Get all current infections
	vector<StrainPtr> strains;
	for(auto infItr = infections.begin(); infItr != infections.end(); infItr++) {
		strains.push_back(infItr->strainPtr);
	}
	
	// Make daughter strains for subset of current infections
	if(strains.size() > 1) {
		cerr << "Multiple infections" << '\n';
		vector<pair<size_t, size_t>> indexPairs = drawMultipleBernoulliIndexPairs(
			rng, strains.size(), pRecombination, false, true
		);
		for(auto & indexPair : indexPairs) {
			cerr << "Drew recombination pair " << indexPair.first << ", " << indexPair.second << '\n';
			
			strains.push_back(popPtr->simPtr->recombineStrains(
				strains[indexPair.first],
				strains[indexPair.second]
			));
		}
	}
	
	// TODO: actually transmit infections
}

void Host::receiveInfection(StrainPtr & strainPtr)
{
	assert(strainPtr->size() > 0);
	
	double time = popPtr->getTime();
	double tLiverStage = popPtr->simPtr->parPtr->tLiverStage;
	size_t initialGeneIndex = tLiverStage == 0 ? 0 : LIVER_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, initialGeneIndex);
	
	// Either start in the liver stage and create a fixed-time transition event,
	// or start in the preactivation stage for the first var gene with a rate-based
	// transition event and a clearance event
	list<Infection>::reverse_iterator reverseItr = infections.rbegin();
	reverseItr++;
	list<Infection>::iterator infectionItr = reverseItr.base();
	if(initialGeneIndex == LIVER_STAGE) {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(infectionItr, time + tLiverStage)
		);
		popPtr->addEvent(infectionItr->transitionEvent.get());
		// Zero-rate clearance event (placeholder)
		infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
			new ClearanceEvent(
				infectionItr, 0.0, time, *(popPtr->rngPtr)
			)
		);
		popPtr->addEvent(infectionItr->clearanceEvent.get());
	}
	else {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(
				infectionItr, infectionItr->activationRate(0), time, *(popPtr->rngPtr)
			)
		);
		popPtr->addEvent(infectionItr->transitionEvent.get());
		infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
			new ClearanceEvent(
				infectionItr, infectionItr->clearanceRate(), time, *(popPtr->rngPtr)
			)
		);
		popPtr->addEvent(infectionItr->clearanceEvent.get());
	}
	
	cerr << time << ": " << popPtr->id << ", " << id << ", " << infectionItr->id << " begun" << '\n';
}

void Host::clearInfection(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	double time = popPtr->getTime();
	cerr << time << ": " << popPtr->id << ", " << id << ", " << infectionItr->id << " clearing" << '\n';
	
	// Remove any events
	if(infectionItr->transitionEvent != nullptr) {
		popPtr->removeEvent(infectionItr->transitionEvent.get());
	}
	if(infectionItr->clearanceEvent != nullptr) {
		popPtr->removeEvent(infectionItr->clearanceEvent.get());
	}
	
	// Remove infection
	infections.erase(infectionItr);
}

void Host::performTransition(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	double time = popPtr->getTime();
	
	size_t geneIndex = infectionItr->currentGeneIndex;
	bool active = infectionItr->active;
	if(geneIndex == infectionItr->strainPtr->size() - 1) {
		clearInfection(infectionItr);
		return;
	}
	else if(geneIndex == LIVER_STAGE) {
		infectionItr->currentGeneIndex = 0;
		assert(!active);
	}
	else {
		if(active) {
			infectionItr->currentGeneIndex++;
			infectionItr->active = false;
		}
		else {
			infectionItr->active = true;
		}
	}
	
	// Update transition rate
	if(infectionItr->active) {
		popPtr->setEventRate(
			infectionItr->transitionEvent.get(),
			infectionItr->deactivationRate(infectionItr->currentGeneIndex)
		);
	}
	else {
		popPtr->setEventRate(
			infectionItr->transitionEvent.get(),
			infectionItr->activationRate(infectionItr->currentGeneIndex)
		);
	}
	
	// Update clearance rate
	popPtr->setEventRate(
		infectionItr->clearanceEvent.get(),
		infectionItr->clearanceRate()
	);
	
	cerr << time << ": " << popPtr->id << ", " << id << ", " << infectionItr->id << " transitioned to " <<
		infectionItr->currentGeneIndex << ", "
		<< (infectionItr->active ? "active" : "not yet active") << '\n';
}

Infection::Infection(Host * hostPtr, size_t id, StrainPtr & strainPtr, size_t initialGeneIndex) :
	hostPtr(hostPtr), id(id), strainPtr(strainPtr),
	currentGeneIndex(initialGeneIndex), active(false)
{
}
	
double Infection::activationRate(size_t geneIndex)
{
	return 1.0;
}

double Infection::deactivationRate(size_t geneIndex)
{
	return 1.0;
}

double Infection::clearanceRate()
{
	return 0.1;
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->die();
}

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
	infectionItr->hostPtr->performTransition(infectionItr);
}

void ClearanceEvent::performEvent(zppsim::EventQueue &queue)
{
	infectionItr->hostPtr->clearInfection(infectionItr);
}
