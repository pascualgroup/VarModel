#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;
using namespace zppsim;

Host::Host(Population * popPtr, size_t id, double deathTime)
	: popPtr(popPtr), id(id), deathTime(deathTime), deathEvent(new DeathEvent(this))
{
	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	popPtr->addEvent(deathEvent.get());
}

void Host::die(zppsim::EventQueue & queue)
{
	cerr << queue.getTime() << ", death event: " << popPtr->id << ", " << id << '\n';
	
	id = popPtr->nextHostId++;
	deathTime = popPtr->simPtr->getTime() + popPtr->simPtr->drawHostLifetime();
	
	cerr << "new id: " << id << '\n';
	cerr << "new death time: " << deathTime << '\n';
	
	deathEvent->setTime(queue, deathTime);
}

void Host::transmitTo(Host & dstHost)
{
}

void Host::receiveInfection(StrainPtr & strainPtr)
{
	assert(strainPtr->size() > 0);
	
	// If there's no infection delay, set initial stage to 0;
	// otherwise set initial stage to INFECTION_STAGE_NULL
	// (infection not yet started)
	double initialDelayTime = calculateTransitionDelay(INFECTION_STAGE_NULL);
	size_t initialStage;
	if(initialDelayTime == 0.0) {
		initialDelayTime = calculateTransitionDelay(0);
		initialStage = 0;
	}
	else {
		initialStage = INFECTION_STAGE_NULL;
	}
	
	infections.emplace_back(
		this,
		strainPtr, initialStage
	);
	
	assert(initialDelayTime > 0.0);
	
	// Create transition event with reference back to list iterator;
	// add event to infection and add event to event queue
	list<Infection>::reverse_iterator rItr = infections.rbegin();
	rItr++;
	list<Infection>::iterator infectionItr = rItr.base();
	assert(infectionItr != infections.end());
	
	infectionItr->nextTransition = unique_ptr<TransitionEvent>(new TransitionEvent(infectionItr, popPtr->getTime() + initialDelayTime));
	assert(infectionItr->nextTransition->infectionItr != infections.end());
	assert(infectionItr->nextTransition->getTime() > 0.0);
	popPtr->addEvent(infectionItr->nextTransition.get());
}

void Host::performTransition(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	StrainPtr strainPtr = infectionItr->strainPtr;
	double time = popPtr->getTime();
	
	// TODO: record immunity due to just-completed stage
	
	cerr << time << ": infection transition pop "
		<< popPtr->id << ", host " << id << ", " << infectionItr->stage
		<< '\n';
	
	// If infection is in the last stage, remove the transition event
	// and remove the infection
	if(infectionItr->stage == strainPtr->size() - 1) {
		popPtr->removeEvent(infectionItr->nextTransition.get());
		infections.erase(infectionItr);
	}
	// Otherwise transition to next stage
	else {
		if(infectionItr->stage == INFECTION_STAGE_NULL) {
			infectionItr->stage = 0;
		}
		else {
			infectionItr->stage++;
		}
		double delay = calculateTransitionDelay(infectionItr->stage);
		popPtr->setEventTime(
			infectionItr->nextTransition.get(),
			time + delay
		);
	}
}

//void Host::pushBackEvents(std::vector<zppsim::Event *> & eventVec)
//{
//	eventVec.push_back(deathEvent.get());
//}

double Host::calculateTransitionDelay(size_t fromStage)
{
	if(fromStage == INFECTION_STAGE_NULL) {
		return 0.1;
	}
	else {
		return 0.5;
	}
}

Infection::Infection(Host * hostPtr, StrainPtr & strainPtr, size_t initialStage) :
	hostPtr(hostPtr), strainPtr(strainPtr), stage(initialStage)
{
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->die(queue);
}

TransitionEvent::TransitionEvent(std::list<Infection>::iterator infectionItr, double time) :
	Event(time), infectionItr(infectionItr)
{
}

void TransitionEvent::performEvent(zppsim::EventQueue &queue)
{
	infectionItr->hostPtr->performTransition(infectionItr);
}
