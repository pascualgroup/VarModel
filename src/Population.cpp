//
//  Population.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 4/21/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "Population.h"
#include "Simulation.h"
#include "SimParameters.h"
#include <iostream>

using namespace std;
using namespace zppsim;

Population::Population(Simulation * simPtr, size_t id) :
	id(id), simPtr(simPtr), parPtr(&(simPtr->parPtr->populations[id])),
	nextHostId(0)
{
	// Create hosts
	hosts.reserve(parPtr->size);
	for(size_t i = 0; i < parPtr->size; i++) {
		hosts.emplace_back(new Host(this, nextHostId++, simPtr->drawHostLifetime()));
	}
	
	// Create biting event
	bitingEvent = unique_ptr<BitingEvent>(new BitingEvent(this, simPtr->rng));
	simPtr->addEvent(bitingEvent.get());
	
	// Create initial infections
	for(size_t i = 0; i < parPtr->nInitialInfections; i++) {
		size_t hostId = drawUniformIndex(simPtr->rng, hosts.size());
		StrainPtr strainPtr = simPtr->generateRandomStrain();
		
		hosts[hostId]->receiveInfection(strainPtr);
	}
}

/*void Population::pushBackEvents(std::vector<Event *> & eventVec)
{
	eventVec.push_back(bitingEvent.get());
	for(auto & hostPtr : hosts) {
		hostPtr->pushBackEvents(eventVec);
	}
}*/

double Population::bitingRate()
{
	return parPtr->bitingRate * hosts.size();
}

size_t Population::size()
{
	return hosts.size();
}

Host * Population::getHost(size_t hostId)
{
	return hosts[hostId].get();
}

double Population::getTime()
{
	return simPtr->getTime();
}

void Population::addEvent(zppsim::Event * event)
{
	simPtr->addEvent(event);
}
	
void Population::removeEvent(zppsim::Event * event)
{
	simPtr->removeEvent(event);
}

void Population::setEventTime(zppsim::Event * event, double time)
{
	simPtr->setEventTime(event, time);
}

void Population::performBitingEvent()
{
	cerr << simPtr->getTime() << ": biting event pop " << id << '\n';
	
	// Choose random host to get bitten
	size_t dstHostId = drawUniformIndex(simPtr->rng, hosts.size());
	Host * dstHostPtr = getHost(dstHostId);
	
	// Choose source host
	Host * srcHostPtr = simPtr->drawSourceHost(id, dstHostId);
	assert(srcHostPtr != dstHostPtr);
	
	cerr << "dst host: " << dstHostId << '\n';
	cerr << "src pop, host: " << srcHostPtr->popPtr->id << ", " << srcHostPtr->id << '\n';
}

/*** BITING EVENT ***/

BitingEvent::BitingEvent(Population * popPtr, zppsim::rng_t & rng) :
	RateEvent(popPtr->bitingRate(), 0.0, rng),
	popPtr(popPtr)
{
}

void BitingEvent::performEvent(zppsim::EventQueue & queue)
{
	popPtr->performBitingEvent();
}
