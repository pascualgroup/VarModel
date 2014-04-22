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
}

void Population::pushBackEvents(std::vector<Event *> & eventVec)
{
	eventVec.push_back(bitingEvent.get());
	for(auto & hostPtr : hosts) {
		hostPtr->pushBackEvents(eventVec);
	}
}

double Population::bitingRate()
{
	return parPtr->bitingRate * hosts.size();
}

void Population::performBitingEvent()
{
	cerr << simPtr->getTime() << ": biting event pop " << id << '\n';
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
