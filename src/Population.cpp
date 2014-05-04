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
	id(id), simPtr(simPtr), rngPtr(&(simPtr->rng)),
	parPtr(&(simPtr->parPtr->populations[id])),
	nextHostId(0)
{
	// Create hosts
	hosts.reserve(parPtr->size);
	for(size_t i = 0; i < parPtr->size; i++) {
		hosts.emplace_back(new Host(this, nextHostId++, simPtr->drawHostLifetime()));
	}
	
	// Create biting event
	bitingEvent = unique_ptr<BitingEvent>(
		new BitingEvent(this, getBitingRate(), simPtr->rng)
	);
	addEvent(bitingEvent.get());
	
	// Create immigration event
	immigrationEvent = unique_ptr<ImmigrationEvent>(
		new ImmigrationEvent(this, getImmigrationRate(), simPtr->rng)
	);
	addEvent(immigrationEvent.get());
	
	// Create initial infections
	for(size_t i = 0; i < parPtr->nInitialInfections; i++) {
		size_t hostId = drawUniformIndex(simPtr->rng, hosts.size());
		StrainPtr strainPtr = simPtr->generateRandomStrain();
		
		hosts[hostId]->receiveInfection(strainPtr);
	}
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

double Population::getBitingRate()
{
	BitingRate brObj = parPtr->bitingRate;
	return brObj.mean + brObj.amplitude * simPtr->getSeasonality();
}

double Population::getImmigrationRate()
{
	return parPtr->immigrationRate;
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

void Population::setEventRate(zppsim::RateEvent * event, double rate)
{
	simPtr->setEventRate(event, rate);
}

void Population::performBitingEvent()
{
	cerr << simPtr->getTime() << ": biting event, src pop " << id << '\n';
	
	size_t srcHostIndex = drawUniformIndex(*rngPtr, hosts.size());
	Host * srcHostPtr = hosts[srcHostIndex].get();
	cerr << "src host: " << srcHostPtr->id << '\n';
	
	Host * dstHostPtr = simPtr->drawDestinationHost(id);
	cerr << "dst pop, host: " << dstHostPtr->popPtr->id << ", " << dstHostPtr->id << '\n';
	
	srcHostPtr->transmitTo(*dstHostPtr, *rngPtr, simPtr->parPtr->pRecombination);
}

void Population::performImmigrationEvent()
{
	cerr << getTime() << ": immigration event, pop " << id << '\n';
}

double Population::getDistance(Population * popPtr)
{
	if(popPtr == this) {
		return parPtr->selfDistance;
	}
	
	double x1 = parPtr->x;
	double y1 = parPtr->y;
	double x2 = popPtr->parPtr->x;
	double y2 = popPtr->parPtr->y;
	
	double xDiff = x1 - x2;
	double yDiff = y1 - y2;
	
	return sqrt(xDiff*xDiff + yDiff*yDiff);
}

void Population::updateRates()
{
	setEventRate(bitingEvent.get(), getBitingRate());
}

std::string Population::toString()
{
	return strprintf("p%u", id);
}

/*** BITING EVENT ***/

BitingEvent::BitingEvent(Population * popPtr, double rate, zppsim::rng_t & rng) :
	RateEvent(rate, 0.0, rng),
	popPtr(popPtr)
{
}

void BitingEvent::performEvent(zppsim::EventQueue & queue)
{
	popPtr->performBitingEvent();
}


/*** IMMIGRATION EVENT ***/

ImmigrationEvent::ImmigrationEvent(Population * popPtr, double rate, zppsim::rng_t & rng) :
	RateEvent(rate, 0.0, rng),
	popPtr(popPtr)
{
}

void ImmigrationEvent::performEvent(zppsim::EventQueue & queue)
{
	popPtr->performImmigrationEvent();
}
