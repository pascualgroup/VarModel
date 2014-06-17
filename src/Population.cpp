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

Population::Population(Simulation * simPtr, int64_t id) :
	id(id), simPtr(simPtr), rngPtr(&(simPtr->rng)),
	parPtr(&(simPtr->parPtr->populations[id]))
{
	// Create hosts
	hosts.reserve(parPtr->size);
	for(int64_t i = 0; i < parPtr->size; i++) {
		int64_t hostId = simPtr->nextHostId++;
		double lifetime = simPtr->drawHostLifetime();
		double birthTime = -uniform_real_distribution<>(0, lifetime)(*rngPtr);
		double deathTime = birthTime + lifetime;
		
		hosts.emplace_back(new Host(this, hostId, birthTime, deathTime, simPtr->hostsTablePtr.get()));
		hostIdIndexMap[hostId] = hosts.size() - 1;
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
	for(int64_t i = 0; i < parPtr->nInitialInfections; i++) {
		int64_t hostId = drawUniformIndex(simPtr->rng, hosts.size());
		StrainPtr strainPtr = simPtr->generateRandomStrain();
		
		hosts[hostId]->receiveInfection(strainPtr);
	}
}

int64_t Population::size()
{
	return hosts.size();
}

Host * Population::getHostAtIndex(int64_t hostIndex)
{
	return hosts[hostIndex].get();
}

void Population::removeHost(Host * hostPtr)
{
	int64_t index = hostIdIndexMap[hostPtr->id];
	hostIdIndexMap.erase(hostPtr->id);
	if(index < hosts.size() - 1) {
		hosts[index] = std::move(hosts.back());
		hostIdIndexMap[hosts[index]->id] = index;
	}
	hosts.pop_back();
	
	// In the future, need to update rates
//	updateRates();
	
	// For now, just create a new one so pop size doesn't change
	createNewHost();
}

Host * Population::createNewHost()
{
	double lifetime = simPtr->drawHostLifetime();
	double birthTime = getTime();
	double deathTime = birthTime + lifetime;
	
	int64_t hostId = simPtr->nextHostId++;
	hosts.emplace_back(new Host(this, hostId, birthTime, deathTime, simPtr->hostsTablePtr.get()));
	hostIdIndexMap[hostId] = hosts.size() - 1;
	
	// In the future, need to update rates
//	updateRates();
	
	return hosts.back().get();
}

double Population::getTime()
{
	return simPtr->getTime();
}

double Population::getBitingRate()
{
	BitingRate brObj = parPtr->bitingRate;
	double perHostBitingRate = brObj.mean + brObj.amplitude * simPtr->getSeasonality();
	return hosts.size() * perHostBitingRate;
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
//	cerr << simPtr->getTime() << ": biting event, src pop " << id << '\n';
	
	int64_t srcHostIndex = drawUniformIndex(*rngPtr, hosts.size());
	Host * srcHostPtr = hosts[srcHostIndex].get();
//	cerr << "src host: " << srcHostPtr->id << '\n';
	
	Host * dstHostPtr = simPtr->drawDestinationHost(id);
//	cerr << "dst pop, host: " << dstHostPtr->popPtr->id << ", " << dstHostPtr->id << '\n';
	
	srcHostPtr->transmitTo(*dstHostPtr);
}

void Population::performImmigrationEvent()
{
	cerr << getTime() << ": immigration event, pop " << id << '\n';
	int64_t hostIndex = drawUniformIndex(*rngPtr, hosts.size());
	StrainPtr strain = simPtr->generateRandomStrain();
	hosts[hostIndex]->receiveInfection(strain);
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

void Population::sampleHosts()
{
	vector<size_t> hostIndices = drawUniformIndices(
		*rngPtr, hosts.size(), size_t(parPtr->sampleSize), true
	);
	for(size_t index : hostIndices) {
		if(simPtr->sampledHostsTablePtr != nullptr) {
			DBRow row;
			row.setReal("time", getTime());
			row.setInteger("hostId", int64_t(hosts[index]->id));
			simPtr->sampledHostsTablePtr->insert(row);
		}
		
		hosts[index]->writeInfections(simPtr->sampledHostInfectionsTablePtr.get());
		hosts[index]->immunity.write(simPtr->sampledHostImmunityTablePtr.get());
		hosts[index]->clinicalImmunity.write(simPtr->sampledHostClinicalImmunityTablePtr.get());
	}
}

std::string Population::toString()
{
	return strprintf("p%u", id);
}

void Population::countTransmission()
{
	simPtr->countTransmission();
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
