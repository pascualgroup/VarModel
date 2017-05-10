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
#include <sstream>
#include <algorithm>

using namespace std;
using namespace zppsim;

Population::Population(Simulation * simPtr, int64_t id) :
	id(id), simPtr(simPtr), rngPtr(&(simPtr->rng)),
	parPtr(&(simPtr->parPtr->populations[id])),
	transmissionCount(0)
{
	// Create hosts
	hosts.reserve(parPtr->size);
	for(int64_t i = 0; i < parPtr->size; i++) {
		int64_t hostId = simPtr->nextHostId++;
		//double lifetime = simPtr->drawHostLifetime();
        double lifetime = exponential_distribution<>(1.0/10800.0)(*rngPtr);
        if (lifetime > 28800.0) lifetime = 28800.0;
        //cout<<lifetime<<endl;
		double birthTime = -uniform_real_distribution<>(0, lifetime)(*rngPtr);
		double deathTime = birthTime + lifetime;
		
		bool writeToDatabase = simPtr->parPtr->outputHosts;
		hosts.emplace_back(new Host(
			this, hostId, birthTime, deathTime,
			writeToDatabase,
			*(simPtr->dbPtr),
			simPtr->hostsTable
		));
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
	
	// Create initial infections, with microsatellites as well
    if (simPtr->parPtr->genes.includeMicrosat) {
        size_t msSampleSize = (int)(parPtr->nInitialInfections + (parPtr->immigrationRate * 360.0))*1.5;
        simPtr->runMSSimCoal(msSampleSize);
        tempMS = simPtr->readMsArray(yearTrack);
        for(int64_t i = 0; i < parPtr->nInitialInfections; i++) {
            int64_t hostId = drawUniformIndex(simPtr->rng, hosts.size());
            StrainPtr strainPtr = simPtr->generateRandomStrain();
            GenePtr msPtr = simPtr->storeMicrosat(tempMS[i]);
            hosts[hostId]->receiveInfection(strainPtr,msPtr);
        }
        immigrationCount = parPtr->nInitialInfections;
    }else{
        for(int64_t i = 0; i < parPtr->nInitialInfections; i++) {
		int64_t hostId = drawUniformIndex(simPtr->rng, hosts.size());
		StrainPtr strainPtr = simPtr->generateRandomStrain();
		hosts[hostId]->receiveInfection(strainPtr);
        }
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
	//double lifetime = simPtr->drawHostLifetime();
    double lifetime = exponential_distribution<>(1.0/10800.0)(*rngPtr);
    if (lifetime > 28800.0) lifetime = 28800.0;
	double birthTime = getTime();
	double deathTime = birthTime + lifetime;
	
	int64_t hostId = simPtr->nextHostId++;
	bool writeToDatabase = simPtr->parPtr->outputHosts;
	hosts.emplace_back(new Host(
		this, hostId, birthTime, deathTime,
		writeToDatabase, *(simPtr->dbPtr), simPtr->hostsTable
	));
	hostIdIndexMap[hostId] = hosts.size() - 1;
	
	// In the future, need to update rates
//	updateRates();
    if (simPtr->parPtr->following.includeHostFollowing) {
        double t = getTime();
        if ((t > simPtr->parPtr->burnIn) && (numberOfHostsFollowed < simPtr->parPtr->following.HostNumber))
        {
            hosts.back()->toTrack = true;
            numberOfHostsFollowed += 1;
            cout << "following hosts " <<hostId<<endl;
        }
    }
	return hosts.back().get();
}

double Population::getTime()
{
	return simPtr->getTime();
}

double Population::getBitingRate()
{
    double t = getTime();
    double perHostBitingRate;
    double bt = parPtr->bitingRate.mean;
	if (monthlyBitingRateDistribution.size()==12) {
        perHostBitingRate = monthlyBitingRateDistribution[(int)(t/30)%12]*bt;
    }else{
        perHostBitingRate = evaluateSinusoid(parPtr->bitingRate, t);
    }
    if ((simPtr->parPtr->intervention.includeIntervention) &&
        (t > simPtr->parPtr->intervention.TimeStart) &&
        (t <= (simPtr->parPtr->intervention.TimeStart + simPtr->parPtr->intervention.duration))) {
        perHostBitingRate = perHostBitingRate * simPtr->parPtr->intervention.amplitude;
        //cout<<"reduced transmission"<<endl;
    }
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
    bool NoMDAflag = true;
    if ((simPtr->parPtr->MDA.includeMDA) &&
    (dstHostPtr->MDAEndTime>(getTime()+simPtr->parPtr->tLiverStage)))
    {
        //express before the MDA is over, then do not perform biting
        bernoulli_distribution flipCoin(1-simPtr->parPtr->MDA.strainFailRate);
        if(flipCoin(*rngPtr)){
            NoMDAflag = false;
        }
    }
    
    if (NoMDAflag) {
        if (simPtr->parPtr->genes.includeMicrosat) {
            srcHostPtr->transmitMSTo(*dstHostPtr);
        }else{
            srcHostPtr->transmitTo(*dstHostPtr);
        }
    }
}

//disable immigration first
void Population::performImmigrationEvent()
{
//	cerr << getTime() << ": immigration event, pop " << id << '\n';
	int64_t hostIndex = drawUniformIndex(*rngPtr, hosts.size());

    
    bool NoMDAflag = true;
    if ((simPtr->parPtr->MDA.includeMDA) &&
        (hosts[hostIndex]->MDAEndTime>(getTime()+simPtr->parPtr->tLiverStage)))
    {
        //express before the MDA is over, then do not perform biting
        bernoulli_distribution flipCoin(1-simPtr->parPtr->MDA.strainFailRate);
        if(flipCoin(*rngPtr)){
            NoMDAflag = false;
        }
    }
    
    
    if(NoMDAflag){
        bernoulli_distribution flipCoin(parPtr->pImmigrationIncludesNewGenes);
        bool includesNewGenes = flipCoin(*rngPtr);
	
        StrainPtr strain;
        if(includesNewGenes) {
            //		cerr << "Generating strain with new genes" << endl;
            strain = simPtr->generateRandomStrain(parPtr->nImmigrationNewGenes);
        }
        else {
            //		cerr << "Generating strain with all old genes" << endl;
            strain = simPtr->generateRandomStrain();
        }
        if (simPtr->parPtr->genes.includeMicrosat) {
            int currentYear = floor(getTime()/360.0);
            if (currentYear > yearTrack) {
                tempMS = simPtr->readMsArray(currentYear);
                yearTrack = currentYear;
                immigrationCount = 0;
            
            }
            GenePtr ms = simPtr->storeMicrosat(tempMS[immigrationCount]);
            hosts[hostIndex]->receiveInfection(strain,ms);
            immigrationCount++;
        }else{
            hosts[hostIndex]->receiveInfection(strain);
        }
        }
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
	Database * dbPtr = simPtr->dbPtr;
	
	/*vector<size_t> hostIndices = drawUniformIndices(
		*rngPtr, hosts.size(), size_t(parPtr->sampleSize), true
	);*/
    //change the samplingHosts to sample enough infected according to sampleSize
    vector<size_t> hostIndices = drawUniformIndices(
         *rngPtr, hosts.size(), hosts.size(), true);
    size_t count = 0;
    size_t moi1count = 0;
    size_t sampledSize = 0;
	for(size_t index : hostIndices) {
		//SampledHostRow row;
		//row.time = getTime();
        //this line records every sampled host
		//row.hostId = hosts[index]->id;
        //dbPtr->insert(simPtr->sampledHostsTable, row);
        //instead, in this version record only how many hosts were sampled before reaching the infected sampleSize
        sampledSize += 1;
        int64_t numActiveInfection = hosts[index]->getActiveInfectionCount();
		if (numActiveInfection>0) {
            hosts[index]->writeInfections(*dbPtr, simPtr->sampledHostInfectionTable, simPtr->strainsTable,simPtr->genesTable,simPtr->lociTable);
            count += 1;
            if (numActiveInfection ==1) {
                moi1count += 1;
            }
		//hosts[index]->immunity.write(*dbPtr, simPtr->sampledHostImmunityTable);
		//hosts[index]->clinicalImmunity.write(*dbPtr, simPtr->sampledHostClinicalImmunityTable);
        }
        if (parPtr->moi1) {
            if (moi1count == size_t(parPtr->sampleSize)) {
                break;
            }
        }else{
            if (count == size_t(parPtr->sampleSize)) {
                break;
            }
        }
       

    }
    //cout<<"sampled count is "<<count<<endl;
    SampledHostRow row;
    row.time = getTime();
    row.sampledNumber = sampledSize;
    dbPtr->insert(simPtr->sampledHostsTable, row);
}

void Population::executeMDA(double time)
{
    std::binomial_distribution<size_t> binoDist(hosts.size(),1-(simPtr->parPtr->MDA.hostFailRate));
    size_t totalHosts = binoDist(*rngPtr);
    vector<size_t> hostIndices = drawUniformIndices(
    *rngPtr, hosts.size(), totalHosts, true);
    for (size_t index : hostIndices) {
        //only give MDA to hosts whose age is older than 3 months
        if ((time - simPtr->parPtr->MDA.drugEffDuration - hosts[index]->birthTime)>90) {
            hosts[index]->MDAEndTime = time;
            hosts[index]->MDAClearInfection();
        }
    }
}


std::string Population::toString()
{
	stringstream ss;
	ss << "p" << id;
	return ss.str();
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
