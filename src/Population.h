//
//  Population.h
//  malariamodel
//
//  Created by Ed Baskerville on 4/21/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__Population__
#define __malariamodel__Population__

#include <iostream>
#include <vector>
#include "EventQueue.hpp"
#include "SimParameters.h"
#include "Host.h"

class Population;
class Simulation;
class SimParameters;

class BitingEvent : public zppsim::RateEvent
{
public:
    BitingEvent(Population * popPtr, double rate, zppsim::rng_t & rng);
    virtual void performEvent(zppsim::EventQueue & queue);
    
private:
    Population * popPtr;
};

class ImmigrationEvent : public zppsim::RateEvent
{
public:
    ImmigrationEvent(Population * popPtr, double rate, zppsim::rng_t & rng);
    virtual void performEvent(zppsim::EventQueue & queue);
private:
    Population * popPtr;
};

class Population
{
friend class Simulation;
friend class Host;
public:
    int64_t const id;
    
    Population(Simulation * simPtr, int64_t id);
    
    double bitingRate();
    int64_t size();
    
    Host * getHostAtIndex(int64_t hostIndex);
    void removeHost(Host * hostPtr);
    Host * createNewHost();
    
    double getTime();
    double getBitingRate();
    double getImmigrationRate();
    
    void addEvent(zppsim::Event * event);
    void removeEvent(zppsim::Event * event);
    void setEventTime(zppsim::Event * event, double time);
    void setEventRate(zppsim::RateEvent * event, double rate);
    
    void performBitingEvent();
    void performImmigrationEvent();
    
    double getDistance(Population * popPtr);
    
    void updateRates();
    void sampleHosts();
    
    std::string toString();
    
    void recordTransmission(Host & srcHost, Host & dstHost, Strain & strain);
    
private:
    Simulation * simPtr;
    zppsim::rng_t * rngPtr;
    PopulationParameters * parPtr;
    
    std::vector<std::unique_ptr<Host>> hosts;
    std::unordered_map<int64_t, int64_t> hostIdIndexMap;
    
    std::unique_ptr<BitingEvent> bitingEvent;
    std::unique_ptr<ImmigrationEvent> immigrationEvent;
    
    int64_t drawSourcePopulation();
    
    int64_t transmissionCount;
};

#endif /* defined(__malariamodel__Population__) */
