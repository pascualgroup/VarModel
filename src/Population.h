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
	BitingEvent(Population * popPtr, zppsim::rng_t & rng);
	virtual void performEvent(zppsim::EventQueue & queue);
	
private:
	Population * popPtr;
};

class Population
{
friend class Host;
public:
	size_t const id;
	
	Population(Simulation * simPtr, size_t id);
	void pushBackEvents(std::vector<zppsim::Event *> & eventVec);
	
	double bitingRate();
	
	void performBitingEvent();
private:
	Simulation * simPtr;
	PopulationParameters * parPtr;
	
	size_t nextHostId;
	
	std::vector<std::unique_ptr<Host>> hosts;
	
	std::unique_ptr<BitingEvent> bitingEvent;
};

#endif /* defined(__malariamodel__Population__) */
