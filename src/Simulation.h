#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__

#include <iostream>

#include <boost/array.hpp>

#include "random.h"
#include "strutils.h"
#include "SimParameters.h"
#include "Host.h"

//#include "EventSampler.h"

/*class Simulation;

class SimulationEvent : public RateEvent
{
public:
	SimulationEvent(Simulation * simPtr, double rate, double initTime, rng_t & rng);
protected:
	Simulation * simPtr;
};

class BitingEvent : public SimulationEvent
{
public:
	BitingEvent(Simulation * simPtr);
	virtual void performEvent(EventSampler & sampler);
};

class IntroductionEvent : public SimulationEvent
{
public:
	IntroductionEvent(Simulation * sim);
	virtual void performEvent(EventSampler & sampler);
};

class Simulation
{
//friend class BitingEvent;
//friend class IntroductionEvent;
public:
	Simulation(SimParameters & params, Database & db);
	void run();
private:
	SimParameters * p;
//	Database * dbPtr;
	rng_t rng;
	
	std::vector<std::unique_ptr<Host>> hosts;
	
//	BitingEvent bitingEvent;
//	IntroductionEvent introductionEvent;
	
	uint64_t nextHostId;
	
	std::unique_ptr<EventSampler> samplerPtr;
};*/

#endif /* defined(__malariamodel__Simulation__) */
