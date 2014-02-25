#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__

#include <iostream>

#include <boost/array.hpp>

#include "Database.h"
#include "random.h"
#include "strutils.h"
#include "SimParameters.h"
#include "Host.h"

#include "EventSampler.h"

class Simulation;

class SimulationEvent : public Event
{
public:
	SimulationEvent(Simulation * sim) : sim(sim) {}
	virtual std::string toJsonString() = 0;
protected:
	Simulation * sim;
};

class BitingEvent : public SimulationEvent
{
public:
	BitingEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual std::string toJsonString();
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class IntroductionEvent : public SimulationEvent
{
public:
	IntroductionEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual std::string toJsonString();
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class Simulation
{
friend class BirthEvent;
friend class DeathEvent;
friend class BitingEvent;
friend class IntroductionEvent;
public:
	Simulation(SimParameters & params, Database & db);
	void run();
private:
	SimParameters * p;
	Database * dbPtr;
	rng_t rng;
	
	std::vector<std::unique_ptr<Host>> hosts;
	
	BitingEvent bitingEvent;
	IntroductionEvent introductionEvent;
	
	uint64_t nextHostId;
	
	std::unique_ptr<EventSampler> samplerPtr;
};

#endif /* defined(__malariamodel__Simulation__) */
