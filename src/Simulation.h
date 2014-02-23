#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__

#include <iostream>

#include <boost/array.hpp>

#include "Database.h"
#include "types.h"
#include "strutils.h"
#include "SimParameters.h"
#include "Host.h"

#include "EventSampler.h"

using namespace std;

class Simulation;

class SimulationEvent : public Event
{
public:
	SimulationEvent(Simulation * sim) : sim(sim) {}
protected:
	Simulation * sim;
};

class BirthEvent : public SimulationEvent
{
public:
	BirthEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class DeathEvent : public SimulationEvent
{
public:
	DeathEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class BitingEvent : public SimulationEvent
{
public:
	BitingEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class ImmigrationEvent : public SimulationEvent
{
public:
	ImmigrationEvent(Simulation * sim) : SimulationEvent(sim) {}
	
	virtual double getRate();
	virtual std::vector<Event *> performEvent(double time);
};

class Simulation
{
friend class BirthEvent;
friend class DeathEvent;
friend class BitingEvent;
friend class ImmigrationEvent;
public:
	Simulation(SimParameters & params, Database & db);
	void run();
private:
	SimParameters * p;
	Database * dbPtr;
	rng_t rng;
	
	std::vector<std::unique_ptr<Host>> hosts;
	
	BirthEvent birthEvent;
	DeathEvent deathEvent;
	BitingEvent bitingEvent;
	ImmigrationEvent immigrationEvent;
	
	std::unique_ptr<EventSampler> samplerPtr;
};

#endif /* defined(__malariamodel__Simulation__) */
