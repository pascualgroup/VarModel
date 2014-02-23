#include "Simulation.h"
#include "vecutils.h"

using namespace std;

Simulation::Simulation(SimParameters & params, Database & db)
:
	p(&params),
	dbPtr(&db),
	rng(p->randomSeed),
	hosts(p->initialPopulationSize),
	birthEvent(this),
	deathEvent(this),
	bitingEvent(this),
	immigrationEvent(this)
{
	// Create hosts
	for(uint32_t i = 0; i < p->initialPopulationSize; i++) {
		hosts[i] = unique_ptr<Host>(new Host());
	}
	
	vector<Event *> events;
	events.push_back((Event *)&birthEvent);
	events.push_back((Event *)&deathEvent);
	events.push_back((Event *)&bitingEvent);
	events.push_back((Event *)&immigrationEvent);
	
	samplerPtr = unique_ptr<EventSampler>(new EventSampler(events, rng));
}



double BirthEvent::getRate()
{
	return 0.0;
}

std::vector<Event *> BirthEvent::performEvent(double time)
{
	return vector<Event *>(0);
}

double DeathEvent::getRate()
{
	return 0.0;
}

std::vector<Event *> DeathEvent::performEvent(double time)
{
	return vector<Event *>(0);
}

double BitingEvent::getRate()
{
	return 0.0;
}

std::vector<Event *> BitingEvent::performEvent(double time)
{
	return vector<Event *>(0);
}

double ImmigrationEvent::getRate()
{
	return 0.0;
}

std::vector<Event *> ImmigrationEvent::performEvent(double time)
{
	return vector<Event *>(0);
}

void Simulation::run()
{
	cerr << "Running" << endl;
	
	while(samplerPtr->getTime() < p->tEnd) {
		Event * event;
		size_t eventId = numeric_limits<size_t>::max();
		double dt;
		if(samplerPtr->performNextEvent(event, eventId, dt))
		{
			cerr << "Event happened!" << endl;
			cerr << "t = " << samplerPtr->getTime() << endl;
			cerr << "event id " << eventId << endl;
		}
		else
		{
			cerr << "No event happened." << endl;
			assert(isinf(samplerPtr->getTime()));
			assert(eventId == numeric_limits<size_t>::max());
		}
	}
}
