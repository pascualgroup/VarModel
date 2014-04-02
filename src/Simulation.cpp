#include "Simulation.h"
#include "vecutils.h"

using namespace std;

/*** SIMULATION INITIALIZATION ***/

//Simulation::Simulation(SimParameters & params, Database & db)
//:
//	p(&params),
//	dbPtr(&db),
//	rng(p->randomSeed),
//	bitingEvent(this),
//	introductionEvent(this),
//	nextHostId(0)
//{
//	// Create hosts
//	hosts.reserve(p->initialPopulationSize);
//	gamma_distribution<> lifetimeDist(p->lifetimeShape, p->lifetimeMean / p->lifetimeShape);
//	for(uint32_t i = 0; i < p->initialPopulationSize; i++) {
//		hosts.push_back(unique_ptr<Host>(new Host(nextHostId++, 0.0, lifetimeDist(rng))));
//	}
//	
//	vector<Event *> events;
//	events.push_back((Event *)&bitingEvent);
//	events.push_back((Event *)&introductionEvent);
//	
//	samplerPtr = unique_ptr<EventSampler>(new EventSampler(events, rng));
//}

/*** MAIN EVENT LOOP ***/

//void Simulation::run()
//{
//	cerr << "Running" << endl;
//	
//	while(samplerPtr->getTime() < p->tEnd) {
//		Event * event;
//		size_t eventId = numeric_limits<size_t>::max();
//		double dt;
//		if(samplerPtr->performNextEvent(event, eventId, dt))
//		{
//			cerr << "Event happened!" << endl;
//			cerr << "t = " << samplerPtr->getTime() << endl;
//			cerr << ((SimulationEvent *)event)->toJsonString() << endl;
//		}
//		else
//		{
//			cerr << "No event happened." << endl;
//			assert(isinf(samplerPtr->getTime()));
//			assert(eventId == numeric_limits<size_t>::max());
//		}
//	}
//}

/*** BITING EVENT ***/

//std::vector<Event *> BitingEvent::performEvent(double time)
//{
//	return vector<Event *>(0);
//}
//
//double BitingEvent::getRate()
//{
//	return sim->hosts.size() * sim->p->bitingRate;
//}

/*** IMMIGRATION EVENT ***/

//string IntroductionEvent::toJsonString()
//{
//	return "{\"name\" : \"introduction\"}";
//}
//
//double IntroductionEvent::getRate()
//{
//	return sim->hosts.size() * sim->p->introductionRate;
//}
//
//std::vector<Event *> IntroductionEvent::performEvent(double time)
//{
//	return vector<Event *>(0);
//}
