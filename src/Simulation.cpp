#include "Simulation.h"
#include "zppsim_util.hpp"

using namespace std;
using namespace zppdata;
using namespace zppsim;

Simulation::Simulation(SimParameters & params, Database & db)
:
	p(&params),
	dbPtr(&db),
	rng(p->randomSeed),
	hosts(p->nPopulations),
	hostToDie(nullptr)
{
	// Create hosts
	for(uint32_t popId = 0; popId < p->nPopulations; popId++) {
		hosts[popId] = vector<unique_ptr<Host>>(p->populationSize[popId]);
		for(uint32_t i = 0; i < p->populationSize[popId]; i++) {
			uint64_t hostId = nextHostId++;
			hosts[popId][i] = unique_ptr<Host>(new Host(hostId, drawHostLifetime()));
			hostIndexMap[hostId] = std::make_pair(popId, i);
		}
	}
	
	// Vector of initial events to populate event queue with
	vector<Event *> initEvents;
	
	// Biting and introduction events
	for(uint32_t popId = 0; popId < p->nPopulations; popId++) {
		bitingEvents.emplace_back(new BitingEvent(this, popId, rng));
		initEvents.push_back((Event *)bitingEvents[popId].get());
		
		introductionEvents.emplace_back(new IntroductionEvent(this, popId, rng));
		initEvents.push_back((Event *)introductionEvents[popId].get());
	}
	
	// Death events
	for(uint32_t popId = 0; popId < p->nPopulations; popId++) {
		for(uint32_t i = 0; i < p->populationSize[popId]; i++) {
			initEvents.push_back((Event *)hosts[popId][i]->deathEvent.get());
		}
	}
	
	// Initialize event queue with list of initial events
	queuePtr = unique_ptr<EventQueue>(new EventQueue(rng, initEvents, this));
}

void Simulation::run()
{
	cout << "Starting run..." << '\n';
	runUntil(p->tEnd);
	cout << "Total event count: " << queuePtr->getEventCount() << endl;
}

void Simulation::runUntil(double time)
{
	while(queuePtr->getNextTime() <= time) {
		runOneEvent();
	}
}

void Simulation::runOneEvent()
{
	Event * event;
	double dt;
	queuePtr->performNextEvent(event, dt);
	double time = queuePtr->getTime();
	
	// Need to handle death outside event queue to avoid premature deletion
	// Shared pointers would avoid this, but this is the only place where they're needed
	// so the overhead was deemed not worth it.
	if(hostToDie != nullptr) {
		uint64_t oldHostId = hostToDie->id;
		cerr << "Killing " << oldHostId << " at " << time << endl;
		
		auto index = hostIndexMap[oldHostId];
		size_t popId = index.first;
		size_t hostIndex = index.second;
		
		// For now, just replace with a new host
		uint64_t newHostId = nextHostId++;
		hosts[popId][hostIndex] = unique_ptr<Host>(new Host(newHostId, time + drawHostLifetime()));
		queuePtr->addEvent((Event *)hosts[popId][hostIndex]->deathEvent.get());
		
		hostIndexMap.erase(oldHostId);
		hostIndexMap[newHostId] = index;
		
		hostToDie = nullptr;
		
		// TODO: when no longer keeping population size, need to adjust biting & introduction rate here.
	}
}

uint64_t Simulation::getNextHostId()
{
	return nextHostId++;
}

double Simulation::drawHostLifetime()
{
	uniform_real_distribution<double> ud(0.0, 20.0);
	return ud(rng);
}

double Simulation::totalBitingRate(uint32_t popId)
{
	return hosts[popId].size() * p->bitingRate[popId];
}

double Simulation::totalIntroductionRate(uint32_t popId)
{
	return hosts[popId].size() * p->introductionRate[popId];
}

void Simulation::bite(uint32_t dstPopId)
{
	cerr << "Biting " << dstPopId << " at " << queuePtr->getTime() << endl;
	
	// Choose source population ID
	uint32_t srcPopId;
	if(p->nPopulations == 1) {
		assert(dstPopId == 0);
		srcPopId = dstPopId;
	}
	else {
		vector<double> popWeights;
		for(uint32_t j = 0; j < p->nPopulations; j++) {
			double weight = p->contactWeight[dstPopId][j];
			if(j == dstPopId) {
				popWeights.push_back((hosts[j].size() - 1) * weight);
			}
			else {
				popWeights.push_back(hosts[j].size() * weight);
			}
			cerr << "weight for pop " << j << " = " << popWeights[j] << endl;
		}
		
		srcPopId = (uint32_t)sampleDiscreteLinearSearch(rng, popWeights);
		cerr << "srcPopId: " << srcPopId << endl;
	}
	
	// Choose hosts
	size_t dstHostIndex = drawUniformIndex(rng, hosts[dstPopId].size());
	size_t srcHostIndex;
	if(srcPopId == dstPopId) {
		srcHostIndex = drawUniformIndexExcept(rng, hosts[dstPopId].size(), dstHostIndex);
	}
	else {
		srcHostIndex = drawUniformIndex(rng, hosts[srcPopId].size());
	}
	if(srcHostIndex == dstHostIndex) {
		assert(srcPopId != dstPopId);
	}
	
	hosts[srcPopId][srcHostIndex]->transmitTo(*hosts[dstPopId][dstHostIndex]);
	
	cerr << "srcHostIndex, dstHostIndex = " << srcHostIndex << ", " << dstHostIndex << endl;
}

void Simulation::introduce(uint32_t popId)
{
	cerr << "Introducing " << popId << " at " << queuePtr->getTime() << endl;
}

void Simulation::markForDeath(Host & host)
{
	cerr << "Host " << host.id << " marked for death" << endl;
	hostToDie = &host;
}


/*** BITING EVENTS ***/

BitingEvent::BitingEvent(Simulation * simPtr, uint32_t popId, rng_t & rng) :
	RateEvent(simPtr->totalBitingRate(popId), 0.0, rng),
	simPtr(simPtr), popId(popId)
{
}

void BitingEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->bite(popId);
}

/*** INTRODUCTION EVENTS ***/

IntroductionEvent::IntroductionEvent(Simulation * simPtr, uint32_t popId, rng_t & rng) :
	RateEvent(simPtr->totalIntroductionRate(popId), 0.0, rng),
	simPtr(simPtr), popId(popId)
{
}

void IntroductionEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->introduce(popId);
}
