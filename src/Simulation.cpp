#include "Simulation.h"
#include "zppsim_util.hpp"

using namespace std;
using namespace zppdata;
using namespace zppsim;

Simulation::Simulation(SimParameters & params, Database & db)
:
	parPtr(&params),
	dbPtr(&db),
	rng(parPtr->randomSeed)
{
	// Create populations
	popPtrs.reserve(parPtr->populations.size());
	for(size_t popId = 0; popId < parPtr->populations.size(); popId++) {
		popPtrs.emplace_back(new Population(this, popId));
	}
	
	// Vector of initial events to populate event queue with
	vector<Event *> initEvents;
	for(auto & popPtr : popPtrs) {
		popPtr->pushBackEvents(initEvents);
	}
	cerr << "# events: " << initEvents.size() << '\n';
	
	// Initialize event queue with list of initial events
	queuePtr = unique_ptr<EventQueue>(new EventQueue(rng, initEvents, this));
}

void Simulation::run()
{
	cout << "Starting run..." << '\n';
	runUntil(parPtr->tEnd);
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
	
	/*// Need to handle death outside event queue to avoid premature deletion
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
	}*/
}

double Simulation::getTime()
{
	return queuePtr->getTime();
}

/*uint64_t Simulation::getNextHostId()
{
	return nextHostId++;
}*/

double Simulation::drawHostLifetime()
{
	return parPtr->hostLifetimeDistribution.draw(rng);
}

/*double Simulation::totalBitingRate(uint32_t popId)
{
	return hosts[popId].size() * parPtr->bitingRate[popId];
}

double Simulation::totalIntroductionRate(uint32_t popId)
{
	return hosts[popId].size() * parPtr->introductionRate[popId];
}*/

/*void Simulation::bite(uint32_t dstPopId)
{
	cerr << "Biting " << dstPopId << " at " << queuePtr->getTime() << endl;
	
	// Choose source population ID
	uint32_t srcPopId;
	if(parPtr->nPopulations == 1) {
		assert(dstPopId == 0);
		srcPopId = dstPopId;
	}
	else {
		vector<double> popWeights;
		for(uint32_t j = 0; j < parPtr->nPopulations; j++) {
			double weight = parPtr->contactWeight[dstPopId][j];
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
}*/


/*** BITING EVENTS ***/

/*BitingEvent::BitingEvent(Simulation * simPtr, uint32_t popId, rng_t & rng) :
	RateEvent(simPtr->totalBitingRate(popId), 0.0, rng),
	simPtr(simPtr), popId(popId)
{
}

void BitingEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->bite(popId);
}*/

/*** INTRODUCTION EVENTS ***/

/*IntroductionEvent::IntroductionEvent(Simulation * simPtr, uint32_t popId, rng_t & rng) :
	RateEvent(simPtr->totalIntroductionRate(popId), 0.0, rng),
	simPtr(simPtr), popId(popId)
{
}

void IntroductionEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->introduce(popId);
}
*/