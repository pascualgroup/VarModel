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
	// Create gene pool
	genes.reserve(params.genePoolSize);
	for(size_t i = 0; i < params.genePoolSize; i++) {
		genes.emplace_back(new Gene());
	}
	
	// Create populations
	popPtrs.reserve(params.populations.size());
	for(size_t popId = 0; popId < params.populations.size(); popId++) {
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
}

double Simulation::getTime()
{
	return queuePtr->getTime();
}

double Simulation::drawHostLifetime()
{
	return parPtr->hostLifetimeDistribution.draw(rng);
}

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
}*/
