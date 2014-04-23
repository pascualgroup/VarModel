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
	// Initialize empty event queue
	queuePtr = unique_ptr<EventQueue>(new EventQueue(rng));
	
	// Create gene pool
	genes.reserve(params.genePoolSize);
	for(size_t i = 0; i < params.genePoolSize; i++) {
		genes.emplace_back(new Gene(i));
	}
	
	// Create populations
	popPtrs.reserve(params.populations.size());
	for(size_t popId = 0; popId < params.populations.size(); popId++) {
		popPtrs.emplace_back(new Population(this, popId));
	}
	
	// Vector of initial events to populate event queue with
//	vector<Event *> initEvents;
//	for(auto & popPtr : popPtrs) {
//		popPtr->pushBackEvents(initEvents);
//	}
	cerr << "# events: " << queuePtr->size() << '\n';
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
}

double Simulation::getTime()
{
	return queuePtr->getTime();
}

double Simulation::drawHostLifetime()
{
	return parPtr->hostLifetimeDistribution.draw(rng);
}

Host * Simulation::drawSourceHost(size_t dstPopId, size_t dstHostId)
{
	size_t srcPopId;
	size_t nPops = popPtrs.size();
	
	if(nPops == 1) {
		srcPopId = 0;
	}
	else {
		std::vector<double> weights(nPops);
		for(size_t i = 0; i < nPops; i++) {
			size_t popSize = popPtrs[i]->size();
			if(i == dstPopId) {
				popSize -= 1;
			}
			weights[i] = parPtr->populations[dstPopId].contactWeight[i] * popSize;
		}
		srcPopId = sampleDiscreteLinearSearch(rng, weights);
	}
	
	size_t srcHostId;
	if(srcPopId == dstPopId) {
		srcHostId = drawUniformIndexExcept(
			rng, popPtrs[srcPopId]->size(), dstHostId
		);
		assert(srcHostId != dstHostId);
	}
	else {
		srcHostId = drawUniformIndex(rng, popPtrs[srcPopId]->size());
	}
	return popPtrs[srcPopId]->getHost(srcHostId);
}

void Simulation::addEvent(Event * event)
{
	queuePtr->addEvent(event);
}

void Simulation::removeEvent(Event * event)
{
	queuePtr->removeEvent(event);
}

void Simulation::setEventTime(zppsim::Event * event, double time)
{
	event->setTime(*queuePtr, time);
}

StrainPtr Simulation::generateRandomStrain()
{
	size_t strainSize = parPtr->strainSize;
	
	std::vector<GenePtr> strainGenes(strainSize);
	for(size_t i = 0; i < strainSize; i++) {
		size_t geneIndex = drawUniformIndex(rng, strainSize);
		strainGenes[i] = genes[geneIndex];
	}
	
	return getStrain(strainGenes);
}

StrainPtr Simulation::getStrain(std::vector<GenePtr> const & strainGenes)
{
	StrainPtr strainPtr;
	auto strainItr = geneVecToStrainIndexMap.find(strainGenes);
	if(strainItr == geneVecToStrainIndexMap.end()) {
		strains.emplace_back(new Strain(strainGenes));
		strainPtr = strains.back();
		geneVecToStrainIndexMap[strainGenes] = strains.size() - 1;
	}
	else {
		strainPtr = strains[strainItr->second];
	}
	return strainPtr;
}
