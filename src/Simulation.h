#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__


#include "SimParameters.h"
#include "Host.h"
#include "Population.h"
#include "Strain.h"
#include "Gene.h"

#include "Database.hpp"
#include "zppsim_random.hpp"

#include <boost/array.hpp>
#include "zppdata_util.hpp"
#include "EventQueue.hpp"

class Simulation;

class Simulation
{
friend class Population;
friend class BitingEvent;
friend class Host;
public:
	Simulation(SimParameters & params, zppdata::Database & db);
	
	void run();
	void runUntil(double time);
	void runOneEvent();
	
	double getTime();
	double drawHostLifetime();
	
	Host * drawSourceHost(size_t dstPopId, size_t dstHostId);
	
	void addEvent(Event * event);
	void removeEvent(Event * event);
	void setEventTime(zppsim::Event * event, double time);
	
	StrainPtr getStrain(std::vector<GenePtr> const & strainGenes);
	StrainPtr generateRandomStrain();
	StrainPtr recombineStrains(StrainPtr const & s1, StrainPtr const & s2);
	
	bool verifyState();
private:
	SimParameters * parPtr;
	zppdata::Database * dbPtr;
	zppsim::rng_t rng;
	
	std::unique_ptr<zppsim::EventQueue> queuePtr;
	
	std::vector<std::unique_ptr<Population>> popPtrs;
	
	// Strain tracking: one strain object for each unique strain
	std::vector<StrainPtr> strains;
	std::unordered_map<StrainPtr, size_t> strainPtrToIndexMap;
	zppsim::unordered_map_bh<std::vector<GenePtr>, size_t> geneVecToStrainIndexMap;
	
	// Gene tracking: right now genes comprise a fixed pool, so no complicated
	// tracking to perform
	std::vector<GenePtr> genes;
};

#endif /* defined(__malariamodel__Simulation__) */
