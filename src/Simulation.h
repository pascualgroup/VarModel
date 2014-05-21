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

#define DECLARE_TABLE(tableName) std::unique_ptr<zppdata::DBTable> tableName ## TablePtr

class Simulation;

class RateUpdateEvent : public zppsim::PeriodicEvent
{
public:
	RateUpdateEvent(Simulation * simPtr, double initialTime, double period);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	Simulation * simPtr;
};

class HostStateSamplingEvent : public zppsim::PeriodicEvent
{
public:
	HostStateSamplingEvent(Simulation * simPtr, double initialTime, double period);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	Simulation * simPtr;
};

class Simulation
{
friend class Population;
friend class BitingEvent;
friend class Host;
public:
	Simulation(SimParameters * parPtr, zppdata::Database * dbPtr);
	
	void run();
	void runUntil(double time);
	void runOneEvent();
	
	double getTime();
	double drawHostLifetime();
	
	void addEvent(Event * event);
	void removeEvent(Event * event);
	void setEventTime(zppsim::Event * event, double time);
	void setEventRate(zppsim::RateEvent * event, double rate);
	
	double getSeasonality();
	double distanceWeightFunction(double d);
	
	GenePtr drawRandomGene();
	GenePtr mutateGene(GenePtr const & srcGene);
	
	StrainPtr getStrain(std::vector<GenePtr> const & strainGenes);
	StrainPtr generateRandomStrain();
	StrainPtr mutateStrain(StrainPtr & strain);
	StrainPtr recombineStrains(StrainPtr const & s1, StrainPtr const & s2);
	
	Host * drawDestinationHost(size_t srcPopId);
	
	void updateRates();
	void sampleHosts();
	
	void countTransmission();
	
	bool verifyState();
private:
	SimParameters * parPtr;
	zppdata::Database * dbPtr;
	zppsim::rng_t rng;
	
	std::unique_ptr<zppsim::EventQueue> queuePtr;
	RateUpdateEvent rateUpdateEvent;
	HostStateSamplingEvent hostStateSamplingEvent;
	
	size_t nextHostId;
	std::vector<std::unique_ptr<Population>> popPtrs;
	
	// Strain tracking: one strain object for each unique strain
	size_t nextStrainId;
	std::vector<StrainPtr> strains;
	std::unordered_map<StrainPtr, size_t> strainPtrToIndexMap;
	zppsim::unordered_map_bh<std::vector<GenePtr>, size_t> geneVecToStrainIndexMap;
	
	// Gene tracking: right now genes comprise a fixed pool, so no complicated
	// tracking to perform
	std::vector<GenePtr> genes;
	std::vector<std::discrete_distribution<>> mutationDistributions;
	
	size_t transmissionCount;
	
	// DECLARE_TABLE is a macro mapping e.g.,
	// DECLARE_TABLE(genes);
	// -->
	// std::unique_ptr<zppdata::DBTable> genesTablePtr;
	DECLARE_TABLE(genes);
	DECLARE_TABLE(strains);
	DECLARE_TABLE(hosts);
	
	DECLARE_TABLE(sampledHosts);
	DECLARE_TABLE(sampledHostInfections);
	DECLARE_TABLE(sampledHostImmunity);
	DECLARE_TABLE(sampledHostClinicalImmunity);
	
	DECLARE_TABLE(sampledTransmissions);
	DECLARE_TABLE(sampledTransmissionInfections);
	DECLARE_TABLE(sampledTransmissionImmunity);
	DECLARE_TABLE(sampledTransmissionClinicalImmunity);
	
	void initializeDatabaseTables();
	void commitDatabase();
};

#endif /* defined(__malariamodel__Simulation__) */
