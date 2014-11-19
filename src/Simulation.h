#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__


#include "SimParameters.h"
#include "Host.h"
#include "Population.h"
#include "Strain.h"
#include "Gene.h"
#include "DiscretizedDistribution.h"

#include "zppdb.hpp"
#include "zppsim_random.hpp"

#include "EventQueue.hpp"


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

class HashGenePtrVec
{
public:
	size_t operator()(std::vector<GenePtr> const & genePtrVec) const
	{
		if(genePtrVec.size() == 0) {
			return 0;
		}
		size_t hashVal = _hash(genePtrVec[0]);
		for(size_t i = 0; i < genePtrVec.size(); i++) {
			hashVal ^= _hash(genePtrVec[i]) + 0x9e3779b9 + (hashVal << 6) + (hashVal >> 2);
		}
		return hashVal;
	}
private:
	std::hash<GenePtr> _hash;
};

class Simulation
{
friend class Population;
friend class BitingEvent;
friend class Host;
public:
	Simulation(SimParameters * parPtr, zppdb::Database * dbPtr);
	
	void run();
	void runUntil(double time);
	void runOneEvent();
	
	double getTime();
	double drawHostLifetime();
	
	void addEvent(zppsim::Event * event);
	void removeEvent(zppsim::Event * event);
	void setEventTime(zppsim::Event * event, double time);
	void setEventRate(zppsim::RateEvent * event, double rate);
	
	double distanceWeightFunction(double d);
	
	GenePtr drawRandomGene();
	GenePtr mutateGene(GenePtr const & srcGene);
	
	StrainPtr getStrain(std::vector<GenePtr> const & strainGenes);
	StrainPtr generateRandomStrain();
	StrainPtr mutateStrain(StrainPtr & strain);
	StrainPtr recombineStrains(StrainPtr const & s1, StrainPtr const & s2);
	
	Host * drawDestinationHost(int64_t srcPopId);
	
	void updateRates();
	void sampleHosts();
	
	void countTransmission();
	
	bool verifyState();
private:
	SimParameters * parPtr;
	zppdb::Database * dbPtr;
	zppsim::rng_t rng;
	
	DiscretizedDistribution hostLifetimeDist;
	
	std::unique_ptr<zppsim::EventQueue> queuePtr;
	RateUpdateEvent rateUpdateEvent;
	HostStateSamplingEvent hostStateSamplingEvent;
	
	int64_t nextHostId;
	std::vector<std::unique_ptr<Population>> popPtrs;
	
	// Strain tracking: one strain object for each unique strain
	int64_t nextStrainId;
	std::vector<StrainPtr> strains;
	std::unordered_map<StrainPtr, int64_t> strainPtrToIndexMap;
	std::unordered_map<std::vector<GenePtr>, int64_t, HashGenePtrVec> geneVecToStrainIndexMap;
	
	// Gene tracking: right now genes comprise a fixed pool, so no complicated
	// tracking to perform
	std::vector<GenePtr> genes;
	std::vector<std::discrete_distribution<>> mutationDistributions;
	
	int64_t transmissionCount;
	
	// Database tables
	zppdb::Table<GeneRow> genesTable;
	zppdb::Table<StrainRow> strainsTable;
	zppdb::Table<HostRow> hostsTable;
	
	zppdb::Table<SampledHostRow> sampledHostsTable;
	zppdb::Table<InfectionRow> sampledHostInfectionTable;
	zppdb::Table<ImmunityRow> sampledHostImmunityTable;
	zppdb::Table<ImmunityRow> sampledHostClinicalImmunityTable;
	
	zppdb::Table<SampledTransmissionRow> sampledTransmissionTable;
	zppdb::Table<InfectionRow> sampledTransmissionInfectionTable;
	zppdb::Table<ImmunityRow> sampledTransmissionImmunityTable;
	zppdb::Table<ImmunityRow> sampledTransmissionClinicalImmunityTable;
	
	void initializeDatabaseTables();
};

#endif /* defined(__malariamodel__Simulation__) */
