#ifndef __malariamodel__Simulation__
#define __malariamodel__Simulation__


#include "SimParameters.h"
#include "DatabaseTypes.h"
#include "CheckpointDatabaseTypes.h"
#include "Host.h"
#include "Population.h"
#include "Strain.h"
#include "Gene.h"
#include "DiscretizedDistribution.h"
#include "random"
#include "zppdb.hpp"
#include "zppsim_random.hpp"
#include <fstream>

#include "EventQueue.hpp"
#include <iterator>


class Simulation;

class CheckpointEvent : public zppsim::PeriodicEvent
{
public:
    CheckpointEvent(Simulation * simPtr, double initialTime, double period);
    virtual void performEvent(zppsim::EventQueue & queue);
private:
    Simulation * simPtr;
};

class VerificationEvent : public zppsim::PeriodicEvent
{
public:
    VerificationEvent(Simulation * simPtr, double initialTime, double period);
    virtual void performEvent(zppsim::EventQueue & queue);
private:
    Simulation * simPtr;
};

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

//add MDA events to simulate giving drugs to all hosts
class MDAEvent : public zppsim::PeriodicEvent
{
public:
    MDAEvent(Simulation * simPtr, double initialTime, double period);
    virtual void performEvent(zppsim::EventQueue & queue);
private:
    Simulation * simPtr;
};

class IRSEvent : public zppsim::OneTimeEvent
{
public:
    IRSEvent(Simulation * simPtr, double time);
    virtual void performEvent(zppsim::EventQueue & queue);
private:
    Simulation * simPtr;
};

class RemoveIRSEvent : public zppsim::OneTimeEvent
{
public:
    RemoveIRSEvent(Simulation * simPtr, double time);
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
friend class ImmuneHistory;
friend class CheckpointEvent;
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
	GenePtr drawRandomGeneExcept(int64_t geneId);
	GenePtr mutateGene(GenePtr const & srcGene, int64_t const source);
	//GenePtr mutateGene2(GenePtr const & srcGene);
    //int64_t recLociId(std::vector<int64_t> & recGeneAlleles);
    int64_t recLociId(std::vector<int64_t> & recGeneAlleles, std::vector<GenePtr> & searchSet);
    double parentsSimilarity(GenePtr const & pGene1, GenePtr const & pGene2, int64_t breakPoint);
    std::vector<GenePtr> ectopicRecomb(GenePtr const & pGene1, GenePtr const & pGene2, bool isConversion);
	
    GenePtr mutateMS(GenePtr const & srcMS);
    
	StrainPtr getStrain(std::vector<GenePtr> const & strainGenes);
	StrainPtr generateRandomStrain();
	StrainPtr generateRandomStrain(int64_t nNewGenes);
	std::vector<GenePtr> mutateStrain(StrainPtr & strain);
    std::vector<GenePtr> ectopicRecStrain(StrainPtr & strain);
	StrainPtr recombineStrains(StrainPtr const & s1, StrainPtr const & s2);
	GenePtr generateRandomMicrosat();
    GenePtr storeMicrosat(std::vector<int64_t> Alleles);
    GenePtr recombineMS(GenePtr const & ms1, GenePtr const & ms2);
	Host * drawDestinationHost(int64_t srcPopId);
	
	void updateRates();
	void sampleHosts();
    void MDA();
    void IRS();
    void RemoveIRS();
    
    void recordImmunity(Host & host, int64_t locusIndex, int64_t alleleId);
	void recordTransmission(Host & srcHost, Host & dstHost, std::vector<StrainPtr> & strains);
    void writeDuration(std::list<Infection>::iterator infectionItr);
    void writeEIR(double time, int64_t infectious);
    void writeFollowedHostInfection(std::list<Infection>::iterator infectionItr);
	void verifyState();
private:
	SimParameters * parPtr;
	zppdb::Database * dbPtr;
	zppsim::rng_t rng;
	
	DiscretizedDistribution hostLifetimeDist;
	
	// MAIN EVENT QUEUE
	std::unique_ptr<zppsim::EventQueue> queuePtr;
	
    CheckpointEvent checkpointEvent;
    VerificationEvent verificationEvent;
    
	RateUpdateEvent rateUpdateEvent;
	HostStateSamplingEvent hostStateSamplingEvent;
    MDAEvent mdaEvent;
    IRSEvent irsEvent;
    RemoveIRSEvent removeirsEvent;
    int64_t mdaCounts = 0;
	
	int64_t nextHostId;
	std::vector<std::unique_ptr<Population>> popPtrs;
	
	// Strain tracking: one strain object for each unique strain
	int64_t nextStrainId;
	std::vector<StrainPtr> strains;
	std::unordered_map<StrainPtr, int64_t> strainPtrToIndexMap;
	std::unordered_map<std::vector<GenePtr>, int64_t, HashGenePtrVec> geneVecToStrainIndexMap;
	
	// Gene tracking
	std::vector<GenePtr> genes;
	//std::vector<std::discrete_distribution<>> mutationDistributions; hqx change
    // get mutation weight distribution of each locus within a gene
    std::vector<double> mutationDistributions = parPtr->genes.mutationWeights.toDoubleVector();
	
    // Loci tracking
    size_t locusNumber = parPtr->genes.locusNumber;
    std::vector<int64_t> alleleNumber;
    
    // microsat tracking, if required
    // read microsat array from files generated by fastsimcoal
    std::vector<GenePtr> microsats;
    std::vector<std::vector<int64_t>> readMsArray(int year);
    std::vector<int64_t> microsatAlleles;
    size_t microsatNumber = parPtr->genes.microsatNumber;
    
	int64_t transmissionCount;
    int64_t mutationCount;
	
	// Database tables
	zppdb::Table<GeneRow> genesTable;
	zppdb::Table<LociRow> lociTable;
	zppdb::Table<StrainRow> strainsTable;
	zppdb::Table<HostRow> hostsTable;
    zppdb::Table<LociRow> microsatTable;
    zppdb::Table<AlleleImmunityRow> alleleImmunityTable;
	
	zppdb::Table<SampledHostRow> sampledHostsTable;
	zppdb::Table<InfectionRow> sampledHostInfectionTable;
	zppdb::Table<ImmunityRow> sampledHostImmunityTable;
    
	zppdb::Table<InfectionDurationRow> InfectionDurationTable;
    zppdb::Table<recordEIRRow> recordEIRTable;
	zppdb::Table<TransmissionRow> sampledTransmissionTable;
	zppdb::Table<TransmissionStrainRow> sampledTransmissionStrainTable;
	zppdb::Table<TransmissionInfectionRow> sampledTransmissionInfectionTable;
	zppdb::Table<TransmissionImmunityRow> sampledTransmissionImmunityTable;
    zppdb::Table<followedHostsRow> followedHostsTable;
	
	GenePtr createGene(std::vector<int64_t> Alleles,bool const functionality, int64_t const source);
	GenePtr createMicrosat(std::vector<int64_t> Alleles);
    void runMSSimCoal(size_t msSampleSize);
	void initializeDatabaseTables();
    
    void saveCheckpoint();
    void writeMetaToCheckpoint(Database & cpdb, Table<CheckpointMetaRow> & metaTable);
    void writeGenesAndStrainsToCheckpoint(Database & cpdb,
        Table<CheckpointAlleleCountRow> & alleleCountsTable,
        Table<GeneRow> & genesTable,
        Table<LociRow> & geneAllelesTable,
        Table<StrainRow> & strainsTable
    );
    void writeMicrosatsToCheckpoint(Database & cpdb,
        Table<CheckpointAlleleCountRow> & msAlleleCountsTable,
        Table<GeneRow> & msTable,
        Table<LociRow> & msAllelesTable
    );
    void writePopulationsToCheckpoint(Database & cpdb,
        Table<CheckpointHostRow> & hostsTable,
        Table<CheckpointInfectionRow> & infectionsTable,
        Table<CheckpointExpressionOrderRow> & expOrderTable,
        Table<CheckpointAlleleImmunityRow> & alleleImmunityTable,
        Table<CheckpointImmunityRow> & immunityTable
    );
    
    void initialize();
    
    void loadCheckpoint();
    void loadAlleleCounts(
        std::vector<CheckpointAlleleCountRow> & alleleCountRows,
        std::vector<int64_t> & alleleCounts,
        int64_t nLoci
    );
    void loadGenes(
        std::vector<GeneRow> & geneRows,
        std::vector<LociRow> & geneAlleleRows,
        std::vector<GenePtr> & genes,
        int64_t nLoci
    );
    void loadStrains(std::vector<StrainRow> & strainRows);
    void loadStrain(int64_t strainId, std::vector<GenePtr> strainGenes);
    
    void loadPopulations(
        Database & cpdb,
        double time,
        Table<CheckpointHostRow> & hostTable,
        Table<CheckpointInfectionRow> & infectionTable,
        Table<CheckpointExpressionOrderRow> & expOrderTable,
        Table<CheckpointAlleleImmunityRow> & alleleImmunityTable,
        Table<CheckpointImmunityRow> & immunityTable
    );
};

#endif /* defined(__malariamodel__Simulation__) */
