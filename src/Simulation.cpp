#include "Simulation.h"
#include "zppsim_util.hpp"

// 100-millisecond delay between database commit retries
#define DB_RETRY_DELAY 100000

// 10-second maximum database retry before failure
#define DB_TIMEOUT 10000000

#define INITIALIZE_TABLE(tableName, columns) \
if(parPtr->dbTablesEnabled[#tableName]) { \
	tableName ## TablePtr = unique_ptr<DBTable>(new DBTable( \
		dbPtr, \
		#tableName, \
		columns \
	)); \
	(tableName ## TablePtr)->create(); \
}

using namespace std;
using namespace zppdata;
using namespace zppsim;

static float elapsed(clock_t clockStart, clock_t clockEnd)
{
	return float(clockEnd - clockStart) / CLOCKS_PER_SEC;
}

template<typename T>
T getEntry(vector<T> vals, size_t index, size_t size)
{
	if(vals.size() == 1) {
		return vals[0];
	}
	else {
		assert(vals.size() == size);
		return vals[index];
	}
}

Simulation::Simulation(SimParameters * parPtr, Database * dbPtr) :
	parPtr(parPtr),
	dbPtr(dbPtr),
	rng(parPtr->randomSeed),
	queuePtr(new EventQueue(rng)),
	rateUpdateEvent(this, 0.0, parPtr->seasonalUpdateEvery),
	hostStateSamplingEvent(this, 0.0, parPtr->sampleHostsEvery),
	nextHostId(0),
	nextStrainId(0),
	transmissionCount(0)
{
	// Construct transition probability distributions for genes
	if(parPtr->genes.mutationWeights.size() > 1) {
		assert(parPtr->genes.mutationWeights.size() == parPtr->genePoolSize);
		for(int64_t i = 0; i < parPtr->genePoolSize; i++) {
			assert(parPtr->genes.mutationWeights[i].size() == parPtr->genePoolSize);
			mutationDistributions.emplace_back(
				parPtr->genes.mutationWeights[i].begin(),
				parPtr->genes.mutationWeights[i].end()
			);
		}
	}
	
	dbPtr->beginTransaction();
	
	initializeDatabaseTables();
	
	queuePtr->addEvent(&rateUpdateEvent);
	queuePtr->addEvent(&hostStateSamplingEvent);
	
	// Create gene pool
	genes.reserve(parPtr->genePoolSize);
	for(int64_t i = 0; i < parPtr->genePoolSize; i++) {
		double transmissibility = getEntry(
			parPtr->genes.transmissibility, i, parPtr->genePoolSize
		);
		double immunityLossRate = getEntry(
			parPtr->genes.immunityLossRate, i, parPtr->genePoolSize
		);
		double clinicalImmunityLossRate = getEntry(
			parPtr->genes.clinicalImmunityLossRate, i, parPtr->genePoolSize
		);
		
		genes.emplace_back(new Gene(
			i,
			transmissibility,
			immunityLossRate,
			clinicalImmunityLossRate,
			genesTablePtr.get()
		));
	}
	
	// Create populations
	popPtrs.reserve(parPtr->populations.size());
	for(int64_t popId = 0; popId < parPtr->populations.size(); popId++) {
		popPtrs.emplace_back(new Population(this, popId));
	}
	
	dbPtr->commitWithRetry(DB_RETRY_DELAY, DB_TIMEOUT, cerr);
	
	cerr << "# events: " << queuePtr->size() << '\n';
}

void Simulation::initializeDatabaseTables()
{
	// INITIALIZE_TABLE is a macro that checks if the table is enabled;
	// and if it is creates a table object and calls xTablePtr->create(),
	// and otherwise leaves xTablePtr == nullptr.
	
	// GENES
	vector<DBColumn> genesColumns = {
		{"geneId", DBType::INTEGER},
		{"transmissibility", DBType::REAL},
		{"immunityLossRate", DBType::REAL},
		{"clinicalImmunityLossRate", DBType::REAL}
	};
	INITIALIZE_TABLE(genes, genesColumns);
	
	// STRAINS
	vector<DBColumn> strainsColumns = {
		{"strainId", DBType::INTEGER},
		{"geneIndex", DBType::INTEGER},
		{"geneId", DBType::INTEGER}
	};
	INITIALIZE_TABLE(strains, strainsColumns);
	
	// HOSTS
	vector<DBColumn> hostsColumns = {
		{"hostId", DBType::INTEGER},
		{"birthTime", DBType::REAL},
		{"deathTime", DBType::REAL}
	};
	INITIALIZE_TABLE(hosts, hostsColumns);
	
	// SAMPLED HOSTS
	vector<DBColumn> sampledHostsColumns = {
		{"time", DBType::REAL},
		{"hostId", DBType::INTEGER}
	};
	INITIALIZE_TABLE(sampledHosts, sampledHostsColumns);
	
	// SAMPLED TRANSMISSION EVENTS
	vector<DBColumn> sampledTransmissionsColumns = {
		{"time", DBType::REAL},
		{"transmissionId", DBType::INTEGER},
		{"sourceHostId", DBType::INTEGER},
		{"targetHostId", DBType::INTEGER}
	};
	INITIALIZE_TABLE(sampledTransmissions, sampledTransmissionsColumns);
	
	// SAMPLED HOSTS & TRANSMISSIONS: INFECTIONS
	vector<DBColumn> infectionColumns = {
		{"time", DBType::REAL},
		{"hostId", DBType::INTEGER},
		{"infectionId", DBType::INTEGER},
		{"strainId", DBType::INTEGER},
		{"geneIndex", DBType::INTEGER},
		{"active", DBType::INTEGER}
	};
	INITIALIZE_TABLE(sampledHostInfections, infectionColumns);
	INITIALIZE_TABLE(sampledTransmissionInfections, infectionColumns);
	
	// SAMPLED HOSTS AND TRANSMISSIONS: IMMUNITY
	vector<DBColumn> immunityColumns = {
		{"time", DBType::REAL},
		{"hostId", DBType::INTEGER},
		{"geneId", DBType::INTEGER},
		{"lossRate", DBType::REAL}
	};
	INITIALIZE_TABLE(sampledHostImmunity, immunityColumns);
	INITIALIZE_TABLE(sampledHostClinicalImmunity, immunityColumns);
	INITIALIZE_TABLE(sampledTransmissionImmunity, immunityColumns);
	INITIALIZE_TABLE(sampledTransmissionClinicalImmunity, immunityColumns);
}

void Simulation::run()
{
	clock_t startClock = clock();
	time_t startTime = time(nullptr);
	fprintf(stderr, "Starting at %s", ctime(&startTime));
	
	bool done = false;
	for(int64_t i = 0; !done; i++) {
		double tNextCommit = i * parPtr->dbCommitPeriod;
		
		dbPtr->beginTransaction();
		if(tNextCommit < parPtr->tEnd) {
			runUntil(tNextCommit);
		}
		else {
			runUntil(parPtr->tEnd);
			done = true;
		}
		dbPtr->commitWithRetry(DB_RETRY_DELAY, DB_TIMEOUT, cerr);
		cerr << "Committed at t = " << getTime() << '\n';
	}
	
	cout << "Total event count: " << queuePtr->getEventCount() << '\n';
	cout << "Transmission count: " << transmissionCount << '\n';
	
	time_t endTime = time(nullptr);
	clock_t endClock = clock();
	fprintf(stderr, "Ending at %s", ctime(&endTime));
	fprintf(stderr, "Total elapsed time: %f\n", elapsed(startClock, endClock));
	
	rusage resourceUsage;
	getrusage(RUSAGE_SELF, &resourceUsage);
	fprintf(stderr, "Memory usage: %ld\n", resourceUsage.ru_maxrss);
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

void Simulation::setEventRate(zppsim::RateEvent * event, double rate)
{
	event->setRate(*queuePtr, rate);
}

double Simulation::getSeasonality()
{
	return sin(2 * M_PI * getTime() / parPtr->tYear);
}



double Simulation::distanceWeightFunction(double d)
{
	assert(d > 0.0);
	return pow(d, -parPtr->distanceFunction.power);
}

Host * Simulation::drawDestinationHost(int64_t srcPopId)
{
	Population * srcPopPtr = popPtrs[srcPopId].get();
	
	vector<double> weights;
	for(auto & popPtr : popPtrs) {
		double dist = srcPopPtr->getDistance(popPtr.get());
		weights.push_back(
			distanceWeightFunction(dist)
			* popPtr->getBitingRate()
			* popPtr->size()
		);
	}
	int64_t dstPopId = sampleDiscreteLinearSearch(rng, weights);
	Population * dstPopPtr = popPtrs[dstPopId].get();
	int64_t dstHostIndex = drawUniformIndex(rng, dstPopPtr->size());
	return dstPopPtr->getHostAtIndex(dstHostIndex);
}

StrainPtr Simulation::generateRandomStrain()
{
	int64_t genesPerStrain = parPtr->genesPerStrain;
	
	// Uniformly randomly draw genes from pool
	std::vector<GenePtr> strainGenes(genesPerStrain);
	for(int64_t i = 0; i < genesPerStrain; i++) {
		strainGenes[i] = drawRandomGene();
	}
	
	return getStrain(strainGenes);
}

StrainPtr Simulation::recombineStrains(StrainPtr const & s1, StrainPtr const & s2)
{
	assert(s1->size() == s2->size());
	
	// Draw random subset of two strains
	vector<GenePtr> allGenes;
	allGenes.reserve(s1->size() + s2->size());
	copy(s1->genes.begin(), s1->genes.end(), std::back_inserter(allGenes));
	copy(s2->genes.begin(), s2->genes.end(), std::back_inserter(allGenes));
	assert(allGenes.size() == s1->size() + s2->size());
	
	vector<size_t> daughterIndices = drawUniformIndices(rng, allGenes.size(), size_t(s1->size()), false);
	vector<GenePtr> daughterGenes(daughterIndices.size());
	for(size_t i = 0; i < daughterIndices.size(); i++) {
		daughterGenes[i] = allGenes[daughterIndices[i]];
	}
	
	return getStrain(daughterGenes);
}

StrainPtr Simulation::mutateStrain(StrainPtr & strain)
{
	vector<int64_t> indices = drawMultipleBernoulli(rng, strain->size(), parPtr->pMutation);
	if(indices.size() == 0) {
		return strain;
	}
	else {
		vector<GenePtr> genes = strain->getGenes();
		for(int64_t index : indices) {
//			cerr << "start gene: " << genes[index]->id << '\n';
			genes[index] = mutateGene(genes[index]);
//			cerr << "end gene: " << genes[index]->id << '\n';
		}
		return getStrain(genes);
	}
}

void Simulation::updateRates()
{
//	cerr << getTime() << ": updating rates" << '\n';
	for(auto & popPtr : popPtrs) {
		popPtr->updateRates();
	}
}

void Simulation::sampleHosts()
{
	cerr << getTime() << ": sampling hosts" << '\n';
	for(auto & popPtr : popPtrs) {
		popPtr->sampleHosts();
	}
}

void Simulation::countTransmission()
{
	transmissionCount++;
}

GenePtr Simulation::drawRandomGene()
{
	int64_t geneIndex = drawUniformIndex(rng, genes.size());
	return genes[geneIndex];
}

GenePtr Simulation::mutateGene(GenePtr const & srcGenePtr) {
	int64_t srcGeneId = srcGenePtr->id;
	
	if(mutationDistributions.size() == 0) {
//		cerr << "not using mutation distributions" << '\n';
		return drawRandomGene();
	}
	else {
		assert(mutationDistributions.size() == genes.size());
		assert(mutationDistributions[srcGeneId].min() == 0);
		assert(mutationDistributions[srcGeneId].max() == genes.size() - 1);
//		cerr << "using mutation distributions" << '\n';
		return genes[mutationDistributions[srcGeneId](rng)];
	}
}


StrainPtr Simulation::getStrain(std::vector<GenePtr> const & strainGenes)
{
	StrainPtr strainPtr;
	auto strainItr = geneVecToStrainIndexMap.find(strainGenes);
	if(strainItr == geneVecToStrainIndexMap.end()) {
		strains.emplace_back(new Strain(nextStrainId++, strainGenes));
		strainPtr = strains.back();
		geneVecToStrainIndexMap[strainGenes] = strains.size() - 1;
		
		if(strainsTablePtr != nullptr) {
			DBRow row;
			row.setInteger("strainId", int64_t(strainPtr->id));
			for(int64_t i = 0; i < strainPtr->size(); i++) {
				row.setInteger("geneIndex", int64_t(i));
				row.setInteger("geneId", int64_t(strainPtr->getGene(i)->id));
				strainsTablePtr->insert(row);
			}
		}
	}
	else {
		strainPtr = strains[strainItr->second];
	}
	return strainPtr;
}

RateUpdateEvent::RateUpdateEvent(Simulation * simPtr, double initialTime, double period) :
	PeriodicEvent(initialTime, period), simPtr(simPtr)
{
}

void RateUpdateEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->updateRates();
}

HostStateSamplingEvent::HostStateSamplingEvent(Simulation * simPtr, double initialTime, double period) :
	PeriodicEvent(initialTime, period), simPtr(simPtr)
{
}

void HostStateSamplingEvent::performEvent(zppsim::EventQueue & queue)
{
	simPtr->sampleHosts();
}
