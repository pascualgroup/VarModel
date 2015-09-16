#include "Simulation.h"
#include "zppdb.hpp"
#include "zppjson.hpp"
#include "zppsim_util.hpp"
#include <sys/time.h>
#include <sys/resource.h>

// 100-millisecond delay between database commit retries
#define DB_RETRY_DELAY 100000

// 10-second maximum database retry before failure
#define DB_TIMEOUT 10000000

using namespace std;
using namespace zppjson;
using namespace zppdb;
using namespace zppsim;

static float elapsed(clock_t clockStart, clock_t clockEnd)
{
	return float(clockEnd - clockStart) / CLOCKS_PER_SEC;
}

static double getEntry(Array<Double> & vals, size_t index, size_t size)
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
	rng(int64_t(parPtr->randomSeed)),
	hostLifetimeDist(
		parPtr->hostLifetimeDistribution.pdf.toDoubleVector(),
		parPtr->hostLifetimeDistribution.x0,
		parPtr->hostLifetimeDistribution.dx.toDoubleVector()
	),
	queuePtr(new EventQueue(rng)),
	rateUpdateEvent(this, 0.0, parPtr->seasonalUpdateEvery),
	hostStateSamplingEvent(this, 0.0, parPtr->sampleHostsEvery),
	nextHostId(0),
	nextStrainId(0),
	transmissionCount(0),
    mutationCount(0),
	genesTable("genes"),
    lociTable("loci"),
	strainsTable("strains"),
	hostsTable("hosts"),
	sampledHostsTable("sampledHosts"),
	sampledHostInfectionTable("sampledHostInfections"),
	sampledHostImmunityTable("sampledHostImmunity"),
	sampledHostClinicalImmunityTable("sampledHostClinicalImmunity"),
	sampledTransmissionTable("sampledTransmissions"),
	sampledTransmissionStrainTable("sampledTransmissionStrains"),
	sampledTransmissionInfectionTable("sampledTransmissionInfections"),
	sampledTransmissionImmunityTable("sampledTransmissionImmunity"),
	sampledTransmissionClinicalImmunityTable("sampledTransmissionClinicalImmunity")
{
	// Construct transition probability distributions for genes
    /** turned off by hqx, new way of mutation, transition prob dist no use
	if(parPtr->genes.mutationWeights.size() > 1) {
		assert(parPtr->genes.mutationWeights.size() == parPtr->genePoolSize);
		for(int64_t i = 0; i < parPtr->genePoolSize; i++) {
			assert(parPtr->genes.mutationWeights[i].size() == parPtr->genePoolSize);
			vector<double> mwi = parPtr->genes.mutationWeights[i].toDoubleVector();
			mwi[i] = 0.0;
			mutationDistributions.emplace_back(mwi.begin(), mwi.end());
		}
	}*/
	
   
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
			parPtr->outputGenes,
            *dbPtr,
			genesTable
		));
	}


	
    // Create loci for each gene in the gene pool, hqx
    lociVec.reserve(parPtr->genePoolSize);
    Array<Double> vals = parPtr->genes.alleleNumber;
	if(vals.size() == 1) {
        for(size_t i =0; i<locusNumber; ++i) {
            alleleNumber.push_back((int64_t)(vals[0]));
        }
	}
	else {
		assert(vals.size() == locusNumber);
        for(size_t i =0; i<locusNumber; ++i) {
            alleleNumber.push_back((int64_t)(vals[i]));
        }
	}
    assert(alleleNumber.size()==size_t(locusNumber));
    std::vector<int64_t> Alleles(locusNumber);
    for(int64_t j = 0; j < locusNumber; j++) {
        std::uniform_int_distribution<int64_t> unif(0, alleleNumber[j]-1);
        Alleles[j] = unif(rng);
        //cout<<Alleles[j]<<" ";
    }
    //cout<<"\n";
    lociVec.emplace_back(new Loci(int64_t(0),
                                  Alleles,true,
                                  parPtr->outputLoci,
                                  *dbPtr,
                                  lociTable));
    int64_t i =0;
    while(i<parPtr->genePoolSize-1) {

        std::vector<int64_t> Alleles(locusNumber);
        for(int64_t j = 0; j < locusNumber; j++) {
            std::uniform_int_distribution<int64_t> unif(0, alleleNumber[j]-1);
            Alleles[j] = unif(rng);
            //cout<<Alleles[j]<<" ";
        }
        //cout<<"\n";
        int64_t checkId = recLociId(Alleles);
        if (checkId > i) {
            i++;
            lociVec.emplace_back(new Loci(i,
                                          Alleles,true,
                                          parPtr->outputLoci,
                                          *dbPtr,
                                          lociTable));
            
        }
    }

    
                
	// Create populations
	popPtrs.reserve(parPtr->populations.size());
	for(int64_t popId = 0; popId < parPtr->populations.size(); popId++) {
		popPtrs.emplace_back(new Population(this, popId));
	}
	
	dbPtr->commitWithRetry(DB_RETRY_DELAY, DB_TIMEOUT, cerr);
	
	cerr << "# events: " << queuePtr->size() << '\n';
}

GenePtr Simulation::createGene()
{
	assert(parPtr->genes.transmissibility.size() == 1);
	double transmissibility = parPtr->genes.transmissibility[0];
	
	assert(parPtr->genes.immunityLossRate.size() == 1);
	double immunityLossRate = parPtr->genes.immunityLossRate[0];
	
	assert(parPtr->genes.clinicalImmunityLossRate.size() == 1);
	double clinicalImmunityLossRate = parPtr->genes.clinicalImmunityLossRate[0];
	
	int64_t index = genes.size();
	genes.emplace_back(new Gene(
		index,
		transmissibility,
		immunityLossRate,
		clinicalImmunityLossRate,
		dbPtr->tableExists(genesTable),
		*dbPtr,
		genesTable
	));
	return genes.back();
}

void Simulation::initializeDatabaseTables()
{
	dbPtr->createTable(genesTable);
	dbPtr->createTable(lociTable);
	dbPtr->createTable(strainsTable);
	dbPtr->createTable(hostsTable);
	dbPtr->createTable(sampledHostsTable);
	dbPtr->createTable(sampledHostInfectionTable);
	dbPtr->createTable(sampledHostImmunityTable);
	dbPtr->createTable(sampledHostClinicalImmunityTable);
	dbPtr->createTable(sampledTransmissionTable);
	dbPtr->createTable(sampledTransmissionStrainTable);
	dbPtr->createTable(sampledTransmissionInfectionTable);
	dbPtr->createTable(sampledTransmissionImmunityTable);
	dbPtr->createTable(sampledTransmissionClinicalImmunityTable);
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
	return hostLifetimeDist.draw(rng);
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



double Simulation::distanceWeightFunction(double d)
{
	assert(d > 0.0);
	return pow(d, -parPtr->distancePower);
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

StrainPtr Simulation::generateRandomStrain(int64_t nNewGenes)
{
	cerr << "generating random strain with new genes " << endl;
	
	int64_t genesPerStrain = parPtr->genesPerStrain;
	
	// Generate shuffled list of the locations of the new genes within the strain
	vector<bool> newGeneLocations(genesPerStrain);
	for(int64_t i = 0; i < nNewGenes; i++) {
		newGeneLocations[i] = true;
	}
	for(int64_t i = nNewGenes; i < genesPerStrain; i++) {
		newGeneLocations[i] = false;
	}
	shuffle(newGeneLocations.begin(), newGeneLocations.end(), rng);
	
	// Assemble strain: old genes at false entries; new genes at true entries
	std::vector<GenePtr> strainGenes(genesPerStrain);
	int64_t nNewGenesCheck = 0;
	for(int64_t i = 0; i < genesPerStrain; i++) {
		if(newGeneLocations[i]) {
			nNewGenesCheck++;
			strainGenes[i] = createGene();
		}
		else {
			strainGenes[i] = drawRandomGene();
		}
	}
	assert(nNewGenes == nNewGenesCheck);
	
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

//new mutation mode -> change function "mutateGene"
//select one gene from the strain to mutate
StrainPtr Simulation::mutateStrain(StrainPtr & strain)
{
	int64_t index = drawUniformIndex(rng,strain->size());
    vector<GenePtr> genes = strain->getGenes();
    genes[index] = mutateGene(genes[index]);
    mutationCount++;
    return getStrain(genes);
}

// randomly select two genes in the strain, and recombine
StrainPtr Simulation::ectopicRecStrain(StrainPtr & strain)
{
    vector<int64_t> indices = drawUniformIndices(rng,strain->size(),int64_t(2),false);
    std::vector<GenePtr> curStrainGenes = strain->getGenes();
    bernoulli_distribution flipCoin(parPtr->percConversion);
    bool isConversion = flipCoin(rng);
    //if(!isConversion) {
    //    cout<<"not conversion"<<endl;
    //}
    vector<GenePtr> genesPtrAfterEctopicRecomb = ectopicRecomb(curStrainGenes[indices[0]], curStrainGenes[indices[1]],isConversion);
    curStrainGenes[indices[0]] = genesPtrAfterEctopicRecomb[0];
    curStrainGenes[indices[1]] = genesPtrAfterEctopicRecomb[1];
    return getStrain(curStrainGenes);
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

void Simulation::recordTransmission(Host &srcHost, Host &dstHost, std::vector<StrainPtr> &strains)
{
	if((transmissionCount + 1) % parPtr->sampleTransmissionEventEvery == 0) {
		TransmissionRow row;
		row.time = getTime();
		row.transmissionId = transmissionCount;
		row.sourceHostId = srcHost.id;
		row.targetHostId = dstHost.id;
		dbPtr->insert(sampledTransmissionTable, row);
		
		for(auto & strainPtr : strains) {
			TransmissionStrainRow row;
			row.transmissionId = transmissionCount;
			row.strainId = strainPtr->id;
			dbPtr->insert(sampledTransmissionStrainTable, row);
		}
		
		srcHost.writeInfections(transmissionCount, *dbPtr, sampledTransmissionInfectionTable);
		dstHost.writeInfections(transmissionCount, *dbPtr, sampledTransmissionInfectionTable);
		
		srcHost.immunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		srcHost.clinicalImmunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		dstHost.immunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		dstHost.clinicalImmunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
	}
	
	transmissionCount++;
}

GenePtr Simulation::drawRandomGene()
{
    bool func = false;
    int64_t geneIndex;
    while(!func) {
        geneIndex = drawUniformIndex(rng, genes.size());
        func = lociVec[geneIndex]->functionality;
    }
	return genes[geneIndex];
}

GenePtr Simulation::drawRandomGeneExcept(int64_t geneId)
{
	int64_t newGeneId = drawUniformIndexExcept(rng, (int64_t)genes.size(), geneId);
	return genes[newGeneId];
}

GenePtr Simulation::mutateGene(GenePtr const & srcGenePtr) {
    //cout<<"mutate gene\n";
    int64_t srcLociId = srcGenePtr->id;
    std::vector<int64_t> srcLociAlleles = lociVec[srcLociId]->Alleles;
    //not using mutationDistributions anymore, mutation weights are locus specific weight, sum(weight) = 1
    //to do: redefine mutationDistributions
    size_t mutateLocusId;
    if(mutationDistributions.size() == 0) {
        mutateLocusId = drawUniformIndex(rng, srcLociAlleles.size());
    }else{
        assert(mutationDistributions.size()==srcLociAlleles.size());
        mutateLocusId = sampleDiscreteLinearSearch(rng, mutationDistributions);
        //cout<<mutateLocusId<<endl;
    }
    //mutate allele
    alleleNumber[mutateLocusId]++;
    std::vector<int64_t> newLoci = srcLociAlleles;
    newLoci[mutateLocusId] = alleleNumber[mutateLocusId]-1;
    //cout<<newLoci[mutateLocusId]<<endl;
    int64_t newGeneId = genes.size();
    lociVec.emplace_back(new Loci(newGeneId,
                                  newLoci,true,parPtr->outputLoci,*dbPtr,lociTable));
    return createGene();
}

//test whether a new allele vector already exist in the lociVec
//return the id number of the vector, all the new gene id
int64_t Simulation::recLociId(std::vector<int64_t> & recGeneAlleles) {
    for (LociPtr j : lociVec) {
        if (recGeneAlleles == j->Alleles) {
            return j->id;
        }
    }
    return lociVec.size();
}

double Simulation::parentsSimilarity(GenePtr const & pGene1, GenePtr const & pGene2) {
    double pSim = 0;
    std::vector<int64_t> pGene1Alleles = lociVec[pGene1->id]->Alleles;
    std::vector<int64_t> pGene2Alleles = lociVec[pGene2->id]->Alleles;
    for (int64_t i=0; i<locusNumber; ++i) {
        if (pGene1Alleles[i] == pGene2Alleles[i]) {
            pSim += 1./locusNumber;
        }
    }
    //cout<<pSim<<endl;
    return pSim>0 ? pSim:0.01;
}

//recombine or conversion of two parent genes, return gene pointers of the two new genes
std::vector<GenePtr> Simulation::ectopicRecomb(GenePtr const & pGene1, GenePtr const & pGene2, bool isConversion) {
    //cout<<"ectopic recombine gene\n";
    std::vector<GenePtr> returnGenes;
    returnGenes.push_back(pGene1);
    returnGenes.push_back(pGene2);
    int64_t breakPoint = 1 + drawUniformIndex(rng,locusNumber-2);
    std::vector<int64_t> pGene1Alleles = lociVec[pGene1->id]->Alleles;
    std::vector<int64_t> pGene2Alleles = lociVec[pGene2->id]->Alleles;
    bernoulli_distribution flipCoin(parentsSimilarity(pGene1,pGene2));
    bool recFunction[] = {flipCoin(rng),flipCoin(rng)};//whether functional for the two recombinants;
    std::vector<int64_t> recGene1Alleles(locusNumber);
    std::vector<int64_t> recGene2Alleles(locusNumber);
    for (int64_t i=0; i<locusNumber; ++i) {
        if (i<breakPoint) {
            recGene1Alleles[i] = pGene1Alleles[i];
            recGene2Alleles[i] = pGene2Alleles[i];
        } else {
            recGene1Alleles[i] = pGene2Alleles[i];
            recGene2Alleles[i] = pGene1Alleles[i];
        }
    }
    //test whether the new recombinant is already in lociVec matrix
    int64_t recId[] = {recLociId(recGene1Alleles),recLociId(recGene2Alleles)};
    if (recId[0] == lociVec.size()) {
    //get whether the new recombinant is functional
        lociVec.emplace_back(new Loci(recId[0],
                                      recGene1Alleles,recFunction[0],parPtr->outputLoci,*dbPtr,lociTable));
        GenePtr recGenePtr = createGene();
        recId[1]++;
        if (recFunction[0]) {
           returnGenes[0] = recGenePtr;
        }
    } else {
        if (lociVec[recId[0]]->functionality) {
            returnGenes[0] = genes[recId[0]];
        }
    }
    if (isConversion) {
        return returnGenes;
        //cout<<returnGenes[0]->id<<endl;
        //cout<<returnGenes[1]->id<<endl;
    }else{
        if(recId[1] == lociVec.size()) {
            lociVec.emplace_back(new Loci(recId[1],
                                          recGene2Alleles,recFunction[1],parPtr->outputLoci,*dbPtr,lociTable));
            GenePtr recGenePtr = createGene();
            if (recFunction[1]) {
                returnGenes[1] = recGenePtr;
            }
        } else {
            if (lociVec[recId[1]]->functionality) {
                returnGenes[1] = genes[recId[1]];
            }
        }
        //cout<<returnGenes[0]->id<<endl;
        //cout<<returnGenes[1]->id<<endl;
        
        return returnGenes;
    }
}
/**
GenePtr Simulation::mutateGene(GenePtr const & srcGenePtr) {
	int64_t srcGeneId = srcGenePtr->id;
	
	if(mutationDistributions.size() == 0) {
//		cerr << "not using mutation distributions" << '\n';
		return drawRandomGeneExcept(srcGeneId);
	}
	else {
		assert(mutationDistributions.size() == genes.size());
		assert(mutationDistributions[srcGeneId].min() == 0);
		assert(mutationDistributions[srcGeneId].max() == genes.size() - 1);
//		cerr << "using mutation distributions" << '\n';
		return genes[mutationDistributions[srcGeneId](rng)];
	}
}
*/

StrainPtr Simulation::getStrain(std::vector<GenePtr> const & strainGenes)
{
	StrainPtr strainPtr;
	auto strainItr = geneVecToStrainIndexMap.find(strainGenes);
	if(strainItr == geneVecToStrainIndexMap.end()) {
		strains.emplace_back(new Strain(nextStrainId++, strainGenes));
		strainPtr = strains.back();
		geneVecToStrainIndexMap[strainGenes] = strains.size() - 1;
		
		if(parPtr->outputStrains) {
			StrainRow row;
			row.strainId = strainPtr->id;
			for(int64_t i = 0; i < strainPtr->size(); i++) {
				row.geneIndex = i;
				row.geneId = strainPtr->getGene(i)->id;
				dbPtr->insert(strainsTable, row);
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
