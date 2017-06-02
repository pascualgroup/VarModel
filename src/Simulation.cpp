#include "Simulation.h"
#include "zppdb.hpp"
#include "zppjson.hpp"
#include "zppsim_util.hpp"
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#include <algorithm>
#include <string>
#include <fstream>
#include <iterator>

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
	hostStateSamplingEvent(this, parPtr->burnIn, parPtr->sampleHostsEvery),
    mdaEvent(this,parPtr->MDA.TimeStartMDA, parPtr->MDA.interval),
    irsEvent(this,parPtr->intervention.TimeStart),
    removeirsEvent(this,parPtr->intervention.TimeStart+parPtr->intervention.duration),
	nextHostId(0),
	nextStrainId(0),
	transmissionCount(0),
    mutationCount(0),
	genesTable("genes"),
    lociTable("loci"),
	strainsTable("strains"),
	hostsTable("hosts"),
    microsatTable("microsats"),
    alleleImmunityTable("hostsAlleleImmunityHistory"),
    InfectionDurationTable("InfectionDuration"),
    recordEIRTable("recordEIR"),
    sampledHostsTable("sampledHosts"),
	sampledHostInfectionTable("sampledHostInfections"),
	sampledHostImmunityTable("sampledHostImmunity"),
	sampledHostClinicalImmunityTable("sampledHostClinicalImmunity"),
	sampledTransmissionTable("sampledTransmissions"),
	sampledTransmissionStrainTable("sampledTransmissionStrains"),
	sampledTransmissionInfectionTable("sampledTransmissionInfections"),
	sampledTransmissionImmunityTable("sampledTransmissionImmunity"),
	sampledTransmissionClinicalImmunityTable("sampledTransmissionClinicalImmunity"),
    followedHostsTable("followedHosts")
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
	
    //make sure burnin time is smaller than end time
    assert(parPtr->burnIn<parPtr->tEnd);
    
	queuePtr->addEvent(&rateUpdateEvent);
	queuePtr->addEvent(&hostStateSamplingEvent);
    if (parPtr->MDA.includeMDA) queuePtr->addEvent(&mdaEvent);
    if (parPtr->intervention.includeIntervention) {
        queuePtr->addEvent(&irsEvent);
        queuePtr->addEvent(&removeirsEvent);
    };
    
    //set migrationPoolSize equal to the initial gene pool size
    migrationPoolSize = parPtr->genePoolSize;
    
    //create variant size for each locus
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
    
	// Create gene pool
	genes.reserve(parPtr->genePoolSize);
	for(int64_t i = 0; i < parPtr->genePoolSize; i++) {
		double transmissibility = getEntry(
			parPtr->genes.transmissibility, i, parPtr->genePoolSize
		);
		double immunityLossRate = getEntry(
			parPtr->genes.immunityLossRate, i, parPtr->genePoolSize
		);
        std::vector<int64_t> Alleles(locusNumber);
        if (i== 0) {
            for(int64_t j = 0; j < locusNumber; j++) {
                std::uniform_int_distribution<int64_t> unif(0, alleleNumber[j]-1);
                Alleles[j] = unif(rng);
                //cout<<Alleles[j]<<" ";
            }
        }else{
            int64_t checkId = 0;
            while(checkId < i) {
                for(int64_t j = 0; j < locusNumber; j++) {
                    std::uniform_int_distribution<int64_t> unif(0, alleleNumber[j]-1);
                    Alleles[j] = unif(rng);
                    //cout<<Alleles[j]<<" ";
                }
                //cout<<"\n";
                checkId = recLociId(Alleles,genes);
            }
        }
        

		genes.emplace_back(new Gene(
			i,
			transmissibility,
			immunityLossRate,
			0,
            true,
            Alleles,
			parPtr->outputGenes,
            parPtr->outputLoci,
            *dbPtr,
			genesTable,
            lociTable
		));
	}
    
    
    //create allele size range for microsatellites, if required
    if (microsatNumber>0){
        Array<Double> msvals = parPtr->genes.microsatAlleles;
        if(msvals.size() == 1) {
            for(size_t i =0; i<microsatNumber; ++i) {
                microsatAlleles.push_back((int64_t)(msvals[0]));
            }
        }
        else {
            assert(msvals.size() == microsatNumber);
            for(size_t i =0; i<microsatNumber; ++i) {
                microsatAlleles.push_back((int64_t)(msvals[i]));
            }
        }
        assert(microsatAlleles.size()==size_t(microsatNumber));
    }
                
	// Create populations
	popPtrs.reserve(parPtr->populations.size());
	for(int64_t popId = 0; popId < parPtr->populations.size(); popId++) {
		popPtrs.emplace_back(new Population(this, popId));
	}
	
	dbPtr->commitWithRetry(DB_RETRY_DELAY, DB_TIMEOUT, cerr);
	
	cerr << "# events: " << queuePtr->size() << '\n';
}

//5/17 new edits, microsats created from real distributions
// create microsats pool

std::vector<std::vector<int64_t>> Simulation::readMsArray(int year) {
    int ttyear = (int)floor(parPtr->tEnd/360.0);
    char temp[512];
    sprintf(temp, "python readinMS.py %d %d", year, ttyear);
    system((char *)temp);
    std::string msFilename("tempMS.txt");
    std::ifstream in(msFilename.c_str());
	if(!in) {
		cerr << "ms file " << msFilename << " cannot be opened.";
	}
    std::vector<std::vector<int64_t> > v;
	
    std::string line;
    while ( getline( in, line ) ) {
        std::istringstream is( line );
        v.push_back(std::vector<int64_t>( std::istream_iterator<int64_t>(is),
                                     std::istream_iterator<int64_t>() ) );
    }
    in.close();
    return (v);
}


//function to run simcoal to generate microsatellites
void Simulation::runMSSimCoal(size_t msSampleSize) {
    char temp[512];
    sprintf(temp, "python generateMS.py %zu %.1f %d %.10f", msSampleSize, (double)parPtr->tEnd, (int)parPtr->genes.microsatNumber, (double)parPtr->pMsMutate);
    system((char *)temp);
}

GenePtr Simulation::createGene(std::vector<int64_t> Alleles,bool const functionality,int64_t const source)
{
	assert(parPtr->genes.transmissibility.size() == 1);
	double transmissibility = parPtr->genes.transmissibility[0];
	
	assert(parPtr->genes.immunityLossRate.size() == 1);
	double immunityLossRate = parPtr->genes.immunityLossRate[0];
	
	int64_t index = genes.size();
    genes.emplace_back(new Gene(
                                index,
                                transmissibility,
                                immunityLossRate,
                                source,
                                functionality,
                                Alleles,
                                parPtr->outputGenes,
                                parPtr->outputLoci,
                                *dbPtr,
                                genesTable,
                                lociTable
                                ));
	return genes.back();
}

void Simulation::initializeDatabaseTables()
{
	dbPtr->createTable(genesTable);
	dbPtr->createTable(lociTable);
	dbPtr->createTable(strainsTable);
	dbPtr->createTable(hostsTable);
	dbPtr->createTable(alleleImmunityTable);
	dbPtr->createTable(sampledHostsTable);
	dbPtr->createTable(sampledHostInfectionTable);
	dbPtr->createTable(sampledHostImmunityTable);
	dbPtr->createTable(sampledHostClinicalImmunityTable);
    dbPtr->createTable(InfectionDurationTable);
    dbPtr->createTable(recordEIRTable);
    dbPtr->createTable(microsatTable);
	dbPtr->createTable(sampledTransmissionTable);
	dbPtr->createTable(sampledTransmissionStrainTable);
	dbPtr->createTable(sampledTransmissionInfectionTable);
	dbPtr->createTable(sampledTransmissionImmunityTable);
	dbPtr->createTable(sampledTransmissionClinicalImmunityTable);
    dbPtr->createTable(followedHostsTable);
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
	//cerr << "generating random strain with new genes " << endl;
	
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
			strainGenes[i] = mutateGene(drawRandomGene(),3);
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

GenePtr Simulation::recombineMS(GenePtr const & ms1, GenePtr const & ms2)
{
    vector<int64_t> newMS = ms1->Alleles;
    for (size_t i=0; i<ms1->Alleles.size(); ++i) {
        if (drawUniformIndex(rng,2) == 1) {
            newMS[i] = ms2->Alleles[i];
            //cout<<"ms2Allele is "<<ms2->Alleles[i]<<endl;
        }
    }
    int64_t index = recLociId(newMS,microsats);
    if(index == microsats.size()) {
        return createMicrosat(newMS);
    }else{
        return microsats[index];
    }
}


//new mutation mode -> change function "mutateGene"
//select one gene from the strain to mutate
std::vector<GenePtr> Simulation::mutateStrain(StrainPtr & strain)
{
	int64_t index = drawUniformIndex(rng,strain->size());
    vector<GenePtr> genes = strain->getGenes();
    genes[index] = mutateGene(genes[index],2);
    mutationCount++;
    return genes;
}

// randomly select two genes in the strain, and recombine
std::vector<GenePtr> Simulation::ectopicRecStrain(StrainPtr & strain)
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
    return curStrainGenes;
}

//generate random microsatellite alleles
GenePtr Simulation::generateRandomMicrosat()
{
    std::vector<int64_t> Alleles(microsatNumber);
    for(int64_t j = 0; j < microsatNumber; j++) {
        std::uniform_int_distribution<int64_t> unif(0, microsatAlleles[j]-1);
        Alleles[j] = unif(rng);
        //cout<<Alleles[j]<<" ";
    }
    //cout<<"\n";
    int64_t checkId = recLociId(Alleles,microsats);
    if (checkId == microsats.size()) {
        return createMicrosat(Alleles);
    }else{
        return microsats[checkId];
    }
}

//check if the microsatellite haplotype exist, if not store in microSat vector
GenePtr Simulation::storeMicrosat(std::vector<int64_t> Alleles){
    int64_t checkId = recLociId(Alleles,microsats);
    if (checkId == microsats.size()) {
        return createMicrosat(Alleles);
    }else{
        return microsats[checkId];
    }
};

//create new microsatellite alleles and store in microsats
GenePtr Simulation::createMicrosat(std::vector<int64_t> Alleles)
{
    int64_t index = microsats.size();
    microsats.emplace_back(new Gene(
                                    index,
                                    0,
                                    0,
                                    0,
                                    true,
                                    Alleles,
                                    false,
                                    parPtr->outputLoci,
                                    *dbPtr,
                                    genesTable,
                                    microsatTable
                                    ));
    return microsats.back();
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
    double t = getTime();
    //if (t > parPtr->burnIn) {
	   cerr << t << ": sampling hosts" << '\n';
        for(auto & popPtr : popPtrs) {
            popPtr->sampleHosts();
        }
    //}
}

void Simulation::MDA()
{
    double t = getTime();
    //change the global pool size for immigration
    migrationPoolSize = parPtr->genePoolSize * parPtr->MDA.MDAGenePoolChange;
    if (mdaCounts > parPtr->MDA.totalNumber) {
        cout<<"remove MDAs at "<<t<<endl;
        removeEvent(&mdaEvent);
        //restore the original migration rate
        //shai preferred the migration and global pool to be kept low
        /*for(auto & popPtr : popPtrs) {
            popPtr->setEventRate(popPtr->immigrationEvent.get(),popPtr->getImmigrationRate());
        }*/
    }else{
        cerr << t << ": start MDA" << '\n';
        for(auto & popPtr : popPtrs) {
            popPtr->executeMDA(t+parPtr->MDA.drugEffDuration);
        }
        mdaCounts++;
        cerr << "totalMDA " <<mdaCounts<< '\n';
    }
}

void Simulation::IRS()
{
    //set biting rate amplitude to IRS biting rate amplitude
    //reduce immigration rate
    migrationPoolSize = parPtr->genePoolSize * parPtr->intervention.IRSGenePoolChange;
    for(auto & popPtr : popPtrs) {
        popPtr->IRSBitingAmplitude = parPtr->intervention.amplitude;
        popPtr->setEventRate(popPtr->immigrationEvent.get(),popPtr->getImmigrationRate()*parPtr->intervention.IRSMRateAmplitude);
    }
    
}

void Simulation::RemoveIRS()
{
    //set biting rate amplitude to back to 1
    //not changing migration rate, nor the global pool size
    for(auto & popPtr : popPtrs) {
        popPtr->IRSBitingAmplitude = 1;
        //popPtr->setEventRate(popPtr->immigrationEvent.get(),popPtr->getImmigrationRate());
    }
}

void Simulation::recordImmunity(Host & host, int64_t locusIndex, int64_t alleleId) {
    AlleleImmunityRow row;
    row.time = getTime();
    row.hostId = host.id;
    row.locusIndex = locusIndex;
    row.alleleId = alleleId;
    dbPtr->insert(alleleImmunityTable,row);
}

void Simulation::recordTransmission(Host &srcHost, Host &dstHost, std::vector<StrainPtr> &strains)
{

	if((transmissionCount + 1) % parPtr->sampleTransmissionEventEvery == 0) {
        for(auto & strainPtr : strains) {
            TransmissionStrainRow row;
            row.time = getTime();
            row.transmissionId = transmissionCount;
            row.strainId = strainPtr->id;
            dbPtr->insert(sampledTransmissionStrainTable, row);
            strainPtr->writeToDatabaseStrain(*dbPtr, strainsTable, genesTable, lociTable);
        }
		TransmissionRow row;
		row.transmissionId = transmissionCount;
		row.sourceHostId = srcHost.id;
		row.targetHostId = dstHost.id;
		dbPtr->insert(sampledTransmissionTable, row);
		
		/**
		srcHost.writeInfectionks(transmissionCount, *dbPtr, sampledTransmissionInfectionTable);
		dstHost.writeInfections(transmissionCount, *dbPtr, sampledTransmissionInfectionTable);
		
		srcHost.immunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		srcHost.clinicalImmunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		dstHost.immunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
		dstHost.clinicalImmunity.write(transmissionCount, *dbPtr, sampledTransmissionImmunityTable);
        */
	}
	
	transmissionCount++;
}

void Simulation::writeDuration(std::list<Infection>::iterator infectionItr)
{
    bernoulli_distribution flipCoin(0.001);
    if(flipCoin(rng)) {
        InfectionDurationRow row;
        row.time = infectionItr->initialTime;
        row.duration = getTime() - infectionItr->initialTime;
        row.hostId = infectionItr->hostPtr->id;
        row.infectionId = infectionItr->id;
        dbPtr->insert(InfectionDurationTable, row);
    }
}

void Simulation::writeFollowedHostInfection(std::list<Infection>::iterator infectionItr)
{
    if (parPtr->following.includeHostFollowing) {
        if (infectionItr->hostPtr->toTrack) {
            followedHostsRow row;
            row.time = infectionItr->initialTime;
            row.duration = getTime()-infectionItr->initialTime;
            row.hostId = infectionItr->hostPtr->id;
            row.infectionId = infectionItr->id;
            row.popId = infectionItr->hostPtr->popPtr->id;
            row.strainId = infectionItr->strainPtr->id;
            double t = getTime() - infectionItr->hostPtr->birthTime;
            
            //if time surpass the tracking period
            if (t > parPtr->following.followDuration) {
               // cout<<"it's larger"<<t<<endl;
                
                infectionItr->hostPtr->toTrack = false;
                
            }else{
                //cout<<"it's within duration"<<t<<endl;
                dbPtr->insert(followedHostsTable, row);
                infectionItr->strainPtr->writeToDatabaseStrain(*dbPtr, strainsTable, genesTable, lociTable);
            }
            

        }
    }
}

void Simulation::writeEIR(double time, int64_t infectious)
{
    bernoulli_distribution flipCoin(0.001);
    if(flipCoin(rng)) {
        recordEIRRow row;
        row.time = time;
        row.infectious = infectious;
        dbPtr->insert(recordEIRTable, row);
    }
    
}



//for immigration events, only sample from the large pool that exists
GenePtr Simulation::drawRandomGene()
{
    bool func = false;
    int64_t geneIndex;
    while(!func) {
        geneIndex = drawUniformIndex(rng, migrationPoolSize);
        func = genes[geneIndex]->functionality;
    }
	return genes[geneIndex];
}

GenePtr Simulation::drawRandomGeneExcept(int64_t geneId)
{
	int64_t newGeneId = drawUniformIndexExcept(rng, (int64_t)genes.size(), geneId);
	return genes[newGeneId];
}

GenePtr Simulation::mutateGene(GenePtr const & srcGenePtr, int64_t const source) {
    //cout<<"mutate gene\n";
    std::vector<int64_t> srcLociAlleles = srcGenePtr->Alleles;
    //not using mutationDistributions anymore, mutation weights are locus specific weight, sum(weight) = 1
    size_t mutateLocusId;
    if(mutationDistributions.size() == 1) {
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
    return createGene(newLoci,true,source);
}

GenePtr Simulation::mutateMS(GenePtr const & srcMS) {
    std::vector<int64_t> srcLociAlleles = srcMS->Alleles;
    size_t mutateLocusId;
    //stepwise mutation of microsatellite allele, by adding 1 or minus 1 length
    mutateLocusId = drawUniformIndex(rng, srcLociAlleles.size());
    std::vector<int64_t> newLoci = srcLociAlleles;
    bernoulli_distribution flipCoin(0.5);
    int incDecrease = (flipCoin(rng)-0.5)*2;//whether adding 1 or minus 1
    //cout<<"incerase "<<incDecrease<<endl;
    newLoci[mutateLocusId] += incDecrease;
    return storeMicrosat(newLoci);
}

//test whether a new allele vector already exist in the genes allele vectors
//return the id number of the vector, all the new gene id
int64_t Simulation::recLociId(std::vector<int64_t> & recGeneAlleles, std::vector<GenePtr> & searchSet) {
    for (GenePtr j : searchSet) {
        if (recGeneAlleles == j->Alleles) {
            return j->id;
        }
    }
    return searchSet.size();
}

double Simulation::parentsSimilarity(GenePtr const & pGene1, GenePtr const & pGene2, int64_t breakPoint) {
    double pDiv = 0;
    double childDiv = 0;
    double rho = 0.8; //recombination tolerance;
    double averageMutation = 5; //average number of mutations per epitope
    std::vector<int64_t> pGene1Alleles = pGene1->Alleles;
    std::vector<int64_t> pGene2Alleles = pGene2->Alleles;
    for (int64_t i=0; i<locusNumber; ++i) {
        if (pGene1Alleles[i] != pGene2Alleles[i]) {
            pDiv += 1;
            if (i<breakPoint) {
                childDiv += 1;
            }
        }
    }
    double rhoPower = childDiv * averageMutation * (pDiv-childDiv)* averageMutation /(pDiv * averageMutation -1);
    double survProb = pow(rho, rhoPower);
    //cout<<survProb<<endl;
    return survProb;
}

//recombine or conversion of two parent genes, return gene pointers of the two new genes
//11-20 change the recombination probability to incorporate child's similarity to the parents.
std::vector<GenePtr> Simulation::ectopicRecomb(GenePtr const & pGene1, GenePtr const & pGene2, bool isConversion) {
    //cout<<"ectopic recombine gene\n";
    std::vector<GenePtr> returnGenes;
    returnGenes.push_back(pGene1);
    returnGenes.push_back(pGene2);
    int64_t breakPoint = drawUniformIndex(rng,locusNumber);
    if ((breakPoint == 0) || (breakPoint == (locusNumber-1))) {
        if(isConversion) {
            if(breakPoint>0) {
                returnGenes[0] = pGene2;
            }else{
                returnGenes[1] = pGene1;
            }
        }
        return returnGenes;
    }else {
        std::vector<int64_t> pGene1Alleles = pGene1->Alleles;
        std::vector<int64_t> pGene2Alleles = pGene2->Alleles;
        bernoulli_distribution flipCoin(parentsSimilarity(pGene1,pGene2, breakPoint));
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
        //test whether the new recombinant is already in genes matrix
        int64_t recId = recLociId(recGene1Alleles,genes);
        if (recId == genes.size()) {
            //get whether the new recombinant is functional
            GenePtr recGenePtr = createGene(recGene1Alleles,recFunction[0],1);
            if (recFunction[0]) {
                returnGenes[0] = recGenePtr;
            }
        } else {
            if (genes[recId]->functionality) {
                returnGenes[0] = genes[recId];
            }
        }
        if (isConversion) {
            return returnGenes;
            //cout<<returnGenes[0]->id<<endl;
            //cout<<returnGenes[1]->id<<endl;
        }else{
            recId = recLociId(recGene2Alleles,genes);
            if(recId == genes.size()) {
                GenePtr recGenePtr = createGene(                                           recGene2Alleles,recFunction[1],1);
                if (recFunction[1]) {
                    returnGenes[1] = recGenePtr;
                }
            } else {
                if (genes[recId]->functionality) {
                returnGenes[1] = genes[recId];
                }
            }
            //cout<<returnGenes[0]->id<<endl;
            //cout<<returnGenes[1]->id<<endl;
        
        return returnGenes;
        }
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

StrainPtr Simulation::getStrain(std::vector<GenePtr> const & oriStrainGenes)
{
	StrainPtr strainPtr;
    std::vector<GenePtr> strainGenes = oriStrainGenes;
    std::sort (strainGenes.begin(),strainGenes.end());

	auto strainItr = geneVecToStrainIndexMap.find(strainGenes);
	if(strainItr == geneVecToStrainIndexMap.end()) {
		strains.emplace_back(new Strain(nextStrainId++, strainGenes,parPtr->outputStrains,*dbPtr, strainsTable));
		strainPtr = strains.back();
		geneVecToStrainIndexMap[strainGenes] = strains.size() - 1;
		
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

//add MDA events to simulate giving drugs to all hosts
MDAEvent::MDAEvent(Simulation * simPtr, double initialTime, double period):
    PeriodicEvent(initialTime,period),simPtr(simPtr)
{
}

void MDAEvent::performEvent(zppsim::EventQueue & queue)
{
    simPtr->MDA();
}

//add IRS event
IRSEvent::IRSEvent(Simulation * simPtr, double time):
    OneTimeEvent(time),simPtr(simPtr)
{
}

void IRSEvent::performEvent(zppsim::EventQueue & queue)
{
    simPtr->IRS();
}

//remove IRS event
RemoveIRSEvent::RemoveIRSEvent(Simulation * simPtr, double time):
OneTimeEvent(time),simPtr(simPtr)
{
}

void RemoveIRSEvent::performEvent(zppsim::EventQueue & queue)
{
    simPtr->RemoveIRS();
}
