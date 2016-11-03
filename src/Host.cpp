#include "Host.h"
#include "Simulation.h"

#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;
using namespace zppsim;

template <typename T>
std::unordered_map<size_t,size_t> ordered(std::vector<T> const& values) {
    std::vector<size_t> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<size_t>(0));
    
    std::sort(
              begin(indices), end(indices),
              [&](size_t a, size_t b) { return values[a] < values[b]; }
              );
    std::unordered_map<size_t,size_t> reorderIndexMap;
    for(size_t i = 0; i<values.size(); i++) {
        reorderIndexMap.insert({indices[i],i});
    }
    return reorderIndexMap;
}

Host::Host(
	Population * popPtr, int64_t id, double birthTime, double deathTime,
	bool writeToDatabase,
	Database & db,
	zppdb::Table<HostRow> & table
):
	id(id), popPtr(popPtr),
	birthTime(birthTime), deathTime(deathTime), nextInfectionId(0),
	deathEvent(new DeathEvent(this)),
	immunity(this, false,popPtr->simPtr->locusNumber,popPtr->simPtr->parPtr->withinHost.infectionTimesToImmune),
	clinicalImmunity(this, true,popPtr->simPtr->locusNumber,popPtr->simPtr->parPtr->withinHost.infectionTimesToImmune)
{
//	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
	
	if(writeToDatabase) {
		HostRow row;
		row.hostId = id;
        row.popId = popPtr->id;
		row.birthTime = birthTime;
		row.deathTime = deathTime;
		db.insert(table, row);
	}
	addEvent(deathEvent.get());
}

double Host::getAge()
{
	return getTime() - birthTime;
}

int64_t Host::getActiveInfectionCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			count++;
		}
	}
	return count;
}

std::vector<GenePtr> Host::getActiveInfectionGenes()
{
	vector<GenePtr> genes;
	genes.reserve(infections.size());
	for(auto & infection : infections) {
		if(infection.isActive()) {
			genes.push_back(infection.getCurrentGene());
		}
	}
	return genes;
}

std::vector<int64_t> Host::getActiveInfectionGeneIds()
{
	vector<int64_t> ids;
	ids.reserve(infections.size());
	for(auto & infection : infections) {
		if(infection.isActive()) {
			ids.push_back(infection.getCurrentGeneId());
		}
	}
	return ids;
}

int64_t Host::getActiveInfectionImmunityCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			if(immunity.isImmune(infection.getCurrentGene())) {
				count++;
			}
		}
	}
	return count;
}

int64_t Host::getActiveInfectionClinicalImmunityCount()
{
	int64_t count = 0;
	for(auto & infection : infections) {
		if(infection.isActive()) {
			if(clinicalImmunity.isImmune(infection.getCurrentGene())) {
				count++;
			}
		}
	}
	return count;
}

void Host::gainAlleleImmunity(GenePtr genePtr) {
    immunity.gainAlleleImmunity(genePtr,false,*popPtr->simPtr->dbPtr,popPtr->simPtr->alleleImmunityTable);
}

void Host::prepareToDie()
{
//	cerr << popPtr->getTime() << ", host going to die: " << toString() << '\n';
	
	// Remove all events
	for(auto & infection : infections) {
		infection.prepareToEnd();
	}
	immunity.prepareToDie();
	clinicalImmunity.prepareToDie();
	
	removeEvent(deathEvent.get());
}

double Host::moiRegulate(Host & dstHost){
    int64_t maxMOI = popPtr->simPtr->parPtr->withinHost.maxMOI;
    int64_t t = dstHost.getActiveInfectionCount();
    double p =  1-(1/(1+exp(-(t-maxMOI))));
    return(p);
}

void Host::transmitTo(Host & dstHost)
{
	rng_t * rngPtr = getRngPtr();
	
	if(infections.size() == 0) {
	//	cerr << "No infections to transmit" << endl;
        popPtr->simPtr->writeEIR(popPtr->getTime(),0);
		return;
	}
    
    popPtr->simPtr->writeEIR(popPtr->getTime(),1);
    
//	cerr << "Transmitting to " << dstHost.popPtr->id << ", " << dstHost.id << '\n';
	double moiControl = moiRegulate(dstHost);
	// Get some current infections according to transmission probability
	vector<StrainPtr> originalStrains;
    for(auto infItr = infections.begin(); infItr != infections.end(); infItr++) {
		if(infItr->isActive()) {
			bernoulli_distribution flipCoin(infItr->transmissionProbability()*moiControl);
			if(flipCoin(*rngPtr)) {
				originalStrains.push_back(
					infItr->strainPtr);
			}
		}
	}
	
	if(originalStrains.size() == 0) {
		return;
	}
	
	vector<StrainPtr> strainsToTransmit;
	if(originalStrains.size() > 1) {
		// Take each original strain with probability 1.0 - pRecombinant
        // disable pRecombinant parameter, use rule of combinations
        // given n strains, the probability of strains of recombined is 1-(1/n);
		double pRecombinant = 1.0-(1.0/double(originalStrains.size()));
		assert(pRecombinant >= 0.0 && pRecombinant <= 1.0);
		strainsToTransmit = drawMultipleBernoulliRandomSubset(
			*rngPtr, originalStrains, 1.0 - pRecombinant, false
		);
		
		// Complete a set of size originalStrains.size() using recombinants
		int64_t nRecombinants = originalStrains.size() - strainsToTransmit.size();
		for(int64_t i = 0; i < nRecombinants; i++) {
			uniform_int_distribution<int64_t> indDist(0, originalStrains.size() - 1);
			int64_t ind1 = indDist(*rngPtr);
			int64_t ind2 = indDist(*rngPtr);
			strainsToTransmit.push_back(
				popPtr->simPtr->recombineStrains(
					originalStrains[ind1],
					originalStrains[ind2]
				)
			);
		}
	}
	else {
		strainsToTransmit = originalStrains;
	}
	
	// Transmit all strains
	for(auto & strainPtr : strainsToTransmit) {
		dstHost.receiveInfection(strainPtr);
	}
	popPtr->simPtr->recordTransmission(*this, dstHost, strainsToTransmit);
}

void Host::transmitMSTo(Host & dstHost)
{
	rng_t * rngPtr = getRngPtr();
	
	if(infections.size() == 0) {
        //	cerr << "No infections to transmit" << endl;
        popPtr->simPtr->writeEIR(popPtr->getTime(),0);
		return;
	}
    
    popPtr->simPtr->writeEIR(popPtr->getTime(),1);
    
    //	cerr << "Transmitting to " << dstHost.popPtr->id << ", " << dstHost.id << '\n';
	
	// Get some current infections according to transmission probability
	vector<StrainPtr> originalStrains;
    vector<GenePtr> originalMS;
	for(auto infItr = infections.begin(); infItr != infections.end(); infItr++) {
		if(infItr->isActive()) {
			bernoulli_distribution flipCoin(infItr->transmissionProbability());
			if(flipCoin(*rngPtr)) {
				originalStrains.push_back(infItr->strainPtr);
                originalMS.push_back(infItr->msPtr);
			}
		}
	}
	
    
	if(originalStrains.size() == 0) {
		return;
	}
	
	vector<StrainPtr> strainsToTransmit;
    vector<GenePtr> msToTransmit;
	if(originalStrains.size() > 1) {
		// Take each original strain with probability 1.0 - pRecombinant
		//double pRecombinant = popPtr->simPtr->parPtr->pRecombinant;
        double pRecombinant = 1.0-(1.0/double(originalStrains.size()));
		assert(pRecombinant >= 0.0 && pRecombinant <= 1.0);
        std::vector<size_t> indices = drawMultipleBernoulli(*rngPtr, originalStrains.size(), 1.0 - pRecombinant, false);
        strainsToTransmit.reserve(indices.size());
        msToTransmit.reserve(indices.size());
        for(size_t index : indices) {
            strainsToTransmit.push_back(originalStrains[index]);
            msToTransmit.push_back(originalMS[index]);
        }
        
		//strainsToTransmit = drawMultipleBernoulliRandomSubset(
		//	*rngPtr, originalStrains, 1.0 - pRecombinant, false
		//);
		
		// Complete a set of size originalStrains.size() using recombinants
		int64_t nRecombinants = originalStrains.size() - strainsToTransmit.size();
		for(int64_t i = 0; i < nRecombinants; i++) {
			uniform_int_distribution<int64_t> indDist(0, originalStrains.size() - 1);
			int64_t ind1 = indDist(*rngPtr);
			int64_t ind2 = indDist(*rngPtr);
			strainsToTransmit.push_back(
                    popPtr->simPtr->recombineStrains(
                    originalStrains[ind1],
                    originalStrains[ind2]
                    )
            );
            msToTransmit.push_back(
            popPtr->simPtr->recombineMS(                                                                    originalMS[ind1],
                originalMS[ind2])
            );
		}
	}
	else {
		strainsToTransmit = originalStrains;
        msToTransmit = originalMS;
	}
	
	// Transmit all strains
    for(size_t i=0; i<strainsToTransmit.size(); ++i) {
        dstHost.receiveInfection(strainsToTransmit[i],msToTransmit[i]);
    }
    /*
     for(auto & strainPtr : strainsToTransmit) {
     dstHost.receiveInfection(strainPtr,infItr->msPtr);
     }
     */
	popPtr->simPtr->recordTransmission(*this, dstHost, strainsToTransmit);
}

void Host::receiveInfection(StrainPtr & strainPtr)
{
	assert(strainPtr->size() > 0);
	
	rng_t * rngPtr = getRngPtr();
	double time = popPtr->getTime();
	
	double tLiverStage = popPtr->simPtr->parPtr->tLiverStage;
	int64_t initialGeneIndex = tLiverStage == 0 ? 0 : WAITING_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, initialGeneIndex, time);
	
	
	// Get an iterator to the infection so that events can
	// point back to their infection
	list<Infection>::reverse_iterator reverseItr = infections.rbegin();
	reverseItr++;
	list<Infection>::iterator infectionItr = reverseItr.base();
	
	// If starting in liver stage, create a fixed-time transition event
	// (liver stage -> first gene not yet active)
	if(initialGeneIndex == WAITING_STAGE) {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(infectionItr, time + tLiverStage)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	// Otherwise create a rate-based transition event
	// (first gene not yet active -> first gene active)
	else {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
			new TransitionEvent(
				infectionItr,
				infectionItr->transitionRate(),
				time,
				*rngPtr
			)
		);
		addEvent(infectionItr->transitionEvent.get());
	}
	
    //Create a mutation event, rate equals pMutation * genesPerStrain
    infectionItr->mutationEvent = unique_ptr<MutationEvent>(new MutationEvent(infectionItr,popPtr->simPtr->parPtr->pMutation * popPtr->simPtr->parPtr->genesPerStrain * popPtr->simPtr->locusNumber,time,*rngPtr));
    addEvent(infectionItr->mutationEvent.get());

    //Create a ectopic recombination event, rate equals pIntraRecomb * C(genesPerStrain,2)
    int64_t numGenes = popPtr->simPtr->parPtr->genesPerStrain;
    infectionItr->recombinationEvent = unique_ptr<RecombinationEvent>(new RecombinationEvent(infectionItr,popPtr->simPtr->parPtr->pIntraRecomb * numGenes * (numGenes -1)/2,time,*rngPtr));
    addEvent(infectionItr->recombinationEvent.get());

	// Create a clearance event
	// (rate will depend on state as determined in clearanceRate()
	// and may be zero)
    infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
        new ClearanceEvent(
		infectionItr, infectionItr->clearanceRate(), time, *rngPtr
		)
    );
	addEvent(infectionItr->clearanceEvent.get());

	
//	cerr << time << ": " << infectionItr->toString() << " begun" << '\n';
}

void Host::receiveInfection(StrainPtr & strainPtr, GenePtr & msPtr)
{
	assert(strainPtr->size() > 0);
	
	rng_t * rngPtr = getRngPtr();
	double time = popPtr->getTime();
	
	double tLiverStage = popPtr->simPtr->parPtr->tLiverStage;
	int64_t initialGeneIndex = tLiverStage == 0 ? 0 : WAITING_STAGE;
	infections.emplace_back(this, nextInfectionId++, strainPtr, msPtr, initialGeneIndex, time);
	
	
	// Get an iterator to the infection so that events can
	// point back to their infection
	list<Infection>::reverse_iterator reverseItr = infections.rbegin();
	reverseItr++;
	list<Infection>::iterator infectionItr = reverseItr.base();
	
	// If starting in liver stage, create a fixed-time transition event
	// (liver stage -> first gene not yet active)
	if(initialGeneIndex == WAITING_STAGE) {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
               new TransitionEvent(infectionItr, time + tLiverStage)
               );
		addEvent(infectionItr->transitionEvent.get());
	}
	// Otherwise create a rate-based transition event
	// (first gene not yet active -> first gene active)
	else {
		infectionItr->transitionEvent = unique_ptr<TransitionEvent>(
             new TransitionEvent(
              infectionItr,
              infectionItr->transitionRate(),
              time,
              *rngPtr
              )
        );
		addEvent(infectionItr->transitionEvent.get());
	}
	
    //Create a mutation event, rate equals pMutation * genesPerStrain
    infectionItr->mutationEvent = unique_ptr<MutationEvent>(new MutationEvent(infectionItr,popPtr->simPtr->parPtr->pMutation * popPtr->simPtr->parPtr->genesPerStrain * popPtr->simPtr->locusNumber,time,*rngPtr));
    addEvent(infectionItr->mutationEvent.get());
    
    //Create a ectopic recombination event, rate equals pIntraRecomb * C(genesPerStrain,2)
    int64_t numGenes = popPtr->simPtr->parPtr->genesPerStrain;
    infectionItr->recombinationEvent = unique_ptr<RecombinationEvent>(new RecombinationEvent(infectionItr,popPtr->simPtr->parPtr->pIntraRecomb * numGenes * (numGenes -1)/2,time,*rngPtr));
    addEvent(infectionItr->recombinationEvent.get());
    
    // Create an ms mutation event, change rate
    infectionItr->msMutationEvent = unique_ptr<MSmutationEvent>(new MSmutationEvent(infectionItr,popPtr->simPtr->parPtr->pMsMutate * popPtr->simPtr->microsatNumber,time,*rngPtr));
    addEvent(infectionItr->msMutationEvent.get());
    
    
	// Create a clearance event
	// (rate will depend on state as determined in clearanceRate()
	// and may be zero)
    infectionItr->clearanceEvent = unique_ptr<ClearanceEvent>(
        new ClearanceEvent(
        infectionItr, infectionItr->clearanceRate(), time, *rngPtr
        )
    );
	addEvent(infectionItr->clearanceEvent.get());
    
	
    //	cerr << time << ": " << infectionItr->toString() << " begun" << '\n';
}

void Host::updateInfectionRates()
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		if(itr->geneIndex != WAITING_STAGE) {
			itr->updateClearanceRate();
			itr->updateTransitionRate();
		}
	}
}

void Host::clearInfection(std::list<Infection>::iterator infectionItr)
{
	assert(infectionItr != infections.end());
	
	bool shouldUpdateAllRates = infectionItr->transitionAffectsAllInfections();
	
//	double time = popPtr->getTime();
//	cerr << time << ": " << infectionItr->toString() << " clearing" << '\n';
	
	// Gain immunity to active gene
	if(infectionItr->active) {
		GenePtr genePtr = infectionItr->getCurrentGene();
        
        getSelectionMode(genePtr, true);
	}
	
    // add table to record duration of infection
    double durationTime = getTime()-infectionItr->initialTime;
    popPtr->simPtr->writeDuration(infectionItr, durationTime);
    
	// Remove infection
	infectionItr->prepareToEnd();
	infections.erase(infectionItr);
	
	if(shouldUpdateAllRates) {
		updateInfectionRates();
	}
}

void Host::getSelectionMode(GenePtr genePtr, bool clearInfection) {
    int64_t selmode = popPtr->simPtr->parPtr->selectionMode;
    //cout<<"selmode is "<<selmode<<endl;
    if (selmode == 1) {
        // specific immunity mode
        if(!popPtr->simPtr->parPtr->withinHost.useAlleleImmunity) {
            immunity.gainImmunity(genePtr);
            if(getSimulationParametersPtr()->trackClinicalImmunity) {
                clinicalImmunity.gainImmunity(genePtr);
            }
        }else{
            //cout<<"gainAlleleImmunity"<<endl;
            gainAlleleImmunity(genePtr);
            
        }
        
    }else if ((selmode == 2) && clearInfection) {
        // general immunity mode
        immunity.gainGeneralImmunity();
    }
    // if selmode == 3, do nothing
    
}

void Host::hstMutateStrain(std::list<Infection>::iterator infectionItr)
{
    std::vector<GenePtr> newStrainGenes = popPtr->simPtr->mutateStrain(infectionItr->strainPtr);
    std::unordered_map<size_t,size_t> reorderIndexMap = ordered(newStrainGenes);
    for (size_t i=0; i<newStrainGenes.size();i++) {
        infectionItr->expressionOrder[i] = reorderIndexMap[infectionItr->expressionOrder[i]];
    }
    infectionItr->strainPtr = popPtr->simPtr->getStrain(newStrainGenes);
}

void Host::RecombineStrain(std::list<Infection>::iterator infectionItr) {
    std::vector<GenePtr> newStrainGenes = popPtr->simPtr->ectopicRecStrain(infectionItr->strainPtr);
    std::unordered_map<size_t,size_t> reorderIndexMap = ordered(newStrainGenes);
    for (size_t i=0; i<newStrainGenes.size();i++) {
        infectionItr->expressionOrder[i] = reorderIndexMap[infectionItr->expressionOrder[i]];
    }
    infectionItr->strainPtr = popPtr->simPtr->getStrain(newStrainGenes);
}

void Host::microsatMutate(std::list<Infection>::iterator infectionItr) {
    infectionItr->msPtr = popPtr->simPtr->mutateMS(infectionItr->msPtr);
}

double Host::getTime()
{
	return popPtr->getTime();
}

rng_t * Host::getRngPtr()
{
	return popPtr->rngPtr;
}

SimParameters * Host::getSimulationParametersPtr()
{
	return popPtr->simPtr->parPtr;
}

PopulationParameters * Host::getPopulationParametersPtr()
{
	return popPtr->parPtr;
}

void Host::addEvent(Event * event)
{
	popPtr->addEvent(event);
}

void Host::removeEvent(Event * event)
{
	popPtr->removeEvent(event);
}

void Host::setEventRate(RateEvent * event, double rate)
{
	popPtr->setEventRate(event, rate);
}

void Host::writeInfections(Database & db, Table<InfectionRow> & table,Table<StrainRow> & strainsTable,Table<GeneRow> & GeneTable,Table<LociRow> & LociTable)
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		itr->write(db, table,strainsTable,GeneTable,LociTable);
        
	}
}

void Host::writeInfections(int64_t transmissionId, Database & db, Table<TransmissionInfectionRow> & table)
{
	for(auto itr = infections.begin(); itr != infections.end(); itr++) {
		itr->write(transmissionId, db, table);
	}
}

std::string Host::toString()
{
	stringstream ss;
	ss << popPtr->toString() << ".h" << id;
	return ss.str();
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->prepareToDie();
	hostPtr->popPtr->removeHost(hostPtr);
}
