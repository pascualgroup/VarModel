//
//  Infection.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "Infection.h"
#include "SimParameters.h"
#include "Host.h"
#include <sstream>
#include <algorithm>

using namespace std;

Infection::Infection(Host * hostPtr, int64_t id, StrainPtr & strainPtr, int64_t initialGeneIndex, double initialTime) :
	hostPtr(hostPtr), id(id), strainPtr(strainPtr),
	geneIndex(initialGeneIndex), active(false),
	transitionTime(initialTime)
{
    for (int64_t i=0; i<strainPtr->size(); i++) {
        expressionOrder.push_back(i);
    }
    std::random_shuffle(expressionOrder.begin(),expressionOrder.end());
    expressionIndex = 0;
}

void Infection::prepareToEnd()
{
	hostPtr->removeEvent(transitionEvent.get());
	hostPtr->removeEvent(clearanceEvent.get());
    //now also remove the mutation Events
    hostPtr->removeEvent(mutationEvent.get());
    //now also remove the recombination events
    hostPtr->removeEvent(recombinationEvent.get());
}

GenePtr Infection::getCurrentGene()
{
	assert(geneIndex != WAITING_STAGE);
	return strainPtr->getGene(geneIndex);
}

int64_t Infection::getCurrentGeneId()
{
	return getCurrentGene()->id;
}

bool Infection::isImmune()
{
	return hostPtr->immunity.isImmune(getCurrentGene());
}

bool Infection::isClinicallyImmune()
{
	return hostPtr->clinicalImmunity.isImmune(getCurrentGene());
}

double Infection::getTransitionTime()
{
	return transitionTime;
}

double Infection::getAgeAtTransitionTime()
{
	return transitionTime - hostPtr->birthTime;
}

void Infection::performTransition()
{
	//SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	transitionTime = hostPtr->getTime();
	
	bool shouldUpdateAllInfections = transitionAffectsAllInfections();
	
	if(geneIndex == WAITING_STAGE) {
		assert(!active);
		geneIndex = expressionOrder[expressionIndex];
	}
	else if(active) {
		assert(expressionIndex != strainPtr->size() - 1);
		GenePtr genePtr = strainPtr->getGene(geneIndex);
        //change it to only gain immunity through the time person gets infected
         /*
        if(!simParPtr->withinHost.useAlleleImmunity) {
            hostPtr->immunity.gainImmunity(genePtr);
            if(hostPtr->getSimulationParametersPtr()->trackClinicalImmunity) {
                hostPtr->clinicalImmunity.gainImmunity(genePtr);
            }
        }else{
            hostPtr->gainAlleleImmunity(genePtr);
        }
         */
        hostPtr->immunity.gainGeneralImmunity();
		expressionIndex++;
		active = false;
	}
	else {
		active = true;
	}
	
	if(shouldUpdateAllInfections) {
		// Every infection's rates need to be updated
		hostPtr->updateInfectionRates();
	}
	else {
		updateTransitionRate();
		updateClearanceRate();
	}
//	cerr << transitionTime << ": " << toString() << " transitioned to " << geneIndex << "(" << getCurrentGene()->toString() << "), " << (active ? "active" : "not yet active") << '\n';
}

void Infection::updateTransitionRate()
{
	hostPtr->setEventRate(transitionEvent.get(), transitionRate());
}

bool Infection::transitionAffectsAllInfections()
{
	if(geneIndex == WAITING_STAGE) {
		return false;
	}
	else {
		return true;
	}
}

double Infection::transitionRate()
{
	assert(geneIndex != WAITING_STAGE);
	
	if(!active) {
		return activationRate();
	}
	else {
		return deactivationRate();
	}
}

void Infection::updateClearanceRate()
{
	hostPtr->setEventRate(clearanceEvent.get(), clearanceRate());
}

double Infection::activationRate()
{
	assert(!active);
//	GenePtr genePtr = getCurrentGene();
//	int64_t geneId = getCurrentGeneId();
//	bool immune = isImmune();
//	bool clinicallyImmune = isClinicallyImmune();
//	double age = getAgeAtTransitionTime();
//	int64_t activeCount = hostPtr->getActiveInfectionCount();
//	vector<GenePtr> activeGenes = hostPtr->getActiveInfectionGenes();
//	vector<int64_t> activeGeneIds = hostPtr->getActiveInfectionGeneIds();
//	int64_t immunityCount = hostPtr->getActiveInfectionImmunityCount();
//	int64_t clinicalImmunityCount = hostPtr->getActiveInfectionClinicalImmunityCount();
	
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	double constant = simParPtr->withinHost.activationRateConstant;
	assert(!std::isnan(constant));
	assert(!std::isinf(constant));
	assert(constant > 0.0);
	double power = simParPtr->withinHost.activationRatePower;
	assert(power <= 0.0);
	assert(!std::isnan(power));
	assert(!std::isinf(power));
	
	double nActiveInfections = hostPtr->getActiveInfectionCount();
	
	if(power == 0.0) {
		return constant;
	}
	
	if(nActiveInfections == 0) {
		return std::numeric_limits<double>::infinity();
	}
	return constant * std::pow(nActiveInfections, power);
}

double Infection::deactivationRate()
{
	assert(active);
	
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	double constant = simParPtr->withinHost.deactivationRateConstant;
	assert(!std::isnan(constant));
	assert(!std::isinf(constant));
	assert(constant > 0.0);
	double power = simParPtr->withinHost.deactivationRatePower;
	assert(!std::isnan(power));
	assert(!std::isinf(power));
	
	double nActiveInfections = hostPtr->getActiveInfectionCount();
	
	return constant * std::pow(nActiveInfections, power);
}

double Infection::clearanceRate()
{
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	// Liver stage
	if(geneIndex == WAITING_STAGE) {
		return 0.0;
	}
	// Gene active: clearance rate depends on immunity
    // ??Ed's bug??
	//else if(!active) {
    else if(active) {
		double clearanceRatePower = simParPtr->withinHost.clearanceRatePower;
		assert(!std::isnan(clearanceRatePower));
		assert(!std::isinf(clearanceRatePower));
		
		double nActiveInfections = hostPtr->getActiveInfectionCount();
		double clearanceRateConstant;
        /*
        if(!simParPtr->withinHost.useAlleleImmunity){
            if(isImmune()) {
                clearanceRateConstant = simParPtr->withinHost.clearanceRateConstantImmune;
            }
        else {
                clearanceRateConstant = simParPtr->withinHost.clearanceRateConstantNotImmune;
            }
        }else{
            double geneImmuneLevel = hostPtr->immunity.checkGeneImmunity(strainPtr->getGene(geneIndex));
            double r1 = simParPtr->withinHost.clearanceRateConstantImmune;
            double r2 = simParPtr->withinHost.clearanceRateConstantNotImmune / (1-geneImmuneLevel);
            if(geneImmuneLevel==1.0) {
                clearanceRateConstant = r1;
            }else if (geneImmuneLevel==0) {
                clearanceRateConstant = simParPtr->withinHost.clearanceRateConstantNotImmune;
            }else{
                clearanceRateConstant = r2;
            }
            //cout<<"clearRate "<<clearanceRateConstant<<endl;
        }
         */
        double geneImmuneLevel = hostPtr->immunity.checkGeneralImmunity();
        double r1 = simParPtr->withinHost.clearanceRateConstantImmune;
        double r2 = simParPtr->withinHost.clearanceRateConstantNotImmune / (1-geneImmuneLevel);
        if(geneImmuneLevel==1.0) {
            clearanceRateConstant = r1;
        }else if (geneImmuneLevel==0) {
            clearanceRateConstant = simParPtr->withinHost.clearanceRateConstantNotImmune;
        }else{
            clearanceRateConstant = r2;
        }
        
		assert(!std::isnan(clearanceRateConstant));
		assert(!std::isinf(clearanceRateConstant));
		assert(clearanceRateConstant > 0.0);
		
		return clearanceRateConstant * std::pow(nActiveInfections, clearanceRatePower);
	}
	// Gene inactive: clearance rate = 0
	else {
		return 0.0;
	}
}

bool Infection::isActive()
{
	return active;
}

double Infection::transmissionProbability()
{
	assert(active);
	
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	GenePtr genePtr = getCurrentGene();
	
	double p = genePtr->transmissibility;
	assert(p > 0.0 && p < 1.0);
	
	if(simParPtr->coinfectionReducesTransmission) {
		p /= hostPtr->getActiveInfectionCount();
	}
	
	return p;
}

std::string Infection::toString()
{
	stringstream ss;
	ss << hostPtr->toString() << ".i" << id;
	return ss.str();
}

void Infection::write(Database & db, Table<InfectionRow> & table)
{
	InfectionRow row;
	row.time = hostPtr->getTime();
	row.hostId = hostPtr->id;
	row.infectionId = id;
	row.strainId = strainPtr->id;
	if(geneIndex == WAITING_STAGE) {
		row.geneIndex.setNull();
		row.active.setNull();
	}
	else {
		row.geneIndex = geneIndex;
		row.active = active;
	}
	db.insert(table, row);
}

void Infection::write(int64_t transmissionId, Database & db, Table<TransmissionInfectionRow> & table)
{
	TransmissionInfectionRow row;
	row.transmissionId = transmissionId;
	row.hostId = hostPtr->id;
	row.infectionId = id;
	row.strainId = strainPtr->id;
	if(geneIndex == WAITING_STAGE) {
		row.geneIndex.setNull();
		row.active.setNull();
	}
	else {
		row.geneIndex = geneIndex;
		row.active = active;
	}
	db.insert(table, row);
}

/*** INFECTION PROCESS EVENTS ***/

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	RateEvent(time), infectionItr(infectionItr)
{
}

InfectionProcessEvent::InfectionProcessEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, zppsim::rng_t & rng
) :
	RateEvent(rate, time, rng), infectionItr(infectionItr)
{
}

TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double time
) :
	InfectionProcessEvent(infectionItr, time)
{
}


TransitionEvent::TransitionEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, zppsim::rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

ClearanceEvent::ClearanceEvent(
	std::list<Infection>::iterator infectionItr, double rate, double time, zppsim::rng_t & rng
) :
	InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

MutationEvent::MutationEvent(
                               std::list<Infection>::iterator infectionItr, double rate, double time, zppsim::rng_t & rng
                               ) :
InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

RecombinationEvent::RecombinationEvent(
                             std::list<Infection>::iterator infectionItr, double rate, double time, zppsim::rng_t & rng
                             ) :
InfectionProcessEvent(infectionItr, rate, time, rng)
{
}

void TransitionEvent::performEvent(zppsim::EventQueue &queue)
{
	// If this is the final deactivation, then it's equivalent to clearing
	if(infectionItr->active
		&& infectionItr->expressionIndex == infectionItr->strainPtr->size() - 1
	) {
		infectionItr->hostPtr->clearInfection(infectionItr);
	}
	// Otherwise, actually perform a transition
	else {
		infectionItr->performTransition();
	}
}

void ClearanceEvent::performEvent(zppsim::EventQueue &queue)
{
	infectionItr->hostPtr->clearInfection(infectionItr);
}

void MutationEvent::performEvent(zppsim::EventQueue &queue)
{
    infectionItr->hostPtr->hstMutateStrain(infectionItr);
    //cout<<"mutated"<<endl;
}

void RecombinationEvent::performEvent(zppsim::EventQueue &queue)
{
    infectionItr->hostPtr->RecombineStrain(infectionItr);
    //cout<<"recomb"<<endl;
}
