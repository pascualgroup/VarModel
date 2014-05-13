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
#include "zppdata_util.hpp"

Infection::Infection(Host * hostPtr, size_t id, StrainPtr & strainPtr, size_t initialGeneIndex) :
	hostPtr(hostPtr), id(id), strainPtr(strainPtr),
	geneIndex(initialGeneIndex), active(false)
{
}

void Infection::prepareToEnd()
{
	hostPtr->removeEvent(transitionEvent.get());
	hostPtr->removeEvent(clearanceEvent.get());
}

GenePtr Infection::getCurrentGene()
{
	assert(geneIndex != WAITING_STAGE);
	return strainPtr->getGene(geneIndex);
}

void Infection::performTransition()
{
	if(geneIndex == WAITING_STAGE) {
		assert(!active);
		geneIndex = 0;
	}
	else if(active) {
		assert(geneIndex != strainPtr->size() - 1);
		GenePtr genePtr = strainPtr->getGene(geneIndex);
		hostPtr->immunity.gainImmunity(genePtr);
		if(hostPtr->getSimulationParametersPtr()->trackClinicalImmunity) {
			hostPtr->clinicalImmunity.gainImmunity(genePtr);
		}
		geneIndex++;
		active = false;
	}
	else {
		active = true;
	}
	
	updateTransitionRate();
	updateClearanceRate();
	
	double time = hostPtr->getTime();
//	cerr << time << ": " << toString() << " transitioned to " << geneIndex << "(" << getCurrentGene()->toString() << "), " << (active ? "active" : "not yet active") << '\n';
}

void Infection::updateTransitionRate()
{
	hostPtr->setEventRate(transitionEvent.get(), transitionRate());
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
//	GenePtr gene = getCurrentGene();
//	bool immune = hostPtr->isImmune(gene);
//	double age = hostPtr->getAge();
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	return simParPtr->withinHost.activationRate;
}

double Infection::deactivationRate()
{
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	GenePtr gene = getCurrentGene();
	bool immune = hostPtr->immunity.isImmune(gene);
	if(immune) {
		return simParPtr->withinHost.deactivationRateImmune;
	}
	else {
		return simParPtr->withinHost.deactivationRateNotImmune;
	}
}

double Infection::clearanceRate()
{
	SimParameters * simParPtr = hostPtr->getSimulationParametersPtr();
	
	// Liver stage
	if(geneIndex == WAITING_STAGE) {
		return 0.0;
	}
	// Gene active
	else if(!active) {
		// Before first activation
		if(geneIndex == 0) {
			return simParPtr->withinHost.clearanceRateInitial;
		}
		// Between activations
		else {
			return simParPtr->withinHost.clearanceRateMidCourse;
		}
	}
	// Gene active
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
	
	// TODO: reduce by active infections or all infections?
	if(simParPtr->transmission.coinfectionReducesTransmission) {
		p /= hostPtr->infections.size();
	}
	return p;
}

std::string Infection::toString()
{
	return hostPtr->toString() + ".i" + strprintf("%u", id);
}

void Infection::write(DBTable * table)
{
	if(table != nullptr) {
		DBRow row;
		row.set("time", hostPtr->getTime());
		row.set("hostId", int64_t(hostPtr->id));
		row.set("infectionId", int64_t(id));
		row.set("strainId", int64_t(strainPtr->id));
		if(geneIndex == WAITING_STAGE) {
			row.set_null("geneIndex");
			row.set_null("active");
		}
		else {
			row.set("geneIndex", int64_t(geneIndex));
			row.set("active", int64_t(active));
		}
		table->insert(row);
	}
}

