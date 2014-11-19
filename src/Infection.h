//
//  Infection.h
//  malariamodel
//
//  Created by Ed Baskerville on 5/5/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__Infection__
#define __malariamodel__Infection__

#include <list>

#include "EventQueue.hpp"
#include "zppsim_util.hpp"
#include "Strain.h"
#include "zppdb.hpp"
#include "DatabaseTypes.h"

class Infection;
class Host;

class InfectionProcessEvent : public zppsim::RateEvent
{
public:
	InfectionProcessEvent(std::list<Infection>::iterator infectionItr, double time);
	InfectionProcessEvent(std::list<Infection>::iterator infectionItr,
		double rate, double time, zppsim::rng_t & rng);
	
	std::list<Infection>::iterator infectionItr;
};

class TransitionEvent : public InfectionProcessEvent
{
public:
	TransitionEvent(std::list<Infection>::iterator infectionItr, double time);
	TransitionEvent(std::list<Infection>::iterator infectionItr,
		double rate, double initTime, zppsim::rng_t & rng);
	virtual void performEvent(zppsim::EventQueue & queue);
};

class ClearanceEvent : public InfectionProcessEvent
{
public:
	ClearanceEvent(std::list<Infection>::iterator infectionItr,
		double rate, double initTime, zppsim::rng_t & rng);
	virtual void performEvent(zppsim::EventQueue & queue);
};

class Infection
{
public:
	Infection(Host * hostPtr, int64_t id, StrainPtr & strainPtr, int64_t initialGeneIndex, double initialTime);
	
	void prepareToEnd();
	
	Host * hostPtr;
	StrainPtr strainPtr;
	
	int64_t id;
	
	int64_t geneIndex;
	bool active;
	double transitionTime;
	
	std::unique_ptr<TransitionEvent> transitionEvent;
	std::unique_ptr<ClearanceEvent> clearanceEvent;
	
	bool isActive();
	GenePtr getCurrentGene();
	int64_t getCurrentGeneId();
	bool isImmune();
	bool isClinicallyImmune();
	double getTransitionTime();
	double getAgeAtTransitionTime();
	
	void performTransition();
	
	void updateTransitionRate();
	bool transitionAffectsAllInfections();
	
	double transitionRate();
	
	void updateClearanceRate();
	double clearanceRate();
	
	double transmissionProbability();
	
	std::string toString();
	
	void write(Database & db, Table<InfectionRow> & table);
	
private:
	double activationRate();
	double deactivationRate();
};

#endif /* defined(__malariamodel__Infection__) */
