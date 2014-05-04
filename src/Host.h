#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include <list>
#include "EventQueue.hpp"
#include "Strain.h"
#include "Gene.h"
#include "zppsim_util.hpp"
#include "SimParameters.h"

#define LIVER_STAGE (std::numeric_limits<size_t>::max())

class Host;
class Population;
class Infection;

class DeathEvent : public zppsim::Event
{
public:
	DeathEvent(Host * hostPtr);
	virtual void performEvent(zppsim::EventQueue & queue);
	
	Host * hostPtr;
};

class ImmunityLossEvent : public zppsim::RateEvent
{
public:
	ImmunityLossEvent(Host * hostPtr, GenePtr genePtr, double rate, double initTime);
	virtual void performEvent(zppsim::EventQueue & queue);
	
	Host * hostPtr;
	GenePtr genePtr;
};

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
	Infection(Host * hostPtr, size_t id, StrainPtr & strainPtr, size_t initialGeneIndex);
	
	void prepareToEnd();
	
	Host * hostPtr;
	StrainPtr strainPtr;
	
	size_t id;
	size_t geneIndex;
	bool active;
	
	std::unique_ptr<TransitionEvent> transitionEvent;
	std::unique_ptr<ClearanceEvent> clearanceEvent;
	
	bool isActive();
	GenePtr getCurrentGene();
	
	void performTransition();
	
	void updateTransitionRate();
	double transitionRate();
	
	void updateClearanceRate();
	double clearanceRate();
	
	double transmissionProbability();
	
	std::string toString();
	
private:
	double activationRate();
	double deactivationRate();
};

class Host
{
friend class Population;
friend class DeathEvent;
public:
	Host(Population * popPtr, size_t id, double birthTime, double deathTime);
	
	void prepareToDie();
	
	void transmitTo(Host & dstHost);
	
	void receiveInfection(StrainPtr & strain);
	
	void gainImmunity(GenePtr genePtr);
	void loseImmunity(GenePtr genePtr);
	
	void updateInfectionRates();
	
	bool isImmune(GenePtr gene);
	double getAge();
	size_t infectionCount();
	
	void performTransition(std::list<Infection>::iterator infectionItr);
	void clearInfection(std::list<Infection>::iterator infectionItr);
	
	double getTime();
	zppsim::rng_t * getRngPtr();
	
	SimParameters * getSimulationParametersPtr();
	PopulationParameters * getPopulationParametersPtr();
	
	void addEvent(zppsim::Event * event);
	void removeEvent(zppsim::Event * event);
	void setEventRate(zppsim::RateEvent * event, double rate);
	
	std::string toString();
private:
	Population * popPtr;
	
	size_t const id;
	double const birthTime;
	double const deathTime;
	
	size_t nextInfectionId;
	
	// Linked list of current infections
	std::list<Infection> infections;
	
	// Hash set of genes that this host has immunity to
	zppsim::unordered_set_bh<GenePtr> immunity;
	zppsim::unordered_map_bh<GenePtr, std::unique_ptr<ImmunityLossEvent>> immunityLossEvents;
	
	std::unique_ptr<DeathEvent> deathEvent;
};

#endif
