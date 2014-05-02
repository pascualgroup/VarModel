#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include <list>
#include "EventQueue.hpp"
#include "Strain.h"
#include "Gene.h"
#include "zppsim_util.hpp"

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
	
	Host * hostPtr;
	StrainPtr strainPtr;
	
	size_t id;
	size_t currentGeneIndex;
	bool active;
	
	std::unique_ptr<TransitionEvent> transitionEvent;
	std::unique_ptr<ClearanceEvent> clearanceEvent;
	
	double activationRate(size_t geneIndex);
	double deactivationRate(size_t geneIndex);
	double clearanceRate();
};

class Host
{
friend class Population;
friend class DeathEvent;
public:
	Host(Population * popPtr, size_t id, double deathTime);
	
	void die();
	void transmitTo(Host & dstHost, zppsim::rng_t & rng, double pRecombination);
	
	void receiveInfection(StrainPtr & strain);
	
	void performTransition(std::list<Infection>::iterator infectionItr);
	void clearInfection(std::list<Infection>::iterator infectionItr);
private:
	Population * popPtr;
	size_t id;
	double deathTime;
	
	size_t nextInfectionId;
	
	// Linked list of current infections
	std::list<Infection> infections;
	
	// Hash set of genes that this host has immunity to
	zppsim::unordered_set_bh<GenePtr> immunity;
	
	std::unique_ptr<DeathEvent> deathEvent;
};

#endif
