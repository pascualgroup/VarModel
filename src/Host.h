#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include <list>
#include "EventQueue.hpp"
#include "Strain.h"
#include "Gene.h"
#include "zppsim_util.hpp"

#define INFECTION_STAGE_NULL (std::numeric_limits<size_t>::max())

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

class TransitionEvent : public zppsim::Event
{
public:
	TransitionEvent(std::list<Infection>::iterator infectionItr, double time);
	virtual void performEvent(zppsim::EventQueue & queue);
	
	std::list<Infection>::iterator infectionItr;
};

class Infection
{
public:
	Infection(Host * hostPtr, StrainPtr & strainPtr, size_t initialStage);
	
	Host * hostPtr;
	StrainPtr strainPtr;
	size_t stage;
	std::unique_ptr<TransitionEvent> nextTransition;
};

class Host
{
friend class Population;
friend class DeathEvent;
public:
	Host(Population * popPtr, size_t id, double deathTime);
	
	void die(zppsim::EventQueue & queue);
	void transmitTo(Host & dstHost, zppsim::rng_t & rng, double pRecombination);
	
	void receiveInfection(StrainPtr & strain);
	void performTransition(std::list<Infection>::iterator infectionItr);
private:
	Population * popPtr;
	size_t id;
	double deathTime;
	
	// Linked list of current infections
	std::list<Infection> infections;
	
	// Hash set of genes that this host has immunity to
	zppsim::unordered_set_bh<GenePtr> immunity;
	
	std::unique_ptr<DeathEvent> deathEvent;
	
	double calculateTransitionDelay(size_t fromStage);
};

#endif
