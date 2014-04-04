#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include "EventQueue.hpp"

class Host;
class Simulation;

class DeathEvent : public zppsim::OneTimeEvent
{
public:
	DeathEvent(Host * hostPtr);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	Host * hostPtr;
};

class Host
{
friend class Simulation;
public:
	Host(uint64_t id, double deathTime);
	uint64_t const id;
	double const deathTime;
	
	void die(Simulation & sim);
	void transmitTo(Host & dstHost);
private:
	
	// Hash set of genes that this host has immunity to
	std::unordered_set<uint32_t> immunity;
	
	std::unique_ptr<DeathEvent> deathEvent;
};

#endif
