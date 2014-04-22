#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>
#include "EventQueue.hpp"
#include "zppsim_util.hpp"

class Host;
class Population;

class DeathEvent : public zppsim::Event
{
public:
	DeathEvent(Host * hostPtr);
	virtual void performEvent(zppsim::EventQueue & queue);
private:
	Host * hostPtr;
};

class Host
{
friend class Population;
friend class DeathEvent;
public:
	Host(Population * popPtr, size_t id, double deathTime);
	
	void die(zppsim::EventQueue & queue);
	void transmitTo(Host & dstHost);
	
	void pushBackEvents(std::vector<zppsim::Event *> & eventVec);
private:
	Population * popPtr;
	size_t id;
	double deathTime;
	
	// Hash set of genes that this host has immunity to
	std::unordered_set<uint32_t> immunity;
	
	std::unique_ptr<DeathEvent> deathEvent;
};

#endif
