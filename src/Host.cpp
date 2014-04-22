#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;
using namespace zppsim;

Host::Host(Population * popPtr, size_t id, double deathTime)
	: popPtr(popPtr), id(id), deathTime(deathTime), deathEvent(new DeathEvent(this))
{
	cerr << "Created host " << id << ", deathTime " << deathTime << '\n';
}

void Host::die(zppsim::EventQueue & queue)
{
	cerr << queue.getTime() << ", death event: " << popPtr->id << ", " << id << '\n';
	
	id = popPtr->nextHostId++;
	deathTime = popPtr->simPtr->getTime() + popPtr->simPtr->drawHostLifetime();
	
	cerr << "new id: " << id << '\n';
	cerr << "new death time: " << deathTime << '\n';
	
	deathEvent->setTime(queue, deathTime);
}

void Host::transmitTo(Host & dstHost)
{
}

void Host::pushBackEvents(std::vector<zppsim::Event *> & eventVec)
{
	eventVec.push_back(deathEvent.get());
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), Event(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	hostPtr->die(queue);
}
