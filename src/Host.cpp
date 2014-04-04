#include "Host.h"
#include "Simulation.h"

#include <iostream>

using namespace std;

Host::Host(uint64_t id, double deathTime)
	: id(id), deathTime(deathTime), deathEvent(new DeathEvent(this))
{
	cerr << "Created host " << id << ", deathTime " << deathTime << endl;
}

void Host::transmitTo(Host & dstHost)
{
}

DeathEvent::DeathEvent(Host * hostPtr):
	hostPtr(hostPtr), OneTimeEvent(hostPtr->deathTime)
{
}

void DeathEvent::performEvent(zppsim::EventQueue & queue)
{
	Simulation * simPtr = (Simulation *)queue.getContext();
	simPtr->markForDeath(*hostPtr);
}
