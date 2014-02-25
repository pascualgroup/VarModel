#include "Host.h"

#include <iostream>

using namespace std;

Host::Host(uint64_t id, double tBirth, double tDeath)
	: id(id), tBirth(tBirth), tDeath(tDeath)
{
	cerr << "Created host " << id << ", " << tBirth << ", " << tDeath << endl;
}
