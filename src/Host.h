#ifndef __malariamodel__Host__
#define __malariamodel__Host__

#include <unordered_set>

class Host
{
public:
	Host(uint64_t id, double tBirth, double tDeath);
	uint64_t const id;
	double const tBirth;
	double const tDeath;
private:
	// Hash set of genes that this host has immunity to
	std::unordered_set<uint32_t> immunity;
};

#endif
