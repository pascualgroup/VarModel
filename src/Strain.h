#ifndef __malariamodel__Strain__
#define __malariamodel__Strain__

#include "random.h"

class Strain
{
public:
	Strain(uint32_t genePoolSize, uint32_t nGenes, rng_t rng);
	std::string toJsonString();
private:
	uint32_t const poolSize;
	std::vector<uint32_t> const genes;
};

#endif
