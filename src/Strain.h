#ifndef __malariamodel__Strain__
#define __malariamodel__Strain__

#include "zppsim_random.hpp"

class Strain
{
public:
	Strain(uint32_t genePoolSize, uint32_t nGenes, zppsim::rng_t rng);
private:
	std::vector<uint32_t> const genes;
};

#endif
