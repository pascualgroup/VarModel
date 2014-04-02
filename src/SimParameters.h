#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#define DEFINE(x) define(#x, x)

#include "PtreeObject.hpp"

using namespace zppdata;

class SimParameters : public PtreeObject
{
public:
	uint32_t randomSeed;
	
	// Initial population size
	uint32_t initialPopulationSize;
	
	// Number of var genes in global pool
	uint32_t genePoolSize;
	
	// Simulation end time
	double tEnd;
	
	// Per-capita biting rate (rate of paired bites)
	double bitingRate;
	
	// Per-capita introduction rate
	double introductionRate;
	
	// Lifetime gamma-distribution parameters
	double lifetimeMean;
	double lifetimeShape;
	
	SimParameters()
	{
		DEFINE(randomSeed);
		DEFINE(initialPopulationSize);
		DEFINE(genePoolSize);
		DEFINE(tEnd);
		DEFINE(bitingRate);
		DEFINE(lifetimeMean);
		DEFINE(lifetimeShape);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
