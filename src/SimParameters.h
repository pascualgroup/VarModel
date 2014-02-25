#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#include "Parameters.h"

class SimParameters : Parameters
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
	
	// Lifetime gamma-distribution parameters
	double lifetimeMean;
	double lifetimeShape;
	
	SimParameters(std::istream & paramStream, Database & db)
		: Parameters(paramStream, db, "parameters")
	{
		LOAD_VALUE(randomSeed);
		
		LOAD_VALUE(initialPopulationSize);
		LOAD_VALUE(genePoolSize);
		LOAD_VALUE(tEnd);
		
		LOAD_VALUE(bitingRate);
		LOAD_VALUE(lifetimeMean);
		LOAD_VALUE(lifetimeShape);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
