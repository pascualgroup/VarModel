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
	uint32_t nVarGenes;
	
	// Simulation time
	double tEnd;
	
	// Per capita biting rate (rate of paired bites)
	double bitingRate;
	
	// Per-capita birth and death rates
	double birthRate;
	double deathRate;
	
	SimParameters(std::istream & paramStream, Database & db)
		: Parameters(paramStream, db, "parameters")
	{
		LOAD_VALUE(randomSeed);
		
		LOAD_VALUE(initialPopulationSize);
		LOAD_VALUE(nVarGenes);
		LOAD_VALUE(tEnd);
		
		LOAD_VALUE(bitingRate);
		LOAD_VALUE(birthRate);
		LOAD_VALUE(deathRate);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
