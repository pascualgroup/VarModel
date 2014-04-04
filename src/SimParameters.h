#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#define DEFINE(x) define(#x, x)

#include "PtreeObject.hpp"

using namespace zppdata;

class SimParameters : public PtreeObject
{
public:
	std::string dbFilename = "output.sqlite";
	
	uint32_t randomSeed = 0;
	
	// Number of populations
	uint32_t nPopulations;
	
	// Initial population size for each population
	std::vector<uint32_t> populationSize;
	
	// Number of var genes in global pool
	uint32_t genePoolSize;
	
	// Simulation end time
	double tEnd;
	
	// Per-capita biting rate, for each population
	std::vector<double> bitingRate = {1.0};
	
	// Relative connectivity between populations,
	// so that contactWeights[i] is a vector that determines
	// the relative probability that individuals
	// in different populations will be the source of a bite;
	// the probability of choosing a source in population j
	// will be populationSize[j] * contactWeight[i][j]
	std::vector<std::vector<double>> contactWeight = {{1.0}};
	
	// Per-capita introduction rate for each population
	std::vector<double> introductionRate;
	
	SimParameters()
	{
		DEFINE(dbFilename);
		DEFINE(randomSeed);
		DEFINE(nPopulations);
		DEFINE(populationSize);
		DEFINE(genePoolSize);
		DEFINE(tEnd);
		DEFINE(bitingRate);
		DEFINE(contactWeight);
		DEFINE(introductionRate);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
