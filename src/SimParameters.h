#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#define DEFINE(x) define(#x, x)

#include "PtreeObject.hpp"
#include "zppsim_random.hpp"
#include <random>

using namespace zppdata;
using namespace zppsim;

class PopulationParameters : public PtreeObject
{
public:
	size_t size;
	size_t nInitialInfections;
	double bitingRate;
	std::vector<double> contactWeight;
	
	PopulationParameters()
	{
		DEFINE(size);
		DEFINE(nInitialInfections);
		DEFINE(bitingRate);
		DEFINE(contactWeight);
	}
};

class DiscreteDistribution : public PtreeObject
{
public:
	double dt;
	std::vector<double> pdf;
	
	std::unique_ptr<std::discrete_distribution<>> discDist;
	std::unique_ptr<std::uniform_real_distribution<>> realDist;
	
	DiscreteDistribution() :
		discDist(nullptr),
		realDist(nullptr)
	{
		DEFINE(dt);
		DEFINE(pdf);
	}
	
	double draw(rng_t & rng)
	{
		if(discDist == nullptr) {
			discDist = std::unique_ptr<std::discrete_distribution<>>(
				new std::discrete_distribution<>(pdf.begin(), pdf.end())
			);
			realDist = std::unique_ptr<std::uniform_real_distribution<>>(
				new std::uniform_real_distribution<>(0.0, 1.0)
			);
		}
		return ((*discDist)(rng) + (*realDist)(rng)) * dt;
	}
};

class SimParameters : public PtreeObject
{
public:
	std::string dbFilename = "output.sqlite";
	
	uint32_t randomSeed = 0;
	
	// Simulation end time
	double tEnd;
	
	// Vector of population parameters
	std::vector<PopulationParameters> populations;
	
	// Host-lifetime distribution
	DiscreteDistribution hostLifetimeDistribution;
	
	// Number of var genes in global pool
	size_t genePoolSize;
	
	// Number of genes in a strain
	size_t strainSize;
	
	// Relative connectivity between populations,
	// so that contactWeights[i] is a vector that determines
	// the relative probability that individuals
	// in different populations will be the source of a bite;
	// the probability of choosing a source in population j
	// will be proportional to populationSize[j] * contactWeight[i][j]
	// This matrix need not be symmetrical, although it probably
	// makes sense if it is.
//	std::vector<std::vector<double>> contactWeight = {{1.0}};
	
	// Per-capita introduction rate for each population
	std::vector<double> introductionRate;
	
	SimParameters()
	{
		DEFINE(dbFilename);
		DEFINE(randomSeed);
		DEFINE(tEnd);
		
		DEFINE(populations);
		DEFINE(hostLifetimeDistribution);
		DEFINE(genePoolSize);
		DEFINE(strainSize);
		
		DEFINE(introductionRate);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
