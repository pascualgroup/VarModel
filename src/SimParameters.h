#ifndef __malariamodel__SimParameters__
#define __malariamodel__SimParameters__

#define DEFINE_PARAM(x) define(#x, x)

#include "PtreeObject.hpp"
#include "zppsim_random.hpp"
#include <random>

using namespace zppdata;
using namespace zppsim;

class BitingRate : public PtreeObject
{
public:
	double mean;
	double amplitude;
	
	BitingRate()
	{
		DEFINE_PARAM(mean);
		DEFINE_PARAM(amplitude);
	}
};

class PopulationParameters : public PtreeObject
{
public:
	size_t size;
	size_t sampleSize;
	
	size_t nInitialInfections;
	BitingRate bitingRate;
	double immigrationRate;
	
	// Location of population in 2D
	double x;
	double y;
	
	// "Distance" to self for use in contact-weight calculations
	double selfDistance;
	
	PopulationParameters()
	{
		DEFINE_PARAM(size);
		DEFINE_PARAM(sampleSize);
		DEFINE_PARAM(nInitialInfections);
		DEFINE_PARAM(bitingRate);
		DEFINE_PARAM(immigrationRate);
		DEFINE_PARAM(x);
		DEFINE_PARAM(y);
		DEFINE_PARAM(selfDistance);
	}
};

class DistanceFunction : public PtreeObject
{
public:
	double alpha;
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
		DEFINE_PARAM(dt);
		DEFINE_PARAM(pdf);
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

class TransmissionParameters : public PtreeObject
{
public:
	bool coinfectionReducesTransmission;
	
	TransmissionParameters()
	{
		DEFINE_PARAM(coinfectionReducesTransmission);
	}
};

class WithinHostParameters : public PtreeObject
{
public:
	double clearanceRateInitial;
	double clearanceRateMidCourse;
	double activationRate;
	double deactivationRateImmune;
	double deactivationRateNotImmune;
	
	WithinHostParameters()
	{
		DEFINE_PARAM(clearanceRateInitial);
		DEFINE_PARAM(clearanceRateMidCourse);
		DEFINE_PARAM(activationRate);
		DEFINE_PARAM(deactivationRateImmune);
		DEFINE_PARAM(deactivationRateNotImmune);
	}
};

class GeneParameters : public PtreeObject
{
public:
	std::vector<double> transmissibility;
	std::vector<double> immunityLossRate;
	std::vector<double> clinicalImmunityLossRate;
	
	GeneParameters()
	{
		DEFINE_PARAM(transmissibility);
		DEFINE_PARAM(immunityLossRate);
		DEFINE_PARAM(clinicalImmunityLossRate);
	}
};

class SimParameters : public PtreeObject
{
public:
	std::string dbFilename = "output.sqlite";
	double dbCommitPeriod;
	std::map<std::string, bool> dbTablesEnabled;
	
	uint32_t randomSeed = 0;
	
	// Simulation year time scale (used for seasonality period, by default in days)
	double tYear;
	
	// Simulation end time
	double tEnd;
	
	// How often to update the seasonal rates
	double seasonalUpdateEvery;
	
	// How often to sample hosts
	double sampleHostsEvery;
	
	// Number of var genes in global pool
	size_t genePoolSize;
	
	// Number of genes in a strain
	size_t genesPerStrain;
	
	// Probability per gene of a mutation when a strain
	// is picked up
	double pMutation;
	
	// Probability that a transmitted strain will be a recombinant
	double pRecombinant;
	
	// Length of liver stage
	double tLiverStage;
	
	// Parameters controlling distance function
	DistanceFunction distanceFunction;
	
	// Host-lifetime distribution, specified as a discrete PDF,
	// with uniform density within each discrete chunk
	DiscreteDistribution hostLifetimeDistribution;
	
	// Vector of population parameters
	std::vector<PopulationParameters> populations;
	
	// Transmission process parameters
	TransmissionParameters transmission;
	
	// Within-host parameters
	WithinHostParameters withinHost;
	
	// Gene parameters
	GeneParameters genes;
	
	// Switch on/off clinical immunity tracking
	bool trackClinicalImmunity;
	
	SimParameters()
	{
		DEFINE_PARAM(dbFilename);
		DEFINE_PARAM(dbCommitPeriod);
		DEFINE_PARAM(dbTablesEnabled);
		
		DEFINE_PARAM(randomSeed);
		
		DEFINE_PARAM(tYear);
		DEFINE_PARAM(tEnd);
		
		DEFINE_PARAM(seasonalUpdateEvery);
		DEFINE_PARAM(sampleHostsEvery);
		
		DEFINE_PARAM(genePoolSize);
		DEFINE_PARAM(genesPerStrain);
		DEFINE_PARAM(pRecombinant);
		DEFINE_PARAM(tLiverStage);
		
		DEFINE_PARAM(hostLifetimeDistribution);
		
		DEFINE_PARAM(populations);
		
		DEFINE_PARAM(transmission);
		DEFINE_PARAM(withinHost);
		DEFINE_PARAM(genes);
		
		DEFINE_PARAM(trackClinicalImmunity);
	}
};

#endif /* defined(__malariamodel__SimParameters__) */
