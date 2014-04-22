//
//  Gene.h
//  malariamodel
//
//  Created by Ed Baskerville on 4/4/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#ifndef __malariamodel__Gene__
#define __malariamodel__Gene__

#include <memory>
#include <vector>
#include <unordered_map>
#include "zppsim_random.hpp"

class Gene;
typedef std::shared_ptr<Gene> GenePtr;

class Gene
{
public:
	Gene();
	Gene(double transmissibility, double meanInfectionDuration, double meanImmunityDuration);
	
	double const transmissibility;
	double const meanInfectionDuration;
	double const meanImmunityDuration;
	
//	static GenePtr makeGene(double transmissibility, double meanInfectionDuration, double meanImmunityDuration);
//	static GenePtr getRandomGene(zppsim::rng_t & rng);
//	static size_t geneCount();
	~Gene();
private:
	static std::unordered_map<Gene *, size_t> ptrToIndexMap;
	static std::vector<std::weak_ptr<Gene>> genes;
};

#endif /* defined(__malariamodel__Gene__) */
