//
//  Gene.cpp
//  malariamodel
//
//  Created by Ed Baskerville on 4/4/14.
//  Copyright (c) 2014 Ed Baskerville. All rights reserved.
//

#include "Gene.h"
#include "zppsim_random.hpp"
#include <cassert>
#include <iostream>

using namespace std;
using namespace zppsim;

//std::vector<std::weak_ptr<Gene>> Gene::genes;
//std::unordered_map<Gene *, size_t> Gene::ptrToIndexMap;

/*GenePtr Gene::makeGene(double transmissibility, double meanInfectionDuration, double meanImmunityDuration)
{
	auto gp = std::make_shared<Gene>(transmissibility, meanInfectionDuration, meanImmunityDuration);
	genes.push_back(std::weak_ptr<Gene>(gp));
	ptrToIndexMap[gp.get()] = genes.size() - 1;
	return gp;
}*/

/*static GenePtr getRandomGene(rng_t & rng)
{
	return GenePtr(genes[drawUniformIndex(rng, genes.size())]);
}

size_t Gene::geneCount()
{
	return genes.size();
}*/

Gene::Gene() :
	transmissibility(numeric_limits<double>::signaling_NaN()),
	meanInfectionDuration(numeric_limits<double>::signaling_NaN()),
	meanImmunityDuration(numeric_limits<double>::signaling_NaN())
{
}

Gene::Gene(double transmissibility, double meanInfectionDuration, double meanImmunityDuration) :
	transmissibility(transmissibility),
	meanInfectionDuration(meanInfectionDuration),
	meanImmunityDuration(meanImmunityDuration)
{
}

/*Gene::~Gene()
{
	size_t index = ptrToIndexMap[this];
	ptrToIndexMap.erase(id);
	if(index != genes.size() - 1) {
		genes[index] = genes.back();
	}
	genes.pop_back();
}*/