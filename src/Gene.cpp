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

Gene::Gene(size_t id) :
	id(id)
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
