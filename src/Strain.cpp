#include "Strain.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace zppsim;
	
//size_t Strain::nextId(0);
//std::vector<std::weak_ptr<Strain>> Strain::strains;
//zppsim::unordered_map_bh<std::vector<GenePtr>, size_t> Strain::geneVecToIndexMap;

Strain::Strain(std::vector<GenePtr> const & genes) :
	genes(genes)
{
}

size_t Strain::size()
{
	return genes.size();
}

std::vector<GenePtr> Strain::getGenes()
{
	return genes;
}

GenePtr Strain::getGene(size_t index)
{
	return genes[index];
}

/*StrainPtr makeRandomStrain(size_t nGenes, zppsim::rng_t rng)
{
	std::vector<GenePtr> genes;
	for(size_t i = 0; i < nGenes; i++) {
//		genes.push_back(Gene::getRandomGene());
	}
	
	auto sp = std::make_shared<Strain>();
//	genes.push_back(std::weak_ptr<Gene>(gp));
//	idToIndexMap[id] = genes.size() - 1;
	return gp;
}
*/
