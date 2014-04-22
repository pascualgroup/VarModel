#ifndef __malariamodel__Strain__
#define __malariamodel__Strain__

#include "zppsim_random.hpp"
#include "zppsim_util.hpp"
#include "Gene.h"
#include <memory>
#include <unordered_map>
#include <vector>

class Strain;
typedef std::shared_ptr<Strain> StrainPtr;

class Strain
{
public:
	Strain(std::vector<GenePtr> const & genes);
private:
	std::vector<GenePtr> genes;
	
//	static std::vector<std::weak_ptr<Strain>> strains;
//	static zppsim::unordered_map_bh<std::vector<GenePtr>, size_t> geneVecToIndexMap;
//	
//	static StrainPtr makeRandomStrain(size_t nGenes, zppsim::rng_t rng);
};

#endif
