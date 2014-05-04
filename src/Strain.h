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
typedef std::weak_ptr<Strain> StrainPtrW;

class Strain
{
friend class Simulation;
public:
	Strain(std::vector<GenePtr> const & genes);
	size_t size();
	std::vector<GenePtr> getGenes();
	GenePtr getGene(size_t index);
private:
	std::vector<GenePtr> genes;
};

#endif
