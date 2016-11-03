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
	int64_t const id;
	
	Strain(int64_t id, std::vector<GenePtr> const & genes, bool writeToDatabase,Database & db, Table<StrainRow> & strainsTable);
	int64_t size();
	std::vector<GenePtr> getGenes();
	GenePtr getGene(int64_t index);
    bool recorded;
    void writeToDatabaseStrain(Database & db, Table<StrainRow> & strainsTable,Table<GeneRow> & GeneTable,Table<LociRow> & LociTable);
private:
	std::vector<GenePtr> genes;
};

#endif
