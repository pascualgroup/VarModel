#include "Strain.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace zppsim;

Strain::Strain(int64_t id, std::vector<GenePtr> const & genes,bool writeToDatabase,Database & db, Table<StrainRow> & strainsTable) :
	id(id), genes(genes),recorded(writeToDatabase)
{
    if (writeToDatabase) {
        StrainRow row;
        row.strainId = id;
        for(int64_t i = 0; i < size(); i++) {
            row.geneIndex = i;
            row.geneId = getGene(i)->id;
            db.insert(strainsTable, row);
        }
        
    }
}

int64_t Strain::size()
{
	return genes.size();
}

std::vector<GenePtr> Strain::getGenes()
{
	return genes;
}

GenePtr Strain::getGene(int64_t index)
{
	return genes[index];
}

void Strain::writeToDatabaseStrain(Database & db, Table<StrainRow> & strainsTable,Table<GeneRow> & GeneTable,Table<LociRow> & LociTable) {
    if (! recorded) {
        StrainRow row;
        row.strainId = id;
        for(int64_t i = 0; i < size(); i++) {
            row.geneIndex = i;
            row.geneId = getGene(i)->id;
            db.insert(strainsTable, row);
            genes[i]->writeToDatabaseGene(db, GeneTable, LociTable);
        }
        recorded = true;
    }

}
