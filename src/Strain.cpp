#include "Strain.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace zppsim;

Strain::Strain(int64_t id, std::vector<GenePtr> const & genes) :
    id(id), genes(genes)
{
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
