#include "Strain.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;
using namespace zppsim;

Strain::Strain(uint32_t genePoolSize, uint32_t nGenes, rng_t rng) :
	genes(drawUniformIndices(rng, genePoolSize, nGenes))
{
}
