#include "Strain.h"
#include <iostream>
#include <sstream>
#include <cassert>

using namespace std;

Strain::Strain(uint32_t genePoolSize, uint32_t nGenes, rng_t rng) :
	poolSize(genePoolSize), genes(drawUniformIndices(rng, genePoolSize, nGenes))
{
	cerr << "Created strain" << toJsonString() << endl;
}

std::string Strain::toJsonString()
{
	stringstream ss;
	ss << '[';
	for(size_t i = 0; i < genes.size(); i++) {
		ss << genes[i];
		if(i < genes.size() - 1) {
			ss << ", ";
		}
	}
	ss << ']';
	return ss.str();
}
