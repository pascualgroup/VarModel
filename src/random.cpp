#include "random.h"

using namespace std;

#include <vector>
#include <iostream>
#include <cassert>

uint32_t drawUniformIndexExcept(rng_t & rng, uint32_t size, uint32_t except) {
	uniform_int_distribution<uint32_t> unif(0, size - 2);
	uint32_t value = unif(rng);
	if(value >= except) {
		value++;
	}
	return value;
}

vector<uint32_t> drawUniformIndicesExcept(rng_t & rng, uint32_t size, uint32_t count, uint32_t except) {
	uniform_int_distribution<uint32_t> unif(0, size - 2);
	vector<uint32_t> indexes(count);
	for(uint32_t i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			if(indexes[i] >= except) {
				indexes[i]++;
			}
			done = true;
			for(uint32_t j = 0; j < i; j++) {
				if(indexes[i] == indexes[j]) {
					done = false;
					break;
				}
			}
		}
	}
	sort(indexes.begin(), indexes.end());
	return indexes;
}

uint32_t drawUniformIndexExcept(rng_t & rng, uint32_t size, std::vector<uint32_t> except) {
	sort(except.begin(), except.end());
	
	uniform_int_distribution<uint32_t> unif(0, size - 1 - uint32_t(except.size()));
	uint32_t value = unif(rng);
	for(uint32_t exceptVal : except) {
		if(value >= exceptVal) {
			value++;
		}
	}
	return value;
}

vector<uint32_t> drawUniformIndicesExcept(rng_t & rng, uint32_t size, uint32_t count, vector<uint32_t> except) {
	sort(except.begin(), except.end());
	
	uniform_int_distribution<uint32_t> unif(0, size - 2);
	vector<uint32_t> indexes(count);
	for(uint32_t i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			for(uint32_t exceptVal : except) {
				if(indexes[i] >= exceptVal) {
					indexes[i]++;
				}
			}
			done = true;
			for(uint32_t j = 0; j < i; j++) {
				if(indexes[i] == indexes[j]) {
					done = false;
					break;
				}
			}
		}
	}
	sort(indexes.begin(), indexes.end());
	return indexes;
}

vector<uint32_t> drawUniformIndices(rng_t & rng, uint32_t size, uint32_t count) {
	uniform_int_distribution<uint32_t> unif(0, size - 1);
	vector<uint32_t> indexes(count);
	for(uint32_t i = 0; i < count; i++) {
		bool done = false;
		while(!done) {
			indexes[i] = unif(rng);
			done = true;
			for(uint32_t j = 0; j < i; j++) {
				if(indexes[i] == indexes[j]) {
					done = false;
					break;
				}
			}
		}
	}
	sort(indexes.begin(), indexes.end());
	return indexes;
}

std::vector<uint32_t> drawMultipleBernoulli(rng_t & rng, uint32_t size, double p) {
	binomial_distribution<int64_t> countDist(size, p);
	int64_t count = countDist(rng);
	std::vector<uint32_t> indices = drawUniformIndices(rng, size, uint32_t(count));
	assert(indices.size() == count);
	return indices;
}
