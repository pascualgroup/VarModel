#ifndef __malariamodel__random__
#define __malariamodel__random__

#include <vector>
#include <random>

typedef std::mt19937_64 rng_t;

uint32_t drawUniformIndexExcept(rng_t & rng, uint32_t size, uint32_t except);
std::vector<uint32_t> drawUniformIndicesExcept(rng_t & rng, uint32_t size, uint32_t count, uint32_t except);
uint32_t drawUniformIndexExcept(rng_t & rng, uint32_t size, std::vector<uint32_t> except);
std::vector<uint32_t> drawUniformIndicesExcept(rng_t & rng, uint32_t size, uint32_t count, std::vector<uint32_t> except);

std::vector<uint32_t> drawUniformIndices(rng_t & rng, uint32_t size, uint32_t count);
std::vector<uint32_t> drawMultipleBernoulli(rng_t & rng, uint32_t size, double p);

#endif
