#ifndef __malariamodel__DiscretizedDistribution__
#define __malariamodel__DiscretizedDistribution__

#include <vector>
#include <random>
#include "zppsim_random.hpp"

class DiscretizedDistribution
{
public:
	DiscretizedDistribution(std::vector<double> const & pdf, double dx);
	double draw(zppsim::rng_t & rng);
private:
	std::discrete_distribution<> discDist;
	std::uniform_real_distribution<> realDist;
	double dx;
};

#endif // #ifndef __malariamodel__DiscretizedDistribution__
