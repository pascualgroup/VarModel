#ifndef __malariamodel__DiscretizedDistribution__
#define __malariamodel__DiscretizedDistribution__

#include <vector>
#include <random>
#include "zppsim_random.hpp"

class DiscretizedDistribution
{
public:
	DiscretizedDistribution(std::vector<double> const & pdf, double x0, std::vector<double> const & dx);
	double draw(zppsim::rng_t & rng);
private:
	//std::uniform_real_distribution<> realDist;
	std::piecewise_constant_distribution<double> pcDist;
    double x0;
};

#endif // #ifndef __malariamodel__DiscretizedDistribution__
