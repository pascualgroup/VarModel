#include "DiscretizedDistribution.h"

//fix several bugs, change dx into vectors of dx chunks, hqx

DiscretizedDistribution::DiscretizedDistribution(std::vector<double> const & pdf, double x0, std::vector<double> const & dx) :
    pcDist(dx.begin(),dx.end(),pdf.begin()),
	x0(x0)
{
}

double DiscretizedDistribution::draw(zppsim::rng_t & rng)
{
	//return x0 + discDist(rng) + realDist(rng) * dx; original, wrong because discDist(rng) only returns 0,1,2,3...., should be timed with dx
    return x0 + pcDist(rng);
}
