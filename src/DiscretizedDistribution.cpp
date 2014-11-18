#include "DiscretizedDistribution.h"

DiscretizedDistribution::DiscretizedDistribution(std::vector<double> const & pdf, double dx) :
	discDist(pdf.begin(), pdf.end()),
	realDist(0.0, 1.0),
	dx(dx)
{
}

double DiscretizedDistribution::draw(zppsim::rng_t & rng)
{
	return discDist(rng) + realDist(rng) * dx;
}
