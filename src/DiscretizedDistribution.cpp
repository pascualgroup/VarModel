#include "DiscretizedDistribution.h"

DiscretizedDistribution::DiscretizedDistribution(std::vector<double> const & pdf, double x0, double dx) :
    discDist(pdf.begin(), pdf.end()),
    realDist(0.0, 1.0),
    x0(x0),
    dx(dx)
{
}

double DiscretizedDistribution::draw(zppsim::rng_t & rng)
{
    return x0 + (discDist(rng) + realDist(rng)) * dx;
}
