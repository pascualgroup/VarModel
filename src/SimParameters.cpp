#include "SimParameters.h"
#include <cmath>

using namespace std;

double evaluateSinusoid(Sinusoid & s, double t)
{
    return s.mean * (
        1.0 + s.relativeAmplitude * sin(
            2 * M_PI * ((t / s.period) - s.phase)
        )
    );
}
