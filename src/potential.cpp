#include "algo.h"

double harmonic_oscillator(double x, double k)
{
    // harmonic oscillator: return 0.5 * k * x * x;
    return 0.5 * k * x * x;  // example: harmonic
}

double barrier_potential(double x, double V0, double a, double b)
{
    // barrier: return (x > a && x < b) ? V0 : 0.0;
    return (x > a && x < b) ? V0 : 0.0;
}

double potential(double x)
{
    switch (potentialParams.type) {
        case PotentialType::None:
            return 0.0;
        case PotentialType::Barrier:
            return barrier_potential(
                x,
                potentialParams.barrier_V0,
                potentialParams.barrier_a,
                potentialParams.barrier_b
            );
        case PotentialType::Harmonic:
            return harmonic_oscillator(x, potentialParams.harmonic_k);
        default:
            return 0.0;
    }
}
