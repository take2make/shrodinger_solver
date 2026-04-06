#include "algo.h"

double harmonic_oscillator(double x, double k)
{
    return 0.5 * k * x * x;  // example: harmonic
}

double barrier_potential(double x, double V0, double a, double b)
{
    return (x > a && x < b) ? V0 : 0.0;
}

double barrier_potential_step(double x, double V0, double a)
{
    return (x <= a) ? V0 : 0.0;
}

double barrier_potential_2D(double x, double y, double V0, double a, double b)
{
    return (x >= a && x <= b && y >= a && y <= b) ? 0.0 : V0;
}

double double_slit_potential_2D(double x, double y, double pos, double L, double shift, double width)
{
    if (abs(y - pos) < 0.1) {
        if (x >= L/2 - shift/2 && x <= L/2 - shift/2 + width) return 0.0;
        if (x >= L/2 + shift/2 - width && x <= L/2 + shift/2) return 0.0;
        return 300.0;
    }
    return 0.0;
}

double potential(double x)
{
    switch (potentialParams.type) {
        case PotentialType::None:
            return 0.0;
        case PotentialType::Barrier:
            return barrier_potential(x, potentialParams.barrier_V0, potentialParams.barrier_a, potentialParams.barrier_b);
        case PotentialType::BarrierStep:
            return barrier_potential_step(x, potentialParams.barrier_V0, potentialParams.barrier_a);
        case PotentialType::Harmonic:
            return harmonic_oscillator(x, potentialParams.harmonic_k);
        default:
            return 0.0;
    }
}

double potential_2D(double x, double y)
{
    switch (potentialParams.type) {
        case PotentialType::None:
            return 0.0;
        case PotentialType::Barrier2D:
            return barrier_potential_2D(x, y, potentialParams.barrier_V0, potentialParams.barrier_a, potentialParams.barrier_b);
        case PotentialType::DoubleSlit2D:
            return double_slit_potential_2D(x, y, potentialParams.barrier_a, gridParams.L, 2.0, 2.0 / gaussParams.k);
        default:
            return 0.0;
    }
}
