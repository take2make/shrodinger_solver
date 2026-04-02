#include "algo.h"

double gauss(double x, double m, double s)
{
    double coeff = 1 / (sqrt(2 * M_PI) * s);
    return coeff * exp(-pow(x - m, 2)/(2 * pow(s, 2)));
}


double gauss_with_velocity_real(double x, double m, double s, double k)
{
    return gauss(x, m, s) * cos(k * x);
}

double gauss_with_velocity_imag(double x, double m, double s, double k)
{
    return gauss(x, m, s) * sin(k * x);
}