#include "algo.h"

double gauss(double x, double m, double s)
{
    double coeff = 1 / (sqrt(2 * M_PI) * s);
    return coeff * exp(-pow(x - m, 2)/(2 * pow(s, 2)));
}

double gauss_2D(double x, double y, double x0, double y0, double s)
{
    return gauss(x, x0, s) * gauss(y, y0, s);
}

double gauss_2D_with_velocity_real(double x, double y, double x0, double y0, double s, double k)
{
    return gauss_2D(x, y, x0, y0, s) * cos(k * y);
}

double gauss_2D_with_velocity_imag(double x, double y, double x0, double y0, double s, double k)
{
    return gauss_2D(x, y, x0, y0, s) * sin(k * y);
}

double gauss_with_velocity_real(double x, double x0, double s, double k)
{
    return gauss(x, x0, s) * cos(k * x);
}

double gauss_with_velocity_imag(double x, double x0, double s, double k)
{
    return gauss(x, x0, s) * sin(k * x);
}