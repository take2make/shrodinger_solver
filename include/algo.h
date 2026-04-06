#ifndef ALGO_H
#define ALGO_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <iomanip>
#include <complex>
#include <math.h>

struct GridParams {
    int N = 200;
    int Nt = 800;
    double L = 10.0;
    double T = 10.0;
    double dx = L / N;
    double dt = 0.001;
    int stride = 50;
    double m = 10.0;
};
inline GridParams gridParams;

struct GaussParams {
    double A = 1.0;
    double x0 = 5.0;
    double y0 = 5.0;
    double sigma = 1.0;
    double k = 0.5;
};
inline GaussParams gaussParams;

enum class PotentialType {
    None,
    Barrier,
    Barrier2D,
    DoubleSlit2D,
    BarrierStep,
    Harmonic,
};

struct PotentialParams {
    PotentialType type = PotentialType::Barrier;
    double barrier_V0 = 30.0;
    double barrier_a = 4.5;
    double barrier_b = 5.5;
    double harmonic_k = 0.5;
};
inline PotentialParams potentialParams;

// utils for write data:
void        write_data_psi(const std::vector<double>& rho, const std::vector<double>& ci_half);
void        ecrit_new(std::string filename, Eigen::VectorXcd u);
void        ecrit(std::string filename, Eigen::VectorXcd u);

// Initial conditions:
double      gauss(double x, double m, double s);
double      gauss_with_velocity_real(double x, double x0, double s, double k);
double      gauss_with_velocity_imag(double x, double x0, double s, double k);
double      gauss_2D(double x, double y, double x0, double y0, double s);
double      gauss_2D_with_velocity_real(double x, double y, double x0, double y0, double s, double k);
double      gauss_2D_with_velocity_imag(double x, double y, double x0, double y0, double s, double k);

// Type of potentials:
double      barrier_potential(double x, double V0, double a, double b);
double      harmonic_oscillator(double x, double k);
double      barrier_potential_2D(double x, double y, double V0, double a, double b);
double      potential(double x);
double      potential_2D(double x, double y);

// Shrodinger solvers:
void        cranck_nickolson_solver();
void        leapfrog_solver();
void        leapfrog_solver_2D();

#endif