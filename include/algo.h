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
    int N_x = 200;
    int N_steps = 800;
    double L = 10.0;
    double T = 10.0;
    double dx = L / N_x;
    double dt = 0.001;
    double m = 10.0;
};
inline GridParams gridParams;

struct GaussParams {
    double A = 1.0;
    double x0 = 5.0;
    double sigma = 0.5;
    double k = 0.5;
};
inline GaussParams gaussParams;

enum class PotentialType {
    None,
    Barrier,
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
void        ecrit_energy(std::string filename, const Eigen::MatrixXcd& H, const Eigen::VectorXcd& u, bool nouveau);

// Initial conditions:
double      gauss(double x, double m, double s);

// Type of potentials:
double      barrier_potential(double x, double V0, double a, double b);
double      harmonic_oscillator(double x, double k);
double      potential(double x);

// Shrodinger solvers:
void        cranck_nickolson_solver();
void        leapfrog_solver();
void        leapfrog_solver2();

#endif