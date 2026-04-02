#include "algo.h"

int test_with_same_dt()
{
    gridParams.m = 0.6;
    gridParams.dt = 0.001;
    gridParams.T = 0.5;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    
    std::cout << "Running leapfrog and Crank-Nicholson solvers with the same dt...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    cranck_nickolson_solver();
    leapfrog_solver();
    return 0;
}

void test_gauss_with_velocity_with_barrier_step() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.x0 = 7.0;
    gaussParams.sigma = 0.6;
    gaussParams.k = -4.0;
    potentialParams.type = PotentialType::BarrierStep;
    potentialParams.barrier_a = 5.5;

    std::cout << "Running leapfrog solver with Gaussian initial condition with step potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver();
}

void test_gauss_without_potential_leapfrog() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::None;

    std::cout << "Running leapfrog solver with Gaussian initial condition without potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver();
}

void test_gauss_without_potential_cranck_nicholson() {
    gridParams.m = 0.6;
    gridParams.dt = 0.01;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::None;

    std::cout << "Running Crank-Nicholson solver with Gaussian initial condition without potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    cranck_nickolson_solver();
}

void test_leapfrog_2D() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::None;

    std::cout << "Running leapfrog solver in 2D with Gaussian initial condition without potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver_2D();
}

int main()
{
    //test_gauss_without_potential_leapfrog();
    //test_gauss_without_potential_cranck_nicholson();
    //test_leapfrog_2D();
    //test_with_same_dt();
    //test_gauss_with_velocity();
    return 0;
}