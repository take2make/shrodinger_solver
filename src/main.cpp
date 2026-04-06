#include "algo.h"

int test_with_same_dt()
{
    gridParams.m = 0.6;
    gridParams.dt = 0.001;
    gridParams.T = 5.0;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::None;
    
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
    potentialParams.barrier_V0 = 10.0;
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
    gridParams.dt = 0.001;
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
    gridParams.T = 2.0;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.sigma = 0.5;
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::None;

    std::cout << "Running leapfrog solver in 2D with Gaussian initial condition without potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver_2D();
}

void test_leapfrog_2D_with_velocity() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.T = 2.0;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.sigma = 1.0;
    gaussParams.k = 5.0;
    potentialParams.type = PotentialType::None;

    std::cout << "Running leapfrog solver in 2D with Gaussian initial condition without potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver_2D();
}

void test_leapfrog_2D_with_potential() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.T = 3.0;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.sigma = 0.8;
    gaussParams.k = 0.0;
    potentialParams.type = PotentialType::Barrier2D;
    potentialParams.barrier_V0 = 2.0;
    potentialParams.barrier_a = 3.0;
    potentialParams.barrier_b = 7.0;

    std::cout << "Running leapfrog solver in 2D with Gaussian initial condition with barrier potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver_2D();
}

void test_leapfrog_double_slit_potential() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.T = 3.0;
    gridParams.Nt = int(gridParams.T / gridParams.dt);
    gaussParams.x0 = 5.0;
    gaussParams.y0 = 2.0;
    gaussParams.sigma = 1.0;
    gaussParams.k = 10.0;
    potentialParams.type = PotentialType::DoubleSlit2D;
    potentialParams.barrier_a = 4.0;

    std::cout << "Running leapfrog solver in 2D with Gaussian initial condition with barrier potential...\n" << std::endl;
    std::cout << "dt: " << gridParams.dt << ", Nt: " << gridParams.Nt << std::endl;
    std::cout << "T real: " << gridParams.Nt * gridParams.dt << std::endl;
    leapfrog_solver_2D();
}

int main()
{
    //test_gauss_without_potential_leapfrog();
    //test_gauss_without_potential_cranck_nicholson();
    //test_gauss_without_potential_cranck_nicholson();
    //test_with_same_dt();
    //test_gauss_with_velocity_with_barrier_step();

    //test_leapfrog_2D();
    //test_leapfrog_2D_with_velocity();
    //test_leapfrog_2D_with_potential();
    test_leapfrog_double_slit_potential();
    return 0;
}