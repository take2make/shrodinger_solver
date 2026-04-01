#include "algo.h"

int main()
{
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.dt = 0.01;
    gridParams.N_steps = int(gridParams.T / gridParams.dt);
    potentialParams.type = PotentialType::Barrier;
    potentialParams.barrier_V0 = 30.0;
    potentialParams.barrier_a = 4.5;
    potentialParams.barrier_b = 5.5;
    
    std::cout << "dt: " << gridParams.dt << ", N_steps: " << gridParams.N_steps << std::endl;
    std::cout << "T real: " << gridParams.N_steps * gridParams.dt << std::endl;
    cranck_nickolson_solver();
    //leapfrog_solver();
    return 0;
}

int main2() {
    gridParams.m = 0.6;
    gridParams.dt = 0.4 * gridParams.m * gridParams.dx * gridParams.dx;
    gridParams.N_steps = int(gridParams.T / gridParams.dt / 100);
    leapfrog_solver2();
    return 0;
}