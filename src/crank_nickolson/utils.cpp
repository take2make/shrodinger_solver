#include "algo.h"

void ecrit_new(std::string filename, Eigen::VectorXcd u){
    std::ofstream file;
    file.open(filename, std::ios::out);
    if (! file){
        std::cerr << "Impossible d'ouvrir " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    file << std::setiosflags(std::ios::scientific) << std::setprecision(7);
    for (int i = 0; i < u.size(); i++){
        file << std::norm(u(i)) << " ";
    }
    file << std::endl;
    file.close();
}

void ecrit(std::string filename, Eigen::VectorXcd u){
    std::ofstream file(filename, std::ios_base::app | std::ios_base::out);
    if (! file){
        std::cerr << "Impossible d'ouvrir " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    file << std::setiosflags(std::ios::scientific) << std::setprecision(7);
    for (int i = 0; i < u.size(); i++){
        file << std::norm(u(i)) << " ";
    }
    file << std::endl;
    file.close();
}

void ecrit_energy(std::string filename, std::complex<double> a, bool nouveau)
{   
    std::ofstream file(filename, nouveau ?std::ios_base::out: std::ios_base::app);
    if (! file){
        std::cerr << "Impossible d'ouvrir " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    file << std::setiosflags(std::ios::scientific) << std::setprecision(7);
    file << a << std::endl; 
    file << std::endl;
    file.close();
}

