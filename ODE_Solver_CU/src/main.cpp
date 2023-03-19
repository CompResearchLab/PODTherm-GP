#include <iostream>

#include "POD_ODE_Solver.hpp"

int main(int argc, char** argv) {
    Thermal::POD::ODE_Solver solver;
    std::cout << "Here is Okay" <<std::endl;
    bool status = solver.init(argc, argv);
    if( !status ) {
        std::cout << "COULDN'T INIT SOLVER! EXITING" << std::endl;
        exit(1);
    }
    solver.run();


}
