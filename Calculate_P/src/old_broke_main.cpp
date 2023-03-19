#include <dolfin.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>
#include <sstream>
#include <fstream>

#include "helpers.hpp"      // load_txt, display_setup
#include "Space.h"          // compiled UFL file with functions
#include "expressions.hpp"  // kappa, sc, DS1, source_term
#include "boundary.hpp"     // BoundaryX0, BoundaryX1,...
#include "Matrix.hpp"

using namespace dolfin;
#include "load_bin.hpp"

// Uncomment this line to have all of the parameters be
// printed.
//#define DEBUG

int main(int argc, char** argv) {

    dolfin::init(argc,argv);
    // Prints certain settings that are available
    parameters("krylov_solver")["absolute_tolerance"] = 1e-25;
    parameters("krylov_solver")["relative_tolerance"] = 1e-23;
    //std::cout << parameters.str(true) << "\n\n";
    //list_lu_solver_methods();
    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;

    // This can be changed to a config file later, maybe
    double l = 0.014, w = 0.012, h = 0.00065;
    unsigned ls = 129, ws = 129, hs = 14;
    unsigned num_modes = 5;  
    bool status;                    // flag for checking status of things

    std::string floorplan_file = "Floorplan.txt";
    std::vector<std::vector<double>> flp;
    status = Helpers::load_txt(floorplan_file, flp);
    if( !status ) {
        std::cout << "ERROR loading floorplan file: \"" << floorplan_file << "\"" << std::endl;
        return 1;
    }

    // END OF INITIALIZATION PHASE

    // set up geometric model, mesh, and function space
    std::shared_ptr<BoxMesh> mesh = 
        std::make_shared<BoxMesh>(
                BoxMesh(Point(0,0,0), Point(l,w,h), ls-1,ws-1,hs-1)
            );

    unsigned counter = 0;
    double tx, ty, tz;

    auto V = std::make_shared<Space::FunctionSpace>( mesh );

    bool res;
    std::vector<std::vector<std::shared_ptr<Function>>> modes;
    if(mpi_rank == 0) {std::cout << "LOADING:";}
    for( unsigned i = 0; i < num_modes; i++ ) {
        if( mpi_rank == 0 ) {std::cout << "," << i << std::flush;}
        modes.push_back(std::vector<std::shared_ptr<Function>>());
        for(unsigned j = 0; j < flp.size(); j++ ) {
            if(mpi_rank == 0 ) {std::cout << ":" << j << std::flush;}
            auto u = std::make_shared<Function>(V);
            std::stringstream ss;
            ss << "modes/mode_" << i << ".bin.0.bin" << j << ".bin";
            res = load_function_from_bin(ss.str(),"MODE_COORDS.bin.0.bin",
                *u,mpi_rank, comm);
            if( !res ) {
                std::cout << "COULDN'T LOAD " << ss.str() << " EXITING!" << std::endl;
                exit(1);
            }
            modes[i].push_back(std::move(u));
        }
        if(mpi_rank == 0) { std::cout << std::endl; }
    }
    if(mpi_rank == 0) { std::cout << std::endl; }

    ::Matrix P_mat(flp.size(),num_modes);
    Space::Form_P_test P(mesh);
    for( unsigned i = 0; i < flp.size(); i++ ) {
        for( unsigned j = 0; j < num_modes; j++ ) {
            P.u1=modes[j][i];
            if(mpi_rank == 0) {std::cerr << i << ":" << j << ",";}
            P_mat[i][j] = assemble(P);
        }
    }
    if( mpi_rank == 0 ) {P_mat.write_file("P_mat.csv",false);}

    return 0;

}
