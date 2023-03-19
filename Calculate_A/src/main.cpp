#include <dolfin.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>
#include <sstream>
#include <fstream>

using namespace dolfin;
#include "load_bin.hpp"
#include "helpers.hpp"      // load_txt, display_setup
#include "Space.h"          // compiled UFL file with functions
#include "Matrix.hpp"
#include "decompressor.hpp"
#include "/usr/include/pugixml.hpp"

int main(int argc, char** argv) {

    dolfin::init(argc,argv);
    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;
    std::ifstream input_file("../../POD_Para.xml");
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(input_file);
    if (!result) {
	    std::cerr << "Failed to parse XML: " << result.description() << std::endl;
	    return 1;
    }
    std::map<std::string, std::string> variables; 
    pugi::xml_node root_node = doc.child("variables"); 
    for (pugi::xml_node var_node = root_node.child("variable"); var_node; var_node = var_node.next_sibling("variable"))
    {
	    std::string name = var_node.attribute("name").value(); 
	    std::string value = var_node.child_value();
	    variables[name] = value;
    }
    // Use these variable
    unsigned N_STATE = std::stoi(variables["N_STATE_in"]);
    double l = std::stod(variables["l_in"]);
    double w = std::stod(variables["w_in"]); 
    double h = std::stod(variables["h_in"]);
    unsigned ls = std::stoi(variables["ls_in"]);
    unsigned ws = std::stoi(variables["ws_in"]);
    unsigned hs = std::stoi(variables["hs_in"]);
    unsigned num_steps = std::stoi(variables["num_steps_in"]);
    double Ta = std::stod(variables["Ta_in"]);
    unsigned Nu = std::stoi(variables["Nu_in"]);
    double thick_actl = std::stod(variables["thick_actl_in"]);
    double thick_Sio2 = std::stod(variables["thick_Sio2_in"]); 
    double chip_area = l*w;
    bool status; 
    double tol = std::stod(variables["tol_in"]);

    // set up geometric model, mesh, and function space
    std::shared_ptr<BoxMesh> mesh = 
        std::make_shared<BoxMesh>(
                BoxMesh(Point(0,0,0), Point(l,w,h), ls-1,ws-1,hs-1)
            );
    std::cout << "NUM CELLS IN MESH: " << mesh->num_cells() << std::endl;

    unsigned counter = 0;
    double tx, ty, tz;

    auto V = std::make_shared<Space::FunctionSpace>( mesh );
    std::cout << "FUNCTION SPACE DIMS: " << V->dim() << std::endl;
    Space::Form_A A(mesh); // Form object for the equation for each cell of the A Matrix

    // Pointers to the functions that will be used for the A
    // equation.
    auto u1 = std::make_shared<Function>(V); 
    auto u2 = std::make_shared<Function>(V);

    // Matrix object to store the A Matrix.
    // Dimensions are Nt*Nt
    ::Matrix A_mat(num_steps,num_steps);
    double cellval; // Temporary value to store each cell of A Matrix

    // The following for loop loads each temperature
    // solution into the us vector.
    // The us vector stores shared pointers to function
    // objects that represent the temperature solutions.
    std::vector<std::shared_ptr<Function>> us;
    if(mpi_rank == 0) {std::cout << "LOADING:";}
    Decompressor dc(mpi_rank, mpi_size);
    for( unsigned i = 0; i < num_steps; i++ ) {
        if( mpi_rank == 0 ) {std::cout << i << "," << std::flush;}
        // First, a blank Function object is created.
        auto u = std::make_shared<Function>(V);
        std::stringstream ss; // ss is used to create the filename from the index of the time step
        ss << "../../solution/file_" << i << "h5";
        // The next three lines load the temperature
        // solution from the HDF5 file.
        //if( mpi_rank == 0 ) {dc.decompress(ss.str());}
        //MPI::barrier(comm);
        //auto infile = HDF5File(mesh->mpi_comm(), "/tmp/temph5", "r");
        auto infile = HDF5File(mesh->mpi_comm(), ss.str(), "r");
        infile.read(*u, "solution");
        infile.close();
        us.push_back(std::move(u));
    }
    std::cout << std::endl;

    // Here a pair of for loops construct the A Matrix.
    for( unsigned i = 0; i < num_steps; i++ ) {
        if( mpi_rank == 0 ) {std::cout << i << std::endl;}
        // u1 is set to the ith time step temperature solution
        A.u1=us[i];

        for( unsigned j = 0; j <= i; j++ ) {
            if( mpi_rank == 0 ) {std::cout << j << "," << std::flush;}
            // u2 is set to the jth time step temperature solution
            A.u2=us[j];
            // The A equation is solved, and the result is
            // divided by Nt to give us the value of A[i,j]
            // and A[j,i], as the A Matrix is symmetric.
            cellval = assemble(A)/num_steps;
            A_mat[i][j] = cellval;
            A_mat[j][i] = cellval;
        }
        if( mpi_rank == 0 ) {std::cout << std::endl;}
       // if( i % 100 == 0 ) {
            if( mpi_rank == 0 ) {
                std::cout << "SAVING A MATRIX AT STEP: " << i << std::endl;
                A_mat.write_file("A.csv",false);
            }
            MPI::barrier(comm);
       // }
    }
    // Store the A Matrix to CSV
    A_mat.write_file("../../A.csv",false);

    return 0;
}
