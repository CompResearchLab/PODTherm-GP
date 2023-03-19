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
#include "/usr/include/pugixml.hpp"
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
    //Parse the POD parameters from POD_para.xml file
    std::ifstream input_file("../../POD_Para.xml");
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load(input_file);
    if (!result) {
	    std::cerr << "Failed to parse XML: " << result.description() << std::endl;
	    return 1;
    }
    std::map<std::string, std::string> variables;
    pugi::xml_node root_node = doc.child("variables");
    for (pugi::xml_node var_node = root_node.child("variable"); var_node; var_node = var_node.next_sibling("variable")){
	    std::string name = var_node.attribute("name").value();
	    std::string value = var_node.child_value();
	    variables[name] = value;
    }
    double l = std::stod(variables["l_in"]);
    double w = std::stod(variables["w_in"]);
    double h = std::stod(variables["h_in"]);
    unsigned ls = std::stoi(variables["ls_in"]);
    unsigned ws = std::stoi(variables["ws_in"]);
    unsigned hs = std::stoi(variables["hs_in"]);
    unsigned num_modes = std::stoi(variables["num_modes_in"]);
    double thick_actl = std::stod(variables["thick_actl_in"]);
    double thick_Sio2 = std::stod(variables["thick_Sio2_in"]);
    bool status;

    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;
    std::string floorplan_file = variables["Flp_path_in"];
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

    bool res = true;
    std::vector<std::shared_ptr<Function>> modes;
    if(mpi_rank == 0) {std::cout << "LOADING";}
    for( unsigned i = 0; i < num_modes; i++ ) {
        if( mpi_rank == 0 ) {std::cout << "," << i << std::flush;}
        //modes.push_back(std::vector<std::shared_ptr<Function>>());
        auto u = std::make_shared<Function>(V);
        std::stringstream ss;
        ss << "../../POD_mode/mode_" << i << "h5";
        //res = load_function_from_bin(ss.str(),"MODE_COORDS_0.0.bin",
            //*u,mpi_rank, comm);
        auto infile = HDF5File(mesh->mpi_comm(), ss.str(), "r");
        infile.read(*u,"solution");
        infile.close();
        if( !res ) {
            std::cout << "COULDN'T LOAD " << ss.str() << " EXITING!" << std::endl;
            exit(1);
        }
        modes.push_back(std::move(u));
        //if(mpi_rank == 0) { std::cout << "," << std::flush; }
    }
    if(mpi_rank == 0) { std::cout << std::endl; }
    ::Matrix P_mat(num_modes,flp.size());
    Space::Form_P P(mesh);
    double E_t;
    if(mpi_rank == 0) {std::cout << "MODE: " << std::flush;}
    //auto z = std::make_shared<Zeroer>();
    std::shared_ptr<Zeroer> z = std::make_shared<Zeroer>();
    z->setParams(h,thick_Sio2,thick_actl);
    //auto z = std::make_shared<Zeroer>();
    z->setFlp(flp);
    for( unsigned i = 0; i < num_modes; i++ ) { // for i = mode idx
        P.u1 = modes[i];
        if(mpi_rank == 0) {std::cout << i << "," << std::flush;}
        for( unsigned j = 0; j < flp.size(); j++ ) {    // for j = flp idx
            z->dont_zero = j;
            P.z = z;
            E_t = assemble(P);
            P_mat[i][j] = E_t;
	   // std::cout << "the P(0,0) is :"<<P_mat[0][0]<<std::endl;
	   //exit(1);
        }
    }
    if( mpi_rank == 0 ) {
        std::cout << "\nSaving!" << std::endl;
        P_mat.write_file("../../P_mat.csv",false);
    }

    return 0;

}
