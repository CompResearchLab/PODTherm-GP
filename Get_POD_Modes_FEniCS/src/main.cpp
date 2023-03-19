#include <dolfin.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>
#include <sstream>
#include <fstream>
#include <cmath>

using namespace dolfin;
#include "load_bin.hpp"
#include "write_bin.hpp"
#include "helpers.hpp"      // load_txt, display_setup
#include "Space.h"          // compiled UFL file with functions
#include "Matrix.hpp"
#include "boundary.hpp"
#include "decompressor.hpp"
#include "/usr/include/pugixml.hpp" 
//#define BEAM
//#define COMPRESSION

int main(int argc, char** argv) {

    dolfin::init(argc,argv);
    
    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;
    // Parse the POD parameters from POD_para.xml file
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
    // Use these variables 
    unsigned N_STATE = std::stoi(variables["N_STATE_in"]);
    double l = std::stod(variables["l_in"]);
    double w = std::stod(variables["w_in"]);
    double h = std::stod(variables["h_in"]);
    unsigned ls = std::stoi(variables["ls_in"]);
    unsigned ws = std::stoi(variables["ws_in"]);
    unsigned hs = std::stoi(variables["hs_in"]);
    unsigned num_modes = std::stoi(variables["num_modes_in"]);
    unsigned num_steps = std::stoi(variables["num_steps_in"]);
    double Ta = std::stod(variables["Ta_in"]);
    unsigned Nu = std::stoi(variables["Nu_in"]);
    double thick_actl = std::stod(variables["thick_actl_in"]);
    double thick_Sio2 = std::stod(variables["thick_Sio2_in"]); 
    double tol = std::stod(variables["tol_in"]);
    double chip_area = l*w;
    bool status;

    // This can be changed to a config file later, maybe
#ifdef BEAM
    num_steps = 0;
    ::Matrix t_steps;
    t_steps.read_file("training_steps.csv",false);
    num_steps = t_steps.getNRows();
    std::cout << "NUM BEAM STEPS: " << num_steps << std::endl;
#endif
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
    double cellval;

    // The temperature solutions are loaded into a vector of
    // shared pointers to function objects.
    std::vector<std::shared_ptr<Function>> temps;
#ifdef COMPRESSION
    Decompressor dc(mpi_rank,mpi_size);
#endif
    if(mpi_rank == 0) {std::cout << "LOADING:";}
    unsigned step_num;
    for( unsigned i = 0; i < num_steps; i++ ) {

#ifdef BEAM
        step_num = t_steps[i][0];
#else
        step_num = i;
#endif
        std::stringstream ss;
        ss << "../../solution/file_" << step_num << "h5";
        if( mpi_rank == 0 ) {std::cout << i << ":" << step_num << "|" << ss.str() << "," << std::endl;}
        auto u = std::make_shared<Function>(V);
#ifdef COMPRESSION
        if( mpi_rank == 0 ) {dc.decompress(ss.str());}
        auto infile = HDF5File(mesh->mpi_comm(), "temph5", "r");
#endif
        auto infile = HDF5File(mesh->mpi_comm(), ss.str(), "r");
        infile.read(*u, "solution");
        infile.close();
        temps.push_back(std::move(u));
        //temps.push_back(u);
    }
    // The eigenvalues and eigenvectors are read into
    // Matrix objects.
    std::cout << "Loaded: " << temps.size() << " Temperature files!" << std::endl;
    ::Matrix eigenvals;
    eigenvals.read_file("../../eigenvals.csv",false);
    std::cout << "Eigenvals: " << eigenvals.getNRows() << " x " << eigenvals.getNCols() << std::endl;
    ::Matrix eigenmat;
    eigenmat.read_file("../../eigenmat.csv",false);
    std::cout << "Eigenmat: " << eigenmat.getNRows() << " x " << eigenmat.getNCols() << std::endl;


    std::shared_ptr<MeshFunction<size_t>> boundary_markers =
        std::make_shared<MeshFunction<size_t>>(mesh, 2);
    boundary_markers->set_all(9999);
    
    // setup boundary objects
    std::shared_ptr<BoundaryX0> bx0 = std::make_shared<BoundaryX0>();
    bx0->setLwh(l,w,h);
    std::shared_ptr<BoundaryX1> bx1 = std::make_shared<BoundaryX1>();
    bx1->setLwh(l,w,h);
    std::shared_ptr<BoundaryY0> by0 = std::make_shared<BoundaryY0>();
    by0->setLwh(l,w,h);
    std::shared_ptr<BoundaryY1> by1 = std::make_shared<BoundaryY1>();
    by1->setLwh(l,w,h);
    std::shared_ptr<BoundaryZ0> bz0 = std::make_shared<BoundaryZ0>();
    bz0->setLwh(l,w,h);
    std::shared_ptr<BoundaryZ1> bz1 = std::make_shared<BoundaryZ1>();
    bz1->setLwh(l,w,h);

    bx0->mark(*boundary_markers, 0);
    bx1->mark(*boundary_markers, 1);
    by0->mark(*boundary_markers, 2);
    by1->mark(*boundary_markers, 3);
    bz0->mark(*boundary_markers, 4);
    bz1->mark(*boundary_markers, 5);

    // A matrix for storing the modes is created.
    // Each row corresponds to a mode, and there should be
    // ls*ws*hs columns for each value within the mode.
    ::Matrix mode_mat(num_modes,ls*hs*ws);

    std::cout << num_steps << std::endl;
    // This is performed for each mode.
    auto v2d = vertex_to_dof_map(*V);
    for( unsigned i = 0; i < num_modes; i++ ) {
        // Initialize a black function for the mode.
        auto a = std::make_shared<Function>(V);
        a->interpolate(Constant(0));
        // Multiply each jth element of the ith eigenvector by the jth temperature solution
        // and add them to a
        for( unsigned j = 0; j < num_steps; j++ ) {
            a->vector()->axpy(eigenmat[i][j], *(temps[j]->vector()));
        }
        // divide each element of a by the ith eigenvalue*Nt
        *(a->vector()) /= (eigenvals[i][0]*num_steps);
        // This sets up Form N for normalizing the modes.
        Space::Form_N N(mesh);
        N.ds = boundary_markers;
        // u1 here corresponds to the mode.
        N.u1 = a;
        // Form N is assembled, and the square root of that
        // is taken.
        double sqrtval = sqrt(assemble(N));
        // Each element of the mode is divided by the square
        // root value.
        *a->vector() /= sqrtval;

        // Here, the mode values are copied into the mode
        // matrix.
        unsigned j =0;
        for( unsigned z = 0; z < hs; z++ ) {
            for( unsigned y = 0; y < ws; y++ ) {
                for( unsigned x = 0; x < ls; x++ ) {
                    unsigned idx = x+(ls*y)+(ls*ws*z);
                    mode_mat[i][j] = (*a->vector())[v2d[idx]];
                    ++j;
                }
            }
        }

        // save mode to HDF5
        std::stringstream ss;
        ss << "../../POD_mode/mode_" << i << "h5";
        auto mode_file = HDF5File(mesh->mpi_comm(), ss.str(), "w");
        mode_file.write(*a, "solution");
        mode_file.close(); 
    }
    // save mode matrix to CSV
    mode_mat.write_file("../../modes.csv",false);

    return 0;
}
