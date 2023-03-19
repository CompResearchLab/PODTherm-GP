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
#include "boundary.hpp"
#include "/usr/include/pugixml.hpp"
class KExpression : public Expression {
    public:
        double k0=0.0, k1=0.0, tol=0.0,h=0.0, thick_Sio2=0.0;

        void setParams(double tolin, double k0in, double k1in,double h_in, double thick_Sio2_in) {
            tol = tolin;
            k0 = k0in;
            k1 = k1in;
	    h   = h_in;
	    thick_Sio2 = thick_Sio2_in;
        }

        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = (( x[2] <= h - thick_Sio2 + tol) ? k0 : k1);
        }
};


int main(int argc, char** argv) {

    dolfin::init(argc,argv);
    
    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;
    // set up geometric model, mesh, and function space
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
    double Ta = std::stod(variables["Ta_in"]);
    unsigned Nu = std::stoi(variables["Nu_in"]);
    double thick_actl = std::stod(variables["thick_actl_in"]);
    double thick_Sio2 = std::stod(variables["thick_Sio2_in"]);
    double tol = std::stod(variables["tol_in"]);
    double k_0 = std::stod(variables["k_0_in"]);
    double k_1 = std::stod(variables["k_1_in"]);
    double h_c = std::stod(variables["h_c_in"]);
    double chip_area = l*w;
    bool status; 

    std::shared_ptr<BoxMesh> mesh = 
        std::make_shared<BoxMesh>(
                BoxMesh(Point(0,0,0), Point(l,w,h), ls-1,ws-1,hs-1)
            );
    std::cout << "NUM CELLS IN MESH: " << mesh->num_cells() << std::endl;

    auto V = std::make_shared<Space::FunctionSpace>( mesh );
    std::cout << "FUNCTION SPACE DIMS: " << V->dim() << std::endl;
    std::shared_ptr<MeshFunction<size_t>> boundary_markers =std::make_shared<MeshFunction<size_t>>(mesh, 2);
    boundary_markers->set_all(9999);
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

    
    // here, the modes are loaded into a vector of Function
    // objects to be used.
    std::vector<std::shared_ptr<Function>> modes;
    if(mpi_rank == 0) {std::cout << "LOADING:";}
    for( unsigned i = 0; i < num_modes; i++ ) {
        if( mpi_rank == 0 ) {std::cout << i << "," << std::flush;}
        auto u = std::make_shared<Function>(V);
        std::stringstream ss;
        ss << "../../POD_mode/mode_" << i << "h5";
        auto infile = HDF5File(mesh->mpi_comm(), ss.str(), "r");
        infile.read(*u,"solution");
        infile.close();
        modes.push_back(std::move(u));
    }

    // The Form for G is initialized
    Space::Form_G G(mesh);

    // the Expression object for kappa is set up
    std::shared_ptr<KExpression> kappa = std::make_shared<KExpression>();
    kappa->setParams(tol,k_0,k_1,h,thick_Sio2);
    
    double cellval;

    // The G Matrix is initialized to be Nmodes*Nmodes
    auto h_coeff = std::make_shared<Constant>(h_c);
    ::Matrix G_mat(num_modes,num_modes);
    G.kappa = kappa;
    G.h_c_in = h_coeff;
    G.ds = boundary_markers;
    // This loop runs to calculate the G Matrix
    // symmetrically
    for( unsigned i = 0; i < modes.size(); i++ ) {
        if( mpi_rank == 0 ) {std::cout << i << std::endl;}
        G.u1=modes[i];
        for( unsigned j = i; j < num_modes; j++ ) {
            if( mpi_rank == 0 ) {std::cout << j << "," << std::flush;}
            G.u2=modes[j];
            cellval = assemble(G);
            G_mat[i][j] = cellval;
            G_mat[j][i] = cellval;
        }
        if( mpi_rank == 0 ) {std::cout << std::endl;}
    }
    // Save G to CSV
    G_mat.write_file("../../G.csv",false);


    return 0;
}
