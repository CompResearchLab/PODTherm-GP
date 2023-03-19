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
#include "/usr/include/pugixml.hpp"

// These are expression classes for the form equations
class DS1Expression : public Expression {
    public:
        double D0=0.0, D1=0.0, tol=0.0, h=0, thick_Sio2=0;

        void setParams(double tolin, double D0in, double D1in,double h_in, double thick_Sio2_in) { 
            tol = tolin; 
            D0 = D0in; 
            D1 = D1in; 
	    h   = h_in;
	    thick_Sio2 = thick_Sio2_in;
        }
        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = ( (x[2] <= h - thick_Sio2 + tol) ? D0 : D1);
        }    
};

class SCExpression : public Expression {
    public:
        double c0=0.0, c1=0.0, tol=0.0,h=0, thick_Sio2=0;

        void setParams(double tolin, double c0in, double c1in, double h_in, double thick_Sio2_in) { 
            tol = tolin; 
            c0 = c0in; 
            c1 = c1in; 
	    h   = h_in;
	    thick_Sio2 = thick_Sio2_in;
        }
        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = ( (x[2] <= h - thick_Sio2 + tol) ? c0 : c1);
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
     double D0 = std::stod(variables["D0_in"]);
     double D1 = std::stod(variables["D1_in"]);
     double c0 = std::stod(variables["c0_in"]);
     double c1 = std::stod(variables["c1_in"]);
     double chip_area = l*w; 
     bool status;
    std::shared_ptr<BoxMesh> mesh = 
        std::make_shared<BoxMesh>(
                BoxMesh(Point(0,0,0), Point(l,w,h), ls-1,ws-1,hs-1)
            );
    std::cout << "NUM CELLS IN MESH: " << mesh->num_cells() << std::endl;

    unsigned counter = 0;
    double tx, ty, tz;

    auto V = std::make_shared<Space::FunctionSpace>( mesh );
    std::cout << "FUNCTION SPACE DIMS: " << V->dim() << std::endl;
    

    // Here, the modes are all loaded in from HDF5 files
    // into a vector of Function objects.
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
        std::cout << (*u->vector())[0] << std::endl;
        modes.push_back(std::move(u));
    }

    std::shared_ptr<DS1Expression> DS1 = std::make_shared<DS1Expression>();
    DS1->setParams(tol,D0,D1,h,thick_Sio2);
    std::shared_ptr<SCExpression> sc = std::make_shared<SCExpression>();
    sc->setParams(tol,c0,c1,h,thick_Sio2);


    auto u1 = std::make_shared<Function>(V);
    auto u2 = std::make_shared<Function>(V);
    // The Form for the C equation is set up, and the
    // density and specific heat parameters are set.
    Space::Form_C C(mesh);
    ::Matrix C_mat(num_modes,num_modes);
    C.sc = sc;
    C.DS1 = DS1;
    double cellval;
    // This pair of loops runs for each mode, and generates
    // the symmetric C Matrix.
    for( unsigned i = 0; i < modes.size(); i++ ) {
        if( mpi_rank == 0 ) {std::cout << i << std::endl;}
        C.u1=modes[i];

        for( unsigned j = i; j < num_modes; j++ ) {
            if( mpi_rank == 0 ) {std::cout << j << "," << std::flush;}
            C.u2=modes[j];
            cellval = assemble(C);
            C_mat[i][j] = cellval;
            C_mat[j][i] = cellval;
        }
        if( mpi_rank == 0 ) {std::cout << std::endl;}
    }
    // Here, the C Matrix is stored to CSV.
    C_mat.write_file("../../C.csv",false);


    return 0;
}
