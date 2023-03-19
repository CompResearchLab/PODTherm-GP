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
#include "/usr/include/pugixml.hpp"

using namespace dolfin;

// Uncomment this line to have all of the parameters be
// printed.
//#define DEBUG

int main(int argc, char** argv) {

    dolfin::init(argc,argv);
    parameters("krylov_solver")["absolute_tolerance"]=1e-25;
    parameters("krylov_solver")["relative_tolerance"]=1e-23;
    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);
    //if(mpi_rank == 0){
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
    for (pugi::xml_node var_node = root_node.child("variable"); var_node; var_node = var_node.next_sibling("variable")) {
	    std::string name = var_node.attribute("name").value();
	    std::string value = var_node.child_value();
	    variables[name] = value;
    }
    // Parse done
    // Use these variables
    unsigned N_STATE = std::stoi(variables["N_STATE_in"]);
    double l = std::stod(variables["l_in"]);
    double w = std::stod(variables["w_in"]);
    double h = std::stod(variables["h_in"]);
    unsigned ls = std::stoi(variables["ls_in"]);
    unsigned ws = std::stoi(variables["ws_in"]); 
    unsigned hs = std::stoi(variables["hs_in"]); 
    double T = std::stod(variables["T_in"]);
    unsigned num_steps = std::stoi(variables["num_steps_in"]); 
    double t = std::stod(variables["t_in"]);
    double delt = T/num_steps;
    double h_c = std::stod(variables["h_c_in"]);
    double Ta = std::stod(variables["Ta_in"]);
    unsigned Nu = std::stoi(variables["Nu_in"]); 
    double thick_actl = std::stod(variables["thick_actl_in"]);
    double thick_Sio2 = std::stod(variables["thick_Sio2_in"]);
    double tol = std::stod(variables["tol_in"]);
    double k_0 = std::stod(variables["k_0_in"]);
    double k_1 = std::stod(variables["k_1_in"]);
    double D0 = std::stod(variables["D0_in"]);
    double D1 = std::stod(variables["D1_in"]);
    double c0 = std::stod(variables["c0_in"]);
    double c1 = std::stod(variables["c1_in"]);
    double chip_area = l*w;
    bool status;

    //std::cout << "RANK: " << mpi_rank << "\tSIZE: " << mpi_size << std::endl;
    // Load files as arrays of rows of doubles
    //std::cout << "LOADING DATA" << std::endl;
    std::string ptrace_file = variables["Power_path_in"];
    std::string floorplan_file = variables["Flp_path_in"];
   // }
   // MPI::barrier(comm);
    std::vector<std::vector<double>> pd, flp;
    status = Helpers::load_txt(ptrace_file, pd);
    if( !status ) {
        std::cout << "ERROR loading ptrace file: \"" << ptrace_file << "\"" << std::endl;
        return 1;
    }
    status = Helpers::load_txt(floorplan_file, flp);
    if( !status ) {
        std::cout << "ERROR loading floorplan file: \"" << floorplan_file << "\"" << std::endl;
        return 1;
    }

#ifdef DEBUG
    std::cout << std::fixed << std::setprecision(16);
    std::cout << "l\t\t" << l << std::endl;
    std::cout << "w\t\t" << w << std::endl;
    std::cout << "h\t\t" << h << std::endl;
    std::cout << "ls\t\t" << ls << std::endl;
    std::cout << "ws\t\t" << ws << std::endl;
    std::cout << "hs\t\t" << hs << std::endl;
    std::cout << "T\t\t" << T << std::endl;
    std::cout << "num_steps\t" << num_steps << std::endl;
    std::cout << "t\t\t" << t << std::endl;
    std::cout << "dt\t\t" << delt << std::endl;
    std::cout << "h_c\t\t" << h_c << std::endl;
    std::cout << "Nu\t\t" << Nu << std::endl;
    std::cout << "thick_actl\t" << thick_actl << std::endl;
    std::cout << "chip_area\t" << chip_area << std::endl;
#endif

    // Compute Power Density
    for( unsigned i = 0; i < num_steps; i++ ) {
        for( unsigned j = 0; j < Nu; j++ ) {
            pd[i][j] = pd[i][j]/(flp[j][0]*flp[j][1]*thick_actl);
        }
    }

    // END OF INITIALIZATION PHASE
   //std::cout << "here is okay"<<std::endl;

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
    
    Helpers::display_setup(tol,k_0,k_1,D0,D1,c0,c1);
    
    std::shared_ptr<KappaExpression> kappa = std::make_shared<KappaExpression>();
    kappa->setParams(tol,k_0,k_1,h,thick_Sio2);
    std::shared_ptr<DS1Expression> DS1 = std::make_shared<DS1Expression>();
    DS1->setParams(tol,D0,D1,h,thick_Sio2);
    std::shared_ptr<SCExpression> sc = std::make_shared<SCExpression>();
    sc->setParams(tol,c0,c1,h,thick_Sio2);

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

    // END OF BOUNDARY PHASE


    std::shared_ptr<SourceTerm> f = std::make_shared<SourceTerm>();
    f->setFlp(flp); f->setPd(pd);
    f->setParams(0,h,thick_actl,Nu);
    
    auto u0 = std::make_shared<Constant>(0.0);
    auto u_n = std::make_shared<Function>(V);
    u_n->interpolate(*u0);

    // for using Dirichlet BCs    
    //auto dcval = std::make_shared<Constant>(10.0);
    //DirichletBoundary boundary;
    //DirichletBC bc(V, dcval, bz1);

    std::string solution_file_name = "";

    auto dt_ptr = std::make_shared<Constant>(delt);
    auto r = std::make_shared<Constant>(h_c);
    auto s = std::make_shared<Constant>(Ta);
    auto g0 = std::make_shared<Constant>(0.0);
    auto g1 = std::make_shared<Constant>(0.0);
    auto g2 = std::make_shared<Constant>(0.0);
    auto g3 = std::make_shared<Constant>(0.0);
    auto g4 = std::make_shared<Constant>(0.0);

    Space::LinearForm L(V);
    Space::BilinearForm a(V,V);

    L.sc = sc; L.DS1 = DS1; L.u_n = u_n;
    L.dt_in = dt_ptr; 
    L.f = f;
    //L.g0 = g0; L.g1 = g1; L.g2 = g2; L.g3 = g3; L.g4 = g4;
    //L.s = s; L.r = r;
    L.ds = boundary_markers;
    //L.dx = boundary_markers_dx;
    //L.set_cell_domains(boundary_markers_dx);

    a.sc = sc; a.DS1 = DS1; 
    a.dt_in = dt_ptr;
    a.kappa = kappa; 
    a.r = r;
    a.ds = boundary_markers;
    //a.dx = boundary_markers_dx;
    //a.set_cell_domains(boundary_markers_dx);

    Function u(V);

    Parameters solver_params;
    solver_params.add("linear_solver","gmres");
#ifdef DEBUG
    solver_params.add("print_rhs",true);
    solver_params.add("print_matrix",true);
#endif

    std::stringstream ss;
    unsigned dd = 0;
    for(unsigned i = 0; i < num_steps; i++ ) {
        std::cout << i << std::endl;
	    //if(i>99){break;}
#ifdef DEBUG
        f->printParams();
        sc->printParams();
        kappa->printParams();
        DS1->printParams();
        std::cout << "OTHER VALUES-----------" << std::endl;
        std::cout << "t:\t\t" << t << std::endl;
        std::cout << "r:\t\t" << double(*r) << std::endl;
        std::cout << "dt:\t\t" << double(*dt_ptr) << std::endl;
#endif
        auto solve_start = std::chrono::high_resolution_clock::now();
        solve( a == L,u, solver_params); 
        auto solve_stop = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> solve_elapsed = solve_stop-solve_start;
        std::cout << solve_elapsed.count() << std::endl;
        // update values
        // std::cout << "UPDATING VALUES FOR NEXT RUN" << std::endl;
        // save solution file 
        ss << i;
        solution_file_name = "../../solution/file_" + ss.str() + "h5";
        auto solution_file = HDF5File(mesh->mpi_comm(), solution_file_name, "w");
        solution_file.write(u, "solution");
        solution_file.close();

        t += delt;
        dd += 1;
        f->updateDD(dd);
        L.f = f;     

        *u_n = u;

        ss.str("");

    } // end solver loop

    return 0;

}
