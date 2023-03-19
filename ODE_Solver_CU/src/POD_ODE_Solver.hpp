#ifndef _THERMAL_POD_ODE_SOLVER_HPP_
#define _THERMAL_POD_ODE_SOLVER_HPP_

#include<iostream>
#include "/usr/include/pugixml.hpp"
#include "Matrix.hpp"

namespace Thermal {
namespace POD {

class ODE_Solver {
    public:
        ODE_Solver() {
            CU.reset(NPOD,NT);
	    std::cout << "THE NPOD" << NPOD <<"THE NT "<<NT<<std::endl;
        }
	bool init(int argc, char** argv);
        int run();

        ::Matrix CU;
    private:
        unsigned CU_rows;
        unsigned CU_cols;
        bool CU_alloced;
        bool realloc_CU(const unsigned& new_rows, const unsigned& new_cols);
        ::Matrix C, Gmatrix, P_lib;
        double** MFloorplan;
        double** Pmatrix1;

        bool load_file(
            const std::string& ifname, 
            double**& mat, 
            const unsigned& nrows,
            const unsigned& ncols
        );

        /**define parameter of ODE solver
        * NT: is the total time step
        * NU: the number of functional unit
        * NPOD: the maximum of podmode used to predict temperature.
        * NC: the size of C and G matrix*/
	int NT =0;
	int NU =0;
	int NPOD=0;
	int NC =0;
	int N_p =0;
	double T =0.0;
	unsigned num_steps = 0;
	double sampling_interval = 0.0;
	double thick_actl =0.0;
	std::string ptrace_file = "";
	std::string floorplan_file = "";
	void read_xml_config(const std::string& filename) {
		pugi::xml_document doc;
		pugi::xml_parse_result result = doc.load_file(filename.c_str());
		if (!result) {
			std::cerr << "Error loading XML config file: " << result.description() << std::endl;
			return;
		}
		pugi::xml_node variables = doc.child("variables");
		NT = variables.find_child_by_attribute("variable", "name", "num_steps_pre_in").text().as_int();
		num_steps = variables.find_child_by_attribute("variable", "name", "num_steps_in").text().as_int();
		NU = variables.find_child_by_attribute("variable", "name", "Nu_in").text().as_int();
		NPOD = variables.find_child_by_attribute("variable", "name", "num_modes_pre_in").text().as_int();
		NC = variables.find_child_by_attribute("variable", "name", "num_modes_in").text().as_int();
		N_p = variables.find_child_by_attribute("variable", "name", "num_modes_in").text().as_int();
		T= variables.find_child_by_attribute("variable", "name", "T_in").text().as_double();
		thick_actl = variables.find_child_by_attribute("variable", "name", "thick_actl_in").text().as_double();
		ptrace_file = variables.find_child_by_attribute("variable", "name", "Power_path_in").text().as_string();
		floorplan_file = variables.find_child_by_attribute("variable", "name", "Flp_path_in").text().as_string();
		sampling_interval = T/num_steps;
		CU.reset(NPOD,NT);
	}




        /* sampling_interval: the size of samling */
        // sampling interval should be the DT from FEM, same
        // with thick_actl

};

}
}

#endif
