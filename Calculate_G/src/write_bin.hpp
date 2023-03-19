#ifndef _DOLFIN_WRITE_CSV_HPP_
#define _DOLFIN_WRITE_CSV_HPP_

#include <iostream>
#include <dolfin.h>
#include <fstream>
#include <vector>
#include <sstream>

#include "Matrix.hpp"

using namespace dolfin;

bool write_solution_as_bin(
        const std::string& ofname,
        const std::string& coords_pref,
        Function& f,
        const unsigned& mpi_rank, MPI_Comm& comm,
        const bool save_coords=false
) {
    
    auto nodal_vals = f.vector();
    std::vector<double> arr;
    nodal_vals->get_local(arr); 

    if( save_coords ) {
        //auto coor = V->tabulate_dof_coordinates();
        auto coor = f.function_space()->tabulate_dof_coordinates();
        ::Matrix coors(arr.size(),3);
        unsigned idx = 0;
        for( unsigned i = 0; i < arr.size(); i++ ) {
            idx = 3*i;
            coors[i][0] = coor[idx];
            coors[i][1] = coor[idx+1];
            coors[i][2] = coor[idx+2];
        }
        std::stringstream out_ofname;
        out_ofname << coords_pref << "." << mpi_rank << ".bin";
        coors.write_file(out_ofname.str());
    }

    std::stringstream ss;
    ss << ofname << "." << mpi_rank << ".bin";
    ::Matrix vals(arr.size(),1);
    for( unsigned i = 0; i < arr.size(); i++ ) {
        vals[i][0] = arr[i];
    }
    vals.write_file(ss.str());

    MPI_Barrier(comm);

    return true;
}

#endif
