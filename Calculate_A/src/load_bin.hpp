#ifndef _DOLFIN_LOAD_CSV_HPP_
#define _DOLFIN_LOAD_CSV_HPP_

#include <dolfin.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "Matrix.hpp"

bool load_function_from_bin(
    const std::string& ifname, const std::string& coords_file,
    Function& f, 
    const unsigned& mpi_rank, MPI_Comm& comm
) {
    ::Matrix vals;
    if(!vals.read_file(ifname)) {
        return false;
    }
    ::Matrix coords_in;
    if(!coords_in.read_file(coords_file)) {
        return false;
    }

    std::map<double,std::map<double,std::map<double,double>>> tree;
    double x,y,z;
    for(unsigned i = 0; i < coords_in.getNRows(); i++ ) {
        x = coords_in[i][0];
        y = coords_in[i][1];
        z = coords_in[i][2];
        if( tree.count(x) == 0 ) {
            tree[x] = std::map<double,std::map<double,double>>();
        }
        if( tree[x].count(y) == 0 ) {
            tree[x][y] = std::map<double,double>();
        }
        tree[x][y][z] = vals[i][0];
    }

    std::vector<double> f_vals;
    f.vector()->get_local(f_vals);
    auto coor = f.function_space()->tabulate_dof_coordinates();

    bool found;
    double x_needed, y_needed, z_needed;
    double x_test, y_test, z_test;
    double x_diff, y_diff, z_diff;
    double tol = 1e-7;
    for( unsigned i = 0; i < coor.size(); i+=3 ) {
        x_needed = coor[i];
        y_needed = coor[i+1];
        z_needed = coor[i+2];

        found = false;
        if( tree.find(x_needed) != tree.end() ) {
            if( tree[x_needed].find(y_needed) != tree[x_needed].end() ) {
                if( tree[x_needed][y_needed].find(z_needed) != tree[x_needed][y_needed].end() ) {
                    f_vals[i/3] = tree[x_needed][y_needed][z_needed];
                    continue;
                }
            }
        }

        // THIS IS KEPT IN CASE THE ELEMENT IS NOT FOUND IN
        // THE TREE
        for( unsigned j = 0; j < coords_in.getNRows(); j++ ) {
            x_test = coords_in[j][0];
            y_test = coords_in[j][1];
            z_test = coords_in[j][2];

            x_diff = x_test-x_needed;
            y_diff = y_test-y_needed;
            z_diff = z_test-z_needed;

            // abs_value
            if( x_diff < 0.0 ) { x_diff = -1*x_diff; }
            if( y_diff < 0.0 ) { y_diff = -1*y_diff; }
            if( z_diff < 0.0 ) { z_diff = -1*z_diff; }

            if( ( x_diff < tol ) &&
                ( y_diff < tol ) &&
                ( z_diff < tol )
            ) {
                f_vals[i/3] = vals[j][0];
                found = true;
                break;
            }
        }
        if( !found ) {
            return false;
        }
    }

    f.vector()->set_local(f_vals);
    MPI_Barrier(comm);

    return true;
}

#endif
