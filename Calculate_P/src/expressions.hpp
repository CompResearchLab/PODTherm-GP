#ifndef _EXPRESSIONS_HPP_
#define _EXPRESSIONS_HPP_

#include<vector>
#include<dolfin.h>
#include<fstream>
#include<iostream>
#include<iomanip>
using namespace dolfin;

class Zeroer : public Expression {
    public:
        std::vector<std::vector<double>> flp;
        const double tol = 1e-14;
        unsigned dont_zero;
	double  h=0, thick_Sio2=0, thick_actl = 0;
	void setParams(double h_in, double thick_Sio2_in, double thick_actl_in){
	
		h   = h_in;
		thick_Sio2 = thick_Sio2_in;
		thick_actl = thick_actl_in;
	} 
        void eval(Array<double>& values, const Array<double>& x) const {
            unsigned n = dont_zero;
	    if((x[2] > (h - thick_Sio2 - thick_actl -tol)) && (x[2] < (h - thick_Sio2 + tol))){
	    	
		    if( (x[0] > flp[n][2] -tol) && (x[0] < (flp[n][0] + flp[n][2] + tol) )&& (x[1] > flp[n][3]-tol) && (x[1] < (flp[n][3] + flp[n][1] + tol))){
		    
			    values[0] = 1.0;
		    }else{
		    
			    values[0] = 0.0;
		    }
	    } else{
	    
		    values[0] = 0.0;
	    }

            /*if( (x[0] > flp[n][2]-tol) && (x[0] < (flp[n][0] + flp[n][2] + tol)) 
                && (x[1] > flp[n][3]-tol) && (x[1] < (flp[n][3] + flp[n][1] + tol)) 
                && (x[2] > (h - thick_Sio2 - thick_actl+tol)) && (x[2] < (h - thick_Sio2 +tol))) {
                values[0] = 1.0;
            } else {
                values[0] = 0.0;
            }*/
        }

        void setFlp( std::vector<std::vector<double>>& flpin ) {
            for( unsigned i = 0; i < flpin.size(); i++) {
                flp.push_back(std::vector<double>());
                for( unsigned j = 0; j < flpin[i].size(); j++ ) {
                    flp[i].push_back( flpin[i][j] );
                }
            }
        }
};

#endif
