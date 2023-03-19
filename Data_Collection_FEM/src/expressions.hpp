#ifndef _EXPRESSIONS_HPP_
#define _EXPRESSIONS_HPP_

#include<vector>
#include<dolfin.h>
#include<fstream>
#include<iostream>
#include<iomanip>
using namespace dolfin;

class KappaExpression : public Expression {
    public:
        // make sure to set these, or there will be issues
        double k_0=0.0, k_1=0.0, tol=0.0, h=0, thick_Sio2=0;

        void setParams(double tolin, double k_0in, double k_1in, double h_in, double thick_Sio2_in) {
            tol = tolin;
            k_0 = k_0in;
            k_1 = k_1in;
	    h	= h_in;
	    thick_Sio2 = thick_Sio2_in;
        }
        
        void printParams() {
            std::cout << "kappa Parameters----------------" << std::endl;
            std::cout << "tol:\t\t" << tol << std::endl;
            std::cout << "k_0:\t\t" << k_0 << std::endl;
            std::cout << "k_1:\t\t" << k_1 << std::endl;
        }

        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = ( (x[2] <= h - thick_Sio2 + tol) ? k_0 : k_1);
        }    
};

class DS1Expression : public Expression {
    public:
        // make sure to set these, or there will be issues
        double D0=0.0, D1=0.0, tol=0.0, h=0, thick_Sio2=0;
        
        void setParams(double tolin, double D0in, double D1in,double h_in, double thick_Sio2_in) {
            tol = tolin;
            D0 = D0in;
            D1 = D1in;
	    h   = h_in;
	    thick_Sio2 = thick_Sio2_in;
        }
        
        void printParams() {
            std::cout << "DS1 Parameters----------------" << std::endl;
            std::cout << "tol:\t\t" << tol << std::endl;
            std::cout << "D0:\t\t" << D0 << std::endl;
            std::cout << "D1:\t\t" << D1 << std::endl;
        }

        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = ( (x[2] <= h - thick_Sio2 + tol) ? D0 : D1);
        }    
};

class SCExpression : public Expression {
    public:
        // make sure to set these, or there will be issues
        double c0=0.0, c1=0.0, tol=0.0, h=0, thick_Sio2=0;
        
        void setParams(double tolin, double c0in, double c1in, double h_in, double thick_Sio2_in) {
            tol = tolin;
            c0 = c0in;
            c1 = c1in;
	    h   = h_in;
	    thick_Sio2 = thick_Sio2_in;
        }

        void printParams() {
            std::cout << "SC Parameters----------------" << std::endl;
            std::cout << "tol:\t\t" << tol << std::endl;
            std::cout << "c0:\t\t" << c0 << std::endl;
            std::cout << "c1:\t\t" << c1 << std::endl;
        }

        void eval(Array<double>& values, const Array<double>& x) const {
            values[0] = ( (x[2] <= h - thick_Sio2 + tol) ? c0 : c1);
        }    
};

class SourceTerm : public Expression {
    public:
        double h, thick_actl;
        std::vector<std::vector<double>> flp, pd;
        unsigned dd, Nu;
        const double tol = 1e-14;
        
        void eval(Array<double>& values, const Array<double>& x) const {
            if( (x[2] > (h-2*thick_actl)) && (x[2] < h-thick_actl)) {
                for( unsigned n = 0; n < Nu; n++ ) {
                    if( (x[0] >= flp[n][2]) && (x[0] < (flp[n][0] + flp[n][2] + tol)) 
                     && (x[1] >= flp[n][3]) && (x[1] < (flp[n][3] + flp[n][1] + tol)) ) {
                        values[0] = pd[dd][n];
                    }
                }
            } else {
                values[0] = 0.0;
            }
        }

        void setFlp( std::vector<std::vector<double>>& flpin ) {
            for( unsigned i = 0; i < flpin.size(); i++) {
                flp.push_back(std::vector<double>());
                for( unsigned j = 0; j < flpin[i].size(); j++ ) {
                    flp[i].push_back( flpin[i][j] );
                }
            }
        }
        
        void setPd( std::vector<std::vector<double>>& pdin ) {
            for( unsigned i = 0; i < pdin.size(); i++) {
                pd.push_back(std::vector<double>());
                for( unsigned j = 0; j < pdin[i].size(); j++ ) {
                    pd[i].push_back( pdin[i][j] );
                }
            }

        }

        void setParams(
                unsigned ddin, double hin, 
                double thick_actlin, unsigned Nuin
        ) {
            dd = ddin;
            h = hin;
            thick_actl = thick_actlin;
            Nu = Nuin;
        }

        void updateDD(unsigned& ddin) {
            dd = ddin;
        }

        void printParams() {
            std::cout << "SOURCE TERM VALUES---" << std::endl;
            std::cout << "DD:\t\t" << dd << std::endl;
            std::cout << "h:\t\t" << h << std::endl;
            std::cout << "thick_actl:\t" << thick_actl << std::endl;
            std::cout << "Nu:\t" << Nu << std::endl;
        }

        bool store_pd(const std::string& filename) {
            if( filename.size() < 1 ) {
                return false;
            }

            std::ofstream ofh;
            ofh.open(filename.c_str());
            ofh << std::fixed << std::setprecision(16);

            for( unsigned row = 0; row < pd.size(); ++row ) {
                for( unsigned col = 0; col < pd[row].size(); ++col ) {
                    ofh << pd[row][col];
                    if( col < pd[row].size()-1 ) {
                        ofh << ",";
                    }
                }
                ofh << std::endl;
            }
            ofh.close();
        }
};

#endif
