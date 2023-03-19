#include "helpers.hpp"

#include <fstream>
#include <sstream>
#include <string>

namespace Helpers {

    bool load_txt(
        const std::string& fname, 
        std::vector<std::vector<double>>& mat
    ) {

        // check file, make sure it's good
        std::ifstream ifh;
        ifh.open(fname);
        if( !ifh ) {
            ifh.close();
            return false;
        }

        // clear mat vectors
        unsigned i = 0;
        for( i=0; i < mat.size(); i++ ) {
            mat[i].clear();
        }
        mat.clear();
   
        std::string temp;
        double dtemp;
        i = 0;
        while( std::getline(ifh, temp) ) { // while we can read a line from the file
            mat.push_back(std::vector<double>()); // set up new row in the matrix
            std::stringstream ss(temp); // load row into stringstream
            temp = "";
            while( ss >> dtemp ) { // read values from stringstream
                mat[i].push_back(dtemp); // put values into matrix
            }
            i++; // increment row counter

        }
        return true;

    } // load_txt
    
    void display_setup(
            double tol, 
            double k0, double k1, 
            double d0, double d1, 
            double c0, double c1
    ) {
        using std::cout;
        using std::endl;

        cout << "TOLERANCE:\t" << tol << endl;
        cout << "Thermal Conductivity:\nk_0:\t\t"
            << k0 << "\nk_1:\t\t" << k1 << endl;
        cout << "Density:\nD0:\t\t" 
            << d0 << "\nD1:\t\t" << d1 << endl;
        cout << "Specific Heat:\nc0:\t\t" 
            << c0 << "\nc1:\t\t" << c1 << endl;
    }

} // namespace Helpers
