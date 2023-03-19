#ifndef _HELPERS_HPP_
#define _HELPERS_HPP_

#include <vector>
#include <iostream>

namespace Helpers {

    bool load_txt(const std::string& fname, std::vector<std::vector<double>>& mat);

    void display_setup(
            double tol, 
            double k0, double k1, 
            double d0, double d1, 
            double c0, double c1);

}


#endif
