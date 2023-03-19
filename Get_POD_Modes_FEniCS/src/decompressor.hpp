#ifndef _DECOMPRESS_HPP_
#define _DECOMPRESS_HPP_

#include <iostream>

class Decompressor {
    public:
        Decompressor(const unsigned rank, const unsigned size) {
            mpi_rank = rank;
            mpi_size = size;
        }
        std::string decompress(const std::string& fn) {
            system(std::string("unxz --stdout "+fn+".xz > " + temp_file).c_str());
            return temp_file;
        }
    private:
        std::string temp_file = "temph5";
        unsigned mpi_rank;
        unsigned mpi_size;

};

#endif
