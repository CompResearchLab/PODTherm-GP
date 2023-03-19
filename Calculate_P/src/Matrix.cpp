#include "Matrix.hpp"

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<iomanip>

void Matrix::init() {
    if( mat.size() != 0 ) {
        for( unsigned i = 0; i < mat.size(); i++ ) {
            mat[i].clear();
        }
        mat.clear();
    }
    for( unsigned i = 0; i < nrows; i++ ) {
        mat.push_back(std::vector<double>());
        for(unsigned j = 0; j < ncols; j++ ) {
            mat[i].push_back(0.0);
        }
    }
}

Matrix Matrix::times(Matrix &other) {
    if(ncols != other.getNRows()) {
        return Matrix(0,0);
    }
    Matrix out(nrows, other.ncols);

    // for every row of the left matrix
    for(unsigned row = 0; row < nrows; row++) {
        // for every col of the right matrix
        for(unsigned col = 0; col < other.ncols; col++) {
            double temp = 0.0;
            // for every cell in those row/col
            for(unsigned i = 0; i < ncols; i++) {
                temp += mat[row][i]*other.mat[i][col];
            }
            out.mat[row][col]=temp;
        }
    }

    return out;
}

void Matrix::print() const {
    for(unsigned i = 0; i < nrows; i++) {
        for(unsigned j = 0; j < ncols; j++) {
            std::cout << mat[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}

double Matrix::getMax() {
    double max = mat[0][0];
    for(unsigned i = 0; i < nrows; i++ ) {
        for(unsigned j = 0; j < ncols; j++ ) {
            if( mat[i][j] > max ) {
                max = mat[i][j];
            }
        }
    }
    return max;
}

bool Matrix::write_file(
        const std::string& ofname,
        bool use_binary,
        char delim,
        unsigned precision
) {
    if( use_binary ) { // binary format
        std::ofstream ofh;
        ofh.open(ofname, std::ios::out | std::ios::binary );
        if( !ofh.good() ) { return false; }
        // metadata
        ofh.write( reinterpret_cast<char*>(&nrows), sizeof nrows );
        ofh.write( reinterpret_cast<char*>(&ncols), sizeof ncols );
        // matrix
        for( unsigned i = 0; i < nrows; ++i ) {
            for( unsigned j = 0; j < ncols; ++j ) {
                ofh.write( reinterpret_cast<char*>(&mat[i][j]), sizeof mat[i][j]);
            }
        }
        ofh.close();
    } else { // CSV
        std::ofstream ofh;
        ofh.open(ofname);
        if( !ofh.good() ) { return false; }
        ofh << std::setprecision(precision);
        for(unsigned i = 0; i < nrows; i++) {
            for(unsigned j = 0; j < ncols; j++) {
                ofh << mat[i][j];
                if( j != ncols-1 ) { ofh << delim; }
            }
            ofh << std::endl;
        }
        ofh.close();
    }

    return true;
}

bool Matrix::read_file(
        const std::string& ifname,
        bool use_binary,
        char delim
) {
    if( use_binary ) {
        std::ifstream ifh;
        ifh.open(ifname, std::ios::in | std::ios::binary );
        if( !ifh.good() ) { return false; }
        // metadata
        ifh.read( reinterpret_cast<char*>(&nrows), sizeof nrows );
        ifh.read( reinterpret_cast<char*>(&ncols), sizeof ncols );
        init();
        // matrix
        for( unsigned i = 0; i < nrows; ++i ) {
            for( unsigned j = 0; j < ncols; ++j ) {
                ifh.read( reinterpret_cast<char*>(&mat[i][j]), sizeof mat[i][j]);
            }
        }
    } else {
        std::ifstream ifh;
        ifh.open(ifname);
        if( !ifh.good() ) { return false; }
        std::vector<std::vector<double>> tmat;
        std::string temp, token;
        double dtemp;
        unsigned i = 0;
        while( std::getline( ifh, temp ) ) {
            tmat.push_back( std::vector<double>() );
            std::stringstream ss(temp);
            while( std::getline(ss, token, delim) ) {
                std::stringstream conv(token);
                conv >> dtemp;
                tmat[i].push_back(dtemp);
            }
            ++i;
        }

        nrows = tmat.size();
        if( nrows == 0 ) { return false; }
        unsigned tncols = tmat[0].size();
        for( i = 1; i < nrows; ++i ) { // error check
            if( tmat[i].size() != tncols ) {
                return false;
            }
        }
        ncols = tncols;
        init();
        for( unsigned i = 0; i < nrows; ++i ) {
            for( unsigned j = 0; j < ncols; ++j ) {
                mat[i][j] = tmat[i][j];
            }
        }
    }

    return true;
}

