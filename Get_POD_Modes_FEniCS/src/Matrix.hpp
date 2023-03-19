#ifndef _MATH_MATRIX_HPP_
#define _MATH_MATRIX_HPP_

#include<vector>
#include<iostream>

class Matrix {
public:
    Matrix(const unsigned& nrows_in, const unsigned& ncols_in) {
        nrows = nrows_in;
        ncols = ncols_in;
        init();
    }

    Matrix() :
        Matrix(0,0)
    {}

    Matrix(const Matrix& other) {
        nrows = other.getNRows();
        ncols = other.getNCols();
        init();
        for( unsigned i = 0; i < nrows; i++ ) {
            for( unsigned j = 0; j < ncols; j++ ) {
                mat[i][j] = other.mat[i][j];
            }
        }
    }

    Matrix times(Matrix& other);

    void print() const;

    unsigned getNRows() const { return nrows; }
    unsigned getNCols() const { return ncols; }

    bool write_file( 
        const std::string& ofname, 
        bool use_binary=true,
        char delim=',',
        unsigned precision=16
    );

    bool read_file( 
        const std::string& ifname, 
        bool use_binary=true,
        char delim=','
    );

    std::vector<double>& operator[](const unsigned& i) {
        return mat[i];
    }

    void reset( const unsigned& nnrow, const unsigned& nncol ) {
        nrows = nnrow;
        ncols = nncol;
        init();
    }
    double getMax();

    void add_row() {
        mat.push_back(std::vector<double>());
        for( unsigned i = 0; i < ncols; ++i ) {
            mat[nrows].push_back(0.0);
        }
        ++nrows;
    }

    void remove_row() {
        if( nrows == 0 ) { return; }
        mat.pop_back();
    }

    std::vector<std::vector<double>> mat;

private:
    void init();

    unsigned nrows;
    unsigned ncols;

};

#endif
