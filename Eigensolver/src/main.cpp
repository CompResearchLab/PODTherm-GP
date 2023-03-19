#include <dolfin.h>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <chrono>
#include <sstream>
#include <fstream>

#include <petscmat.h>
#include <petscsys.h>

using namespace dolfin;
#include "Matrix.hpp"

// Uncomment this line to have all of the parameters be
// printed.
//#define DEBUG

int main(int argc, char** argv) {

    dolfin::init(argc,argv);

    if( argc != 2 ) {
        std::cout << "NEED MATRIX TO EIGENSOLVE" << std::endl;
        exit(1);
    }

    auto comm = MPI_COMM_WORLD;
    int mpi_size,mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    // Here the A Matrix is loaded from the CSV file.
    ::Matrix A_mat;
    bool res = A_mat.read_file(argv[1],false);
    if( !res ) {
        std::cout << "COULDN'T READ FILE!" << std::endl;
        exit(1);
    }
    std::cout << "ROWS: " << A_mat.getNRows() << " COLS: " << A_mat.getNCols() << std::endl;

    // The following block, up to line 59, sets up A PETSc Matrix, and copies the A Matrix into it.
    Mat A_petsc;
    //PetscMalloc2(A_mat.getNRows(),&dnnz,A_mat.getNCols(),&onnz); // allocate 2 arrays of matrix
    MatCreateDense(comm,PETSC_DECIDE,PETSC_DECIDE,A_mat.getNRows(),A_mat.getNCols(),NULL,&A_petsc); // RC is inverse matrix of C matrix
    MatSetUp(A_petsc);
    //PetscFree2(dnnz,onnz);
    MatSetOption(A_petsc,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
    for (unsigned i=0; i<A_mat.getNRows(); i++) {
        for(unsigned j=0; j<A_mat.getNCols();j++){
             MatSetValue(A_petsc,i,j,A_mat[i][j],INSERT_VALUES);
        }
    }
    MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY);

    // A shared pointer to the PETSc A Matrix is created
    // to be passed to the Eigensolver.
    auto A = std::make_shared<PETScMatrix>(A_petsc);
    std::cout << "SOLVING!" << std::endl;

    // Parameters for the Eigensolver are set.
    PETScOptions::set("st_ksp_type", "preonly");// Sets the KSP Solver to only do preconditioning
    PETScOptions::set("st_pc_type", "lu"); // Sets the preconditioner type to only do LU
    PETScOptions::set("-pc_factor_mat_solver_type", "mumps"); // Sets matrix solver to MUMPS
    // More docs: https://slepc.upv.es/documentation/slepc.pdf
    SLEPcEigenSolver esolver(A);
    esolver.parameters["spectrum"] = "largest real"; // sets the eigenvalue spectrum to get the largest real eigenvalues.
    esolver.parameters["problem_type"] = "hermitian"; 
    //esolver.parameters["solver"] = "krylov-schur"; 
    esolver.parameters["tolerance"] = 1e-14; 
    esolver.solve();

    // r, c, rx, and cx are temporary variables for getting
    // the eigenvalues and eigenvectors from the solver.
    double r,c;
    PETScVector rx, cx;
    // r_local is used to gather the eigenvector into when
    // MPI is used.
    std::vector<double> r_local;
    // The following block initializes the eigenmatrix and
    // eigenvalues matrix
    ::Matrix eigenmat, eigenvals;
    if( mpi_rank == 0) {
        // The dimensions of A are used, as we are only
        // interested in that many eigenvalues/vectors
        eigenmat.reset(A_mat.getNRows(), A_mat.getNCols());
        // The eigenvals matrix only needs one column.
        eigenvals.reset(A_mat.getNRows(),1);
    }
    std::cout << std::setprecision(16);
    for( unsigned i = 0; i < eigenmat.getNCols(); i++ ) {
        esolver.get_eigenpair(r,c,rx,cx,i); // get eigenpair from solver
        if(mpi_rank == 0) {std::cout << r << std::endl;}
        rx.gather_on_zero(r_local); // gather MPI values into process 0
        if( mpi_rank == 0 ) {
            eigenvals[i][0] = r; // copy eigenvalue into matrix.
            for( unsigned j = 0; j < r_local.size(); j++ ) {
                eigenmat[i][j] = r_local[j]; // copy eigenvector into matrix
            }
        }
    }

    if( mpi_rank == 0 ) {
        // stores the eigenvalues and eigenvectors as CSV
        // files.
        eigenmat.write_file("../../eigenmat.csv",false);
        eigenvals.write_file("../../eigenvals.csv",false);
    }

    std::cout << "LAST LINE" << std::endl;
    MatDestroy(&A_petsc);

    exit(0);

    //return 0;
}
