#include "POD_ODE_Solver.hpp"

#include <petscmat.h>
#include <stdio.h>
#include <petscts.h>
#include <petscksp.h>

#include<sstream>

namespace Thermal {
namespace POD {
        
/* import C, G, P_vector and temperature gradient to solve ODE.  And then save the coefficient of ppodmode*/

typedef struct {
    Mat A;                 /* RHS mat, used with IFunction interface */
    Vec vp;
} AppCtx;

// include user defined routine
extern PetscErrorCode RHSFunctionHeat(TS ts,PetscReal t,Vec X,Vec r,void *ctx);

bool ODE_Solver::init( int argc, char** argv ) {
    
    PetscInitialize(&argc,&argv,0,0);
    read_xml_config("../../POD_Para.xml");
    bool res = true;
    res &= C.read_file("../../C.csv",false);
    C.print_info("C");
    if(!res) { printf("After loading C!\n"); }
    res &= Gmatrix.read_file("../../G.csv",false);
    Gmatrix.print_info("G");
    if(!res) { printf("After loading G!\n"); }
    printf("FLP\n");
    res = res && load_file(floorplan_file ,MFloorplan,NU,4);
    if(!res) { printf("After loading FLP!\n"); }
    printf("power\n");
    res = res && load_file(ptrace_file,Pmatrix1,NT,NU); // power to be predicted
    if(!res) { printf("After loading power!\n"); }
    res &= P_lib.read_file("../../P_mat.csv",false);
    P_lib.print_info("P");
    if(!res) { printf("After loading P!\n"); }
   /* res &= temp_grad_lib.read_file("TG_mat.csv",false);
    temp_grad_lib.print_info("TG");
    if(!res) { printf("After loading TG!\n"); }
    res &= times_per_step.read_file("timings.csv",false);
    if(!res) { printf("After loading timings!\n"); }*/
    
    return res;
}

int ODE_Solver::run()
{
    // define the petsc parameter
    Mat   CP,GP,PP,PSP,RC,CG,CPSP; // context of eps  
    Vec              u,u0; // eigenvector
    MPI_Comm        comm;
    PetscInt        *dnnz,*onnz,i,j; // variable will be used in slepc
    /*****************************************************************************/
    IS             perm,iperm;
    MatFactorInfo  info;
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    comm = MPI_COMM_WORLD;
    
    /***********************************************************
     * compute the power density. Hereafter, Pmatrix is power density.
     * ***********************************************************/
    printf("PD\n");
    for(int i = 0; i < NT; i++) {
        for(int j =0; j < NU;  j++) {
            Pmatrix1[i][j] = Pmatrix1[i][j]/(MFloorplan[j][0]*MFloorplan[j][1]*thick_actl);             
            //printf("%lf\n",Pmatrix1[i][j]);
        }
    }
    /***********************************************************
     * compute the P vector.
     * ***********************************************************/
    double **PSmatrix = (double**)malloc(NT*sizeof(double*));
    if (!PSmatrix)
        return -1;
    for(int i = 0; i < NT; i++) {
        PSmatrix[i] = (double*)malloc(N_p*sizeof(double));
        if (!PSmatrix[i])
            return -1;
    }

    /**
    ::Matrix new_PS;
    new_PS.read_file("PS_4s.csv",false);
    **/
    double sum_p =0;
    for (int i = 0; i < NT; i++) {
        for (int j = 0; j < N_p; j++) {
            sum_p =0;
            for(int k =0; k < NU; k++) {
                sum_p =  sum_p +P_lib[j][k]*Pmatrix1[i][k];
            }
	    //std::cout << temp_grad_lib.getNRows() << "," << temp_grad_lib.getNCols()
	//	    << ":" << i << "," << j << std::endl;
            //PSmatrix[i][j] = sum_p + temp_grad_lib[j][i]; // BAD
            PSmatrix[i][j] = sum_p; // GOOD
        }
    }
    //check PS
    FILE *fileps = fopen("PS.txt","w");
    for(int i=0; i < NT; i++){
    	for(int j= 0; j<N_p; j++){
		fprintf(fileps,"%.16lg\t",PSmatrix[i][j]);
	}
	fprintf(fileps,"\n");
    }
    fclose(fileps);


    printf("Starting Solver!\n");
    for (int num_pod = 1; num_pod <= NPOD; num_pod++) {
	   // std::cout << "The NPOD is "<<NPOD <<std::endl;
        //std::cout << "STEP: " << num_pod << std::endl;
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,num_pod,num_pod,num_pod,num_pod,NULL,&CP);
        MatSetUp(CP);
        PetscFree2(dnnz,onnz);
        MatSetOption(CP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                MatSetValue(CP,i,j,C[i][j],INSERT_VALUES);
            }
        }
        MatAssemblyBegin(CP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(CP,MAT_FINAL_ASSEMBLY);

        MatGetOrdering(CP,MATORDERINGRCM,&perm,&iperm);
        MatFactorInfoInitialize(&info);
        info.fill          = 1.0;
        MatLUFactor(CP,perm,iperm,&info);

        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,num_pod,num_pod,num_pod,num_pod,NULL,&RC); // RC is inverse matrix of C matrix
        MatSetUp(RC);
        PetscFree2(dnnz,onnz);
        MatSetOption(RC,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                if(i==j)
                    MatSetValue(RC,i,j,1,INSERT_VALUES);
                else
                    MatSetValue(RC,i,j,0,INSERT_VALUES);
            }
        }
        MatAssemblyBegin(RC,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(RC,MAT_FINAL_ASSEMBLY);
        ////////////////////////////////////////////////////////
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,num_pod,num_pod,num_pod,num_pod,NULL,&PP);
        MatSetUp(PP);
        PetscFree2(dnnz,onnz);
        MatSetOption(PP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                if(i==j)
                    MatSetValue(PP,i,j,1,INSERT_VALUES);
                else
                    MatSetValue(PP,i,j,0,INSERT_VALUES);
            }
        }
        MatAssemblyBegin(PP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PP,MAT_FINAL_ASSEMBLY);

        MatMatSolve(CP,PP,RC);

        /*-----------------------------------------------------------------------------------------------------
          save G matrix in petsc format GP
          -------------------------------------------------------------------------------------------------------*/
        PetscMalloc2(num_pod,&dnnz,num_pod,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,num_pod,num_pod,num_pod,num_pod,NULL,&GP);
        MatSetUp(GP);
        PetscFree2(dnnz,onnz);
        MatSetOption(GP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);
        for (i=0; i<num_pod; i++) {
            for(j=0; j<num_pod;j++){
                MatSetValue(GP,i,j,-Gmatrix[i][j],INSERT_VALUES);
            }
        }
        MatAssemblyBegin(GP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(GP,MAT_FINAL_ASSEMBLY);
        //free Gmatrix

        /*------------------------------------------------------------------------------------------------------
          compute C inverse matrix times G matrix
          ---------------------------------------------------------------------------------------------------------*/
        MatMatMult(RC,GP,MAT_INITIAL_MATRIX,1,&CG);

        /*------------------------------------------------------------------------------------------------------
          save source term in petsc format matrix PSP
          ---------------------------------------------------------------------------------------------------------*/

        PetscMalloc2(num_pod,&dnnz,NT,&onnz); // allocate 2 arrays of matrix
        MatCreateDense(comm,num_pod,NT,num_pod,NT,NULL,&PSP);
        MatSetUp(PSP);
        PetscFree2(dnnz,onnz);
        MatSetOption(PSP,MAT_NEW_NONZERO_ALLOCATION_ERR,PETSC_TRUE);

        for (i=0; i<num_pod; i++) {
            for(j=0; j<NT;j++){
                MatSetValue(PSP,i,j,PSmatrix[j][i],INSERT_VALUES); // psp is transpose matrix of PSmatrix
            }
        }
        MatAssemblyBegin(PSP,MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(PSP,MAT_FINAL_ASSEMBLY);

        MatMatMult(RC,PSP,MAT_INITIAL_MATRIX,1,&CPSP);

        /*-----------------------------------------------------------------------------------------------------
          save mat4 matrix in petsc format PHI
          -------------------------------------------------------------------------------------------------------*/
        /*-------------------------------------------------------------------------------------------
          SLOVE THE EQUATION
          -------------------------------------------------------------------------------------------------*/   
        AppCtx         appctx;                 /* user-defined application context */
        TS             ts;                     /* timestepping context */
        // Vec            u,u0;                      /* approximate solution vector */
        PetscInt       time_steps_max = 100;   /* default max timesteps */

        PetscReal      dt;
        VecCreate(PETSC_COMM_WORLD,&u);
        VecSetSizes(u,PETSC_DECIDE,num_pod);
        VecCreate(PETSC_COMM_WORLD,&u0);
        VecSetSizes(u0,PETSC_DECIDE,num_pod);
        VecCreate(PETSC_COMM_WORLD,&appctx.vp);
        VecSetSizes(appctx.vp,PETSC_DECIDE,num_pod);
        VecSetType(appctx.vp, VECMPI);
        VecSetType(u, VECMPI);
        VecSetType(u0, VECMPI);

        VecZeroEntries(u);
        VecZeroEntries(u0);
        VecZeroEntries(appctx.vp);

        MatCreate(PETSC_COMM_SELF,&appctx.A);
        MatSetSizes(appctx.A,PETSC_DECIDE,PETSC_DECIDE,num_pod,num_pod);
        MatSetFromOptions(appctx.A);
        MatSetUp(appctx.A);
        MatDuplicate(CG, MAT_COPY_VALUES,&appctx.A);
	//std::cout <<"The NT"<< NT<<std::endl;

        for (int i = 0; i < NT; i ++){
	  //      sampling_interval = times_per_step[i][0];

            TSCreate(PETSC_COMM_SELF,&ts);
            TSSetProblemType(ts,TS_LINEAR);
            TSSetType(ts,TSRK);
            TSRKSetType( ts,TSRK4);
            TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP );
            TSSetRHSJacobian(ts,CG,CG,TSComputeRHSJacobianConstant,&appctx);

            TSSetMaxSteps(ts,time_steps_max);
            TSSetMaxTime(ts,sampling_interval);

            // set solution vector and initial timestepping
            dt = sampling_interval/20;
            TSSetTimeStep(ts,dt);

            MatGetColumnVector(CPSP,appctx.vp,i);

            //VecCopy(vp,u);
            TSSetRHSFunction(ts,NULL,RHSFunctionHeat,&appctx);
            TSSetSolution(ts,u0);
            TSSolve(ts,u);
            VecZeroEntries(u0);
            VecCopy(u,u0);
	    //std::cout << num_pod <<std::endl;
            for( int k = 0; k < num_pod; k++) {
		    //std::cout << "Num_mode " << num_pod<<"the K " <<k<<std::endl;
                VecGetValues(u,1,&k,&CU[k][i]);
            }
        }
        std::cout << "WRITING CU " << num_pod << std::endl;
        char namecu[200];
        sprintf(namecu, "../../pod_result/CU%d.txt", num_pod);
        FILE *fileMAT4 = fopen(namecu,"w");

        for(int j = 0; j < NT; j++){
            for (int i = 0; i < num_pod; i++){
                fprintf(fileMAT4,"%.16lg\t",CU[i][j]);
            }
            fprintf(fileMAT4,"\n");
        }
        fclose(fileMAT4); 
        TSDestroy(&ts);
        MatDestroy(&CP);
        MatDestroy(&RC);
        MatDestroy(&PP);
        MatDestroy(&GP);
        MatDestroy(&CG);
        MatDestroy(&PSP);
        MatDestroy(&CPSP);
        MatDestroy(&appctx.A);
        VecDestroy(&u);
        VecDestroy(&u0);
        VecDestroy(&appctx.vp);
    }

    // free memory
    // for (int i= 0 ; i < NC; i++){
    //     free(C[i]);
    // }
    // free(C);

    // for (int i= 0 ; i < NC; i++){
    //     free(Gmatrix[i]);
    // }
    // free(Gmatrix);

    for (int i= 0 ; i < NU; i++){
        free(MFloorplan[i]);
    }
    free(MFloorplan);

    for (int i= 0 ; i < NT; i++){
        free(Pmatrix1[i]);
    }
    free(Pmatrix1);

    // for (int i= 0 ; i < N_p; i++){
    //     free(P_lib[i]);
    // }
    // free(P_lib);

    // for (int i= 0 ; i < NT; i++){
    //     free(temp_grad_lib[i]);
    // }
    // free(temp_grad_lib);

    for (int i= 0 ; i < NT; i++){
        free(PSmatrix[i]);
    }
    free(PSmatrix);

    return 0;
}

PetscErrorCode  RHSFunctionHeat(TS ts,PetscReal t,Vec X,Vec r,void *ctx)
{
    AppCtx *appctx = (AppCtx*)ctx;     /* user-defined application context */

    MatMult(appctx -> A,X,r);
    VecAXPY(r,1.0,appctx -> vp);

    return 0;
}

bool ODE_Solver::load_file(
    const std::string& ifname,
    double**& mat,
    const unsigned& nrows,
    const unsigned& ncols
) {
    mat = (double**)malloc(nrows*sizeof(double*));
    if(!mat) { return false; }
    for(int i = 0; i < nrows; i++){
        mat[i] = (double*)malloc(ncols*sizeof(double));
        if (!mat[i]) { return false; }
    }

    FILE* ifh = fopen(ifname.c_str(),"r");
    for(int i = 0; i < nrows; i++){
        for(int j = 0; j < ncols; j++){
            mat[i][j]=0;
            if(!fscanf(ifh , "%lf", &mat[i][j])) { break; }
            //printf("%lf\n",mat[i][j]);
        }
    }
    fclose(ifh);
    return true;
}


} // ns POD
} // ns Thermal
