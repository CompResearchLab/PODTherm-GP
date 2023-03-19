#!/bin/bash
echo "Data Collection!"
cd ./Data_Collection_FEM/build
mpirun -n 20 ./Therm_FEM 
cd ../..
echo "Computing Correlation Matrix!"
cd ./Calculate_A/build
mpirun -n 20 ./calculate_A_Matrix
cd ../../
echo "Solve Eigen Pairs!"
cd ./Eigensolver/build
./Eigensolver ../../A.csv 
cd ../..
echo "Generate POD Modes!"
cd ./Get_POD_Modes_FEniCS/build
./get_POD_Modes
cd ../..
echo "Compute Conductance Matrix!"
cd ./Calculate_G/build/
./calculate_G_Matrix
cd ../..
echo "Compute Capacitance Matrix!"
cd ./Calculate_C/build
./calculate_C_Matrix
cd ../..
echo "Compute Mode Integration"
cd ./Calculate_P/build
./calculate_P_Matrix
cd ../..
