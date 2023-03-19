#!/bin/bash
echo "Predicting results from FEM!"
cd ./Data_Collection_FEM/build
#mpirun -n 20 ./Therm_FEM
cd ../..
mkdir pod_result
echo "Solve ODE Equations!"
cd ./ODE_Solver_CU/build
./Calculate_CU
cd ../..
echo "Post_processing!"
python3 Post_Processing.py
