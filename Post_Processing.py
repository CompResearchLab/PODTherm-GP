# fenics code for the thermal simulation of chip for the publication  for frequency is 3.5 GHZ
import sys
import meshio
from fenics import *
from dolfin import * 
from mshr import *
from numpy import loadtxt
from petsc4py import PETSc
import numpy as np
import time
import xml.etree.ElementTree as ET
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
tree = ET.parse('POD_Para.xml')
root = tree.getroot()
variables = {}
for var in root.findall('variable'):
     name = var.get('name')
     value = var.text
     variables[name] = value

num_steps = int(variables['num_steps_pre_in'])                                              # num_steps is the number of time steps
t = float(variables['t_in'])
Train_steps = num_steps
Ta = float(variables['Ta_in'])    
N_mode =int(variables['num_modes_pre_in'])
l = float(variables['l_in'])
w = float(variables['w_in'])
h = float(variables['h_in'])
ls = int(variables['ls_in']) 
ws = int(variables['ws_in'])
hs = int(variables['hs_in'])
mesh = BoxMesh(Point(0,0,0), Point(l,w,h),ls-1,ws-1,hs-1)
coor1 = mesh.coordinates()
lr = coor1[:,0].max()
ll = coor1[:,0].min()
wb = coor1[:,1].min()
wt = coor1[:,1].max()
hmax = coor1[:,2].max()
hmin = coor1[:,2].min()
V = FunctionSpace(mesh, 'P', 1)

Num_nodes = mesh.num_vertices()
thick_actl =float(variables['thick_actl_in'])
thick_Sio2 = float(variables['thick_Sio2_in'])
#number_mode = loadtxt('config_block.txt')
#define thermal conductivity 
tol = 1E-14
#define initial value 
u0 = Constant(Ta)                                                                   # Ta is initial temperature 
u_n = interpolate(u0,V)
# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
solution = []
for n in range(0,Train_steps):
	solution.append('u1')
# Collect Neumann integrals
#a = DS1*sc*u*v*dx + dt*kappa*dot(grad(u), grad(v))*dx + sum(integrals_R_a)
coor = mesh.coordinates()
v2d = vertex_to_dof_map(V)
h = coor[:,2].max()
#CU = loadtxt('pod_result/CU.txt')
#print(len(CU[0]))
# if we have the solutiuon, we just need to  Load solution# #######################read the solution from a solution file ######################
u = Function(V)
for n in range(0,Train_steps):
    solution_load_file_name = "./solution/file_" + str(n) + "h5"
    solution_file = HDF5File(mesh.mpi_comm(), solution_load_file_name, "r")
    solution_file.read(u, "solution")
    solution[n] = interpolate(u0,V)
    solution[n].assign(u)
    solution_file.close()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~generate podmode~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ###########################################################################
podmode = []
for n in range(0,N_mode):
    podmode.append('u1')
for n in range(0,N_mode):
    podmode_load_file_name = "./POD_mode/mode_" + str(n) + "h5"
    podmode_file = HDF5File(mesh.mpi_comm(), podmode_load_file_name, "r")
    podmode_file.read(u, "solution")
    podmode[n] = interpolate(u0,V)
    podmode[n].assign(u)
    podmode_file.close()

E_T = 0
for i in range (0,Train_steps):
    E_T += assemble(dot(solution[i],solution[i])*dx)
E_diff = 0

for n_t in range(0, N_mode):
    print(n_t)
    CU_filename = 'pod_result/CU'+str(n_t +1)+".txt"
    CU = loadtxt(CU_filename)
    E_diff = 0
    for i in range (0,Train_steps):
        a_init = Constant(0)
        a = interpolate(a_init,V)
        for j in range (0, n_t +1):
            if n_t ==0:
                a.vector().axpy(CU[i], podmode[j].vector())
            else:
                a.vector().axpy(CU[i][j], podmode[j].vector())
        E_diff += assemble(dot((a - solution[i]),(a - solution[i]))*dx)
    LS_name = "pod_result/LS_error.txt"
    LS_file = open(LS_name,"a")
    LS_file.write('%d\t' % (n_t))
    LS_file.write('%.16g\n' % (100*sqrt(E_diff/E_T)))
    LS_file.close()
    tol = 1e-14 # avoid hitting points outside the domain
    y = np.linspace(wb+ tol, wt- tol, 100)
    points = [(0.00361, y_,0.000242) for y_ in y]  # 2D points
    w_line = np.array([a(point) for point in points])
    p_line = np.array([solution[239](point) for point in points])
    if n_t ==0:
        header_data = 'pod_result/y_block'
        data_file_name = header_data + '.txt'
        data_file = open(data_file_name,'w')
        for ii in range(0,len(y)):
            data_file.write('%.16g\n' % (y[ii]))
        data_file.close()
    header_data = 'pod_result/T_Y_'+str(n_t + 1)+'_mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for ii in range(0,len(y)):
        data_file.write('%.16g\t' % (w_line[ii]))
        data_file.write('%.16g\n' % (p_line[ii]))
    data_file.close()
    y = np.linspace(ll + tol, lr- tol, 100)
    points = [(y_,0.009205,0.000242) for y_ in y]
    w_line = np.array([a(point) for point in points])
    p_line = np.array([solution[239](point) for point in points])
    if n_t == 0:
        header_data = 'pod_result/x_block'
        data_file_name = header_data + '.txt'
        data_file = open(data_file_name,'w')
        for ii in range(0,len(y)):
            data_file.write('%.16g\n' % (y[ii]))
        data_file.close()
    header_data = 'pod_result/T_X_'+str(n_t +1)+'_mode'
    data_file_name = header_data + '.txt'
    data_file = open(data_file_name,'w')
    for ii in range(0,len(y)):
        data_file.write('%.16g\t' % (w_line[ii]))
        data_file.write('%.16g\n' % (p_line[ii]))
    data_file.close()

