#ifndef POD_PAR_H
#define POD_PAR_H
unsigned N_STATE = 0; 		// 0: Perform thermal simulation using FEM; 1: Read the existing FEM result.
double l = 0.014;		// The Length
double w = 0.012;		// Width
double h = 0.0002977976;		// Thickness in total.
unsigned ls = 200;		//Number of grids in length
unsigned ws = 200;		//Number of grids in width
unsigned hs = 17;		//Number of grids in thickness
unsigned num_steps = 240;	//Number of time step
double T = 240*3.125e-6;	//total simulation time
double delt = T/num_steps;	//the time step size
double t = 0;			//the initial time step
double h_c = 2.40598e4;		//The heat transfer coefficient applied on the bottome of substrate
double Ta = 0.0;		//The initial temperature value
unsigned Nu = 11; 		//The number of functional units of processor floorplan
double thick_actl = 0.0000557976;	//The thickness of device layer
double thick_Sio2 = thick_actl;	//The thickness of the Sio2 layer
double chip_area = l*w;		//The area of processor
bool status;
double tol = 1E-14;		//The tolerance used in the simulation
double k_0 = 100.0, k_1 = 1.2;	//The thermal conductivity k_0 (Silicon) k_1 (Sio2)
double D0 = 2.33e3, D1 = 2.65e3;//The density; Do(Silicon) D1(Sio2)
double c0 = 751.1, c1 = 680;	//Specific heat; c0 (Silicon) c1(Sio2)
#endif

