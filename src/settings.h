#ifndef SETTINGS_H
#define SETTINGS_H

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
/*    Computation Settings        */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
#define LX 1                // w (width)
#define LY 2.5              // h (height)
#define FACTOR 320          // Cell numbers per unit of length
#define NI int(LX * FACTOR) // Mesh size in x-direction
#define NJ int(LY * FACTOR) // Mesh size in y-direction
#define MAXSIZE (((NI > NJ) ? NI : NJ) + 2)
#define dx (double(LX) / NI) // Delta x
#define dy (double(LY) / NJ) // Delta y
#define TOL 1e-9             // Tolerance of L2 norm of change in solution
#define MAX_ITERS 50000      // Max number of iterations allowed
#define IMPLICIT_DT 0.11     // dt for implicit scheme
#define OMEGA 1.25           // Overrelaxation parameter

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
/*          Case Setup            */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */

#define RE 300.   // Reynolds number
#define beta 0.75 // Artificial compressibility parameter

/* Initial Condition Parameters */
#define P0 1.
#define u0 1.
#define v0 1.

/* Boundary Condotion Paramters */
#define U_TOP 1.

#endif