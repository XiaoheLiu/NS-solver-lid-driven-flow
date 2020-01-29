#ifndef NAVIER_STOKES_H
#define NAVIER_STOKES_H

#include "VectorField.h"
#include "field_functions.h"
#include "settings.h"
#include "jacobians.h"  // Flux jacobians, and helper functions for linear algebra
#include "tri_thomas.h" // Block Thomas algorithm

void initialize(VectorField &U);
/* Function: set initial condition, and set ghost cell values 
   Parameters: 
    - "U": a blank VectorField to be initialized
   Returns: Null
*/

void calculate_flux_integral(const VectorField &U, VectorField &FI);
/* Function: calculate the flux integral using 2nd-order centered scheme
   Parameters: 
    - "U": the current solution field
    - "FI": the flux integral field to be set through this function
   Returns: Null
*/

void assemble_jacobians(const VectorField &U, const VectorField &dU, VectorField &RHS);
/* Function: assemble the flux jacobians according to the RHS of Eq. (1) in section 1.2
   Parameters: 
    - "U": U^n, current solution
    - "dU": delta U, change in solution
    - "RHS": the VectorField to store the result of this function
   Returns: Null
*/

void update_ghost_cells(VectorField &U);
/* Function: update ghost cells according to section 3.1
   Parameters: 
    - "U": U^n, current solution to be updated
   Returns: Null
*/

void approximate_factorization(VectorField &U, const double dt, double (&L2norms)[3]);
/* Function: perform one iteration of Implicit Euler time advance using approximate factorization.
   Parameters: 
    - "U": the current field U^n. It will be updated to U^{n+1} after the function ends.
    - "dt": time step size
    - "L2norms": the L2 norms of the change in P, u, v will be populated into this size 3 array
   Returns: Null
*/

void run_to_convergence(VectorField &U, const string log_name);
/* Function: iterate using the Implicit Euler time scheme until the solution reaches convergence or maximum iterations allowed.
   Parameters: 
    - "U": the initialized field U^n. It will be updated to the converged solution when the function ends.
    - "log_name": file name to save the convergence log
   Returns: Null
*/

#endif