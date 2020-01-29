#include "settings.h"        // Settings for the computation domain and solver
#include "VectorField.h"     // Custom class for a vector field
#include "field_functions.h" // Functions that define the initial and boundary fields
#include "navier-stokes.h"   // Functions for the 2D incompressible Navier-Stokes solver

/* Validation Suites */
void test_residual();  // 1.1 test correctness of residual
void test_jacobians(); // 1.2 test correctness implicit discretization - flux jacobian

int main()
{
  // test_residual();
  // test_jacobians();

  VectorField U(LX, LY, NI, NJ);

  initialize(U);

  string fileName = "4_mesh" + to_string(FACTOR);

  run_to_convergence(U, fileName);

  U.write_csv(fileName);

  U.output_vorticity(fileName);
  // U.output_midline_u();

  return 0;
}

void test_residual()
{
  cout << "Test correctness of residual...\n";
  double (*initial_U[])(double, double) = {initial_P, initial_u, initial_v};
  VectorField U(LX, LY, NI, NJ, initial_U);

  VectorField FI(LX, LY, NI, NJ);
  calculate_flux_integral(U, FI);

  double (*exact_FI[])(double, double) = {exact_FI_P, exact_FI_u, exact_FI_v};
  VectorField FI_e(LX, LY, NI, NJ, exact_FI);

  FI.calculate_L2_error(FI_e);
}

void test_jacobians()
{
  cout << "Test correctness of flux Jacobian...\n";
  double (*initial_U[])(double, double) = {initial_P, initial_u, initial_v};
  VectorField U(LX, LY, NI, NJ, initial_U);

  VectorField FI_n(LX, LY, NI, NJ);
  calculate_flux_integral(U, FI_n);

  VectorField dU(LX, LY, NI, NJ);
  double delta[3]{1e-6, 1e-6, 1e-6};
  // double delta[3]{0, 1e-6, 0};
  // double delta[3]{0, 0, 1e-6};
  dU.set(10, 10, delta);

  U += dU;

  VectorField FI_p(LX, LY, NI, NJ);
  calculate_flux_integral(U, FI_p);
  FI_p -= FI_n; // Now FI_p is (-LHS) of Eqn. 1

  VectorField RHS(LX, LY, NI, NJ);
  assemble_jacobians(U, dU, RHS); // Now RHS is the RHS of Eqn. 1

  RHS += FI_p; // Calculate the error
  RHS.print_partial(0);
  RHS.print_partial(1);
  RHS.print_partial(2);
}