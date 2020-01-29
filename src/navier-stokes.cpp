#include "navier-stokes.h"

void initialize(VectorField &U)
{
  double (*initial_U[])(double, double) = {initial_P, initial_u, initial_v};
  U.set_fields(initial_U);
  // U.write_dat("initial_U"); // verified.
  update_ghost_cells(U);
}

void calculate_flux_integral(const VectorField &U, VectorField &FI)
{
  for (int i = 1; i < U.Ni + 1; i++)
  {
    for (int j = 1; j < U.Nj + 1; j++)
    {
      double fi[3];
      // F.I for P:
      fi[0] = -(((U.u(i + 1, j) - U.u(i - 1, j)) / (2. * beta)) / dx + ((U.v(i, j + 1) - U.v(i, j - 1)) / (2. * beta)) / dy);

      // F.I. for u:
      fi[1] = -((pow(U.u(i + 1, j) + U.u(i, j), 2) - pow(U.u(i, j) + U.u(i - 1, j), 2)) / 4 + (U.P(i + 1, j) - U.P(i - 1, j)) / 2 - (U.u(i + 1, j) - 2 * U.u(i, j) + U.u(i - 1, j)) / (RE * dx)) / dx;

      fi[1] -= ((((U.u(i, j + 1) + U.u(i, j)) * (U.v(i, j + 1) + U.v(i, j))) - (U.u(i, j) + U.u(i, j - 1)) * (U.v(i, j) + U.v(i, j - 1))) / 4 - (U.u(i, j + 1) - 2 * U.u(i, j) + U.u(i, j - 1)) / (RE * dy)) / dy;

      // F.I. for v:
      fi[2] = -((((U.u(i + 1, j) + U.u(i, j)) * (U.v(i + 1, j) + U.v(i, j))) - (U.u(i, j) + U.u(i - 1, j)) * (U.v(i, j) + U.v(i - 1, j))) / 4 - (U.v(i + 1, j) - 2 * U.v(i, j) + U.v(i - 1, j)) / (RE * dx)) / dx;

      fi[2] -= ((pow(U.v(i, j + 1) + U.v(i, j), 2) - pow(U.v(i, j) + U.v(i, j - 1), 2)) / 4 + (U.P(i, j + 1) - U.P(i, j - 1)) / 2 - (U.v(i, j + 1) - 2 * U.v(i, j) + U.v(i, j - 1)) / (RE * dy)) / dy;

      FI.set(i, j, fi);
    }
  }
}

void assemble_jacobians(const VectorField &U, const VectorField &dU, VectorField &RHS)
{
  for (int i = 1; i < dU.Ni + 1; i++)
  {
    for (int j = 1; j < dU.Nj + 1; j++)
    {
      double Ax[3][3], Bx[3][3], Cx[3][3];
      double Ay[3][3], By[3][3], Cy[3][3];
      set_x_jacobians(U, i, j, Ax, Bx, Cx);
      set_y_jacobians(U, i, j, Ay, By, Cy);

      double dUim[3]{dU.P(i - 1, j), dU.u(i - 1, j), dU.v(i - 1, j)},
          dUij[3]{dU.P(i, j), dU.u(i, j), dU.v(i, j)},
          dUip[3]{dU.P(i + 1, j), dU.u(i + 1, j), dU.v(i + 1, j)},
          dUjm[3]{dU.P(i, j - 1), dU.u(i, j - 1), dU.v(i, j - 1)},
          dUjp[3]{dU.P(i, j + 1), dU.u(i, j + 1), dU.v(i, j + 1)};

      double v[6][3]; // six vectors to be summed
      mat_dot_vec(Ax, dUim, v[0]);
      mat_dot_vec(Bx, dUij, v[1]);
      mat_dot_vec(Cx, dUip, v[2]);
      mat_dot_vec(Ay, dUjm, v[3]);
      mat_dot_vec(By, dUij, v[4]);
      mat_dot_vec(Cy, dUjp, v[5]);

      double new_values[3];
      sum_six_vectors(v, new_values);

      RHS.set(i, j, new_values);
    }
  }
}

void update_ghost_cells(VectorField &U)
{
  // Dirichlet conditions for u and v, Neumann for P
  double U0 = 0., U_top = U_TOP, bc[3];
  for (int j = 0; j < U.Nj + 2; j++)
  {
    // Left: dP=0, u=v=0
    bc[0] = U.P(1, j);
    bc[1] = 2. * U0 - U.u(1, j);
    bc[2] = 2. * U0 - U.v(1, j);
    U.set(0, j, bc);

    // Right: dP=0, u=v=0
    bc[0] = U.P(U.Ni, j);
    bc[1] = 2. * U0 - U.u(U.Ni, j);
    bc[2] = 2. * U0 - U.v(U.Ni, j);
    U.set(U.Ni + 1, j, bc);
  }

  for (int i = 0; i < U.Ni + 2; i++)
  {
    // Top: dP=0, u=U_top, v=0
    bc[0] = U.P(i, U.Nj);
    bc[1] = 2. * U_top - U.u(i, U.Nj);
    bc[2] = 2. * U0 - U.v(i, U.Nj);
    U.set(i, U.Nj + 1, bc);

    // Bottom: dP=0, u=v=0
    bc[0] = U.P(i, 1);
    bc[1] = 2. * U0 - U.u(i, 1);
    bc[2] = 2. * U0 - U.v(i, 1);
    U.set(i, 0, bc);
  }
}

void approximate_factorization(VectorField &U, const double dt, double (&L2norms)[3])
{
  /* (1) First pass (solving for lines of constant j):
           ([I] + dt * [Dx]) dU_star = dt * [R]
  */
  double LHS[MAXSIZE][3][3][3], RHS[MAXSIZE][3];
  VectorField FI(LX, LY, NI, NJ);
  calculate_flux_integral(U, FI);
  FI *= dt;
  // FI.write_dat("before_1st_pass"); // verified.

  for (int j = 1; j < U.Nj + 1; j++)
  {
    /* (1.1) Set up interior matrices for LHS
                LHS[i][0] = dt * [Ax;i,j]
                LHS[i][1] = [I] + dt * [Bx;i,j]
                LHS[i][2] = dt * [Cx;i,j]
    */
    for (int i = 1; i < U.Ni + 1; i++)
    {
      double Ax[3][3], Bx[3][3], Cx[3][3];
      set_x_jacobians(U, i, j, Ax, Bx, Cx);
      for (int m = 0; m < 3; m++)
      {
        for (int n = 0; n < 3; n++)
        {
          LHS[i][0][m][n] = dt * Ax[m][n];
          LHS[i][1][m][n] = dt * Bx[m][n] + ((m == n ? 1. : 0.));
          LHS[i][2][m][n] = dt * Cx[m][n];
        }
      }

      // (1.2) Set up RHS = dt * [R]
      RHS[i][0] = FI.P(i, j);
      RHS[i][1] = FI.u(i, j);
      RHS[i][2] = FI.v(i, j);
    }
    // (1.3) Set up implicit boundary conditions
    // i = 0, Left
    set_diag_mat(LHS[0][1], 1., 1., 1.);
    set_diag_mat(LHS[0][2], -1., 1., 1.);
    set_zero_vector(RHS[0]);
    // i = NI + 1, Right
    set_diag_mat(LHS[U.Ni + 1][0], 1., 1., 1.);
    set_diag_mat(LHS[U.Ni + 1][1], -1., 1., 1.);
    set_zero_vector(RHS[U.Ni + 1]);

    // (1.4) Solve the tri-diagonal system
    SolveBlockTri(LHS, RHS, U.Ni + 2);

    // (1.5) Save the intermediate solution dU_star(:,j) = RHS
    for (int i = 0; i < U.Ni + 2; i++)
    {
      FI.set(i, j, RHS[i]); // Overwrite FI with the intermediate solution dU_star
    }
  }
  // FI.write_dat("intermediate_dU");

  /* (2) Second pass (solving for lines of constant i):
           ([I] + dt * [Dy]) dU = dU_star
  */

  for (int i = 1; i < U.Ni + 1; i++)
  {
    /* (2.1) Set up interior matrices for LHS
                LHS[j][0] = dt * [Ay;i,j]
                LHS[j][1] = [I] + dt * [By;i,j]
                LHS[j][2] = dt * [Cy;i,j]
    */
    for (int j = 1; j < U.Nj + 1; j++)
    {
      double Ay[3][3], By[3][3], Cy[3][3];
      set_y_jacobians(U, i, j, Ay, By, Cy);
      for (int m = 0; m < 3; m++)
      {
        for (int n = 0; n < 3; n++)
        {
          LHS[j][0][m][n] = dt * Ay[m][n];
          LHS[j][1][m][n] = dt * By[m][n] + ((m == n ? 1. : 0.));
          LHS[j][2][m][n] = dt * Cy[m][n];
        }
      }

      // (2.2) Set up RHS = dU_star
      RHS[j][0] = FI.P(i, j);
      RHS[j][1] = FI.u(i, j);
      RHS[j][2] = FI.v(i, j);
    }
    // (2.3) Set up implicit boundary conditions
    // j = 0, Bottom
    set_diag_mat(LHS[0][1], 1., 1., 1.);
    set_diag_mat(LHS[0][2], -1., 1., 1.);
    set_zero_vector(RHS[0]);
    // j = NJ + 1, Top
    set_diag_mat(LHS[U.Nj + 1][0], 1., 1., 1.);
    set_diag_mat(LHS[U.Nj + 1][1], -1., 1., 1.);
    set_zero_vector(RHS[U.Nj + 1]);

    // (2.4) Solve the tri-diagonal system
    SolveBlockTri(LHS, RHS, U.Nj + 2);

    // (2.5) Save the final solution dU(i,:) = RHS
    for (int j = 0; j < U.Nj + 2; j++)
    {
      FI.set(i, j, RHS[j]); // Overwrite FI with the solution dU
    }
  }
  // FI.write_dat("final_dU");

  // (3) Update interior: T = T + dT
  FI *= OMEGA; // Over relaxation
  U += FI;
  // U.write_dat("final_U");

  L2norms[0] = FI.calculate_L2_norm(0);
  L2norms[1] = FI.calculate_L2_norm(1);
  L2norms[2] = FI.calculate_L2_norm(2);
}

void run_to_convergence(VectorField &U, const string log_name)
{
  double L2norms[3];
  int iter = 1;
  ofstream log;
  log.open("./outputs/log-" + log_name + ".csv");
  log << "Iter,dP,du,dv\n";
  cout << "Start iterating...\n";
  double max_change = 1.;

  while (iter <= MAX_ITERS && max_change > TOL)
  {
    approximate_factorization(U, IMPLICIT_DT, L2norms);
    log << iter << "," << L2norms[0] << "," << L2norms[1] << "," << L2norms[2] << "\n";
    if (iter % 10 == 0)
    {
      cout << "[" << iter << "] \t dP: " << L2norms[0] << "\t, |du|: " << L2norms[1] << "\t, |dv|: " << L2norms[2] << "\n";
    }
    update_ghost_cells(U);
    iter += 1;
    max_change = max_element_arr(L2norms);
  }
}