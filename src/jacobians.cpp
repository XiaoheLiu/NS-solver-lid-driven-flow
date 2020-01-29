#include "jacobians.h"

void set_Ax(const VectorField &U, int i, int j, double (&J)[3][3])
{
  J[0][0] = 0;
  J[0][1] = (-1. / dx) / (2 * beta);
  J[0][2] = 0;
  J[1][0] = (-1. / dx) / 2;
  J[1][1] = (-1. / dx) * ((U.u(i, j) + U.u(i - 1, j)) / 2 + 1. / (dx * RE));
  J[1][2] = 0;
  J[2][0] = 0;
  J[2][1] = (-1. / dx) * ((U.v(i, j) + U.v(i - 1, j)) / 4);
  J[2][2] = (-1. / dx) * ((U.u(i, j) + U.u(i - 1, j)) / 4 + 1. / (dx * RE));
}

void set_Cx(const VectorField &U, int i, int j, double (&J)[3][3])
{
  J[0][0] = 0;
  J[0][1] = (1. / dx) / (2 * beta);
  J[0][2] = 0;
  J[1][0] = (1. / dx) / 2;
  J[1][1] = (1. / dx) * ((U.u(i, j) + U.u(i + 1, j)) / 2 - 1 / (dx * RE));
  J[1][2] = 0;
  J[2][0] = 0;
  J[2][1] = (1. / dx) * ((U.v(i, j) + U.v(i + 1, j)) / 4);
  J[2][2] = (1. / dx) * ((U.u(i, j) + U.u(i + 1, j)) / 4 - 1 / (dx * RE));
}

void set_Bx(const VectorField &U, int i, int j, double (&J)[3][3])
{
  double Ax[3][3], Cx[3][3];
  set_Ax(U, i + 1, j, Ax);
  set_Cx(U, i - 1, j, Cx);
  for (int m = 0; m < 3; m++)
    for (int n = 0; n < 3; n++)
      J[m][n] = -Ax[m][n] - Cx[m][n];
}

void set_Ay(const VectorField &U, int i, int j, double (&J)[3][3])
{
  J[0][0] = 0;
  J[0][1] = 0;
  J[0][2] = (-1. / dy) / (2 * beta);
  J[1][0] = 0;
  J[1][1] = (-1. / dy) * ((U.v(i, j) + U.v(i, j - 1)) / 4 + 1 / (dy * RE));
  J[1][2] = (-1. / dy) * ((U.u(i, j) + U.u(i, j - 1)) / 4);
  J[2][0] = (-1. / dy) / 2;
  J[2][1] = 0;
  J[2][2] = (-1. / dy) * ((U.v(i, j) + U.v(i, j - 1)) / 2 + 1 / (dy * RE));
}

void set_Cy(const VectorField &U, int i, int j, double (&J)[3][3])
{
  J[0][0] = 0;
  J[0][1] = 0;
  J[0][2] = (1. / dy) / (2 * beta);
  J[1][0] = 0;
  J[1][1] = (1. / dy) * ((U.v(i, j) + U.v(i, j + 1)) / 4 - 1 / (dy * RE));
  J[1][2] = (1. / dy) * ((U.u(i, j) + U.u(i, j + 1)) / 4);
  J[2][0] = (1. / dy) / 2;
  J[2][1] = 0;
  J[2][2] = (1. / dy) * ((U.v(i, j) + U.v(i, j + 1)) / 2 - 1 / (dy * RE));
}

void set_By(const VectorField &U, int i, int j, double (&J)[3][3])
{
  double Ay[3][3], Cy[3][3];
  set_Ay(U, i, j + 1, Ay);
  set_Cy(U, i, j - 1, Cy);
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      J[m][n] = -Ay[m][n] - Cy[m][n];
    }
  }
}

void set_x_jacobians(const VectorField &U, int i, int j, double (&Ax)[3][3], double (&Bx)[3][3], double (&Cx)[3][3])
{
  set_Ax(U, i, j, Ax);
  set_Bx(U, i, j, Bx);
  set_Cx(U, i, j, Cx);
}

void set_y_jacobians(const VectorField &U, int i, int j, double (&Ay)[3][3], double (&By)[3][3], double (&Cy)[3][3])
{
  set_Ay(U, i, j, Ay);
  set_By(U, i, j, By);
  set_Cy(U, i, j, Cy);
}

void mat_dot_vec(const double (&M)[3][3], const double (&V)[3], double(Result)[3])
{
  Result[0] = M[0][0] * V[0] + M[0][1] * V[1] + M[0][2] * V[2];
  Result[1] = M[1][0] * V[0] + M[1][1] * V[1] + M[1][2] * V[2];
  Result[2] = M[2][0] * V[0] + M[2][1] * V[1] + M[2][2] * V[2];
}

void sum_six_vectors(const double (&vecs)[6][3], double (&result)[3])
{
  for (int k = 0; k < 3; k++)
  {
    result[k] = 0.;
    for (int i = 0; i < 6; i++)
    {
      result[k] += vecs[i][k];
    }
  }
}

void set_diag_mat(double (&M)[3][3], double a, double b, double c)
{
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      M[i][j] = 0.;
    }
  }
  M[0][0] = a;
  M[1][1] = b;
  M[2][2] = c;
}

void set_zero_vector(double (&V)[3])
{
  V[0] = 0.;
  V[1] = 0.;
  V[2] = 0.;
}

double max_element_arr(double (&a)[3])
{
  double max = a[0];
  for (int i = 1; i < 3; i++)
  {
    if (a[i] > max)
      max = a[i];
  }
  return max;
}