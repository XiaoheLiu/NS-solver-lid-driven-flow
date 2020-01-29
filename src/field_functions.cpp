#include "field_functions.h"

#define cpx cos(M_PI *x)
#define cpy cos(M_PI *y)
#define spx sin(M_PI *x)
#define spy sin(M_PI *y)
#define c2px cos(2 * M_PI * x)
#define c2py cos(2 * M_PI * y)
#define s2px sin(2 * M_PI * x)
#define s2py sin(2 * M_PI * y)
#define pi M_PI

double initial_P(double x, double y)
{
  return P0 * cpx * cpy;
}

double initial_u(double x, double y)
{
  return u0 * spx * s2py;
}

double initial_v(double x, double y)
{
  return v0 * s2px * spy;
}

double exact_FI_P(double x, double y)
{
  return -pi / beta * (u0 * cpx * s2py + v0 * s2px * cpy);
}

double exact_FI_u(double x, double y)
{
  double u = P0 * pi * spx * cpy - pow(u0, 2) * pi * s2px * pow(s2py, 2);
  u -= u0 * v0 * pi * spx * s2px * (cpy * s2py + 2. * c2py * spy);
  u -= u0 * 5 * pi * pi * spx * s2py / RE;
  return u;
}

double exact_FI_v(double x, double y)
{
  double v = P0 * pi * cpx * spy - pow(v0, 2) * pi * s2py * pow(s2px, 2);
  v -= u0 * v0 * pi * spy * s2py * (cpx * s2px + 2. * c2px * spx);
  v -= v0 * 5 * pi * pi * spy * s2px / RE;
  return v;
}