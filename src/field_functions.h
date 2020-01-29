#ifndef FIELD_FUNCTIONS_H
#define FIELD_FUNCTIONS_H

#include <math.h>
#include "settings.h"

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
/*      Field Functions           */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */

// Initial fields:
double initial_P(double x, double y);

double initial_u(double x, double y);

double initial_v(double x, double y);

// Exact flux integral fields for problem 1.1:
double exact_FI_P(double x, double y);

double exact_FI_u(double x, double y);

double exact_FI_v(double x, double y);

#endif