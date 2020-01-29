#ifndef JACOBIANS_H
#define JACOBIANS_H

#include "VectorField.h"
#include "settings.h"

/* Flux Jacobians */

void set_Ax(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_Cx(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_Bx(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_Ay(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_Cy(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_By(const VectorField &U, int i, int j, double (&J)[3][3]);

void set_x_jacobians(const VectorField &U, int i, int j, double (&Ax)[3][3], double (&Bx)[3][3], double (&Cx)[3][3]);

void set_y_jacobians(const VectorField &U, int i, int j, double (&Ay)[3][3], double (&By)[3][3], double (&Cy)[3][3]);

/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */
/*         Helper Functions       */
/* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */

void sum_six_vectors(const double (&vecs)[6][3], double (&result)[3]);
/* Sum six vectors and store result in "result" */

void mat_dot_vec(const double (&M)[3][3], const double (&V)[3], double(Result)[3]);
/* Matrix * Vector = Result */

void set_diag_mat(double (&M)[3][3], double a, double b, double c);
/* Set 3*3 matrix M to be: 
    a, 0, 0,
    0, b, 0,
    0, 0, c,
*/

void set_zero_vector(double (&V)[3]);
/* Set the vector V to be all zeros */

double max_element_arr(double (&a)[3]);
/* Returns max element in an array */

#endif