#ifndef VECTORFIELD_H
#define VECTORFIELD_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>

using namespace std;

/* Custom Class for a size-3 vector field
  - Supports storing, reading and modifying the field data on a uniform 2D mesh.
  - Includes ghost cells
*/
class VectorField
{
private:
  double *array; // pointer to the array

public:
  double Lx;           // length of domain in the x-direction
  double Ly;           // height of domain in the y-direction
  int Ni;              // number of cells on the x-direction (excluding ghost cells)
  int Nj;              // number of cells on the y-direction (excluding ghost cells)
  int Nk = 3;          // length of vectors
  double Dx = Lx / Ni; // Delta x
  double Dy = Ly / Nj; // Delta y

  VectorField(double lx, double ly, int ni, int nj);
  /* Constructor: initiate a field with all zero data.
    Parameters: 
      - lx: length of domain in the x-direction
      - ly: height of domain in the y-direction
      - ni: number of cells on the x-direction (excluding ghost cells)
      - nj: number of cells on the y-direction (excluding ghost cells)
    Returns: Null
  */

  VectorField(double lx, double ly, int ni, int nj, double (*f[3])(double, double));
  /* Constructor: initiate a field by the given field functions.
    Parameters: 
      - lx, ly, ni, nj: same as above
      - (*f[3]): an array of field functions to set the vector field
    Returns: Null
  */

  ~VectorField();
  /* Destructor: delete the VectorField instance from memory.
    Parameters: 
      None
    Returns: Null
  */

  VectorField &operator+=(const VectorField &dU);
  /* Element-wise "-=" operation on every data values */

  VectorField &operator-=(const VectorField &dU);
  /* Element-wise "+=" operation on every data values */

  VectorField &operator*=(const double scalar);
  /* Element-wise "*=" operation on every data values */

  double at(int k, int i, int j) const;
  /* Method to access data at certain location.
    Parameters:
      - k: index of the vector component (0, 1, or 2) 
      - i: index of the x direction (from 0 to ni+1)
      - j: index of the y direction (from 0 to nj+1)
    Returns: data of component k at the location (i, j)
  */

  double P(int i, int j) const;
  /* Shorthand for calling at(1, i, j);
  */

  double u(int i, int j) const;
  /* Shorthand for calling at(1, i, j);
  */

  double v(int i, int j) const;
  /* Shorthand for calling at(1, i, j);
  */

  void set(int k, int i, int j, double value);
  /* Method to update the data at certain location to a new value.
    Parameters: 
      - i, j, k: same as in at(k, i, j)
      - value: new value for the data
    Returns: Null
  */

  void set(int i, int j, double values[3]);
  /* Method to update the data at certain location to a new value.
    Parameters: 
      - i, j: same as in at(k, i, j)
      - values: an array holding the values of the updated vector
    Returns: Null
  */

  void set_fields(double (*f[3])(double, double));
  /* Method to update the data in the mesh cells according to given function fields.
    Parameters: 
      - f: (*f[3]): an array of field functions to set the vector field
    Returns: Null
  */

  void calculate_L2_error(VectorField &Exact);
  /* Method to calculate the L2 error norms of each component compared to a given exact solution, and print them to screen (the ghost cells are not included, only the inner field values).
    Parameters: 
      - Exact: VectorField instance that holds the exact solution
    Returns: Null
  */

  double calculate_L2_norm(int k);
  /* Method to calculate the L2 norm of selected component (the ghost cells are not included, only the inner field values).
    Parameters: 
      - k: index of the component (0, 1, or 2)
    Returns: Null
  */

  void print(int k);
  /* Method to print the field of one component to screen - for debugging only.
    Parameters: 
      - k: index of the vector component (0, 1, or 2) 
    Returns: Null
  */

  void print_partial(int k);
  /* Method part of the field to screen - for debugging problem 1.2.
    Parameters: 
      None 
    Returns: Null
  */

  void write_dat(string fileName);
  /* Method to write all data in the format similar to the given comparison data.
    Parameters:  
      - fileName: file name to save the data.
    Returns: Null
  */

  void write_csv(string fileName);
  /* Method to write each component of the VectorField into a separate csv file.
    Parameters:  
      - fileName: file name to save the data. "P", "u", "v" will be appended to the fileName to indicate which component it is.
    Returns: Null
  */

  void output_midline_u();
  /* Method to output the solution of u at x=1/2 into csv file named "u_midline-[Ni]by[Nj].csv"
    Parameters:  
      None
    Returns: Null
  */

  void output_vorticity(string fileName);
  /* Method to output the vorticity field.
    Parameters:  
      - fileName: the output will will be named "vorticity-[fileName].csv"
    Returns: Null
  */

  double x(int i);
  /* Method to find the x position given index i.
    Parameters: 
      - i: the given index
    Returns: position x
  */

  double y(int j);
  /* Similar as x(i) but for y direction */

protected:
  int index(int k, int i, int j) const;
};

#endif