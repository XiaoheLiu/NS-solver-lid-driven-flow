# CFD Project: Lid Driven Cavity

Athena Liu

December, 2019

## Project Description
In this project, I wrote a program in C++ to solve the two-dimensional incompressible laminar Navier-Stokes equations on regular Cartesion grids using finite volume method. The solver is applied to simulate the lid driven cavity flow (a rectangular cavity with horizontally-moving top). See `./Problem_Description.pdf` for a detailed description of the problem.

## Results
The steady-state flow fields are output as `csv` files in `/src/outputs`, and post-processing is done in `/src/Plotting.ipynb` using Python / Jupyter Notebook.

The analysis of the results and the flow visualizations are presented in `./Report.pdf`.


## Files in `/src`:
- __main.cpp__: Entry point of the code; runs the main program.
- __settings.h__: Sets up parameters for the computation and physical problem.
- __VectorField.h/cpp__: Custom class for a vector field.
- __field_functions.h/cpp__: Functions that define the initial and boundary fields.
- __navier-stokes.h/cpp__: Functions for the 2D incompressible Navier-Stokes solver.
- __jacobians.h/cpp__ : Functions to compute the flux Jacobians, and a few helper functions for matrix manipulations.
- __tri_thomas__ : Function to solve a block tri-diagonal system using the Thomas algorithm.

On a Mac, compile and run the program using
```
source ./src/run.sh
```

## Further Note
This project was done as the final project for the course _MECH 510:  Computational Methods in Transport Phenomena I_ at the University of the British Columbia, taught by the awesome Dr. Carl Ollivier-Gooch.