//  file: diffeq_routines.cpp
//
//  Routines for Euler's and 4th order Runge-Kutta diff. eq. routines
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//   30-Jan-2004 --- original version, translated from diffeq_routines.c,
//                    which was based on rk4.c from
//                    "Computational Physics" by Landau and Paez
//   14-Feb-2004 --- added 2nd order Runge-Kutta routine
//   30-Jan-2005 --- comments improved and function names changed
//   22-Jan-2006 --- made i local to loops
//
//   * Based on the discussion of differential equations in Chap. 9
//      of "Computational Physics" by Landau and Paez
//   * As a convention (advocated in "Practical C++"), we'll append
//      "_ptr" to all pointers.
//
//  To do:
//
//************************************************************************

// include files
#include <iostream>		// note that .h is omitted
#include <iomanip>		// note that .h is omitted
#include <fstream>		// note that .h is omitted
#include <cmath>
#include "diffeq_routines.h"	// diffeq routine prototypes

const int NMAX=5;   // maximum of NMAX coupled equations



double V (double x, params_ptr)
{ if(abs(x)<10.) return (16.); else return(0);}


double f (double x,double y[],double F[],double E, params_ptr)
{
  F[0]=y[1]; F[1]=-(0.4829)*(E-V(x))*y[0];
 return (0);
}

//************************************************************************
//
//   4th Order Runge-Kutta Differential Equation Solver
//
// This routine takes all of the y's one step, from t to t+h.
//  The original values of y[0], y[1], etc. are lost.
//
// inputs:
//   N --- number of y(t)'s
//   t --- independent variable
//   y[] --- vector of y(t)'s
//   h --- step size
//   f --- function for the right hand sides
//   *params_ptr --- pointer to parameters for rhs function f
//
// outputs:
//   y[] --- predictions for the values of y(t+h)
//
// Notes:
//   * The algorithm is from Eq.(9.46) in Landau and Paez.
//
//***********************************************************************
int
 runge4 (const int N, double t, double y[], double h,
    double (*f) (double t, double y[], int i, void *quantum_params_ptr),
            void *quantum_params_ptr);
{

  double NMAX = ((quantum_parameters *) quantum_params_ptr)->Neqs;
  double t = ((quantum_parameters *) quantum_params_ptr)->t;
  double y[] = ((quantum_parameters *) quantum_params_ptr)->y[];
  double h = ((quantum_parameters *) quantum_params_ptr)->h;
  double i = ((quantum_parameters *) quantum_params_ptr)->i;//loading parameters into the function
  double (*f) = ((quantum_parameters *)quantum_params_ptr)-> f;
  double y1[NMAX], y2[NMAX], y3[NMAX];	// intermediate y values
  double k1[NMAX], k2[NMAX], k3[NMAX], k4[NMAX]; // Runge-Kutta notation
  double f(t,y,F,E)
  for (int i = 0; i < N; i++)
    {
      k1[i] = h * f (t, y, i, params_ptr);
      y1[i] = y[i] + k1[i] / 2.;	// argument for k2
    }

  for (int i = 0; i < N; i++)
    {
      k2[i] = h * f (t + h / 2., y1, i, params_ptr);
      y2[i] = y[i] + k2[i] / 2.;	// argument for k3
    }

  for (int i = 0; i < N; i++)
    {
      k3[i] = h * f (t + h / 2., y2, i, params_ptr);
      y3[i] = y[i] + k3[i];	// argument for k4
    }

  for (int i = 0; i < N; i++)
    {
      k4[i] = h * f (t + h, y3, i, params_ptr);
    }

  for (int i = 0; i < N; i++)
    {
      y[i] += (k1[i] + 2. * k2[i] + 2. * k3[i] + k4[i]) / 6.0;
    }

  return (0);			// successful completion
}
