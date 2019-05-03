#include "FINAL_PROJECT.h"
#include <cmath>

double n_steps;
double left, right, x;
int ix, nL, nR, i_match;


double diff (double E, double h, void *params_ptr) {

    double n_steps;
    double left, right, x;
    int ix, nL, nR, i_match;
    double y[]=y[2];
    double NMAX = ((quantum_parameters *) params_ptr)->Neqs;
    double t = ((quantum_parameters *) params_ptr)->t;
    double y[] = ((quantum_parameters *) params_ptr)->y[];
    double h = ((quantum_parameters *) params_ptr)->h;
    double i = ((quantum_parameters *) params_ptr)->i; //loading parameters into the function

    i_match = n_steps/3;
    nL =i_match +1;
    y[0] = 1.e-15;
    y[1] =y[0]* sqrt(-E*0.4829);

    for ( ix = 0; ix < nL+1 ; ix ++) {
          x= h * (ix - n_steps/2);}
          left = y[1]/y[0];
          y[0] =1.e-15;
          y[1] = -y[0]* sqrt (-E* 0.4829);

    for ( ix = n_steps ; ix >nL +1; ix --)
    {
        x = h* (ix+1 - n_steps/2);
        runge4 (x,y,-h ,2 ,E,);
     }
        right = y[1]/y[0];
        return( (left - right)/(left + right));}
