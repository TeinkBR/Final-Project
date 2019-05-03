#include "FINAL_PROJECT.h"
#include <cmath>
int plot ( double E, double h){

     ofstream my_out;    // now open a stream to a file for output
     my_out.open(my_stringstream.str().c_str());


     // print initial values and column headings
     ostringstream my_stringstream;  // declare a stringstream object
     my_stringstream << "Quantum"  << setprecision(8)
                     << vdp_rhs.get_mu()
                     << "_L _" << setprecision(8) << L
                     << "_R_" << setprecision(8) << R
                     << "_V_" << setprecision(8) << Vx <<".dat";
     ofstream my_out;    // now open a stream to a file for output
     //my_out.open(my_stringstream.str().c_str());

double left, right, normL, x= 0. ;
int ix,nL, nR, i_match;
double y[2]; double yL[2][505];
int n_steps =1501;     // Total No. integration steps
i_match = 500 ;       // Matching point
nL = i_match +1 ;
y[0] =1.e-40;           // initial wf on the left
y[1] = - sqrt (-E*0.4829)*y[0];
for (ix =0; ix <= nL; ix ++)
{yL[0][ix]=y[0];
yL[1][ix]=y[1];
x= h* (ix -n_steps/2);
rk4(x,y,h,2,E);}  // integrate to the left
