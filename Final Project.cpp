  #include <iostream>
  #include <iomanip>
  #include <fstream>
  #include <cmath>
  #include "FinalProject.h"
  #include <stdlib.h>
  #include "Gnuplotpipe.h"

  using namespace std;

  double diff(double E, double h);
  inline double sqr (double x) {return x*x;};	// square a double
  class QuantumEigen      // define a type to hold parameters
  {
    public     //
          //
    const double eps = 1e-6;      //
    const double n_steps = 501;    //
  };    //
   /*private:
    double count;   //
    double count_max;  //
    double E; */
      // now we can define a structure of this type
          //   using the keyword "force_parameters"
  int main (void){


    double E =-17., h=0.04;
    // Initial E in MeV, step size in fm
    double  Emax, Emin, Diff;

    int count , count_max= 100;

   Emax= 1.1*E ;Emin=E/1.1;
   // Iteration loop
   for ( count =0; count <= count_max;count++){
     E =(Emax+Emin)/2.;              // Bisection
       Diff = diff (E,h);
       cout << "E =" << E << ", L-R Log deriv(E) =" << Diff << endl;
   };
     if (diff(Emax,h)*Diff > 0) Emax = E;
     else Emin = E ;
   if ( abs(Diff)<eps) break ;};
   plot(E,h);
    cout << "Final eigenvalue E = " << E << endl;
    cout << "iterations, max =" << count << "," << count_max << endl;
    cout << "WF in QuantumL/R.dat, V in QuantumV.dat";
  } // end main
  double diff(double E, double h) {// L-R log deriv

  double left, right, x;
  int ix, nL, nR, i_match;
  double y[] = new double[2];
  i_match = n_steps/3;
  nL =i_match +1;
  y[0] = 1.e-15;
  y[1] =y[0]* sqrt(-E*0.4829);

  for ( ix = 0; ix < nL+1 ; ix ++) {
    x= h* (ix -n_steps/2);
    void rk4 (x,y,h,2,E);
  }
  left = y[1]/y[0];
  y[0] =1.e-15;
  y[1] = -y[0]* sqrt (-E* 0.4829);

  for ( ix = n_steps ; ix >nL +1; ix --)
  {
    x = h* (ix+1 - n_steps/2);
    rk4 (x,y, -h ,2 , E);
    }
    right = y[1]/y[0];
    return( (left - right)/(left + right));}
    // Repeat integrations for plot, can't integrate out decaying wf
     void plot ( double E, double h)
    //ofstream.out
    /* PrintWriter L =
          new PrintWriter (new FileOutputStream("QuantumL.dat"),true);
    PrintWriter R=
          new PrintWriter (new FileOutputStream ("QuantuR, dat"), true);
    PrintWriter Vx =
          new PrintWriter (new FileOutputStream ("QuantumV.dat"), true); */
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
    double y[]= new double[2]; double yL[][]= new double[2][505];
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
y[0]= -1.e-15; // -slope : even ; reverse for odd
y[1] = - sqrt(-E*0.4829)* y[0]; // Initialise R wf
for ( ix = n_steps -1 ; ix >= nL+1 ; ix --) { // Integrate in
  x = h * (ix +1 - n_steps/2);
  /* R.println (x + " " + y[0] +" " + y[1]);  // File print
  Vx.println(x+"  " + 1.e 34*V(x));  // Scaled V
  rk4 (x,y, -h, 2,E);

}
  x = x - h ;
  R.println(x+ "  " + y[0]+ " "+y[1]);  */// File print
    my_out << "# x        y[0]     y[1]"<<setprecision(12)<< R << endl;
    my_out << "# x        1.e*34*V(x)"  <<setprecision(12)<< Vx <<endl;
    rk (x,y,-h, 2, E); }
    x= x - h;
    my_out << "# x        y[0]     y[1]"<<setprecision(12)<< R <<endl;

  normL = y[0]/ yL[0][nL];              //Renormalize L wf and derivative
  for (ix = 0 ; ix < = nL ; ix ++) {
    x = h* (ix -n_steps/2 +1);
    y[0]= yL[0][ix] *normL;
    y[1]=yL[1][ix]*normL;
    my_out << "# x        y[0]      y[1]"<< setprecision(12) << L << endl;
    my_out << "# x        y[0]      y[1]"<< setprecision(12) << Vx << endl;
    my_out << "# x        1.e 34 * V(x)" << setprecision(12) << Vx << endl;
     /* L.println(x+ " " +y[0]+ " " + y[1]);
    Vx.println(x+ "  "+ y[0] + "  "+ y[1]);
     Vx.println(x + "  "+ 1.e 34 *V(x));  } */
     return;
   }
    void f(double x, double y[], double F[], double E)
   { F[0] = y[1]; F[1] = -(0.4829)*(E- V(x)*y[0]; }

   double V(double x)
   { if (abs(x) < 10.)} return ( -16.); else return(0.);}
                                            //rk4 algorithm
     void rk4(double t, double y[], double h,)
                 int Neqs, double E) {
                   int i;
                   double F[] = new double [Neqs];
                   double ydumb[] = new double [Neqs] ;
                   double k1[] = new double [Neqs]; double k2[] = new double[Neqs];
                   double k3[] = new double[Neqs] ; double k4[]; new double[Neqs];
                   f(t,y,F, E);
                   for (i=0; i< Neqs; i++)
                   {k1[i]=h*F[i];
                   ydumb[i]= y[i]+ k1[i]/2;}

                   f(t+h/2,ydumb, F, E);
                   for (i=0; i< Neqs; i++)
                   { k2[i] = h*F[i];
                     ydumb[i] = y[i] +k2[i]/2;}

            f(t+h/2,ydumb, F, E);
            for (i=0 ; i <Neqs, i++)
            {k3[i]=h*F[i];
            ydumb[i]= y[i]+ k3[i];

          f(t+h, ydumb,F,E);
          for (i=0; i < Neqs ; i++)
          { k4[i]=h*F[i];
            y[i]=y[i]+ (k1[i]+2*(k2[i]+k3[i])+k4[i])/6;}
          }

                   }
