  #include <iostream>
  #include <iomanip>
  #include <fstream>
  #include <cmath>
  #include <string>
  using namespace std;
  #include "GnuplotPipe.h"
  #include "diffeq_routines.h"


  double runge4 (double t, double y[], double h,
             int Neqs, double E, void *params_ptr);
  double diff(double E, double h, void *params_ptr);

  int plot(double E,double h, void *params_ptr);


  const double pi = M_PI;  // use the system-defined 3.14159...
  inline double sqr (double x) {return x*x;};	// square a double

  // structures
  typedef struct      // define a type to hold parameters
  {
    double Emax;    // natural frequency
    double Emin;      // coefficient of friction
    double E;      // amplitude of external force
    double count;    // frequency of external force
    double count_max;    // phase angle for external force
    double eps;
    double n_steps;
    double Diff;
    double Diff_max;
    double Neqs;
    double h;
    double t;
    double y;
    int i;

    double left, right, x;
    int ix, nL, nR, i_match;
  }
  quantum_parameters;    // now we can define a structure of this type
          //   using the keyword "force_parameters"

//**********************main program***********************************

int main (void){


    void *quantum_params_ptr;
    quantum_parameters ode_parameters;
    //int i;
    double n_steps = 501.;
    double eps = 1e-6;
    double E =-17., h=0.04;
    //double t;
    double  Emax, Emin, Diff, Diff_max;
    int count , count_max= 100;
    double Neqs=2;
    double y[2];
    Emax= 1.1*E ;Emin=E/1.1;




    ode_parameters.Emax = Emax;
    ode_parameters.Emin= Emin;
    ode_parameters.E = E;
    ode_parameters.count = count;
    ode_parameters.count_max = count_max;
    ode_parameters.eps=eps;
    ode_parameters.n_steps=n_steps;
    ode_parameters.Diff= Diff;
    ode_parameters.h=h;
    ode_parameters.y=y[2];
    //ode_parameters.t=t;
    //ode_parameters.i=i;
    ode_parameters.Neqs= Neqs;
    quantum_params_ptr = & ode_parameters;  // structure to pass to function



      for (count = 0; count <= count_max; count++){
       E =(Emax+Emin)/2.;              // Bisection

         Diff_max = diff(E ,h,quantum_params_ptr);
        ode_parameters.Diff_max= Diff_max;
       cout << "E =" << E << ", L-R Log deriv(E) =" << Diff << endl;

     if (Diff_max*Diff > 0) Emax = E;
     else Emin = E ;
     //QuantumEigen myQuantumEigen;
     double temp_eps = ode_parameters.eps;
   if ( abs(Diff) < temp_eps)
   break ;
    }
   plot(E,h, quantum_params_ptr);
    cout << "Final eigenvalue E = " << E << endl;
    cout << "iterations, max =" << count << "," << count_max << endl;
    cout << "WF in QuantumL/R.dat, V in QuantumV.dat";
  return 0;
  }
//*************************diff.cpp*****************************************


double diff (void *params_ptr) {




    double E = ((quantum_parameters *) params_ptr)-> E;

    double y[2] = ((quantum_parameters *) params_ptr)->y[];
    double h = ((quantum_parameters *) params_ptr)->h;
    int nL = ((quantum_parameters *) params_ptr)->nL; //loading parameters into the function
    double i_match = ((quantum_parameters *) params_ptr)->i_match;
    double n_steps =((quantum_parameters *) params_ptr)->n_steps;
    double ix = ((quantum_parameters*)params_ptr)->ix;
    int x = ((quantum_parameters*)params_ptr)->x;

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
        runge4 (x,y,-h ,2 ,E, params_ptr);
     }
        right = y[1]/y[0];
        return( (left - right)/(left + right));}

//****************************rk4*******************************************************
int
 runge4 (const int N, double (*f) (void *params_ptr),
            void *params_ptr)
{

  double Neqs = ((quantum_parameters *) params_ptr)->NMAX;
  double t = ((quantum_parameters *) params_ptr)->t;
  double y[] = ((quantum_parameters *) params_ptr)->y[];
  double h = ((quantum_parameters *) params_ptr)->h;
  double i = ((quantum_parameters *) params_ptr)->i;//loading parameters into the function
  double (*f) = ((quantum_parameters *)params_ptr)-> f;
  double y1[NMAX], y2[NMAX], y3[NMAX];	// intermediate y values
  double k1[NMAX], k2[NMAX], k3[NMAX], k4[NMAX]; // Runge-Kutta notation
  double f(t,y,F,E,params_ptr)
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
//************************plot****************************************
int plot ( double E, double h, void *params_ptr){


      h = ((quantum_parameters *) params_ptr)->h;
      E = ((quantum_parameters *) params_ptr)->E;
 ofstream my_outL("QuantumL.dat");
 ofstream my_outR("QuantumR.dat");
 ofstream my_outVx("QuantumV.dat");


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
runge4(x,y,h,2,E);}  // integrate to the left
y[]=1.e-15;
y[1]=sqrt(-e*0.4829)*y[0];
for(ix=n_steps-1;ix>=nL+1;ix--){
  x=h*(ix+1-n_steps/2);
  my_outR<<"#x     y[0]     y[1]"<<endl;
  my_outR<<setprecision(12)<<x<<setprecision(12)<<y[0]<<setprecision(12)<<y[1]<<endl;
  my_outVx<<"#x    1.e34*V(x)"<<setprecision(12)<<x<<setprecision(12)<<1.e34*Vx<<endl;
  runge4 (x,y,-h,2,E,params_ptr);

}
x=x-h;
my_outR<<"#x     y[0]     y[1]"<<endl;
my_outR<<setprecision(12)<<x<<setprecision(12)<<y[0]<<setprecision(12)<<y[1]<<endl;
normL = y[0]/yL[0][nL];
for (ix=0; ix <= nL; ix++){
x= h*(ix-n_steps/2+1);
y[0]=yL[0][ix]*normL;
y[1]=yL[1][ix]*normL;
my_outL<<"#x     y[0]      y[1]"<<setprecision(12)<<x<<setprecision(12)<<y[0]<<setprecision(12)<<y[1]<<endl;
my_outVx<<"#x    1.e34*V(x)"<<setprecision(12)<<x<<setprecision(12)<<1.e34*Vx<<endl;


}
return;
}
int f (double x, double y[],double F[],double E)
{F[0]=y[1]; F[1]=-(0.4829)*(E-V(x))*y[0];}
double V(double x)
{if(abs(x)<10.) return(16.); else return(0.);}
