//  file: diffeq_routines.h
// 
//  Header file for diffeq_routines.cpp
//
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision History:
//    02/07/04 --- original version converted from C version, 
//                  based on rk4.cpp from "Computational
//                  Physics" and derivative_test.cpp
//    02/14/04 --- added 2nd order Runge-Kutta routine       
//
//  To do:
//
//************************************************************************

// function prototypes 
/* extern int euler ( const int N, double t, double y[], double h,
	    double (*f) (double t, double y[], int i, void *params_ptr), 
            void *params_ptr ); */
 
extern int runge4 ( const int N, double t, double y[], double h,
	    double (*f) (double t, double y[], int i, void *params_ptr), 
            void *params_ptr );
 
/* extern int runge2 ( const int N, double t, double y[], double h,
	    double (*f) (double t, double y[], int i, void *params_ptr), 
            void *params_ptr ); */
 
