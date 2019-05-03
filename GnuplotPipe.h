//  file: GnuplotPipe.h
//
//  Header file for the GnuplotPipe C++ class.
//
//  Programmer:  Dick Furnstahl  furnstahl.1@osu.edu
//
//  Revision history:
//      02/06/06  original version, based on gnuplot_pipe 
//      02/07/09  minor upgrades
//
//*****************************************************************
#ifndef GNUPLOTPIPE_H
#define GNUPLOTPIPE_H

// include files
#include <iostream>
#include <fstream>
#include <string>
using std::string;

class GnuplotPipe
{
 public:
  GnuplotPipe ( );   // constructor
  ~GnuplotPipe ( );  // destructor
  
  // accessor functions --- these set private variables
  void set_filename (const string &t_filename) {filename = t_filename;};
  void set_filename2 (const string &t_filename2) {filename2 = t_filename2;};
  void set_title (const string &t_title) {title = t_title;};
  void set_xlabel (const string &t_xlabel) {xlabel = t_xlabel;};
  void set_ylabel (const string &t_ylabel) {ylabel = t_ylabel;};
  void set_plot_title (const string &t_plot_title) {plot_title = t_plot_title;};
  void set_plot_title2 (const string &t_plot_title2) 
                            {plot_title2 = t_plot_title2;};
  void set_xmin (const double t_xmin) {xmin = t_xmin;};
  void set_xmax (const double t_xmax) {xmax = t_xmax;};
  void set_ymin (const double t_ymin) {ymin = t_ymin;};
  void set_ymax (const double t_ymax) {ymax = t_ymax;};
  void set_delay (const int t_delay) {delay = t_delay;};

  int init ();                // initialize a plot
  int plot (const double x, const double y);  // plot a point on first line
  int plot2 (const double x, const double y);  // plot a point on second line
  int finish ();              // close down a gnuplot session
  int gnuplot_cmd (const string &plot_cmd);   // send a command to gnuplot

 private:
  FILE *gp_cmd;         // file handle for gnuplot process
  FILE *fileout;        // file handle for temporary output file 1
  string filename;      // file name for temporary output file 1
  FILE *fileout2;       // file handle for temporary output file 2
  string filename2;  // file name for temporary output file 2
  string title;      // overall plot title
  string xlabel;     // x-axis label
  string ylabel;     // y-axis label  
  string plot_title;  // title for first line
  string plot_title2;  // title for second line
  double xmin;    // lower x-limit of graph
  double xmax;    // upper x-limit of graph
  double ymin;    // lower y-limit of graph
  double ymax;    // upper y-limit of graph
  int delay;      // plotting delay (in microseconds)
  int plot2_flag;       // flag to indicate whether 2nd plot has started
  string plot_cmd;      // command to plot only the first plot
  string plot_cmd2;     // command to plot both plots
};

#endif
