#ifndef HISTOGRAM1D_H
#define HISTOGRAM1D_H 

#include <string>
#include "locate_array.h"
using std::string;

// default values
#define DEFAULT_NAME "histogram"


// This is copied from the spectrum.hh class, but simplified. What you should really do is let spectrum inherit from this class.
// Doesn't have a mu, phi, or time dependence
// Something else that might be a little confusing: in the original spectrum class, there were two arrays, count and click. Oddly, the old click array seemed more like it performed a counting function to me (the old count array accumulated total photon energies). So I have the count array here doing the sort of counting that the old click array used to do.

class Histogram1d {
 
private:

  // histogram name
  string name;
  string datadir;

  // number of elements and index helpers
  int n_elements;

  // bin array
  LOCATE_ARRAY bin_array;
  
  // counting array
  //  int *count;
  double *count; //using a double because sometimes you'll have weighted counts

  // Used for MPI
  //  void MPI_Allreduce_Array(double *arr);
    
public:

  // constructors
  Histogram1d();
  
  // Initialize
  void Init(double, double, double);
  void Set_Name(string);
  void Set_Data_Directory(string f) {  datadir = f; }
  // MPI functions
  //  void MPI_Sum_All();
  //  void MPI_Average_All();
  
  // Count and Normalize counted packets
  void Count(double);
  void Count(double,double); // with a weight
  //  void Normalize();
  //  void Rescale(double);
  //  void Get_Magnitude(int t, int m, double *result);
  //  void Bolometric_Luminosity(int t, int m, double *result);

  // Print out
  //  void Print();
  void Print_Histogram();

  //  double& operator()(int m, int p, int l, int s);
  //  double& operator()(int m, int p, int l);
  
  void Wipe();
  //  double N_Escaped() {return n_escaped;}
};

#endif
