#ifndef HISTOGRAMND_H
#define HISTOGRAMND_H 

#include <string>
#include "locate_array.h"

using std::string;
using std::vector;

// default define
#define DEFAULT_NAME "histogramNd"

class HistogramNd {
 
private:

  string base_name; // histogram name
  vector<string> dimension_names;
  string datadir;

  int dimension;

  vector<int> n_elements;

  int n_bins_total; // total number of bins in the multidimensional grid. 

  // bin array
  vector<LOCATE_ARRAY> bin_array; // one bin array for each dimension
  
  vector<double> count; // The actual count data. Flattened. Using doubles because sometimes you'll have weighted counts

  // Used for MPI
  //  void MPI_Allreduce_Array(double *arr);

    
public:

  // constructors
  HistogramNd();
  
  HistogramNd(vector<vector<double> >, vector<string>, string);
  
  // Initialize
  void Init(vector<vector<double> >, vector<string>, string);
  void InitFrom(vector<vector<double> >, vector<string>, string);
  void Set_Dimension_Names(vector<string>);
  void Set_Dimension_Name(int, string);
  void Set_Base_Name(string);
  void Set_Data_Directory(string f) {  datadir = f; }

  int Get_Dimension();
  int Get_Num_Bins(int);
  int Get_Num_Bins_Total();
  vector<string> Get_Dimension_Names();
  string Get_Dimension_Name(int);
  string Get_Base_Name();

  double Get_Span(int);
  double Get_Low_Edge(int);
  double Get_High_Edge(int);
  double Get_Bin_Width(int, int);
  double Get_Bin_Center(int, int);
  double Get_Bin_Count(vector<int>); // use bin coordinates
  double Get_Bin_Count(int); // use flattened index

  int Flatten_Index(vector<int>);
  vector<int> Unflatten_Index(int);
  
  // MPI functions
  //  void MPI_Sum_All();
  //  void MPI_Average_All();

 
  // Count and Normalize counted packets
  void Count(vector<double>);
  void Count(vector<double>,double); // with a weight
  void Count(vector<int>,double); // based on coordinates with a weight
  //  void Normalize();
  //  void Rescale(double);

  void Set_All_Counts(vector<double>); // when passed a flattened counts array

  vector<double> Get_All_Counts();

  //Projections
  HistogramNd Create_Projected_Histogram(vector<int>);
  
  // Print out
  void Print_Histogram_1D(int axis);  // not sure if specifying axes is the way to go, might want to allow arbitrary projections
  void Print_Histogram_2D(int axis1, int axis2); 

  //  double& operator()(int m, int p, int l, int s);
  //  double& operator()(int m, int p, int l);

  
  void Wipe();
  //  double N_Escaped() {return n_escaped;}
 
};

#endif
