#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "physical_constants.h"
#include "histogram1d.h"


//***************************************************************
// Constructors
//***************************************************************

Histogram1d::Histogram1d()
{
  name = DEFAULT_NAME;
  n_elements = 0;
}

void Histogram1d::Init(double lb, double ub, double db)
{
  
  // initalize grids
  bin_array.New(lb,ub,db);
  n_elements = bin_array.Size();

  // allocate
  //  count  = new int[n_elements];
  count  = new double[n_elements]; // using a double because sometimes you'll have weighted counts


  // clear 
  Wipe();
}


void Histogram1d::Set_Name(string n)
{
  name = n;
}

//***************************************************************
// Functional procedure: Wipe
//***************************************************************
void Histogram1d::Wipe()
{

  for (int i=0;i<n_elements;i++) {
    //    count[i]  = 0;
    count[i]  = 0.;

  }
}


//***************************************************************
//--------------------------------------------------------------
//********************** COUNT FUNCTIONS ************************
//--------------------------------------------------------------
//***************************************************************

void Histogram1d::Count(double data)
{

  // locate bin number
  int bin_index = bin_array.Locate_Bounded(data);

  count[bin_index] += 1.;
}

void Histogram1d::Count(double data, double weight)
{

  // locate bin number
  int bin_index = bin_array.Locate_Bounded(data);

  count[bin_index] += weight;
}


//***************************************************************
//--------------------------------------------------------------
//********************** PRINT FUNCTIONS *************************
//--------------------------------------------------------------
//***************************************************************
void Histogram1d::Print_Histogram()
{
  string filename = name + ".hist";
  FILE *out = fopen(filename.c_str(),"w");

  fprintf(out, "%-16s %-16s %-16s\n","#LOW", "HIGH", "COUNT");


  // NEED TO THINK MORE ABOUT HOW TO DO THIS FOR MULTI-D CASE!!
  // PROBABLY SHOULD PROJECT FIRST INTO Histogram1d
  
  for (int j=0;j<n_elements;j++) 
    {
      fprintf(out,"%-16.6e %-16.6e %-16.6e\n",bin_array.Low(j),bin_array.High(j), count[j]);
    }
  fclose(out);
}


//***************************************************************
//--------------------------------------------------------------
//********************** MPI FUNCTIONS **************************
//--------------------------------------------------------------
//***************************************************************

/*
void Histogram1d::Normalize()
{
  // renormalize flux
  double norm_factor =  1.0/(1.0*mu_grid.Size()*phi_grid.Size());
  for (int i=0;i<time_grid.Size();i++) 
    for (int j=0;j<lambda_grid.Size();j++) 
      for (int m=0;m<mu_grid.Size();m++) 
	for (int p=0;p<phi_grid.Size();p++) 
	  {
	      count[Index(i,j,m,p)] /= 
		time_grid.Delta(i)*lambda_grid.Delta(j)*norm_factor;

	  }
}


void Histogram1d::Rescale(double r)
{
  for (int i=0;i<n_elements;i++) count[i] *= r;
}
*/

/*
void Histogram1d::MPI_Sum_All()
{
  MPI_Allreduce_Array(count);
  MPI_Allreduce_Array(click);
}

void Histogram1d::MPI_Average_All()
{
  MPI_Sum_All();

  int mpi_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  for (int i=0;i<n_elements;i++) {
    count[i] /= mpi_procs; }
}


void Histogram1d::MPI_Allreduce_Array(double *arr)
{
  double *new_ptr,*this_ptr;
  int j;          

  // do it in lambda chunks
  int chunk = a1;
  // allocate the memory for new pointer
  new_ptr = new double[chunk];
         
  for (int i=0;i<time_grid.Size();i++)
  {
    // zero out array
    for (j=0;j<chunk;j++) new_ptr[j] = 0;
    // reduce the stuff
    this_ptr = &(arr[i*chunk]);
    MPI_Allreduce(this_ptr,new_ptr,chunk,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    // put back into place
    for (j=0;j<chunk;j++) arr[i*chunk + j] = new_ptr[j];
  }
  // free up the memory
  delete new_ptr;
}

*/
