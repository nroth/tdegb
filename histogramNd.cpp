#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <vector>
#include "physical_constants.h"
#include "histogramNd.h"

using std::vector;
using std::string;

//***************************************************************
// Constructors
//***************************************************************

HistogramNd::HistogramNd()
{
  dimension = 1;
  n_bins_total = 0;
}

HistogramNd::HistogramNd(vector<vector<double> > bin_specs, vector<string> d_names, string b_name)
{
  Init(bin_specs, d_names, b_name);
}

void HistogramNd::Init(vector<vector<double> > bin_specs, vector<string> d_names, string b_name)
{
  
  // initalize grids
  dimension = bin_specs.size();
  bin_array.resize(dimension);
  n_elements.resize(dimension);
  dimension_names.resize(dimension);


  for (int i = 0; i < dimension; i++)
    {
      bin_array[i].New(bin_specs[i][0],bin_specs[i][1],bin_specs[i][2]);
      n_elements[i] = bin_array[i].Size(); // Size here is a locate_array method

      dimension_names[i] = d_names[i];

      if (i == 0)
	n_bins_total = n_elements[i];
      else
	n_bins_total *= n_elements[i];
    }

  base_name = b_name;

  
  // allocate
  //  count  = new int[n_elements];
  count.resize(n_bins_total);

  // clear 
  Wipe();
}

void HistogramNd::InitFrom(vector<vector<double> > bin_edges, vector<string> d_names, string b_name)
{
  
  // initalize grids
  dimension = bin_edges.size();
  bin_array.resize(dimension);
  n_elements.resize(dimension);
  dimension_names.resize(dimension);

  for (int i = 0; i < dimension; i++)
    {
      bin_array[i].New(bin_edges[i]);
      n_elements[i] = bin_array[i].Size(); // Size here is a locate_array method

      dimension_names[i] = d_names[i];

      if (i == 0)
	n_bins_total = n_elements[i];
      else
	n_bins_total *= n_elements[i];
    }

  base_name = b_name;

  count.resize(n_bins_total);

  Wipe();
}


void HistogramNd::Set_Dimension_Names(vector<string> d_names)
{
  if (d_names.size() != dimension)
    {
      printf("ERROR: the provided list of names does not match the histogram dimension in Set_Names()!\n");
    }
    for (int i = 0; i < dimension; i++)
    {
      dimension_names[i] = d_names[i];
    }
}

void HistogramNd::Set_Base_Name(string b_name)
{
  base_name = b_name;
}

void HistogramNd::Set_Dimension_Name(int axis, string d_name)
{
  if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Set_Dimension_Name()\n");
    }
  
  dimension_names[axis] = d_name;

}

string HistogramNd::Get_Dimension_Name(int axis)
{
  if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_Dimension_Name()\n");
    }
  return dimension_names[axis];
}

double HistogramNd::Get_Span(int axis)
{
  if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_Span()\n");
    }
  return bin_array[axis].Span();
}

double HistogramNd::Get_Low_Edge(int axis)
{
    if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_Low_Edge()\n");
    }
  return bin_array[axis].Min();
}

double HistogramNd::Get_High_Edge(int axis)
{
    if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_High_Edge()\n");
    }
  return bin_array[axis].Max();
}

double HistogramNd::Get_Bin_Width(int axis, int index)
{
  if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_Bin_Width()\n");
    }
  if (index < 0 || index >= n_elements[axis])
    {
      printf("ERROR: Invalid index in Get_Bin_Width()\n");
    }
  
  return bin_array[axis].Delta(index);
}

double HistogramNd::Get_Bin_Center(int axis, int index)
{

  if (axis < 0 || axis >= dimension )
    {
      printf("ERROR: Invalid axis in Get_Bin_Center()\n");
    }
  if (index < 0 || index >= n_elements[axis])
    {
      printf("ERROR: Invalid index in Get_Bin_Center()\n");
    }
  
  return bin_array[axis].Center(index);
}

double HistogramNd::Get_Bin_Count(vector<int> coord_indices)
{
  if (coord_indices.size() != dimension)
    {
      printf("ERROR: Invalid number of coordinate indices in Get_Bin_Count\n");
    }
  
  int flattened_index = Flatten_Index(coord_indices);
  return count[flattened_index];
}

double HistogramNd::Get_Bin_Count(int flattened_index)
{
  if (flattened_index < 0 || flattened_index >= n_bins_total)
    {
      printf("ERROR: Invalid index in Get_Bin_Count() \n");
    }

  return count[flattened_index];
}

int HistogramNd::Flatten_Index(vector<int> coord_indices)
{

  if (coord_indices.size() != dimension)
    {
      printf("ERROR: The provided list of coordinate indices does not match the histogram dimension in Flatten_Index()!\n");
    }

  int current_value = 0;
  for (int i = 0; i < coord_indices.size(); i++)
    {
      current_value = coord_indices[i] + n_elements[i] * current_value;
    }
  return current_value;
}

vector<int> HistogramNd::Unflatten_Index(int index)
{

  if (index < 0 || index >= n_bins_total)
    {
      printf("ERROR: The provided index is invalid in Unflatten_Index()\n");
    }
  int current_divisor = 1;
  vector<int> coordinates(dimension,0);
  for (int j = dimension - 1; j > -1; j--)
    {
      if (j == dimension - 1)
	{
	  coordinates[j] = index % n_elements[j];
	}
      else
	{
	  current_divisor *= n_elements[j+1];
	  coordinates[j] = (int)floor(index / current_divisor) % n_elements[j];
	}
    }

  return coordinates;
}

int HistogramNd::Get_Dimension()
{
  return dimension;
}

vector<string> HistogramNd::Get_Dimension_Names()
{
  return dimension_names;
}

string HistogramNd::Get_Base_Name()
{
  return base_name;
}

int HistogramNd::Get_Num_Bins(int axis)
{
  if (axis < 0 || axis >= dimension)
    {
      printf("ERROR: The provided axis is invalid in Get_Num_Bins()\n");
    }
  return n_elements[axis];
}

int HistogramNd::Get_Num_Bins_Total()
{
  return n_bins_total;
}


//***************************************************************
// Functional procedure: Wipe
//***************************************************************
void HistogramNd::Wipe()
{

  for (int i=0;i<n_bins_total;i++) {
    count[i]  = 0.;
  }
}



//***************************************************************
//--------------------------------------------------------------
//********************** COUNT FUNCTIONS ************************
//--------------------------------------------------------------
//***************************************************************


void HistogramNd::Count(vector<double> data)
{
  if (data.size() != dimension)
    {
      printf("ERROR: The provided data is invalid in Count()\n");
    }

  vector<int> bin_coordinates(dimension);
  for (int i = 0; i < dimension; i++)
    {
      bin_coordinates[i] = bin_array[i].Locate_Bounded(data[i]);
    }

  int bin_index = Flatten_Index(bin_coordinates);
  count[bin_index] += 1.;
}

void HistogramNd::Count(vector<double> data, double weight)
{
    if (data.size() != dimension)
    {
      printf("ERROR: The provided data is invalid in Count() with weight\n");
    }

  vector<int> bin_coordinates(dimension);
  for (int i = 0; i < dimension; i++)
    {
      bin_coordinates[i] = bin_array[i].Locate_Bounded(data[i]);
    }

  int bin_index = Flatten_Index(bin_coordinates);
  count[bin_index] += weight;
}

void HistogramNd::Count(vector<int> coordinates, double weight)
{

  if (coordinates.size() != dimension)
    {
      printf("ERROR: The provided data is invalid in Count() based on coordinates with weight\n");
    }
  
  int bin_index = Flatten_Index(coordinates);
  count[bin_index] += weight;
}

void HistogramNd::Set_All_Counts(vector<double> flattened_counts)
{

  count = flattened_counts;

}

vector<double> HistogramNd::Get_All_Counts()
{

  return count;

}


//***************************************************************
//--------------------------------------------------------------
//********************** PROJECTIONS  *************************
//--------------------------------------------------------------
//***************************************************************


HistogramNd HistogramNd::Create_Projected_Histogram(vector<int> kept_axes)
{

  int old_dimension = dimension;
  int new_dimension = kept_axes.size();

  if (new_dimension > old_dimension)
    {
      printf("ERROR: cannot project into more dimensions than before");
    }
  else if (new_dimension == old_dimension)
    {
      printf("WARNING: projecting into same size as before; histogram will just be copied");
    }
  else
    {
      printf("Projecting from dimension %d to dimension %d\n",old_dimension,new_dimension);
    }
  
  // First create the projected histogram based on the axes that will be kept
  vector<vector<double> > kept_bin_edges;
  vector<string> kept_names;
  for (int i =0; i < new_dimension; i++)
    {
      kept_bin_edges.push_back(bin_array[kept_axes[i]].Get_Edges());
      kept_names.push_back(dimension_names[kept_axes[i]]);
    }
  
  HistogramNd hist_new;
  hist_new.InitFrom(kept_bin_edges, kept_names, base_name);

  // Loop through all bins in original histogram and increment in new projected histogram
  vector<int> old_bin_coordinates(old_dimension);
  vector<int> new_bin_coordinates(new_dimension);

  for (int i = 0; i < n_bins_total; i++)
    {

      old_bin_coordinates = Unflatten_Index(i);
      for (int j = 0; j < new_dimension; j++)
	{
	  // the coordinates of this bin are a subset of the old coordinates
	  new_bin_coordinates[j] = old_bin_coordinates[kept_axes[j]];
	}

      hist_new.Count(new_bin_coordinates,count[i]);
      
    }

  return hist_new;

  
}



//***************************************************************
//--------------------------------------------------------------
//********************** PRINT FUNCTIONS *************************
//--------------------------------------------------------------
//***************************************************************

// not sure if specifying axes is the way to go here, might want to allow for arbitrary projections
void HistogramNd::Print_Histogram_1D(int axis)
{
  string filename = base_name + "_" + dimension_names[axis] + ".hist";
  FILE *out = fopen(filename.c_str(),"w");

  fprintf(out, "%-16s %-16s %-16s\n","#LOW", "HIGH", "COUNT");


  // NEED TO THINK MORE ABOUT HOW TO DO THIS FOR MULTI-D CASE!!
  // PROBABLY SHOULD PROJECT FIRST INTO HISTOGRAM1D
  
  for (int j=0;j<n_elements[axis];j++) 
    {
      fprintf(out,"%-16.6e %-16.6e %-16.6e\n",bin_array[axis].Low(j),bin_array[axis].High(j), count[j]);
    }
  fclose(out);
}


void HistogramNd::Print_Histogram_2D(int axis1, int axis2)
{
  string filename = base_name + "_" + dimension_names[axis1] + "_" + dimension_names[axis2] + ".hist";
  FILE *out = fopen(filename.c_str(),"w");

  // write out bin edges

  vector<double> edges1 = bin_array[axis1].Get_Edges();
  vector<double> edges2 = bin_array[axis2].Get_Edges();
  
  for (int i=0; i < edges1.size(); i++) 
      {
	fprintf(out," %-16.6e", edges1[i]);
      }
  fprintf(out,"\n");
  for (int i=0; i < edges2.size(); i++) 
      {
	fprintf(out," %-16.6e", edges2[i]);
      }
  fprintf(out,"\n");


  //write out counts as a matrix
  vector<int> output_coordinates(2);
  for (int i = 0; i < n_elements[axis1]; i++)
    {
      for (int j=0;j<n_elements[axis2]; j++) 
	{
	  int coords[] = {i,j};
	  output_coordinates.assign(coords, coords + 2);

	  int flattened_index = Flatten_Index(output_coordinates);
	  fprintf(out,"%-16.6e ",count[flattened_index]);
	}
      fprintf(out,"\n");
    }
  fclose(out);
}



//***************************************************************
//--------------------------------------------------------------
//********************** MPI FUNCTIONS **************************
//--------------------------------------------------------------
//***************************************************************

/*
void HistogramNd::Normalize()
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


void HistogramNd::Rescale(double r)
{
  for (int i=0;i<n_elements;i++) count[i] *= r;
}
*/

/*
void HistogramNd::MPI_Sum_All()
{
  MPI_Allreduce_Array(count);
  MPI_Allreduce_Array(click);
}

void HistogramNd::MPI_Average_All()
{
  MPI_Sum_All();

  int mpi_procs;
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_procs );
  for (int i=0;i<n_elements;i++) {
    count[i] /= mpi_procs; }
}


void HistogramNd::MPI_Allreduce_Array(double *arr)
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

