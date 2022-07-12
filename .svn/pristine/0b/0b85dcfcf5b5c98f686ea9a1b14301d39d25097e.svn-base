#include "Double3DArray.h"
#include <cstdlib>
// #include <iostream>
// using std::cout;
// using std::endl;

Double3DArray::Double3DArray()
{
  this->dim1=0;
  this->dim2=0;
  this->dim3=0;
  this->dim23=0;
  this->dim=0;
  this->Flag.Initialize();
  this->data=NULL;
}


Double3DArray::Double3DArray(int dim1, int dim2, int dim3)
{
  this->dim1=dim1;
  this->dim2=dim2;
  this->dim3=dim3;
  this->dim=dim1*dim2*dim3;
  this->dim23=dim2*dim3;
  this->Flag.Initialize();
  this->data=new double[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}

Double3DArray::~Double3DArray()
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete [] this->data;
}

// copy constructor
//
// array = array to copy
// duplicateFlag = true if datas have to be duplicated
Double3DArray::Double3DArray(const Double3DArray& array, bool duplicateFlag) 
{
  this->dim1 = array.dim1;
  this->dim2 = array.dim2;
  this->dim3 = array.dim3;
  this->dim = array.dim;
  this->dim23 = array.dim23;
  if (array.dim == 0)
    {
      this->Flag.Initialize();
      this->data = NULL;
    }
  else
    if (duplicateFlag == false)
      {
	this->Flag = array.Flag;
	this->data = array.data;
      }
    else
      {
	this->Flag.Initialize();
	this->data = new double [dim];
	for (int i = 0; i < this->dim; i++)
	  this->data[i] = array.data[i];
      }
}

// assignment
//
// array = array to assign
// return value = reference on current array
Double3DArray& Double3DArray::operator = (const Double3DArray& array) 
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete[] this->data;
  this->Flag = array.Flag;
  this->data = array.data;
  this->dim1 = array.dim1;
  this->dim2 = array.dim2;
  this->dim3 = array.dim3;
  this->dim = array.dim;
  this->dim23 = array.dim23;
  return *this;
}


void Double3DArray::Redefine(int dim1, int dim2, int dim3)
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete [] this->data;
  this->Flag.Initialize();
  this->dim1 = dim1;
  this->dim2 = dim2;
  this->dim3 = dim3;
  this->dim=dim1*dim2*dim3;
  this->dim23=dim2*dim3;
  this->data=new double[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}

void Double3DArray::Set(int i, int j, int k, double val)
{
  data[i*dim23+j*dim3+k]=val;
}

void Double3DArray::Add(int i, int j, int k, double val)
{
  data[i*dim23+j*dim3+k]+=val;
}

void Double3DArray::Multiply(int i, int j, int k, double val)
{
  data[i*dim23+j*dim3+k]*=val;
}

