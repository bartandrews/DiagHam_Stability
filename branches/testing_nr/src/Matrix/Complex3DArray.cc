#include "Complex3DArray.h"
// #include <iostream>
// using std::cout;
// using std::endl;

Complex3DArray::Complex3DArray()
{
  this->dim1=0;
  this->dim2=0;
  this->dim3=0;
  this->dim23=0;
  this->dim=0;
  this->Flag.Initialize();
  this->data=NULL;
}


Complex3DArray::Complex3DArray(int dim1, int dim2, int dim3)
{
  this->dim1=dim1;
  this->dim2=dim2;
  this->dim3=dim3;
  this->dim=dim1*dim2*dim3;
  this->dim23=dim2*dim3;
  this->Flag.Initialize();
  data = new Complex[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}

Complex3DArray::~Complex3DArray()
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete [] this->data;
}

// copy constructor
//
// array = array to copy
// duplicateFlag = true if datas have to be duplicated
Complex3DArray::Complex3DArray(const Complex3DArray& array, bool duplicateFlag) 
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
	this->data = new Complex [dim];
	for (int i = 0; i < this->dim; i++)
	  this->data[i] = array.data[i];
      }
}

// assignment
//
// array = array to assign
// return value = reference on current array
Complex3DArray& Complex3DArray::operator = (const Complex3DArray& array) 
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


void Complex3DArray::Redefine(int dim1, int dim2, int dim3)
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
  data = new Complex[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}


void Complex3DArray::Set(int i, int j, int k, const Complex &val)
{
  int pos = i*dim23+j*dim3+k;
  data[pos]=val;
}

void Complex3DArray::Set(int i, int j, int k, double val)
{
  int pos = i*dim23+j*dim3+k;
  data[pos]=val;
}

void Complex3DArray::Add(int i, int j, int k, const Complex &val)
{
  int pos = i*dim23+j*dim3+k;
  data[pos]+=val;
}

void Complex3DArray::Multiply(int i, int j, int k, const Complex &val)
{
  int pos = i*dim23+j*dim3+k;
  data[pos]*=val;
}

void Complex3DArray::Multiply(int i, int j, int k, double val)
{
  int pos = i*dim23+j*dim3+k;
  data[pos]*=val;
}

