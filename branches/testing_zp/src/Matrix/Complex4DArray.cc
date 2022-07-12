#include "Complex4DArray.h"
// #include <cstdlib>
// #include <iostream>
// using std::cout;
// using std::endl;

Complex4DArray::Complex4DArray()
{
  this->dim1=0;
  this->dim2=0;
  this->dim3=0;
  this->dim4=0;
  this->dim234=0;
  this->dim34=0;
  this->dim=0;
  this->Flag.Initialize();
  this->data=NULL;
}


Complex4DArray::Complex4DArray(int dim1, int dim2, int dim3, int dim4)
{
  this->dim1=dim1;
  this->dim2=dim2;
  this->dim3=dim3;
  this->dim4=dim4;
  this->dim=dim1*dim2*dim3*dim4;
  this->dim234=dim2*dim3*dim4;
  this->dim34=dim3*dim4;
  this->Flag.Initialize();
  data = new Complex[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}


// destructor
//
Complex4DArray::~Complex4DArray()
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete [] this->data;
}


// copy constructor
//
// array = array to copy
// duplicateFlag = true if datas have to be duplicated
Complex4DArray::Complex4DArray(const Complex4DArray& array, bool duplicateFlag) 
{
  this->dim1 = array.dim1;
  this->dim2 = array.dim2;
  this->dim3 = array.dim3;
  this->dim4 = array.dim4;
  this->dim = array.dim;
  this->dim234 = array.dim234;
  this->dim34 = array.dim34;
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
Complex4DArray& Complex4DArray::operator = (const Complex4DArray& array) 
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete[] this->data;
  this->Flag = array.Flag;
  this->data = array.data;
  this->dim1 = array.dim1;
  this->dim2 = array.dim2;
  this->dim3 = array.dim3;
  this->dim4 = array.dim4;
  this->dim = array.dim;
  this->dim234 = array.dim234;
  this->dim34 = array.dim34;
  return *this;
}

void Complex4DArray::Redefine(int dim1, int dim2, int dim3, int dim4)
{
  if ((this->data != NULL) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      delete [] this->data;
  this->Flag.Initialize();
  this->dim1 = dim1;
  this->dim2 = dim2;
  this->dim3 = dim3;
  this->dim4 = dim4;
  this->dim=dim1*dim2*dim3*dim4;
  this->dim234=dim2*dim3*dim4;
  this->dim34=dim3*dim4;
  data = new Complex[dim];
  for(int i=0;i<dim;++i) data[i]=0.0;
}

void Complex4DArray::Set(int i, int j, int k, int l, const Complex &val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  data[pos]=val;
}

void Complex4DArray::Set(int i, int j, int k, int l, double val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;  
  data[pos]=val;
}

void Complex4DArray::Add(int i, int j, int k, int l, const Complex &val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  data[pos]+=(val);
}

void Complex4DArray::Add(int i, int j, int k, int l, double val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  data[pos]+=val;
}

void Complex4DArray::Multiply(int i, int j, int k, int l, const Complex &val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  data[pos]*=val;
}

void Complex4DArray::Multiply(int i, int j, int k, int l, double val)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  data[pos]*=val;
}
