#ifndef COMPLEX_4D_ARRAY_H
#define COMPLEX_4D_ARRAY_H

#include "MathTools/Complex.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/GarbageFlag.h"

class Complex4DArray
{
public:
  Complex4DArray();
  Complex4DArray(int dim1, int dim2, int dim3, int dim4);
  Complex4DArray(const Complex4DArray& array, bool duplicateFlag=false);
  ~Complex4DArray();

  Complex4DArray& operator = (const Complex4DArray& array); 
  
  // redefine the array, loosing all contents
  void Redefine(int dim1, int dim2, int dim3, int dim4);

  // access a matrix element
  Complex& operator () (int i, int j, int k, int l);

  // access with index-within-bounds testing
  Complex Get(int i, int j, int k, int l);

  // access a whole vector of elements, where the first three indices coincide:
  Complex* GetVector(int i, int j, int k);
  
  void Set(int i, int j, int k, int l, const Complex &val);
  void Set(int i, int j, int k, int l, double val);
  void Add(int i, int j, int k, int l, const Complex &val);
  void Add(int i, int j, int k, int l, double val);
  void Multiply(int i, int j, int k, int l, const Complex &val);
  void Multiply(int i, int j, int k, int l, double val);

private:
  int dim1, dim2, dim3, dim4;
  int dim234, dim34, dim;
  Complex* data;
  GarbageFlag Flag;
};

 
inline Complex& Complex4DArray::operator () (int i, int j, int k, int l)
{
  int pos = i*dim234+j*dim34+k*dim4+l;
  return data[pos];
}

inline Complex Complex4DArray::Get(int i, int j, int k, int l)
{
  if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2) && (k>-1) && (k<dim3) && (l>-1) && (l<dim4))
    {
      int pos = i*dim234+j*dim34+k*dim4+l;
      return data[pos];
    }
  else return Complex();
}

inline Complex* Complex4DArray::GetVector(int i, int j, int k)
{
  if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2) && (k>-1) && (k<dim3) )
    {
      int pos = i*dim234+j*dim34+k*dim4  ;
      return &(data[pos]);
    }
  else return NULL;
}



#endif //COMPLEX_4D_ARRAY_H
