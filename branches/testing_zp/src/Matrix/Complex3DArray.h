#ifndef COMPLEX_3D_ARRAY_H
#define COMPLEX_3D_ARRAY_H

#include "MathTools/Complex.h"
#include "GeneralTools/GarbageFlag.h"

class Complex3DArray
{
public:
  Complex3DArray();
  Complex3DArray(int dim1, int dim2, int dim3);
  Complex3DArray(const Complex3DArray& array, bool duplicateFlag=false);
  ~Complex3DArray();
  
  Complex3DArray& operator = (const Complex3DArray& array); 
  void Redefine(int dim1, int dim2, int dim3);
  Complex& operator () (int i, int j, int k);
  Complex Get(int i, int j, int k);

  // access a whole vector of elements, where the first three indices coincide:
  Complex* GetVector(int i, int j);
  void Set(int i, int j, int k, const Complex &val);
  void Set(int i, int j, int k, double val);
  void Add(int i, int j, int k, const Complex &val);
  void Add(int i, int j, int k, double val);
  void Multiply(int i, int j, int k, const Complex &val);
  void Multiply(int i, int j, int k, double val);

private:
  int dim1, dim2, dim3;
  int dim23, dim;
  Complex* data;
  GarbageFlag Flag;
};


inline Complex& Complex3DArray::operator () (int i, int j, int k)
{
  int pos = i*dim23+j*dim2+k;
  return data[pos];
}

inline Complex Complex3DArray::Get(int i, int j, int k)
{
  if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2) && (k>-1) && (k<dim3))
    {
      int pos = i*dim23+j*dim2+k;
      return data[pos];
    }
  else return Complex();
}


inline Complex* Complex3DArray::GetVector(int i, int j)
{
  if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2) )
    {
      int pos = i*dim23+j*dim2;
      return &(data[pos]);
    }
  else return NULL;
}


#endif //COMPLEX_3D_ARRAY_H
