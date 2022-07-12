#ifndef DOUBLE_3D_ARRAY_H
#define DOUBLE_3D_ARRAY_H

#include "GeneralTools/GarbageFlag.h"

class Double3DArray
{
public:
  Double3DArray();
  Double3DArray(int dim1, int dim2, int dim3);
  Double3DArray(const Double3DArray& array, bool duplicateFlag=false);
  ~Double3DArray();

  Double3DArray& operator = (const Double3DArray& array); 

  void Redefine(int dim1, int dim2, int dim3);
  double& operator () (int i, int j, int k);
  double Get(int i, int j, int k);
  double* GetVector(int i, int j);
  void Set(int i, int j, int k, double val);
  void Add(int i, int j, int k, double val);
  void Multiply(int i, int j, int k, double val);

private:
  int dim1, dim2, dim3;
  int dim, dim23;
  double *data;
  GarbageFlag Flag;
};


inline double& Double3DArray::operator () (int i, int j, int k)
{
  return (data[i*dim23+j*dim3+k]);
}

inline double Double3DArray::Get(int i, int j, int k)
{
  if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2) && (k>-1) && (k<dim3))
    return (data[i*dim23+j*dim3+k]);
  else return 0.0;
}

inline double* Double3DArray::GetVector(int i, int j)
{
   if ((i>-1) && (i<dim1) && (j>-1) && (j<dim2))
     return &(data[i*dim23+j*dim3]);
   else return 0;
}


#endif //DOUBLE_3D_ARRAY_H
