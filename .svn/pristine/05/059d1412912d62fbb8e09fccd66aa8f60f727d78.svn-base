#include "ComplexLapackDeterminant.h"

#ifdef HAVE_LAPACK

#include <iostream>
using std::cout;
using std::endl;

// binding to the LAPACK zgetrf routine for LU decomposition
//
extern "C" void FORTRAN_NAME(zgetrf)(const int* dimensionM, const int* dimensionN, const doublecomplex* matrixA,
				     const int* leadingDimensionA, const int *ipiv, const int *info);



// default constructor
//
ComplexLapackDeterminant::ComplexLapackDeterminant()
{
  this->Dimension=0;
  this->Components=0;
}
   

// constructor for an empty matrix
//
// dimension = number of columns and rows
// zero = tue if matrix elements have to be set to zero
ComplexLapackDeterminant::ComplexLapackDeterminant(int dimension, bool zero)
{
  this->Dimension=dimension;
  this->Components=new doublecomplex[dimension*dimension];
  this->Permutation = new int[dimension];
  if (zero)
    for (int i=0;i<dimension*dimension;++i) {Components[i].r=0.0;Components[i].i=0.0;}
}

// destructor
ComplexLapackDeterminant::~ComplexLapackDeterminant()
{
  if (this->Dimension>0)
    {
      delete [] this->Components;
      delete [] this->Permutation;
    }
}

// set a matrix element
//
// i = line position
// j = column position
// x = new value for matrix element
void ComplexLapackDeterminant::SetMatrixElement(int i, int j, const Complex& x)
{
  int index=i+j*Dimension;
  this->Components[index].r = Real(x);
  this->Components[index].i = Imag(x);
}

// set a matrix element
//
// i = line position
// j = column position
// (re,im) = new value for matrix element
void ComplexLapackDeterminant::SetMatrixElement(int i, int j, const double re, const double im)
{
  int index=i+j*Dimension;
  this->Components[index].r = re;
  this->Components[index].i = im;
}

// calculate the determinant, loosing information about the matrix elements
Complex ComplexLapackDeterminant::Determinant()
{
  int Information = 0;
  FORTRAN_NAME(zgetrf)(&Dimension, &Dimension, Components, &Dimension , Permutation, &Information);
  if (Information < 0)
    {
      cout << "Illegal argument " << -Information << " in LAPACK function call in ComplexLapackDeterminant.cc, line "<< __LINE__<<endl;
      exit(1);
    }

  int sign=0;
  Complex Result(1.0,0.0);
  
  for (int i=0; i<Dimension; ++i)
    {
      if (Permutation[i]!=i+1) sign ^= 1;
      Result *= Complex(Components[i+Dimension*i].r,Components[i+Dimension*i].i);
    }
  if (sign & 1) Result*=-1.0;

  return Result;
}

// resize determinant dimensions
void ComplexLapackDeterminant::Resize(int newDimension)
{
  if ((newDimension > 0) && (newDimension <= this->Dimension))
    {
      this->Dimension=newDimension;
    }
  else
    {
      if (Dimension > 0)
	{
	  delete [] this->Components;
	  delete [] this->Permutation;
	}
      if (newDimension>0)
	{
	  this->Dimension=newDimension;
	  this->Components=new doublecomplex[Dimension*Dimension];
	  this->Permutation = new int[Dimension];
	}
      else
	this->Dimension=0;	
    }
}

#endif
