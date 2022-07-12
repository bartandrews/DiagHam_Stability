#include "NSphereParameters.h"

#include <cmath>

using std::cos;
using std::sin;

// default constructor 
//
NSphereParameters::NSphereParameters()
{
  this->Dimension=0;
  this->IsComplex=false;
  this->NbrParameters=0;
}

// constructor from an integer
//
// dim = dimension of sphere
// isComplex = flag indicating whether complex N-Sphere is used
NSphereParameters::NSphereParameters(int dim, bool isComplex)
{
  this->Dimension=dim;
  this->IsComplex=isComplex;
  if (IsComplex)
    {
      this->NbrParameters=2*dim-2;
      this->RealCoordinates.Resize(dim);
      this->ComplexCoordinates.Resize(dim);
    }
  else
    {
      this->NbrParameters=dim-1;
      this->RealCoordinates.Resize(dim);
    }
  this->Parameters.Resize(NbrParameters);
  this->Parameters.ClearVector();
  this->CosTable = new double[NbrParameters];
  this->SinTable = new double[NbrParameters];
  for (int i=0; i<NbrParameters; ++i)
    {
      this->CosTable[i]=1.0;
      this->SinTable[i]=0.0;
    }
}

// destructor
//
NSphereParameters::~NSphereParameters()
{
  if (this->Dimension!=0)
    {
      delete [] CosTable;
      delete [] SinTable;
    }
}

// assignement
//
// factorial = factorial coefficient to assign
// return value = reference on current factorial coefficient
NSphereParameters& NSphereParameters::operator = (const NSphereParameters& sphere)
{
  if (this->NbrParameters!=sphere.NbrParameters)
    {
      if (this->Dimension!=0)
	{
	  delete [] CosTable;
	  delete [] SinTable;
	}
      this->CosTable = new double[sphere.NbrParameters];
      this->SinTable = new double[sphere.NbrParameters];
    }
  this->Dimension=sphere.Dimension;
  this->IsComplex=sphere.IsComplex;
  this->NbrParameters=sphere.NbrParameters;
  this->ComplexCoordinates=ComplexVector(sphere.ComplexCoordinates,true);
  this->RealCoordinates= RealVector(sphere.RealCoordinates,true);
  this->Parameters=RealVector(sphere.Parameters,true);
  for (int i=0; i<NbrParameters; ++i)
    {
      this->CosTable[i] = sphere.CosTable[i];
      this->SinTable[i] = sphere.SinTable[i];
    }
  return *this;
}

// set parameters
void NSphereParameters::SetParameters(double *parameters)
{
  for (int i=0; i<NbrParameters; ++i)
    {
      if (this->Parameters[i]!=parameters[i])
	{
	  this->CosTable[i]=cos(parameters[i]);
	  this->SinTable[i]=sin(parameters[i]);
	}
    }
  for (int i=0; i<Dimension-1; ++i)
    this->RealCoordinates[i]=CosTable[i];
  this->RealCoordinates[Dimension-1]=1.0;
  for (int i=1; i<Dimension; ++i)
    for (int j=0; j<i; ++j)
      this->RealCoordinates[i]*=SinTable[j];
  if (this->IsComplex)
    {
      this->ComplexCoordinates[0]=this->RealCoordinates[0];
      for (int i=0; i<Dimension-1; ++i)
	{
	  this->ComplexCoordinates[1+i].Re=this->RealCoordinates[1+i]*CosTable[Dimension-1+i];
	  this->ComplexCoordinates[1+i].Im=this->RealCoordinates[1+i]*SinTable[Dimension-1+i];
	}
    }
}
  
