////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere density operator               //
//                                                                            //
//                        last modification : 11/12/2002                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;



// default constructor
//

ParticleOnSphereDensityOperator::ParticleOnSphereDensityOperator()
{
}

// constructor from default datas
//
// particle = hilbert space associated to the particles
// index = index of the density operator

ParticleOnSphereDensityOperator::ParticleOnSphereDensityOperator(ParticleOnSphere* particle, int index)
{
  this->Particle= (ParticleOnSphere*) (particle->Clone());
  this->OperatorIndex = index;
  this->OperatorIndexDagger = index;
}

// constructor when dealing with two different Hilbert spaces
//
// particle = hilbert space associated to the right hand state (target space has to be fixed to the hilbert space associated to the left hand state)
// indexDagger = index of the creation operator that is part of the density operator
// index = index of the annihilation operator that is part of the density operator
 
ParticleOnSphereDensityOperator::ParticleOnSphereDensityOperator(ParticleOnSphere* particle, int indexDagger, int index)
{
  this->Particle= (ParticleOnSphere*) (particle->Clone());
  this->OperatorIndexDagger = indexDagger;
  this->OperatorIndex = index;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnSphereDensityOperator::ParticleOnSphereDensityOperator(ParticleOnSphereDensityOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->OperatorIndexDagger = oper.OperatorIndexDagger;
  this->OperatorIndex = oper.OperatorIndex;
}

// destructor
//

ParticleOnSphereDensityOperator::~ParticleOnSphereDensityOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereDensityOperator::Clone ()
{
  return new ParticleOnSphereDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  double Element = 0.0;
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      Element += V1[i] * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += V1[TmpIndex] * V2[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      Element += V1[i] * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += V1[TmpIndex] * V2[i] * TmpCoefficient;
	    }
	}
     }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								 int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] +=  vSource[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] +=  vSource[i] * TmpCoefficient;
	    }
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnSphereDensityOperator::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									 int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient;
		    }
		}
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient;
		    }
		}
	    }
	}
    }
  return vDestinations;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  Complex Element = 0.0;
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      Element += Conj(V1[i]) * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  int FullDim = this->Particle->GetTargetHilbertSpaceDimension();
	  for (int i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if ((TmpIndex != FullDim) && (TmpCoefficient != 0.0))
		Element += Conj(V1[TmpIndex]) * V2[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      Element += Conj(V1[i]) * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Dim; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		Element += Conj(V1[TmpIndex]) * V2[i] * TmpCoefficient;
	    }
	}
     }
  return Element;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  int Last = firstComponent + nbrComponent;;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] += vSource[i] * TmpCoefficient;
	    }
	}
    }
  else
    {
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      vDestination[i] += vSource[i] * this->Particle->AdA(i, this->OperatorIndex);
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  long Last = ((long) firstComponent) + ((long) nbrComponent);
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		vDestination[TmpIndex] += vSource[i] * TmpCoefficient;
	    }
	}
    }
  return vDestination;
}
  
// multiply a set of vectors by the current operator for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnSphereDensityOperator::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									    int firstComponent, int nbrComponent)
{
  if (((int) this->Particle->GetLargeHilbertSpaceDimension()) == this->Particle->GetHilbertSpaceDimension())
    {
      int Last = firstComponent + nbrComponent;;
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  int TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (int i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient;
		    }
		}
	    }
	}
    }
  else
    {
      long Last = ((long) firstComponent) + ((long) nbrComponent);
      if (this->OperatorIndexDagger == this->OperatorIndex)
	{
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      double TmpCoefficient = this->Particle->AdA(i, this->OperatorIndex);
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  vDestinations[k][i] += vSources[k][i] * TmpCoefficient;
		}
	    }
	}
      else
	{
	  long TmpIndex;
	  double TmpCoefficient = 0.0;
	  for (long i = firstComponent; i < Last; ++i)
	    {
	      TmpIndex =  this->Particle->AdA(i, this->OperatorIndexDagger, this->OperatorIndex, TmpCoefficient);
	      if (TmpCoefficient != 0.0)
		{
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][TmpIndex] +=  vSources[k][i] * TmpCoefficient;
		    }
		}
	    }
	}
    }
  return vDestinations;
}



