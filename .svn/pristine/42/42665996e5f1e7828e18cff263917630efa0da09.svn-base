////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2013 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class Spin operator for particle with spin                 //
//                                                                            //
//                        last modification : 29/10/2013                      //
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
#include "Operator/ParticleOnSphereSpinOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// coordinate = Coordinate index : 0 for x, 1 for y and 2 for z
// lzMax = maximum single particle momentum

ParticleOnSphereSpinOperator::ParticleOnSphereSpinOperator(ParticleOnSphereWithSpin* particle, int coordinate, int lzMax)
{
  this->Particle= particle;
  if((coordinate<0)||(coordinate>2))
  {
    cout << "Error: ParticleOnSphereSpinOperator coordinate index must be 0 (for x), 1 (y) or 2 (z). Aborting." << endl;	 
    exit(1);
  }
  this->CoordinateIndex=coordinate;
  this->LMax=lzMax;
}

// destructor
//

ParticleOnSphereSpinOperator::~ParticleOnSphereSpinOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereSpinOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereSpinOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereSpinOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereSpinOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
 
// evaluate the SxSquare matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSpinOperator::MatrixElementSquare (RealVector& V1, RealVector& V2)
{
  return this->PartialMatrixElementSquare(V1,V2,0l, V2.GetLargeVectorDimension());
}

// evaluate part of the matrix element, within a given range of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereSpinOperator::PartialMatrixElementSquare (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  double Element = 0.0; double Sum; double lim= this->LMax;
  double Coefficient=0.0;
  int TmpIndex;
  RealVector TmpState(Dim, true); // the second argument sets every coeff to zero
  RealVector TmpState2(Dim, true);

  switch(this->CoordinateIndex)
  {
//================================== Sx case =========================================================
    case 0:
      for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int m1 = 0; m1 <= lim; ++m1) 
	{
	  int Index=this->Particle->AduAd(i,m1,m1,Coefficient);
	  if ( (Index<Dim) && (Coefficient != 0.0))
	    TmpState[Index] += V2[i]*Coefficient;
	  Index=this->Particle->AddAu(i,m1,m1,Coefficient);
	  if ( (Index<Dim) && (Coefficient != 0.0))
	    TmpState[Index] += V2[i]*Coefficient;
	}
      }
      
      for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int m1 = 0; m1 <= lim; ++m1) 
	{
	  int Index=this->Particle->AduAd(i,m1,m1,Coefficient);
	  if ( (Index<Dim) && (Coefficient != 0.0))
	    TmpState2[Index] += TmpState[i]*Coefficient;
	  Index=this->Particle->AddAu(i,m1,m1,Coefficient);
	  if ( (Index<Dim) && (Coefficient != 0.0))
	    TmpState2[Index] += TmpState[i]*Coefficient;
	}
      }
      for (int i = (int) firstComponent; i < Dim; ++i)
	Element += V1[i]*TmpState2[i];

      break;
//================================= Sy case ==========================================================>
    case 1:
      for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int m1 = 0; m1 <= lim; ++m1) 
	  {
	    int Index=this->Particle->AduAd(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState[Index] += V2[i]*Coefficient;
	    Index=this->Particle->AddAu(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState[Index] -= V2[i]*Coefficient;
	  }
	}

	for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int m1 = 0; m1 <= lim; ++m1) 
	  {
	    int Index=this->Particle->AduAd(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState2[Index] += TmpState[i]*Coefficient;
	    Index=this->Particle->AddAu(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState2[Index] -= TmpState[i]*Coefficient;
	  }
	}

	for (int i = (int) firstComponent; i < Dim; ++i)
	  Element += V1[i]*TmpState2[i];

	Element *= -1.;

	break;	
//================================= Sz case ==========================================================>
    case 2:
      for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int m1 = 0; m1 <= lim; ++m1) 
	  {
	    int Index=this->Particle->AduAu(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState[Index] += V2[i]*Coefficient;
	    Index=this->Particle->AddAd(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState[Index] -= V2[i]*Coefficient;
	  }
	}

	for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int m1 = 0; m1 <= lim; ++m1) 
	  {
	    int Index=this->Particle->AduAu(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState2[Index] += TmpState[i]*Coefficient;
	    Index=this->Particle->AddAd(i,m1,m1,Coefficient);
	    if ( (Index<Dim) && (Coefficient != 0.0))
	      TmpState2[Index] -= TmpState[i]*Coefficient;
	  }
	}

	for (int i = (int) firstComponent; i < Dim; ++i)
	  Element += V1[i]*TmpState2[i];

	break;
  }

  return 0.25*Complex(Element);
} 

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereSpinOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  double Element = 0.0; double Sum; double lim= this->LMax;
  double coeff, coeffBis;
  int StateIndex, StateIndexBis;

  switch(this->CoordinateIndex)
  {
//================================== Sx case =========================================================
    case 0:
      for (int i = (int) firstComponent; i < Dim; ++i)
      {
	Sum = 0;
	for (int m = 0; m <= lim; ++m) 
	{
	  StateIndex = this->Particle->AduAd(i, m, m, coeff);
	  StateIndexBis = this->Particle->AddAu(i, m, m, coeffBis);
	  Sum += V2[StateIndex]+V2[StateIndexBis]; 
	}
	Element += V1[i] * Sum;
      }
      break;
//================================= Sy case ==========================================================>
    case 1:
      for (int i = (int) firstComponent; i < Dim; ++i)
      {
	Sum = 0;
	for (int m = 0; m <= lim; ++m) 
	{
	  StateIndex = this->Particle->AduAd(i, m, m, coeff);
	  StateIndexBis = this->Particle->AddAu(i, m, m, coeffBis);
	  Sum += -V2[StateIndex]+V2[StateIndexBis]; 
	}
	Element += V1[i] * Sum;
      }
      break;
//================================= Sz case ==========================================================>
    case 2:
      for (int i = (int) firstComponent; i < Dim; ++i)
	{
	   Sum = 0;
	   for (int m = 0; m <= lim; ++m) Sum += (this->Particle->AduAu(i,m)-this->Particle->AddAd(i,m));
	   Element += V2[i] * V1[i] * Sum ;      
	}
      break;
  }

  return 0.5*Complex(Element);
}


// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereSpinOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									  int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  double coeff=0., coeffBis=0.;
  int StateIndex, StateIndexBis;

  switch(this->CoordinateIndex)
  {
//================================== Sx case =========================================================
    case 0:
      for (int i = firstComponent; i < Last; ++i)
	{
	  for (int m = 0; m <= this->LMax; ++m) 
	  {
	    StateIndex = this->Particle->AduAd(i, m, m, coeff);
	    StateIndexBis = this->Particle->AddAu(i, m, m, coeffBis);
	    vDestination[StateIndex] += vSource[i] * coeff;      
	    vDestination[StateIndexBis] += vSource[i] * coeffBis;      
	  }
	}
      break;

//================================== Sy case =========================================================
    case 1:
      for (int i = firstComponent; i < Last; ++i)
	{
	  for (int m = 0; m <= this->LMax; ++m) 
	  {
	    StateIndex = this->Particle->AduAd(i, m, m, coeff);
	    StateIndexBis = this->Particle->AddAu(i, m, m, coeffBis);
	    vDestination[StateIndex] += vSource[i] * coeff;      
	    vDestination[StateIndexBis] -= vSource[i] * coeffBis;      
	  }
	}
      break;

//================================== Sz case =========================================================
    case 2:
      for (int i = firstComponent; i < Last; ++i)
	{
	  int Index = this->Particle->SzToMinusSz(i, coeff);
	  vDestination[Index] += vSource[i] * coeff;      
	}
      break;
  }

  return vDestination;
}



