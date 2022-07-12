////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on sphere density-density operator           //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex1 = index of the leftmost creation operator
// creationIndex2 = index of the rightmost creation operator
// annihilationIndex1 = index of the leftmost annihilation operator
// annihilationIndex2 = index of the rightmost annihilation operator

ParticleOnSphereDensityDensityOperator::ParticleOnSphereDensityDensityOperator(ParticleOnSphere* particle, int creationIndex1, int creationIndex2,
									       int annihilationIndex1, int annihilationIndex2)
{
  this->Particle= particle;
  this->CreationIndex1 = creationIndex1;
  this->CreationIndex2 = creationIndex2;
  this->AnnihilationIndex1 = annihilationIndex1;
  this->AnnihilationIndex2 = annihilationIndex2;
}

// destructor
//

ParticleOnSphereDensityDensityOperator::~ParticleOnSphereDensityDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereDensityDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereDensityDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereDensityDensityOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      Element += V1[Index] * V2[i] * Coefficient;
    }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereDensityDensityOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& ParticleOnSphereDensityDensityOperator::Multiply(RealVector& vSource, RealVector& vDestination)
{
  return this->Multiply(vSource, vDestination, 0, this->Particle->GetHilbertSpaceDimension());
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereDensityDensityOperator::Multiply(RealVector& vSource, RealVector& vDestination, 
							     int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      vDestination[Index] = vSource[i] * Coefficient;
    }
  return vDestination;
}
  
// multiply a vector by the current operator and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereDensityDensityOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return this->Multiply(vSource, vDestination, 0, this->Particle->GetHilbertSpaceDimension());
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereDensityDensityOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  return vDestination;
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnSphereDensityDensityOperator& O)
{
  RealVector TmpV2 (O.Particle->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.Particle->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.Particle->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < O.Particle->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i];
	  Str << "   ";
	}
      Str << endl;
    }
  return Str;
}
 
// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnSphereDensityDensityOperator& O)
{
  RealVector TmpV2 (O.Particle->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [O.Particle->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(O.Particle->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2[i - 1] = 0.0;
	}
      TmpV2[i] = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (O.Particle->GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (O.Particle->GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j][i];
	  Str << ",";
	}
      Str << TmpV[O.Particle->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (O.Particle->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][O.Particle->GetHilbertSpaceDimension() - 1];
      Str << ",";
    }
  Str << TmpV[O.Particle->GetHilbertSpaceDimension() - 1][O.Particle->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}
