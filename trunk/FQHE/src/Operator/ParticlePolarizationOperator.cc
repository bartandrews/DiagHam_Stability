////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of particle polarization operator                 //
//                                                                            //
//                        last modification : 16/09/2002                      //
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


#include "Operator/ParticlePolarizationOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// nbrParticle = number of particles

ParticlePolarizationOperator::ParticlePolarizationOperator(ParticleOnTorusWithSpin* particle, int nbrParticle)
{
  this->Particle = particle;
  this->NbrParticle = nbrParticle;
}

// destructor
//

ParticlePolarizationOperator::~ParticlePolarizationOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticlePolarizationOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticlePolarizationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnTorusWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticlePolarizationOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticlePolarizationOperator::GetHilbertSpaceDimension ()
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

Complex ParticlePolarizationOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  double Factor = 2.0 / ((double) this->NbrParticle);
  double Coefficient = 0.0;
  double Element = 0.0;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      this->Particle->SumAudAu(i, Coefficient);
      Element += V1[i] * V2[i] * (Factor * Coefficient - 1.0);
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

RealVector& ParticlePolarizationOperator::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double Factor = 2.0 / ((double) this->NbrParticle);
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      this->Particle->SumAudAu(i, Coefficient);
      vDestination[i] = vSource[i] * (Factor * Coefficient - 1.0);
    }
  return vDestination;
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticlePolarizationOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  double Factor = 2.0 / ((double) this->NbrParticle);
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      this->Particle->SumAudAu(i, Coefficient);
      vDestination[i] += vSource[i] * (Factor * Coefficient - 1.0);
    }
  return vDestination;
}							       

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticlePolarizationOperator& O)
{
  ComplexVector TmpV2 (O.Particle->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [O.Particle->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(O.Particle->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < O.Particle->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0.0)
	    {
	      Str << ((TmpV[j].Im(i))) << "i";
	    }
	  else
	    if (TmpV[j].Im(i) > 0.0)
	      {
		Str << "+" << TmpV[j].Im(i) << "i";
	      }
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

MathematicaOutput& operator << (MathematicaOutput& Str, ParticlePolarizationOperator& O)
{
  ComplexVector TmpV2 (O.Particle->GetHilbertSpaceDimension(), true);
  ComplexVector* TmpV = new ComplexVector [O.Particle->GetHilbertSpaceDimension()];
  for (int i = 0; i < O.Particle->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = ComplexVector(O.Particle->GetHilbertSpaceDimension());
      if (i > 0)
	{
	  TmpV2.Re(i - 1) = 0.0;
	}
      TmpV2.Re(i) = 1.0;
      O.Multiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (O.Particle->GetHilbertSpaceDimension() - 1); ++i)
    {
      Str << "{";
      for (int j = 0; j < (O.Particle->GetHilbertSpaceDimension() - 1); ++j)
	{
	  Str << TmpV[j].Re(i);
	  if (TmpV[j].Im(i) < 0)
	    {
	      Str << ((TmpV[j].Im(i))) << "I";
	    }
	  else
	    if (TmpV[j].Im(i) > 0)
	      {
		Str << "+" << ((TmpV[j].Im(i))) << "I";
	      }	  
	  Str << ",";
	}
      Str << TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Re(i);
      if (TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(i) < 0)
	{
	  Str << ((TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	}
      else
	if (TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(i) > 0)
	  {
	    Str << "+" << ((TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(i))) << "I";
	  }	  
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (O.Particle->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j].Re(O.Particle->GetHilbertSpaceDimension() - 1);
      if (TmpV[j].Im(O.Particle->GetHilbertSpaceDimension() - 1) < 0)
	{
	  Str << ((TmpV[j].Im(O.Particle->GetHilbertSpaceDimension() - 1))) << "I";
	}
      else
	if (TmpV[j].Im(O.Particle->GetHilbertSpaceDimension() - 1) > 0)
	  {
	    Str << "+" << ((TmpV[j].Im(O.Particle->GetHilbertSpaceDimension() - 1))) << "I";
	  }	  
      Str << ",";
    }
  Str << TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Re(O.Particle->GetHilbertSpaceDimension() - 1);
  if (TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(O.Particle->GetHilbertSpaceDimension() - 1) < 0)
    {
      Str << ((TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(O.Particle->GetHilbertSpaceDimension() - 1))) << "I";
    }
  else
    if (TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(O.Particle->GetHilbertSpaceDimension() - 1) > 0)
      {
	Str << "+" << ((TmpV[O.Particle->GetHilbertSpaceDimension() - 1].Im(O.Particle->GetHilbertSpaceDimension() - 1))) << "I";
      }	  
  Str << "}}";
  return Str;
}


