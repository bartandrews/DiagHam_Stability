////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class density-density operator for particle with spin          //
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
#include "Operator/ParticleOnSphereWithSpinAllSzDensityOddChannel.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <stdlib.h>


using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationMomentumIndex2 = momentum index of the leftmost creation operator (from 0 to 2S)
// annihilationMomentumIndex1 = momentum index of the leftmost annihilation operator (from 0 to 2S)
// annihilationMomentumIndex2 = momentum index of the rightmost annihilation operator(from 0 to 2S)
ParticleOnSphereWithSpinAllSzDensityOddChannel::ParticleOnSphereWithSpinAllSzDensityOddChannel(FermionOnSphereWithSpinAllSz* particle, int creationMomentumIndex1, int annihilationMomentumIndex1)
{
  this->Particle= particle;
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->AnnihilationMomentumIndex1 = annihilationMomentumIndex1;
}

// destructor
//

ParticleOnSphereWithSpinAllSzDensityOddChannel::~ParticleOnSphereWithSpinAllSzDensityOddChannel()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinAllSzDensityOddChannel::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinAllSzDensityOddChannel::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnSphereWithSpinAllSz*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinAllSzDensityOddChannel::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinAllSzDensityOddChannel::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinAllSzDensityOddChannel::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int Index;
  double Element = 0.0;
  int m = this->AnnihilationMomentumIndex1;
  if (m != this->CreationMomentumIndex1)
   {
    cout << "error: difference in indices " << endl;
   }

  for (int i = 0; i < Dim; ++i)
   {
     Coefficient = this->Particle->AduAu(i, m);
     Element += 0.5 * V1[i] * V2[i] * Coefficient;
      
     Coefficient = this->Particle->AddAd(i, m);
     Element += 0.5 * V1[i] * V2[i] * Coefficient;

     Index = this->Particle->AduAd(i, m, m, Coefficient); 	
     if (Index < Dim)
      {
	Element += (-0.5) * V1[Index] * V2[i] * Coefficient; 
      }

     Index = this->Particle->AddAu(i, m, m, Coefficient); 	
     if (Index < Dim)
      {
	Element += (-0.5) * V1[Index] * V2[i] * Coefficient; 
      }

   }

  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinAllSzDensityOddChannel::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinAllSzDensityOddChannel::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;
  int Index;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int Sign;
  int m = this->AnnihilationMomentumIndex1;

   if (m != this->CreationMomentumIndex1)
    cout << "error: difference in indices" << endl;

   for (int i = 0; i < Dim; ++i)
   {
     Coefficient = this->Particle->AduAu(i, m);
     vDestination[i] += 0.5 * vSource[i] * Coefficient;
      
     Coefficient = this->Particle->AddAd(i, m);
     vDestination[i] += 0.5 * vSource[i] * Coefficient;

     Index = this->Particle->AduAd(i, m, m, Coefficient); 	
     if (Index < Dim)
      {
	vDestination[Index] += (-0.5) * vSource[i] * Coefficient; 
      }

     Index = this->Particle->AddAu(i, m, m, Coefficient); 	
     if (Index < Dim)
      {
	vDestination[Index] += (-0.5) * vSource[i] * Coefficient; 
      }

   }

 return vDestination;

}

