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
#include "Operator/ParticleOnSphereWithSpinAllSzExcitonOrder.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex1 = index of the leftmost creation operator
// creationIndex2 = index of the rightmost creation operator
// annihilationIndex1 = index of the leftmost annihilation operator
// annihilationIndex2 = index of the rightmost annihilation operator

ParticleOnSphereWithSpinAllSzExcitonOrder::ParticleOnSphereWithSpinAllSzExcitonOrder(FermionOnSphereWithSpinAllSz* particle, int creationIndex1, int creationIndex2, int annihilationIndex1, int annihilationIndex2)
{
  this->Particle= particle;
  this->CreationIndex1 = creationIndex1;
  this->CreationIndex2 = creationIndex2;
  this->AnnihilationIndex1 = annihilationIndex1;
  this->AnnihilationIndex2 = annihilationIndex2;
}

// destructor
//

ParticleOnSphereWithSpinAllSzExcitonOrder::~ParticleOnSphereWithSpinAllSzExcitonOrder()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinAllSzExcitonOrder::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinAllSzExcitonOrder::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnSphereWithSpinAllSz*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinAllSzExcitonOrder::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinAllSzExcitonOrder::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
 
// Get exciton order parameter <c^*_up,n c^*_down,0 c_up,o c_down,n>
//
// vSource = groundstate vector
// firstIndex = index of the first exciton orbital
// secondIndex = index of the second exciton orbital
double ParticleOnSphereWithSpinAllSzExcitonOrder::GetExcitonOrderLevel2(RealVector& vSource)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();

  if (vSource.GetVectorDimension() != Dim)
   {
     cout << "Dimension mismatch "<<endl;
     exit(1);
   }
  RealVector FinalState(Dim, true);
  double TmpCoefficient;

  for (int i = 0; i < Dim; ++i)
   {
     double Coefficient = vSource[i] * this->Particle->AuAd(i, this->AnnihilationIndex1, this->AnnihilationIndex2);
     if (Coefficient != 0.0)
      {
	int Index = this->Particle->AduAdd(this->CreationIndex1, this->CreationIndex2, TmpCoefficient);
	if (Index < Dim)
	 FinalState[Index] += Coefficient * TmpCoefficient;
      }
    }
  return (FinalState * vSource);
}
