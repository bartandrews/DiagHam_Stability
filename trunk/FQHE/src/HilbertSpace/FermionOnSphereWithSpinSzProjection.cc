////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of fermions on sphere                       //
//                                                                            //
//                        last modification : 24/06/2002                      //
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


#include "HilbertSpace/FermionOnSphereWithSpinSzProjection.h"

#include <iostream>
#include <cstdlib>

using std::cout;
using std::endl;

FermionOnSphereWithSpinSzProjection::FermionOnSphereWithSpinSzProjection()
{
  this->FullSpace=0;
  this->HilbertSpaceDimension=0;
}

  
// basic construct^or
// ^
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// szProjectionValue = -1 / +1
// memory = amount of memory granted for precalculations

FermionOnSphereWithSpinSzProjection::FermionOnSphereWithSpinSzProjection (int nbrFermions, int totalLz, int lzMax, int totalSpin, int szProjectionValue, unsigned long memory)
{
  this->TargetSpace = this;
  if (szProjectionValue==0)
    {
      std::cout << "Spin Projection cannot be zero, please select +/- 1" << std::endl;
      exit(1);
    }
  this->SzProjectionValue = szProjectionValue;
  this->FullSpace = new FermionOnSphereWithSpin(nbrFermions, totalLz, lzMax, totalSpin, memory);
  this->HilbertSpaceDimension = FullSpace->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  this->Flag.Initialize();
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
FermionOnSphereWithSpinSzProjection::FermionOnSphereWithSpinSzProjection(const FermionOnSphereWithSpinSzProjection& fermions)
{
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->Flag=fermions.Flag;
  this->FullSpace = fermions.FullSpace;
  this->SzProjectionValue = fermions.SzProjectionValue;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//
FermionOnSphereWithSpinSzProjection::~FermionOnSphereWithSpinSzProjection ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete this->FullSpace;
    }
}
  

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space
FermionOnSphereWithSpinSzProjection& FermionOnSphereWithSpinSzProjection::operator = (const FermionOnSphereWithSpinSzProjection& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete this->FullSpace;
    }
  this->Flag=fermions.Flag;
  this->FullSpace = fermions.FullSpace;
  this->SzProjectionValue = fermions.SzProjectionValue;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  return *this;
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space
AbstractHilbertSpace* FermionOnSphereWithSpinSzProjection::Clone()
{
  return new FermionOnSphereWithSpinSzProjection(*this);
}



// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number
List<AbstractQuantumNumber*> FermionOnSphereWithSpinSzProjection::GetQuantumNumbers ()
{
  return this->FullSpace->GetQuantumNumbers();
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number
AbstractQuantumNumber* FermionOnSphereWithSpinSzProjection::GetQuantumNumber (int index)
{
  return this->FullSpace->GetQuantumNumber(index);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace
AbstractHilbertSpace* FermionOnSphereWithSpinSzProjection::ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter)
{
  return 0;
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension
int FermionOnSphereWithSpinSzProjection::GetHilbertSpaceDimension()
{
  return FullSpace->GetHilbertSpaceDimension();
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereWithSpinSzProjection::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereWithSpinSzProjection*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereWithSpinSzProjection::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}



// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinSzProjection::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  if (SzProjectionValue>0)
    return FullSpace->AduAduAuAu(index, m1, m2, n1, n2, coefficient);
  else
    return FullSpace->AddAddAdAd(index, m1, m2, n1, n2, coefficient);
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 
double FermionOnSphereWithSpinSzProjection::AA (int index, int n1, int n2)
{
  if (SzProjectionValue>0)
    return FullSpace->AuAu(index, n1, n2);
  else
    return FullSpace->AdAd(index, n1, n2);
}
  


// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereWithSpinSzProjection::AdAd (int m1, int m2, double& coefficient)
{
  if (SzProjectionValue>0)
    return FullSpace->AduAdu( m1, m2, coefficient);
  else
    return FullSpace->AddAdd( m1, m2, coefficient);
}
  

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnSphereWithSpinSzProjection::AdA (int index, int m)
{
  if (SzProjectionValue>0)
    return FullSpace->AduAu( index, m);
  else
    return FullSpace->AddAd( index, m);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnSphereWithSpinSzProjection::AdA (long index, int m)
{
  if (SzProjectionValue>0)
    return FullSpace->AduAu( index, m);
  else
    return FullSpace->AddAd( index, m);
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& FermionOnSphereWithSpinSzProjection::PrintState (ostream& Str, int state)
{
  return FullSpace->PrintState (Str, state);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex FermionOnSphereWithSpinSzProjection::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent)
{
  return Complex();
}
  
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used
void FermionOnSphereWithSpinSzProjection::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
