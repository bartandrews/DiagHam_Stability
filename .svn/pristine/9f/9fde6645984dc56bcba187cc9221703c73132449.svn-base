////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of particle on sphere with spin                  //
//                                                                            //
//                        last modification : 12/12/2005                      //
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
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"


#include <iostream>
using std::cout;
using std::endl;

// virtual destructor
//

ParticleOnSphereWithSpin::~ParticleOnSphereWithSpin ()
{
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space
void ParticleOnSphereWithSpin::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->SetTargetSpace((ParticleOnSphereWithSpin*) targetSpace);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space
void ParticleOnSphereWithSpin::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension
int ParticleOnSphereWithSpin::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}
  

// apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_d a^+_m2_u a_n1_d a_n2_u operator to a given state (with m1+m2=n1+n2, one spin up an one spin own)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin up)
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// m = Lz value of particle to be added
// s = spin index of particle to be added (0=down, 1=up)
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long ParticleOnSphereWithSpin::Ad (unsigned long state, int m, int s, double& coefficient)
{
  cout << "Attention: calling placeholder function ParticleOnSphereWithSpin::Ad - please override in inherited class!" <<endl;
  return 0x0l;
}


// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AdA (int index, int m)
{
  return (this->AddAd(index, m) + this->AduAu(index, m));
}

// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AdA (long index, int m)
{
  return (this->AddAd(index, m) + this->AduAu(index, m));
}


// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::AuAu (int index, int n1, int n2)
{
  return 0.0;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin:: AdAd (int index, int n1, int n2)
{
  return 0.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::AuAd (int index, int n1, int n2)
{
  return 0.0;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAu (int index, int m, int n, double& coefficient)
{
  if (m == n)
    {
      coefficient = this->AduAu(index, m);
      return index;
    }
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAd (int index, int m, int n, double& coefficient)
{
  if (m == n)
    {
      coefficient = this->AddAd(index, m);
      return index;
    }
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d a_n_d operator to a given state  and check if resulting state belongs to Hilbert space
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, double& coefficient)
{
  return this->AduAd(index, m, m, coefficient);
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AduAu(index, m, m, coefficient);
}
  
// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AddAd(index, m, m, coefficient);
}
  

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AduAd(index, m, m, coefficient, nbrTranslation);
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AduAd(index, m, n, coefficient);
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AduAu(index, m, m, coefficient, nbrTranslationX);
}
  
// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AddAd(index, m, m, coefficient, nbrTranslationX);
}
  

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AduAd(index, m, m, coefficient, nbrTranslationX);
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AduAd(index, m, n, coefficient, nbrTranslationX);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAu (int index, int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AddAu(index, m, m, coefficient, nbrTranslationX);
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AddAu(index, m, n, coefficient, nbrTranslationX);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int ParticleOnSphereWithSpin::AddAu (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAu (int index, int m, double& coefficient)
{
  return this->AddAu(index, m, m, coefficient);
}


// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAu (int index, int m, double& coefficient, int& nbrTranslation)
{
  return this->AddAu(index, m, m, coefficient, nbrTranslation);
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AddAu(index, m, n, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AduAdu (m1, m2, coefficient);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AddAdd (m1, m2, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AduAdd (m1, m2, coefficient);
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AduAdu (m1, m2, coefficient, nbrTranslationX);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AddAdd (m1, m2, coefficient, nbrTranslationX);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AduAdd (m1, m2, coefficient, nbrTranslationX);
}

// apply a^+_m1_u a^+_m2_u operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdu (int index, int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}
   
// apply a^+_m1_d a^+_m2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdd (int index, int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}
  
// apply a^+_m1_u a^+_m2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdd (int index, int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a_n1_u a_n2_d operator to a state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AuAd (int index, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}
  
  
// apply a^+_m1_u a^+_m2_u operator to the state, assuming a different target space
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_d a^+_m2_d operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_d operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AduAdd (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_d a^+_m2_u operator to the state, assuming a different target space 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AddAdu (int index, int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->HilbertSpaceDimension;
}

  

// apply a_n1_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Adu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::Au (int index, int n1)
{
  return 0.0;
}

// apply a_n1_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next Adu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::Ad (int index, int n1)
{
  return 0.0; 
}

// apply a^+_m1_u operator to the state produced using Au method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::Adu (int m1, double& coefficient)
{
  return this->HilbertSpaceDimension; 
}

// apply a^+_m1_d operator to the state produced using Au method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::Add (int m1, double& coefficient)
{
  return this->HilbertSpaceDimension; 
}

// apply a^+_m1_d operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored	
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::Add (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  return this->Add (m, coefficient);
}

// apply a^+_m1_d operator to the state produced using Au method (without destroying it)
//
// m = first index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored	
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::Adu (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationX = 0;
  nbrTranslationY = 0;
  return this->Add (m, coefficient);
}

// apply a^+_m_u  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphereWithSpin::Adu (int index, int m, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphereWithSpin::Add (int index, int m, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a_m_u  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphereWithSpin::Au (int index, int m, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a_m_d  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// m = index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphereWithSpin::Ad (int index, int m, double& coefficient)
{
  return this->HilbertSpaceDimension;
}
  
// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  if (sigma1 == 0)
    {
      if (sigma2 == 0)
	return this->AuAu(index, n1, n2);
      else
	return this->AuAd(index, n1, n2);
    }
  if (sigma2 == 1)
    return this->AdAd(index, n1, n2);
  else
    {
      // warning : can be statistic sensitive
      return this->AuAd(index, n2, n1);
    }
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(2) index for the first creation operator
// sigma2 = SU(2) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  if (sigma1 == 0)
    {
      if (sigma2 == 0)
	return this->AduAdu(m1, m2, coefficient);
      else
	return this->AduAdd(m1, m2, coefficient);
    }
  if (sigma2 == 1)
    return this->AddAdd(m1, m2, coefficient);
  else
    {
      // warning : can be statistic sensitive
      return this->AduAdd(m2, m1, coefficient);
    }
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  return 0.0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnSphereWithSpin::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  return 0.0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  cout << "warning using dummy method ParticleOnSphereWithSpin::ProdAd" << endl;
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  cout << "warning using dummy method ParticleOnSphereWithSpin::ProdAd" << endl;
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning using dummy method ParticleOnSphereWithSpin::ProdAd" << endl;
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphereWithSpin::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning using dummy method ParticleOnSphereWithSpin::ProdAd" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AddAd (int index, int m)
{
  return 0.0;
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AduAu (int index, int m)
{
  return 0.0;
}

// apply a^+_m_s a_m_s operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AdsigmaAsigma (int index, int m, int sigma)
{
  return 0.0;
}

// apply a^+_m_s a_m_s operator to a given state)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// sigma = internal degree of freedom label of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSpin::AdsigmaAsigma (long index, int m, int sigma)
{
  return 0.0;
}

// apply a^+_m1_s1 a_m2_s2 operator to a given state
//
// index = index of the state on which the operator has to be applied
// m1 = index of the creation operator
// sigma1 = internal degree of freedom label of the creation operator
// m2 = index of the annihilation operator
// sigma2 = internal degree of freedom label of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int ParticleOnSphereWithSpin::AdsigmaAsigma (int index, int m1, int sigma1, int m2, int sigma2, double& coefficient)
{
  cout << "warning using dummy method ParticleOnSphereWithSpin::AdsigmaAsigma" << endl;
  return this->HilbertSpaceDimension;
}

// carefully test whether state is in Hilbert-space and find corresponding state index
//
// stateDescription = unsigned integer describing the state
// highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
// return value = corresponding index, or dimension of space, if not found
int ParticleOnSphereWithSpin::CarefulFindStateIndex(unsigned long stateDescription, int highestBit)
{
  return this->HilbertSpaceDimension;
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex ParticleOnSphereWithSpin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSpin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
									 AbstractFunctionBasis& basis, 
									 int nextCoordinates, int firstComponent, 
									 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnSphereWithSpin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
// Evaluate the Density Matrix of the spin up fermions in a sector with a fixed lzUp 
//
// lzUp = twice total momentum of up fermions.
// groundstate = reference on the total system groundstate
// return value = density matrix of the subsystem of spins up fermions.

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixSpinSeparation (int lzUp, RealVector & groundstate)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealMatrix ParticleOnSphereWithSpin::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{
  RealMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrFermionSector, int lzSector, int szSector, RealVector& groundState)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrNUpSector = number of spin up  that belong to the subsytem 
// nbrNDownSector = number of spin down  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, int nbrNUpSector, int nbrNDownSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}



// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm   
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix ParticleOnSphereWithSpin::EvaluatePartialEntanglementMatrixParticlePartition (int nbrFermionSector, int lzSector, int szSector, RealVector& groundState, 
											 bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  RealMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz,Sz=0, Sz parity sectors.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated. It should be equal to zero
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix ParticleOnSphereWithSpin::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, int szParity, RealVector& groundState, 
											 bool removeBinomialCoefficient, AbstractArchitecture* architecture)
{
  RealMatrix TmpMatrix;
  return TmpMatrix;
}
   

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrNUpSector = number of spin up  that belong to the subsytem 
// nbrNDownSector = number of spin down  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector,   ComplexVector& groundState, AbstractArchitecture* architecture)
{
  HermitianMatrix TmpMatrix;
  return TmpMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphereWithSpin::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, int szSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix)
{
  return entanglementMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSpin::EvaluatePartialDensityMatrixGenericRealSpacePartition (int nbrFermionSector, int lzSector, int szSector, int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
												     int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealVector& groundState, 
												     AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpMatrix;
  return TmpMatrix;
}


// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphereWithSpin::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, int szSector, 
															 int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
															 int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix)
{
  return entanglementMatrix;
}
  
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsystem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// szParitySector = parity sector for the discrete symmetry Sz<->-Sz. Can be either -1 or +1
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphereWithSpin::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector, int szParity,
															 int nbrOrbitalA, double* weightOrbitalAUp, 
															 double* weightOrbitalADown, 
															 int nbrOrbitalB, double* weightOrbitalBUp, 
															 double* weightOrbitalBDown, 
															 RealMatrix& entanglementMatrix)
{
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition breaking the momentum conservation. 
// The entanglement matrix is computed from precalculated particle entanglement matrices in each momentum sector
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// nbrConnectedOrbitalAUp = number of orbitals connected to a given one by the A part real space cut (for the spin up)
// nbrConnectedOrbitalADown = number of orbitals connected to a given one by the A part real space cut (for the spin down)
// connectedOrbitalAUp = orbitals taht connected to a given one by the A part real space cut (for the spin up)
// connectedOrbitalADown = orbitals taht connected to a given one by the A part real space cut (for the spin down)
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// nbrConnectedOrbitalBUp = number of orbitals connected to a given one by the B part real space cut (for the spin up)
// nbrConnectedOrbitalBDown = number of orbitals connected to a given one by the B part real space cut (for the spin down)
// connectedOrbitalBUp = orbitals taht connected to a given one by the B part real space cut (for the spin up)
// connectedOrbitalBDown = orbitals taht connected to a given one by the B part real space cut (for the spin down)
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// nbrEntanglementMatrices = number of available entanglement matrices with a fixed moementum
// entanglementMatrixLzSectors = momentum sector of each entanglement matrix
// entanglementMatrices = array containing the entanglement matrices with a fixed moementum
// return value = real space entanglement matrix

RealMatrix ParticleOnSphereWithSpin::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int szSector,
															int nbrOrbitalA, int* nbrConnectedOrbitalAUp, int* nbrConnectedOrbitalADown,
															int** connectedOrbitalAUp, int** connectedOrbitalADown, 
															double** weightOrbitalAUp, double** weightOrbitalADown, 
															int nbrOrbitalB, int* nbrConnectedOrbitalBUp, int* nbrConnectedOrbitalBDown, 
															int** connectedOrbitalBUp, int** connectedOrbitalBDown, 
															double** weightOrbitalBUp, double** weightOrbitalBDown, 
															int nbrEntanglementMatrices, int* entanglementMatrixLzSectors,
															RealMatrix* entanglementMatrices)
{
  RealMatrix TmpEntanglementMatrix;
  return TmpEntanglementMatrix;
}

// flip all spins of a given state
// 
// index = index of the state on which the operator has to be applied
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSpin::SzToMinusSz (int index, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// convert a given state from a generic basis to the current Sz subspace basis
//
// state = reference on the vector to convert
// basis = reference on the basis associated to state
// return value = converted vector

ComplexVector ParticleOnSphereWithSpin::ConvertFromNbodyBasis(ComplexVector& state, ParticleOnSphereWithSpin& basis)
{
  cout << "using dummy ParticleOnSphereWithSpin::ConvertFromNbodyBasis" << endl;
  ComplexVector TmpVector;
  return TmpVector;
}


// get the total spin
//
//return value: total spin of the Hilbert space
int ParticleOnSphereWithSpin::GetTotalSpin()
{
  cout << "using dummy ParticleOnSphereWithSpin::GetTotalSpin" << endl;
  return 0; 
}

// convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

RealVector ParticleOnSphereWithSpin::SU2ToSU2AllSz(RealVector& state, ParticleOnSphereWithSpin* su2space)
{
  cout << "using dummy ParticleOnSphereWithSpin::SU2ToSU2AllSz" << endl;
  ComplexVector Dummy;
  return Dummy;
}

// convert state of a SU(2) Hilbert space with fixed Sz to a SU(2) space with all sz sectors
//
// state = state that needs to be projected
// su2space = SU(2) space with fixed sz of the input state
// return value = input state expression in the SU(2) basis

ComplexVector ParticleOnSphereWithSpin::SU2ToSU2AllSz(ComplexVector& state, ParticleOnSphereWithSpin* su2space)
{
  cout << "using dummy ParticleOnSphereWithSpin::SU2ToSU2AllSz" << endl;
  ComplexVector Dummy;
  return Dummy;
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void ParticleOnSphereWithSpin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  cout << "using dummy ParticleOnSphereWithSpin::TransformOneBodyBasis" << endl;
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void ParticleOnSphereWithSpin::TransformOneBodyBasis(RealVector& initialState, RealVector& targetState, RealMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  cout << "using dummy ParticleOnSphereWithSpin::TransformOneBodyBasis" << endl;
}

// symmetrized a product of two decoupled states 
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

RealVector ParticleOnSphereWithSpin::SymmetrizeSU2SU2State (RealVector& leftVector, RealVector& rightVector, ParticleOnSphereWithSpin* leftSpace, 
							    ParticleOnSphereWithSpin* rightSpace, bool unnormalizedBasisFlag, 
							    AbstractArchitecture* architecture)
{
  RealVector TmpVector(this->GetHilbertSpaceDimension(), true);
  FQHESphereSymmetrizeU1U1StateOperation TmpOperation(this, leftSpace, rightSpace, &TmpVector, &leftVector, &rightVector, unnormalizedBasisFlag);
  TmpOperation.ApplyOperation(architecture);
//   this->SymmetrizeSU2SU2StateCore(TmpVector, leftVector, rightVector, leftSpace, rightSpace, unnormalizedBasisFlag, 
// 				  0, leftSpace->GetHilbertSpaceDimension());
  TmpVector.Normalize();
  return TmpVector;
}
  
// symmetrized a product of two decoupled states, core part
//
// outputVector = reference on the vector which will contain the symmetrized state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// firstComponent = index of the first component
// nbrComponents = number of components to symmetrize
// return value = symmetrized state

void ParticleOnSphereWithSpin::SymmetrizeSU2SU2StateCore (RealVector& symmetrizedVector, RealVector& leftVector, RealVector& rightVector, 
							  ParticleOnSphereWithSpin* leftSpace, ParticleOnSphereWithSpin* rightSpace, 
							  bool unnormalizedBasisFlag, unsigned long firstComponent, unsigned long nbrComponents)
{
  cout << "warning, using dummy method ParticleOnSphereWithSpin::SymmetrizeSU2SU2StateCore" << endl;
}


// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector ParticleOnSphereWithSpin::ForgeSU2FromU1(RealVector& upState, ParticleOnSphere* upStateSpace, RealVector& downState, ParticleOnSphere* downStateSpace)
{
  cout << "warning, using dummy method ParticleOnSphereWithSpin::ForgeSU2FromU1" << endl;
  return RealVector();
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

ComplexVector ParticleOnSphereWithSpin::ForgeSU2FromU1(ComplexVector& upState, ParticleOnSphere* upStateSpace, ComplexVector& downState, ParticleOnSphere* downStateSpace)
{
  cout << "warning, using dummy method ParticleOnSphereWithSpin::ForgeSU2FromU1" << endl;
  return ComplexVector();
}

