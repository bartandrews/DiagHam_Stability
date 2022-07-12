////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on sphere                        //
//                                                                            //
//                        last modification : 15/07/2002                      //
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
#include "HilbertSpace/ParticleOnSphere.h"
#include "MathTools/LongRational.h"
#include "Vector/LongRationalVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/ComplexVector.h"

#include <iostream>

using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnSphere::~ParticleOnSphere ()
{
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void ParticleOnSphere::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  printf("Attention - trying to call SetTargetSpace on generic class ParticleOnSphere\n");
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int ParticleOnSphere::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}


// get the number of particles in the target space
//
// return value = number of particles in the target space
int ParticleOnSphere::GetTargetNbrParticles()
{
  return 0;
}  


//return dimension of the subspace of the target space whose elements are related by the tz<->-tz and Z3 symmetry
//
//return value = dimension of the subspace
int ParticleOnSphere::GetSymmetryDimension(int i)
{
 return this->GetSymmetryDimension(i); 
}


// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int ParticleOnSphere::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnSphere::NoSymmetry;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2), safe version i.e. works with any numbers of particles
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAASafe (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->AdAdAA(index, m1, m2, n1, n2, coefficient);
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

int ParticleOnSphere::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  double Coefficient = this->AA(index, n1, n2);
  int Index = this->AdAd(m1, m2, coefficient);
  coefficient *= Coefficient;
  return Index;
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

long ParticleOnSphere::AdAdAA (long index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return this->LargeHilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 a^+_m3 a_n1 a_n2 a_n3 operator to a given state (with m1+m2+m3=n1+n2+n3)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// m3 = third index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// n3 = third index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAdAAA (int index, int m1, int m2, int m3, int n1, int n2, int n3, double& coefficient)
{
  int TmpM[3];
  int TmpN[3];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  return ProdAdProdA(index, TmpM, TmpN, 3, coefficient);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 a^+_m3 a^+_m4 a_n1 a_n2 a_n3 a_n4 operator to a given state (with m1+m2+m3+m4=n1+n2+n3+n4)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// m3 = third index for creation operator
// m4 = fourth index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// n3 = third index for annihilation operator
// n4 = fourth index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAdAdAdAAAA (int index, int m1, int m2, int m3, int m4, int n1, int n2, int n3, int n4, double& coefficient)
{
  int TmpM[4];
  int TmpN[4];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpM[3] = m4;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  TmpN[3] = n4;
  return ProdAdProdA(index, TmpM, TmpN, 4, coefficient);
}

// apply creation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// m = Lz value of particle to be added
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long ParticleOnSphere::Ad (unsigned long state, int m, double& coefficient)
{
  cout << "Attention: calling placeholder function ParticleOnSphereWithSpin::Ad - please override in inherited class!" <<endl;
  return 0x0l;
}


// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphere::AA (int index, int n1, int n2)
{
  return 0.0;
}

// apply a_n1 a_n2 operator to a given state without keeping it in cache
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AA (int index, int n1, int n2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}


// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnSphere::ProdA (int index, int* n, int nbrIndices)
{
  return 0.0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
// use double when calculating normalization factors to avoid overflow
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnSphere::ProdAL (int index, int* n, int nbrIndices)
{
  return 0.0;
}

// apply Prod_i a_ni operator to a given state
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdA (int index, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAd (int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 operator to the state 
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdAd (int index, int m1, int m2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphere::AdAd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AdAd(m1, m2, coefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AdAd(m1, m2, coefficient, nbrTranslationX);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->ProdAd(m, nbrIndices, coefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->ProdAd(m, nbrIndices, coefficient, nbrTranslationX);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
// use double when calculating normalization factors to avoid overflow
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAdL (int* m, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_ni operator to a given state
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::ProdAd (int index, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphere::AdA (int index, int m)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::AdA (int index, int m, int n, double& coefficient)
{
  if (m != n)
    return this->HilbertSpaceDimension;
  else
    {
      coefficient = this->AdA(index, m);
      return index;
    }
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphere::AdA (long index, int m)
{
  return this->LargeHilbertSpaceDimension;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long ParticleOnSphere::AdA (long index, int m, int n, double& coefficient)
{
  if (m != n)
    return this->HilbertSpaceDimension;
  else
    {
      coefficient = this->AdA(index, m);
      return index;
    }
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long ParticleOnSphere::AdA (long index, int m, int n, Complex& coefficient)
{
  return this->AdA(index, m, n, coefficient.Re);
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::AdA (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  nbrTranslation = 0;
  return this->AdA(index, m, n, coefficient);
}
    
// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  nbrTranslationY = 0;
  return this->AdA(index, m, n, coefficient, nbrTranslationX);
}
 
// apply a_n  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphere::A (int index, int nn)
{
  cout << "Warning: using the default operator ParticleOnSphere::A" << endl;
  return 0.0;
}

// apply a^+_m  operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad or A call
//
// index = index of the state on which the operator has to be applied
// m = index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphere::Ad (int index, int m)
{
  cout << "Warning: using the default operator ParticleOnSphere::Ad" << endl;
  return 0.0;
}

// apply a_n operator to the state produced using the A or Ad method (without destroying it)
//
// n = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::A (int n, double& coefficient)
{
  cout << "Warning: using the default operator ParticleOnSphere::A" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m operator to the state produced using the A or Ad method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphere::Ad (int m, double& coefficient)
{
  cout << "Warning: using the default operator ParticleOnSphere::Ad" << endl;
  return this->HilbertSpaceDimension;
}
 

// apply a^+_m operator to the state produced using AuAu method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::Ad (int m, double& coefficient, int& nbrTranslationX)
{
  cout << "using int ParticleOnSphere::Ad (int m, double& coefficient, int& nbrTranslationX)"<< endl;
  return this->HilbertSpaceDimension;
}

 
// apply a^+_m operator to the state produced using AuAu method (without destroying it)
//
// m = first index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int ParticleOnSphere::Ad (int m, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  return this->HilbertSpaceDimension;
}

// apply a_n  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// n = index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphere::A (int index, int n, double& coefficient)
{
  cout << "Warning: using the default operator ParticleOnSphere::A" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_n  operator to a given state. 
//
// index = index of the state on which the operator has to be applied
// n = index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value =  index of the resulting state 

int ParticleOnSphere::Ad (int index, int n, double& coefficient)
{
  cout << "Warning: using the default operator ParticleOnSphere::Ad" << endl;
  return this->HilbertSpaceDimension;
}

// check whether HilbertSpace implements ordering of operators
//

bool ParticleOnSphere::HaveOrder ()
{
  return false;
}

// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value
int ParticleOnSphere::CheckOrder (int* m, int* n, int nbrIndices)
{
  return 0;
}


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex ParticleOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
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
Complex ParticleOnSphere::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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

Complex ParticleOnSphere::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, 
								 int nextCoordinates, int firstComponent, 
								 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnSphere::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, RealVector& groundState)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

HermitianMatrix ParticleOnSphere::EvaluatePartialDensityMatrix (int subsytemSize, int nbrFermionSector, int lzSector, ComplexVector& groundState)
{
  HermitianMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem

RealMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{ 
  RealMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = entanglement matrix of the subsytem
LongRationalMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrix (int subsytemSize, int nbrFermionSector, int lzSector, LongRationalVector& groundState)
{
  LongRationalMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;
}
  

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnSphere::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, ComplexVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  HermitianMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									  RealVector& groundState, RealSymmetricMatrix* densityMatrix)
{
  return 0l;
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									  ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  return 0l;
}

// core part of the evaluation density matrix particle partition calculation involving a sum of projetors 
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// nbrGroundStates = number of projectors
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									  int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  HermitianMatrix TmpDensityMatrix(densityMatrix->GetNbrRow());
  long Tmp = 0l;
  for (int i = 0; i < nbrGroundStates; ++i)
    {
      TmpDensityMatrix.ClearMatrix();
      Tmp += this->EvaluatePartialDensityMatrixParticlePartitionCore(minIndex, nbrIndex, complementaryHilbertSpace, destinationHilbertSpace, groundStates[i], &TmpDensityMatrix);
      TmpDensityMatrix *= weights[i];
      (*densityMatrix) += TmpDensityMatrix;
    }
  return Tmp;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// perimeter = cylinder perimeter
// height = height of a cylinder (from -H/2 to H/2) 
// xcut = x-coordinate of a cylinder cut
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartitionCylinder (int nbrBosonSector, int lzSector, double perimeter, double height, double xcut, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixGenericRealSpacePartition (int nbrFermionSector, int lzSector, int nbrOrbitalA, double* weightOrbitalA, 
											     int nbrOrbitalB, double* weightOrbitalB, RealVector& groundState, 
											     AbstractArchitecture* architecture)
{
  cout << "warning, EvaluatePartialDensityMatrixGenericRealSpacePartition not implemented" << endl;
  RealSymmetricMatrix PartialDensityMatrix;
  return PartialDensityMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphere::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
														 int nbrOrbitalA, double* weightOrbitalA, 
														 int nbrOrbitalB, double* weightOrbitalB, RealMatrix& entanglementMatrix)
{
  cout << "warning, EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix not implemented" << endl;
  return entanglementMatrix;
}


// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = pointer to array of entanglement matrix (will be overwritten)
// nbrEntanglementMatrices = number of entanglement matrices to be computed
// return value = reference on the entanglement matrix
RealMatrix* ParticleOnSphere::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
														 int nbrOrbitalA, double* weightOrbitalA, 
														 int nbrOrbitalB, double* weightOrbitalB, RealMatrix* entanglementMatrix, int nbrEntanglementMatrices)
{
  cout << "warning, EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix not implemented" << endl;
  return entanglementMatrix;
}
													 
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

ComplexMatrix& ParticleOnSphere::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
														 int nbrOrbitalA, double* weightOrbitalA, 
														 int nbrOrbitalB, double* weightOrbitalB, ComplexMatrix& entanglementMatrix)
{
  cout << "warning, EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix not implemented" << endl;
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalA = weight of each orbital in the A part (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalB = weight of each orbital in the B part (starting from the leftmost orbital)
// entanglementMatrix = pointer to array of entanglement matrix (will be overwritten)
// nbrEntanglementMatrices = number of entanglement matrices to be computed
// return value = reference on the entanglement matrix

ComplexMatrix* ParticleOnSphere::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrFermionSector, int lzSector, 
														 int nbrOrbitalA, double* weightOrbitalA, 
														 int nbrOrbitalB, double* weightOrbitalB, ComplexMatrix* entanglementMatrix, int nbrEntanglementMatrices)
{
  cout << "warning, EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix not implemented" << endl;
  return entanglementMatrix;
}
// core part of the evaluation density matrix real space partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
// incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// return value = number of components that have been added to the density matrix

long ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									   RealVector& groundState,  RealSymmetricMatrix* densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange)
{
  return 0l;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  RealMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// core part of the entanglement matrix evaluation for the particle partition
// 
// minIndex = first index to consider in the complementary Hilbert space
// nbrIndex = number of indices to consider in the complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// entanglementMatrix = pointer to entanglement matrix
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = number of components that have been added to the entanglement matrix

long ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace, 
									       ParticleOnSphere* destinationHilbertSpace, RealVector& groundState, RealMatrix* entanglementMatrix, 
									       bool removeBinomialCoefficient)
{
  cout << "warning, using generic method ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartitionCore" << endl;
  return 0l;
}
   
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
// nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB,
										RealVector& groundState, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  RealMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix* ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  RealMatrix* PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
// nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
// groundState = reference on the total system ground state
// nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix* ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB,
										RealVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  RealMatrix* PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}


// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, ComplexVector& groundState, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  ComplexMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
// nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB,
										ComplexVector& groundState, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  ComplexMatrix PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix* ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, ComplexVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  ComplexMatrix* PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// The entanglement matrix is only evaluated in a given Lz sector and both A and B are resticted to a given number of orbitals
// 
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// nbrOrbitalA = number of orbitals that have to be kept for the A part (starting from the leftmost orbital)
// nbrOrbitalA = number of orbitals that have to be kept for the B part (starting from the rightmost orbital)
// groundState = reference on the total system ground state
// nbrEntanglementMatrices = number of vectors whose entanglement matrix has to be computed
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix* ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition(int nbrFermionSector, int lzSector, int nbrOrbitalA, int nbrOrbitalB,
										ComplexVector* groundState, int nbrEntanglementMatrices, bool removeBinomialCoefficient)
{
  cout << "calling non defined function ParticleOnSphere::EvaluatePartialEntanglementMatrixParticlePartition" << endl;
  ComplexMatrix* PartialEntanglementMatrix;
  return PartialEntanglementMatrix;  
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector
// and computed from precalculated entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// entanglementMatrix = reference on the entanglement matrix
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphere::EvaluatePartialDensityMatrixParticlePartitionFromEntanglementMatrix (int nbrBosonSector, int lzSector, RealMatrix& entanglementMatrix)
{
  RealSymmetricMatrix PartialDensityMatrix(entanglementMatrix.GetNbrRow(), true);
  for (int i = 0; i < entanglementMatrix.GetNbrRow(); ++i)
    {
      for (int j = i; j < entanglementMatrix.GetNbrRow(); ++j)
	{
	  for (int k = 0; k < entanglementMatrix.GetNbrColumn(); ++k)
	    {
	      PartialDensityMatrix.AddToMatrixElement(i, j, entanglementMatrix(i, k) * entanglementMatrix(j, k));
	    }
	}
    }
  return PartialDensityMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphere::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix)
{
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

ComplexMatrix& ParticleOnSphere::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, ComplexMatrix& entanglementMatrix)
{
  cout << "calling non defined function ParticleOnSphere::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix" << endl;
  return entanglementMatrix;
}

// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using real space partition on a cylinder. The entanglement matrix is only evaluated in a given Lz sector.
// and computed from precalculated particle entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// perimeter = cylinder perimeter
// height = height of a cylinder (from -H/2 to H/2) 
// xcut = x-coordinate of the cut   // entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& ParticleOnSphere::EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrixCylinder (int nbrBosonSector, int lzSector, double perimeter, double height, double xcut, RealMatrix& entanglementMatrix)
{
  return entanglementMatrix;
}

// evaluate coeffecicents requested to compute the real space partition
//
// lzMax = twice the maximum angular momentum
// thetaTop = inclination angle defining one edge of the cut in radians
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in radians
// incompleteBetaThetaTop = reference on the pointer to the array (allocation done by the method) where the top part coefficients will be stored
// incompleteBetaThetaBotton = reference on the pointer to the array (allocation done by the method) where the bottom part coefficients will be stored

void ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(int lzMax, double thetaTop, double thetaBottom, double*& incompleteBetaThetaTop, double*& incompleteBetaThetaBottom)
{
  double * LogFactorialsNphi = new double[lzMax +2];
  LogFactorialsNphi[0] = 0.0;
  LogFactorialsNphi[1] = 0.0;
  for (int i = 2 ; i < lzMax+2; ++i)
    LogFactorialsNphi[i] = LogFactorialsNphi[i - 1] + log((double) i); 
  
  
  double LogSinThetaTop = 0;
  
  if (thetaTop > 1e-10)
    LogSinThetaTop = 2.0*log(sin(thetaTop/2.));
  
  double LogCosThetaTop = 2.0*log(cos(thetaTop/2.));
  	
  double LogSinThetaBot = 2.0*log(sin(thetaBottom/2.));
  double LogCosThetaBot = 0;
  if (thetaBottom < (179 + 1e-10))
    LogCosThetaBot = 2.0*log(cos(thetaBottom/2.));
  
  // Compute the incomplete beta function for x=sin^2(thetaTop/2) and x=sin^2(thetaBottom/2)
  incompleteBetaThetaTop = new double[lzMax + 1];
  incompleteBetaThetaBottom = new double[lzMax + 1];
  if ( thetaTop < 1e-10)
    incompleteBetaThetaTop [lzMax] = 0.0;
  else
    incompleteBetaThetaTop [lzMax] = exp((lzMax + 1.0) * LogSinThetaTop);
  if( thetaBottom < (179 + 1e-10) )
    incompleteBetaThetaBottom [lzMax] =  exp((lzMax + 1.0) * LogSinThetaBot);
  else
    incompleteBetaThetaBottom [lzMax] = 1.0;  
  
  for (int i = lzMax - 1; i >=  0; --i)
    {
      if (thetaTop < 1e-10)
	incompleteBetaThetaTop[i] = 0.0;
      else
	incompleteBetaThetaTop[i] = exp(LogFactorialsNphi[lzMax+1] - LogFactorialsNphi[i+1] - LogFactorialsNphi[lzMax-i] + (i + 1.0)*LogSinThetaTop + (lzMax - i)*LogCosThetaTop) + incompleteBetaThetaTop[i+1] ;
      if (thetaBottom < (179 + 1e-10) )
	incompleteBetaThetaBottom[i] = exp(LogFactorialsNphi[lzMax+1] - LogFactorialsNphi[i+1] - LogFactorialsNphi[lzMax-i] + (i + 1.0)*LogSinThetaBot + (lzMax - i)*LogCosThetaBot) + incompleteBetaThetaBottom[i+1] ;
      else
	incompleteBetaThetaBottom[i] = 1.0;
    }
  delete [] LogFactorialsNphi;
}

// evaluate coeffecicents requested to compute the real space partition (cylinder geometry)
//
// lzMax = twice the maximum angular momentum
// perimeter = cylinder perimeter
// xcut = x-coordinate of the cut
// incompleteBetaThetaTop = reference on the pointer to the array (allocation done by the method) where the top part coefficients will be stored

void ParticleOnSphere::EvaluatePartialDensityMatrixRealSpacePartitionCoefficientCylinder(int lzMax, double perimeter, double xcut, double*& incompleteBetaThetaTop)
{
  incompleteBetaThetaTop = new double[lzMax + 1];
//  cout<<"Real space coefficients: "<<endl;
  for (int i = lzMax; i >=0; --i)
    {
        double Xm = (i - 0.5 * lzMax) * 2.0 * M_PI/perimeter;
	incompleteBetaThetaTop[i] = 0.5 * (1.0 + erf(xcut - Xm));
        //cout<<incompleteBetaThetaTop[i]<<" ";
    }
//  cout<<endl;
}

// compute part of the Schmidt decomposition, allowing cut in the reduced denisty matrix eigenvalue space
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrFermionSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
// rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
// diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
// transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
// return value = reference on rebuiltSchmidtGroundState

RealVector& ParticleOnSphere::EvaluatePartialSchmidtDecomposition(int subsytemSize, int nbrFermionSector, int lzSector, double eigenvalueCut,
								  RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
								  RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  return rebuiltSchmidtGroundState;
}

// compute part of the Schmidt decomposition for the particle partition, allowing cut in the reduced denisty matrix eigenvalue space
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// eigenvalueCut = discard all contribution from the reduced density matrix whose eigenvalues is lower than eigenvalueCut
// rebuiltSchmidtGroundState = reference on the state to whose current sector contribution to the Schmidt decomposition will be added 
// diagonalizedDensityMatrix = reference on the diagonalized reduced density matrix associated to the current sector (with down ordered diagonal elements)
// transformationMatrix =  reference on the transformation matric that diagonalizes the reduced density matrix
// return value = reference on rebuiltSchmidtGroundState

RealVector& ParticleOnSphere::EvaluatePartialSchmidtDecompositionParticlePartition(int nbrParticleSector, int lzSector, double eigenvalueCut,
										   RealVector& groundState, RealVector& rebuiltSchmidtGroundState,
										   RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  return rebuiltSchmidtGroundState;
}
  
// rebuild a state from its Schmidt decomposition for the particle partition
// 
// nbrParticleSector = number of particles that belong to the subsytem (i.e. part A)
// lzSector = Lz sector in which the density matrix has to be evaluated  (i.e. part A)
// schmidtDecomposedState = reference on the vector to which the rebuild state will be added
// nbrSingularValues = number of singular values (can be lower than the actual number of ingular values to perform a truncation)
// singularValues = array containing the singular values
// aVectors = matrix than contains the singular vectors of the part A
// bVectors = matrix than contains the singular vectors of the part B

void ParticleOnSphere::RebuildStateFromSchmidtDecompositionParticlePartition(int nbrParticleSector, int lzSector, RealVector& schmidtDecomposedState, 
									     int nbrSingularValues, double* singularValues, RealMatrix& aVectors, RealMatrix& bVectors)
{
}

// convert a state to its occupation number representation
//
// index = index of the state
// finalState = reference on the array where the occupation number representation has to be stored

void ParticleOnSphere::GetOccupationNumber(long index, unsigned long*& finalState)
{
}

// get the list of occupied orbitals in a given state
//
// state = ID of the state
// orbitals = list of orbitals to be filled

void ParticleOnSphere::GetOccupied(int state, int* orbitals)
{
}

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int ParticleOnSphere::FindStateIndex(char* stateDescription)
{
  return -1;
}

// find state index from an array of occupied orbitals
//
// stateDescription = array describing the state (stored as k1,k2,k3,...)
// return value = corresponding index, -1 if an error occured
int ParticleOnSphere::FindStateIndex(int* stateDescription)
{
  return -1;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component has to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& ParticleOnSphere::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& ParticleOnSphere::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  return state;
}

// convert a state such that its components are now expressed in the normalized basis, without applying the global normalization to the final state
//
// state = reference to the state to convert
// return value = converted state

RealVector& ParticleOnSphere::ConvertFromUnnormalizedMonomialNoGlobalNormalization(RealVector& state)
{
  cout << "warning, using dummy method ParticleOnSphere::ConvertFromUnnormalizedMonomialNoGlobalNormalization" << endl;
  return state;
}
 

// convert a state such that its components are now expressed in the normalized basis, shifting all orbitals
//
// state = reference to the state to convert
// shift = shift to apply to each orbitals
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& ParticleOnSphere::ShiftedConvertFromUnnormalizedMonomial(RealVector& state, int shift, long reference)
{
  return state;
}

// convert a state such that its components, given in the conformal limit,  are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& ParticleOnSphere::ConvertFromConformalLimit(RealVector& state, long reference)
{
  return this->ConvertFromUnnormalizedMonomial(state, reference, false);
}


// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ParticleOnSphere::PrintStateMonomial (ostream& Str, long state)
{
  return Str;
}

// print a given State using the monomial notation, with one column per particle (using space as a seperator)
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ParticleOnSphere::PrintColumnFormattedStateMonomial (ostream& Str, long state)
{
  return Str;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& ParticleOnSphere::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
					  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag, double coefficient)
{
  return outputVector;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

LongRationalVector& ParticleOnSphere::FuseStates (LongRationalVector& outputVector, LongRationalVector& leftVector, LongRationalVector& rightVector, int padding, 
						  ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace, bool symmetrizedFlag, LongRational& coefficient)
{
  return outputVector;
}

// fuse multiple states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& ParticleOnSphere::FuseMultipleStates (RealVector& outputVector, int nbrInputVectors, RealVector* inputVectors, int* paddings, 
						  ParticleOnSphere** inputSpaces, bool symmetrizedFlag, double coefficient)
{
  return outputVector;
}

// fuse multiple states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// nbrInputVectors = number of input vectors
// inputVectors = input vectors whose Hilbert space will be fuse from  left to right
// paddings = number of unoccupied one body states that have to be inserted between two consecutive fused spaces
// inputSpaces = point to the Hilbert space that will be fuse to the left
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

LongRationalVector& ParticleOnSphere::FuseMultipleStates (LongRationalVector& outputVector, int nbrInputVectors, LongRationalVector* inputVectors, int* paddings, 
							  ParticleOnSphere** inputSpaces, bool symmetrizedFlag, LongRational& coefficient)
{
  return outputVector;
}


// use product rule to produce part of the components of a system from a smaller one
//
// outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
// inputVector = reference on the vector associated to the smaller system
// inputSpace = pointer to the Hilbert space of the smaller system
// commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
// commonPatterSize = number of elements in the commonPattern array
// addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
// addedPatterSize = number of elements in the addedPattern array
// coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the product rule state

RealVector& ParticleOnSphere::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
					    int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
					    double coefficient, bool symmetrizedFlag)
{
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double ParticleOnSphere::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  return 0.0;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational ParticleOnSphere::JackSqrNormalization (LongRationalVector& outputVector, long minIndex, long nbrComponents)
{
  return LongRational();
}

// compute part of the Jack polynomial scalar product in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double ParticleOnSphere::JackScalarProduct (RealVector& state1, RealVector& state2, long minIndex, long nbrComponents)
{
  return 0.0;
}


// compute part of the Jack polynomial square normalization in a given range of indices
//
// state1 = reference on the first unnormalized Jack polynomial
// state2 = reference on the second unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational ParticleOnSphere::JackScalarProduct (LongRationalVector& state1, LongRationalVector& state2, long minIndex, long nbrComponents)
{
  return LongRational();
}
  

// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector ParticleOnSphere::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnSphere* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  RealVector TmpVector;
  return TmpVector;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals
bool ParticleOnSphere::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  return false;
}


// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the Lz component
int ParticleOnSphere::GetLzValue(int j)
{
  return 0;
}

// get Sz component of the spin
//
// j = index of the vector in Hilbert space
// return value = Sz component

int ParticleOnSphere::GetSzValue(int j)
{
  return 0;
} 

// transform a vector belonging to this vector space in the lz->-lz
//
// finalSpace = the space obtained after the lz->-lz operation
// initialVector = vector on which the operation will be apply
// return value = vector resulting of the operation

RealVector ParticleOnSphere::GetLzSymmetricVector(ParticleOnSphere* finalSpace, RealVector& initialVector)
{
  RealVector Tmp (this->LargeHilbertSpaceDimension, true);
  return Tmp;
}

// transform a vector belonging to this vector space in the lz->-lz
//
// finalSpace = the space obtained after the lz->-lz operation
// initialVector = vector on which the operation will be apply
// return value = vector resulting of the operation

LongRationalVector ParticleOnSphere::GetLzSymmetricVector(ParticleOnSphere* finalSpace, LongRationalVector& initialVector)
{
  LongRationalVector Tmp (this->LargeHilbertSpaceDimension, true);
  return Tmp;
}


// compute the number of particles in each Landau level
//
// state = ID of the state to handle
// lLOccupationConfiguration = array where the decomposition will be store

void ParticleOnSphere::LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration)
{
}

void ParticleOnSphere::EvaluatePartialDensityMatrixMultipartiteParticlePartition(ParticleOnSphere * spaceA, ParticleOnSphere * spaceB,ParticleOnSphere * spaceC,  RealVector groundstate, RealSymmetricMatrix* densityMatrix, AbstractArchitecture* architecture)
{
  cout <<"calling non defined function ParticleOnSphere::EvaluatePartialDensityMatrixMultipartiteParticlePartition"<<endl;
  return ;  
}

// print a given state using the most compact notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& ParticleOnSphere::PrintCompactState (ostream& Str, long state)
{
  return this->PrintState(Str, (int) state);
}

// convert the vector with a given Lz to the full space (all Lz components)
// inputState = input vector
// inputSpace = input Hilbert space with given Lz
// return value = vector in the full Hilbert space

void ParticleOnSphere::ConvertToAllLz (ComplexVector& inputState, ParticleOnSphere* inputSpace, ComplexVector& outputState)
{
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& ParticleOnSphere::NormalizeJackToCylinder(RealVector& state, double aspect)
{
  cout << "warning using dummy ParticleOnSphere::NormalizeJackToCylinder" << endl;
  return state;
}

// normalize from the cylinder geometry to the Jack normalization
//
// state = reference to the state to unnormalize
// aspect = cylinder aspect ratio
// reference = set which component as to be normalized to 1
// return value = unnormalized state

RealVector& ParticleOnSphere::NormalizeCylinderToJack(RealVector& state, double aspect, long reference)
{
  cout << "warning using dummy ParticleOnSphere::NormalizeCylinderToJack" << endl;
  return state;
}

// normalize Jack with respect to cylinder basis
//
// state = reference to the Jack state to normalize
// aspect = aspect ratio of cylinder
// return value = normalized state

RealVector& ParticleOnSphere::NormalizeSphereToCylinder(RealVector& state, double aspect)
{
  cout << "warning using dummy ParticleOnSphere::NormalizeSphereToCylinder" << endl;
  return state;
}


// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void ParticleOnSphere::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
						      long memory, long initialIndex, long nbrComponents)
{
  cout << "warning using dummy ParticleOnSphere::CreateStateFromMPSDescription" << endl;
}

// create a state from its MPS description, inclusing additional quasihole matrices
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// memory = amount of memory that can be use to precompute matrix multiplications  
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void  ParticleOnSphere::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholeBMatrices,
						      ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
						      long memory, long initialIndex, long nbrComponents)
{
  cout << "warning using dummy ParticleOnSphere::CreateStateFromMPSDescription" << endl;
}

// create a state from its site-dependent MPS description
//
// bMatrices = array that gives the site-dependent MPS
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void ParticleOnSphere::CreateStateFromSiteDependentMPSDescription (SparseRealMatrix** bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
								   long initialIndex, long nbrComponents)
{
  cout << "warning using dummy ParticleOnSphere::CreateStateFromSiteDependentMPSDescription" << endl;
}

// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnSphere::ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  ComplexVector TmpVector;
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnSphere::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  ComplexVector TmpVector;
  return TmpVector;
}

// apply a Gutzwiller projection (in the orbital space) to a given state
//
// state = reference on the state to project
// space = pointer to the Hilbert space where state is defined
// return value = Gutzwiller projected state

ComplexVector ParticleOnSphere::GutzwillerProjection(ComplexVector& state, ParticleOnSphere* space)
{
  ComplexVector TmpVector;
  return TmpVector;
}

// Compute the overlap of two states made from different one-body wavefunction
//
//  firstVector = reference on the first vector 
//  secondVector = reference on the second vector 
// overlapMatrix = pointer to the table with the overlap between the one-body states
Complex ParticleOnSphere::ComputeOverlapWaveFunctionsWithDifferentGamma (ComplexVector& firstVector, ComplexVector& secondVector, Complex * overlapMatrix)
{
  cout <<"using dummy function ParticleOnSphere::ComputeOverlapWaveFunctionsWithDifferentGamma"<<endl;
  return Complex();
}


// compute sum of positions in the x and y direction for lattice class
//
// index = index of the state in the basis whose position sums are to be computed
// positionX = reference on the sum of positions in the x direction
// positionY = reference on the sum of positions in the y direction
void ParticleOnSphere::GetPositionSum(int index,int & positionX, int & positionY)
{
  cout <<"using dummy function ParticleOnSphere::GetPositionSum"<<endl;
}

void ParticleOnSphere::GetPositionSum(unsigned long * monomial, int & positionX, int & positionY)
{  cout <<"using dummy function ParticleOnSphere::GetPositionSum"<<endl;}
