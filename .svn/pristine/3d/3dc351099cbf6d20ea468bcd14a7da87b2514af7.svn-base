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
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int ParticleOnSphere::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int ParticleOnSphere::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnSphere::NoSymmetry;
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
  int TmpM[2];
  int TmpN[2];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpN[0] = n1;
  TmpN[1] = n2;
  return ProdAdProdA(index, TmpM, TmpN, 2, coefficient);
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

// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int ParticleOnSphere::FindStateIndex(char* stateDescription)
{
  return -1;
}

