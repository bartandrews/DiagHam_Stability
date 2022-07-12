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

#include <iostream>
using std::cout;
using std::endl;

// virtual destructor
//

ParticleOnSphereWithSpin::~ParticleOnSphereWithSpin ()
{
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
