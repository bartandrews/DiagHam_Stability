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
#include "HilbertSpace/ParticleOnLattice.h"

#include <iostream>
using std::cout;
using std::endl;

// virtual destructor
//

ParticleOnLattice::~ParticleOnLattice ()
{
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void ParticleOnLattice::SetTargetSpace(ParticleOnLattice* targetSpace)
{
  cout << "Attention: ParticleOnLattice::SetTargetSpace not overloaded!"<<endl;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int ParticleOnLattice::GetTargetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// get information about any additional symmetry of the Hilbert space
//
// return value = symmetry id

int ParticleOnLattice::GetHilbertSpaceAdditionalSymmetry()
{
  return ParticleOnLattice::NoSymmetry;
}

// apply annihilation operator to a word, using the conventions
// for state-coding and quantum numbers of this space
// state = word to be acted upon
// q = quantum number of boson to be added
// coefficient = reference on the double where the multiplicative factor has to be stored
unsigned long ParticleOnLattice::A (unsigned long state, int q, double &coefficient)
{
  cout << "ParticleOnLattice::A"<<endl;
  return 0x0ul;
}

// apply a_n1 / a^\dagger_n1 operator to a given state and search in target space
//
// index = index of the state on which the operator has to be applied
// q = index for annihilation operator
// coefficient = prefactor
// return value =  index in target space
int ParticleOnLattice::A (int index, int q, double &coefficient)
{
  cout << "ParticleOnLattice::A"<<endl;
  return 0x0ul;
}

int ParticleOnLattice::Ad (int index, int q, double &coefficient)
{
  cout << "ParticleOnLattice::Ad"<<endl;
  return 0x0ul;
}


// apply a^+_q1 a^+_q2 a_r1 a_r2 operator to a given state (with q1+q2=r1+r2)
//
// index = index of the state on which the operator has to be applied
// q1 = first index for creation operator
// q2 = second index for creation operator
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnLattice::AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient)
{
  cout << "ParticleOnLattice::AdAdAA"<<endl;
  return this->HilbertSpaceDimension;
}


// apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// r1 = first index for annihilation operator
// r2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnLattice::AA (int index, int r1, int r2)
{
  return 0.0;
}


// apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
//
// q1 = first index for creation operator
// q2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnLattice::AdAd (int q1, int q2, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor   

double ParticleOnLattice::ProdA (int index, int* n, int nbrIndices)
{
  return 0.0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnLattice::ProdAd (int* m, int nbrIndices, double& coefficient)
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

int ParticleOnLattice::AdA (int index, int m, int n, double& coefficient)
{
  if (m != n)
    return this->HilbertSpaceDimension;
  else
    {
      coefficient = this->AdA(index, m);
      return index;
    }
}

// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// NbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// QValues = array of quantum numbers where an interaction is present
// InteractionPerQ = coefficient U_q of the interaction
//
double ParticleOnLattice::AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  double coefficient, result=0.0;
  for (int i=0; i<nbrInteraction; ++i)
    {
      this->AdAdAA(i,qValues[i],qValues[i],qValues[i],qValues[i],coefficient);
      result+=coefficient*interactionPerQ[i];
    }
  return result;
}

// calculate (possibly non-local) density-density interactions \sum q V_{q1,q2} : n_q1 n_q2 :
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum
// interactionPerQ12 = coefficient V_(q1, q2) of the interaction
// q12Values = array of quantum numbers of the orbitals in tuples (q1, q2), 2*nbrInteraction entries in total
//
double ParticleOnLattice::RhoRhoDiagonal(int index, int nbrInteraction, double *interactionPerQ12, int *q12Values)
{
  double coefficient, result=0.0;
  for (int i=0; i<nbrInteraction; ++i)
    {
      this->AdAdAA(index, q12Values[i<<1], q12Values[(i<<1)+1], q12Values[(i<<1)+1], q12Values[i<<1], coefficient);
      result+=coefficient*interactionPerQ12[i];
    }
  return result;
}


// apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
// index = index of the state on which the operator has to be applied
// nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
// qValues = array of quantum numbers where an interaction is present
// interactionPerQ = coefficient U_q of the interaction
//
double ParticleOnLattice::ProdAdProdADiagonal(int index,int nbrBody, int nbrInteraction, double *interactionPerQ, int *qValues)
{
  cout << "ParticleOnLattice::AdAdAA"<<endl;
  return this->HilbertSpaceDimension;
}

// check whether HilbertSpace implements ordering of operators
//
bool ParticleOnLattice::HaveOrder ()
{
  return false;
}

// check whether a given operator \prod c^\dagger_m \prod c_n increases or decreases the index of a state
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value = 1, if created state is of higher value, 0 if equal, and -1 if lesser value

int ParticleOnLattice::CheckOrder (int* m, int* n, int nbrIndices)
{
  return 0;
}



// obtain a list of quantum numbers in state
// quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles

void ParticleOnLattice::ListQuantumNumbers(int index, int *quantumNumbers)
{
  double tmp;
  this->ListQuantumNumbers(index, quantumNumbers, tmp);
}


// apply a gauge transformation
// phases = phases in array ordered according to the quantum number q
// input = vector that has to be transformed according to that gauge

ComplexVector& ParticleOnLattice::GaugeTransformVector(double *phases, ComplexVector& input)
{
  cout << "Ignoring gauge transform: method needs XXXOnLattice::GaugeTransformVector needs to be implemented!"<<endl;
  return input;
}

// conversion to generic (full) many-body representation in real-space basis
// state: many-body state in Ky-momentum basis
// nbodyBasis: full Hilbert-space in real-space representation
// returns: vector in many-body basis of targetSpace

ComplexVector& ParticleOnLattice::ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis)
{
  return this->ConvertToNbodyBasis(state, nbodyBasis, 0, this->GetHilbertSpaceDimension());
}

// conversion to generic (full) many-body representation in real-space basis
// state: many-body state in Ky-momentum basis
// nbodyBasis: full Hilbert-space in real-space representation
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// returns: vector in many-body basis of targetSpace

ComplexVector& ParticleOnLattice::ConvertToNbodyBasis(ComplexVector& state, ParticleOnLattice &nbodyBasis, int firstComponent, int nbrComponent)
{
  return state;
}


// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnLattice::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex ParticleOnLattice::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex ParticleOnLattice::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int firstComponent, int nbrComponent)
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

Complex ParticleOnLattice::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, int nextCoordinates, int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnLattice::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnLattice::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  cout << "Attention: calling dummy routine ParticleOnLattice::EvaluatePartialDensityMatrixParticlePartition"<<endl;
  return HermitianMatrix();
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. 
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnLattice::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  cout << "Attention: calling dummy routine ParticleOnLattice::EvaluatePartialDensityMatrixParticlePartition"<<endl;
  return HermitianMatrix();
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

long ParticleOnLattice::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnLattice* complementaryHilbertSpace,  ParticleOnLattice* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  return 0l;
}


// get maximum possible momentum for this geometry
// return = maximum value of Ky
int ParticleOnLattice::GetMaximumKy()
{
  cout <<"Calling ill defined function ParticleOnLattice::GetMaximumKy()"<<endl; 
  return 0;
}
