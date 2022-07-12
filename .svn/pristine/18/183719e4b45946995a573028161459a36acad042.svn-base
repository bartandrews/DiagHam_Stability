////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere with SU(4) spin                //
//                                                                            //
//                        last modification : 29/11/2006                      //
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
#include "HilbertSpace/ParticleOnSphereWithSU4Spin.h"


// virtual destructor
//

ParticleOnSphereWithSU4Spin::~ParticleOnSphereWithSU4Spin ()
{
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU4Spin::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSU4Spin::AdA (int index, int m)
{
  return (this->AdupAup(index, m) + this->AdumAum(index, m) + this->AddpAdp(index, m) + this->AddmAdm(index, m));
}

// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSU4Spin::AdA (long index, int m)
{
  return (this->AdupAup(index, m) + this->AdumAum(index, m) + this->AddpAdp(index, m) + this->AddmAdm(index, m));
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSU4Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
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

Complex ParticleOnSphereWithSU4Spin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
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
Complex ParticleOnSphereWithSU4Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
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

Complex ParticleOnSphereWithSU4Spin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
									 AbstractFunctionBasis& basis, 
									 int nextCoordinates, int firstComponent, 
									 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnSphereWithSU4Spin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrNUpPlusSector = number of particles with quantum number up-plus that belong to the subsytem 
// nbrNUpMinusSector = number of particles with quantum number up-minus that belong to the subsytem 
// nbrNDownPlusSector = number of particles with quantum number down-plus that belong to the subsytem 
// nbrNDownMinusSector = number of particles with quantum number down-plus that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSU4Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, 
												int nbrNUpPlusSector, int nbrNUpMinusSector, int nbrNDownPlusSector, int nbrNDownMinusSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}
