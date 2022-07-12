////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of particle on a torus with spin                  //
//                                                                            //
//                        last modification : 10/09/2002                      //
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
#include "HilbertSpace/ParticleOnTorusWithSpinAndMagneticTranslations.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnTorusWithSpinAndMagneticTranslations::~ParticleOnTorusWithSpinAndMagneticTranslations ()
{
}



// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double ParticleOnTorusWithSpinAndMagneticTranslations::AuAu (int index, int n1, int n2)
{
  cout << "Warning: using dummy function ParticleOnTorusWithSpinAndMagneticTranslations::AuAu"<<endl;
  return 0.0;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnTorusWithSpinAndMagneticTranslations:: AdAd (int index, int n1, int n2)
{
  cout << "Warning: using dummy function ParticleOnTorusWithSpinAndMagneticTranslations::AdAd"<<endl;
  return 0.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnTorusWithSpinAndMagneticTranslations::AuAd (int index, int n1, int n2)
{
  cout << "Warning: using dummy function ParticleOnTorusWithSpinAndMagneticTranslations::AuAd"<<endl;
  return 0.0;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double ParticleOnTorusWithSpinAndMagneticTranslations::AdAu (int index, int n1, int n2)
{
  cout << "Warning: using dummy function ParticleOnTorusWithSpinAndMagneticTranslations::AdAu"<<endl;
  return 0.0;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_d a_m_u operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::AddAu (int index, int m, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_u a_m_d operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::AduAd (int index, int m, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to be applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnTorusWithSpinAndMagneticTranslations::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnTorusWithSpinAndMagneticTranslations::ConvertToKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  ComplexVector TmpVector;
  return TmpVector;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnTorusWithSpinAndMagneticTranslations::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  ComplexVector TmpVector;
  return TmpVector;
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

long ParticleOnTorusWithSpinAndMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, 
												 ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
												 ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
												 ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  cout << "warning : EvaluatePartialDensityMatrixParticlePartitionCore not implemented" << endl; 
  return 0l; 
}

// core part of the evaluation density matrix particle partition calculation involving a sum of projectors
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// nbrGroundStates = number of projectors
// groundStates = array of degenerate groundstates associated to each projector
// weights = array of weights in front of each projector
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long ParticleOnTorusWithSpinAndMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithSpinAndMagneticTranslations* complementaryHilbertSpace,  
													ParticleOnTorusWithSpinAndMagneticTranslations* destinationHilbertSpace,
													int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  cout << "warning : EvaluatePartialDensityMatrixParticlePartitionCore not implemented" << endl; 
  return 0l; 
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

int ParticleOnTorusWithSpinAndMagneticTranslations::GetKxMomentum()
{
  return 0;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

int ParticleOnTorusWithSpinAndMagneticTranslations::GetKyMomentum()
{
  return 0;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

int ParticleOnTorusWithSpinAndMagneticTranslations::GetMaxXMomentum()
{
  return 1;
}
  
// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

int ParticleOnTorusWithSpinAndMagneticTranslations::GetMaxYMomentum()
{
  return 1;
}
  
