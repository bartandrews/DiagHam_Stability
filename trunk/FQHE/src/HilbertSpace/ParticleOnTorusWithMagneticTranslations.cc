////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of particle on a torus taking into account magnetic translations   //
//                                                                            //
//                        last modification : 10/09/2003                      //
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
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHETorusApplyCNRotationOperation.h"

#include <iostream>

using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnTorusWithMagneticTranslations::~ParticleOnTorusWithMagneticTranslations ()
{
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnTorusWithMagneticTranslations::AA (int index, int n1, int n2)
{
  return 0.0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double ParticleOnTorusWithMagneticTranslations::ProdA (int index, int* n, int nbrIndices)
{
  return 0.0;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnTorusWithMagneticTranslations::AdAd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  coefficient = 0.0;
  nbrTranslation = 0;
  return this->HilbertSpaceDimension;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnTorusWithMagneticTranslations::ProdAd (int* m, int nbrIndices, double& coefficient, int& nbrTranslation)
{
  coefficient = 0.0;
  nbrTranslation = 0;
  return this->HilbertSpaceDimension;
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& ParticleOnTorusWithMagneticTranslations::A (int i, Matrix& M)
{
  return M;
}

// return matrix representation ofthw creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& ParticleOnTorusWithMagneticTranslations::Ad (int i, Matrix& M)
{
  return M;
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

long ParticleOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
												 ComplexVector& groundState, HermitianMatrix* densityMatrix)
{
  cout << "warning, using dummy ParticleOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore" << endl;
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

long ParticleOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnTorusWithMagneticTranslations* complementaryHilbertSpace,  ParticleOnTorusWithMagneticTranslations* destinationHilbertSpace,
												 int nbrGroundStates, ComplexVector* groundStates, double* weights, HermitianMatrix* densityMatrix)
{
  cout << "warning, using dummy ParticleOnTorusWithMagneticTranslations::EvaluatePartialDensityMatrixParticlePartitionCore" << endl;
  return 0l;
}

// convert a state defined in the Ky basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnTorusWithMagneticTranslations::ConvertToKxKyBasis(ComplexVector& state, ParticleOnTorus* space)
{
  return ComplexVector();
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnTorusWithMagneticTranslations::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnTorus* space)
{
  return ComplexVector();
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// oldKx = vnew value of the relative quantum number kx
// return value = state in the (Kx,Ky) basis

ComplexVector ParticleOnTorusWithMagneticTranslations::ConvertFromKxKyBasisAndModifyKx(ComplexVector& state, ParticleOnTorus* space, int oldKx)
{
  return ComplexVector();
}

// get the C2 symmetric state of a given state 
//
// index = index of the state whose symmetric counterpart has to be computed
// nbrTranslation = number of translations that has to be applied to C2 symmetric state to be in the canonical form

int ParticleOnTorusWithMagneticTranslations::GetC2SymmetricState (int index, int& nbrTranslation)
{
  cout << "warning : GetC2SymmetricState not implemented" << endl;
  return 0;
}

// apply the C4 rotation to a given state assumin it is an eigenstate of both kx and ky
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// architecture = pointer to the architecture
// clockwise = the rotation is done clockwise
// return value = rotated state

ComplexVector ParticleOnTorusWithMagneticTranslations::C4Rotation (ComplexVector& inputState, ParticleOnTorusWithMagneticTranslations* inputSpace, bool clockwise, AbstractArchitecture* architecture)
{
  ComplexVector OutputState (this->HilbertSpaceDimension, true);
  if (architecture != 0)
    {
      if (clockwise == false)
	{
	  FQHETorusApplyCNRotationOperation Operation(4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
      else
	{
	  FQHETorusApplyCNRotationOperation Operation(-4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
    }
  else
    {
      this->CoreC4Rotation(inputState, inputSpace, OutputState, 0 , this->HilbertSpaceDimension, clockwise);
    }
  return OutputState;
}

// core part of the C4 rotation
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// outputState = reference on the rotated state
// minIndex = minimum index that has to be computed
// nbrIndices = number of indices that have to be computed
// clockwise = the rotation is done clockwise
// return value = reference on the rotated state

ComplexVector& ParticleOnTorusWithMagneticTranslations::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorusWithMagneticTranslations* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise)
{
  cout << "warning : CoreC4Rotation not implemented" << endl;
  return outputState;
}

// get the momentum along the x axis
// 
// return avlue = momentum along the x axis

int ParticleOnTorusWithMagneticTranslations::GetKxMomentum()
{
  return 0;
}

// get the momentum along the y axis
// 
// return avlue = momentum along the y axis

int ParticleOnTorusWithMagneticTranslations::GetKyMomentum()
{
  return 0;
}

// get the maximum momentum along the x axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the x axis

int ParticleOnTorusWithMagneticTranslations::GetMaxXMomentum()
{
  return 1;
}
  
// get the maximum momentum along the y axis (i.e. the number of momentum sectors)
// 
// return avlue = maximum momentum along the y axis

int ParticleOnTorusWithMagneticTranslations::GetMaxYMomentum()
{
  return 1;
}
  
