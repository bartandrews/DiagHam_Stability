////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class superconductor order parameter operator              // 
//             for particle with spin on a lattice using translations         //
//                                                                            //
//                        last modification : 09/05/2015                      //
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
#include "Operator/ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for up, 1 for plus)
// creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
// creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for up, 1 for plus)

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle, 
															   int creationMomentumIndex1, int creationSymmetryIndex1,
															   int creationMomentumIndex2, int creationSymmetryIndex2)
{
  this->Particle = (ParticleOnSphereWithSpin*) (particle->Clone());
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationSymmetryIndex1 = creationSymmetryIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->CreationSymmetryIndex2 = creationSymmetryIndex2;
  this->CombinationFlag = false;
  this->CreationSymmetryIndex1SecondTerm = creationSymmetryIndex1;
  this->CreationSymmetryIndex2SecondTerm = creationSymmetryIndex2;
  this->CombinationSign = 0.0;  
}

// constructor for operator such as a+_{sigma_1,i_1} a+_{sigma_2,i_2} +/- a+_{sigma_3,i_1} a+_{sigma_4,i_2}
//
// particle = hilbert space associated to the particles with the small number of particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationSymmetryIndex1 = symmetry index of the leftmost creation operator (0 for down, 1 for up)
// creationMomentumIndex2 = momentum index of the rightmost creation operator (from 0 to 2S)
// creationSymmetryIndex2 = symmetry index of the rightmost creation operator (0 for down, 1 for up)
// creationSymmetryIndex1SecondTerm = symmetry index of the rightmost creation operator for the second term
// creationSymmetryIndex2SecondTerm = symmetry index of the leftmost creation operator for the second term
// sign = sign in front of the second term

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* particle,  int creationMomentumIndex1, int creationSymmetryIndex1,
																					       int creationMomentumIndex2, int creationSymmetryIndex2,
																					       int creationSymmetryIndex1SecondTerm, int creationSymmetryIndex2SecondTerm, 
																					       double sign)
{
  this->Particle = (ParticleOnSphereWithSpin*) (particle->Clone());
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationSymmetryIndex1 = creationSymmetryIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->CreationSymmetryIndex2 = creationSymmetryIndex2;
  this->CombinationFlag = true;
  this->CreationSymmetryIndex1SecondTerm = creationSymmetryIndex1SecondTerm;
  this->CreationSymmetryIndex2SecondTerm = creationSymmetryIndex2SecondTerm;
  this->CombinationSign = sign;
}

// copy constructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator& oper)
{
  this->Particle = (ParticleOnSphereWithSpin*) (oper.Particle->Clone());
  this->CreationMomentumIndex1 = oper.CreationMomentumIndex1;
  this->CreationSymmetryIndex1 = oper.CreationSymmetryIndex1;
  this->CreationMomentumIndex2 = oper.CreationMomentumIndex2;
  this->CreationSymmetryIndex2 = oper.CreationSymmetryIndex2;
  this->CombinationFlag = oper.CombinationFlag;
  this->CombinationSign = oper.CombinationSign;  
  this->CreationSymmetryIndex1SecondTerm = oper.CreationSymmetryIndex1SecondTerm;
  this->CreationSymmetryIndex2SecondTerm = oper.CreationSymmetryIndex2SecondTerm;
  
}

// destructor
//

ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::~ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::Clone ()
{
  return new ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) this->Particle;
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = TmpHilbertSpace->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (TmpHilbertSpace->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;

  int CreationOrbitalIndex1;
  int CreationXPosition1;
  int CreationYPosition1;
  TmpHilbertSpace->GetLinearizedIndex(this->CreationMomentumIndex1, CreationXPosition1, CreationYPosition1, CreationOrbitalIndex1);

  int CreationOrbitalIndex2;
  int CreationXPosition2;
  int CreationYPosition2;
  TmpHilbertSpace->GetLinearizedIndex(this->CreationMomentumIndex2, CreationXPosition2, CreationYPosition2, CreationOrbitalIndex2);

  int TmpNbrUnitCells = TmpHilbertSpace->GetMaxXMomentum() * TmpHilbertSpace->GetMaxYMomentum();
  int* ShiftedCreationIndices1 = new int [TmpNbrUnitCells];
  int* ShiftedCreationIndices2 = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  Complex** LeftPhases = new Complex* [TmpHilbertSpace->GetMaxXMomentum()];
  int Tmp = 0;  
  for (int i = 0; i < TmpHilbertSpace->GetMaxXMomentum(); ++i)
    {
      LeftPhases[i] = new Complex[TmpHilbertSpace->GetMaxYMomentum()];
      for (int j = 0; j < TmpHilbertSpace->GetMaxYMomentum(); ++j)
	{
	  ShiftedCreationIndices1[Tmp] = TmpHilbertSpace->GetLinearizedIndexSafe(CreationXPosition1 - i, CreationYPosition1 - j, CreationOrbitalIndex1);
	  ShiftedCreationIndices2[Tmp] = TmpHilbertSpace->GetLinearizedIndexSafe(CreationXPosition2 - i, CreationYPosition2 - j, CreationOrbitalIndex2);
	  RightPhases[Tmp] = Phase (2.0 * M_PI * ((((double) (i * TmpHilbertSpace->GetKxMomentum())) / ((double) TmpHilbertSpace->GetMaxXMomentum()))
						  + (((double) (j * TmpHilbertSpace->GetKyMomentum())) / ((double) TmpHilbertSpace->GetMaxYMomentum())))) / ((double) TmpNbrUnitCells);
	  LeftPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * TmpHilbertSpace->TargetSpace->GetKxMomentum())) / ((double) TmpHilbertSpace->GetMaxXMomentum())
						   + ((double) (j * TmpHilbertSpace->TargetSpace->GetKyMomentum())) / ((double) TmpHilbertSpace->GetMaxYMomentum())));
	  ++Tmp;
	}
    }

  int NbrTranslationX;
  int NbrTranslationY;  

  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AddAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != FullDim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
      }
      break;
    case 3 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdu(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != FullDim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
      }
      break;
    case 1 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices2[j], ShiftedCreationIndices1[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != FullDim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (Sign * Coefficient);      
		  }
	      }
	  }
      }
      break;
    case 2 :
      {
	for (int i = (int) firstComponent; i < Dim; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != FullDim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
      }
      break;
   }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm; 
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AddAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != FullDim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (this->CombinationSign * Coefficient);      
		      }
		  }
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdu(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != FullDim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (this->CombinationSign * Coefficient);      
		      }
		  }
	      }
	  }
	  break;
	case 1 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices2[j], ShiftedCreationIndices1[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != FullDim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (this->CombinationSign * Sign * Coefficient);      
		      }
		  }
	      }
	  }
	  break;
	case 2 :
	  {
	    for (int i = (int) firstComponent; i < Dim; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != FullDim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			Element += Conj(V1[Index]) * V2[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (this->CombinationSign * Coefficient);      
		      }
		  }
	      }
	  }
	  break;
	}
   }
  delete[] ShiftedCreationIndices1;
  delete[] ShiftedCreationIndices2;
  delete[] RightPhases;
  for (int i = 0; i < TmpHilbertSpace->GetMaxXMomentum(); ++i)
    delete[] LeftPhases[i];
  delete[] LeftPhases;
  return Element;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeRealSpaceWithSpinAnd2DTranslationSuperconductorOrderParameterOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
															   int firstComponent, int nbrComponent)
{
  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* TmpHilbertSpace = (FermionOnLatticeWithSpinRealSpaceAnd2DTranslation*) this->Particle;
  int Last = firstComponent + nbrComponent;;
  int Dim = TmpHilbertSpace->GetTargetHilbertSpaceDimension();
  double Coefficient = 0.0;
  int SymmetryIndex = (this->CreationSymmetryIndex1 << 1) | this->CreationSymmetryIndex2;
  double Sign = 1.0;
  if (TmpHilbertSpace->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    Sign = -1.0;
  int CreationOrbitalIndex1;
  int CreationXPosition1;
  int CreationYPosition1;
  TmpHilbertSpace->GetLinearizedIndex(this->CreationMomentumIndex1, CreationXPosition1, CreationYPosition1, CreationOrbitalIndex1);

  int CreationOrbitalIndex2;
  int CreationXPosition2;
  int CreationYPosition2;
  TmpHilbertSpace->GetLinearizedIndex(this->CreationMomentumIndex2, CreationXPosition2, CreationYPosition2, CreationOrbitalIndex2);

  int TmpNbrUnitCells = TmpHilbertSpace->GetMaxXMomentum() * TmpHilbertSpace->GetMaxYMomentum();
  int* ShiftedCreationIndices1 = new int [TmpNbrUnitCells];
  int* ShiftedCreationIndices2 = new int [TmpNbrUnitCells];
  Complex* RightPhases = new Complex [TmpNbrUnitCells];
  Complex** LeftPhases = new Complex* [TmpHilbertSpace->GetMaxXMomentum()];
  int Tmp = 0;  
  for (int i = 0; i < TmpHilbertSpace->GetMaxXMomentum(); ++i)
    {
      LeftPhases[i] = new Complex[TmpHilbertSpace->GetMaxYMomentum()];
      for (int j = 0; j < TmpHilbertSpace->GetMaxYMomentum(); ++j)
	{
	  ShiftedCreationIndices1[Tmp] = TmpHilbertSpace->GetLinearizedIndexSafe(CreationXPosition1 - i, CreationYPosition1 - j, CreationOrbitalIndex1);
	  ShiftedCreationIndices2[Tmp] = TmpHilbertSpace->GetLinearizedIndexSafe(CreationXPosition2 - i, CreationYPosition2 - j, CreationOrbitalIndex2);
	  RightPhases[Tmp] = Phase (2.0 * M_PI * (((double) (i * TmpHilbertSpace->GetKxMomentum())) / ((double) TmpHilbertSpace->GetMaxXMomentum())
						   + ((double) (j * TmpHilbertSpace->GetKyMomentum())) / ((double) TmpHilbertSpace->GetMaxYMomentum()))) / ((double) TmpNbrUnitCells);
	  LeftPhases[i][j] =  Phase (2.0 * M_PI * (((double) (i * TmpHilbertSpace->TargetSpace->GetKxMomentum())) / ((double) TmpHilbertSpace->GetMaxXMomentum())
						   + ((double) (j * TmpHilbertSpace->TargetSpace->GetKyMomentum())) / ((double) TmpHilbertSpace->GetMaxYMomentum())));
	  ++Tmp;
	}
    }

  int NbrTranslationX;
  int NbrTranslationY;  
  switch (SymmetryIndex)
    {
    case 0 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AddAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices1[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != Dim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    vDestination[Index] += vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
      }
      break;
    case 3 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdu(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != Dim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    vDestination[Index] += vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
	break;
      }
    case 1 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices2[j], ShiftedCreationIndices1[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != Dim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    vDestination[Index] += vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * (Sign * Coefficient);      
		  }
	      }
	  }
 	break;
     }
    case 2 :
      {
	for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int j = 0; j < TmpNbrUnitCells; ++j)
	      {
		int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		if (Index != Dim)
		  {
		    NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationX < 0)
		      NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
		    NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
		    if (NbrTranslationY < 0)
		      NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
		    vDestination[Index] += vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		  }
	      }
	  }
	break;
      }
    }
  if (this->CombinationFlag == true)
    {
      SymmetryIndex = (this->CreationSymmetryIndex1SecondTerm << 1) | this->CreationSymmetryIndex2SecondTerm;
      switch (SymmetryIndex)
	{
	case 0 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AddAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != Dim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			vDestination[Index] += this->CombinationSign * vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		      }
		  }
	      }
	  }
	  break;
	case 3 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdu(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != Dim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			vDestination[Index] += this->CombinationSign * vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		      }
		  }
	      }
	    break;
	  }
	case 1 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices2[j], ShiftedCreationIndices1[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != Dim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			vDestination[Index] += this->CombinationSign * vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Sign * Coefficient;      
		      }
		  }
	      }
	    break;
	  }
	case 2 :
	  {
	    for (int i = firstComponent; i < Last; ++i)
	      {
		for (int j = 0; j < TmpNbrUnitCells; ++j)
		  {
		    int Index = TmpHilbertSpace->AduAdd(i, ShiftedCreationIndices1[j], ShiftedCreationIndices2[j], Coefficient, NbrTranslationX, NbrTranslationY);
		    if (Index != Dim)
		      {
			NbrTranslationX -= j / TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationX < 0)
			  NbrTranslationX += TmpHilbertSpace->GetMaxXMomentum();
			NbrTranslationY -= j % TmpHilbertSpace->GetMaxYMomentum();
			if (NbrTranslationY < 0)
			  NbrTranslationY += TmpHilbertSpace->GetMaxYMomentum();
			vDestination[Index] += this->CombinationSign * vSource[i] * RightPhases[j] * LeftPhases[NbrTranslationX][NbrTranslationY] * Coefficient;      
		      }
		  }
	      }
	break;
	  }
	}
    }
  delete[] ShiftedCreationIndices1;
  delete[] ShiftedCreationIndices2;
  delete[] RightPhases;
  for (int i = 0; i < TmpHilbertSpace->GetMaxXMomentum(); ++i)
    delete[] LeftPhases[i];
  delete[] LeftPhases;
  return vDestination;
}

