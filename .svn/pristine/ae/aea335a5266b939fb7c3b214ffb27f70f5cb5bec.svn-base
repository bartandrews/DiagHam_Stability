////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of hamiltonian associated to periodic DMRG algorithm       //
//                                                                            //
//                        last modification : 10/05/2002                      //
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


#include "Hamiltonian/DMRGHamiltonian/PeriodicDMRGHamiltonian.h"
#include "TensorProduct/TensorProductRealVector.h"
#include "TensorProduct/FullTensorProductStructure.h"
#include "MathTools/Complex.h"

  
// constructor from default datas
//
// leftBlockLeftInteractionPart = tensor associated to hamiltonian part describing 
// left block and interaction between left block and new left added block 
// rightBlockRightInteractionPart = tensor associated to hamiltonian part describing 
// right block and interaction between right block and new right added block
// interaction = reference on matrix describing interaction between the two blocks  
// spaceStructure = tensor product struture associated to the Hilbert space where 
// Hamiltonian has to be applied
// leftSpaceFineStructureArray = array of fine structure describing total left space
// rightSpaceFineStructureArray = array of fine structure describing total right space
// blockTotalDimension = total dimension of Hilbert space associated to left or right block
// interactionBlockTotalDimension = total dimension of Hilbert space associated to left or right 
// interaction block

PeriodicDMRGHamiltonian::PeriodicDMRGHamiltonian(OneSpaceTensor leftBlockLeftInteractionPart, 
						 OneSpaceTensor rightBlockRightInteractionPart, 
						 Matrix& interaction,
						 CompositeTensorProductStructure* spaceStructure,
						 FullTensorProductStructure** 
						 leftSpaceFineStructureArray, 
						 FullTensorProductStructure** 
						 rightSpaceFineStructureArray, 
						 int blockTotalDimension,
						 int interactionBlockTotalDimension)
{
//  cout << interaction << endl;
//  cout << leftBlockLeftInteractionPart << endl;
  this->LeftBlockLeftInteractionPart = leftBlockLeftInteractionPart;
  this->RightBlockRightInteractionPart = rightBlockRightInteractionPart;
  this->SpaceStructure = spaceStructure;
  this->Interaction = interaction.Clone();
  this->TestInteractingStateFactor = 0;
  this->BlockTotalDimension = blockTotalDimension;
  this->InteractionBlockTotalDimension = interactionBlockTotalDimension;

  this->InteractionPrecision = 1e-14;
  int TmpFactor = 1;
  while (TmpFactor < this->Interaction->GetNbrRow())
    {
       TmpFactor <<= 1;
       this->TestInteractingStateFactor++;
    }
  //  cout << this->TestInteractingStateFactor << endl;
  this->TestInteractingStateFlags = new bool [this->Interaction->GetNbrColumn() 
					      << this->TestInteractingStateFactor];
  for (int i = 0; i < this->Interaction->GetNbrColumn(); i++)
    {
      TmpFactor = i << this->TestInteractingStateFactor;
      for (int j = 0; j < this->Interaction->GetNbrRow(); j++)
	if (fabs((*(this->Interaction))(i, j)) < this->InteractionPrecision)
	  this->TestInteractingStateFlags[TmpFactor + j] = false;
	else
	  this->TestInteractingStateFlags[TmpFactor + j] = true;	  
    }
  this->InitializeNonNullInteractions (leftSpaceFineStructureArray, rightSpaceFineStructureArray);
}

// destructor
//

PeriodicDMRGHamiltonian::~PeriodicDMRGHamiltonian() 
{
  int TotalDim = this->SpaceStructure->GetTotalDimension();
  for (int i = 0; i < TotalDim; ++i)
    {
      delete[] this->InteractionFactorBulk[i];
      delete[] this->InteractionIndexBulk[i];
    }
  delete[] this->InteractionFactorBulk;
  delete[] this->InteractionIndexBulk;
  delete[] this->GroupPositionBulk;
  delete[] this->GroupSizeBulk;
  delete[] this->LeftRightAddedBlockGlobalIndexBulk;
  delete[] this->LeftRightAddedBlockIndexBulk;
  delete[] this->NbrInteractionBulk;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void PeriodicDMRGHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* PeriodicDMRGHamiltonian::GetHilbertSpace () 
{
  return 0;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int PeriodicDMRGHamiltonian::GetHilbertSpaceDimension () 
{
  return this->SpaceStructure->GetTotalDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void PeriodicDMRGHamiltonian::ShiftHamiltonian (double shift) 
{
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicDMRGHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  TensorProductRealVector TmpVectorSource (this->SpaceStructure, V2);
  RealVector Tmp (this->SpaceStructure->GetTotalDimension());
  TensorProductRealVector TmpVectorDestination (this->SpaceStructure, Tmp);
  for (int i = 0; i < Tmp.GetVectorDimension(); i++)
    Tmp[i] = 0.0;
  TmpVectorDestination.Multiply(this->LeftBlockLeftInteractionPart, TmpVectorSource);
  TmpVectorDestination.AddMultiply(this->RightBlockRightInteractionPart, TmpVectorSource);
  this->InteractionAddMultiply(V2, Tmp, 0, this->SpaceStructure->GetTotalDimension());
  return Complex (V1 * Tmp, 0.0);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex PeriodicDMRGHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& PeriodicDMRGHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  TensorProductRealVector TmpVectorSource (this->SpaceStructure, vSource);
  TensorProductRealVector TmpVectorDestination (this->SpaceStructure, vDestination);
  for (int i = 0; i < vDestination.GetVectorDimension(); i++)
    vDestination[i] = 0.0;
  TmpVectorDestination.Multiply(this->LeftBlockLeftInteractionPart, TmpVectorSource);
  TmpVectorDestination.AddMultiply(this->RightBlockRightInteractionPart, TmpVectorSource);
  this->InteractionAddMultiply(vSource, vDestination, 0, this->SpaceStructure->GetTotalDimension());
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PeriodicDMRGHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
						      int firstComponent, int nbrComponent) 
{
  TensorProductRealVector TmpVectorSource (this->SpaceStructure, vSource);
  TensorProductRealVector TmpVectorDestination (this->SpaceStructure, vDestination);
  TmpVectorDestination.Multiply(this->LeftBlockLeftInteractionPart, TmpVectorSource, firstComponent, 
				nbrComponent);
  TmpVectorDestination.AddMultiply(this->RightBlockRightInteractionPart, TmpVectorSource, firstComponent,
				   nbrComponent);
  this->InteractionAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& PeriodicDMRGHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  TensorProductRealVector TmpVectorSource (this->SpaceStructure, vSource);
  TensorProductRealVector TmpVectorDestination (this->SpaceStructure, vDestination);
  for (int i = 0; i < vDestination.GetVectorDimension(); i++)
    vDestination[i] = 0.0;
  TmpVectorDestination.AddMultiply(this->LeftBlockLeftInteractionPart, TmpVectorSource);
  TmpVectorDestination.AddMultiply(this->RightBlockRightInteractionPart, TmpVectorSource);
  this->InteractionAddMultiply(vSource, vDestination, 0, this->SpaceStructure->GetTotalDimension());
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& PeriodicDMRGHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
							    int firstComponent, int nbrComponent)
{
  TensorProductRealVector TmpVectorSource (this->SpaceStructure, vSource);
  TensorProductRealVector TmpVectorDestination (this->SpaceStructure, vDestination);
  TmpVectorDestination.AddMultiply(this->LeftBlockLeftInteractionPart, TmpVectorSource, firstComponent, 
				   nbrComponent);
  TmpVectorDestination.AddMultiply(this->RightBlockRightInteractionPart, TmpVectorSource, firstComponent,
				   nbrComponent);
  this->InteractionAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicDMRGHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& PeriodicDMRGHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& PeriodicDMRGHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							    int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& PeriodicDMRGHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							 int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> PeriodicDMRGHamiltonian::LeftInteractionOperators() 
{
  return List<Matrix*>();
}

// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> PeriodicDMRGHamiltonian::RightInteractionOperators() 
{
  return List<Matrix*>();
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, PeriodicDMRGHamiltonian& H) 
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicDMRGHamiltonian& H) 
{
  return Str;
}

// find all components of a state that are coupled by interaction  
//
// leftSpaceFineStructureArray = array of fine structure describing total left space
// rightSpaceFineStructureArray = array of fine structure describing total right space

void PeriodicDMRGHamiltonian::InitializeNonNullInteractions(FullTensorProductStructure** 
							    leftSpaceFineStructureArray, 
							    FullTensorProductStructure** 
							    rightSpaceFineStructureArray)
{
  int TotalDim = this->SpaceStructure->GetTotalDimension();
  int Index = 0;
  this->LeftRightAddedBlockGlobalIndexBoundary = new int[TotalDim];
  this->LeftRightAddedBlockIndexBoundary = new int[TotalDim];
  this->LeftRightAddedBlockGlobalIndexBulk = new int[TotalDim];
  this->LeftRightAddedBlockIndexBulk = new int[TotalDim];
  int* TmpLeftRightAddedBlockIndexBoundary = new int[TotalDim];
  int* TmpLeftRightAddedBlockIndexBulk = new int[TotalDim];

  int NbrSubspace = this->SpaceStructure->GetNbrSubspace();
  int Index2;
  int Index3;
  int Index4;
  int Index5;
  int Index6;
  int Index7;
  int RightFactor = 0;
  while ((1 << RightFactor) < this->BlockTotalDimension)
    RightFactor++;
  int* SortTableBulk = new int [this->BlockTotalDimension << (RightFactor + this->TestInteractingStateFactor)];
  int* SortTableBoundary = new int [this->BlockTotalDimension << (RightFactor + this->TestInteractingStateFactor)];
  int ReducedSortTableDimension = this->BlockTotalDimension << RightFactor;
  int* TmpNbrInteractionBulk = new int [ReducedSortTableDimension];
  int* TmpNbrInteractionBoundary = new int [ReducedSortTableDimension];
  for (int i = 0; i < ReducedSortTableDimension; i++)
    {
      TmpNbrInteractionBulk[i] = 0;
      TmpNbrInteractionBoundary[i] = 0;
    }
  this->NbrInteractingStateGroupBulk = 0;
  this->NbrInteractingStateGroupBoundary = 0;

  for (int i = 0 ; i < NbrSubspace; i++)
    {
      for (int j2 = 0 ; j2 < rightSpaceFineStructureArray[i]->GetTotalDimension(); j2++)
	{
	  Index2 = (*rightSpaceFineStructureArray[i])(j2, 0);
	  Index5 = (*rightSpaceFineStructureArray[i])(j2, 1);
	  Index3 = Index2 * this->InteractionBlockTotalDimension;
	  for (int j1 = 0 ; j1 < leftSpaceFineStructureArray[i]->GetTotalDimension(); j1++)
	    {
	      Index6 = (*leftSpaceFineStructureArray[i])(j1, 0);
	      Index7 = (*leftSpaceFineStructureArray[i])(j1, 1);

	      Index4 = (Index6 << RightFactor) | Index5;
	      SortTableBulk[(Index4 << this->TestInteractingStateFactor) + TmpNbrInteractionBulk[Index4]] = Index;
	      if (TmpNbrInteractionBulk[Index4] == 0)
		++this->NbrInteractingStateGroupBulk;
	      ++TmpNbrInteractionBulk[Index4];
	      TmpLeftRightAddedBlockIndexBulk[Index] = Index3 + Index7;

	      Index4 = (Index2 << RightFactor) | Index7;
	      SortTableBoundary[(Index4 << this->TestInteractingStateFactor) + TmpNbrInteractionBoundary[Index4]] = Index;
	      if (TmpNbrInteractionBoundary[Index4] == 0)
		++this->NbrInteractingStateGroupBoundary;
	      ++TmpNbrInteractionBoundary[Index4];
	      TmpLeftRightAddedBlockIndexBoundary[Index] = Index6 * this->InteractionBlockTotalDimension + Index5;

	      ++Index;
	    }
	}
    }

  Index = 0;
  Index3 = 0;
  Index4 = 0;
  Index5 = 0;
  this->GroupPositionBulk = new int [this->NbrInteractingStateGroupBulk];
  this->GroupSizeBulk = new int [this->NbrInteractingStateGroupBulk];
  this->GroupPositionBoundary = new int [this->NbrInteractingStateGroupBoundary];
  this->GroupSizeBoundary = new int [this->NbrInteractingStateGroupBoundary];
  for (int i = 0; i < ReducedSortTableDimension; i++)
    {
      if (TmpNbrInteractionBulk[i] != 0)
	{
	  this->GroupPositionBulk[Index3] = Index;
	  this->GroupSizeBulk[Index3++] = TmpNbrInteractionBulk[i];
	  Index2 = i << this->TestInteractingStateFactor;
	  for (int j = 0; j < TmpNbrInteractionBulk[i]; j++)
	    {
	      this->LeftRightAddedBlockGlobalIndexBulk[Index] =  SortTableBulk[Index2 + j];
	      this->LeftRightAddedBlockIndexBulk[Index] = TmpLeftRightAddedBlockIndexBulk
		[this->LeftRightAddedBlockGlobalIndexBulk[Index]];
	      Index++;
	    }
	}
      if (TmpNbrInteractionBoundary[i] != 0)
	{
	  this->GroupPositionBoundary[Index4] = Index5;
	  this->GroupSizeBoundary[Index4++] = TmpNbrInteractionBoundary[i];
	  Index2 = i << this->TestInteractingStateFactor;
	  for (int j = 0; j < TmpNbrInteractionBoundary[i]; j++)
	    {
	      this->LeftRightAddedBlockGlobalIndexBoundary[Index5] =  SortTableBoundary[Index2 + j];
	      this->LeftRightAddedBlockIndexBoundary[Index5] = TmpLeftRightAddedBlockIndexBoundary
		[this->LeftRightAddedBlockGlobalIndexBoundary[Index5]];
	      Index5++;
	    }
	}
    }
  delete[] TmpLeftRightAddedBlockIndexBulk;
  delete[] TmpLeftRightAddedBlockIndexBoundary;
  delete[] TmpNbrInteractionBulk;
  delete[] TmpNbrInteractionBoundary;
  delete[] SortTableBulk;
  delete[] SortTableBoundary;

  int k2;
  this->NbrInteractionBulk = new int [TotalDim];
  this->InteractionIndexBulk = new int* [TotalDim];
  this->InteractionFactorBulk = new double* [TotalDim];
  for (int g = 0; g < this->NbrInteractingStateGroupBulk; g++)
    {
      Index = this->GroupPositionBulk[g];
      for (int k1 = 0; k1 < this->GroupSizeBulk[g]; k1++)
	{
	  this->NbrInteractionBulk[Index + k1] = 0;
	  this->InteractionIndexBulk[Index + k1] = new int [this->GroupSizeBulk[g]];
	  this->InteractionFactorBulk[Index + k1] = new double [this->GroupSizeBulk[g]];
	}
      for (int k1 = 0; k1 < this->GroupSizeBulk[g]; k1++)
	{
	  k2 = k1 + 1;
	  Index2 = this->LeftRightAddedBlockIndexBulk[Index + k1] << this->TestInteractingStateFactor;
	  if (this->TestInteractingStateFlags[Index2 | 
					     this->LeftRightAddedBlockIndexBulk[Index + k1]] == true)
	    {
	      this->InteractionIndexBulk[Index + k1][this->NbrInteractionBulk[Index + k1]] = 
		this->LeftRightAddedBlockGlobalIndexBulk[Index + k1];
	      this->InteractionFactorBulk[Index + k1][this->NbrInteractionBulk[Index + k1]] = 
		(*(this->Interaction))(this->LeftRightAddedBlockIndexBulk[Index + k1], 
				       this->LeftRightAddedBlockIndexBulk[Index + k1]);
	      this->NbrInteractionBulk[Index + k1]++;
	    }
	  while (k2 < this->GroupSizeBulk[g])
	    {
	      if (this->TestInteractingStateFlags[Index2 | 
						 this->LeftRightAddedBlockIndexBulk[Index + k2]] == true)
		{
		  this->InteractionIndexBulk[Index + k1][this->NbrInteractionBulk[Index + k1]] = 
		    this->LeftRightAddedBlockGlobalIndexBulk[Index + k2];
		  this->InteractionIndexBulk[Index + k2][this->NbrInteractionBulk[Index + k2]] = 
		    this->LeftRightAddedBlockGlobalIndexBulk[Index + k1];
		  this->InteractionFactorBulk[Index + k1][this->NbrInteractionBulk[Index + k1]] = 
		    (*(this->Interaction))(this->LeftRightAddedBlockIndexBulk[Index + k1], 
					   this->LeftRightAddedBlockIndexBulk[Index + k2]);
		  this->InteractionFactorBulk[Index + k2][this->NbrInteractionBulk[Index + k2]] = 
		    this->InteractionFactorBulk[Index + k1][this->NbrInteractionBulk[Index + k1]];
		  this->NbrInteractionBulk[Index + k1]++;
		  this->NbrInteractionBulk[Index + k2]++;
		}
	      k2++;
	    }
	}
    }

  this->NbrInteractionBoundary = new int [TotalDim];
  this->InteractionIndexBoundary = new int* [TotalDim];
  this->InteractionFactorBoundary = new double* [TotalDim];
  for (int g = 0; g < this->NbrInteractingStateGroupBoundary; g++)
    {
      Index = this->GroupPositionBoundary[g];
      for (int k1 = 0; k1 < this->GroupSizeBoundary[g]; k1++)
	{
	  this->NbrInteractionBoundary[Index + k1] = 0;
	  this->InteractionIndexBoundary[Index + k1] = new int [this->GroupSizeBoundary[g]];
	  this->InteractionFactorBoundary[Index + k1] = new double [this->GroupSizeBoundary[g]];
	}
      for (int k1 = 0; k1 < this->GroupSizeBoundary[g]; k1++)
	{
	  k2 = k1 + 1;
	  Index2 = this->LeftRightAddedBlockIndexBoundary[Index + k1] << this->TestInteractingStateFactor;
	  if (this->TestInteractingStateFlags[Index2 | 
					     this->LeftRightAddedBlockIndexBoundary[Index + k1]] == true)
	    {
	      this->InteractionIndexBoundary[Index + k1][this->NbrInteractionBoundary[Index + k1]] = 
		this->LeftRightAddedBlockGlobalIndexBoundary[Index + k1];
	      this->InteractionFactorBoundary[Index + k1][this->NbrInteractionBoundary[Index + k1]] = 
		(*(this->Interaction))(this->LeftRightAddedBlockIndexBoundary[Index + k1], 
				       this->LeftRightAddedBlockIndexBoundary[Index + k1]);
	      this->NbrInteractionBoundary[Index + k1]++;
	    }
	  while (k2 < this->GroupSizeBoundary[g])
	    {
	      if (this->TestInteractingStateFlags[Index2 | 
						 this->LeftRightAddedBlockIndexBoundary[Index + k2]] == true)
		{
		  this->InteractionIndexBoundary[Index + k1][this->NbrInteractionBoundary[Index + k1]] = 
		    this->LeftRightAddedBlockGlobalIndexBoundary[Index + k2];
		  this->InteractionIndexBoundary[Index + k2][this->NbrInteractionBoundary[Index + k2]] = 
		    this->LeftRightAddedBlockGlobalIndexBoundary[Index + k1];
		  this->InteractionFactorBoundary[Index + k1][this->NbrInteractionBoundary[Index + k1]] = 
		    (*(this->Interaction))(this->LeftRightAddedBlockIndexBoundary[Index + k1], 
					   this->LeftRightAddedBlockIndexBoundary[Index + k2]);
		  this->InteractionFactorBoundary[Index + k2][this->NbrInteractionBoundary[Index + k2]] = 
		    this->InteractionFactorBoundary[Index + k1][this->NbrInteractionBoundary[Index + k1]];
		  this->NbrInteractionBoundary[Index + k1]++;
		  this->NbrInteractionBoundary[Index + k2]++;
		}
	      k2++;
	    }
	}
    }
}

// multiply a vector by the current hamiltonian part corresponding to interaction between left and 
// right blocks and add result to another vector, for a given range of components
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vector where result has been stored

RealVector& PeriodicDMRGHamiltonian::InteractionAddMultiply(RealVector& vSource, 
							    RealVector& vDestination,
							    int firstComponent, int NbrComponent) 
{
  int TotalDim = firstComponent + NbrComponent;
  double* TmpInteractionFactorBulk;
  int* TmpInteractionIndexBulk;
  int LocalNbrInteractionBulk;
  double CoefBulk;
  double* TmpInteractionFactorBoundary;
  int* TmpInteractionIndexBoundary;
  int LocalNbrInteractionBoundary;
  double CoefBoundary;
  int j;
  for (int i = firstComponent; i < TotalDim; ++i)
    {
      TmpInteractionFactorBoundary = this->InteractionFactorBoundary[i];
      TmpInteractionIndexBoundary = this->InteractionIndexBoundary[i];
      LocalNbrInteractionBoundary = this->NbrInteractionBoundary[i];
      CoefBoundary = 0.0;
      TmpInteractionFactorBulk = this->InteractionFactorBulk[i];
      TmpInteractionIndexBulk = this->InteractionIndexBulk[i];
      LocalNbrInteractionBulk = this->NbrInteractionBulk[i];
      CoefBulk = 0.0;

      if (LocalNbrInteractionBulk <= LocalNbrInteractionBoundary)
	{
	  for (j = 0; j < LocalNbrInteractionBulk; ++j)
	    {
	      CoefBoundary += TmpInteractionFactorBoundary[j] * vSource[TmpInteractionIndexBoundary[j]];
	      CoefBulk += TmpInteractionFactorBulk[j] * vSource[TmpInteractionIndexBulk[j]];
	    }
	  for (; j < LocalNbrInteractionBoundary; ++j)
	    CoefBoundary += TmpInteractionFactorBoundary[j] * vSource[TmpInteractionIndexBoundary[j]];	  
	}
      else
	{
	  for (j = 0; j < LocalNbrInteractionBoundary; ++j)
	    {
	      CoefBoundary += TmpInteractionFactorBoundary[j] * vSource[TmpInteractionIndexBoundary[j]];
	      CoefBulk += TmpInteractionFactorBulk[j] * vSource[TmpInteractionIndexBulk[j]];
	    }
	  for (; j < LocalNbrInteractionBulk; ++j)
	    CoefBulk += TmpInteractionFactorBulk[j] * vSource[TmpInteractionIndexBulk[j]];
	}
      vDestination[this->LeftRightAddedBlockGlobalIndexBoundary[i]] += CoefBoundary;
      vDestination[this->LeftRightAddedBlockGlobalIndexBulk[i]] += CoefBulk;
    }
  return vDestination;
}

