////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to  non-periodic DMRG algorithm     //
//                                                                            //
//                        last modification : 08/06/2001                      //
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


#include "Hamiltonian/DMRGHamiltonian/NonPeriodicDMRGHamiltonian.h"
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

NonPeriodicDMRGHamiltonian::NonPeriodicDMRGHamiltonian(OneSpaceTensor leftBlockLeftInteractionPart, 
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

NonPeriodicDMRGHamiltonian::~NonPeriodicDMRGHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void NonPeriodicDMRGHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* NonPeriodicDMRGHamiltonian::GetHilbertSpace () 
{
  return 0;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int NonPeriodicDMRGHamiltonian::GetHilbertSpaceDimension () 
{
  return this->SpaceStructure->GetTotalDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void NonPeriodicDMRGHamiltonian::ShiftHamiltonian (double shift) 
{
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
bool NonPeriodicDMRGHamiltonian::SavePrecalculation (char* fileName)
{
  return false;
}


// return matrix representation of current Hamiltonian
//
// return value = reference to representation

Matrix* NonPeriodicDMRGHamiltonian::GetHamiltonian () 
{
  return 0;
}

// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex NonPeriodicDMRGHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex NonPeriodicDMRGHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& NonPeriodicDMRGHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
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

RealVector& NonPeriodicDMRGHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
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

RealVector& NonPeriodicDMRGHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
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

RealVector& NonPeriodicDMRGHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
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

ComplexVector& NonPeriodicDMRGHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& NonPeriodicDMRGHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
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
ComplexVector& NonPeriodicDMRGHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
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

ComplexVector& NonPeriodicDMRGHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  return vDestination;
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> NonPeriodicDMRGHamiltonian::LeftInteractionOperators() 
{
  return List<Matrix*>();
}

// return a list of right interaction operators 
//
// return value = list of right interaction operators

List<Matrix*> NonPeriodicDMRGHamiltonian::RightInteractionOperators() 
{
  return List<Matrix*>();
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, NonPeriodicDMRGHamiltonian& H) 
{
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, NonPeriodicDMRGHamiltonian& H) 
{
  return Str;
}

// find all components of a state that are coupled by interaction  
//
// leftSpaceFineStructureArray = array of fine structure describing total left space
// rightSpaceFineStructureArray = array of fine structure describing total right space

void NonPeriodicDMRGHamiltonian::InitializeNonNullInteractions(FullTensorProductStructure** 
							       leftSpaceFineStructureArray, 
							       FullTensorProductStructure** 
							       rightSpaceFineStructureArray)
{
  int TotalDim = this->SpaceStructure->GetTotalDimension();
  int Index = 0;
  this->LeftRightAddedBlockGlobalIndex = new int[TotalDim];
  this->LeftRightAddedBlockIndex = new int[TotalDim];
  int* TmpLeftRightAddedBlockIndex = new int[TotalDim];
  int NbrSubspace = this->SpaceStructure->GetNbrSubspace();
  int Index2;
  int Index3;
  int Index4;
  int RightFactor = 0;
  while ((1 << RightFactor) < this->BlockTotalDimension)
    RightFactor++;
  int* SortTable = new int [this->BlockTotalDimension << (RightFactor + this->TestInteractingStateFactor)];
  int ReducedSortTableDimension = this->BlockTotalDimension << RightFactor;
  int* TmpNbrInteraction = new int [ReducedSortTableDimension];
  for (int i = 0; i < ReducedSortTableDimension; i++)
    TmpNbrInteraction[i] = 0;
  this->NbrInteractingStateGroup = 0;
  for (int i = 0 ; i < NbrSubspace; i++)
    {
      for (int j2 = 0 ; j2 < rightSpaceFineStructureArray[i]->GetTotalDimension(); j2++)
	{
	  Index2 = (*rightSpaceFineStructureArray[i])(j2, 0) << RightFactor;
	  Index3 = (*rightSpaceFineStructureArray[i])(j2, 1) * this->InteractionBlockTotalDimension;
	  for (int j1 = 0 ; j1 < leftSpaceFineStructureArray[i]->GetTotalDimension(); j1++)
	    {
	      Index4 = Index2 | (*leftSpaceFineStructureArray[i])(j1, 0);
	      SortTable[(Index4 << this->TestInteractingStateFactor) + TmpNbrInteraction[Index4]] = Index;
	      if (TmpNbrInteraction[Index4] == 0)
		this->NbrInteractingStateGroup++;
	      TmpNbrInteraction[Index4]++;
	      TmpLeftRightAddedBlockIndex[Index++] = Index3 + (*leftSpaceFineStructureArray[i])(j1, 1);
	    }
	}
    }
  Index = 0;
  Index3 = 0;
  this->GroupPosition = new int [this->NbrInteractingStateGroup];
  this->GroupSize = new int [this->NbrInteractingStateGroup];
  for (int i = 0; i < ReducedSortTableDimension; i++)
    {
      if (TmpNbrInteraction[i] != 0)
	{
	  this->GroupPosition[Index3] = Index;
	  this->GroupSize[Index3++] = TmpNbrInteraction[i];
	  Index2 = i << this->TestInteractingStateFactor;
	  for (int j = 0; j < TmpNbrInteraction[i]; j++)
	    {
	      this->LeftRightAddedBlockGlobalIndex[Index] =  SortTable[Index2 + j];
	      this->LeftRightAddedBlockIndex[Index] = TmpLeftRightAddedBlockIndex
		[this->LeftRightAddedBlockGlobalIndex[Index]];
	      Index++;
	    }
	}
    }
  delete[] TmpLeftRightAddedBlockIndex;
  delete[] TmpNbrInteraction;
  delete[] SortTable;

  int k2;
  this->NbrInteraction = new int [TotalDim];
  this->InteractionIndex = new int* [TotalDim];
  this->InteractionFactor = new double* [TotalDim];
  for (int g = 0; g < this->NbrInteractingStateGroup; g++)
    {
      Index = this->GroupPosition[g];
      for (int k1 = 0; k1 < this->GroupSize[g]; k1++)
	{
	  this->NbrInteraction[Index + k1] = 0;
	  this->InteractionIndex[Index + k1] = new int [this->GroupSize[g]];
	  this->InteractionFactor[Index + k1] = new double [this->GroupSize[g]];
	}
      for (int k1 = 0; k1 < this->GroupSize[g]; k1++)
	{
	  k2 = k1 + 1;
	  Index2 = this->LeftRightAddedBlockIndex[Index + k1] << this->TestInteractingStateFactor;
	  if (this->TestInteractingStateFlags[Index2 | 
					     this->LeftRightAddedBlockIndex[Index + k1]] == true)
	    {
	      this->InteractionIndex[Index + k1][this->NbrInteraction[Index + k1]] = 
		this->LeftRightAddedBlockGlobalIndex[Index + k1];
	      this->InteractionFactor[Index + k1][this->NbrInteraction[Index + k1]] = 
		(*(this->Interaction))(this->LeftRightAddedBlockIndex[Index + k1], 
				       this->LeftRightAddedBlockIndex[Index + k1]);
	      this->NbrInteraction[Index + k1]++;
	    }
	  while (k2 < this->GroupSize[g])
	    {
	      if (this->TestInteractingStateFlags[Index2 | 
						 this->LeftRightAddedBlockIndex[Index + k2]] == true)
		{
		  this->InteractionIndex[Index + k1][this->NbrInteraction[Index + k1]] = 
		    this->LeftRightAddedBlockGlobalIndex[Index + k2];
		  this->InteractionIndex[Index + k2][this->NbrInteraction[Index + k2]] = 
		    this->LeftRightAddedBlockGlobalIndex[Index + k1];
		  this->InteractionFactor[Index + k1][this->NbrInteraction[Index + k1]] = 
		    (*(this->Interaction))(this->LeftRightAddedBlockIndex[Index + k1], 
					   this->LeftRightAddedBlockIndex[Index + k2]);
		  this->InteractionFactor[Index + k2][this->NbrInteraction[Index + k2]] = 
		    this->InteractionFactor[Index + k1][this->NbrInteraction[Index + k1]];
		  this->NbrInteraction[Index + k1]++;
		  this->NbrInteraction[Index + k2]++;
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

RealVector& NonPeriodicDMRGHamiltonian::InteractionAddMultiply(RealVector& vSource, 
							       RealVector& vDestination,
							       int firstComponent, int NbrComponent) 
{
  int TotalDim = firstComponent + NbrComponent;
  double* TmpInteractionFactor;
  int* TmpInteractionIndex;
  double Coef;
  int LocalNbrInteraction;
  for (int i = firstComponent; i < TotalDim; i++)
    {
      TmpInteractionFactor = this->InteractionFactor[i];
      TmpInteractionIndex = this->InteractionIndex[i];
      LocalNbrInteraction = this->NbrInteraction[i];
      Coef = 0.0;
      for (int j = 0; j < LocalNbrInteraction; j++)
	{
	  Coef += TmpInteractionFactor[j] * vSource[TmpInteractionIndex[j]];
	}
      vDestination[this->LeftRightAddedBlockGlobalIndex[i]] += Coef;
    }
  return vDestination;
}

