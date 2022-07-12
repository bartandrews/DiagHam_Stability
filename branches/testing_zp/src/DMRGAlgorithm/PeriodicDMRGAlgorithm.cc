////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of periodic DMRG algorithm                     //
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


#include "DMRGAlgorithm/PeriodicDMRGAlgorithm.h"
#include "DMRGAlgorithm/DMRGBlock.h"
#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/DMRGHamiltonian/PeriodicDMRGHamiltonian.h"
#include "Matrix/Matrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace//DMRGHilbertSpace/DMRGPartialHilbertSpace.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "TensorProduct/FullTensorProductStructure.h"
#include "TensorProduct/CompositeTensorProductStructure.h"
#include "Tensor/TwoSpaceTensor.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "Interaction/BasicInteraction.h"

#include <math.h>


using std::cout;
using std::endl;


void TestPeriodicDiagonalize(TwoSpaceTensor* T);

// constructor from datas
//
// blockHamiltonian = Hamiltonian associated to left and right blocks
// interactionBlockHamiltonian = Hamiltonian associated to interaction blocks 
// leftInteraction = interaction between left block and left interaction block
// rightInteraction = interaction between right block and right interaction block
// lanczosAlgorithm = Lanczos algorithm to use
// hilbertSpaceSize = number of states kept for a block

PeriodicDMRGAlgorithm::PeriodicDMRGAlgorithm (AbstractHamiltonian* blockHamiltonian, 
					      AbstractHamiltonian* interactionBlockHamiltonian,
					      AbstractInteraction* leftInteraction,
					      AbstractInteraction* rightInteraction,
					      AbstractLanczosAlgorithm* lanczosAlgorithm,
					      int hilbertSpaceSize)
{
  this->BlockHamiltonian = blockHamiltonian;
  this->InteractionBlockHamiltonian = interactionBlockHamiltonian;
  this->LeftInteraction = leftInteraction;
  this->RightInteraction = rightInteraction;
  this->LanczosAlgorithm = lanczosAlgorithm;
  this->HilbertSpaceSize = hilbertSpaceSize;
  this->Blocks += new DMRGBlock (this->ReduceHamiltonian(this->BlockHamiltonian, 
							 this->BlockHamiltonian->
							 GetHilbertSpaceDimension()), 
				 0);
  this->ExplicitInteractionBlockHamiltonian = this->ReduceHamiltonian(this->InteractionBlockHamiltonian,
								      this->InteractionBlockHamiltonian->GetHilbertSpaceDimension());

  this->GlobalQuantumNumberConstraint = false;
}

// destructor
//

PeriodicDMRGAlgorithm::~PeriodicDMRGAlgorithm () 
{
}

// force constraint on global quantum number
//
// quantumNumber = pointer to the global quantum number to use

void PeriodicDMRGAlgorithm::Constraint(AbstractQuantumNumber* quantumNumber)
{
  this->GlobalQuantumNumber = quantumNumber;
  this->GlobalQuantumNumberConstraint = true;
}

// run DMRG algorithm
//

void PeriodicDMRGAlgorithm::RunDMRG(int currentBlockIndex) 
{
  // construct explicit interaction matrix between left added block and right previous block

/*  TensorProductStructure* InteractionTensorStructure = new TensorProductStructure(2);
  InteractionTensorStructure->SetDimension(0, this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpaceDimension());
  InteractionTensorStructure->SetDimension(1, this->InteractionBlockHamiltonian->GetHilbertSpaceDimension());
  this->RightInteraction->SetTensorProductStructure(InteractionTensorStructure);
  this->RightInteraction->SetLeftSpaceIndex(0);
  this->RightInteraction->SetRightSpaceIndex(1);
  List<Matrix*> TmpLeftInteraction = this->ExplicitInteractionBlockHamiltonian->RightInteractionOperators();
  List<Matrix*> TmpRightInteraction = this->Blocks[currentBlockIndex]->Hamiltonian->LeftInteractionOperators();
  TwoSpaceTensor TmpInteraction = this->RightInteraction->Interaction(TmpRightInteraction, TmpLeftInteraction);
  this->InteractionRightLeftBlocks = TmpInteraction.ElementaryMatrix->Clone();*/

 // this->HilbertSpaceSize = this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpaceDimension() 
//    * this->InteractionBlockHamiltonian->GetHilbertSpaceDimension();

  TensorProductStructure* InteractionTensorStructure = new TensorProductStructure(2);
  InteractionTensorStructure->SetDimension(1, this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpaceDimension());
  InteractionTensorStructure->SetDimension(0, this->InteractionBlockHamiltonian->GetHilbertSpaceDimension());
  this->RightInteraction->SetTensorProductStructure(InteractionTensorStructure);
  this->RightInteraction->SetLeftSpaceIndex(0);
  this->RightInteraction->SetRightSpaceIndex(1);
  List<Matrix*> TmpLeftInteraction = this->ExplicitInteractionBlockHamiltonian->RightInteractionOperators();
  List<Matrix*> TmpRightInteraction = this->Blocks[currentBlockIndex]->Hamiltonian->LeftInteractionOperators();
  TwoSpaceTensor TmpInteraction = this->RightInteraction->Interaction(TmpLeftInteraction, TmpRightInteraction);
  this->InteractionRightLeftBlocks = TmpInteraction.ElementaryMatrix->Clone();

  // construction of Hamiltonian corresponding to left main block and left added block
  //

  // find all possible quantum numbers for the total left block
  //

  List<AbstractQuantumNumber*> TmpListQ1 = this->Blocks[currentBlockIndex]->Hamiltonian->
    GetHilbertSpace()->GetQuantumNumbers();
  List<AbstractQuantumNumber*> TmpListQ2 = this->ExplicitInteractionBlockHamiltonian->
    GetHilbertSpace()->GetQuantumNumbers();
  ListIterator<AbstractQuantumNumber*> IterListQ1(TmpListQ1);
  int Lim = (TmpListQ1.GetNbrElement() * TmpListQ2.GetNbrElement());
  AbstractQuantumNumber** SumQ = new AbstractQuantumNumber* [TmpListQ1.GetNbrElement() * 
							    TmpListQ2.GetNbrElement()];
  AbstractQuantumNumber** TmpQ1;
  AbstractQuantumNumber** TmpQ2;
  int Pos = 0;
  while ((TmpQ1 = IterListQ1()))
    {
      ListIterator<AbstractQuantumNumber*> IterListQ2(TmpListQ2);
      while ((TmpQ2 = IterListQ2()))
	{
	  SumQ[Pos++] = (*TmpQ1)->Add(**TmpQ2);
	}      
    }

  //find decomposition into subspace with respect to the previously found quantum numbers
  //

  int NbrSubspace = 0;
  int* NbrDirectSumSubspace = new int [Lim];
  int** DirectSumSubspaces = new int* [Lim];
  int* DimensionDirectSumSubspace = new int [Lim];
  int NbrSumDeleted = 0;
  for (int i = 0; i < Lim; i++)
    {
      int NbrDirectSum = 0;
      DirectSumSubspaces[NbrSubspace] = new int [2 * (Lim - NbrSumDeleted)];
      NbrDirectSumSubspace[NbrSubspace] = 0;
      if (SumQ[i] != 0)
	{
	  DimensionDirectSumSubspace[NbrSubspace] = 0; 
	  for (int j = i + 1; j < Lim; j++)
	    if ((SumQ[j] != 0) && (SumQ[j]->IsEqual(*(SumQ[i]))))
	      {
		delete SumQ[j];
		SumQ[j] = 0;
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = j / TmpListQ2.GetNbrElement();
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
		  j - TmpListQ2.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
		DimensionDirectSumSubspace[NbrSubspace] += 
		  (((DMRGPartialHilbertSpace*) this->Blocks[currentBlockIndex]->
		    Hamiltonian->GetHilbertSpace())->GetSubspaceDescription
		   (DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum])).GetSubspaceDimension() *
		  (((DMRGPartialHilbertSpace*) this->ExplicitInteractionBlockHamiltonian->
		    GetHilbertSpace())->GetSubspaceDescription
		   (DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1])).GetSubspaceDimension();
		NbrDirectSumSubspace[NbrSubspace]++;
		NbrDirectSum++;
	      }
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = i / TmpListQ2.GetNbrElement();
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
	    i - TmpListQ2.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
	  SumQ[NbrSubspace] = SumQ[i];
	  if (i != NbrSubspace)
	    SumQ[i] = 0;
	  DimensionDirectSumSubspace[NbrSubspace] += (((DMRGPartialHilbertSpace*) 
						       this->Blocks[currentBlockIndex]->
						       Hamiltonian->GetHilbertSpace())->
						      GetSubspaceDescription
						      (DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum])).
	    GetSubspaceDimension() * (((DMRGPartialHilbertSpace*) 
				       this->ExplicitInteractionBlockHamiltonian->
				       GetHilbertSpace())->GetSubspaceDescription
				      (DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1])).
	    GetSubspaceDimension();
	  NbrDirectSumSubspace[NbrSubspace++]++;
	  NbrDirectSum++;
	  NbrSumDeleted += NbrDirectSum;
	}
    }

  // construct hamiltonian on each subspace of total left block
  //

  List<ExplicitHamiltonian*> LeftHamiltonians;
  int** SpacePermutation =  new int* [this->Blocks[currentBlockIndex]->Hamiltonian->
				     GetHilbertSpaceDimension()];
  for (int i = 0; i < this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpaceDimension(); i++)
    SpacePermutation[i] = new int [this->ExplicitInteractionBlockHamiltonian->GetHilbertSpaceDimension()];
  int GlobalPosition = 0;
  int* GlobalSubspacePosition = new int [NbrSubspace];
  GlobalSubspacePosition[0] = 0;
  for (int i = 0; i < NbrSubspace; i++)
    {
      List<SubspaceSpaceConverter> Space1Subspaces;
      List<SubspaceSpaceConverter> Space2Subspaces;
      int* Space1DecompositionArray = new int [NbrDirectSumSubspace[i]];
      int* Space2DecompositionArray = new int [NbrDirectSumSubspace[i]];
      Space1DecompositionArray[0] = 0;
      Space2DecompositionArray[0] = 0;
      FullTensorProductStructure* FineStructure = new FullTensorProductStructure(2, DimensionDirectSumSubspace[i]);
      int TmpState = 0;
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  // find subspace decomposition
	  Space1Subspaces += ((DMRGPartialHilbertSpace*) this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpace())->
	    GetSubspaceDescription(DirectSumSubspaces[i][2 * j]);
	  Space2Subspaces += ((DMRGPartialHilbertSpace*) this->ExplicitInteractionBlockHamiltonian->GetHilbertSpace())->
	    GetSubspaceDescription(DirectSumSubspaces[i][2 * j + 1]);
	  if (j < (NbrDirectSumSubspace[i] - 1))
	    {
	      Space1DecompositionArray[j + 1] = Space1DecompositionArray[j] + 
		((DMRGPartialHilbertSpace*) this->Blocks[currentBlockIndex]->Hamiltonian->GetHilbertSpace())->
		GetSubspaceDescription(DirectSumSubspaces[i][2 * j]).GetSubspaceDimension();
	      Space2DecompositionArray[j + 1] = Space2DecompositionArray[j] + 
		((DMRGPartialHilbertSpace*) this->ExplicitInteractionBlockHamiltonian->GetHilbertSpace())->
		GetSubspaceDescription(DirectSumSubspaces[i][2 * j + 1]).GetSubspaceDimension();
	    }
	  // find permutation of total space indices
	  int k1Lim = ((DMRGPartialHilbertSpace*) this->Blocks[currentBlockIndex]->Hamiltonian->
		       GetHilbertSpace())->GetSubspaceDescription(DirectSumSubspaces[i][2 * j]).
	    GetSubspaceDimension();
	  int k2Lim = ((DMRGPartialHilbertSpace*) this->ExplicitInteractionBlockHamiltonian->
		       GetHilbertSpace())->GetSubspaceDescription(DirectSumSubspaces[i][2 * j + 1]).
	    GetSubspaceDimension();
	  for (int k2 = 0; k2 < k2Lim; k2++)
	    { 
	      int StateIndex2 = ((DMRGPartialHilbertSpace*) this->ExplicitInteractionBlockHamiltonian->
				 GetHilbertSpace())->GetSubspaceDescription(DirectSumSubspaces[i]
									    [2 * j + 1]).GetSpaceIndex(k2);
	      for (int k1 = 0; k1 < k1Lim; k1++)
		{
		  SpacePermutation[((DMRGPartialHilbertSpace*) this->Blocks[currentBlockIndex]->
				    Hamiltonian->GetHilbertSpace())->
				  GetSubspaceDescription(DirectSumSubspaces[i][2 * j]).GetSpaceIndex(k1)]
		    [StateIndex2] = GlobalPosition++;
		  (*FineStructure)(TmpState, 0) = ((DMRGPartialHilbertSpace*) 
						   this->Blocks[currentBlockIndex]->Hamiltonian->
						   GetHilbertSpace())->
		    GetSubspaceDescription(DirectSumSubspaces[i][2 * j]).GetSpaceIndex(k1);
		  (*FineStructure)(TmpState++, 1) = StateIndex2;
		}	      
	    }
	}
      SubspaceSpaceConverter Space1SubspaceDescription (Space1Subspaces);      
      SubspaceSpaceConverter Space2SubspaceDescription (Space2Subspaces);
      TensorProductStructure* Structure = new TensorProductStructure(2);
      Structure->SetDimension(0, Space1SubspaceDescription.GetSubspaceDimension());
      Structure->SetDimension(1, Space2SubspaceDescription.GetSubspaceDimension());
      FineStructure->SetDimension(0, Space1SubspaceDescription.GetSubspaceDimension());
      FineStructure->SetDimension(1, Space2SubspaceDescription.GetSubspaceDimension());
      if (i < (NbrSubspace - 1))  
	{
	  GlobalSubspacePosition[i + 1] = GlobalSubspacePosition[i] + 
	    Space1SubspaceDescription.GetSubspaceDimension();
	}

      // transform boundary operators for total left block
      List<Matrix*> InteractionOperators1;
      List<Matrix*> InteractionOperators2;
      List<Matrix*> TmpInteractionOperators = (this->Blocks[currentBlockIndex]->Hamiltonian->
					       RightInteractionOperators()); 
      ListIterator<Matrix*> IterInteractionOperators(TmpInteractionOperators);
      Matrix** TmpInteraction;
      while ((TmpInteraction = IterInteractionOperators()))
	{
	  InteractionOperators1 += (*TmpInteraction)->Project(Space1SubspaceDescription);
	}
      TmpInteractionOperators = (this->ExplicitInteractionBlockHamiltonian->
				 LeftInteractionOperators());
      IterInteractionOperators.DefineList(TmpInteractionOperators);
      while ((TmpInteraction = IterInteractionOperators()))
	{
	  InteractionOperators2 += (*TmpInteraction)->Project(Space2SubspaceDescription);
	}

      // evaluate interaction part of the total left block
      SpaceDecomposition Space1Decomposition (InteractionOperators1[0]->GetNbrRow(), 
					      NbrDirectSumSubspace[i], Space1DecompositionArray);
      SpaceDecomposition Space2Decomposition (InteractionOperators2[0]->GetNbrRow(), 
					      NbrDirectSumSubspace[i], Space2DecompositionArray);
      this->LeftInteraction->SetTensorProductStructure(FineStructure);
      this->LeftInteraction->SetLeftSpaceIndex(0);
      this->LeftInteraction->SetRightSpaceIndex(1);
      TwoSpaceTensor* Tensor2Left = (TwoSpaceTensor*) this->LeftInteraction->
	Interaction(InteractionOperators1, InteractionOperators2, 
		    Space1Decomposition, Space2Decomposition).Clone();

      Matrix* TmpBlockHamiltonian = (this->Blocks[currentBlockIndex]->Hamiltonian->
				     GetHamiltonian()->
				     Project(Space1SubspaceDescription));
      Matrix* TmpInteractionHamiltonian = (this->ExplicitInteractionBlockHamiltonian->
					   GetHamiltonian()->
					   Project(Space2SubspaceDescription));
      
      // add hamiltonian contributions of initial left block and added block
      Tensor2Left->AddToFirstSpace(*(RealDiagonalMatrix*)(TmpBlockHamiltonian), 
				   Space1Decomposition, Space2Decomposition);
      Tensor2Left->AddToSecondSpace(*(RealDiagonalMatrix*)(TmpInteractionHamiltonian), 
				    Space1Decomposition, Space2Decomposition);
      

      delete TmpBlockHamiltonian;
      delete TmpInteractionHamiltonian;

      // store new hamiltonian
      UndescribedHilbertSpace* TmpHilbertSpace = 
	new UndescribedHilbertSpace(((Matrix*) Tensor2Left)->GetNbrRow(), (*SumQ[i]));
//      TestPeriodicDiagonalize(Tensor2Left);
      LeftHamiltonians += new ExplicitHamiltonian (TmpHilbertSpace, Tensor2Left, InteractionOperators1, InteractionOperators2);
    }

//  return;

  // evaluate total left space decomposition
  SpaceDecomposition TotalLeftSpaceDecomposition(this->Blocks[currentBlockIndex]->Hamiltonian->
						 GetHilbertSpaceDimension() * 
						 this->ExplicitInteractionBlockHamiltonian->
						 GetHilbertSpaceDimension(), 
						 NbrSubspace, GlobalSubspacePosition);

  // free memory  
  for (int i = 0; i < Lim; i++)
    if (SumQ[i] != 0)
      delete SumQ[i];
  delete[] SumQ;
  delete[] NbrDirectSumSubspace;
  delete[] DimensionDirectSumSubspace;
  for (int i = 0; i < NbrSubspace; i++)
    delete[] DirectSumSubspaces[i];
  delete[] DirectSumSubspaces;

  // construct total system
  //

  // find all possible quantun numbers for the total system
  //
  Lim = (LeftHamiltonians.GetNbrElement() * LeftHamiltonians.GetNbrElement());
  SumQ = new AbstractQuantumNumber* [Lim];
  Pos = 0;
  ListIterator<ExplicitHamiltonian*> IterListH1(LeftHamiltonians);
  ExplicitHamiltonian** TmpH1;
  ExplicitHamiltonian** TmpH2;
  while ((TmpH1 = IterListH1()))
    {
      ListIterator<ExplicitHamiltonian*> IterListH2(LeftHamiltonians);
      while ((TmpH2 = IterListH2()))
	{
	  SumQ[Pos++] = (((*TmpH1)->GetHilbertSpace())->GetQuantumNumber(0))->Add(*(((*TmpH2)->GetHilbertSpace())->GetQuantumNumber(0)));
	}      
    }
  NbrSubspace = 0;
  NbrDirectSumSubspace = new int [Lim];
  DirectSumSubspaces = new int* [Lim];
  NbrSumDeleted = 0;

  for (int i = 0; i < Lim; i++)
    {
      if ((SumQ[i] != 0) && ((GlobalQuantumNumberConstraint == false) || 
			     ((GlobalQuantumNumberConstraint == true) &&
			      (this->GlobalQuantumNumber->IsEqual(*(SumQ[i]))))))
	{	  
	  int NbrDirectSum = 0;
	  DirectSumSubspaces[NbrSubspace] = new int [2 * (Lim - NbrSumDeleted)];
	  NbrDirectSumSubspace[NbrSubspace] = 0;
	  for (int j = i + 1; j < Lim; j++)
	    if ((SumQ[j] != 0) && (SumQ[j]->IsEqual(*(SumQ[i]))))
	      {
		delete SumQ[j];
		SumQ[j] = 0;
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = j / LeftHamiltonians.GetNbrElement();
		DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
		  j - LeftHamiltonians.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
		NbrDirectSumSubspace[NbrSubspace]++;
		NbrDirectSum++;
	      }
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum] = i / LeftHamiltonians.GetNbrElement();
	  DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum + 1] = 
	    i - LeftHamiltonians.GetNbrElement() * DirectSumSubspaces[NbrSubspace][2 * NbrDirectSum];
	  NbrDirectSumSubspace[NbrSubspace++]++;
	  NbrDirectSum++;
	  NbrSumDeleted += NbrDirectSum;
	}
    }

//  return;

  List<ExplicitHamiltonian*> Hamiltonians2;
  for (int i = 0; i < NbrSubspace; i++)
    {
      List<Matrix*> RightHamiltonianRepresentationList;
      List<Matrix*> LeftHamiltonianRepresentationList;
      List<TensorProductStructure> SpaceStructure;
      List<SubspaceSpaceConverter> RightSpaceSubspaces;
      List<SubspaceSpaceConverter> LeftSpaceSubspaces;
      int* LeftSpaceDecompositionArray = new int [NbrDirectSumSubspace[i]];
      int* RightSpaceDecompositionArray = new int [NbrDirectSumSubspace[i]];
      FullTensorProductStructure** LeftSpaceFineStructureArray = 
	new FullTensorProductStructure* [NbrDirectSumSubspace[i]];
      FullTensorProductStructure** RightSpaceFineStructureArray = 
	new FullTensorProductStructure* [NbrDirectSumSubspace[i]];
      LeftSpaceDecompositionArray[0] = 0;
      RightSpaceDecompositionArray[0] = 0;
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  LeftSpaceFineStructureArray[j] = (FullTensorProductStructure*) 
	    ((TwoSpaceTensor*) (LeftHamiltonians[DirectSumSubspaces[i][2 * j]]->
				GetHamiltonian()))->GetTensorProductStructure();
	  RightSpaceFineStructureArray[j] = (FullTensorProductStructure*) 
	    ((TwoSpaceTensor*) (LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->
				GetHamiltonian()))->GetTensorProductStructure();
	  RightSpaceSubspaces += SubspaceSpaceConverter(TotalLeftSpaceDecomposition.GetSpaceDimension(), 
							LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->
							GetHilbertSpaceDimension(),
							TotalLeftSpaceDecomposition.
							GetSubspacePosition(DirectSumSubspaces[i]
									    [2 * j + 1]));
	  LeftSpaceSubspaces += SubspaceSpaceConverter(TotalLeftSpaceDecomposition.GetSpaceDimension(), 
						       LeftHamiltonians[DirectSumSubspaces[i][2 * j]]->
						       GetHilbertSpaceDimension(),
						       TotalLeftSpaceDecomposition.
						       GetSubspacePosition(DirectSumSubspaces[i]
									   [2 * j]));
	  RightHamiltonianRepresentationList += ((TwoSpaceTensor*) (LeftHamiltonians
								    [DirectSumSubspaces[i][2 * j + 1]]->
								    GetHamiltonian()))->
	    ElementaryMatrix->Clone();
	  LeftHamiltonianRepresentationList += 
	    ((TwoSpaceTensor*) (LeftHamiltonians[DirectSumSubspaces[i][2 * j]]->GetHamiltonian()))->
	    ElementaryMatrix->Clone();

	  if (j < (NbrDirectSumSubspace[i] - 1))
	    {
	      RightSpaceDecompositionArray[j + 1] = RightSpaceDecompositionArray[j] +
		LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->GetHilbertSpaceDimension();
	      LeftSpaceDecompositionArray[j + 1] = LeftSpaceDecompositionArray[j] + 
		LeftHamiltonians[DirectSumSubspaces[i][2 * j]]->GetHilbertSpaceDimension();
	    }
	  TensorProductStructure TmpStructure (2);
	  TmpStructure.SetDimension(0, LeftHamiltonians[DirectSumSubspaces[i][2 * j]]->
				    GetHamiltonian()->GetNbrRow());
	  TmpStructure.SetDimension(1, LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->
				    GetHamiltonian()->GetNbrRow());
	  SpaceStructure += TmpStructure;
	}
      SubspaceSpaceConverter LeftSpaceSubspaceDescription (LeftSpaceSubspaces);      
      SubspaceSpaceConverter RightSpaceSubspaceDescription (RightSpaceSubspaces);

      TensorProductStructure* Structure = new TensorProductStructure(2);
      Structure->SetDimension(0, LeftSpaceSubspaceDescription.GetSubspaceDimension());
      Structure->SetDimension(1, RightSpaceSubspaceDescription.GetSubspaceDimension());
      SpaceDecomposition LeftSpaceDecomposition (LeftSpaceSubspaceDescription.GetSubspaceDimension(), 
						 NbrDirectSumSubspace[i], LeftSpaceDecompositionArray);
      SpaceDecomposition RightSpaceDecomposition (RightSpaceSubspaceDescription.GetSubspaceDimension(), 
						   NbrDirectSumSubspace[i], RightSpaceDecompositionArray);

      CompositeTensorProductStructure* SpaceTensorProductStructure = 
	new CompositeTensorProductStructure (SpaceStructure);
      BlockDiagonalMatrix LeftMatrix (LeftHamiltonianRepresentationList);
      OneSpaceTensor LeftTensor(SpaceTensorProductStructure, LeftMatrix, 0);
      BlockDiagonalMatrix RightMatrix (RightHamiltonianRepresentationList);
      OneSpaceTensor RightTensor(SpaceTensorProductStructure, RightMatrix, 1);
      PeriodicDMRGHamiltonian TotalHamiltonian (LeftTensor, RightTensor,
						*(this->InteractionRightLeftBlocks),
						SpaceTensorProductStructure, LeftSpaceFineStructureArray,
						RightSpaceFineStructureArray, 
						this->Blocks[currentBlockIndex]->Hamiltonian->
						GetHilbertSpaceDimension(), 
						this->ExplicitInteractionBlockHamiltonian->
						GetHilbertSpaceDimension());

/*      RealSymmetricMatrix HRep (TotalHamiltonian.GetHilbertSpaceDimension());
      RealTriDiagonalSymmetricMatrix TmpTriDiag (TotalHamiltonian.GetHilbertSpaceDimension());
      TotalHamiltonian.GetHamiltonian(HRep);
      cout << HRep << endl;
      HRep.Householder(TmpTriDiag, 1e-14);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      for (int j = 0; j < TotalHamiltonian.GetHilbertSpaceDimension(); j++)
	cout << TmpTriDiag.DiagonalElement(j) << " ";
      cout << endl;*/
      

      this->LanczosAlgorithm->SetHamiltonian(&TotalHamiltonian);
      this->LanczosAlgorithm->InitializeLanczosAlgorithm();
      int MaxNbrIterLanczos = 200; 
      this->LanczosAlgorithm->RunLanczosAlgorithm(4);
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      this->GroundStateEnergy = PreviousLowest;
      this->NbrLanczosIteration = 4;
      while ((Precision > MACHINE_PRECISION) && (this->NbrLanczosIteration++ < MaxNbrIterLanczos))
	{
	  this->LanczosAlgorithm->RunLanczosAlgorithm(1);
	  this->GroundStateEnergy = this->LanczosAlgorithm->GetGroundStateEnergy();
	  Precision = fabs((PreviousLowest - this->GroundStateEnergy) / PreviousLowest);
//	  cout << this->NbrLanczosIteration << " " << this->GroundStateEnergy << " " << Precision << endl;
	  PreviousLowest = this->GroundStateEnergy;
	}
//      cout << "end Lanczos" << endl;      
      // evaluate density matrix
      RealVector& GroundState = (RealVector&) this->LanczosAlgorithm->GetGroundState();
      RealMatrix* DensityMatrixEigenvectors = new RealMatrix [NbrDirectSumSubspace[i]];
      RealTriDiagonalSymmetricMatrix* DensityMatrixEigenvalues = 
	new RealTriDiagonalSymmetricMatrix [NbrDirectSumSubspace[i]];
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  // evaluate density matrix on a subspace
	  int TmpDim1 = SpaceTensorProductStructure->GetDimension(1, j);
	  int TmpDim2 = SpaceTensorProductStructure->GetDimension(0, j);
	  int Pos11 = SpaceTensorProductStructure->GetSubspaceIncrement(j);
	  int Pos12;
	  int Pos21 = SpaceTensorProductStructure->GetSubspaceIncrement(j);
	  int Pos22;
	  double Coef;
	  RealSymmetricMatrix DensityMatrix(TmpDim1);
	  for (int k1 = 0; k1 < TmpDim1; k1++)
	    {		  
	      Pos22 = Pos21 + k1 * TmpDim2;
	      for (int k2 = k1; k2 < TmpDim1; k2++)
		{
		  Pos12 = Pos11;
		  Coef = 0.0;
		  for (int k3 = 0; k3 < TmpDim2; k3++)
		    Coef += GroundState[Pos12++] * GroundState[Pos22++];
		  DensityMatrix(k1, k2) = Coef;
		}
	      Pos11 += TmpDim2;
	    }
	  // diagonalize density matrix
	  if (TmpDim1 > 1)
	    {
	      DensityMatrixEigenvalues[j] = RealTriDiagonalSymmetricMatrix (TmpDim1);
	      DensityMatrixEigenvectors[j] = RealMatrix (TmpDim1, TmpDim1);
	      DensityMatrix.Householder(DensityMatrixEigenvalues[j], MACHINE_PRECISION, 
					DensityMatrixEigenvectors[j]);
	      DensityMatrixEigenvalues[j].Diagonalize(DensityMatrixEigenvectors[j]);
	      DensityMatrixEigenvalues[j].SortMatrixDownOrder(DensityMatrixEigenvectors[j]);
	    }
	  else
	    {
	      DensityMatrixEigenvalues[j] = RealTriDiagonalSymmetricMatrix (1);
	      DensityMatrixEigenvalues[j].DiagonalElement(0) = DensityMatrix(0, 0);
	      DensityMatrixEigenvectors[j] = RealMatrix (1, 1);
	      DensityMatrixEigenvectors[j](0, 0) = 1.0;
	    }		
	}
//      cout << "end density matrix diagonalization" << endl;  

      // sort and truncate density matrix
      int* DensityMatrixKeptEigenvalues = new int [NbrDirectSumSubspace[i]];
      bool* DensityMatrixEigenvaluesCompleteFlag = new bool [NbrDirectSumSubspace[i]];
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  DensityMatrixKeptEigenvalues[j] = 0;
	  DensityMatrixEigenvaluesCompleteFlag[j] = false;
	}
      int CurrentTruncationDimension = 0;
      int MinPos;
      double MinValue;	  
      while (CurrentTruncationDimension < this->HilbertSpaceSize)
	{
	  MinPos = 0;
	  while (DensityMatrixEigenvaluesCompleteFlag[MinPos] == true)
	    MinPos++;
	  MinValue =  DensityMatrixEigenvalues[MinPos].
	    DiagonalElement(DensityMatrixKeptEigenvalues[MinPos]);
	  for (int j = MinPos + 1; j < NbrDirectSumSubspace[i]; j++)
	    if ((DensityMatrixEigenvaluesCompleteFlag[j] == false) && 
		(MinValue < DensityMatrixEigenvalues[j].DiagonalElement(DensityMatrixKeptEigenvalues[j])))
	      {
		MinValue = DensityMatrixEigenvalues[j].DiagonalElement(DensityMatrixKeptEigenvalues[j]);
		MinPos = j;						   
	      }	      
	  DensityMatrixKeptEigenvalues[MinPos]++;
	  if (DensityMatrixKeptEigenvalues[MinPos] == DensityMatrixEigenvalues[MinPos].GetNbrRow())
	    DensityMatrixEigenvaluesCompleteFlag[MinPos] = true;
	  CurrentTruncationDimension++;
	}
      this->TruncationError = 0.0;
      MinPos = 0;
      int NbrLeftSpace = 0;
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  if (DensityMatrixKeptEigenvalues[j] > 0)
	    {
	      for (int k = 0; k < DensityMatrixKeptEigenvalues[j]; k++)
		{
		  this->TruncationError += DensityMatrixEigenvalues[j].DiagonalElement(k);
		  MinPos++;		      
		}
	      DensityMatrixEigenvectors[j].Resize(DensityMatrixEigenvectors[j].GetNbrRow(), 
						  DensityMatrixKeptEigenvalues[j]);
	      NbrLeftSpace++;
	    }
	}
      this->TruncationError = 1.0 - this->TruncationError; 
      
//       cout << "end truncation" << endl;      
     // conjugate total left block hamiltonian and associated interactions
      int* TruncatedSubspacePosition = new int [NbrLeftSpace];
      int ColumnPos = 0;
      List<AbstractQuantumNumber*> ListLeftQuantumNumber;
      List<Matrix*> ListTransformationDensityMatrix;
      int BlockPos = 0;
      Pos = 0;
      RealDiagonalMatrix* TruncatedHamiltonian = new RealDiagonalMatrix(this->HilbertSpaceSize);
      for (int j = 0; j < NbrDirectSumSubspace[i]; j++)
	{
	  if (DensityMatrixKeptEigenvalues[j] > 0)
	    {		  
	      TruncatedSubspacePosition[BlockPos++] = ColumnPos;
	      ColumnPos += DensityMatrixKeptEigenvalues[j];
	      ListLeftQuantumNumber += LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->
		GetHilbertSpace()->GetQuantumNumber(0);
	      Matrix* TmpHamiltonian = 
		((RealSymmetricMatrix*) ((TwoSpaceTensor*) LeftHamiltonians[DirectSumSubspaces[i][2 * j + 1]]->
					 GetHamiltonian())->ElementaryMatrix)->
		Conjugate(DensityMatrixEigenvectors[j]);
	      if (TmpHamiltonian->GetNbrRow() > 1)
		{
		  RealMatrix TmpEigenvectors (TmpHamiltonian->GetNbrRow(), TmpHamiltonian->GetNbrColumn());
		  
		  RealTriDiagonalSymmetricMatrix TmpEigenvalues (TmpHamiltonian->GetNbrRow());
		  ((RealSymmetricMatrix*) TmpHamiltonian)->Householder(TmpEigenvalues, MACHINE_PRECISION, 
								       TmpEigenvectors);
		  TmpEigenvalues.Diagonalize(TmpEigenvectors);
		  DensityMatrixEigenvectors[j] = DensityMatrixEigenvectors[j] * TmpEigenvectors;
		  ListTransformationDensityMatrix += (Matrix*) &(DensityMatrixEigenvectors[j]);
		  for (int k = 0; k < TmpEigenvalues.GetNbrRow(); k++)
		    {
		      (*TruncatedHamiltonian)[Pos++] = TmpEigenvalues.DiagonalElement(k);
		    }
		}
	      else
		{
		  RealMatrix TmpEigenvectors (1, 1);
		  TmpEigenvectors(0, 0) = 1.0;
		  ListTransformationDensityMatrix += (Matrix*) &(DensityMatrixEigenvectors[j]);
		  (*TruncatedHamiltonian)[Pos++] = (*TmpHamiltonian)(0,0);
		}
	      delete TmpHamiltonian;
	    }
	}
//      cout << "end transformation hamiltonian left block" << endl;      

      BlockDiagonalMatrix TransformationMatrix  (ListTransformationDensityMatrix);

      List<Matrix*> TmpListInteraction = this->ExplicitInteractionBlockHamiltonian->RightInteractionOperators();
      List<Matrix*> ConjugatedTotalLeftBlockLeftInteractions = 
	this->TransformOperator(TmpListInteraction, true, TransformationMatrix, 
				NbrDirectSumSubspace[i], DensityMatrixKeptEigenvalues, 
				RightSpaceFineStructureArray, this->Blocks[currentBlockIndex]
				->Hamiltonian->GetHilbertSpaceDimension(), 
				this->ExplicitInteractionBlockHamiltonian->
				GetHilbertSpaceDimension());
      TmpListInteraction = this->Blocks[currentBlockIndex]->Hamiltonian->LeftInteractionOperators();
      List<Matrix*> ConjugatedTotalLeftBlockRightInteractions = 
	this->TransformOperator(TmpListInteraction, false, TransformationMatrix, 
				NbrDirectSumSubspace[i], DensityMatrixKeptEigenvalues, 
				RightSpaceFineStructureArray, this->Blocks[currentBlockIndex]
				->Hamiltonian->GetHilbertSpaceDimension(), 
				this->ExplicitInteractionBlockHamiltonian->
				GetHilbertSpaceDimension());
//      cout << "end transformation left block" << endl;      

      // store truncated hilbert space and hamiltonian of the new left block 
      DMRGPartialHilbertSpace* LeftSpace = 
	new DMRGPartialHilbertSpace(*TruncatedHamiltonian, ListLeftQuantumNumber,
				    SpaceDecomposition(this->HilbertSpaceSize, 
						       ListLeftQuantumNumber.GetNbrElement(),
						       TruncatedSubspacePosition));
      ExplicitHamiltonian* LeftHamiltonian = 
	new ExplicitHamiltonian ((AbstractHilbertSpace*) LeftSpace, TruncatedHamiltonian, 				 
				 ConjugatedTotalLeftBlockRightInteractions,
				 ConjugatedTotalLeftBlockLeftInteractions);
      
      this->Blocks += new DMRGBlock (LeftHamiltonian, currentBlockIndex + 1);


//      cout << "end new block" << endl;      

      // free memory
      delete[] DensityMatrixEigenvectors;
      delete[] DensityMatrixKeptEigenvalues;
      delete[] DensityMatrixEigenvaluesCompleteFlag;
      delete[] DensityMatrixEigenvalues;
      delete[] LeftSpaceFineStructureArray;
      delete[] RightSpaceFineStructureArray;
//      cout << "end first free" << endl;      
    }

  ListIterator<ExplicitHamiltonian*> IterHamiltonians(LeftHamiltonians);
  ExplicitHamiltonian** TmpHamiltonian;
  while ((TmpHamiltonian = IterHamiltonians()))
    delete *TmpHamiltonian;
//  cout << "end second free" << endl;      
  for (int i = 0; i < Lim; i++)
    if (SumQ[i] != 0)
      delete SumQ[i];
  delete[] SumQ;
}

// evaluate an operator after density matrix reduction
//
// observables = reference on a list of operators to transform
// spaceFlag = true if belong to the added block
// transformationMatrix = reference on transformation matrix to use
// nbrSubspace = number of subspaces
// keptSubspace = array containing number of states kept in each subspace
// spaceFineStructureArray = array of fine structure describing total space
// initialBlockSize = size of Hilbert space associated to initial block
// addedBlockSize = size of Hilbert space associated to added block
// return value = list of transformed operator

List<Matrix*> PeriodicDMRGAlgorithm::TransformOperator(List<Matrix*>& observables, bool spaceFlag, 
						       BlockDiagonalMatrix& transformationMatrix, 
						       int nbrSubspace, int* keptSubspace, 
						       FullTensorProductStructure** spaceFineStructureArray,
						       int initialBlockSize, int addedBlockSize)
{
  if (spaceFlag == true)
    {
      int* SubspaceSize = new int [initialBlockSize];
      int** StateIndex = new int* [initialBlockSize];
      int** StateGlobalIndex = new int* [initialBlockSize];
      int GlobalIndexShift = 0;
      int TmpIndex0;
      int TmpIndex1;
      FullTensorProductStructure* TmpStructure;
      for (int i = 0; i < initialBlockSize; i++)
	{
	  SubspaceSize[i] = 0;
	  StateIndex[i] = new int [addedBlockSize];
	  StateGlobalIndex[i] = new int [addedBlockSize];
	}
      for (int i = 0; i < nbrSubspace; i++)
	{
	  if (keptSubspace[i] > 0)
	    {
	      TmpStructure = spaceFineStructureArray[i];
	      for (int j = 0; j < TmpStructure->GetTotalDimension(); j++)
		{
		  TmpIndex0 = (*TmpStructure)(j, 0);
		  TmpIndex1 = (*TmpStructure)(j, 1);
		  StateIndex[TmpIndex0][SubspaceSize[TmpIndex0]] = TmpIndex1;
		  StateGlobalIndex[TmpIndex0][SubspaceSize[TmpIndex0]] = GlobalIndexShift + j;
		  SubspaceSize[TmpIndex0]++;
		}
	      GlobalIndexShift += TmpStructure->GetTotalDimension();
	    }
	}
      int MinValue;
      int MinPos;
      int* TmpStateGlobalIndex;
      int* TmpStateIndex;
      for (int i = 0; i < initialBlockSize; i++)
	if (SubspaceSize[i] > 0)
	  {
	    TmpStateGlobalIndex = StateGlobalIndex[i];
	    TmpStateIndex = StateIndex[i];
	    TmpIndex0 = SubspaceSize[i];
	    TmpIndex1 = TmpIndex0 - 1;
	    for (int j1 = 0; j1 < TmpIndex0; j1++)
	      {
		MinPos = TmpIndex1;
		MinValue = TmpStateIndex[MinPos];
		for (int j2 = TmpIndex1 - 1; j2 >= j1; j2--)
		  if (MinValue > TmpStateIndex[j2])
		    {
		      MinPos = j2;
		      MinValue = TmpStateIndex[j2];
		    }
		TmpStateIndex[MinPos] = TmpStateIndex[j1];
		TmpStateIndex[j1] = MinValue;
		MinValue = TmpStateGlobalIndex[MinPos];
		TmpStateGlobalIndex[MinPos] = TmpStateGlobalIndex[j1];
		TmpStateGlobalIndex[j1] = MinValue;
	      }
	  }
      TensorProductStructure* TmpStructure2 = new TensorProductStructure(2);
      TmpStructure2->SetDimension(0, GlobalIndexShift);
      TmpStructure2->SetDimension(1, GlobalIndexShift);
      List<Matrix*> TransformedObservables;
      ListIterator<Matrix*> IterObservables (observables);
      Matrix** TmpObservable;
      while ((TmpObservable = IterObservables()))
	{
	  TwoSpaceTensor TmpTensor (TmpStructure2, **TmpObservable, GlobalIndexShift, initialBlockSize, 
				    SubspaceSize, StateIndex, StateGlobalIndex, 0);
	  TransformedObservables += TmpTensor.ElementaryMatrix->Conjugate(transformationMatrix);
	}
      delete[] SubspaceSize;
      for (int i = 0; i < initialBlockSize; i++)
	{
	  delete[] StateIndex[i];
	  delete[] StateGlobalIndex[i];
	}
      delete[] StateIndex;
      delete[] StateGlobalIndex;
      return TransformedObservables;
    }
  else
    {
      int* SubspaceSize = new int [addedBlockSize];
      int** StateIndex = new int* [addedBlockSize];
      int** StateGlobalIndex = new int* [addedBlockSize];
      int GlobalIndexShift = 0;
      int TmpIndex0;
      int TmpIndex1;
      FullTensorProductStructure* TmpStructure;
      for (int i = 0; i < addedBlockSize; i++)
	{
	  SubspaceSize[i] = 0;
	  StateIndex[i] = new int [initialBlockSize];
	  StateGlobalIndex[i] = new int [initialBlockSize];
	}
      for (int i = 0; i < nbrSubspace; i++)
	{
	  if (keptSubspace[i] > 0)
	    {
	      TmpStructure = spaceFineStructureArray[i];
	      for (int j = 0; j < TmpStructure->GetTotalDimension(); j++)
		{
		  TmpIndex0 = (*TmpStructure)(j, 0);
		  TmpIndex1 = (*TmpStructure)(j, 1);
		  StateIndex[TmpIndex1][SubspaceSize[TmpIndex1]] = TmpIndex0;
		  StateGlobalIndex[TmpIndex1][SubspaceSize[TmpIndex1]] = GlobalIndexShift + j;
		  SubspaceSize[TmpIndex1]++;
		}
	      GlobalIndexShift += TmpStructure->GetTotalDimension();
	    }
	}
      int MinValue;
      int MinPos;
      int* TmpStateGlobalIndex;
      int* TmpStateIndex;
      for (int i = 0; i < addedBlockSize; i++)
	if (SubspaceSize[i] > 0)
	  {
	    TmpStateGlobalIndex = StateGlobalIndex[i];
	    TmpStateIndex = StateIndex[i];
	    TmpIndex0 = SubspaceSize[i];
	    TmpIndex1 = TmpIndex0 - 1;
	    for (int j1 = 0; j1 < TmpIndex0; j1++)
	      {
		MinPos = TmpIndex1;
		MinValue = TmpStateIndex[MinPos];
		for (int j2 = TmpIndex1 - 1; j2 >= j1; j2--)
		  if (MinValue > TmpStateIndex[j2])
		    {
		      MinPos = j2;
		      MinValue = TmpStateIndex[j2];
		    }
		TmpStateIndex[MinPos] = TmpStateIndex[j1];
		TmpStateIndex[j1] = MinValue;
		MinValue = TmpStateGlobalIndex[MinPos];
		TmpStateGlobalIndex[MinPos] = TmpStateGlobalIndex[j1];
		TmpStateGlobalIndex[j1] = MinValue;
	      }
	  }
      TensorProductStructure* TmpStructure2 = new TensorProductStructure(2);
      TmpStructure2->SetDimension(0, GlobalIndexShift);
      TmpStructure2->SetDimension(1, GlobalIndexShift);
      List<Matrix*> TransformedObservables;
      ListIterator<Matrix*> IterObservables (observables);
      Matrix** TmpObservable;
      while ((TmpObservable = IterObservables()))
	{
	  TwoSpaceTensor TmpTensor (TmpStructure2, **TmpObservable, GlobalIndexShift, addedBlockSize, 
				    SubspaceSize, StateIndex, StateGlobalIndex, 0);
	  TransformedObservables += TmpTensor.ElementaryMatrix->Conjugate(transformationMatrix);
	}
      delete[] SubspaceSize;
      for (int i = 0; i < addedBlockSize; i++)
	{
	  delete[] StateIndex[i];
	  delete[] StateGlobalIndex[i];
	}
      delete[] StateIndex;
      delete[] StateGlobalIndex;
      return TransformedObservables;
    }
}

// reduce hamiltonian to the basis of its n-first eigenstates
//
// hamiltonian = hamiltonian
// truncatedSize = number of states to keep
// return value = reduced hamiltonian (with associated Hilbert space)

ExplicitHamiltonian* PeriodicDMRGAlgorithm::ReduceHamiltonian (AbstractHamiltonian* hamiltonian, 
							       int truncatedSize) 
{
  AbstractHilbertSpace* TmpSpace = hamiltonian->GetHilbertSpace();
  List<SubspaceSpaceConverter> Converters;
  int SubspaceIndex = 0;
  List<AbstractQuantumNumber*> ListQ = TmpSpace->GetQuantumNumbers();
//  cout << hamiltonian->GetHilbertSpaceDimension() << endl;
  ListIterator<AbstractQuantumNumber*> IterQ (ListQ);
  AbstractQuantumNumber** TmpQ;

  // diagonalize Hamiltonian

  RealTriDiagonalSymmetricMatrix* TmpEigenvalues = 
    new RealTriDiagonalSymmetricMatrix [ListQ.GetNbrElement()];
  RealMatrix* TmpEigenvectors = new RealMatrix [ListQ.GetNbrElement()];
//  cout << hamiltonian->GetHilbertSpaceDimension() << endl;
  while ((TmpQ = IterQ()))
    {
//      cout <<**TmpQ << endl;
      //      cout << **TmpQ << endl;
      SubspaceSpaceConverter Converter2;
      AbstractHilbertSpace* SubChain2 = TmpSpace->ExtractSubspace(**TmpQ, Converter2);
//      cout << "dim=" << SubChain2->GetHilbertSpaceDimension() << endl;
      hamiltonian->SetHilbertSpace(SubChain2);
//      for (int  i = 0; i < SubChain2->GetHilbertSpaceDimension(); i++)
//	SubChain2->PrintState(cout ,i) << (*(SubChain2->GetQuantumNumber(i))) << " " << i << " " << 
//	  ((Spin1_2Chain*) SubChain2)->FindStateIndex(((Spin1_2Chain*) SubChain2)->ChainDescription[i]) << endl;
	
      RealSymmetricMatrix HRep2 (SubChain2->GetHilbertSpaceDimension());
      hamiltonian->GetHamiltonian(HRep2);
      if (SubChain2->GetHilbertSpaceDimension() > 1)
	{
	  TmpEigenvectors[SubspaceIndex] = RealMatrix (SubChain2->GetHilbertSpaceDimension(), 
						       SubChain2->GetHilbertSpaceDimension());
	  TmpEigenvalues[SubspaceIndex] = RealTriDiagonalSymmetricMatrix (SubChain2->
									  GetHilbertSpaceDimension());
	  HRep2.Householder(TmpEigenvalues[SubspaceIndex], MACHINE_PRECISION, 
			    TmpEigenvectors[SubspaceIndex]);
	  TmpEigenvalues[SubspaceIndex].Diagonalize(TmpEigenvectors[SubspaceIndex]);
	  TmpEigenvalues[SubspaceIndex].SortMatrixUpOrder(TmpEigenvectors[SubspaceIndex]);
	}
      else
	{
	  TmpEigenvectors[SubspaceIndex] = RealMatrix (1, 1);
	  TmpEigenvalues[SubspaceIndex] = RealTriDiagonalSymmetricMatrix (1);
	  TmpEigenvalues[SubspaceIndex].DiagonalElement(0) = HRep2(0, 0);
	  TmpEigenvectors[SubspaceIndex](0, 0) = 1.0;
	}
      SubspaceIndex++;
      Converters += Converter2;
      delete SubChain2;
    }
  SubspaceSpaceConverter Converter(Converters);

  int Count = 0;
  for (int i = 0; i < SubspaceIndex; i++)
    for (int j = 0; j < TmpEigenvalues[i].GetNbrRow(); j++)
      {
	Count++;
      }
  // truncate Hilbert space

  double* TruncatedSpectrum = new double [truncatedSize];
  int* KeptStates = new int [SubspaceIndex];
  bool* KeptStatesCompleteFlag = new bool [SubspaceIndex];
  for (int i = 0; i < SubspaceIndex; i++)
    {
      KeptStates[i] = 0;
      KeptStatesCompleteFlag[i] = false;
    }  
  int CurrentTruncatedSize = 0;
  int MinPos;
  double MinValue;	  
  while (CurrentTruncatedSize < truncatedSize)
    {
      MinPos = 0;
      while (KeptStatesCompleteFlag[MinPos] == true)
	MinPos++;
      MinValue = TmpEigenvalues[MinPos].DiagonalElement(KeptStates[MinPos]);
      for (int i = MinPos + 1; i < SubspaceIndex; i++)
	if ((KeptStatesCompleteFlag[i] == false) && 
	    (MinValue > TmpEigenvalues[i].DiagonalElement(KeptStates[i])))
	  {
	    MinPos = i;
	    MinValue = TmpEigenvalues[i].DiagonalElement(KeptStates[i]);
	  }
      KeptStates[MinPos]++;
      if (KeptStates[MinPos] == TmpEigenvalues[MinPos].GetNbrRow())
	KeptStatesCompleteFlag[MinPos] = true;
      CurrentTruncatedSize++;
    }

  // evaluate diagonalized and truncated hamiltonian, find all kept quantum numbers 
  // and prepare tranformation matrix
  
  MinPos = 0;
  List<AbstractQuantumNumber*> TruncatedListQuantumNumber;
  List<Matrix*> TransformationMatrixList;
  int TruncatedNbrSubspace = 0;
  for (int i = 0; i < SubspaceIndex; i++)
    {
      for (int j = 0; j < KeptStates[i]; j++)
	{
	  TruncatedSpectrum[MinPos++] = TmpEigenvalues[i].DiagonalElement(j);
	}
      if (KeptStates[i] > 0)
	{
	  TruncatedListQuantumNumber += ListQ[i];
	  TmpEigenvectors[i].Resize(TmpEigenvectors[i].GetNbrRow(), KeptStates[i]);
	  TransformationMatrixList += &(TmpEigenvectors[i]);
	  TruncatedNbrSubspace++;
	}
    }
  BlockDiagonalMatrix TransformationMatrix(TransformationMatrixList);
  RealDiagonalMatrix* TruncatedHamiltonian = new RealDiagonalMatrix (TruncatedSpectrum, truncatedSize);

  // find truncated space structure with respect to quantum numbers

  int* TruncatedSubspacePosition = new int [TruncatedNbrSubspace];
  TruncatedSubspacePosition[0] = 0;
  MinPos = 1;
  for (int i = 1; i < SubspaceIndex; i++)
    if (KeptStates[i] > 0)
      {
	TruncatedSubspacePosition[MinPos] = TruncatedSubspacePosition[MinPos - 1] + KeptStates[i - 1];
	MinPos++;
      }

  // conjugate interaction operators

  hamiltonian->SetHilbertSpace(TmpSpace);
  List<Matrix*> RightInteractionOperators (hamiltonian->RightInteractionOperators());
  List<Matrix*> RightConjugatedInteractionOperators;
  ListIterator<Matrix*> IterOperators(RightInteractionOperators);
  Matrix** TmpMatrix;
  Matrix* TmpMatrix2;
  int i = 0;
  while ((TmpMatrix = IterOperators()))
    {
      TmpMatrix2 = (*TmpMatrix)->Project(Converter);
      RightConjugatedInteractionOperators += TmpMatrix2->Conjugate(TransformationMatrix);
    }
  List<Matrix*> LeftInteractionOperators (hamiltonian->LeftInteractionOperators());
  List<Matrix*> LeftConjugatedInteractionOperators;
  IterOperators.DefineList(LeftInteractionOperators);
  i = 0;
  while ((TmpMatrix = IterOperators()))
    {
      TmpMatrix2 = (*TmpMatrix)->Project(Converter);
      LeftConjugatedInteractionOperators += TmpMatrix2->Conjugate(TransformationMatrix);
    }

  DMRGPartialHilbertSpace* LeftSpace = 
    new DMRGPartialHilbertSpace(*TruncatedHamiltonian, TruncatedListQuantumNumber,
				SpaceDecomposition(truncatedSize, 
						   TruncatedListQuantumNumber.GetNbrElement(),
						   TruncatedSubspacePosition));


  // delete temporay matrices
  delete[] TmpEigenvalues;
  delete[] TmpEigenvectors;
  delete[] KeptStates;
  delete[] KeptStatesCompleteFlag;
  
  return new ExplicitHamiltonian ((AbstractHilbertSpace*) LeftSpace, TruncatedHamiltonian, 
				  RightConjugatedInteractionOperators, LeftConjugatedInteractionOperators);

}

void TestPeriodicDiagonalize(TwoSpaceTensor* T)
{
  RealSymmetricMatrix* HRep = (RealSymmetricMatrix*) (T->ElementaryMatrix);
  cout << "Dimension = " << HRep->GetNbrRow() << endl;
  if (HRep->GetNbrRow() <= 0)
    return;
  if (HRep->GetNbrRow() == 1)
    {
      cout << (*HRep)(0, 0) << endl;
      return;
    }
  RealTriDiagonalSymmetricMatrix D (HRep->GetNbrRow());
  HRep->Householder (D, MACHINE_PRECISION);
  D.Diagonalize();
  D.SortMatrixUpOrder();
  for (int i = 0; i < HRep->GetNbrRow(); i++)
    cout << D.DiagonalElement(i) << "  ";
  cout << endl << " end " << endl;
}
