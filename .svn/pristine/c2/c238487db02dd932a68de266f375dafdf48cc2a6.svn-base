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


#ifndef PERIODICDMRGALGORITHM_H
#define PERIODICDMRGALGORITHM_H


#include "DMRGAlgorithm/AbstractDMRGAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Interaction/AbstractInteraction.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "GeneralTools/List.h"


class DMRGBlock;
class ExplicitHamiltonian;
class Matrix;
class AbstractQuantumNumber;
class BlockDiagonalMatrix;
class FullTensorProductStructure;


class PeriodicDMRGAlgorithm : public AbstractDMRGAlgorithm
{

 protected:

  AbstractHamiltonian* BlockHamiltonian;
  AbstractHamiltonian* InteractionBlockHamiltonian;

  ExplicitHamiltonian* ExplicitInteractionBlockHamiltonian;

  AbstractInteraction* LeftInteraction;
  AbstractInteraction* RightInteraction;
  Matrix* InteractionRightLeftBlocks;

  int HilbertSpaceSize;

  List<DMRGBlock*> Blocks;

  bool GlobalQuantumNumberConstraint;
  AbstractQuantumNumber* GlobalQuantumNumber;

  AbstractLanczosAlgorithm* LanczosAlgorithm;

  double GroundStateEnergy;
  double TruncationError;
  int NbrLanczosIteration;


 public:

  // constructor from datas
  //
  // blockHamiltonian = Hamiltonian associated to left and right blocks
  // interactionBlockHamiltonian = Hamiltonian associated to interaction blocks 
  // leftInteraction = interaction between left block and left interaction block
  // rightInteraction = interaction between right block and right interaction block
  // lanczosAlgorithm = Lanczos algorithm to use
  // hilbertSpaceSize = number of states kept for a block
  PeriodicDMRGAlgorithm (AbstractHamiltonian* blockHamiltonian, 
			 AbstractHamiltonian* interactionBlockHamiltonian,
			 AbstractInteraction* leftInteraction,
			 AbstractInteraction* rightInteraction,
			 AbstractLanczosAlgorithm* lanczosAlgorithm,
			 int hilbertSpaceSize);

  // destructor
  //
  ~PeriodicDMRGAlgorithm ();

  // force constraint on global quantum number
  //
  // quantumNumber = pointer to the global quantum number to use
  void Constraint(AbstractQuantumNumber* quantumNumber);

  // run DMRG algorithm
  //
  // currentBlockIndex = 
  void RunDMRG(int currentBlockIndex = 0);
 
  // get ground state energy of the last iteration
  // 
  // return value = ground state energy  
  double GetGroundStateEnergy() {return this->GroundStateEnergy;};

  // get truncation error from the last iteration
  // 
  // return value = truncation error
  double GetTruncationError() {return this->TruncationError;};

  // get number of Lanczos iterations needed during the last iteration
  //   
  // return value = number of Lanczos iterations
  int GetNbrLanczosIteration() {return this->NbrLanczosIteration;};

 private:

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
  List<Matrix*> TransformOperator(List<Matrix*>& observables, bool spaceFlag, 
				  BlockDiagonalMatrix& transformationMatrix, 
				  int nbrSubspace, int* keptSubspace, 
				  FullTensorProductStructure** spaceFineStructureArray,
				  int initialBlockSize, int addedBlockSize);
    
  // reduce hamiltonian to the basis of its n-first eigenstates
  //
  // hamiltonian = hamiltonian
  // truncatedSize = number of states to keep
  // return value = reduced hamiltonian (with associated Hilbert space)
  ExplicitHamiltonian* ReduceHamiltonian (AbstractHamiltonian* hamiltonian, int truncatedSize);

};

#endif
