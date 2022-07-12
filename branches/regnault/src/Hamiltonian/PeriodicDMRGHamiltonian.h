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


#ifndef PERIODICDMRGHAMILTONIAN_H
#define PERIODICDMRGHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Tensor/OneSpaceTensor.h"
#include "TensorProduct/CompositeTensorProductStructure.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class Matrix;
class FullTensorProductStructure;


class PeriodicDMRGHamiltonian : public AbstractHamiltonian
{

 protected:
  
  AbstractHilbertSpace* HilbertSpace;

  OneSpaceTensor LeftBlockLeftInteractionPart;
  OneSpaceTensor RightBlockRightInteractionPart;
  Matrix* Interaction;
  CompositeTensorProductStructure* SpaceStructure;

  int BlockTotalDimension;
  int InteractionBlockTotalDimension;

  int NbrInteractingStateGroupBulk;
  int* GroupPositionBulk;
  int* GroupSizeBulk;
  int* LeftRightAddedBlockGlobalIndexBulk;
  int* LeftRightAddedBlockIndexBulk;
  double** InteractionFactorBulk;
  int* NbrInteractionBulk;
  int** InteractionIndexBulk;

  int NbrInteractingStateGroupBoundary;
  int* GroupPositionBoundary;
  int* GroupSizeBoundary;
  int* LeftRightAddedBlockGlobalIndexBoundary;
  int* LeftRightAddedBlockIndexBoundary;
  double** InteractionFactorBoundary;
  int* NbrInteractionBoundary;
  int** InteractionIndexBoundary;

  int TestInteractingStateFactor;
  bool* TestInteractingStateFlags;
  double InteractionPrecision;

 public:

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
  PeriodicDMRGHamiltonian(OneSpaceTensor leftBlockLeftInteractionPart, 
			  OneSpaceTensor rightBlockRightInteractionPart, 
			  Matrix& interaction,
			  CompositeTensorProductStructure* spaceStructure,
			  FullTensorProductStructure** leftSpaceFineStructureArray, 
			  FullTensorProductStructure** rightSpaceFineStructureArray,
			  int blockTotalDimension, int interactionBlockTotalDimension);

  // destructor
  //
  ~PeriodicDMRGHamiltonian();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);

  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (RealVector& V1, RealVector& V2);
  
  // evaluate matrix element
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // return value = corresponding matrix element
  Complex MatrixElement (ComplexVector& V1, ComplexVector& V2);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of idinces 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
		       int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
			  int firstComponent, int nbrComponent);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // return value = reference on vectorwhere result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
				     int firstComponent, int nbrComponent);
 
  // return a list of left interaction operators
  //
  // return value = list of left interaction operators
  List<Matrix*> LeftInteractionOperators();  

  // return a list of right interaction operators 
  //
  // return value = list of right interaction operators
  List<Matrix*> RightInteractionOperators();  

  // Output Stream overload
  //
  // Str = reference on output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend ostream& operator << (ostream& Str, PeriodicDMRGHamiltonian& H);

  // Mathematica Output Stream overload
  //
  // Str = reference on Mathematica output stream
  // H = Hamiltonian to print
  // return value = reference on output stream
  friend MathematicaOutput& operator << (MathematicaOutput& Str, PeriodicDMRGHamiltonian& H);

 private:

  // find all components of a state that are coupled by interaction  
  //
  // leftSpaceFineStructureArray = array of fine structure describing total left space
  // rightSpaceFineStructureArray = array of fine structure describing total right space
  void InitializeNonNullInteractions(FullTensorProductStructure** leftSpaceFineStructureArray, 
				     FullTensorProductStructure** rightSpaceFineStructureArray);

  // multiply a vector by the current hamiltonian part corresponding to interaction between left and 
  // right blocks and add result to another vector, for a given range of components
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  RealVector& InteractionAddMultiply(RealVector& vSource, RealVector& vDestination,
				     int firstComponent, int NbrComponent);

};

#endif
