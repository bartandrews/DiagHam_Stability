////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of diamond interaction                       //
//                                                                            //
//                        last modification : 14/12/2001                      //
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


#include "Interaction/DiamondInteraction.h"
#include "Tensor/TwoSpaceTensor.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SpaceDecomposition.h"

#include <iostream>


using std::dec;
using std::cout;
using std::endl;


// default constructor
//

DiamondInteraction::DiamondInteraction() 
{
  this->UpperCouplingConstant = 1.0;
  this->LowerCouplingConstant = 1.0;
  this->LeftSpaceIndex = 0;
  this->RightSpaceIndex = 0;
  this->Structure = 0;
}

// constructor from partial datas
//
// upperCouplingConstant = couling constant for upper link
// lowerCoupligConstant = couling constant for lower link

DiamondInteraction::DiamondInteraction(double upperCouplingConstant, double lowerCouplingConstant)
{
  this->Structure = 0;
  this->UpperCouplingConstant = upperCouplingConstant;
  this->LowerCouplingConstant = lowerCouplingConstant;
  this->LeftSpaceIndex = 0;
  this->RightSpaceIndex = 0;
}

// constructor from complete datas
//
// upperCouplingConstant = couling constant for upper link
// lowerCoupligConstant = couling constant for lower link
// leftSpaceIndex = 
// rightSpaceIndex = 
// struture = reference on tensor product structure

DiamondInteraction::DiamondInteraction(double upperCouplingConstant, double lowerCouplingConstant, 
				       int leftSpaceIndex, int rightSpaceIndex, 
				       AbstractTensorProductStructure* structure)
{
  this->UpperCouplingConstant = upperCouplingConstant;
  this->LowerCouplingConstant = lowerCouplingConstant;
  this->LeftSpaceIndex = leftSpaceIndex;
  this->RightSpaceIndex = rightSpaceIndex;
  this->Structure = structure;
}

// copy constructor
//
// interaction = reference to interaction to copy

DiamondInteraction::DiamondInteraction (const DiamondInteraction& interaction) 
{
  this->UpperCouplingConstant = interaction.UpperCouplingConstant;
  this->LowerCouplingConstant = interaction.LowerCouplingConstant;
  this->LeftSpaceIndex = interaction.LeftSpaceIndex;
  this->RightSpaceIndex = interaction.RightSpaceIndex;
  this->Structure = interaction.Structure;
}

// destructor
//

DiamondInteraction::~DiamondInteraction() 
{
}

// assignment
//
// interaction = reference to interaction to assign
// return value = reference to current interaction

DiamondInteraction& DiamondInteraction::operator = (const DiamondInteraction& interaction) 
{
  this->UpperCouplingConstant = this->UpperCouplingConstant;
  this->LowerCouplingConstant = this->LowerCouplingConstant;
  this->LeftSpaceIndex = interaction.LeftSpaceIndex;
  this->RightSpaceIndex = interaction.RightSpaceIndex;
  this->Structure = interaction.Structure;
  return *this;
}

// evaluate interaction between two systems from operators of each system
//
// leftSystemOperators = list of operators associated to the left system
// rightSystemOperators = list of operators associated to the right system
// return value = tensor corresponding to the interaction

TwoSpaceTensor DiamondInteraction::Interaction (List<Matrix*>& leftSystemOperators, 
						List<Matrix*>& rightSystemOperators)
{
  if ((leftSystemOperators.GetNbrElement() == 3) && (rightSystemOperators.GetNbrElement() == 6))
    {
      ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
      ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
      Matrix** OperatorSx = IterLeftSystem();
      Matrix** OperatorSy = IterLeftSystem();
      Matrix** OperatorSz = IterLeftSystem();
      Matrix** Operator2 = IterRightSystem();
      int SpaceIndex = this->LeftSpaceIndex;
      if (SpaceIndex > this->RightSpaceIndex)
	SpaceIndex = this->RightSpaceIndex;
      TwoSpaceTensor Interaction(this->Structure, *OperatorSx, *Operator2, this->LeftSpaceIndex,
				 this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(*OperatorSy, *Operator2, -this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(*OperatorSz, *Operator2, this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(*OperatorSx, *Operator2, this->LowerCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(*OperatorSy, *Operator2, -this->LowerCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(*OperatorSz, *Operator2, this->LowerCouplingConstant);
      return Interaction;
    }
  if ((leftSystemOperators.GetNbrElement() == 6) && (rightSystemOperators.GetNbrElement() == 3))
    {
      ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
      ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
      Matrix** OperatorSx = IterRightSystem();
      Matrix** OperatorSy = IterRightSystem();
      Matrix** OperatorSz = IterRightSystem();
      Matrix** Operator2 = IterLeftSystem();
      int SpaceIndex = this->LeftSpaceIndex;
      if (SpaceIndex > this->RightSpaceIndex)
	SpaceIndex = this->RightSpaceIndex;
      TwoSpaceTensor Interaction(this->Structure, *Operator2, *OperatorSx, this->LeftSpaceIndex,
				 this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(*Operator2, *OperatorSy, -this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(*Operator2, *OperatorSz, this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(*Operator2, *OperatorSx, this->LowerCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(*Operator2, *OperatorSy, -this->LowerCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(*Operator2, *OperatorSz, this->LowerCouplingConstant);
      return Interaction;
    }
  cout << "error!!!!!! " << leftSystemOperators.GetNbrElement() << " " << rightSystemOperators.GetNbrElement() << endl;
  return TwoSpaceTensor();
}

// evaluate interaction between two systems from operators of each system with a given
// space decomposition for each space
//
// firstSystemOperators = list of operators associated to the first system
// secondSystemOperators = list of operators associated to the second system
// firstSpaceDecomposition = space decomposition of the space associated to the first system
// secondSpaceDecomposition = space decomposition of the space associated to the second system
// return value = tensor corresponding to the interaction

TwoSpaceTensor DiamondInteraction::Interaction (List<Matrix*>& leftSystemOperators, 
						List<Matrix*>& rightSystemOperators, 
						SpaceDecomposition& leftSpaceDecomposition,
						SpaceDecomposition& rightSpaceDecomposition)
{
  if ((leftSystemOperators.GetNbrElement() == 3) && (rightSystemOperators.GetNbrElement() == 6))
    {
      ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
      ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
      Matrix** OperatorSx = IterLeftSystem();
      Matrix** OperatorSy = IterLeftSystem();
      Matrix** OperatorSz = IterLeftSystem();
      Matrix** Operator2 = IterRightSystem();
      int SpaceIndex = this->LeftSpaceIndex;
      if (SpaceIndex > this->RightSpaceIndex)
	SpaceIndex = this->RightSpaceIndex;
      TwoSpaceTensor Interaction(this->Structure, **OperatorSx, **Operator2, 
				 leftSpaceDecomposition, rightSpaceDecomposition,
				 this->LeftSpaceIndex, this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(**OperatorSy, **Operator2, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   -this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(**OperatorSz, **Operator2, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->UpperCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(**OperatorSx, **Operator2, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->LowerCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(**OperatorSy, **Operator2, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   -this->LowerCouplingConstant);
      Operator2 = IterRightSystem();
      Interaction.AddTensorProductMatrices(**OperatorSz, **Operator2, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->LowerCouplingConstant);
      return Interaction;
    }
  if ((leftSystemOperators.GetNbrElement() == 6) && (rightSystemOperators.GetNbrElement() == 3))
    {
      ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
      ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
      Matrix** OperatorSx = IterRightSystem();
      Matrix** OperatorSy = IterRightSystem();
      Matrix** OperatorSz = IterRightSystem();
      Matrix** Operator2 = IterLeftSystem();
      int SpaceIndex = this->LeftSpaceIndex;
      if (SpaceIndex > this->RightSpaceIndex)
	SpaceIndex = this->RightSpaceIndex;
      TwoSpaceTensor Interaction(this->Structure, **Operator2, **OperatorSx, 
				 leftSpaceDecomposition, rightSpaceDecomposition,
				 this->LeftSpaceIndex, this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(**Operator2, **OperatorSy, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   -this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(**Operator2, **OperatorSz, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->UpperCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(**Operator2, **OperatorSx,
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->LowerCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(**Operator2, **OperatorSy, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   -this->LowerCouplingConstant);
      Operator2 = IterLeftSystem();
      Interaction.AddTensorProductMatrices(**Operator2, **OperatorSz, 
					   leftSpaceDecomposition, rightSpaceDecomposition,
					   this->LowerCouplingConstant);
      return Interaction;
    }
  cout << "error!!!!!!" << leftSystemOperators.GetNbrElement() << " " << rightSystemOperators.GetNbrElement() << endl;
  return TwoSpaceTensor();

}

