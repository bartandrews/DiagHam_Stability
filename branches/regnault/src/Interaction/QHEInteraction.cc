////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quantum hall effect interaction                 //
//                                                                            //
//                        last modification : 25/11/2002                      //
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


#include "Interaction/QHEInteraction.h"
#include "Tensor/TwoSpaceTensor.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SpaceDecomposition.h"
#include "Matrix/RealSymmetricMatrix.h"


using std::cout;
using std::endl;


// constructor
//
// fermion = flag to indicate if particles are fermions
// maxMomentum = maximum momentum that can be reached in the system
// ratio = ratio between the width in the x direction and the width in the y direction
// precision = precision on interaction coefficient 

QHEInteraction::QHEInteraction(double* couplingConstants, int nbrCouplingConstant)
{
  this->Structure = 0;
  this->NbrCouplingConstant = nbrCouplingConstant;
  this->CouplingConstants = couplingConstants;
  this->LeftSpaceIndex = 0;
  this->RightSpaceIndex = 0;
}

// copy constructor
//
// interaction = reference to interaction to copy

QHEInteraction::QHEInteraction (const QHEInteraction& interaction) 
{
}

// destructor
//

QHEInteraction::~QHEInteraction() 
{
}

// assignment
//
// interaction = reference to interaction to assign
// return value = reference to current interaction

QHEInteraction& QHEInteraction::operator = (const QHEInteraction& interaction) 
{
  return *this;
}

// evaluate interaction between two systems from operators of each system
//
// leftSystemOperators = list of operators associated to the left system
// rightSystemOperators = list of operators associated to the right system
// return value = tensor corresponding to the interaction

TwoSpaceTensor QHEInteraction::Interaction (List<Matrix*>& leftSystemOperators, 
					    List<Matrix*>& rightSystemOperators)
{
  if ((this->NbrCouplingConstant == 0) || (leftSystemOperators.GetNbrElement() < this->NbrCouplingConstant)
      || (rightSystemOperators.GetNbrElement() < this->NbrCouplingConstant))
    return TwoSpaceTensor();
  ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
  ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
  Matrix** Operator1 = IterLeftSystem();
  Matrix** Operator2 = IterRightSystem();
  int SpaceIndex = this->LeftSpaceIndex;
  if (SpaceIndex > this->RightSpaceIndex)
    SpaceIndex = this->RightSpaceIndex;
  TwoSpaceTensor Interaction(this->Structure, *Operator1, *Operator2, this->LeftSpaceIndex,
			     this->CouplingConstants[0]);
  for (int i = 1; i < this->NbrCouplingConstant; i++)
    {
      Operator1 = IterLeftSystem();
      Operator2 = IterRightSystem();
      if (this->CouplingConstants[i] != 0)
	Interaction.AddTensorProductMatrices(*Operator1, *Operator2, this->CouplingConstants[i]);
    }
  return Interaction;
}

// evaluate interaction between two systems from operators of each system with a given
// space decomposition for each space
//
// firstSystemOperators = list of operators associated to the first system
// secondSystemOperators = list of operators associated to the second system
// firstSpaceDecomposition = space decomposition of the space associated to the first system
// secondSpaceDecomposition = space decomposition of the space associated to the second system
// return value = tensor corresponding to the interaction

TwoSpaceTensor QHEInteraction::Interaction (List<Matrix*>& leftSystemOperators, 
					    List<Matrix*>& rightSystemOperators, 
					    SpaceDecomposition& leftSpaceDecomposition,
					    SpaceDecomposition& rightSpaceDecomposition)
{
  if ((this->NbrCouplingConstant == 0) || (leftSystemOperators.GetNbrElement() < this->NbrCouplingConstant)
      || (rightSystemOperators.GetNbrElement() < this->NbrCouplingConstant))
    {
      cout << "error " << this->NbrCouplingConstant << " " << leftSystemOperators.GetNbrElement() << " " << rightSystemOperators.GetNbrElement() << endl;
      return TwoSpaceTensor();
    }
  ListIterator<Matrix*> IterLeftSystem(leftSystemOperators);
  ListIterator<Matrix*> IterRightSystem(rightSystemOperators);
  Matrix** Operator1 = IterLeftSystem();
  Matrix** Operator2 = IterRightSystem();
  int i = 0;
  while (this->CouplingConstants[i] == 0.0)
    {
      Operator1 = IterLeftSystem();
      Operator2 = IterRightSystem();
      i++;
    }
  if (i == this->NbrCouplingConstant)
    {
      RealSymmetricMatrix TmpMatrix(this->Structure->GetTotalDimension(), true);
      TwoSpaceTensor Interaction(this->Structure, (Matrix&) TmpMatrix, this->LeftSpaceIndex);
      return Interaction;
    }
  int SpaceIndex = this->LeftSpaceIndex;
  if (SpaceIndex > this->RightSpaceIndex)
    SpaceIndex = this->RightSpaceIndex;
//  cout << "Operators : " << endl << **Operator1 << endl << **Operator2 << endl;
  TwoSpaceTensor Interaction(this->Structure, **Operator1, **Operator2, 
			     leftSpaceDecomposition, rightSpaceDecomposition,
			     this->LeftSpaceIndex, this->CouplingConstants[i++]);
//  cout << Interaction << endl;
  for (; i < this->NbrCouplingConstant; i++)
    {
      Operator1 = IterLeftSystem();
      Operator2 = IterRightSystem();
//      cout << "Operators : " << endl << **Operator1 << endl << **Operator2 << endl;
      if (this->CouplingConstants[i] != 0)
	Interaction.AddTensorProductMatrices(**Operator1, **Operator2, 
					     leftSpaceDecomposition, rightSpaceDecomposition,
					     this->CouplingConstants[i]);
//      cout << Interaction << endl;
    }
  return Interaction;
}

