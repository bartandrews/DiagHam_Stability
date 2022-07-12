////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class of Fermions                             //
//                                                                            //
//                        last modification : 11/05/2001                      //
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


#include "HilbertSpace/DMRGHilbertSpace/DMRGPartialHilbertSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "GeneralTools/ListIterator.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"


// default constructor
//

DMRGPartialHilbertSpace::DMRGPartialHilbertSpace () 
{
  this->HilbertSpaceDimension = 0;
  this->Energies = RealDiagonalMatrix();
  this->ListQuantumNumber = List<AbstractQuantumNumber*>();
}
  
// constructor from datas 
//
// energies = matrix conatining energies
// quantumNumbers = list of quantum numbers
// decomposition = space decompostion with respect to the subspaces

DMRGPartialHilbertSpace::DMRGPartialHilbertSpace (RealDiagonalMatrix& energies, 
						  List<AbstractQuantumNumber*> quantumNumbers, 
						  const SpaceDecomposition& decomposition)
{
  this->Energies = energies;
  this->HilbertSpaceDimension = this->Energies.GetNbrRow();
  ListIterator<AbstractQuantumNumber*> IterQ(quantumNumbers);
  AbstractQuantumNumber** TmpQ;
  while ((TmpQ = IterQ()))
    this->ListQuantumNumber += (*TmpQ)->Clone();
  this->Decomposition = decomposition;
  this->Flag.Initialize();
}
  
// copy constructor (without duplicating datas)
//
// space = reference on DMRG partial Hilbert space to copy

DMRGPartialHilbertSpace::DMRGPartialHilbertSpace (const DMRGPartialHilbertSpace& space) 
{
  this->Energies = space.Energies;
  this->HilbertSpaceDimension = this->Energies.GetNbrRow();
  this->ListQuantumNumber = space.ListQuantumNumber;
  this->Decomposition = space.Decomposition;
  this->Flag = space.Flag;
}

// destructor
//

DMRGPartialHilbertSpace::~DMRGPartialHilbertSpace () 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      ListIterator<AbstractQuantumNumber*> IterQ(this->ListQuantumNumber);
      AbstractQuantumNumber** TmpQ;
      while ((TmpQ = IterQ()))
	delete *TmpQ;      
    }
}

// assignement (without duplicating datas)
//
// space = reference on DMRG partial Hilbert space to copy
// return value = reference on current fermions

DMRGPartialHilbertSpace& DMRGPartialHilbertSpace::operator = (const DMRGPartialHilbertSpace& space) 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      ListIterator<AbstractQuantumNumber*> IterQ(this->ListQuantumNumber);
      AbstractQuantumNumber** TmpQ;
      while ((TmpQ = IterQ()))
	delete *TmpQ;      
    }
  this->Energies = space.Energies;
  this->HilbertSpaceDimension = this->Energies.GetNbrRow();
  this->ListQuantumNumber = space.ListQuantumNumber;
  this->Decomposition = space.Decomposition;
  this->Flag = space.Flag;
  return *this; 
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DMRGPartialHilbertSpace::Clone()
{
  return new DMRGPartialHilbertSpace(*this);
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int DMRGPartialHilbertSpace::GetHilbertSpaceDimension() 
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> DMRGPartialHilbertSpace::GetQuantumNumbers () 
{
  return this->ListQuantumNumber;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* DMRGPartialHilbertSpace::GetQuantumNumber (int index) 
{
  return this->ListQuantumNumber[this->Decomposition.FindComponent(index)];
}

// get energy of a given state
//
// index = index of the state to test
// return value = energy

double DMRGPartialHilbertSpace::GetEnergy (int index) 
{
  return this->Energies[index];
}

// get description of a given subspace
//
// index = index of the subspace to describe
// return value = reference on subspace description

SubspaceSpaceConverter& DMRGPartialHilbertSpace::GetSubspaceDescription (int index)
{
  return this->Decomposition.GetSubspaceDescription(index);
}
  
// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* DMRGPartialHilbertSpace::ExtractSubspace (AbstractQuantumNumber& q, 
								SubspaceSpaceConverter& converter) 
{
  AbstractQuantumNumber** TmpQ = 0;
  ListIterator<AbstractQuantumNumber*> IterQ(this->ListQuantumNumber);
  int Pos = 0;
  while (((TmpQ = IterQ())) && ((*TmpQ)->IsDifferent(q)))
    Pos++;
  if (TmpQ == 0)
    return new DMRGPartialHilbertSpace();
  int Dim =  this->Decomposition.GetSubspaceDimension(Pos);
  RealDiagonalMatrix TmpM(Dim);
  for (int i = 0; i < Dim; i++)
    TmpM[i] = this->Energies[i + this->Decomposition.GetSubspacePosition(Pos)];
  int* TmpSubspace = new int [1];
  TmpSubspace[0] = 0;
  List<AbstractQuantumNumber*> TmpListQ;
  TmpListQ += *TmpQ;
  return new DMRGPartialHilbertSpace (TmpM, TmpListQ, 
				      SpaceDecomposition(this->Decomposition.GetSubspaceDimension(Pos), 
							 1, TmpSubspace));
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DMRGPartialHilbertSpace::PrintState (ostream& Str, int state) 
{
  if ((state < this->HilbertSpaceDimension) && (state >= 0))
    Str << "State " << state << " Energy = " << this->Energies[state];
  return Str; 
}
