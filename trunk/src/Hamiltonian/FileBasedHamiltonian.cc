////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of hamiltonian where matrix elements are stored in a text file     //
//                                                                            //
//                        last modification : 31/03/2010                      //
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


#include "Hamiltonian/FileBasedHamiltonian.h"
#include "MathTools/Complex.h" 
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"


#include <iostream>


using std::cout;
using std::endl;


// contructor from default datas without boundary operators
//
// fileName = name of the file where the matrix is stored 
// elementColumnIndex = index of the column where matrix elements are stored
// symmetricFlag = hamiltonian is stored using only the upper or lower triangular part
// fortranIndices = indicates that indices use fortran convention (i.e. 1 based)
// nbrSkippedLines = number of lines to skip in the input file
// nosortRowIndices = if true, assume that the row indices are sorted from the smallest to the largest

FileBasedHamiltonian::FileBasedHamiltonian(char* fileName, int elementColumnIndex, bool symmetricFlag, bool fortranIndices, int nbrSkippedLines, bool nosortRowIndices)
{
  this->NbrElements = GetFileNbrLines(fileName);
  this->SymmetricStorageFlag = symmetricFlag;
  this->HamiltonianShift = 0.0;
  this->NbrElements -= nbrSkippedLines;
  if (this->NbrElements > 0)
    {
      this->RowIndices = new int [this->NbrElements];
      this->ColumnIndices = new int [this->NbrElements];
      this->MatrixElements = new double [this->NbrElements];
      ifstream File;
      File.open(fileName, ios::binary | ios::in);
      if (nbrSkippedLines > 0)
	{
	  
	  char TmpChar;
	  while (nbrSkippedLines > 0)
	    {
	      File.read(&TmpChar, 1);
	      if (TmpChar == '\n')
		--nbrSkippedLines;
	    }
	}
      switch (elementColumnIndex)
	{
	case 0:
	  for (long i = 0; i < this->NbrElements; ++i)
	    {
	      File >> this->MatrixElements[i] >> this->ColumnIndices[i] >> this->RowIndices[i];
	    }
	  break;
	case 1:
	  for (long i = 0; i < this->NbrElements; ++i)
	    File >> this->ColumnIndices[i] >> this->MatrixElements[i] >> this->RowIndices[i];
	  break;
	case 2:
	  for (long i = 0; i < this->NbrElements; ++i)
	    {
	      File >> this->ColumnIndices[i] >> this->RowIndices[i] >> this->MatrixElements[i];
	    }
	  break;
	}
      File.close();

      cout << "nbr read matrix elements = " << this->NbrElements << endl;
      int HamiltonianDimension = 0;
      if (fortranIndices == true)
	{
	  for (long i = 0; i < this->NbrElements; ++i)
	    {
	      if (this->RowIndices[i] > HamiltonianDimension)
		HamiltonianDimension = this->RowIndices[i];
	      if (this->ColumnIndices[i] > HamiltonianDimension)
		HamiltonianDimension = this->ColumnIndices[i];
	      --this->RowIndices[i];
	      --this->ColumnIndices[i];
	    }
	}
      else
	{
	  for (long i = 0; i < this->NbrElements; ++i)
	    {
	      if (this->RowIndices[i] > HamiltonianDimension)
		HamiltonianDimension = this->RowIndices[i];
	      if (this->ColumnIndices[i] > HamiltonianDimension)
		HamiltonianDimension = this->ColumnIndices[i];
	    }
	  ++HamiltonianDimension;
	}
      this->HilbertSpace = new UndescribedHilbertSpace(HamiltonianDimension);
    }
  if (nosortRowIndices == false)
    {
      SortArrayUpOrdering<double>(this->RowIndices, this->ColumnIndices, this->MatrixElements, this->NbrElements);
      // for (int i= 0; i < this->NbrElements; ++i)
      // 	{
      // 	  cout << i << " : " << this->RowIndices[i] << " " << this->ColumnIndices[i] << " : " << this->MatrixElements[i] << endl;
      // 	}
    }
}

// destructor
//

FileBasedHamiltonian::~FileBasedHamiltonian() 
{
  delete[] this->RowIndices;
  delete[] this->ColumnIndices;
  delete[] this->MatrixElements;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void FileBasedHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace) 
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* FileBasedHamiltonian::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int FileBasedHamiltonian::GetHilbertSpaceDimension () 
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void FileBasedHamiltonian::ShiftHamiltonian (double shift) 
{
  this->HamiltonianShift = shift;
}


// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex FileBasedHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  return Complex();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex FileBasedHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& FileBasedHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  long StartingIndex = 0l;
  long LastIndex = this->NbrElements - 1l;
  long MidIndex = 0l;
  int LastComponent = firstComponent + nbrComponent;
  while ((LastIndex - StartingIndex) > 1l)
    {      
      MidIndex = (LastIndex + StartingIndex) >> 1;
      if (this->RowIndices[MidIndex] <= firstComponent)
	StartingIndex = MidIndex;
      else
	LastIndex = MidIndex;
    }
  if ((this->RowIndices[LastIndex] == firstComponent) || (this->RowIndices[StartingIndex] == firstComponent))
    {
      if (this->RowIndices[LastIndex] == firstComponent)
	StartingIndex = LastIndex;
      while ((StartingIndex >= 0l) && (this->RowIndices[StartingIndex] == firstComponent))
	--StartingIndex;
      ++StartingIndex;
      while ((StartingIndex < this->NbrElements) && (this->RowIndices[StartingIndex] < LastComponent))
	{
	  vDestination[this->ColumnIndices[StartingIndex]] += this->MatrixElements[StartingIndex] * vSource[this->RowIndices[StartingIndex]];
	  ++StartingIndex;
	}
    }
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* FileBasedHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  long StartingIndex = 0l;
  long LastIndex = this->NbrElements - 1l;
  long MidIndex = 0l;
  int LastComponent = firstComponent + nbrComponent;
  while ((LastIndex - StartingIndex) > 1l)
    {      
      MidIndex = (LastIndex + StartingIndex) >> 1;
      if (this->RowIndices[MidIndex] <= firstComponent)
	StartingIndex = MidIndex;
      else
	LastIndex = MidIndex;
    }
  if ((this->RowIndices[LastIndex] == firstComponent) || (this->RowIndices[StartingIndex] == firstComponent))
    {
      if (this->RowIndices[LastIndex] == firstComponent)
	StartingIndex = LastIndex;
      while ((StartingIndex >= 0l) && (this->RowIndices[StartingIndex] == firstComponent))
	--StartingIndex;
      ++StartingIndex;
      double TmpMatrixElement;
      int InputIndex;
      int OutputIndex;
      while ((StartingIndex < this->NbrElements) && (this->RowIndices[StartingIndex] < LastComponent))
	{
	  TmpMatrixElement = this->MatrixElements[StartingIndex];
	  InputIndex = this->RowIndices[StartingIndex];
	  OutputIndex = this->ColumnIndices[StartingIndex];
	  for (int k= 0; k < nbrVectors; ++k)
	    vDestinations[k][OutputIndex] += TmpMatrixElement * vSources[k][InputIndex];
	  ++StartingIndex;
	}
    }
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  RealVector& TmpDestination = vDestinations[k];
	  RealVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}
