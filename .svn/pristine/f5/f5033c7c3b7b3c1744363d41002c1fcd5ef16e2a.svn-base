////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of herimitian hamiltonian where                    //
//                 matrix elements are stored in a text file                  //
//                                                                            //
//                        last modification : 22/11/2011                      //
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


#ifndef FILEBASEDHERMITIANHAMILTONIAN_H
#define FILEBASEDHERMITIANHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::ostream;
class MathematicaOutput;
class Matrix;


class FileBasedHermitianHamiltonian : public AbstractHamiltonian
{

 protected:

  // Hilbert space assocaited to the Hamiltonian 
  AbstractHilbertSpace* HilbertSpace;
  
  // number of non-zero matrix elements
  long NbrElements;
  // array of row indices  
  int* RowIndices;
  // array of column indices
  int* ColumnIndices;
  // array of matrix elements
  Complex* MatrixElements;
  // true if only the upper triangular part of the Hamiltonian is stored
  bool SymmetricStorageFlag;

  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;


 public:

  // contructor from default datas without boundary operators
  //
  // fileName = name of the file where the matrix is stored 
  // elementColumnIndex = index of the column where matrix elements are stored
  // symmetricFlag = hamiltonian is stored using only the upper or lower triangular part
  // fortranIndices = indicates that indices use fortran convention (i.e. 1 based)
  // nbrSkippedLines = number of lines to skip in the input file
  // nosortRowIndices = if true, assume that the row indices are sorted from the smallest to the largest
  FileBasedHermitianHamiltonian(char* fileName, int elementColumnIndex = 0, bool symmetricFlag = false, bool fortranIndices = false, int nbrSkippedLines = 0, bool nosortRowIndices = false);

  // destructor
  //
  ~FileBasedHermitianHamiltonian();

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

  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
					     int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
						     int firstComponent, int nbrComponent);


};

#endif
