////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//             class for calculating a determinant of a complex matrix                    //
//                                                                            //
//                        last modification : 06/07/2007                      //
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



#ifndef COMPLEXLAPACKDETERMINANT_H
#define COMPLEXLAPACKDETERMINANT_H

#include "config.h"

#ifdef HAVE_LAPACK
#include "MathTools/Complex.h"
#include <iostream>

using std::ostream;


class ComplexLapackDeterminant
{
 private:
  int Dimension;
  int *Permutation;
  doublecomplex *Components;
  
 public:
  // default constructor
  //
  ComplexLapackDeterminant();

  // constructor for an empty matrix
  //
  // dimension = number of columns and rows
  // zero = tue if matrix elements have to be set to zero
  ComplexLapackDeterminant(int dimension, bool zero = false);

  // destructor
  ~ComplexLapackDeterminant();

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // (re,im) = new value for matrix element
  void SetMatrixElement(int i, int j, const double re, const double im);

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;


  // calculate the determinant, loosing information about the matrix elements
  Complex Determinant();

  // resize determinant dimensions
  void Resize(int newDimension);

  
};

// get a matrix element
//
// i = line position
// j = column position
// x = reference on the variable where to store the requested matrix element

inline void ComplexLapackDeterminant::GetMatrixElement(int i, int j, Complex& x) const
{
  x.Re = this->Components[i+j*Dimension].r;
  x.Im = this->Components[i+j*Dimension].i;
}


#else

#include "ComplexMatrix.h"
#include <iostream>

using std::ostream;


class ComplexLapackDeterminant
{
 private:
  ComplexMatrix *M;
  
 public:
  // default constructor
  //
  ComplexLapackDeterminant();

  // constructor for an empty matrix
  //
  // dimension = number of columns and rows
  // zero = tue if matrix elements have to be set to zero
  ComplexLapackDeterminant(int dimension, bool zero = false);

  // destructor
  ~ComplexLapackDeterminant();

  // set a matrix element
  //
  // i = line position
  // j = column position
  // x = new value for matrix element
  void SetMatrixElement(int i, int j, const Complex& x);

  // set a matrix element
  //
  // i = line position
  // j = column position
  // (re,im) = new value for matrix element
  void SetMatrixElement(int i, int j, const double re, const double im);

  // get a matrix element
  //
  // i = line position
  // j = column position
  // x = reference on the variable where to store the requested matrix element
  void GetMatrixElement(int i, int j, Complex& x) const;


  // calculate the determinant, loosing information about the matrix elements
  Complex Determinant();

  // resize determinant dimensions
  void Resize(int newDimension);

  
};

#endif // HAVE_LAPACK

#endif // COMPLEXLAPACKDETERMINANT


