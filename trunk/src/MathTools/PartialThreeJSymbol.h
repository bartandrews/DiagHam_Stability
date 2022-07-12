////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class to generate Three-J Symbols                      //
//                                                                            //
//                        last modification : 06/05/2009                      //
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


#include "config.h"
#include "MathTools/ThreeJSymbol.h"

#include <cstdlib>
#include <cmath>

#ifndef PARTIALTHREEJSYMBOL_H
#define PARTIALTHREEJSYMBOL_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;
using std::ofstream;
using std::ifstream;


class PartialThreeJSymbol
{
 protected:

  // first angular momentum (twice the value to avoid half integer value)
  int J1;
  // second angular momentum (twice the value to avoid half integer value)
  int J2;
  // z-component of resulting angular momentum
  int M3;
  // minimum permitted value of j3
  int J3Min;
  // number of valid angular momenta J3
  int NbrJ3;

  // smallest permitted value of M1
  int MinM1;

  // largest permitted value of M1
  int MaxM1;

  // number of M1 values
  int NbrM1;

  // Clebsch Gordan coefficient array accessed as Coefficients[j1 + m1][j3]
  double** Coefficients;
  
  // garbage flag associated to coefficient array
  GarbageFlag Flag;

  // position associated to the projection of first angular momentum for current iterator
  int M1;
  // resulting angular momentum for current iterator  
  int J;
  // position associated to the resulting angular momentum  for current iterator  
  int CurrentPosition;

  
 public:
  // default constructor
  //
  PartialThreeJSymbol();

  // constructor
  //
  // j1 = first angular momentum (twice the value to avoid half integer value)
  // j2 = second angular momentum (twice the value to avoid half integer value)
  // m3 = twice the value of the z-component of the target angular momentum
  PartialThreeJSymbol(int j1, int j2, int m3, ThreeJSymbol *fullSymbol=NULL);
  
  // constructor
  //
  // j1 = first angular momentum (twice the value to avoid half integer value)
  // j2 = second angular momentum (twice the value to avoid half integer value)
  // m3 = twice the value of the z-component of the target angular momentum
  // minM1 = smallest allowed value for m1
  // maxM1 = largest allowed value for m1
  PartialThreeJSymbol(int j1, int j2, int m3, int minM1, int maxM1, ThreeJSymbol *fullSymbol=NULL);

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Three J Symbol to copy
  PartialThreeJSymbol (const PartialThreeJSymbol& symbol);

  // destructor
  //
  ~PartialThreeJSymbol ();

  // assignment (without duplicating datas)
  //
  // coefficients = reference on Three J Symbol to assign
  // return value = reference on current Three J Symbol
  PartialThreeJSymbol& operator = (const PartialThreeJSymbol& symbol);

  // get a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // m1 = projection of first angular momentum 
  // j = resulting angular momentum
  // return value = corresponding Three J symbol
  double GetCoefficient (int m1, int j);

  // initial an iterator on all Clebsch Gordan coefficients for fixed m1 and m2 values
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum   
  void InitializeCoefficientIterator(int m1);

  // return next coefficient associated with current iterator (with increasing j value)
  //
  // j = reference on integer where resulting angular momentum has to be stored
  // coefficient = reference on double where Clebsch Gordan coefficient has to be stored
  // return value = false if no coefficient has been returned
  bool Iterate(int& j, double& coefficient);

  // print a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // str = reference on output stream
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum 
  // j = resulting angular momentum
  // return value = reference on output stream
  ostream& PrintCoefficient (ostream& str, int m1, int j);


  // operator overloads to read from and write to a stream
  friend ofstream& operator << (ofstream& str, const PartialThreeJSymbol& threeJ);
  friend ifstream& operator >> (ifstream& str, PartialThreeJSymbol& threeJ);

  // pointer-friendly I/O functions
  void WriteSymbol(ofstream& str);
  void ReadSymbol(ifstream& str);

 protected:
  // initialize tables, copying values from provided 3J symbol
  // fullSymbol = 3J symbol at j1, j2
  //
  void InitializeTable(ThreeJSymbol *fullSymbol);
  
};

#endif
