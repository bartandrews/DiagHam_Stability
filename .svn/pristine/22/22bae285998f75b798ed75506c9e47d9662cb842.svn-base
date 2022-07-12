////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 19/06/2002                      //
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


#ifndef CLEBSCHGORDANCOEFFICIENTS_H
#define CLEBSCHGORDANCOEFFICIENTS_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients
{

 protected:

  // first angular momentum (twice the value to avoid half integer value)
  int J1;
  // second angular momentum (twice the value to avoid half integer value)
  int J2;

  // Clebsch Gordan coefficient array accessed as Coefficients[j1 + m1][j2 + m2][j3]
  double*** Coefficients;
  // minimum j value in each sector of given m1 and m2
  int** JMin;
  // garbage flag associated to coefficient array
  GarbageFlag Flag;


  // position associated to the projection of first angular momentum for current iterator
  int M1;
  // position associated to the projection of second angular momentum for current iterator
  int M2;
  // resulting angular momentum for current iterator  
  int J;
  // position associated to the resulting angular momentum  for current iterator  
  int CurrentPosition;

 public:

  // default constructor
  //
  ClebschGordanCoefficients();

  // constructor 
  //
  // j1 = first angular momentum (twice the value to avoid half integer value)
  // j2 = second angular momentum (twice the value to avoid half integer value)
  ClebschGordanCoefficients(int j1, int j2);

  // copy constructor (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to copy
  ClebschGordanCoefficients (const ClebschGordanCoefficients& coefficients);

  // destructor
  //
  ~ClebschGordanCoefficients ();

  // assignment (without duplicating datas)
  //
  // coefficients = reference on Clebsch Gordan coefficients to assign
  // return value = reference on current Clebsch Gordan coefficients
  ClebschGordanCoefficients& operator = (const ClebschGordanCoefficients& coefficients);

  // get a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum 
  // j = resulting angular momentum
  // return value = corresponding Clebsch Gordan coefficient
  double GetCoefficient (int m1, int m2, int j);

  // initial an iterator on all Clebsch Gordan coefficients for fixed m1 and m2 values
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum   
  void InitializeCoefficientIterator(int m1, int m2);

  // return next coefficient associated with current iterator (with increasing j value)
  //
  // j = reference on integer where resulting angular momentum has to be stored
  // coefficient = reference on double where Clebsch Gordan coefficient has to be stored
  // return value = false if no coefficient has been returned
  bool Iterate(int& j, double& coefficient);

  // return the values of the angular momenta that have been coupled
  int GetJ1(){return this->J1;}
  int GetJ2(){return this->J2;}

  // print a particular coefficient (without testing if m1, m2 and j are valid)
  //
  // str = reference on output stream
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum 
  // j = resulting angular momentum
  // return value = reference on output stream
  ostream& PrintCoefficient (ostream& str, int m1, int m2, int j);

 protected:
  
  // evaluate all Clebsch Gordan coefficients using Schulten Gordon recursion algorithm
  //
  void EvaluateClebschGordanCoefficients();

  // evaluate C coefficient needed during the recursion
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum   
  // j = resulting angular momentum
  // return value = corresponding D coefficient
  double CCoefficient (int m1, int m2, int j);
  
  // evaluate D coefficient needed during the recursion
  //
  // m1 = projection of first angular momentum 
  // m2 = projection of second angular momentum   
  // j = resulting angular momentum
  // return value = corresponding D coefficient
  double DCoefficient (int m1, int m2, int j);

  // main part of the recursion algorithm
  //
  // m1 = projection of first angular momentum 
  // m2Min = minimum value for the projection of second angular momentum   
  // m2Max = maximum value for the projection of second angular momentum   
  // j = resulting angular momentum
  void RecursionAlgorithm (int m1, int m2Min, int m2Max, int j);

};

#endif


