////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author : Cecile Repellin                      //
//                                                                            //
//              class of SU(3) irreducible representations                    //
//                                                                            //
//                        last modification : 17/02/2013                      //
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


#ifndef SU3IRREDUCIBLEREPRESENTATIONS_H
#define SU3IRREDUCIBLEREPRESENTATIONS_H


#include "config.h"

#include <iostream>


using std::ostream;


class SU3IrreducibleRepresentations
{

  protected:
    //array giving the Q(M) indices for each y, tz, degeneracy
    int*** QIndices;
    //array giving tz, y for each Q(M)
    int** QuantumNumbers;
    //(p,q) = indices that define the SU(3) representation
    int P;
    int Q;
    //representation dimensison
    int RepresentationDimension;
    //array that gives the correspondance between an index Q(M) and the corresponding state (M11, M12, M22)
    int** States;
    //array that gives the degeneracy of a state defined by its quantum numbers (tz, y)
    int** Degeneracy;
    //array the number of accessible values of tz for each value of y
    int* NbrTzValuesPerShiftedY;
    
    
  public:
    
    // constructor 
    //
    // (p,q) = indices that define an irreducible SU(3) representation
    SU3IrreducibleRepresentations(int p, int q);

    
    // destructor
    //
    ~SU3IrreducibleRepresentations ();

    //generate all states in the representation (p,q)
    //
    void IrreducibleRepresentationGenerateStates();
    
    
    //compute the value of tz and y for a given index Q(M)
    //
    //index = value of Q(M)
    //tz = reference on the value of tz
    //y = reference on the value of y
    void GetQuantumNumbers(int index, int& tz, int& y);
    
    //get the degeneracy of a state with quantum numbers (tz, y) in representation (p,q)
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //return value = degeneracy
    int GetDegeneracy(int tz, int y);
    
    //get the Q(M) index corresponding to tz, y, and a given degeneracy
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //degeneracy = degeneracy index
    //return value = Q(M) index
    int GetQIndices(int tz, int y, int degeneracy);
    
    //return the dimension of the representation
    //
    //return value = dimension
    int GetDimension();
    
    
  protected:
    
    // Compute the degeneracy for all values of (tz,y) accessible in representation (p,q)
    //
    void ComputeDegeneracy();
    
    //Computes the value of momenta tz and y given their positive integer value
    //
    //tz = reference on the value of tz
    //y = reference on the value of y
    //shiftedTz = shifted value of tz (positive integer that is used as index in an array)
    //shiftedY = shifted value of y (positive integer that is used as index in an array)
    void GetQuantumNumbersFromShiftedIndices(int& tz, int& y, int shiftedTz, int shiftedY);
    
    //Computes the shifted values of momenta tz and y
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //shiftedTz = reference on the shifted value of tz (positive integer that is used as index in an array)
    //shiftedY = reference to the shifted value of y (positive integer that is used as index in an array)
    void GetShiftedIndices(int tz, int y, int& shiftedTz, int& shiftedY);
    

    
    
};

    //Computes the value of momenta tz and y given their positive integer value
    //
    //tz = reference on the value of tz
    //y = reference on the value of y
    //shiftedTz = shifted value of tz (positive integer that is used as index in an array)
    //shiftedY = shifted value of y (positive integer that is used as index in an array)
 inline void SU3IrreducibleRepresentations::GetQuantumNumbersFromShiftedIndices(int& tz, int& y, int shiftedTz, int shiftedY)
    {
      y = 3*shiftedY - 2*this->P - this->Q;
      int tzMax;
      if (shiftedY < this->P)
	tzMax = shiftedY + this->Q;
      else
	tzMax = 2*this->P + this->Q - shiftedY;
      tz = 2*shiftedTz - tzMax;
    }
    
    
    //Computes the shifted values of momenta tz and y
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //shiftedTz = reference on the shifted value of tz (positive integer that is used as index in an array)
    //shiftedY = reference to the shifted value of y (positive integer that is used as index in an array)
 inline void SU3IrreducibleRepresentations::GetShiftedIndices(int tz, int y, int& shiftedTz, int& shiftedY)
    {
      int tzMax;
      shiftedY = (y + 2*this->P + this->Q)/3;
      if (shiftedY < this->P)
	tzMax = shiftedY + this->Q;
      else
	tzMax = 2*this->P + this->Q - shiftedY;
      shiftedTz = (tzMax + tz)/2;
    }
    
  //compute the value of tz and y for a given index Q(M)
  //
  //index = value of Q(M)
  //tz = reference on the value of tz
  //y = reference on the value of y
  inline void SU3IrreducibleRepresentations::GetQuantumNumbers(int index, int& tz, int& y)
    {
      int lambda1 = 2*this->States[index][2] - (this->States[index][0] + this->States[index][1]);
      int lambda2 = 2*(this->States[index][0] + this->States[index][1]) - (this->States[index][2] + this->P + 2*this->Q);
      tz = lambda1 + lambda2;
      y = lambda1 - lambda2;
    }
    
    //get the degeneracy of a state with quantum numbers (tz, y) in representation (p,q)
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //return value = degeneracy
  inline int SU3IrreducibleRepresentations::GetDegeneracy(int tz, int y)
    {
      int shiftedTz;
      int shiftedY;
      this->GetShiftedIndices(tz, y, shiftedTz, shiftedY);
      return this->Degeneracy[shiftedY][shiftedTz];
    }
    
    //get the Q(M) index corresponding to tz, y, and a given degeneracy
    //
    //tz = value of momentum quantum number tz
    //y = value of momentum quantum number y
    //degeneracy = degeneracy index
    //return value = Q(M) index
  inline int SU3IrreducibleRepresentations::GetQIndices(int tz, int y, int degeneracy)
    {
      int shiftedTz;
      int shiftedY;
      this->GetShiftedIndices(tz, y, shiftedTz, shiftedY);
      return this->QIndices[shiftedY][shiftedTz][degeneracy];
    }
    
    //return the dimension of the representation
    //
    //return value = dimension
  inline int SU3IrreducibleRepresentations::GetDimension()
    {
      return this->RepresentationDimension; 
    }
 
 #endif