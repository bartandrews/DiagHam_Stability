////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                        Class author : Cecile Repellin                      //
//                                                                            //
//                     class of SU(3) clebsch gordan coefficients             //
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


#ifndef SU3CLEBSCHGORDANCOEFFICIENTS_H
#define SU3CLEBSCHGORDANCOEFFICIENTS_H


#include "config.h"

#include <iostream>


using std::ostream;


class SU3ClebschGordanCoefficients 
{
      
  protected:
    //(P1, Q1) representation of the first space
    int P1;
    int Q1;
    //dimension of the (P1,Q1) representation
    int NbrQM1Values;
    //(P2, Q2) representation of the second space
    int P2;
    int Q2;
    //dimension of the (P2,Q2) representation
    int NbrQM2Values;
    //array of all the (p,q) destination spaces
    int** PQ;
    // array of Clebsch Gordan Coefficients
    double**** Coefficients;
    //total number of representations that arise from (p1,q1)*(p2,q2)
    int NbrPQRepresentations;
    //array that gives the degeneracy in every representation (p,q) for every point in the original representation
    int*** Degeneracy;
    
    
    
  public:
    
    // constructor 
    //
    // (p1,q1) = first irrep
    // (p2,q2) = second irrep to be combined with (p1, q1)
    SU3ClebschGordanCoefficients(int p1, int q1, int p2, int q2);

    // copy constructor (without duplicating datas)
    //
    // coefficients = reference on Clebsch Gordan coefficients to copy
    SU3ClebschGordanCoefficients (const SU3ClebschGordanCoefficients& coefficients);
    
    // destructor
    ~SU3ClebschGordanCoefficients ();
    
   
    //get all (p,q) destination representations
    //
    //return value = number of such representations
    int GetAllPQRepresentations();

    //get the degeneracy of a point in the (p,q) representation
    //
    //(p,q) destination representation
    //index = position in the original subspaces index = Q(M1) + dim(p1,q1)*Q(M2)
    //return value = degeneracy
    int GetClebschGordanDegeneracy(int p, int q, int index);

    //get one Clebsch Gordan Coefficient
    //
    //representationIndex = index for the (p,q) destination representation
    //q1Index = index of the first particle in (p1, q1)
    //q2Index = index of the second particle in (p2, q2)
    //qIndex = Q(M) index of the point in the destination representation representation
    //return value = Clebsch Gordan coefficient    
    double GetClebschGordanCoefficient(int representationIndex, int q1Index, int q2Index, int qIndex);
     
    //get the dimension of a (p,q) representation represented by its index index
    //
    //return value  = dimension
    int GetRepresentationDimension(int index);

    // return the degeneracy of a point in representation (p,q)
    //
    // representationIndex = index that defines the (p,q) representation
    // q1Index = index of particle 1 in representation (p1, q1)
    // q2Index = index of particle 2 in representation (p2, q2)
    // return value = degeneracy
    int GetDegeneracy(int representationIndex, int q1Index, int q2Index);
     
    
    //return the number of different (p,q) representations
    //
    //return value = number of representations
    int GetNbrPQRepresentations();
     
  protected:
  
    //evaluate all Clebsch Gordan coefficients
    //
    void EvaluateClebschGordanCoefficients();    
    
};

    //get one Clebsch Gordan Coefficient
    //
    //representationIndex = index for the (p,q) destination representation
    //q1Index = index of the first particle in (p1, q1)
    //q2Index = index of the second particle in (p2, q2)
    //qIndex = Q(M) index of the point in the destination representation representation
    //return value = Clebsch Gordan coefficient    
    inline double SU3ClebschGordanCoefficients::GetClebschGordanCoefficient(int representationIndex, int q1Index, int q2Index, int qIndex)
      {
	return this->Coefficients[representationIndex][q1Index][q2Index][qIndex]; 
      }
      
      //get the dimension of a (p,q) representation represented by its index index
      //
      //return value  = dimension
      inline int SU3ClebschGordanCoefficients::GetRepresentationDimension(int index)
      {
	return (this->PQ[index][0] + 1) * (this->PQ[index][1] + 1) * (this->PQ[index][0] + this->PQ[index][1] + 2)/2;
      }
      
      // return the degeneracy of a point in representation (p,q)
      //
      // representationIndex = index that defines the (p,q) representation
      // q1Index = index of particle 1 in representation (p1, q1)
      // q2Index = index of particle 2 in representation (p2, q2)
      // return value = degeneracy
      inline int SU3ClebschGordanCoefficients::GetDegeneracy(int representationIndex, int q1Index, int q2Index)
	{
	  return this->Degeneracy[representationIndex][q1Index][q2Index]; 
	}
	
	//return the number of different (p,q) representations
	//
	//return value = number of representations
      inline int SU3ClebschGordanCoefficients::GetNbrPQRepresentations()
	{
	  return this->NbrPQRepresentations; 
	}
	
	

#endif