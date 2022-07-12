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


#include "config.h"
#include "SU3ClebschGordanCoefficients.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// constructor 
//
// (p1,q1) = first irrep
// (p2,q2) = second irrep to be combined with (p1, q1)

SU3ClebschGordanCoefficients::SU3ClebschGordanCoefficients(int p1, int q1, int p2, int q2)
{
  this->P1 = p1;
  this->P2 = p2;
  this->Q1 = q1;
  this->Q2 = q2;
  this->NbrQM1Values = (this->P1 + 1)*(this->Q1 + 1)*(this->P1 + this->Q1 + 2)/2;
  this->NbrQM2Values = (this->P2 + 1)*(this->Q2 + 1)*(this->P2 + this->Q2 + 2)/2;
  this->NbrPQRepresentations = this->GetAllPQRepresentations();
  this->Degeneracy = new int**[this->NbrPQRepresentations];
  this->Coefficients = new double***[this->NbrPQRepresentations];
  for (int i = 0; i < this->NbrPQRepresentations; ++i)
  {
    this->Degeneracy[i] = new int*[this->NbrQM1Values];
    this->Coefficients[i] = new double**[this->NbrQM1Values];
    for (int j = 0; j < this->NbrQM1Values; ++j)
    {
      this->Degeneracy[i][j] = new int[this->NbrQM2Values];
      this->Coefficients[i][j] = new double*[this->NbrQM2Values];
      for (int k = 0; k < this->NbrQM2Values; ++k)
      {
	this->Degeneracy[i][j][k] = this->GetClebschGordanDegeneracy(this->PQ[i][0], this->PQ[i][1], j);
	this->Coefficients[i][j][k] = new double[this->GetRepresentationDimension(i)];
	for (int l = 0; l < this->GetRepresentationDimension(i); ++l)
	  this->Coefficients[i][j][k][l] = 0;
      }
    }
  }
  this->EvaluateClebschGordanCoefficients();
}

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy

SU3ClebschGordanCoefficients::SU3ClebschGordanCoefficients (const SU3ClebschGordanCoefficients& coefficients)
{
  this->P1 = coefficients.P1;
  this->Q1 = coefficients.Q1;
  this->P2 = coefficients.P2;
  this->Q2 = coefficients.Q2;
  this->Coefficients = coefficients.Coefficients;
  this->PQ = coefficients.PQ;
  this->NbrPQRepresentations = coefficients.NbrPQRepresentations;
  this->Degeneracy = coefficients.Degeneracy;
//   this->QIndices = coefficients.QIndices;
}

//destructor
SU3ClebschGordanCoefficients::~SU3ClebschGordanCoefficients()
{
 delete[] this->PQ;
 delete[] this->Degeneracy;
 delete[] this->Coefficients;
//  delete[] this->QIndices;
}

void SU3ClebschGordanCoefficients::EvaluateClebschGordanCoefficients()
{
 char* FileName = 0;
 int TmpNbrLinesInFile;
 int* TmpMArray;
 int* TmpM1Array;
 int* TmpM2Array;
 double* TmpClebschArray;
 int p;
 int q;
 for (int i = 0; i < this->NbrPQRepresentations; ++i)
 {
   
   p = this->PQ[i][0];
   q = this->PQ[i][1];
   FileName = new char[64];
   sprintf(FileName, "CGCSU3_p1_%d_q1_%d_p2_%d_q2_%d_p_%d_q_%d.dat", this->P1, this->Q1, this->P2, this->Q2, p, q);
   MultiColumnASCIIFile ClebschFile;
   if (ClebschFile.Parse(FileName) == false)
    {
      ClebschFile.DumpErrors(cout);
    }
   TmpNbrLinesInFile = ClebschFile.GetNbrLines();
   if (ClebschFile.GetNbrColumns() != 4)
   {
    cout << "CGC file should have 4 columns" << endl;
   }
   TmpClebschArray = ClebschFile.GetAsDoubleArray(3); 
   TmpMArray = ClebschFile.GetAsIntegerArray(0); 
   TmpM1Array = ClebschFile.GetAsIntegerArray(1); 
   TmpM2Array = ClebschFile.GetAsIntegerArray(2);
   for (int j = 0; j < this->NbrQM1Values; ++j)
   {
     for (int k = 0; k < this->NbrQM2Values; ++k)
     {
      this->Degeneracy[i][j][k] = 0;
      for (int l = 0; l < TmpNbrLinesInFile; ++l)
      {
	if ((TmpM1Array[l] - 1 == j) && (TmpM2Array[l] - 1 == k))
	{
	  this->Coefficients[i][j][k][TmpMArray[l] - 1] = TmpClebschArray[l];
	  this->Degeneracy[i][j][k] += 1;
       }
      }
     }
   }  
 }
 delete[] FileName;
 delete[] TmpM1Array;
 delete[] TmpM2Array;
 delete[] TmpClebschArray;
}

int SU3ClebschGordanCoefficients::GetAllPQRepresentations()
{
  if ((this->P1 == this->P2) && (this->Q1 == 0) && (this->Q2 == 0))
    {
      this->PQ = new int*[this->P1 + 1];
      for (int i = 0; i <= this->P1; ++i)
	{
	  this->PQ[i] = new int[2];
	  this->PQ[i][0] = 2*(this->P1 - i);
	  this->PQ[i][1] = i;
	}
      return (this->P1 + 1);
    }
  else
    {
      return 0;
    }
}

int SU3ClebschGordanCoefficients::GetClebschGordanDegeneracy(int p, int q, int index)
{
  char* FileName = new char[64];
  int* TmpM1Array;
  int* TmpM2Array;
  sprintf(FileName, "CGCSU3_p1_%d_q1_%d_p2_%d_q2_%d_p_%d_q_%d.dat", this->P1, this->Q1, this->P2, this->Q2, p, q);
  MultiColumnASCIIFile ClebschFile;
  if (ClebschFile.Parse(FileName) == false)
    {
      ClebschFile.DumpErrors(cout);
    }
  if (ClebschFile.GetNbrColumns() != 4)
   {
    cout << "CGC file should have 4 columns" << endl;
   }
  TmpM1Array = ClebschFile.GetAsIntegerArray(1); 
  TmpM2Array = ClebschFile.GetAsIntegerArray(2);
  int NbrLines = ClebschFile.GetNbrLines();
  int Q1 = index % this->NbrQM1Values;
  int Q2 = index / this->NbrQM1Values;
  int degeneracy = 0;
  for (int i = 0; i < NbrLines; ++i)
  {
   if ((Q1 == TmpM1Array[i]) && (Q2 == TmpM2Array[i]))
     degeneracy += 1;
  }
  delete[] FileName;
  delete[] TmpM1Array;
  delete[] TmpM2Array;
  return degeneracy;
}
