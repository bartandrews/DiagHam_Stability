////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                     class of clebsch gordan coefficients                   //
//                                                                            //
//                        last modification : 02/06/2009                      //
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
#include "MathTools/PartialThreeJSymbol.h"

#include <stdlib.h>
#include <math.h>

using std::cout;
using std::endl;

// default constructor
//

PartialThreeJSymbol::PartialThreeJSymbol()
{
  this->J1 = -1;
  this->J2 = -1;
  this->M3 = 0;
  this->MinM1=0;
  this->MaxM1=0;
  this->NbrJ3=0;
  this->CurrentPosition = -1;
  this->J3Min = 0l;
  this->Coefficients = 0l;
}


// constructor
//
// j1 = first angular momentum (twice the value to avoid half integer value)
// j2 = second angular momentum (twice the value to avoid half integer value)
// m3 = twice the value of the z-component of the target angular momentum
PartialThreeJSymbol::PartialThreeJSymbol(int j1, int j2, int m3, ThreeJSymbol *fullSymbol)
{
  this->J1 = j1;
  this->J2 = j2;
  this->M3 = m3;
  this->MinM1=-j1;
  this->MaxM1=j1;  
  this->Flag.Initialize();
  this->CurrentPosition = -1;
  if ((fullSymbol==NULL)||(fullSymbol->GetJ1()!=j1)||(fullSymbol->GetJ2()!=j2))
    {
      ThreeJSymbol *TmpSymbol = new ThreeJSymbol(j1, j2);
      this->InitializeTable(*TmpSymbol);
      delete TmpSymbol;
    }
  else this->InitializeTable(*fullSymbol);
}
  
// constructor
//
// j1 = first angular momentum (twice the value to avoid half integer value)
// j2 = second angular momentum (twice the value to avoid half integer value)
// m3 = twice the value of the z-component of the target angular momentum
// minM1 = smallest allowed value for m1
// maxM1 = largest allowed value for m1
PartialThreeJSymbol::PartialThreeJSymbol(int j1, int j2, int m3, int minM1, int maxM1, ThreeJSymbol *fullSymbol)
{
  this->J1 = j1;
  this->J2 = j2;
  this->M3 = m3;
  this->MinM1=minM1;
  if (MinM1<-j1)
    MinM1=-j1;
  this->MaxM1=maxM1;
  if (MaxM1>j1)
    MaxM1=j1;
  if (minM1>maxM1)
    {
      cout << "It is required that minM1 <= maxM1 in PartialThreeJSymbol"<<endl;
      exit(-1);
    }
  this->Flag.Initialize();
  this->CurrentPosition = -1;
  if ((fullSymbol==NULL)||(fullSymbol->GetJ1()!=j1)||(fullSymbol->GetJ2()!=j2))
    {
      ThreeJSymbol *TmpSymbol = new ThreeJSymbol(j1, j2);
      this->InitializeTable(*TmpSymbol);
      delete TmpSymbol;
    }
  else this->InitializeTable(*fullSymbol);
}



// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy

PartialThreeJSymbol::PartialThreeJSymbol (const PartialThreeJSymbol& coefficients)
{
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->M3 = coefficients.M3;
  this->NbrJ3 = coefficients.NbrJ3;
  this->MinM1 = coefficients.MinM1;
  this->MaxM1 = coefficients.MaxM1;
  this->NbrM1 = coefficients.NbrM1;
  this->Coefficients = coefficients.Coefficients;
  this->J3Min = coefficients.J3Min;
  this->CurrentPosition = -1;
  this->Flag = coefficients.Flag;
}

// destructor
//

PartialThreeJSymbol::~PartialThreeJSymbol ()
{
   if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
     {
       for (int i = 0; i < NbrM1; ++i)
	{
	  delete[] this->Coefficients[i];
	}
      delete[] this->Coefficients;
    }
}


// assignment (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to assign
// return value = reference on current Clebsch Gordan coefficients

PartialThreeJSymbol& PartialThreeJSymbol::operator = (const PartialThreeJSymbol& coefficients)
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < NbrM1; ++i)
	{
	  delete[] this->Coefficients[i];
	}
      delete[] this->Coefficients;
    }
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->M3 = coefficients.M3;
  this->NbrJ3 = coefficients.NbrJ3;
  this->MinM1 = coefficients.MinM1;
  this->MaxM1 = coefficients.MaxM1;
  this->NbrM1 = coefficients.NbrM1;
  this->Coefficients = coefficients.Coefficients;
  this->J3Min = coefficients.J3Min;
  this->CurrentPosition = -1;
  this->Flag = coefficients.Flag;
  return *this;
}

// get a particular coefficient (without testing if m1, m2 and j are valid)
//
// m1 = projection of first angular momentum 
// j = resulting angular momentum
// return value = corresponding Clebsch Gordan coefficient

double PartialThreeJSymbol::GetCoefficient (int m1, int j)
{
  int TmpPos1 = (m1-MinM1) >> 1;
  return this->Coefficients[TmpPos1][(j - this->J3Min) >> 1];
}

// initial an iterator on all Clebsch Gordan coefficients for fixed m1 and m2 values
//
// m1 = projection of first angular momentum 

void PartialThreeJSymbol::InitializeCoefficientIterator(int m1)
{
  this->M1 = (m1-MinM1) >> 1;
  this->CurrentPosition = 0;
  this->J = J3Min;
}

// return next coefficient associated with current iterator (with increasing j value)
//
// j = reference on integer where resulting angular momentum has to be stored
// coefficient = reference on double where Clebsch Gordan coefficient has to be stored
// return value = false if no coefficient has been returned

bool PartialThreeJSymbol::Iterate(int& j, double& coefficient)
{
  if (this->CurrentPosition == -1)
    return false;
  coefficient = this->Coefficients[this->M1][this->CurrentPosition];
  j = this->J;
  this->J += 2;
  if (this->J > (this->J1 + this->J2))
    this->CurrentPosition = - 1;
  else
    ++this->CurrentPosition;
  return true;
}

// print a particular coefficient (without testing if m1, m2 and j are valid)
//
// str = reference on output stream
// m1 = projection of first angular momentum 
// j = resulting angular momentum
// return value = reference on output stream

ostream& PartialThreeJSymbol::PrintCoefficient (ostream& str, int m1, int j)
{
  int TmpPos1 = (m1 - MinM1) >> 1;
  str << " < ";
  if ((this->J1 & 0x1) == 0)
    str << (this->J1 >> 1) << " " << (m1 / 2) << " ; ";
  else
    str << this->J1 << "/2 " << m1 << "/2 ; ";
  if ((this->J2 & 0x1) == 0)
    str << (this->J2 >> 1) << " " << ( (-M3-m1) / 2) << " | ";
  else
    str << this->J2 << "/2 " << -M3-m1 << "/2 | ";
  if ((j & 0x1) == 0)
    str << (j >> 1) << " " << (M3 / 2);
  else
    str << j << "/2 " << M3 << "/2";
  str << " > = " << this->Coefficients[TmpPos1][(j - this->J3Min) >> 1];
  return str;
}

// initialize tables, copying values from provided 3J symbol
// fullSymbol = 3J symbol at j1, j2
//
void PartialThreeJSymbol::InitializeTable(ThreeJSymbol &fullSymbol)
{  
  this->J3Min=(J1-J2>0 ? J1-J2 : J2-J1);
  if (this->J3Min<M3) this->J3Min=M3;
  int MaxJ3=J1+J2;
  if (MaxJ3<M3)
    {
      cout << "The requested total angular momentum "<<M3<<" cannot be formed from j1="<<J1<<", and j2="<<J2<<endl;
      exit(-1);
    }
  this->NbrJ3=((MaxJ3-this->J3Min)>>1)+1;
  if (MinM1<-J2-M3)
    MinM1=-J2-M3;
  if (MaxM1>J2-M3)
    MaxM1=J2-M3;
  this->NbrM1 = ((MaxM1-MinM1)>>1)+1;
  this->Coefficients = new double*[NbrM1];  
  for (int pos=0, m1=MinM1; pos<NbrM1; ++pos, m1+=2)
    {
      this->Coefficients[pos] = new double[NbrJ3];
      for (int pos2=0, j3=this->J3Min; pos2<NbrJ3; ++pos2, j3+=2)
	{
	  this->Coefficients[pos][pos2]=fullSymbol.GetCoefficient(m1, -M3-m1, j3);
	}
    }
}
