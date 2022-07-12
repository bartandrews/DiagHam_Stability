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


#include "config.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <stdlib.h>
#include <math.h>


// constructor 
//
// j1 = first angular momentum (twice the value to avoid half integer value)
// j2 = second angular momentum (twice the value to avoid half integer value)

ClebschGordanCoefficients::ClebschGordanCoefficients(int j1, int j2)
{
  this->J1 = j1;
  this->J2 = j2;
  this->Flag.Initialize();
  this->CurrentPosition = -1;
  this->EvaluateClebschGordanCoefficients();
}

// copy constructor (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to copy

ClebschGordanCoefficients::ClebschGordanCoefficients (const ClebschGordanCoefficients& coefficients)
{
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->Coefficients = coefficients.Coefficients;
  this->JMin = coefficients.JMin;
  this->CurrentPosition = -1;
}

// destructor
//

ClebschGordanCoefficients::~ClebschGordanCoefficients ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int TotalM1 = (this->J1 + 1);
      int TotalM2 = (this->J2 + 1);      
      for (int i = 0; i < TotalM1; ++i)
	{
	  for (int j = 0; j < TotalM2; ++j)
	    delete[] this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	  delete[] this->JMin[i];
	}
      delete[] this->JMin;
    }
}

// assignment (without duplicating datas)
//
// coefficients = reference on Clebsch Gordan coefficients to assign
// return value = reference on current Clebsch Gordan coefficients

ClebschGordanCoefficients& ClebschGordanCoefficients::operator = (const ClebschGordanCoefficients& coefficients)
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      int TotalM1 = (this->J1 + 1);
      int TotalM2 = (this->J2 + 1);      
      for (int i = 0; i < TotalM1; ++i)
	{
	  for (int j = 0; j < TotalM2; ++j)
	    delete[] this->Coefficients[i][j];
	  delete[] this->Coefficients[i];
	  delete[] this->JMin[i];
	}
      delete[] this->JMin;
    }
  this->J1 = coefficients.J1;
  this->J2 = coefficients.J2;
  this->Coefficients = coefficients.Coefficients;
  this->JMin = coefficients.JMin;
  this->CurrentPosition = -1;
  return *this;
}

// get a particular coefficient (without testing if m1, m2 and j are valid)
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = corresponding Clebsch Gordan coefficient

double ClebschGordanCoefficients::GetCoefficient (int m1, int m2, int j)
{
  int TmpPos1 = (this-> J1 + m1) >> 1;
  int TmpPos2 = (this-> J2 + m2) >> 1;
  return this->Coefficients[TmpPos1][TmpPos2][(j - this->JMin[TmpPos1][TmpPos2]) >> 1];
}

// initial an iterator on all Clebsch Gordan coefficients for fixed m1 and m2 values
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum   

void ClebschGordanCoefficients::InitializeCoefficientIterator(int m1, int m2)
{
  this->M1 = (this->J1 + m1) >> 1;
  this->M2 = (this->J2 + m2) >> 1;
  this->CurrentPosition = 0;
  this->J = this->JMin[this->M1][this->M2];
}

// return next coefficient associated with current iterator (with increasing j value)
//
// j = reference on integer where resulting angular momentum has to be stored
// coefficient = reference on double where Clebsch Gordan coefficient has to be stored
// return value = false if no coefficient has been returned

bool ClebschGordanCoefficients::Iterate(int& j, double& coefficient)
{
  if (this->CurrentPosition == -1)
    return false;
  coefficient = this->Coefficients[this->M1][this->M2][this->CurrentPosition];
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
// m2 = projection of second angular momentum 
// j = resulting angular momentum
// return value = reference on output stream

ostream& ClebschGordanCoefficients::PrintCoefficient (ostream& str, int m1, int m2, int j)
{
  int TmpPos1 = (this-> J1 + m1) >> 1;
  int TmpPos2 = (this-> J2 + m2) >> 1;
  str << " < ";
  if ((this->J1 & 0x1) == 0)
    str << (this->J1 >> 1) << " " << (m1 / 2) << " ; ";
  else
    str << this->J1 << "/2 " << m1 << "/2 ; ";
  if ((this->J2 & 0x1) == 0)
    str << (this->J2 >> 1) << " " << (m2 / 2) << " | ";
  else
    str << this->J2 << "/2 " << m2 << "/2 | ";
  if ((j & 0x1) == 0)
    str << (j >> 1) << " " << ((m1 + m2) / 2);
  else
    str << j << "/2 " << (m1 + m2) << "/2";
  str << " > = " << this->Coefficients[TmpPos1][TmpPos2][(j - this->JMin[TmpPos1][TmpPos2]) >> 1];
  return str;
}

// evaluate all Clebsch Gordan coefficients using Schulten Gordon recursion algorithm
//

void ClebschGordanCoefficients::EvaluateClebschGordanCoefficients()
{
  int TotalM1 = (this->J1 + 1);
  int TotalM2 = (this->J2 + 1);      
  int SumM;
  int SumJ = this->J1 + this->J2;
  int DiffJ = abs(this->J1 - this->J2);
  this->JMin = new int* [TotalM1];
  this->Coefficients = new double** [TotalM1];
  for (int i = 0; i < TotalM1; ++i)
    {
      this->JMin[i] = new int [TotalM2];
      this->Coefficients[i] = new double* [TotalM2];
      for (int j = 0; j < TotalM2; ++j)
	{
	  SumM = abs (2 * (i + j) - SumJ);
	  if (SumM > DiffJ)
	    this->JMin[i][j] = SumM;
	  else
	    this->JMin[i][j] = DiffJ;
	  this->Coefficients[i][j] = new double [SumJ - this->JMin[i][j] + 1];
	}
    }

  int NbrValue;
  double* TmpCoefficients;
  int M2Max;
  int M2Min;
  int* TmpJMin;
  for (int m1 = -this->J1; m1 <= this->J1; m1 += 2)
    {
      for (int j = DiffJ; j <= SumJ; j += 2)
	{
	  M2Max = j - m1;
	  if(M2Max > this->J2)
	    M2Max = this->J2;
	  M2Min = -j - m1;
	  if(M2Min < (-this->J2))
	    M2Min = -this->J2;
	  if (M2Max == M2Min)
	    {
	      int M1Pos = (m1 + this->J1) >> 1;
	      int M2Pos = (M2Min + this->J2) >> 1;
	      int Sign = (this->J2 - m1 - j) / 2;
	      if ((Sign & 0x1) == 0)
		{
		  this->Coefficients[M1Pos][M2Pos][(j - this->JMin[M1Pos][M2Pos]) >> 1] = 1.0 / sqrt (this->J1 + 1.0);
		}
	      else
		{
		  this->Coefficients[M1Pos][M2Pos][(j - this->JMin[M1Pos][M2Pos]) >> 1] = -1.0 / sqrt (this->J1 + 1.0);
		}
	    }
	  else
	    if (M2Max > M2Min)
	      {
		this->RecursionAlgorithm(m1, M2Min, M2Max, j);
	      }
	}
    }
  double J;
  double SignCoef = 1.0;
  double SignCoef2;
  if ((this->J2 & 0x1) == 0x1)
    SignCoef = -1;
  for (int i = 0; i < TotalM1; ++i)
    {
      SignCoef2 = SignCoef;
      TmpJMin = this->JMin[i];
      for (int j = 0; j < TotalM2; ++j)
	{
	  TmpCoefficients = this->Coefficients[i][j];
	  NbrValue = TmpJMin[j];
	  J = (double) NbrValue;
	  NbrValue = (SumJ - NbrValue) >> 1;
	  for (int k = 0; k <= NbrValue; ++k)
	    {
	      TmpCoefficients[k] *= SignCoef2 * sqrt(J + 1.0);
	      J += 2.0;
	    }
	  SignCoef2 *= -1.0;
	}
      SignCoef *= -1.0;
    }
}

// evaluate C coefficient needed during the resursion
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum   
// j = resulting angular momentum
// return value = corresponding D coefficient

double ClebschGordanCoefficients::CCoefficient (int m1, int m2, int j)
{
  double C = (double) (this->J2 - m2 + 2);
  C *= (double) (this->J2 + m2);
  C *= (double) (j + (m1 + m2));
  C *= (double) (j - (m1 + m2) + 2);
  return (0.25 * sqrt(C)); 
}

// evaluate D coefficient needed during the resursion
//
// m1 = projection of first angular momentum 
// m2 = projection of second angular momentum   
// j = resulting angular momentum
// return value = corresponding D coefficient

double ClebschGordanCoefficients::DCoefficient (int m1, int m2, int j)
{
  return (0.25 * (double) (this->J2 * (this->J2 + 2) + j * (j + 2) - this->J1 * (this->J1 + 2) - 2 * (m2 * (m1 + m2)))); 
}

// main part of the recursion algorithm
//
// m1 = projection of first angular momentum 
// m2Min = minimum value for the projection of second angular momentum   
// m2Max = maximum value for the projection of second angular momentum   
// j = resulting angular momentum

void ClebschGordanCoefficients::RecursionAlgorithm (int m1, int m2Min, int m2Max, int j)
{
  int M2Mid = (m2Max + m2Min) / 2;
  int M1Pos = (m1 + this->J1) >> 1;
  int M2Pos = (m2Min + this->J2) >> 1;
  double**  TmpCoefficients2 = this->Coefficients[M1Pos];
  int* TmpJMin = this->JMin[M1Pos];
  double PreviousWigner3jCoefficient = 1.0;
  double CurrentWigner3jCoefficient = - (DCoefficient(m1, m2Min, j) / CCoefficient(m1, m2Min + 2, j)); 
  double NextWigner3jCoefficient = 0.0;
  if ((m2Max - m2Min) == 2)
    {
      double NormalizationFactor = 1.0 / sqrt((this->J1 + 1) * (1.0 + CurrentWigner3jCoefficient * CurrentWigner3jCoefficient));
      int Sign = (this->J2 - m1 - j) / 2;
      if ((((Sign & 0x1) == 0) && (CurrentWigner3jCoefficient < 0)) || (((Sign & 0x1) == 1) && (CurrentWigner3jCoefficient > 0)))
	NormalizationFactor *= -1.0;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NormalizationFactor;
      ++M2Pos;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = CurrentWigner3jCoefficient * NormalizationFactor;
      return;
    }
  if ((m2Max - m2Min) == 4)
    {
      NextWigner3jCoefficient = - ((DCoefficient(m1, M2Mid, j) * CurrentWigner3jCoefficient + 
				    CCoefficient(m1, M2Mid, j) * PreviousWigner3jCoefficient) / CCoefficient(m1, m2Max, j));
      double NormalizationFactor = 1.0 / sqrt((this->J1 + 1) * (1.0 + NextWigner3jCoefficient * NextWigner3jCoefficient +
								CurrentWigner3jCoefficient * CurrentWigner3jCoefficient));
      int Sign = (this->J2 - m1 - j) / 2;
      if ((((Sign & 0x1) == 0) && (NextWigner3jCoefficient < 0)) || (((Sign & 0x1) == 1) && (NextWigner3jCoefficient > 0)))
	NormalizationFactor *= -1.0;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NormalizationFactor;
      ++M2Pos;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = CurrentWigner3jCoefficient * NormalizationFactor;
      ++M2Pos;      
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NextWigner3jCoefficient * NormalizationFactor;
      return;
    }
  if ((m2Max - m2Min) == 6)
    {
      NextWigner3jCoefficient = - ((DCoefficient(m1, m2Min + 2, j) * CurrentWigner3jCoefficient + 
				    CCoefficient(m1, m2Min + 2, j) * PreviousWigner3jCoefficient) / CCoefficient(m1, m2Max - 2, j));
      double LastWigner3jCoefficient = - ((DCoefficient(m1, m2Max - 2, j) * NextWigner3jCoefficient + 
					   CCoefficient(m1, m2Max - 2, j) * CurrentWigner3jCoefficient) / CCoefficient(m1, m2Max, j));
      double NormalizationFactor = 1.0 / sqrt((this->J1 + 1) * (1.0 + NextWigner3jCoefficient * NextWigner3jCoefficient +
								CurrentWigner3jCoefficient * CurrentWigner3jCoefficient + 
								LastWigner3jCoefficient * LastWigner3jCoefficient));
      int Sign = (this->J2 - m1 - j) / 2;
      if ((((Sign & 0x1) == 0) && (LastWigner3jCoefficient < 0)) || (((Sign & 0x1) == 1) && (LastWigner3jCoefficient > 0)))
	NormalizationFactor *= -1.0;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NormalizationFactor;
      ++M2Pos;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = CurrentWigner3jCoefficient * NormalizationFactor;
      ++M2Pos;      
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NextWigner3jCoefficient * NormalizationFactor;
      ++M2Pos;            
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = LastWigner3jCoefficient * NormalizationFactor;
      return;
    }
  if ((m2Max & 0x1) != (M2Mid & 0x1))
    ++M2Mid;
  double Normalization = PreviousWigner3jCoefficient * PreviousWigner3jCoefficient;
  TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = PreviousWigner3jCoefficient;
  ++M2Pos;
  Normalization += CurrentWigner3jCoefficient * CurrentWigner3jCoefficient;
  TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = CurrentWigner3jCoefficient;
  ++M2Pos;
  for (int m2 = m2Min + 2; m2 < M2Mid; m2 += 2)
    {	
      NextWigner3jCoefficient = - ((DCoefficient(m1, m2, j) * CurrentWigner3jCoefficient + 
				    CCoefficient(m1, m2, j) * PreviousWigner3jCoefficient) / CCoefficient(m1, m2 + 2, j));
      Normalization += NextWigner3jCoefficient * NextWigner3jCoefficient;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NextWigner3jCoefficient;
      ++M2Pos;
      PreviousWigner3jCoefficient = CurrentWigner3jCoefficient;
      CurrentWigner3jCoefficient = NextWigner3jCoefficient;
    }
  double MidWigner3jCoefficient = CurrentWigner3jCoefficient;
  if (fabs(MidWigner3jCoefficient) > MACHINE_PRECISION)
    {
      M2Pos = (m2Max + this->J2) >> 1;
      PreviousWigner3jCoefficient = 1.0;
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = PreviousWigner3jCoefficient;
      --M2Pos;
      CurrentWigner3jCoefficient = - (DCoefficient(m1, m2Max, j) / CCoefficient(m1, m2Max, j)); 
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = CurrentWigner3jCoefficient;
      --M2Pos;
      int ReducedM2Mid = M2Mid + 2;
      for (int m2 = m2Max - 2; m2 > ReducedM2Mid; m2 -= 2)
	{	
	  NextWigner3jCoefficient = - ((DCoefficient(m1, m2, j) * CurrentWigner3jCoefficient + 
					CCoefficient(m1, m2 + 2, j) * PreviousWigner3jCoefficient) / CCoefficient(m1, m2, j));
	  
	  TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NextWigner3jCoefficient;
	  --M2Pos;
	  PreviousWigner3jCoefficient = CurrentWigner3jCoefficient;
	  CurrentWigner3jCoefficient = NextWigner3jCoefficient;
	}
      NextWigner3jCoefficient = - ((DCoefficient(m1, ReducedM2Mid, j) * CurrentWigner3jCoefficient + 
				    CCoefficient(m1, ReducedM2Mid + 2, j) * PreviousWigner3jCoefficient) 
				   / CCoefficient(m1, ReducedM2Mid, j));
      ++M2Pos;
      CurrentWigner3jCoefficient = MidWigner3jCoefficient / NextWigner3jCoefficient;
      for (int m2 = ReducedM2Mid; m2 <= m2Max;  m2 += 2)
	{
	  TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] *= CurrentWigner3jCoefficient;
	  NextWigner3jCoefficient = TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1];
	  ++M2Pos;
	  Normalization += NextWigner3jCoefficient * NextWigner3jCoefficient;
	}
    }
  else
    {
      for (int m2 = M2Mid; m2 < m2Max; m2 += 2)
	{	
	  NextWigner3jCoefficient = - ((DCoefficient(m1, m2, j) * CurrentWigner3jCoefficient + 
					CCoefficient(m1, m2, j) * PreviousWigner3jCoefficient) / CCoefficient(m1, m2 + 2, j));
	  Normalization += NextWigner3jCoefficient * NextWigner3jCoefficient;
	  TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] = NextWigner3jCoefficient;
	  ++M2Pos;
	  PreviousWigner3jCoefficient = CurrentWigner3jCoefficient;
	  CurrentWigner3jCoefficient = NextWigner3jCoefficient;
	}
    }
  Normalization = 1.0 / sqrt (Normalization * (this->J1 + 1.0));
  int Sign = (this->J2 - m1 - j) / 2;
  if ((((Sign & 0x1) == 0) && (NextWigner3jCoefficient < 0)) || (((Sign & 0x1) == 1) && (NextWigner3jCoefficient > 0)))
    Normalization *= -1.0;
  M2Pos = (m2Min + this->J2) >> 1;
  for (int m2 = m2Min; m2 <= m2Max; m2 += 2)
    {	
      TmpCoefficients2[M2Pos][(j - TmpJMin[M2Pos]) >> 1] *= Normalization;
      ++M2Pos;
    }
}
