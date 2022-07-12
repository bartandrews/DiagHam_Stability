////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of MPS matrix for the Laughlin state                //
//                                                                            //
//                        last modification : 30/10/2012                      //
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
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix()
{
  this->UniformChargeIndexRange = true;
  this->BosonicVersion = false;
  this->TorusFlag = false;
  this->SiteDependentMatrixNbrOrbitals = 0;
  this->SiteDependentMatrices = 0;
  this->SiteDependentMatrixOrbitalIndices = 0;
  this->NbrSiteDependentMatrices = 0;
  this->SiteDependentPhysicalIndices = 0;
}

// constructor 
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->ComplexBMatrices = 0;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->BosonicVersion = bosonicVersion;
  this->TorusFlag = false;
  this->TwistedTorusFlag = false;
  if (this->BosonicVersion == true)
    this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex) + 1;
  else
    this->NbrNValue = (2 * this->PLevel) + this->LaughlinIndex;
  this->NValueGlobalShift = this->PLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  this->SiteDependentMatrixNbrOrbitals = 0;
  this->SiteDependentMatrices = 0;
  this->SiteDependentMatrixOrbitalIndices = 0;
  this->NbrSiteDependentMatrices = 0;
  this->SiteDependentPhysicalIndices = 0;
  if (this->BosonicVersion == true)
    this->AlternateCreateBMatrices();
  else
    this->CreateBMatrices();
  this->ComputeSiteDependentMatrices(-8, 3);
}

// constructor for the torus geometry
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// trimChargeIndices = trim the charge indices
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool trimChargeIndices, 
					     int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion)
{
  this->NbrBMatrices = nbrBMatrices;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->BosonicVersion = bosonicVersion;
  if (this->BosonicVersion == true)
    this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex) + 1;
  else
    this->NbrNValue = (2 * this->PLevel) + this->LaughlinIndex;
  this->NValueGlobalShift = this->PLevel;
  this->CylinderFlag = false;
  this->TorusFlag = true;
  this->TorusNbrFluxQuanta = nbrFluxQuanta;
  this->TorusAngle = angle;
  this->TorusAspectRatio = aspectRatio;
  this->TorusFluxInsertion = fluxInsertion;
  this->Kappa = sqrt(2.0 * M_PI * this->TorusAspectRatio / ((double) this->TorusNbrFluxQuanta));
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->SiteDependentMatrixNbrOrbitals = 0;
  this->SiteDependentMatrices = 0;
  this->SiteDependentMatrixOrbitalIndices = 0;
  this->NbrSiteDependentMatrices = 0;
  this->SiteDependentPhysicalIndices = 0;
  if ((this->TorusAngle != 0.0) || (this->TorusFluxInsertion != 0.0))
    {
      this->TwistedTorusFlag = true;
      this->TauFactor = ((2.0 * M_PI * this->TorusAspectRatio  / ((double) this->TorusNbrFluxQuanta))
			 * Complex(-sin(this->TorusAngle * M_PI), cos(this->TorusAngle * M_PI)));
      this->ComplexBMatrices = new SparseComplexMatrix [this->NbrBMatrices];
      this->RealBMatrices = 0;
    }
  else
    {
      this->TwistedTorusFlag = false;
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      this->ComplexBMatrices = 0;
    }
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  if (this->BosonicVersion == true)
    this->AlternateCreateBMatrices();
  else
    this->CreateBMatrices();
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// trimChargeIndices = trim the charge indices, assuming an iMPS
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinMatrix::FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  this->NValueGlobalShift = this->PLevel;
  this->CylinderFlag = cylinderFlag;
  this->TorusFlag = false;
  this->TwistedTorusFlag = false;
  this->Kappa = kappa;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->SiteDependentMatrixNbrOrbitals = 0;
  this->SiteDependentMatrices = 0;
  this->SiteDependentMatrixOrbitalIndices = 0;
  this->NbrSiteDependentMatrices = 0;
  this->SiteDependentPhysicalIndices = 0;
  this->LoadMatrices(fileName);
}

// destructor
//

FQHEMPSLaughlinMatrix::~FQHEMPSLaughlinMatrix()
{
  delete[] this->TotalStartingIndexPerPLevel;
  delete[] this->NbrIndicesPerPLevel;
  delete[] this->NbrNValuesPerPLevel;
  delete[] this->NInitialValuePerPLevel;
  delete[] this->NLastValuePerPLevel;
  if (this->SiteDependentMatrixNbrOrbitals != 0)
    {
      for (int i = 0; i < this->SiteDependentMatrixNbrOrbitals; ++i)
	{
	  delete[] this->SiteDependentMatrices[i];
	  delete[] this->SiteDependentPhysicalIndices[i];
	}
      delete[] this->SiteDependentMatrices;
      delete[] this->SiteDependentMatrixOrbitalIndices;
      delete[] this->NbrSiteDependentMatrices;
      delete[] this->SiteDependentPhysicalIndices;
    }
}
  
// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSLaughlinMatrix::GetName()
{
  char* TmpName = new char[16];
  sprintf (TmpName, "laughlin%d", this->LaughlinIndex);
  return TmpName;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSLaughlinMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 1;
  denominator = this->LaughlinIndex;
}

// create the B matrices for the laughlin state
//

void FQHEMPSLaughlinMatrix::CreateBMatrices ()
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
    }
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValuesPerPLevel[0];
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValuesPerPLevel[i];
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << "x" << MatrixSize << endl;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    ++TmpNbrElementPerRow[this->GetMatrixIndex(i, k, j - 1)];
	}
    }
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
	      double Tmp = 1.0;
	      if (this->CylinderFlag)
// symmetric normalization
// 		Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
// 							 + ((j - 1.0 - this->NValueGlobalShift) * (j - 1.0 - this->NValueGlobalShift) / (4.0 * (double) this->LaughlinIndex))
// 							 + (((j - this->NValueGlobalShift) * (j - this->NValueGlobalShift)) / (4.0 * (double) this->LaughlinIndex))));
		Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
							 + (((j - 1.0 - this->NValueGlobalShift) * (j - 1.0 - this->NValueGlobalShift)) / (2.0 * (double) this->LaughlinIndex))));
	      BMatrices[0].SetMatrixElement(this->GetMatrixIndex(i, k, j - 1), this->GetMatrixIndex(i, k, j), Tmp);
	    }
	}
    }

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;
  
  for (int m = 1; m < this->NbrBMatrices; ++m)
    {
      for (int i = 0; i < MatrixSize; ++i)
	TmpNbrElementPerRow[i] = 0;
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	      int N2 = (j - i) + this->NValueGlobalShift;
	      int N1 = N2 + (this->LaughlinIndex - 1);
  	      if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
  		  && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		{ 
		  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		    {
		      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
			{
			  ++TmpNbrElementPerRow[this->GetMatrixIndex(i, k1, N1)];
			}
		    }
		}
	    }
	}

      BMatrices[m] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);

      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	      int N2 = (j - i) + this->NValueGlobalShift;
	      int N1 = N2 + (this->LaughlinIndex - 1);
  	      if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
  		  && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		{ 
		  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		    {
		      TmpSpace1->GetOccupationNumber(k1, Partition1);
		      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
			{
			  TmpSpace2->GetOccupationNumber(k2, Partition2);
			  double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex * m * m, 1, Partition1, Partition2, i, j, Coef);
			  if (this->CylinderFlag)
// symmetric normalization
// 			    Tmp *= exp(-this->Kappa * this->Kappa * (( 0.5 *  ((double) (i + j)))
// 								     + ((N1 - this->NValueGlobalShift) * (N1 - this->NValueGlobalShift)  / (4.0 * (double) this->LaughlinIndex))
// 								     + (((N2 - this->NValueGlobalShift) * (N2 - this->NValueGlobalShift))  / (4.0 * (double) this->LaughlinIndex))));
			    Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
								     + (((N1 - this->NValueGlobalShift) * (N1 - this->NValueGlobalShift))  / (2.0 * (double) this->LaughlinIndex))));
			  BMatrices[m].SetMatrixElement(this->GetMatrixIndex(i, k1, N1), this->GetMatrixIndex(j, k2, N2), Tmp);
			}
		    }
		}
	    }
	}
    }
  
  delete[] Partition1;
  delete[] Partition2;

  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;
}

// create the B matrices for the laughlin state
//

void FQHEMPSLaughlinMatrix::AlternateCreateBMatrices ()
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  SparseRealMatrix* BMatrices = 0;
  SparseComplexMatrix* TmpComplexBMatrices = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
    }
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValuesPerPLevel[0];
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValuesPerPLevel[i];
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << "x" << MatrixSize << endl;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    ++TmpNbrElementPerRow[this->GetMatrixIndex(i, k, j - 1)];
	}
    }
  int* TmpNbrElementPerRowInverseBMatrixZero = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRowInverseBMatrixZero[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = this->NInitialValuePerPLevel[i]; j < this->NLastValuePerPLevel[i]; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    ++TmpNbrElementPerRowInverseBMatrixZero[this->GetMatrixIndex(i, k, j)];
	}
    }
  if (this->TwistedTorusFlag == false)
    {
      BMatrices = new SparseRealMatrix[this->NbrBMatrices];
      BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
      InverseBMatrixZero = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
	  for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	    {
	      for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
		{
		  double Tmp = 1.0;
		  if (this->CylinderFlag)
		    {
		      Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
							       + (((j - 1.0 - this->NValueGlobalShift) * (j - 1.0 - this->NValueGlobalShift)) / (2.0 * (double) this->LaughlinIndex))));
		    }
		  else
		    {
		      if (this->TorusFlag)
			{
			  Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
								   + (((j - 1.0 - this->NValueGlobalShift) * (j - 1.0 - this->NValueGlobalShift)) / (2.0 * (double) this->LaughlinIndex))));
			}
		    }
		  BMatrices[0].SetMatrixElement(this->GetMatrixIndex(i, k, j - 1), this->GetMatrixIndex(i, k, j), Tmp);
		  if (Tmp != 0.0)
		    {
		      InverseBMatrixZero.SetMatrixElement(this->GetMatrixIndex(i, k, j), this->GetMatrixIndex(i, k, j - 1), 1.0 / Tmp);
		    }
		}
	    }
	}
    }
  else
    {
      TmpComplexBMatrices = new SparseComplexMatrix[this->NbrBMatrices];
      TmpComplexBMatrices[0] = SparseComplexMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
	  for (int j = this->NInitialValuePerPLevel[i] + 1; j <= this->NLastValuePerPLevel[i]; ++j)
	    {
	      for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
		{
		  Complex Tmp = exp(this->TauFactor * (((double) i)
						       + (((j - 1.0 - this->NValueGlobalShift) * (j - 1.0 - this->NValueGlobalShift)) / (2.0 * (double) this->LaughlinIndex))));
		  TmpComplexBMatrices[0].SetMatrixElement(this->GetMatrixIndex(i, k, j - 1), this->GetMatrixIndex(i, k, j), Tmp);
		}
	    }
	}
    }

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;
  
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	  int N2 = (j - i) + this->NValueGlobalShift;
	  int N1 = N2 + (this->LaughlinIndex);
	  if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
  		  && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
	    { 
	      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		{
		  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		    {
		      ++TmpNbrElementPerRow[this->GetMatrixIndex(i, k1, N1)];
		    }
		}
	    }
	}
    }
  
  SparseRealMatrix V0Matrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	  int N2 = (j - i) + this->NValueGlobalShift;
	  int N1 = N2 + (this->LaughlinIndex);
	  if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
	      && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
	    { 
	      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
		{
		  TmpSpace1->GetOccupationNumber(k1, Partition1);
		  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		    {
		      TmpSpace2->GetOccupationNumber(k2, Partition2);
		      double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex, 1, Partition1, Partition2, i, j, Coef);
		      V0Matrix.SetMatrixElement(this->GetMatrixIndex(i, k1, N1), this->GetMatrixIndex(j, k2, N2), Tmp);
		    }
		}
	    }
	}
    }
  
  delete[] Partition1;
  delete[] Partition2;

  int TmpNbrBMatrices = 0;
  if (this->TwistedTorusFlag == false)
    {
      for (int m = 1; m < this->NbrBMatrices; ++m)
	{
	  BMatrices[m] = MemoryEfficientMultiply(BMatrices[m - 1], V0Matrix);
	  if (BMatrices[m].GetNbrRow() > 0)
	    {
	      ++TmpNbrBMatrices;
	      if ((this->CylinderFlag) || (this->TorusFlag))
		BMatrices[m] /= sqrt((double) m);
	    }
	}
    }
  else
    {
      for (int m = 1; m < this->NbrBMatrices; ++m)
	{
	  TmpComplexBMatrices[m] = MemoryEfficientMultiply(TmpComplexBMatrices[m - 1], V0Matrix);
	  if (TmpComplexBMatrices[m].GetNbrRow() > 0)
	    {
	      ++TmpNbrBMatrices;
	      TmpComplexBMatrices[m] /= sqrt((double) m);
	    }
	}
    }
  TmpNbrBMatrices = 0;
  if (this->TwistedTorusFlag == false)
    {
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (BMatrices[i].GetNbrRow() > 0)
	    this->RealBMatrices[TmpNbrBMatrices++] = BMatrices[i];
	}
      delete[] BMatrices;
    }
  else
    {
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (TmpComplexBMatrices[i].GetNbrRow() > 0)
	    this->ComplexBMatrices[TmpNbrBMatrices++] = TmpComplexBMatrices[i];
	}
      delete[] TmpComplexBMatrices;
    }
  this->NbrBMatrices = TmpNbrBMatrices;
  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;
}

// create the matrix element of the B matrix U(1) part
//
// chargeNumerator = numerator of the charge (in sqrt(q) unit)
// chargeDenominator = denominator of the charge (in sqrt(q) unit)
// partition1 = U(1) partition associated to the left state
// p1Level = length of partition1
// partition2 = U(1) partition associated to the left state
// p1Level = length of partition2
// coef = reference on a temporary factorial coefficient
// return value = matrix element

double FQHEMPSLaughlinMatrix::CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, 
							    unsigned long* partition1, unsigned long* partition2, 
							    int p1Level, int p2Level, FactorialCoefficient& coef)
{
  double Tmp = 1.0;
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      double Tmp2 = 0.0;
      for (int j = 0; j <= partition1[i]; ++j)
	{
	  int k = partition2[i] + j - partition1[i];
	  if ((k >= 0) && (k <= partition2[i]))
	    {
	      int Sum = k + j;
	      coef.SetToOne();
	      coef.PartialFactorialMultiply(partition1[i] - j + 1, partition1[i]);
	      coef.PartialFactorialMultiply(partition2[i] - k + 1, partition2[i]);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.PowerNMultiply(chargeNumerator, Sum);
	      coef.PowerNDivide(chargeDenominator, Sum);
	      coef.PowerNDivide(i, Sum);
	      Tmp2 += ((double) (1 - ((j & 1) << 1))) * sqrt(coef.GetNumericalValue());
	    }
	}
      Tmp *= Tmp2;
    }
  return Tmp;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = tuncation level of the block left indices
// q1 = charge index of the block left indices
// pLevel1 = tuncation level of the block right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

SparseRealMatrix FQHEMPSLaughlinMatrix::ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2)
{
  double Tmp;

  int NbrK1 = this->GetBondIndexRange(pLevel1, q1);
  int NbrK2 = this->GetBondIndexRange(pLevel2, q2);
  SparseRealMatrix TmpMatrix (NbrK1, NbrK2);
  for (int k1 = 0; k1 < NbrK1; ++k1)
    {
      for (int k2 = 0; k2 < NbrK2; ++k2)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargeAndPLevel(k1, pLevel1, q1), 
				  this->GetBondIndexWithFixedChargeAndPLevel(k2, pLevel2, q2), Tmp);
	  if (Tmp != 0.0)
	    TmpMatrix.SetMatrixElement(k1, k2, Tmp);
	}
    }

  return TmpMatrix;
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSLaughlinMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevel[pLevel]) || (qValue > this->NLastValuePerPLevel[pLevel]))
    return 0;
  return this->NbrIndicesPerPLevel[pLevel] / this->NbrNValuesPerPLevel[pLevel];  
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSLaughlinMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{  
  return (this->TotalStartingIndexPerPLevel[pLevel] + (localIndex * this->NbrNValuesPerPLevel[pLevel] + (qValue - this->NInitialValuePerPLevel[pLevel])));
}

// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSLaughlinMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevel[pLevel];
  maxQ = this->NLastValuePerPLevel[pLevel];
  return;
}

// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of particles

int FQHEMPSLaughlinMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  int Numerator;
  int Denominator;
  this->GetFillingFactor(Numerator, Denominator);
  if (this->TorusFlag == false)
    {
      int NbrParticles = ((nbrFluxQuanta + 1) * Numerator);
      if ((NbrParticles % Denominator) == 0)
	return (NbrParticles / Denominator);
      else
	return ((NbrParticles / Denominator) + 1);
    }
  else
    {
      return ((nbrFluxQuanta * Numerator) / Denominator);      
    }
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool FQHEMPSLaughlinMatrix::LoadHeader (ifstream& file)
{
  int HeaderSize = 0;
  ReadLittleEndian(file, HeaderSize);
  ReadLittleEndian(file, this->PLevel);
  ReadLittleEndian(file, this->LaughlinIndex);
  ReadLittleEndian(file, this->NbrNValue);
  ReadLittleEndian(file, this->NValueGlobalShift);
  ReadLittleEndian(file, this->CylinderFlag);
  ReadLittleEndian(file, this->Kappa);
  ReadLittleEndian(file, this->UniformChargeIndexRange);
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NbrNValuesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NInitialValuePerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      ReadLittleEndian(file, this->NLastValuePerPLevel[i]);
    }
  return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool FQHEMPSLaughlinMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = (this->PLevel + 1) * (5 * sizeof(int)) + (sizeof(int) * 5) + (sizeof(bool) * 2) + sizeof(double);
  WriteLittleEndian(file, HeaderSize);
  WriteLittleEndian(file, this->PLevel);
  WriteLittleEndian(file, this->LaughlinIndex);
  WriteLittleEndian(file, this->NbrNValue);
  WriteLittleEndian(file, this->NValueGlobalShift);
  WriteLittleEndian(file, this->CylinderFlag);
  WriteLittleEndian(file, this->Kappa);
  WriteLittleEndian(file, this->UniformChargeIndexRange);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->TotalStartingIndexPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NbrIndicesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NbrNValuesPerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NInitialValuePerPLevel[i]);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      WriteLittleEndian(file, this->NLastValuePerPLevel[i]);
    }
  return true;
}

// 2 times the maximal P value in the 'root partition' with a given charge N at a particular cut
// N takes value in [-Pmax, Pmax + q - 1]
static int twop(int N, int q)
{
    // (N // (q - 1)), python-style downward rounding
    int r = (N - ((N < 0)?(q-2):0)) / (q - 1);
    return r * (2 * N - (q - 1) * (r + 1));
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSLaughlinMatrix::ComputeChargeIndexRange(int pLevel, int& minQ, int& maxQ)
{
  if (this->UniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }
  
  if (this->BosonicVersion == false)
    {
      for (minQ = 0; minQ < this->NbrNValue; ++minQ)
	if (2 * pLevel + twop(minQ - this->NValueGlobalShift, this->LaughlinIndex) <= 2 * this->PLevel)
	  break;
      for (maxQ = this->NbrNValue - 1; maxQ >= 0; --maxQ)
	if (2 * pLevel + twop(maxQ - this->NValueGlobalShift, this->LaughlinIndex) <= 2 * this->PLevel)
	  break;
//      --minQ;
//      ++maxQ;
      cout << "N range at " << pLevel << ": [" << minQ - this->NValueGlobalShift << "," << maxQ - this->NValueGlobalShift << "] (+" << this->NValueGlobalShift << ")" << endl;
      
      cout << "other method" << endl;
      int TmpMinQ = this->NbrNValue - 1;
      int TmpMaxQ = 0;    
      int NValueShift = this->PLevel;
      for (int Q = 0; Q < this->NbrNValue; ++Q)
	{
	  int QPrime = Q;
	  int TmpP = 0;
	  int TmpMaxP = -1;
	  QPrime -= (this->LaughlinIndex - 1);
	  TmpP += QPrime - NValueShift;
	  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
	    {
	      if (TmpP > TmpMaxP)
		TmpMaxP = TmpP;	    
	      QPrime -= (this->LaughlinIndex - 1);
	      TmpP += QPrime - NValueShift;
	    }
	  QPrime = Q;
	  TmpP = 0;
	  int TmpMaxP2 = -1;
	  TmpP -= QPrime - NValueShift;
	  QPrime += (this->LaughlinIndex - 1);
	  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
	    {
	      if (TmpP > TmpMaxP2)
		TmpMaxP2 = TmpP;	    
	      TmpP -= QPrime - NValueShift;
	      QPrime += (this->LaughlinIndex - 1);
	    }
	  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
	    {
	      if (Q < TmpMinQ)
		TmpMinQ = Q;
	      if (Q > TmpMaxQ)
		TmpMaxQ = Q;	    
	    }
	}
      cout << "range at " << pLevel << " : " << (TmpMinQ - this->NValueGlobalShift) << " " << (TmpMaxQ - this->NValueGlobalShift) << " (" << this->NbrNValue << ")" << endl;   
    }
  else
    {
      int TmpMinQ = this->NbrNValue - 1;
      int TmpMaxQ = 0;    
      int NValueShift = this->PLevel;
      for (int Q = 0; Q < this->NbrNValue; ++Q)
	{
	  int QPrime = Q;
	  int TmpP = 0;
	  int TmpMaxP = -1;
	  QPrime -= ((this->LaughlinIndex * this->NbrBMatrices) - 1);
	  TmpP += QPrime - NValueShift;
	  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
	    {
	      if (TmpP > TmpMaxP)
		TmpMaxP = TmpP;	    
	      QPrime -= ((this->LaughlinIndex * this->NbrBMatrices) - 1);
	      TmpP += QPrime - NValueShift;
	    }
	  QPrime = Q;
	  TmpP = 0;
	  int TmpMaxP2 = -1;
	  TmpP -= QPrime - NValueShift;
	  QPrime += ((this->LaughlinIndex * this->NbrBMatrices) - 1);
	  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
	    {
	      if (TmpP > TmpMaxP2)
		TmpMaxP2 = TmpP;	    
	      TmpP -= QPrime - NValueShift;
	      QPrime += ((this->LaughlinIndex * this->NbrBMatrices) - 1);
	    }
	  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
	    {
	      if (Q < TmpMinQ)
		TmpMinQ = Q;
	      if (Q > TmpMaxQ)
		TmpMaxQ = Q;	    
	    }
	}
      TmpMaxQ  += 1;
      if (TmpMaxQ >= this->NbrNValue)
	TmpMaxQ = this->NbrNValue - 1;
      minQ = TmpMinQ;
      maxQ = TmpMaxQ;
      cout << "range at " << pLevel << " : " << (TmpMinQ - this->NValueGlobalShift) << " " << (TmpMaxQ - this->NValueGlobalShift) << " (" << this->NbrNValue << ")" << endl;   
    }
}

// compute the global charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSLaughlinMatrix::ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ)
{
  minQ = 0;
  maxQ = this->NbrNValue - 1;
  return;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSLaughlinMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      rowIndex = this->NValueGlobalShift + (this->LaughlinIndex - 1) / 2 - MinQ;
      columnIndex = rowIndex;
    }
  else
    {
      rowIndex = this->NValueGlobalShift + this->LaughlinIndex - 1 - MinQ;
      columnIndex = this->NValueGlobalShift - MinQ;
    }
}

// get the matrix that into account the Jordan Wigner string on the torus geometry
//
// nbrFermions = number of fermions in the system
// return value = corresponding matrix

SparseRealMatrix FQHEMPSLaughlinMatrix::GetTorusStringMatrix(int nbrFermions)
{
  if ((nbrFermions == 0) || ((nbrFermions & 1) == 1))
    {
      return this->AbstractFQHEMPSMatrix::GetTorusStringMatrix(nbrFermions);
    }
  int TmpDimension = 0;
  if (this->RealBMatrices != 0)
    {
      TmpDimension = this->RealBMatrices[0].GetNbrColumn();
    }
  else
    {
      TmpDimension = this->ComplexBMatrices[0].GetNbrColumn();
    }
  int* TmpNbrElementPerRow =  new int [TmpDimension];
  for (int i = 0; i < TmpDimension; ++i)
    {
      TmpNbrElementPerRow[i] = 1;
    }
  SparseRealMatrix StringMatrix;
  if (this->RealBMatrices != 0)
    StringMatrix = SparseRealMatrix(this->RealBMatrices[0].GetNbrRow(), TmpDimension, TmpNbrElementPerRow);
  else
    StringMatrix = SparseRealMatrix(this->ComplexBMatrices[0].GetNbrRow(), TmpDimension, TmpNbrElementPerRow);    
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      int MinQValue;
      int MaxQValue;
      this->ComputeGlobalChargeIndexRange(CurrentPLevel, MinQValue, MaxQValue);
      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	{
	  int TmpBondIndexRange = this->GetBondIndexRange(CurrentPLevel, CurrentQValue);
	  if ((CurrentQValue & 1) == 0)
	    {
	      for (int i = 0; i < TmpBondIndexRange; ++i)
		{
		  int TmpIndex = this->GetBondIndexWithFixedChargeAndPLevel(i, CurrentPLevel, CurrentQValue);	      
		  StringMatrix.SetMatrixElement(TmpIndex, TmpIndex, 1.0);
		}
	    }
	  else
	    {
	      for (int i = 0; i < TmpBondIndexRange; ++i)
		{
		  int TmpIndex = this->GetBondIndexWithFixedChargeAndPLevel(i, CurrentPLevel, CurrentQValue);	      
		  StringMatrix.SetMatrixElement(TmpIndex, TmpIndex, -1.0);
		}
	    }
	}
    }
  return StringMatrix;
}

// get the auxiliary space indices that are related to a given topological scetor
//
// topologicalSector = index of the topological sector to select
// nbrIndices = reference on the integer that will be set to the number of indices
// return value = array that contains the auxiliary space indices related to the selected topological sector

int* FQHEMPSLaughlinMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
{
  nbrIndices = 0;
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      int MinQValue;
      int MaxQValue;
      this->ComputeGlobalChargeIndexRange(CurrentPLevel, MinQValue, MaxQValue);
      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	{
	  if ((CurrentQValue % this->LaughlinIndex) == topologicalSector)
	    {
	      nbrIndices += this->GetBondIndexRange(CurrentPLevel, CurrentQValue);
	    }
	}
    }
  int* TmpIndices =  new int [nbrIndices];
  nbrIndices = 0;
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      int MinQValue;
      int MaxQValue;
      this->ComputeGlobalChargeIndexRange(CurrentPLevel, MinQValue, MaxQValue);
      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	{
	  if ((CurrentQValue % this->LaughlinIndex) == topologicalSector)
	    {
	      int TmpBondIndexRange = this->GetBondIndexRange(CurrentPLevel, CurrentQValue);
	      for (int i = 0; i < TmpBondIndexRange; ++i)
		{
		  TmpIndices[nbrIndices] = this->GetBondIndexWithFixedChargeAndPLevel(i, CurrentPLevel, CurrentQValue);
		  ++nbrIndices;
		}
	    }
	}
    }
  return TmpIndices;
}
  
// get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
// 
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true if we are dealing with fermions
// return value = minimum ky momentum 

int FQHEMPSLaughlinMatrix::GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics)
{
  if (statistics == false)
    return 0;
  if ((nbrParticles & 1) == 0)
    return (FindGCD(nbrParticles, nbrFluxQuanta) / 2);
  else
    return 0;
}
  
// print a given state of the auxiliary space
//
// str = reference on the output stream
// index = index of the state
// return value = reference on the output stream

ostream& FQHEMPSLaughlinMatrix::PrintAuxiliarySpaceState(ostream& str, int index)
{
  int TmpPLevel;
  int TmpQ;
  this->GetChargeAndPLevelFromMatrixIndex(index, TmpPLevel, TmpQ);
  str << "|" << index << ": Q=" << TmpQ << " P=" << TmpPLevel << ">";
  return str;
}

// get a given physical indiex
//
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void FQHEMPSLaughlinMatrix::GetPhysicalIndex(int index, unsigned long* configuration)
{  
  configuration[0] = this->PhysicalIndices[index];
}

// get the array where the site-dependent matrices are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseRealMatrix** FQHEMPSLaughlinMatrix::GetSiteDependentMatrices(int nbrFluxQuanta)
{
  cout << "FQHEMPSLaughlinMatrix::GetSiteDependentMatrices is not properly implemented" << endl;

  int TmpNbrOrbitals = nbrFluxQuanta;
  if (this->TorusFlag == false)
    {
      ++TmpNbrOrbitals;
    }

  SparseRealMatrix** BMatrices = new SparseRealMatrix*[TmpNbrOrbitals];
  for (int i = 0; i < TmpNbrOrbitals; ++i)
    {
      BMatrices[i] = new SparseRealMatrix[this->NbrBMatrices];
    }

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
    }
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValuesPerPLevel[0];
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValuesPerPLevel[i];
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << "x" << MatrixSize << endl;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;
  
  
  
  for (int OrbitalIndex = 0; OrbitalIndex < TmpNbrOrbitals; ++OrbitalIndex)
    {
      int* TmpNbrElementPerRow = new int[MatrixSize];
      for (int i = 0; i < MatrixSize; ++i)
	TmpNbrElementPerRow[i] = 0;
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      if ((j - i) == OrbitalIndex)
		{
		  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
		  int N2 = (j - i) + this->NValueGlobalShift;
		  int N1 = N2 + (this->LaughlinIndex);
		  if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
		      && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		    { 
		      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
			{
			  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
			    {
			      ++TmpNbrElementPerRow[this->GetMatrixIndex(i, k1, N1)];
			    }
			}
		    }
		}
	    }
	}
      SparseRealMatrix V0Matrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      if ((j - i) == OrbitalIndex)
		{
		  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
		  int N2 = (j - i) + this->NValueGlobalShift;
		  int N1 = N2 + (this->LaughlinIndex);
		  if (((N1 >= this->NInitialValuePerPLevel[i]) && (N1 <= this->NLastValuePerPLevel[i]))
		      && ((N2 >= this->NInitialValuePerPLevel[j]) && (N2 <= this->NLastValuePerPLevel[j])))
		    { 
		      for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
			{
			  TmpSpace1->GetOccupationNumber(k1, Partition1);
			  for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
			    {
			      TmpSpace2->GetOccupationNumber(k2, Partition2);
			      double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex, 1, Partition1, Partition2, i, j, Coef);
			      V0Matrix.SetMatrixElement(this->GetMatrixIndex(i, k1, N1), this->GetMatrixIndex(j, k2, N2), Tmp);
			    }
			}
		    }
		}
	    }
	}
      BMatrices[OrbitalIndex][1] = V0Matrix;
    }
  
  delete[] Partition1;
  delete[] Partition2;

//   int TmpNbrBMatrices = 0;
//   for (int m = 1; m < this->NbrBMatrices; ++m)
//     {
//       BMatrices[m] = MemoryEfficientMultiply(BMatrices[m - 1], V0Matrix);
//       if (BMatrices[m].GetNbrRow() > 0)
// 	{
// 	  ++TmpNbrBMatrices;
// 	  if ((this->CylinderFlag) || (this->TorusFlag))
// 	    BMatrices[m] /= sqrt((double) m);
// 	}
//     }

//   TmpNbrBMatrices = 0;

//   for (int i = 0; i < this->NbrBMatrices; ++i)
//     {
//       if (BMatrices[i].GetNbrRow() > 0)
// 	this->RealBMatrices[TmpNbrBMatrices++] = BMatrices[i];
//     }
//   delete[] BMatrices;

  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;

  return BMatrices;
}

// compute the site-dependent matrices
//
// initialOrbitalIndex = index of the first orbital
// lastOrbitalIndex = index of the last orbital

void FQHEMPSLaughlinMatrix::ComputeSiteDependentMatrices(int initialOrbitalIndex, int lastOrbitalIndex)
{
  if (this->TorusFlag == true)
    {
      cout << "FQHEMPSLaughlinMatrix::ComputeSiteDependentMatrices not implemented for the torus geometry" << endl;
      exit(0);
    }
   if (this->BosonicVersion == false)
     {
       cout << "FQHEMPSLaughlinMatrix::ComputeSiteDependentMatrices only works using the MPS bosonic convention" << endl;
       exit(0);      
     }
   if (this->SiteDependentMatrixNbrOrbitals != 0)
    {
      for (int i = 0; i < this->SiteDependentMatrixNbrOrbitals; ++i)
	{
	  delete[] this->SiteDependentMatrices[i];
	  delete[] this->SiteDependentPhysicalIndices[i];
	}
      delete[] this->SiteDependentMatrices;
      delete[] this->SiteDependentMatrixOrbitalIndices;
      delete[] this->NbrSiteDependentMatrices;
      delete[] this->SiteDependentPhysicalIndices;
    }


  int OrbitalZeroIndex = -1;
   int OrbitalMinusOneIndex = -1;
   if (initialOrbitalIndex < lastOrbitalIndex)
     {
       this->SiteDependentMatrixNbrOrbitals = lastOrbitalIndex - initialOrbitalIndex + 1;
       this->SiteDependentMatrixOrbitalIndices = new int [this->SiteDependentMatrixNbrOrbitals];
       for (int i = initialOrbitalIndex; i <= lastOrbitalIndex; ++i)
	 {
	   this->SiteDependentMatrixOrbitalIndices[i - initialOrbitalIndex] = i;
	   if (i == 0)
	     {
	       OrbitalZeroIndex = i - initialOrbitalIndex;
	     }
	   else
	     {
	       if (i == -1)
		 {
		   OrbitalMinusOneIndex = i - initialOrbitalIndex;
		 }	       
	     }
	 }
     }
   else
     {
       this->SiteDependentMatrixNbrOrbitals = initialOrbitalIndex - lastOrbitalIndex  + 1;
       this->SiteDependentMatrixOrbitalIndices = new int [this->SiteDependentMatrixNbrOrbitals];
       for (int i = initialOrbitalIndex; i >= lastOrbitalIndex; --i)
	 {
	   this->SiteDependentMatrixOrbitalIndices[i - lastOrbitalIndex] = i;
	   if (i == 0)
	     {
	       OrbitalZeroIndex = i - lastOrbitalIndex;
	     }
	   else
	     {
	       if (i == -1)
		 {
		   OrbitalMinusOneIndex = i - lastOrbitalIndex;
		 }	       
	     }
	 }
     }
   if (OrbitalZeroIndex == -1)
     {
       OrbitalZeroIndex = this->SiteDependentMatrixNbrOrbitals;
     }
   this->NbrSiteDependentMatrices = new int [this->SiteDependentMatrixNbrOrbitals];
   this->SiteDependentMatrices = new SparseRealMatrix* [this->SiteDependentMatrixNbrOrbitals];
   this->SiteDependentPhysicalIndices = new unsigned long* [this->SiteDependentMatrixNbrOrbitals];

   //  computing all site-dependent matrices for positive orbital indices
   SparseRealMatrix* TmpMatrices = new SparseRealMatrix [this->NbrBMatrices];
   int TmpNbrMatrices = 0;
   if ((OrbitalZeroIndex != this->SiteDependentMatrixNbrOrbitals) || (this->SiteDependentMatrixOrbitalIndices[0] > 0))
     {
       bool TmpFlag = false;
       TmpNbrMatrices = 0;
       for (int i = 0; (i < this->NbrBMatrices) && (TmpFlag == false); ++i)
	 {
	   TmpMatrices[i] = MemoryEfficientMultiply(this->RealBMatrices[i], this->InverseBMatrixZero);
	   if (TmpMatrices[i].GetNbrRow() > 0)
	     {
	       ++TmpNbrMatrices;
	     }
	   else
	     {
	       TmpFlag = true;
	     }
	 }
       if (OrbitalZeroIndex == this->SiteDependentMatrixNbrOrbitals)
	 {
	   int TmpIndex = 0;
	   while (TmpIndex < this->SiteDependentMatrixOrbitalIndices[0])
	     {
	       TmpFlag = false;
	       SparseRealMatrix* TmpMatrices2 = new SparseRealMatrix [TmpNbrMatrices];
	       TmpMatrices2[0] = TmpMatrices[0];
	       int TmpNbrMatrices2 = 1;
	       for (int i = 1; (i < TmpNbrMatrices) && (TmpFlag == false); ++i)
		 {
		   SparseRealMatrix TmpMatrix = MemoryEfficientMultiply(this->RealBMatrices[0], TmpMatrices[i]);
		   TmpMatrices2[i] = MemoryEfficientMultiply(TmpMatrix, this->InverseBMatrixZero);
		   if (TmpMatrices2[i].GetNbrRow() > 0)
		     {
		       ++TmpNbrMatrices2;
		     }
		   else
		     {
		       TmpFlag = true;
		     }
		 }
	       delete[] TmpMatrices;
	       TmpNbrMatrices = TmpNbrMatrices2;
	       TmpMatrices = TmpMatrices2;
	       ++TmpIndex;
	     }
	   OrbitalZeroIndex = 0;
	 }
       this->NbrSiteDependentMatrices[OrbitalZeroIndex] = TmpNbrMatrices;
       this->SiteDependentMatrices[OrbitalZeroIndex] = new SparseRealMatrix[TmpNbrMatrices];
       for (int i = 0; i < TmpNbrMatrices; ++i)
	 {	   
	   this->SiteDependentMatrices[OrbitalZeroIndex][i] = TmpMatrices[i];
	 }
       OrbitalZeroIndex++;
     }
   for (; OrbitalZeroIndex < this->SiteDependentMatrixNbrOrbitals; ++OrbitalZeroIndex)
     {
       bool TmpFlag = false;
       TmpMatrices[0] = this->SiteDependentMatrices[OrbitalZeroIndex - 1][0];
       TmpNbrMatrices = 1;
       for (int i = 1; (i < this->NbrSiteDependentMatrices[OrbitalZeroIndex - 1]) && (TmpFlag == false); ++i)
	 {
	   SparseRealMatrix TmpMatrix = MemoryEfficientMultiply(this->RealBMatrices[0], this->SiteDependentMatrices[OrbitalZeroIndex - 1][i]);
	   TmpMatrices[i] = MemoryEfficientMultiply(TmpMatrix, this->InverseBMatrixZero);
	   if (TmpMatrices[i].GetNbrRow() > 0)
	     {
	       ++TmpNbrMatrices;
	     }
	   else
	     {
	       TmpFlag = true;
	     }
	 }
       this->NbrSiteDependentMatrices[OrbitalZeroIndex] = TmpNbrMatrices;
       this->SiteDependentMatrices[OrbitalZeroIndex] = new SparseRealMatrix[this->NbrSiteDependentMatrices[OrbitalZeroIndex]];
       for (int i = 0; i < this->NbrSiteDependentMatrices[OrbitalZeroIndex]; ++i)
	 {
	   this->SiteDependentMatrices[OrbitalZeroIndex][i] = TmpMatrices[i];
	 }
     }

   //  computing all site-dependent matrices for negative orbital indices
   if ((OrbitalMinusOneIndex != -1) || (this->SiteDependentMatrixOrbitalIndices[this->SiteDependentMatrixNbrOrbitals - 1] < 0))
     {
       bool TmpFlag = false;
       TmpNbrMatrices = 0;
       for (int i = 0; i < (this->NbrBMatrices) && (TmpFlag == false); ++i)
	 {
	   TmpMatrices[i] = MemoryEfficientMultiply(this->InverseBMatrixZero, this->RealBMatrices[i]);
	   if (TmpMatrices[i].GetNbrRow() > 0)
	     {
	       ++TmpNbrMatrices;
	     }
	   else
	     {
	       TmpFlag = true;
	     }
	 }
       if (OrbitalMinusOneIndex == -1)
	 {
	   TmpFlag = false;
	   int TmpIndex = -1;
	   while (TmpIndex > this->SiteDependentMatrixOrbitalIndices[this->SiteDependentMatrixNbrOrbitals - 1])
	     {
	       SparseRealMatrix* TmpMatrices2 = new SparseRealMatrix [TmpNbrMatrices];
	       TmpMatrices2[0] = TmpMatrices[0];
	       int TmpNbrMatrices2 = 1;	   
	       for (int i = 1; (i < TmpNbrMatrices) && (TmpFlag == false); ++i)
		 {
		   SparseRealMatrix TmpMatrix = MemoryEfficientMultiply(this->InverseBMatrixZero, TmpMatrices[i]);
		   TmpMatrices2[i] = MemoryEfficientMultiply(TmpMatrix, this->RealBMatrices[0]);
		   if (TmpMatrices2[i].GetNbrRow() > 0)
		     {
		       ++TmpNbrMatrices2;
		     }
		   else
		     {
		       TmpFlag = true;
		     }
		 }
	       delete[] TmpMatrices;
	       TmpNbrMatrices = TmpNbrMatrices2;
	       TmpMatrices = TmpMatrices2;
	       --TmpIndex;
	     }
	   OrbitalMinusOneIndex = this->SiteDependentMatrixNbrOrbitals - 1;
	 }
       this->NbrSiteDependentMatrices[OrbitalMinusOneIndex] = TmpNbrMatrices;
       this->SiteDependentMatrices[OrbitalMinusOneIndex] = new SparseRealMatrix[TmpNbrMatrices];
       for (int i = 0; i < TmpNbrMatrices; ++i)
	 {	   
	   this->SiteDependentMatrices[OrbitalMinusOneIndex][i] = TmpMatrices[i];
	 }
       OrbitalMinusOneIndex--;
     }
   for (; OrbitalMinusOneIndex >= 0; --OrbitalMinusOneIndex)
     {
       bool TmpFlag = false;
       TmpMatrices[0] = this->SiteDependentMatrices[OrbitalMinusOneIndex + 1][0];
       TmpNbrMatrices = 1;
       for (int i = 1; (i < this->NbrSiteDependentMatrices[OrbitalMinusOneIndex + 1]) && (TmpFlag == false); ++i)
	 {
	   SparseRealMatrix TmpMatrix = MemoryEfficientMultiply(this->InverseBMatrixZero, this->SiteDependentMatrices[OrbitalMinusOneIndex + 1][i]);
	   TmpMatrices[i] = MemoryEfficientMultiply(TmpMatrix, this->RealBMatrices[0]);
	   if (TmpMatrices[i].GetNbrRow() > 0)
	     {
	       ++TmpNbrMatrices;
	     }
	   else
	     {
	       TmpFlag = true;
	     }
	 }
       this->NbrSiteDependentMatrices[OrbitalMinusOneIndex] = TmpNbrMatrices;
       this->SiteDependentMatrices[OrbitalMinusOneIndex] = new SparseRealMatrix[this->NbrSiteDependentMatrices[OrbitalMinusOneIndex]];
       for (int i = 0; i < this->NbrSiteDependentMatrices[OrbitalMinusOneIndex]; ++i)
	 {
	   this->SiteDependentMatrices[OrbitalMinusOneIndex][i] = TmpMatrices[i];
	 }
     }

   for (int i = 0; i < this->SiteDependentMatrixNbrOrbitals; ++i)
     {
       this->SiteDependentPhysicalIndices[i] = new unsigned long[this->NbrSiteDependentMatrices[i]];
       for (int j = 0; j < this->NbrSiteDependentMatrices[i]; ++j)
	 {
	   this->SiteDependentPhysicalIndices[i][j] = (unsigned long) j;
	 }
     }
   
   delete[] TmpMatrices;
}

// get the site-dependent matrices (real version) computed through ComputeSiteDependentMatrices
//
// siteDependentMatrices = reference on the site-dependent matrices
// nbrSiteDependentMatrices = reference on the array providing the number of site-dependent matrices per orbital
// siteDependentMatrixOrbitalIndices = reference on the array providing the orbital indices 
// siteDependentPhysicalIndices = reference on the array providing the physical indices associated to each site-dependent matrix
// return value = number of orbitals covered by the site-dependent matrices

int FQHEMPSLaughlinMatrix::GetSiteDependentMatrices(SparseRealMatrix**& siteDependentMatrices, int*& nbrSiteDependentMatrices, int*& siteDependentMatrixOrbitalIndices, unsigned long**& siteDependentPhysicalIndices)
{
  siteDependentMatrices = this->SiteDependentMatrices;
  nbrSiteDependentMatrices = this->NbrSiteDependentMatrices;
  siteDependentMatrixOrbitalIndices = this->SiteDependentMatrixOrbitalIndices;
  siteDependentPhysicalIndices = this->SiteDependentPhysicalIndices;
  return this->SiteDependentMatrixNbrOrbitals;
}
  
