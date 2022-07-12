////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of MPS matrix for the PH-Pfaffian state               //
//                                                                            //
//                        last modification : 19/05/2017                      //
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
#include "Tools/FQHEMPS/FQHEMPSPHPfaffianMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSPHPfaffianMatrix::FQHEMPSPHPfaffianMatrix()
{
  this->CFTDirectory = 0;
  this->Architecture = 0;
  this->UseRationalFlag = true;
}

// constructor 
//
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction
// architecture = architecture to use for precalculation

FQHEMPSPHPfaffianMatrix::FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, const char* cftDirectory, bool bosonicVersion, bool useRational, 
						 bool trimChargeIndices, bool cylinderFlag, double kappa, 
						 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
						 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RIndex = 2;
  //  this->BosonicVersion = bosonicVersion;
  this->BosonicVersion = true;
  this->LaughlinIndex = laughlinIndex;
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  this->PLevel = 2 * pLevel;
  this->PLevelShift = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevelShift;
  this->CylinderFlag = cylinderFlag;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->Kappa = kappa;
  this->TorusFlag = torusFlag;
  if (this->TorusFlag == true)
    {
      this->TorusNbrFluxQuanta = nbrFluxQuanta;
      this->TorusAngle = angle;
      this->TorusAspectRatio = aspectRatio;
      this->TorusFluxInsertion = fluxInsertion;
      this->Kappa = sqrt(2.0 * M_PI * this->TorusAspectRatio / ((double) this->TorusNbrFluxQuanta));
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
    }
  else
    {
      this->TwistedTorusFlag = false;
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      this->ComplexBMatrices = 0;
    }
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "phpfaffian");
  this->CreateBMatrices(cftDirectory, architecture);
  if (cftDirectory != 0)
    {
      this->CFTDirectory = new char [strlen(cftDirectory) + 1];
      strcpy (this->CFTDirectory, cftDirectory);
    }
  else
    {
      this->CFTDirectory = 0;
    }
  this->Architecture = architecture;
}


// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction

FQHEMPSPHPfaffianMatrix::FQHEMPSPHPfaffianMatrix(int laughlinIndex, int pLevel, const char* fileName, 
						 bool trimChargeIndices, bool cylinderFlag, double kappa, 
						 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion)
{
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = 2 * pLevel;
  this->PLevelShift = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevelShift;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->TorusFlag = torusFlag;
  if (this->TorusFlag == true)
    {
      this->TorusNbrFluxQuanta = nbrFluxQuanta;
      this->TorusAngle = angle;
      this->TorusAspectRatio = aspectRatio;
      this->TorusFluxInsertion = fluxInsertion;
      this->Kappa = sqrt(2.0 * M_PI * this->TorusAspectRatio / ((double) this->TorusNbrFluxQuanta));
      if ((this->TorusAngle != 0.0) || (this->TorusFluxInsertion != 0.0))
	{
	  this->TwistedTorusFlag = true;
	  this->TauFactor = ((2.0 * M_PI * this->TorusAspectRatio  / ((double) this->TorusNbrFluxQuanta))
			     * Complex(-sin(this->TorusAngle * M_PI), cos(this->TorusAngle * M_PI)));
	}
      else
	{
	  this->TwistedTorusFlag = false;
	}
    }
  else
    {
      this->TwistedTorusFlag = false;
    }
  this->LoadMatrices(fileName);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "phpfaffian");
}

// destructor
//

FQHEMPSPHPfaffianMatrix::~FQHEMPSPHPfaffianMatrix()
{
}
  
// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int FQHEMPSPHPfaffianMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  return ((nbrFluxQuanta + 1) / 2);
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSPHPfaffianMatrix::CreateBMatrices (const char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;

  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();

  double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevelShift + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevelShift + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevelShift + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevelShift + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevelShift + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevelShift + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevelShift + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevelShift + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevelShift + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevelShift + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevelShift + 1];

  LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevelShift + 1];
  double** MultiplicityFactor = new double*[this->PLevelShift + 1];
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevelShift + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevelShift + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevelShift + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevelShift + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevelShift + 1];
      RationalMultiplicityFactor[i] = new LongRational[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      MultiplicityFactor[i] = new double[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      for (int j = 0; j < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
	{
	  U1BosonBasis[i]->GetOccupationNumber(j, TmpPartition);	    
	  RationalMultiplicityFactor[i][j] = 1l;
	  MultiplicityFactor[i][j] = 1.0;
 	  for (int k = 1; k <= i; ++k)
 	    if (TmpPartition[k] > 1ul)
	      {
		RationalMultiplicityFactor[i][j].FactorialDivide(TmpPartition[k]);
		double Tmp = 1.0;
		for (unsigned long l = 2l; l <= TmpPartition[k]; ++l)
		  Tmp *=  (double) l;
		MultiplicityFactor[i][j] /= Tmp;
	      }
	}
    }

  cout << "weight: " <<   this->WeightIdentity << " " << this->WeightPsi << endl;
  char* TmpScalarProductIdentityFileName = 0; 
  char* TmpScalarProductPsiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	}
      this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductIdentityFileName, architecture, RationalScalarProductIdentity, ScalarProductIdentity, i, U1BosonBasis,
					   CentralCharge12, CentralCharge12Numerical, this->WeightIdentity, WeightIdentityNumerical, "identity",
					   OrthogonalBasisIdentityLeft, OrthogonalBasisIdentityRight, RationalMultiplicityFactor, MultiplicityFactor);
      this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductPsiFileName, architecture, RationalScalarProductPsi, ScalarProductPsi, i, U1BosonBasis,
					   CentralCharge12, CentralCharge12Numerical, this->WeightPsi, WeightPsiNumerical, "psi",
					   OrthogonalBasisPsiLeft, OrthogonalBasisPsiRight, RationalMultiplicityFactor, MultiplicityFactor);
      cout << "---------------------------------" << endl;
    }
  this->RescaleFullScalarProductMatrix(RationalScalarProductIdentity, ScalarProductIdentity, RationalMultiplicityFactor, MultiplicityFactor);
  this->RescaleFullScalarProductMatrix(RationalScalarProductPsi, ScalarProductPsi, RationalMultiplicityFactor, MultiplicityFactor);

  this->U1BasisDimension = new int [this->PLevelShift + 1];	
  this->NeutralSectorDimension = new int* [2];
  this->NeutralSectorDimension[0] = new int [this->PLevelShift + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevelShift + 1];
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }

  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  double ExtraCylinderFactor = 1.0;
  if (this->BosonicVersion == false)
    {
      QValue = 1 + (this->RIndex / 2);
      this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
      NValueShift = 2 * this->PLevel - 1;
      QValueDenominator = 1;
    }
  else
    {
      QValue = this->LaughlinIndex - 1 + (this->RIndex / 2);
      this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 2;
      NValueShift = 2 * this->PLevel - 1 - this->RIndex / 2;
      QValueDenominator = 1;
    }
  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  for (int j = 0; j <= this->PLevelShift; ++j)
    {
      for (int i = 0; i <= this->PLevelShift; ++i)
	{
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	      else
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	    }
	  this->ComputeFullMatrixElements(cftDirectory, TmpScalarProductIdentityFileName, architecture, 
					  RationalMatrixPsi01, MatrixPsi01, i, j, U1BosonBasis, 
					  CentralCharge12, CentralCharge12Numerical, 
					  this->WeightIdentity, WeightIdentityNumerical, 
					  this->WeightPsi, WeightPsiNumerical,
					  this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	  this->ComputeFullMatrixElements(cftDirectory, TmpScalarProductPsiFileName, architecture, 
					  RationalMatrixPsi10, MatrixPsi10, i, j, U1BosonBasis, 
					  CentralCharge12, CentralCharge12Numerical, 
					  this->WeightPsi, WeightPsiNumerical, 
					  this->WeightIdentity, WeightIdentityNumerical,
					  this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	}
    }
  this->RescaleFullMatrixElements(RationalMatrixPsi01, MatrixPsi01, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);
  this->RescaleFullMatrixElements(RationalMatrixPsi10, MatrixPsi10, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);

  cout << "building B matrices" << endl;

  SparseRealMatrix* BMatrices = 0;
  SparseComplexMatrix* TmpComplexBMatrices = 0;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;

  // B^[0]  matrix evaluation
  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel)  && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[n];
	      RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[n];
	      RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[n];
	      RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[n];
	      RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[n];
	      for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
		{	      
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}
	    }
	}
    }

  if (this->TwistedTorusFlag == false)
    {
      BMatrices = new SparseRealMatrix[this->NbrBMatrices];
      BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
    }
  else
    {
      TmpComplexBMatrices = new SparseComplexMatrix[this->NbrBMatrices];
      TmpComplexBMatrices[0] = SparseComplexMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
    }
  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel)  && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[n];
	      RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[n];
	      RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[n];
	      RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[n];
	      RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[n];
	      for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
		{	      
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) (p + n))
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) (p + n))
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				    }
				  
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j - 1, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightIdentityNumerical +  ((double) (p + n))
									      + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j - 1, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) (p + n))
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) (p + n))
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j - 1, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightPsiNumerical +  ((double) (p + n))
									      + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j - 1, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  // B^[1]  matrix evaluation
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  
  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel) && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
	      for (int m = 0; m <= this->PLevelShift; ++m)
		{
		  for (int q = 0; q <= this->PLevelShift; ++q)
		    {
		      if (((this->PLevelShift + q - m) <= this->PLevel) && ((this->PLevelShift + q - m) >= 0))
			{
			  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
			  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[m];
			  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[m];
			  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[m];
			  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[n][m];
			  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[n][m];	
			  
			  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			    {	      
			      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
				{	      
				  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
				  int N2;
				  int N1;
				  N2 = (2 * ((q - m) - (p - n)) - this->RIndex + 1 + NValueShift) / 2;
				  N1 = N2 + QValue;
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][1])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					    }
					}
				      
				    }
				  N2 = (2 * ((q - m) - (p - n)) + 1 + NValueShift) / 2;
				  N1 = N2 + QValue;
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, N1, p, ChargedIndex1, NeutralIndex1)];
					    }
					}
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  SparseRealMatrix V0Matrix (MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel) && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
	      for (int m = 0; m <= this->PLevelShift; ++m)
		{
		  for (int q = 0; q <= this->PLevelShift; ++q)
		    {
		      if (((this->PLevelShift + q - m) <= this->PLevel) && ((this->PLevelShift + q - m) >= 0))
			{
			  double TmpProjectionFactor = 1.0;
			  if (this->CylinderFlag)
			    {
			      TmpProjectionFactor = exp (this->Kappa * this->Kappa * ((double) ((n - m) * (n - m))));
			    }
			  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
			  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[m];
			  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[m];
			  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[m];
			  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[n][m];
			  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[n][m];	
			  
			  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			    {	      
			      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
				{	      
				  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
				  int N2;
				  int N1;
				  N2 = (2 * ((q - m) - (p - n)) - this->RIndex + 1 + NValueShift) / 2;
				  N1 = N2 + QValue;
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][1])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      double Tmp = 0.0;
					      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
						{
						  double Tmp1 = 0.0;			      
						  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						    {
						      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisIdentity1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, p, q, Coef);
					      Tmp *= TmpProjectionFactor;
					      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(this->PLevelShift + q - m, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    }
					  
					}
				    }
				  N2 = (2 * ((q - m) - (p - n)) + 1 + NValueShift) / 2;
				  N1 = N2 + QValue;
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      double Tmp = 0.0;
					      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
						{
						  double Tmp1 = 0.0;			      
						  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						    {
						      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentity2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, p, q, Coef);
					      Tmp *= TmpProjectionFactor;
					      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(this->PLevelShift + q - m, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    }
					}
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}

    }

  if (this->TwistedTorusFlag == false)
    {
      for (int m = 1; m < this->NbrBMatrices; ++m)
	{
	  BMatrices[m] = MemoryEfficientMultiply(BMatrices[m - 1], V0Matrix);
	  if (BMatrices[m].GetNbrRow() > 0)
	    {
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
	      TmpComplexBMatrices[m] /= sqrt((double) m);
	    }
	}
    }

  int TmpNbrBMatrices = 0;
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
  
  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete U1BosonBasis[i];
      delete[] RationalMultiplicityFactor[i];
      delete[] MultiplicityFactor[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
  delete[] RationalMultiplicityFactor;
  delete[] MultiplicityFactor;
}

// compute the various arrays required to convert from quantum numbers and local indices to a global linearized index
//
// return value = dimension of the B matrix

int FQHEMPSPHPfaffianMatrix::ComputeLinearizedIndexArrays()
{
  this->NbrNValuesPerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NInitialValuePerPLevelCFTSector = new int* [this->PLevel + 1];
  this->NLastValuePerPLevelCFTSector = new int* [this->PLevel + 1];     
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NbrNValuesPerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      this->NInitialValuePerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      this->NLastValuePerPLevelCFTSector[i] = new int [this->NbrCFTSectors];
      for (int l = 0; l < this->NbrCFTSectors; ++l)
	{
	  this->ComputeChargeIndexRange(i, l, this->NInitialValuePerPLevelCFTSector[i][l], this->NLastValuePerPLevelCFTSector[i][l]);
	  this->NbrNValuesPerPLevelCFTSector[i][l] =  this->NLastValuePerPLevelCFTSector[i][l] - this->NInitialValuePerPLevelCFTSector[i][l] + 1;
	}
    }

  this->StartingIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  this->NbrIndexPerPLevelCFTSectorQValue = new int** [this->PLevel + 1];
  this->StartingIndexPerPLevelCFTSectorQValueU1Sector = new int*** [this->PLevel + 1];
  this->NbrIndexPerPLevelCFTSectorQValueU1Sector = new int*** [this->PLevel + 1];
  int TotalIndex = 0;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->StartingIndexPerPLevelCFTSectorQValue[i] = new int* [this->NbrCFTSectors];
      this->NbrIndexPerPLevelCFTSectorQValue[i] = new int* [this->NbrCFTSectors];
      this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i] = new int** [this->NbrCFTSectors];
      this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i] = new int** [this->NbrCFTSectors];
      for (int l = 0; l < this->NbrCFTSectors; ++l)
	{
	  this->StartingIndexPerPLevelCFTSectorQValue[i][l] = new int [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->NbrIndexPerPLevelCFTSectorQValue[i][l] = new int [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l] = new int* [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l] = new int* [this->NbrNValuesPerPLevelCFTSector[i][l]];
	  for (int j = this->NInitialValuePerPLevelCFTSector[i][l];  j <= this->NLastValuePerPLevelCFTSector[i][l]; ++j)
	    {
	      this->StartingIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = TotalIndex;
	      this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = new int [i + 1];
	      this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = new int [i + 1];
	      int TmpMinChargePSector = i - this->PLevelShift;
	      if (TmpMinChargePSector < 0)
		{
		  TmpMinChargePSector = 0;
		}
	      int TmpMaxChargePSector = i;
	      if (TmpMaxChargePSector > this->PLevelShift)
		{
		  TmpMaxChargePSector = this->PLevelShift;
		}
	      for (int k = TmpMinChargePSector; k <= TmpMaxChargePSector; ++k)
		{
		  this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = TotalIndex;
		  int Tmp = this->U1BasisDimension[k] * this->NeutralSectorDimension[l][k + this->PLevelShift - i];
		  this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = Tmp;
		  TotalIndex += Tmp;	      
		}
	      this->NbrIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = TotalIndex - this->StartingIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]];
	    }
	}
    }
  return TotalIndex;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSPHPfaffianMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(this->PLevelShift, MinQ, MaxQ);
  if (padding == true)
    {
      if (this->BosonicVersion == false)
	{
	  rowIndex = this->StartingIndexPerPLevelCFTSectorQValue[this->PLevelShift][0][(2 * this->PLevelShift + 2) - MinQ];
//	  rowIndex = (2 * this->PLevelShift + 2) - MinQ;
	  columnIndex = rowIndex;
	}
      else
	{
	  rowIndex = this->StartingIndexPerPLevelCFTSectorQValue[this->PLevelShift][0][(2 * this->PLevelShift + 2) - MinQ];
	  columnIndex = rowIndex;
	}
    }
  else
    {
      if (this->BosonicVersion == false)
	{
	  rowIndex = this->PLevel + this->RIndex - MinQ;
	  columnIndex = this->PLevel - MinQ;
	}
      else
	{
	  rowIndex = this->PLevel + this->LaughlinIndex - MinQ;
	  columnIndex = this->PLevel - MinQ;
	}
    }
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSPHPfaffianMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
{
  if (this->UniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }
  minQ = 0;
  maxQ = this->NbrNValue - 1;
  int TmpMinQ = this->NbrNValue - 1;
  int TmpMaxQ = 0;    

  if (this->BosonicVersion == false)
    {
      if ((this->RIndex & 1) == 0)
	{
	  int NValueShift = this->PLevel;
	  int QValue = 1 + (this->RIndex / 2);
	  if (cftSector == 0) 
	    {
	      for (int Q = 0; Q < this->NbrNValue; ++Q)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 1);
		      TmpP += QPrime - (this->RIndex / 2) - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - 1);
			  TmpP += QPrime - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= QPrime - NValueShift;
		      QPrime += (QValue - 1);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= QPrime  - (this->RIndex / 2) - NValueShift;
			  QPrime += (QValue - 1);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	
		    }    
		}
	    }
	  else
	    {
	      for (int Q = 0; Q < this->NbrNValue; ++Q)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP3 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP3)
			TmpMaxP3 = TmpP;	    
		      QPrime -= (QValue - 1);
		      TmpP += QPrime - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP3)
			    TmpMaxP3 = TmpP;	    
			  QPrime -= (QValue - 1);
			  TmpP += QPrime - (this->RIndex / 2) - NValueShift;
			}
		    }
		  
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP4 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP4)
			TmpMaxP4 = TmpP;	    
		      TmpP -= QPrime - NValueShift - (this->RIndex / 2);
		      QPrime += (QValue - 1);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP4)
			    TmpMaxP4 = TmpP;	    
			  TmpP -= QPrime - NValueShift;
			  QPrime += (QValue - 1);
			}
		    }
		  if (((this->PLevel - TmpMaxP3) >= pLevel) && ((this->PLevel - TmpMaxP4) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	    
		    }
		}
	    }
	}
      else
	{
	  int NValueShift = this->PLevel;
	  int QValue = this->RIndex + 2;
	  if (cftSector == 0)
	    {
	      for (int Q = 0; Q < this->NbrNValue; Q += 2)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 2);
		      TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - 2);
			  TmpP += QPrime / 2 - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime / 2) - NValueShift;
		      QPrime += (QValue - 2);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
			  QPrime += (QValue - 2);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	    
		    }
		}
	    }
	  else
	    {
	      for (int Q = 1; Q < this->NbrNValue; Q += 2)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - 2);
		      TmpP += QPrime / 2 - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - 2);
			  TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
		      QPrime += (QValue - 2);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime / 2) - NValueShift;
			  QPrime += (QValue - 2);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	    
		    }
		}
	    }
	}
    }
  else
    {
      if ((this->RIndex & 1) == 0)
	{
	  int NValueShift = this->PLevel;
	  int QValue = 1 + (this->RIndex / 2);
	  if (cftSector == 0) 
	    {
	      for (int Q = 0; Q < this->NbrNValue; ++Q)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += QPrime - (this->RIndex / 2) - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += QPrime - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= QPrime - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= QPrime  - (this->RIndex / 2) - NValueShift;
			  QPrime += (QValue);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	
		    }    
		}
	    }
	  else
	    {
	      for (int Q = 0; Q < this->NbrNValue; ++Q)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP3 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP3)
			TmpMaxP3 = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += QPrime - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP3)
			    TmpMaxP3 = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += QPrime - (this->RIndex / 2) - NValueShift;
			}
		    }
		  
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP4 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP4)
			TmpMaxP4 = TmpP;	    
		      TmpP -= QPrime - NValueShift - (this->RIndex / 2);
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP4)
			    TmpMaxP4 = TmpP;	    
			  TmpP -= QPrime - NValueShift;
			  QPrime += (QValue);
			}
		    }
		  if (((this->PLevel - TmpMaxP3) >= pLevel) && ((this->PLevel - TmpMaxP4) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	    
		    }
		}
	    }
	}
      else
	{
	  int NValueShift = this->PLevel;
	  int QValue = this->RIndex + 2;
	  if (cftSector == 0)
	    {
	      for (int Q = 0; Q < this->NbrNValue; Q += 2)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += QPrime / 2 - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime / 2) - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
			  QPrime += (QValue);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	    
		    }
		}
	    }
	  else
	    {
	      for (int Q = 1; Q < this->NbrNValue; Q += 2)
		{
		  int QPrime = Q;
		  int TmpP = 0;
		  int TmpMaxP = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += QPrime / 2 - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += (QPrime - this->RIndex) / 2 - NValueShift;
			}
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - this->RIndex) / 2 - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime / 2) - NValueShift;
			  QPrime += (QValue);
			}
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
		    TmpMaxQ = Q;	    
		    }
		}
	    }
	}
    }
  minQ = TmpMinQ;
  maxQ = TmpMaxQ;
  cout << "range at p=" << pLevel << ", x=" << cftSector << " : " << minQ << " " << maxQ << " (" << this->NbrNValue << ")" << endl;   
}

// get the array where the site-dependent matrices for the geometry are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseRealMatrix** FQHEMPSPHPfaffianMatrix::GetSphereSiteDependentMatrices(int nbrFluxQuanta)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;

  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();

  double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevelShift + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevelShift + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevelShift + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevelShift + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevelShift + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevelShift + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevelShift + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevelShift + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevelShift + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevelShift + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevelShift + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevelShift + 1];

  LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevelShift + 1];
  double** MultiplicityFactor = new double*[this->PLevelShift + 1];
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevelShift + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevelShift + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevelShift + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevelShift + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevelShift + 1];
      RationalMultiplicityFactor[i] = new LongRational[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      MultiplicityFactor[i] = new double[U1BosonBasis[i]->GetHilbertSpaceDimension()];
      for (int j = 0; j < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
	{
	  U1BosonBasis[i]->GetOccupationNumber(j, TmpPartition);	    
	  RationalMultiplicityFactor[i][j] = 1l;
	  MultiplicityFactor[i][j] = 1.0;
 	  for (int k = 1; k <= i; ++k)
 	    if (TmpPartition[k] > 1ul)
	      {
		RationalMultiplicityFactor[i][j].FactorialDivide(TmpPartition[k]);
		double Tmp = 1.0;
		for (unsigned long l = 2l; l <= TmpPartition[k]; ++l)
		  Tmp *=  (double) l;
		MultiplicityFactor[i][j] /= Tmp;
	      }
	}
    }

  cout << "weight: " <<   this->WeightIdentity << " " << this->WeightPsi << endl;
  char* TmpScalarProductIdentityFileName = 0; 
  char* TmpScalarProductPsiFileName = 0;
  if (this->CFTDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(this->CFTDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(this->CFTDirectory)];
    }
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (this->CFTDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_scalarproducts_psi_level_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_scalarproducts_identity_level_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_scalarproducts_psi_level_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i);
	    }
	}
      this->ComputeFullScalarProductMatrix(this->CFTDirectory, TmpScalarProductIdentityFileName, this->Architecture, RationalScalarProductIdentity, ScalarProductIdentity, i, U1BosonBasis,
					   CentralCharge12, CentralCharge12Numerical, this->WeightIdentity, WeightIdentityNumerical, "identity",
					   OrthogonalBasisIdentityLeft, OrthogonalBasisIdentityRight, RationalMultiplicityFactor, MultiplicityFactor);
      this->ComputeFullScalarProductMatrix(this->CFTDirectory, TmpScalarProductPsiFileName, this->Architecture, RationalScalarProductPsi, ScalarProductPsi, i, U1BosonBasis,
					   CentralCharge12, CentralCharge12Numerical, this->WeightPsi, WeightPsiNumerical, "psi",
					   OrthogonalBasisPsiLeft, OrthogonalBasisPsiRight, RationalMultiplicityFactor, MultiplicityFactor);
      cout << "---------------------------------" << endl;
    }
  this->RescaleFullScalarProductMatrix(RationalScalarProductIdentity, ScalarProductIdentity, RationalMultiplicityFactor, MultiplicityFactor);
  this->RescaleFullScalarProductMatrix(RationalScalarProductPsi, ScalarProductPsi, RationalMultiplicityFactor, MultiplicityFactor);

  this->U1BasisDimension = new int [this->PLevelShift + 1];	
  this->NeutralSectorDimension = new int* [2];
  this->NeutralSectorDimension[0] = new int [this->PLevelShift + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevelShift + 1];
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }

  int NValueShift;
  int QValue;
  int QValueDenominator;
  double ExtraCylinderFactor = 1.0;

  QValue = 2;
  this->NbrNValue = 4 * this->PLevelShift + 5;
  NValueShift = 4 * this->PLevelShift + 4;
  QValueDenominator = 1;

  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  for (int j = 0; j <= this->PLevelShift; ++j)
    {
      for (int i = 0; i <= this->PLevelShift; ++i)
	{
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (this->CFTDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_matrixelement_identitypsi_level_%d_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_matrixelement_psiidentity_level_%d_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i, j);
		}
	      else
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_matrixelement_identitypsi_level_%d_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_matrixelement_psiidentity_level_%d_%d.dat", this->CFTDirectory, this->BMatrixOutputName, i, j);
		}
	    }
	  this->ComputeFullMatrixElements(this->CFTDirectory, TmpScalarProductIdentityFileName, this->Architecture, 
					  RationalMatrixPsi01, MatrixPsi01, i, j, U1BosonBasis, 
					  CentralCharge12, CentralCharge12Numerical, 
					  this->WeightIdentity, WeightIdentityNumerical, 
					  this->WeightPsi, WeightPsiNumerical,
					  this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	  this->ComputeFullMatrixElements(this->CFTDirectory, TmpScalarProductPsiFileName, this->Architecture, 
					  RationalMatrixPsi10, MatrixPsi10, i, j, U1BosonBasis, 
					  CentralCharge12, CentralCharge12Numerical, 
					  this->WeightPsi, WeightPsiNumerical, 
					  this->WeightIdentity, WeightIdentityNumerical,
					  this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	}
    }
  this->RescaleFullMatrixElements(RationalMatrixPsi01, MatrixPsi01, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);
  this->RescaleFullMatrixElements(RationalMatrixPsi10, MatrixPsi10, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);

  cout << "building B matrices" << endl;

  char** TmpLabels = new char* [MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    {
      TmpLabels[i] = this->GetAuxiliarySpaceLabel(i);
    }

  SparseRealMatrix TmpB0Matrix;
  SparseRealMatrix* BMatrices = 0;
  SparseComplexMatrix* TmpComplexBMatrices = 0;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;

  // B^[0]  matrix evaluation
  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel)  && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[n];
	      RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[n];
	      RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[n];
	      RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[n];
	      RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[n];
	      for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
		{	      
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}
	    }
	}
    }

  BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  TmpB0Matrix = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);

  for (int n = 0; n <= this->PLevelShift; ++n)
    {
      for (int p = 0; p <= this->PLevelShift; ++p)
	{
	  if (((this->PLevelShift + p - n) <= this->PLevel)  && ((this->PLevelShift + p - n) >= 0))
	    {
	      BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[n];
	      RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[n];
	      RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[n];
	      RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[n];
	      RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[n];
	      RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[n];
	      RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[n];
	      for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
		{	      
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      TmpB0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j - 1, p, ChargedIndex, NeutralIndex1),
							    this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      TmpB0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j - 1, p, ChargedIndex, NeutralIndex1),
							    this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		}
	    }
	}
    }
  
  // B^[1]  matrix evaluation
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];

  for (int i = 0; i < MatrixSize; ++i)
    {
      this->PrintAuxiliarySpaceState(cout, i) << endl;
    }

  int NbrV0MatrixIndices = (2 * this->PLevelShift) + 2;
  int V0MatrixIndexShift = this->PLevelShift + 1;
  SparseRealMatrix* V0Matrices = new SparseRealMatrix[NbrV0MatrixIndices];

  for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
    {
      for (int i = 0; i < MatrixSize; ++i)
	TmpNbrElementPerRow[i] = 0;
      long TmpTotalNbrElements = 0l;
      for (int n = 0; n <= this->PLevelShift; ++n)
	{
	  for (int p = 0; p <= this->PLevelShift; ++p)
	    {
	      if (((this->PLevelShift + p - n) <= this->PLevel) && ((this->PLevelShift + p - n) >= 0))
		{
		  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
		  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
		  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
		  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
		  for (int m = 0; m <= this->PLevelShift; ++m)
		    {
		      for (int q = 0; q <= this->PLevelShift; ++q)
			{
			  if (((this->PLevelShift + q - m) <= this->PLevel) && ((this->PLevelShift + q - m) >= 0))
			    {
			      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
			      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[m];
			      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[m];
			      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[m];
			      RealMatrix& TmpMatrixPsi01 = MatrixPsi01[n][m];
			      RealMatrix& TmpMatrixPsi10 = MatrixPsi10[n][m];	
			      
			      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
				{	      
				  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
				  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
				    {	      
				      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
				      int N2;
				      int N1;
				      N2 = (2 * ((q - m) - (p - n)) - 2 + NValueShift) / 2;
				      N1 = N2 + QValue;
				      if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]))
					  && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][1]))
					  && ((V0MatrixIndex - V0MatrixIndexShift) == (n - m - 1)))
					{ 
//					  cout << "x1=0 x2=1 q=" << q << " m=" << m << " p=" << p << " n=" << n << " N2=" << N2 << " N1=" << N1 << " " << " NValueShift=" << NValueShift << endl;
					  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
					    {
					      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
						{
						  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, N1, p, ChargedIndex1, NeutralIndex1)];
						  ++TmpTotalNbrElements;
						}
					    }
					  
					}
				      N2 = (2 * ((q - m) - (p - n)) + NValueShift) / 2;
				      N1 = N2 + QValue;
				      if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]))
					  && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][0]))
					  && ((V0MatrixIndex - V0MatrixIndexShift) == (n - m)))
					{ 
//					  cout << "x1=1 x2=0 q=" << q << " m=" << m << " p=" << p << " n=" << n << " N2=" << N2 << " N1=" << N1 << " " << " NValueShift=" << NValueShift << endl;
					  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
					    {
					      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
						{
						  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, N1, p, ChargedIndex1, NeutralIndex1)];
						  ++TmpTotalNbrElements;
						}
					    }
					}
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
      
      if (TmpTotalNbrElements > 0l)
	{
	  SparseRealMatrix V0Matrix (MatrixSize, MatrixSize, TmpNbrElementPerRow);
	  for (int n = 0; n <= this->PLevelShift; ++n)
	    {
	      for (int p = 0; p <= this->PLevelShift; ++p)
		{
		  if (((this->PLevelShift + p - n) <= this->PLevel) && ((this->PLevelShift + p - n) >= 0))
		    {
		      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
		      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
		      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
		      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
		      for (int m = 0; m <= this->PLevelShift; ++m)
			{
			  for (int q = 0; q <= this->PLevelShift; ++q)
			    {
			      if (((this->PLevelShift + q - m) <= this->PLevel) && ((this->PLevelShift + q - m) >= 0))
				{
				  double TmpProjectionFactor = 1.0;
				  if (this->CylinderFlag)
				    {
				      TmpProjectionFactor = exp (this->Kappa * this->Kappa * ((double) ((n - m) * (n - m))));
				    }
				  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
				  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[m];
				  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[m];
				  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[m];
				  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[n][m];
				  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[n][m];	
				  
				  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
				    {	      
				      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
				      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
					{	      
					  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
					  int N2;
					  int N1;
					  N2 = (2 * ((q - m) - (p - n)) - 2 + NValueShift) / 2;
					  N1 = N2 + QValue;
					  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][0]))
					      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][1]))
					      && ((V0MatrixIndex - V0MatrixIndexShift) == (n - m - 1)))
					    { 
					      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
						{
						  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
						    {
						      double Tmp = 0.0;
						      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
							{
							  double Tmp1 = 0.0;			      
							  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
							    {
							      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
							    }
							  Tmp += TmpOrthogonalBasisIdentity1(NeutralIndex3, NeutralIndex1) * Tmp1;
							}
						      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, p, q, Coef);
						      Tmp *= TmpProjectionFactor;
						      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 0, N1, p, ChargedIndex1, NeutralIndex1),
										this->Get2RMatrixIndexV2(this->PLevelShift + q - m, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
						    }
						  
						}
					    }
					  N2 = (2 * ((q - m) - (p - n)) + NValueShift) / 2;
					  N1 = N2 + QValue;
					  if (((N1 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + p - n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + p - n][1]))
					      && ((N2 >= this->NInitialValuePerPLevelCFTSector[this->PLevelShift + q - m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[this->PLevelShift + q - m][0]))
					      && ((V0MatrixIndex - V0MatrixIndexShift) == (n - m)))
					    { 
					      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
						{
						  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
						    {
						      double Tmp = 0.0;
						      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
							{
							  double Tmp1 = 0.0;			      
							  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
							    {
							      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentity2(NeutralIndex4, NeutralIndex2);				  
							    }
							  Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
							}
						      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, p, q, Coef);
						      Tmp *= TmpProjectionFactor;
						      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(this->PLevelShift + p - n, 1, N1, p, ChargedIndex1, NeutralIndex1),
										this->Get2RMatrixIndexV2(this->PLevelShift + q - m, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
						    }
						}
					    }
					}
				    }
				}
			    }
			}	 
		    }     
		}
	    }
 	  V0Matrices[V0MatrixIndex] = V0Matrix;
  	  cout << "V0[" << V0MatrixIndex << "] = " << endl;
  	  V0Matrices[V0MatrixIndex].PrintNonZero(cout, TmpLabels, TmpLabels) << endl;
	}
    }

  cout.precision(14);
  double** TmpProjectorCoefficients = new double*[nbrFluxQuanta + 1];
  if (true)
    {
      BinomialCoefficients TmpCoef (nbrFluxQuanta);
      for (int i = 0; i <= nbrFluxQuanta; ++i)
	{
	  TmpProjectorCoefficients[i] = new double[NbrV0MatrixIndices];
	  for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
	    {
	      if (((i + (V0MatrixIndex - V0MatrixIndexShift)) >= 0) && ((i + (V0MatrixIndex - V0MatrixIndexShift)) <= nbrFluxQuanta))
		{
		  TmpProjectorCoefficients[i][V0MatrixIndex] =  (sqrt(4.0 * M_PI / ((double) (nbrFluxQuanta + 1)) * TmpCoef(nbrFluxQuanta, i)) 
								 / TmpCoef(nbrFluxQuanta, i + (V0MatrixIndex - V0MatrixIndexShift)));
		}
	      else
		{
		  TmpProjectorCoefficients[i][V0MatrixIndex] = 0.0;
		}
	    }
	}
    }
  else
    {
//       FactorialCoefficient TmpCoef;
//       for (int i = 0; i <= nbrFluxQuanta; ++i)
// 	{
// 	  TmpProjectorCoefficients[i] = new double[NbrV0MatrixIndices];
// 	  for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
// 	    {
// 	      if ((i + (V0MatrixIndex - V0MatrixIndexShift)) >= 0)
// 		{
// 		  TmpCoef.SetToOne();
// 		  TmpCoef.FactorialMultiply(i + (V0MatrixIndex - V0MatrixIndexShift));
// 		  TmpCoef.FactorialMultiply(i + (V0MatrixIndex - V0MatrixIndexShift));
// 		  TmpCoef.FactorialDivide(i);
// 		  TmpCoef.Power2Multiply (2 * (i + (V0MatrixIndex - V0MatrixIndexShift)) - i);
// 		  cout << i << " " << V0MatrixIndex << " : " << TmpCoef << endl;
// 		  TmpProjectorCoefficients[i][V0MatrixIndex] =  sqrt(2.0 * M_PI * TmpCoef.GetNumericalValue());
// 		}
// 	      else
// 		{
// 		  TmpProjectorCoefficients[i][V0MatrixIndex] = 0.0;
// 		}
// 	    }
// 	}


      // cylinder
//       double Perimeter = 5.0;
//       this->Kappa =  (2.0 * M_PI) / Perimeter;
//       for (int i = 0; i <= nbrFluxQuanta; ++i)
// 	 {
//  	  TmpProjectorCoefficients[i] = new double[NbrV0MatrixIndices];
//  	  for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
//  	    {
// 	      TmpProjectorCoefficients[i][V0MatrixIndex] =  exp (0.5 * this->Kappa * this->Kappa * ((double) ((2 * (i + V0MatrixIndex - V0MatrixIndexShift) - nbrFluxQuanta) * 
// 													      (2 * (i + V0MatrixIndex - V0MatrixIndexShift) - nbrFluxQuanta)
// 													      - ((2 * i  - nbrFluxQuanta) * (2 * i - nbrFluxQuanta)))));
// 	    }
// 	 }
      
    }

  SparseRealMatrix** TmpMatrices = new SparseRealMatrix*[nbrFluxQuanta + 1];
  SparseRealMatrix* TmpMatrices2 = new SparseRealMatrix[NbrV0MatrixIndices];
  double* TmpCoefficients  = new double[NbrV0MatrixIndices];
  for (int i = 0; i <= nbrFluxQuanta; ++i)
    {
      TmpMatrices[i] = new  SparseRealMatrix[this->NbrBMatrices];
      TmpMatrices[i][0] = TmpB0Matrix;
      int TmpNbrV0Matrices = 0;
      for (int m = 0; m < NbrV0MatrixIndices; ++m)
	{
	  if ((TmpProjectorCoefficients[i][m] != 0.0) && (V0Matrices[m].GetNbrRow() > 0))
	    {
	      TmpMatrices2[TmpNbrV0Matrices] = V0Matrices[m];
	      TmpCoefficients[TmpNbrV0Matrices] = TmpProjectorCoefficients[i][m];
	      ++TmpNbrV0Matrices;
	    }
	}
//       cout << "B[0," << i << "] = " << endl;
//       TmpMatrices[i][0].PrintNonZero(cout) << endl;
      if (TmpNbrV0Matrices > 0)
	{
	  SparseRealMatrix TmpV0Matrix = SparseRealMatrixLinearCombination(TmpNbrV0Matrices, TmpCoefficients, TmpMatrices2);
	  for (int m = 1; m < this->NbrBMatrices; ++m)
	    {
	      TmpMatrices[i][m] = MemoryEfficientMultiply(TmpMatrices[i][m - 1], TmpV0Matrix);
	    }
   	  cout << "B[1," << i << "] = " << endl;
    	  TmpMatrices[i][1].PrintNonZero(cout, TmpLabels, TmpLabels) << endl;
	}
    }
  delete[] TmpMatrices2;
  delete[] TmpCoefficients;
  
  for (int i = 0; i <= nbrFluxQuanta; ++i)
    {
      delete[] TmpProjectorCoefficients[i];
    }
  delete[] TmpProjectorCoefficients;

  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
  for (int i = 0; i <= this->PLevelShift; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete U1BosonBasis[i];
      delete[] RationalMultiplicityFactor[i];
      delete[] MultiplicityFactor[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
  delete[] RationalMultiplicityFactor;
  delete[] MultiplicityFactor;
  return TmpMatrices;
}

