////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                                                                            //
//                        last modification : 27/09/2012                      //
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
#include "Tools/FQHEMPS/FQHEMPSClustered2ROptimizedMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

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

FQHEMPSClustered2ROptimizedMatrix::FQHEMPSClustered2ROptimizedMatrix()
{
  this->UseRationalFlag = true;
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2ROptimizedMatrix::FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool useRational,
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_opt_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation
  
FQHEMPSClustered2ROptimizedMatrix::FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool useRational, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->Kappa = kappa;
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_opt_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(cftDirectory, architecture);
}


// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2ROptimizedMatrix::FQHEMPSClustered2ROptimizedMatrix(int pLevel, int nbrBMatrices, char* fileName, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->PLevel = pLevel;
  
  ConfigurationParser StateDefinition;
  if (StateDefinition.Parse(fileName) == false)
    {
      StateDefinition.DumpErrors(cout) << endl;
    }
  else
    {
      bool ErrorFlag = true;
      ErrorFlag = StateDefinition.GetAsSingleInteger("RIndex", this->RIndex);
      ErrorFlag = StateDefinition.GetAsSingleInteger("LaughlinIndex", this->LaughlinIndex);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightIdentity", this->WeightIdentity);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPsi", this->WeightPsi);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
      if (StateDefinition["UseRational"] != 0)	
	this->UseRationalFlag = false;
      else
	ErrorFlag = StateDefinition.GetAsBoolean("UseRational", this->UseRationalFlag);
      if (StateDefinition["EMatrixDegeneracy"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleInteger("EMatrixDegeneracy", this->TransferMatrixDegeneracy);	  
	}
      else
	{
	  this->TransferMatrixDegeneracy = this->RIndex + 2;
	}
      if (StateDefinition["Name"] != 0)
	{
	  this->BMatrixOutputName = new char[strlen(StateDefinition["Name"]) + 1]; 
	  strcpy(this->BMatrixOutputName, StateDefinition["Name"]);
	}
      else
	{
	  this->BMatrixOutputName = new char[256]; 
	  sprintf(this->BMatrixOutputName, "clustered_opt_k_2_r_%d", this->RIndex);
	}
      if (StateDefinition["PsiSquareMatrixElement"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("PsiSquareMatrixElement", this->SquareMatrixElementNormalization);
	}
      else
	{
	  this->SquareMatrixElementNormalization = LongRational(1, 1);
	}
      this->MatrixElementNormalization = sqrt(fabs(this->SquareMatrixElementNormalization.GetNumericalValue()));
      if (ErrorFlag == true)
	{
	  this->WeightPrimaryFieldMatrixElement = this->WeightPsi;
	  this->CreateBMatrices(StateDefinition["CFTMatrixDirectory"], architecture);
	}
    }
}

// constructor from stored B matrices
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2ROptimizedMatrix::FQHEMPSClustered2ROptimizedMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
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
  sprintf(this->BMatrixOutputName, "clustered_opt_k_2_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2ROptimizedMatrix::~FQHEMPSClustered2ROptimizedMatrix()
{
}
  
// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2ROptimizedMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;

  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();

  double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];

  LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevel + 1];
  double** MultiplicityFactor = new double*[this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
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
  for (int i = 0; i <= this->PLevel; ++i)
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

  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [2];
  this->NeutralSectorDimension[0] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }

  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  double ExtraCylinderFactor = 1.0;
  if ((this->RIndex & 1) == 0)
    {
      QValue = 1 + (this->RIndex / 2);
      this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
      NValueShift = 2 * this->PLevel - 1;
      QValueDenominator = 1;
    }
  else
    {
      QValue = 2 + this->RIndex;
      this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex;
      NValueShift = 4 * this->PLevel - 2;
      QValueDenominator = 2;
      ExtraCylinderFactor = 4.0;
      QValue /= 2;
      this->NbrNValue /= 2;
      ++this->NbrNValue;
      NValueShift /= 2;
      QValueDenominator /= 2;
    }

  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
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

  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;

  // B^[0]  matrix evaluation
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}	      
	    }
	}
    }

  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
			      if (this->CylinderFlag)
				Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									 + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
			      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
 							    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
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

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift);
			      N1 = N2 + QValue - 1;
			    }			  
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
				    }
				}

			    }
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift);
			      N1 = N2 + QValue - 1;
			    }
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
			    { 
			      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				{
				  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				    {
				      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1)];
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

  BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift);
			      N1 = N2 + QValue - 1;
			    }			  
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
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
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
				    }
				  
				}
			    }
			  if ((this->RIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift);
			      N1 = N2 + QValue - 1;
			    }
			  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
			      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
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
				      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				      if (this->CylinderFlag)
					Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
										       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
										       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
				      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
				    }
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
  for (int i = 0; i <= this->PLevel; ++i)
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

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSClustered2ROptimizedMatrix::GetQValueCFTSectorShift(int cftSector)
{
  if ((this->RIndex & 1) == 1)
    return 0;
  if (cftSector == 0)
    return 0;
  return ((this->RIndex + 2) / 2);
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSClustered2ROptimizedMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel))
    return 0;
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][0]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][0]))
    {
      int Tmp = this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
      return Tmp;
    }
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
    return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
  return 0;
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSClustered2ROptimizedMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
      (qValue > this->NLastValuePerPLevelCFTSector[pLevel][cftSector]) || (cftSector > 1) ||  (cftSector < 0))
    return 0;
  return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]];
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSClustered2ROptimizedMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][0]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][0]))
    {
      if (localIndex < this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]])
	{	  
	  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]  + localIndex);
	}
      else
	{
	  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]]  
		  + (localIndex - this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]));
	}
    }
  else
    {
      return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]] + localIndex);
    }
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int FQHEMPSClustered2ROptimizedMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{
  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]]  + localIndex);
}


// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSClustered2ROptimizedMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      if ((this->RIndex & 1) == 0)
	rowIndex = this->PLevel + (this->RIndex / 2) - MinQ;
      else
	rowIndex = this->PLevel + (this->RIndex - 1) / 2 - MinQ;      
      columnIndex = rowIndex;
    }
  else
    {
      rowIndex = this->PLevel + this->RIndex - MinQ;
      columnIndex = this->PLevel - MinQ;
    }
}


// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2ROptimizedMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
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
  minQ = TmpMinQ;
  maxQ = TmpMaxQ;
  cout << "range at p=" << pLevel << ", x=" << cftSector << " : " << minQ << " " << maxQ << " (" << this->NbrNValue << ")" << endl;   
}

// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int FQHEMPSClustered2ROptimizedMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  nbrFluxQuanta += this->RIndex + 1;
  nbrFluxQuanta *= 2;
  return (nbrFluxQuanta / (this->RIndex + 2));
}

