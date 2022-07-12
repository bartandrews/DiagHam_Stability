////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the N=1 superconformal states           //
//                                                                            //
//                        last modification : 11/03/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSN1SuperconformalMatrix.h"
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

FQHEMPSN1SuperconformalMatrix::FQHEMPSN1SuperconformalMatrix()
{
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

FQHEMPSN1SuperconformalMatrix::FQHEMPSN1SuperconformalMatrix(int pLevel, int nbrBMatrices, char* fileName,  bool trimChargeIndices, bool cylinderFlag, double kappa, 
							     AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->PLevel = pLevel;
  this->RIndex = 6;
  this->NbrCFTSectors = 2;

  ConfigurationParser StateDefinition;
  if (StateDefinition.Parse(fileName) == false)
    {
      StateDefinition.DumpErrors(cout) << endl;
    }
  else
    {
      bool ErrorFlag = true;
      ErrorFlag = StateDefinition.GetAsSingleInteger("LaughlinIndex", this->LaughlinIndex);
      if ((StateDefinition["MinimalModelP"] != 0) && (StateDefinition["MinimalModelQ"] != 0))
	{
	  int PValue;
	  int QValue;
	  ErrorFlag = StateDefinition.GetAsSingleInteger("MinimalModelP", PValue);
	  ErrorFlag = StateDefinition.GetAsSingleInteger("MinimalModelQ", QValue);
	  this->CentralCharge = LongRational(3l * ((5l * ((long) PValue) * ((long) QValue)) 
						   - (2l * ((long) PValue) * ((long) PValue)) 
						   - (2l * ((long) QValue) * ((long) QValue))), 
					     2l * ((long) PValue) * ((long) QValue));
	  this->WeightIdentity = 0l;
	}
      else
	{
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightIdentity", this->WeightIdentity);
	}
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
	  char* TmpCentralCharge = this->CentralCharge.GetString('_');
	  this->BMatrixOutputName = new char[256 + strlen(TmpCentralCharge)]; 
	  sprintf(this->BMatrixOutputName, "n1superconformal_c_%s", TmpCentralCharge);
	  delete[] TmpCentralCharge;
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
      this->UseRationalFlag = true;
      if (StateDefinition["UseNonRational"])
	{
	  bool TmpFlag;
	  StateDefinition.GetAsBoolean("UseNonRational", TmpFlag);
	  this->UseRationalFlag = !TmpFlag;
	}
      if (this->UseRationalFlag)
	{
	  cout << "using rational numbers" << endl;
	}
      else
	{
	  cout << "using double numbers" << endl;
	}
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

FQHEMPSN1SuperconformalMatrix::FQHEMPSN1SuperconformalMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName,  bool trimChargeIndices, bool cylinderFlag, double kappa)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->LoadMatrices(fileName);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  char* TmpCentralCharge = this->CentralCharge.GetString('_');
  this->BMatrixOutputName = new char[256 + strlen(TmpCentralCharge)]; 
  sprintf(this->BMatrixOutputName, "n1superconformal_c_%s", TmpCentralCharge);
  delete[] TmpCentralCharge;
}

// destructor
//

FQHEMPSN1SuperconformalMatrix::~FQHEMPSN1SuperconformalMatrix()
{
}
  
// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSN1SuperconformalMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  this->EffectivePLevel = (2 * this->PLevel + 3);
  LongRational InvCentralCharge3 (3l, 1l);
  InvCentralCharge3 /= this->CentralCharge;
  long* Partition = new long[2 * (this->EffectivePLevel + 1) + 1];
  unsigned long* TmpPartition = new unsigned long [this->EffectivePLevel + 2];

  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double InvCentralCharge3Numerical = InvCentralCharge3.GetNumericalValue();
  double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  BosonOnDiskShort** SupersymmetricU1BosonBasis = new BosonOnDiskShort* [this->EffectivePLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->EffectivePLevel + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->EffectivePLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 2];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 2];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 2];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 2];
  char* TmpScalarProductIdentityFileName = 0; 
  char* TmpScalarProductPsiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(cftDirectory)];
    }

  LongRational** RationalMultiplicityFactor = new LongRational*[this->EffectivePLevel + 1];
  double** MultiplicityFactor = new double*[this->EffectivePLevel + 1];
  for (int i = 0; i <= this->EffectivePLevel; ++i)
    {
      BosonOnDiskShort TmpU1BosonBasis (i, i, this->EffectivePLevel + 1);
      int EffectiveDimension = 0;
      bool* U1BosonBasisKeptFlag = new bool[TmpU1BosonBasis.GetHilbertSpaceDimension()];
      for (int n = 0; n < TmpU1BosonBasis.GetHilbertSpaceDimension(); ++n)
	{	  
	  TmpU1BosonBasis.GetOccupationNumber(n, TmpPartition);	    
	  bool Flag = true;
	  for (int k = 1; (k <= i) && (Flag == true); k += 2)
	    if (TmpPartition[k] > 1l)
	      Flag = false;
	  if (Flag == true)
	    ++EffectiveDimension;
	  U1BosonBasisKeptFlag[n] = Flag;
	}
      SupersymmetricU1BosonBasis[i] = new BosonOnDiskShort(TmpU1BosonBasis, EffectiveDimension, U1BosonBasisKeptFlag);
      delete[] U1BosonBasisKeptFlag;
      if (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() > 0)
	{
	  RationalMultiplicityFactor[i] = new LongRational[SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension()];
	  MultiplicityFactor[i] = new double[SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension()];
	  for (int j = 0; j < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
	    {
	      SupersymmetricU1BosonBasis[i]->GetOccupationNumber(j, TmpPartition);	    
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
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
    }

  cout << "Primary field conformal weight: " <<   this->WeightIdentity << endl;
  for (int i = 0; i <= this->EffectivePLevel; ++i)
    {
      if ((i & 1) == 0)
	cout << "Level = " <<  (i / 2) << endl;
      else
	cout << "Level = " <<  i << "/2" << endl;
      RationalScalarProductIdentity[i] = LongRationalMatrix(SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
	{
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
	    }
	  else
	    {
	      ScalarProductIdentity[i].ReadMatrix(TmpScalarProductIdentityFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, SupersymmetricU1BosonBasis, i, CentralCharge12, InvCentralCharge3,
						     this->WeightIdentity,
						     RationalScalarProductIdentity,  i - 1);
	      Operation1.ApplyOperation(architecture);
	      RationalScalarProductIdentity[i] = Operation1.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, SupersymmetricU1BosonBasis, i, CentralCharge12Numerical, InvCentralCharge3Numerical,
						     WeightIdentityNumerical,
						     ScalarProductIdentity,  i - 1);
	      Operation1.ApplyOperation(architecture);
	      ScalarProductIdentity[i] = Operation1.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	}  

      RealSymmetricMatrix TmpMatrix;
      TmpMatrix.Copy(ScalarProductIdentity[i]);
      if (this->UseRationalFlag == true)
 	{
 	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductIdentity[i].GetNbrRow(), RationalScalarProductIdentity[i].GetNbrColumn());
 	  for (int k = 0; k < RationalScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductIdentity[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
       else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductIdentity[i].GetNbrRow(), ScalarProductIdentity[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductIdentity[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductIdentity[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		TmpMatrix.SetMatrixElement(k, l, Tmp);
	      }
	}
      RealMatrix TmpBasis(SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      double Error = 0.0;
      for (int n = 0; n < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      int Count  = 0;
      for (int n = 0; n < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
      cout << "nbr of null vectors identity sector = " << Count << " (" << (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  if ((i & 1) == 0)
	    {
	      OrthogonalBasisIdentityLeft[i / 2] = RealMatrix (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      OrthogonalBasisIdentityRight[i / 2] = RealMatrix (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      Count = 0;
	      for (int n = 0; n < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisIdentityLeft[i / 2][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisIdentityRight[i / 2][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisIdentityLeft[i / 2][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisIdentityRight[i / 2][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisIdentityLeft[i / 2][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisIdentityRight[i / 2][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
		  }
	    }
	  else
	    {
	      if (i >= 3)
		{
		  OrthogonalBasisPsiLeft[(i - 3) / 2] = RealMatrix (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
		  OrthogonalBasisPsiRight[(i - 3) / 2] = RealMatrix (SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(), SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
		  Count = 0;
		  for (int n = 0; n < SupersymmetricU1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		    if (fabs(TmpDiag(n, n)) > Error)
		      {
			OrthogonalBasisPsiLeft[(i - 3)/ 2][Count].Copy(TmpBasis[n]);
			OrthogonalBasisPsiRight[(i - 3) / 2][Count].Copy(TmpBasis[n]);
			if (TmpDiag(n, n) > 0)
			  {
			    OrthogonalBasisPsiLeft[(i - 3) / 2][Count] /=  sqrt(TmpDiag(n, n));
			    OrthogonalBasisPsiRight[(i - 3) / 2][Count] /=  sqrt(TmpDiag(n, n));
			  }
			else
			  {
			    OrthogonalBasisPsiLeft[(i - 3) / 2][Count] /=  sqrt(-TmpDiag(n, n));
			    OrthogonalBasisPsiRight[(i - 3) / 2][Count] /=  -sqrt(-TmpDiag(n, n));
			  }
			++Count;
		      }
		}
	    }
	}
      else
	{
	  if ((i & 1) == 0)
	    {
	      OrthogonalBasisIdentityLeft[i / 2] = RealMatrix();
	      OrthogonalBasisIdentityRight[i / 2] = RealMatrix();
	    }
	  else
	    {
	      if (i >= 3)
		{
		  OrthogonalBasisPsiLeft[(i - 3) / 2] = RealMatrix();
		  OrthogonalBasisPsiRight[(i - 3) / 2] = RealMatrix();
		}
	    }
	}
      cout << "---------------------------------" << endl;
    }
  
   for (int i = 0; i <= this->EffectivePLevel; ++i)
     {
       if (this->UseRationalFlag == true)
  	{
  	  for (int k = 0; k < RationalScalarProductIdentity[i].GetNbrRow(); ++k)
 	    for (int l = 0; l < RationalScalarProductIdentity[i].GetNbrColumn(); ++l)
 	      {
 		RationalScalarProductIdentity[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
 	      }
  	  ScalarProductIdentity[i] = RationalScalarProductIdentity[i];
  	}
       else
 	{
 	  for (int k = 0; k < ScalarProductIdentity[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductIdentity[i].GetNbrColumn(); ++l)
 	      {
 		double Tmp;
 		ScalarProductIdentity[i].GetMatrixElement(k, l, Tmp);
 		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
 		ScalarProductIdentity[i].SetMatrixElement(k, l, Tmp);
 	      }
 	}
     }

  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [this->NbrCFTSectors];
  for (int i = 0; i < this->NbrCFTSectors; ++i)
    this->NeutralSectorDimension[i] = new int [this->PLevel + 1];
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
  QValue = 1 + (this->RIndex / 2);
  this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
  NValueShift = 2 * this->PLevel - 1;
  QValueDenominator = 1;

  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (this->WeightPsi);
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  cout << "Levels = " <<  i << " " << (2 * j + 3) << "/2" << endl;
	  if (cftDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_matrixelement_rrns_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_matrixelement_nsrr_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	      else
		{
		  sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_matrixelement_rrns_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		  sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_matrixelement_nsrr_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductIdentityFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
		}
	      else
		{
		  MatrixPsi01[i][j].ReadMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, SupersymmetricU1BosonBasis, 2 * i, 2 * j + 3, CentralCharge12, InvCentralCharge3,
							 this->WeightIdentity, 
							 RationalMatrixPsi01, 2 * (i - 1), 2 * j + 3);
		  Operation1.ApplyOperation(architecture);
		  RationalMatrixPsi01[i][j] = Operation1.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, SupersymmetricU1BosonBasis, 2 * i, 2 * j + 3, CentralCharge12Numerical, InvCentralCharge3Numerical,
							 WeightIdentityNumerical, 
							 MatrixPsi01, 2 * (i - 1), 2 * j + 3);
		  Operation1.ApplyOperation(architecture);
		  MatrixPsi01[i][j] = Operation1.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi01[i][j].WriteMatrix(TmpScalarProductIdentityFileName);
		    }
		}
	    }
	  cout << "Levels = " <<  (2 * i + 3) << "/2 " << j << endl;
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
		}
	      else
		{
		  MatrixPsi10[i][j].ReadMatrix(TmpScalarProductPsiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, SupersymmetricU1BosonBasis, 2 * i + 3, 2 * j, CentralCharge12, InvCentralCharge3,
							 this->WeightIdentity,
							 RationalMatrixPsi10,  2 * i + 1, 2 * j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi10[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, SupersymmetricU1BosonBasis, 2 * i + 3, 2 * j, CentralCharge12Numerical, InvCentralCharge3Numerical,
							 WeightIdentityNumerical,
							 MatrixPsi10,  2 * i + 1, 2 * j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi10[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi10[i][j].WriteMatrix(TmpScalarProductPsiFileName);
		    }
		}
	    }
	}
    }
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < RationalMatrixPsi01[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi01[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi01[i][j][l][k] *= (RationalMultiplicityFactor[2 * i][k] * RationalMultiplicityFactor[2 * j + 3][l]);
		  }
	      for (int k = 0; k < RationalMatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi10[i][j][l][k] *= (RationalMultiplicityFactor[2 * i + 3][k] * RationalMultiplicityFactor[2 * j][l]);
		  }
	      MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	      MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
	    }
	  else
	    {
	      for (int k = 0; k < MatrixPsi01[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi01[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi01[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[2 * i][k] * MultiplicityFactor[2 * j + 3][l]);
		    MatrixPsi01[i][j].SetMatrixElement(k, l, Tmp);
		  }
	      for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi10[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[2 * i + 3][k] * MultiplicityFactor[2 * j][l]);
		    MatrixPsi10[i][j].SetMatrixElement(k, l, Tmp);
		  }
	    }
	  MatrixPsi01[i][j] *= this->MatrixElementNormalization;
	  MatrixPsi10[i][j] *= this->MatrixElementNormalization;
	}
    }
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
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[2 * (i - p)];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductIdentity[2 * (i - p) + 3];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
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
	    }
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
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

  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = SupersymmetricU1BosonBasis[2 * (i - p)];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[2 * (i - p)];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductIdentity[2 * (i - p) + 3];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		{
		  TmpSpaceNeutral = SupersymmetricU1BosonBasis[2 * (i - p)];
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
		  TmpSpaceNeutral = SupersymmetricU1BosonBasis[2 * (i - p) + 3];
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
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
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
			  N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			  N1 = N2 + QValue - 1;
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
			  N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			  N1 = N2 + QValue - 1;
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
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
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
			  BosonOnDiskShort* TmpSpaceNeutral1 = SupersymmetricU1BosonBasis[2 * (i - p)];
			  BosonOnDiskShort* TmpSpaceNeutral2 = SupersymmetricU1BosonBasis[2 * (j - q) + 3];
			  int N2;
			  int N1;
			  N2 = (2 * (j - i) + this->RIndex + 1 + NValueShift) / 2;
			  N1 = N2 + QValue - 1;
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
			  TmpSpaceNeutral1 = SupersymmetricU1BosonBasis[2 * (i - p) + 3];
			  TmpSpaceNeutral2 = SupersymmetricU1BosonBasis[2 * (j - q)];
			  N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			  N1 = N2 + QValue - 1;
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
  for (int i = 0; i <= this->PLevel; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete U1BosonBasis[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
}

// compute the scalar product matrices of the Virasoro descendant
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// invCentralCharge3 = reference on the value of three divided by the central charge
// weight = weight of the primary field that is considered
// return value = scalar product

LongRational FQHEMPSN1SuperconformalMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
										    LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  if ((partition[0] & 1l) == 0)
	    {
	      LongRational Tmp1 (centralCharge12);
	      Tmp1 *= (partition[0] / 2l) * ((partition[0] / 2l) * (partition[0] / 2l) - 1l);
	      LongRational Tmp2 (weight);
	      Tmp2 *= partition[0];
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	  else
	    {
	      LongRational Tmp1 (invCentralCharge3);
	      Tmp1 *= weight;
	      LongRational Tmp2 (partition[0] * partition[0] - 1l, 8l);
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	}
    }
  LongRational Tmp(0l);
  if ((partition[position - 1] + partition[position]) == 0)
    {
      if ((partition[position] & 1l) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((Store / 2l) * ((Store  / 2l) * (Store / 2l) - 1l)) * centralCharge12)
		   + ((Store / 2l) * (2l * weight - TmpLength))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((LongRational((Store * Store - 1l), 8l)
		   + invCentralCharge3 * (weight - LongRational(TmpLength, 2l))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
    }
  else
    {
      if (((partition[position - 1] & 1l) == 0l) && ((partition[position] & 1l) == 0l))
	{
	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += (((Store1 - Store2) / 2l)
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position, centralCharge12, invCentralCharge3, weight));
	    }
	  else
	    {
	      Tmp += (((Store1 - Store2) / 2l) 
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position - 1, centralCharge12, invCentralCharge3, weight));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      else
	{
	  if (((partition[position - 1] & 1l) == 1l) && ((partition[position] & 1l) == 1l))
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  Tmp += ( invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position, centralCharge12, invCentralCharge3, weight));
		}
	      else
		{
		  Tmp += (invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position - 1, centralCharge12, invCentralCharge3, weight));
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	  else
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += (LongRational((Store1 - 2l * Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight));
		    }
		  else
		    {
		      Tmp += (LongRational((2l * Store1 - Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight));
		    }
		}
	      else
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += (LongRational((Store1 - 2l * Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight));
		    }
		  else
		    {
		      Tmp += (LongRational((2l * Store1 - Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight));
		    }
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	}
    }

  if ((partition[position - 1] & partition[position] & 1l) == 0l)
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  else
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp -= this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// invCentralCharge3 = reference on the value of three divided by the central charge
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

LongRational FQHEMPSN1SuperconformalMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
										    LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight,
										    LongRationalMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
										    BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  if ((partition[0] & 1l) == 0)
	    {
	      LongRational Tmp1 (centralCharge12);
	      Tmp1 *= (partition[0] / 2l) * ((partition[0] / 2l) * (partition[0] / 2l) - 1l);
	      LongRational Tmp2 (weight);
	      Tmp2 *= partition[0];
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	  else
	    {
	      LongRational Tmp1 (invCentralCharge3);
	      Tmp1 *= weight;
	      LongRational Tmp2 (partition[0] * partition[0] - 1l, 8l);
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	}
    }
  int TmpPosition = 0;
  while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
    ++TmpPosition;
  if (TmpPosition == position)
    {
      TmpPosition = 1;
      int TmpPLevel1 = partition[0];
      bool FlagSorted = true;
      if ((partition[0] & 1) == 1)
	{
	  while ((TmpPosition < position) && (FlagSorted == true) && ((partition[TmpPosition] & 1) == 1))
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if (TmpPosition < position)
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      ++TmpPosition;
	    }
	}
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if ((partition[TmpPosition - 1] < partition[TmpPosition]) || ((partition[TmpPosition] & 1) == 1))
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	{     
	  int TmpPLevel2 = -partition[TmpPosition];	  
	  FlagSorted = true;
	  ++TmpPosition;
	  if ((partition[TmpPosition - 1] & 1) == 0)
	    {
	      while ((TmpPosition < partitionLength) && (FlagSorted == true) && ((partition[TmpPosition] & 1) == 0))
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  if (partition[TmpPosition - 1] < partition[TmpPosition])
		    FlagSorted = false;
		  ++TmpPosition;
		}
	      if (TmpPosition < partitionLength)
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  ++TmpPosition;
		}
	    }
	  while ((TmpPosition < partitionLength) && (FlagSorted == true))
	    {
	      TmpPLevel2 -= partition[TmpPosition];
	      if ((partition[TmpPosition - 1] < partition[TmpPosition]) || ((partition[TmpPosition] & 1) == 0))
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	    {
	      for (int k = 0; k <= (this->EffectivePLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		{
		  temporaryOccupationNumber[partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel1 - position;
	      FlagSorted = true;
	      for (int k = 1; k <= (this->EffectivePLevel + 1); k += 2)
		if (temporaryOccupationNumber[k] > 1l)
		  FlagSorted = false;
	      if (FlagSorted == true)
		{
		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  for (int k = 0; k <= (this->EffectivePLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		    {
		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
		  FlagSorted = true;
		  for (int k = 1; k <= (this->EffectivePLevel + 1); k += 2)
		    if (temporaryOccupationNumber[k] > 1l)
		      FlagSorted = false;
		  if (FlagSorted == true)
		    {
		      LongRational Tmp;		      
		      int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		      return Tmp;
		    }
		}
	    }
	}
    }
  LongRational Tmp(0l);
  if ((partition[position - 1] + partition[position]) == 0)
    {
      if ((partition[position] & 1l) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((Store / 2l) * ((Store  / 2l) * (Store / 2l) - 1l)) * centralCharge12)
		   + ((Store / 2l) * (2l * weight - TmpLength))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight,
							       precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((LongRational((Store * Store - 1l), 8l)
		   + invCentralCharge3 * (weight - LongRational(TmpLength, 2l))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight,
							       precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
    }
  else
    {
      if (((partition[position - 1] & 1l) == 0l) && ((partition[position] & 1l) == 0l))
	{
	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += (((Store1 - Store2) / 2l)
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position, centralCharge12, invCentralCharge3, weight,
								     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	    }
	  else
	    {
	      Tmp += (((Store1 - Store2) / 2l) 
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position - 1, centralCharge12, invCentralCharge3, weight,
								     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      else
	{
	  if (((partition[position - 1] & 1l) == 1l) && ((partition[position] & 1l) == 1l))
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  Tmp += ( invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position, centralCharge12, invCentralCharge3, weight,
									 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		}
	      else
		{
		  Tmp += (invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position - 1, centralCharge12, invCentralCharge3, weight,
									 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	  else
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += (LongRational((Store1 - 2l * Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		  else
		    {
		      Tmp += (LongRational((2l * Store1 - Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		}
	      else
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += (LongRational((Store1 - 2l * Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		  else
		    {
		      Tmp += (LongRational((2l * Store1 - Store2), 4l)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	}
    }

  if ((partition[position - 1] & partition[position] & 1l) == 0l)
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight,
							  precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  else
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp -= this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight,
							  precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  return Tmp;
}

// compute the scalar product matrices of the Virasoro descendant
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// invCentralCharge3 = reference on the value of three divided by the central charge
// weight = weight of the primary field that is considered
// return value = scalar product

double FQHEMPSN1SuperconformalMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									      double& centralCharge12, double& invCentralCharge3, double& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  if ((partition[0] & 1l) == 0)
	    {
	      double Tmp1 = centralCharge12;
	      Tmp1 *= (partition[0] / 2l) * ((partition[0] / 2l) * (partition[0] / 2l) - 1l);
	      double Tmp2 = weight;
	      Tmp2 *= partition[0];
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	  else
	    {
	      double Tmp1 = invCentralCharge3;
	      Tmp1 *= weight;
	      double Tmp2 = ((double) partition[0] * partition[0] - 1l) /  8.0;
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	}
    }
  double Tmp = 0.0;
  if ((partition[position - 1] + partition[position]) == 0)
    {
      if ((partition[position] & 1l) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((Store / 2l) * ((Store  / 2l) * (Store / 2l) - 1l)) * centralCharge12)
		   + ((Store / 2l) * (2l * weight - TmpLength))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((double) (Store * Store - 1l))/  8l)
		   + invCentralCharge3 * (weight - (((double) TmpLength) / 2.0))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
    }
  else
    {
      if (((partition[position - 1] & 1l) == 0l) && ((partition[position] & 1l) == 0l))
	{
	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += (((Store1 - Store2) / 2l)
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position, centralCharge12, invCentralCharge3, weight));
	    }
	  else
	    {
	      Tmp += (((Store1 - Store2) / 2l) 
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position - 1, centralCharge12, invCentralCharge3, weight));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      else
	{
	  if (((partition[position - 1] & 1l) == 1l) && ((partition[position] & 1l) == 1l))
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  Tmp += ( invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position, centralCharge12, invCentralCharge3, weight));
		}
	      else
		{
		  Tmp += (invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position - 1, centralCharge12, invCentralCharge3, weight));
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	  else
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += ((((double) (Store1 - 2l * Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight));
		    }
		  else
		    {
		      Tmp += ((((double) (2l * Store1 - Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight));
		    }
		}
	      else
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += ((((double) (Store1 - 2l * Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight));
		    }
		  else
		    {
		      Tmp += ((((double) (2l * Store1 - Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight));
		    }
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	}
    }

  if ((partition[position - 1] & partition[position] & 1l) == 0l)
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  else
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp -= this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// invCentralCharge3 = reference on the value of three divided by the central charge
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

double FQHEMPSN1SuperconformalMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									      double& centralCharge12, double& invCentralCharge3, double& weight,
									      RealSymmetricMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
									      BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  if ((partition[0] & 1l) == 0)
	    {
	      double Tmp1 = centralCharge12;
	      Tmp1 *= (partition[0] / 2l) * ((partition[0] / 2l) * (partition[0] / 2l) - 1l);
	      double Tmp2 = weight;
	      Tmp2 *= partition[0];
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	  else
	    {
	      double Tmp1 = invCentralCharge3 ;
	      Tmp1 *= weight;
	      double Tmp2 = ((double) (partition[0] * partition[0] - 1l)) / 8.0;
	      Tmp1 += Tmp2;
	      return Tmp1;
	    }
	}
    }
  int TmpPosition = 0;
  while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
    ++TmpPosition;
  if (TmpPosition == position)
    {
      TmpPosition = 1;
      int TmpPLevel1 = partition[0];
      bool FlagSorted = true;
      if ((partition[0] & 1) == 1)
	{
	  while ((TmpPosition < position) && (FlagSorted == true) && ((partition[TmpPosition] & 1) == 1))
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if (TmpPosition < position)
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      ++TmpPosition;
	    }
	}
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if ((partition[TmpPosition - 1] < partition[TmpPosition]) || ((partition[TmpPosition] & 1) == 1))
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	{     
	  int TmpPLevel2 = -partition[TmpPosition];	  
	  FlagSorted = true;
	  ++TmpPosition;
	  if ((partition[TmpPosition - 1] & 1) == 0)
	    {
	      while ((TmpPosition < partitionLength) && (FlagSorted == true) && ((partition[TmpPosition] & 1) == 0))
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  if (partition[TmpPosition - 1] < partition[TmpPosition])
		    FlagSorted = false;
		  ++TmpPosition;
		}
	      if (TmpPosition < partitionLength)
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  ++TmpPosition;
		}
	    }
	  while ((TmpPosition < partitionLength) && (FlagSorted == true))
	    {
	      TmpPLevel2 -= partition[TmpPosition];
	      if ((partition[TmpPosition - 1] < partition[TmpPosition]) || ((partition[TmpPosition] & 1) == 0))
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
	    {
	      for (int k = 0; k <= (this->EffectivePLevel + 1); ++k)
		temporaryOccupationNumber[k] = 0x0ul;	  
	      for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		{
		  temporaryOccupationNumber[partition[TmpPosition]]++;	      
		}
	      temporaryOccupationNumber[0] = TmpPLevel1 - position;
	      FlagSorted = true;
	      for (int k = 1; k <= (this->EffectivePLevel + 1); k += 2)
		if (temporaryOccupationNumber[k] > 1l)
		  FlagSorted = false;
	      if (FlagSorted == true)
		{
		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  for (int k = 0; k <= (this->EffectivePLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		    {
		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
		  FlagSorted = true;
		  for (int k = 1; k <= (this->EffectivePLevel + 1); k += 2)
		    if (temporaryOccupationNumber[k] > 1l)
		      FlagSorted = false;
		  if (FlagSorted == true)
		    {
		      double Tmp;		      
		      int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		      return Tmp;
		    }
		}
	    }
	}
    }
  double Tmp = 0.0;
  if ((partition[position - 1] + partition[position]) == 0)
    {
      if ((partition[position] & 1l) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((Store / 2l) * ((Store  / 2l) * (Store / 2l) - 1l)) * centralCharge12)
		   + ((Store / 2l) * (2l * weight - TmpLength))) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight,
							       precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += (((((double) (Store * Store - 1l)) / 8.0)
		   + invCentralCharge3 * (weight - ((double) TmpLength) *0.5)) * 
		  this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, invCentralCharge3, weight,
							       precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
    }
  else
    {
      if (((partition[position - 1] & 1l) == 0l) && ((partition[position] & 1l) == 0l))
	{
	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += (((Store1 - Store2) / 2l)
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position, centralCharge12, invCentralCharge3, weight,
								     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	    }
	  else
	    {
	      Tmp += (((Store1 - Store2) / 2l) 
		      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								     position - 1, centralCharge12, invCentralCharge3, weight,
								     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      else
	{
	  if (((partition[position - 1] & 1l) == 1l) && ((partition[position] & 1l) == 1l))
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  Tmp += ( invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position, centralCharge12, invCentralCharge3, weight,
									 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		}
	      else
		{
		  Tmp += (invCentralCharge3
			  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									 position - 1, centralCharge12, invCentralCharge3, weight,
									 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	  else
	    {
	      long Store1 = partition[position - 1];
	      long Store2 = partition[position];
	      partition[position - 1] += partition[position];
	      for (int i = position + 1; i < partitionLength; ++i)
		partition[i - 1] = partition[i];
	      if ((Store1 + Store2) > 0)
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += ((((double) (Store1 - 2l * Store2)) /  4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		  else
		    {
		      Tmp += ((((double) ((2l * Store1 - Store2))) /  4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		}
	      else
		{
		  if ((Store1 & 1l) == 0l)
		    {
		      Tmp += ((((double) (Store1 - 2l * Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		  else
		    {
		      Tmp += ((((double) (2l * Store1 - Store2)) / 4.0)
			      * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
									     position - 1, centralCharge12, invCentralCharge3, weight,
									     precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
		    }
		}
	      for (int i = partitionLength - 1; i > position; --i)
		partition[i] = partition[i - 1];
	      partition[position] = Store2;
	      partition[position - 1] = Store1;
	    }
	}
    }

  if ((partition[position - 1] & partition[position] & 1l) == 0l)
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight,
							  precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  else
    {
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp -= this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, invCentralCharge3, weight,
							  precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;
    }
  return Tmp;
}

// compute the matrix elements of any primary field in the Virasoro descendant basis
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight1 = weight of the primary field that is considered for the left state
// weight2 = weight of the primary field that is considered for the right state
// weight = weight of the primary field whose matrix elements are computed
// return value = matrix element
  
LongRational FQHEMPSN1SuperconformalMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
								       int descendantPosition, int position, 
								       LongRational& centralCharge12, LongRational& weight1, 
								       LongRational& weight2, LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
	}
      Tmp *= 1l - (2l * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition < position)
    {
      LongRational Tmp(0l);
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
						       position - 1, centralCharge12, weight1, weight2, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
 	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position, centralCharge12, weight1, weight2, weight));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position - 1, centralCharge12, weight1, weight2, weight));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
						  centralCharge12, weight1, weight2, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  LongRational Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight);
  LongRational Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					       centralCharge12, weight1, weight2, weight);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

// compute the matrix elements of any primary field in the Virasoro descendant basis
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// descendantPosition = location of the primary field
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight1 = weight of the primary field that is considered for the left state
// weight2 = weight of the primary field that is considered for the right state
// weight = weight of the primary field whose matrix elements are computed
// precomputedDescendantMatrixElement = matrices where matrix elements computed for previous levels are stored
// precomputedDescendantMatrixElementMaxLeftPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the left entry
// precomputedDescendantMatrixElementMaxRightPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the right entry
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = matrix element
  
LongRational FQHEMPSN1SuperconformalMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
									    int descendantPosition, int position, 
									    LongRational& centralCharge12, LongRational& weight1, 
									    LongRational& weight2, LongRational& weight,
									    LongRationalMatrix** precomputedDescendantMatrixElement, 
									    int precomputedDescendantMatrixElementMaxLeftPLevel, 
									    int precomputedDescendantMatrixElementMaxRightPLevel, 
									    BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
	}
      Tmp *= 1l - (2l * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition == position)
    {
      int TmpPosition = 0;
      while ((TmpPosition < position) && (partition[TmpPosition] >= 0))
	++TmpPosition;
      if (TmpPosition == position)
	{
	  TmpPosition = 1;
	  int TmpPLevel1 = partition[0];
	  bool FlagSorted = true;
	  while (TmpPosition < position)
	    {
	      TmpPLevel1 += partition[TmpPosition];
	      if (partition[TmpPosition - 1] < partition[TmpPosition])
		FlagSorted = false;
	      ++TmpPosition;
	    }
	  if ((TmpPLevel1 <= precomputedDescendantMatrixElementMaxLeftPLevel) && (FlagSorted == true))
	    {
	      int TmpPLevel2 = -partition[TmpPosition];	  
	      FlagSorted = true;
	      ++TmpPosition;
	      while (TmpPosition < partitionLength)
		{
		  TmpPLevel2 -= partition[TmpPosition];
		  if (partition[TmpPosition - 1] < partition[TmpPosition])
		    FlagSorted = false;
		  ++TmpPosition;
		}
	      if ((TmpPLevel2 <= precomputedDescendantMatrixElementMaxRightPLevel) && (FlagSorted == true))
		{
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = 0; TmpPosition < position; ++TmpPosition)
		    {
		      temporaryOccupationNumber[partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel1 - position;
		  int TmpIndex1 = basis[TmpPLevel1]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  for (int k = 0; k <= (this->PLevel + 1); ++k)
		    temporaryOccupationNumber[k] = 0x0ul;	  
		  for (TmpPosition = position; TmpPosition < partitionLength; ++TmpPosition)
		    {
		      temporaryOccupationNumber[-partition[TmpPosition]]++;	      
		    }
		  temporaryOccupationNumber[0] = TmpPLevel2 - partitionLength + position;
		  int TmpIndex2 = basis[TmpPLevel2]->FindStateIndexFromOccupationNumber(temporaryOccupationNumber);
		  LongRational Tmp;
		  precomputedDescendantMatrixElement[TmpPLevel1][TmpPLevel2].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		  return Tmp;
		}
	    }
	}
    }
  if (descendantPosition < position)
    {
      LongRational Tmp(0l);
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  this->ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, 
						       position - 1, centralCharge12, weight1, weight2, weight,
						       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
 	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, 
							     position - 1, centralCharge12, weight1, weight2, weight,
							     precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							     precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, 
						  centralCharge12, weight1, weight2, weight,
						  precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
						  precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  LongRational Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight,
							   precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							   precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
  LongRational Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= this->ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					       centralCharge12, weight1, weight2, weight,
					       precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
					       precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

