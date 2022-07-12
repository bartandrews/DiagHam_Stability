////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
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
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

#include "GeneralTools/FilenameTools.h"


// default constructor 
//

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix()
{
}

// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
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

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool useRational, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion, 
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->BosonicVersion = bosonicVersion;
  this->PLevel = pLevel;
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
    }
  this->TransferMatrixDegeneracy = 5;
  this->NbrCFTSectors = 4;
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction
// architecture = architecture to use for precalculation

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion, 
						   bool useRational, bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion, 
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->BosonicVersion = bosonicVersion;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
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
    }
  this->TransferMatrixDegeneracy = 5;
  this->NbrCFTSectors = 4;
  this->CreateBMatrices(cftDirectory, architecture);
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction

FQHEMPSReadRezayi3Matrix::FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion)
{
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
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
    }
  this->TransferMatrixDegeneracy = 5;
  this->NbrCFTSectors = 4;
  this->LoadMatrices(fileName);
}


// destructor
//

FQHEMPSReadRezayi3Matrix::~FQHEMPSReadRezayi3Matrix()
{
}
  
// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSReadRezayi3Matrix::GetName()
{
  char* TmpName = new char[16];
  sprintf (TmpName, "readrezayi3");
  return TmpName;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSReadRezayi3Matrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 3;
  denominator = 2 + 3 * (this->LaughlinIndex - 1);
}


// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int FQHEMPSReadRezayi3Matrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  if (this->TorusFlag == false)
    {
      nbrFluxQuanta += this->LaughlinIndex;
      nbrFluxQuanta *= 3;
      return (nbrFluxQuanta / (2 + 3 * (this->LaughlinIndex - 1)));
    }
  else
    {
      return ((nbrFluxQuanta * 3) / (2 + 3 * (this->LaughlinIndex - 1)));      
    }
}
// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSReadRezayi3Matrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge (4l, 5l);
  cout << "central charge = " << CentralCharge << endl;
  LongRational CentralCharge12(CentralCharge);
  CentralCharge12 /= 12l;
  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (2l, 3l);
  LongRational WeightW (3l);
  double WeightIdentityNumerical = WeightIdentity.GetNumericalValue();
  double WeightPsiNumerical = WeightPsi.GetNumericalValue();
  double WeightWNumerical = WeightW.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductW = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductW = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi11 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi12 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi21 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi11 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi12 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi21 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisWLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisWRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];
  LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevel + 1];
  double** MultiplicityFactor = new double*[this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi11[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi21[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi12[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi11[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi21[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi12[i] = new LongRationalMatrix[this->PLevel + 1];
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
  
  char* TmpScalarProductIdentityFileName = 0; 
  char* TmpScalarProductPsiFileName = 0;
  char* TmpScalarProductWFileName = 0;
  char* TmpMatrixElementIdentityPsiFileName = 0;
  char* TmpMatrixElementPsiIdentityFileName = 0;
  char* TmpMatrixElementPsiPsiFileName = 0;
  char* TmpMatrixElementWPsiFileName = 0;
  char* TmpMatrixElementPsiWFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductWFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementIdentityPsiFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementPsiIdentityFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementPsiPsiFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementWPsiFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementPsiWFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      RationalScalarProductIdentity[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      RationalScalarProductPsi[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if ((3 + i) <= this->PLevel)
	RationalScalarProductW[3 + i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (i < 3)
	RationalScalarProductW[i] = LongRationalMatrix();
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_readrezayi3_scalarproducts_identity_level_%d.dat", cftDirectory, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_readrezayi3_scalarproducts_psi_level_%d.dat", cftDirectory, i);	  
	      if ((3 + i) <= this->PLevel)
		sprintf (TmpScalarProductWFileName, "%s/cft_readrezayi3_scalarproducts_w_level_%d.dat", cftDirectory, (i + 3));
	    }
	  else
	    {
	      sprintf (TmpScalarProductIdentityFileName, "%s/cft_readrezayi3_num_scalarproducts_identity_level_%d.dat", cftDirectory, i);
	      sprintf (TmpScalarProductPsiFileName, "%s/cft_readrezayi3_num_scalarproducts_psi_level_%d.dat", cftDirectory, i);	  
	      if ((3 + i) <= this->PLevel)
		sprintf (TmpScalarProductWFileName, "%s/cft_readrezayi3_num_scalarproducts_w_level_%d.dat", cftDirectory, (i + 3));
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
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
						     WeightIdentity,
						     RationalScalarProductIdentity,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      RationalScalarProductIdentity[i] = Operation1.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical,
						     WeightIdentityNumerical,
						     ScalarProductIdentity,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      ScalarProductIdentity[i] = Operation1.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductIdentity[i].WriteMatrix(TmpScalarProductIdentityFileName);
		}
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductPsiFileName)))
	{
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
	    }
	  else
	    {
	      ScalarProductPsi[i].ReadMatrix(TmpScalarProductPsiFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
						     WeightPsi,
						     RationalScalarProductPsi,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      RationalScalarProductPsi[i] = Operation2.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12Numerical, 
						     WeightPsiNumerical,
						     ScalarProductPsi,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      ScalarProductPsi[i] = Operation2.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductPsi[i].WriteMatrix(TmpScalarProductPsiFileName);
		}
	    }
	}
      if ((3 + i) <= this->PLevel)
	{
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductWFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalScalarProductW[3 + i].ReadMatrix(TmpScalarProductWFileName);
		}
	      else
		{
		  ScalarProductW[3 + i].ReadMatrix(TmpScalarProductWFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
							 WeightW,
							 RationalScalarProductW + 3,  i - 1);
		  Operation1.ApplyOperation(architecture);
		  RationalScalarProductW[3 + i] = Operation1.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductW[3 + i].WriteMatrix(TmpScalarProductWFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical,
							 WeightWNumerical,
							 ScalarProductW + 3,  i- 1);
		  Operation1.ApplyOperation(architecture);
		  ScalarProductW[3 + i] = Operation1.GetOverlapMatrix();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      ScalarProductW[3 + i].WriteMatrix(TmpScalarProductWFileName);
		    }
		}
	    }
	}
      if (i < 3)
	ScalarProductW[i] = RealSymmetricMatrix();
      
      
      RealSymmetricMatrix TmpMatrix;
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
      RealMatrix TmpBasis(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      double Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      int Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	{
	  if (fabs(TmpDiag(n, n)) < Error)
	    ++Count;
	}
      cout << "nbr of null vectors identity sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisIdentityRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisIdentityLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisIdentityRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix();
	  OrthogonalBasisIdentityRight[i] = RealMatrix();
	}
      
      if ((3 + i) <= this->PLevel)
	{	  
	  if (this->UseRationalFlag == true)
	    {
	      LongRationalMatrix TmpRationalMatrix(RationalScalarProductW[3 + i].GetNbrRow(), RationalScalarProductW[3 + i].GetNbrColumn());
	      for (int k = 0; k < RationalScalarProductW[3 + i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductW[3 + i].GetNbrColumn(); ++l)
		  {
		    TmpRationalMatrix[l][k] = RationalScalarProductW[3 + i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      TmpMatrix = TmpRationalMatrix;
	    }
	  else
	    {
	      TmpMatrix = RealSymmetricMatrix (ScalarProductW[3 + i].GetNbrRow(), ScalarProductW[3 + i].GetNbrColumn());
	      for (int k = 0; k < ScalarProductW[3 + i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductW[3 + i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductW[3 + i].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		    TmpMatrix.SetMatrixElement(k, l, Tmp);
		  }
	    }
	  TmpBasis.SetToIdentity();
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
	  TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
	  double Error = 0.0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      Error = fabs(TmpDiag(n, n));
	  Error *= 1e-14;
	  if (Error < 1e-14)
	    Error = 1e-14;
	  int Count  = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    {
	      if (fabs(TmpDiag(n, n)) < Error)
		++Count;
	    }
	  cout << "nbr of null vectors W_-3 sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      OrthogonalBasisWRight[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisWLeft[3 + i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisWRight[3 + i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
		  }
	    }
	  else
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix();
	      OrthogonalBasisWRight[3 + i] = RealMatrix();
	    }
	}
      if (i < 3)
	{
	  OrthogonalBasisWLeft[i] = RealMatrix();
	  OrthogonalBasisWRight[i] = RealMatrix();
	}
      
      
      if (this->UseRationalFlag == true)
 	{
	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductPsi[i].GetNbrRow(), RationalScalarProductPsi[i].GetNbrColumn());
	  for (int k = 0; k < RationalScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductPsi[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
       else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductPsi[i].GetNbrRow(), ScalarProductPsi[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductPsi[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		TmpMatrix.SetMatrixElement(k, l, Tmp);
	      }
	}

       TmpBasis.SetToIdentity();
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	{
	  if (fabs(TmpDiag(n, n)) < Error)
	    ++Count;
	}
      cout << "nbr of null vectors Psi sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;

      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisPsiRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisPsiLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisPsiRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix();
	  OrthogonalBasisPsiRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  for (int i = 0; i <= this->PLevel; ++i)
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
      if (this->UseRationalFlag == true)
 	{
	  for (int k = 0; k < RationalScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		RationalScalarProductPsi[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  ScalarProductPsi[i] = RationalScalarProductPsi[i];
 	}
       else
	{
	  for (int k = 0; k < ScalarProductPsi[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductPsi[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductPsi[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductPsi[i].SetMatrixElement(k, l, Tmp);
	      }
	}
      if ((3 + i) <= this->PLevel)
	{	  
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < RationalScalarProductW[3 + i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductW[3 + i].GetNbrColumn(); ++l)
		  {
		    RationalScalarProductW[3 + i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      ScalarProductW[3 + i] = RationalScalarProductW[3 + i];
	    }
	  else
	    {
	      for (int k = 0; k < ScalarProductW[3 + i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductW[3 + i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductW[3 + i].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		    ScalarProductW[3 + i].SetMatrixElement(k, l, Tmp);
		  }
	    }
	}
      if (i < 3)
	ScalarProductW[i] = RealSymmetricMatrix();
    }


  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [4];
  this->NeutralSectorDimension[0] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[2] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[3] = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisIdentityLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[2][i] = OrthogonalBasisPsiLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[3][i] = OrthogonalBasisWLeft[i].GetNbrColumn();
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }
  
  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  double WeightNumerical = Weight.GetNumericalValue();
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  sprintf (TmpMatrixElementIdentityPsiFileName, "%s/cft_readrezayi3_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementPsiIdentityFileName, "%s/cft_readrezayi3_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementPsiPsiFileName, "%s/cft_readrezayi3_matrixelement_psipsi_level_%d_%d.dat", cftDirectory, i, j);
		  if ((3 + i) <= this->PLevel)
		    sprintf (TmpMatrixElementWPsiFileName, "%s/cft_readrezayi3_matrixelement_wpsi_level_%d_%d.dat", cftDirectory, (3 + i), j);
		  if ((3 + j) <= this->PLevel)
		    sprintf (TmpMatrixElementPsiWFileName, "%s/cft_readrezayi3_matrixelement_psiw_level_%d_%d.dat", cftDirectory, i, (3 + j));
		}
	      else
		{
		  sprintf (TmpMatrixElementIdentityPsiFileName, "%s/cft_readrezayi3_num_matrixelement_identitypsi_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementPsiIdentityFileName, "%s/cft_readrezayi3_num_matrixelement_psiidentity_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementPsiPsiFileName, "%s/cft_readrezayi3_num_matrixelement_psipsi_level_%d_%d.dat", cftDirectory, i, j);
		  if ((3 + i) <= this->PLevel)
		    sprintf (TmpMatrixElementWPsiFileName, "%s/cft_readrezayi3_num_matrixelement_wpsi_level_%d_%d.dat", cftDirectory, (3 + i), j);
		  if ((3 + j) <= this->PLevel)
		    sprintf (TmpMatrixElementPsiWFileName, "%s/cft_readrezayi3_num_matrixelement_psiw_level_%d_%d.dat", cftDirectory, i, (3 + j));
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementIdentityPsiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi01[i][j].ReadMatrix(TmpMatrixElementIdentityPsiFileName);
		}
	      else
		{
		  MatrixPsi01[i][j].ReadMatrix(TmpMatrixElementIdentityPsiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightIdentity, WeightPsi, Weight,
							 RationalMatrixPsi01,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi01[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpMatrixElementIdentityPsiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightIdentityNumerical, WeightPsiNumerical, WeightNumerical,
							 MatrixPsi01,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi01[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi01[i][j].WriteMatrix(TmpMatrixElementIdentityPsiFileName);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementPsiIdentityFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi10[i][j].ReadMatrix(TmpMatrixElementPsiIdentityFileName);
		}
	      else
		{
		  MatrixPsi10[i][j].ReadMatrix(TmpMatrixElementPsiIdentityFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightPsi, WeightIdentity, Weight,
							 RationalMatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi10[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpMatrixElementPsiIdentityFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightPsiNumerical, WeightIdentityNumerical, WeightNumerical,
							 MatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi10[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi10[i][j].WriteMatrix(TmpMatrixElementPsiIdentityFileName);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementPsiPsiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi11[i][j].ReadMatrix(TmpMatrixElementPsiPsiFileName);
		}
	      else
		{
		  MatrixPsi11[i][j].ReadMatrix(TmpMatrixElementPsiPsiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightPsi, WeightPsi, Weight,
							 RationalMatrixPsi11,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi11[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi11[i][j].WriteMatrix(TmpMatrixElementPsiPsiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightPsiNumerical, WeightPsiNumerical, WeightNumerical,
							 MatrixPsi11,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi11[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi11[i][j].WriteMatrix(TmpMatrixElementPsiPsiFileName);
		    }
		}
	    }
	  if ((3 + i) <= this->PLevel)
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementWPsiFileName)))
		{
		  if (this->UseRationalFlag == true)
		    {
		      RationalMatrixPsi21[i][j].ReadMatrix(TmpMatrixElementWPsiFileName);
		    }
		  else
		    {
		      MatrixPsi21[i][j].ReadMatrix(TmpMatrixElementWPsiFileName);
		    }
		}
	      else
		{
		  if (this->UseRationalFlag == true)
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightW, WeightPsi, Weight,
							     RationalMatrixPsi21,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi21[i][j] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi21[i][j].WriteMatrix(TmpMatrixElementWPsiFileName);
			}
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							     WeightWNumerical, WeightPsiNumerical, WeightNumerical,
							     MatrixPsi21,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      MatrixPsi21[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  MatrixPsi21[i][j].WriteMatrix(TmpMatrixElementWPsiFileName);
			}
		    }
		}
	    }
	  if ((3 + j) <= this->PLevel)
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementPsiWFileName)))
		{
		  if (this->UseRationalFlag == true)
		    {
		      RationalMatrixPsi12[i][j].ReadMatrix(TmpMatrixElementPsiWFileName);
		    }
		  else
		    {
		      MatrixPsi12[i][j].ReadMatrix(TmpMatrixElementPsiWFileName);
		    }
		}
	      else
		{
		  if (this->UseRationalFlag == true)
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightPsi, WeightW, Weight,
							     RationalMatrixPsi12,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi12[i][j] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi12[i][j].WriteMatrix(TmpMatrixElementPsiWFileName);
			}
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							     WeightPsiNumerical, WeightWNumerical, WeightNumerical,
							     MatrixPsi12,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      MatrixPsi12[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  MatrixPsi12[i][j].WriteMatrix(TmpMatrixElementPsiWFileName);
			}
		    }
		}	      
	    }
	}
    }

  if (this->UseRationalFlag == false)
    {
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  for (int i = this->PLevel - 3; i >= 0; --i)
	    {
	      MatrixPsi21[3 + i][j] = MatrixPsi21[i][j];
	    }
	  if (this->PLevel >= 2)
	    MatrixPsi21[2][j] = RealMatrix();
	  if (this->PLevel >= 1)
	    MatrixPsi21[1][j] = RealMatrix();
	  MatrixPsi21[0][j] = RealMatrix();
	}
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int j = this->PLevel - 3; j >= 0; --j)
	    {
	      MatrixPsi12[i][3 + j] = MatrixPsi12[i][j];
	    }
	  if (this->PLevel >= 2)
	    MatrixPsi12[i][2] = RealMatrix();
	  if (this->PLevel >= 1)
	    MatrixPsi12[i][1] = RealMatrix();
	  MatrixPsi12[i][0] = RealMatrix();
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
		    RationalMatrixPsi01[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		  }
	      for (int k = 0; k < RationalMatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi10[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		  }
	      for (int k = 0; k < RationalMatrixPsi11[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < RationalMatrixPsi11[i][j].GetNbrColumn(); ++l)
		  {
		    RationalMatrixPsi11[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		  }
	      if ((3 + j) <= this->PLevel)
		{	
		  for (int k = 0; k < RationalMatrixPsi11[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < RationalMatrixPsi11[i][j].GetNbrColumn(); ++l)
		      {
			RationalMatrixPsi12[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		      }
		}
	      if ((3 + i) <= this->PLevel)
		{	  
		  for (int k = 0; k < RationalMatrixPsi11[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < RationalMatrixPsi11[i][j].GetNbrColumn(); ++l)
		      {
			RationalMatrixPsi21[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		      }
		}
	      MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	      MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
	      MatrixPsi11[i][j] = RationalMatrixPsi11[i][j];	  
	      if ((3 + j) <= this->PLevel)
		{	
		  MatrixPsi12[i][3 + j] = RationalMatrixPsi12[i][j];
		}
	      if ((3 + i) <= this->PLevel)
		{	  
		  MatrixPsi21[3 + i][j] = RationalMatrixPsi21[i][j];
		}
	    }
	  else
	    {
	      for (int k = 0; k < MatrixPsi01[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi01[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi01[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
		    MatrixPsi01[i][j].SetMatrixElement(k, l, Tmp);
		  }
	      for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi10[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
		    MatrixPsi10[i][j].SetMatrixElement(k, l, Tmp);
		  }
	      for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    MatrixPsi11[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
		    MatrixPsi11[i][j].SetMatrixElement(k, l, Tmp);
		  }
	      if ((3 + j) <= this->PLevel)
		{	
		  for (int k = 0; k < MatrixPsi01[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < MatrixPsi01[i][j].GetNbrColumn(); ++l)
		      {
			double Tmp;
			MatrixPsi12[i][3 + j].GetMatrixElement(k, l, Tmp);
			Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
			MatrixPsi12[i][3 + j].SetMatrixElement(k, l, Tmp);
		  }
		}
	      if ((3 + i) <= this->PLevel)
		{	  
		  for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		      {
			double Tmp;
			MatrixPsi21[3 + i][j].GetMatrixElement(k, l, Tmp);
			Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
			MatrixPsi21[3 + i][j].SetMatrixElement(k, l, Tmp);
		      }
		}
	    }
	  MatrixPsi11[i][j] *= 2.0 / sqrt(3.0);	  
	  if ((3 + j) <= this->PLevel)
	    {	
	      MatrixPsi12[i][3 + j] *= -sqrt(26.0) / 9.0;
	    }
	  if ((3 + i) <= this->PLevel)
	    {	  
	      MatrixPsi21[3 + i][j] *= sqrt(26.0) / 9.0;
	    }
	}
    }

  cout << "building B matrices" << endl;

  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  if (this->BosonicVersion == false)
    {
      QValue = 5;
      QValueDenominator = 3;
      this->NbrNValue = QValueDenominator * (2 * this->PLevel) + 4 + 2 + 1;
      NValueShift = QValueDenominator * this->PLevel;
    }
  else
    {
      QValue = 2 +  3 * (this->LaughlinIndex - 1);
      QValueDenominator = 3;
      this->NbrNValue = QValueDenominator * (2 * this->PLevel) + QValue + 4 + 2 + 1;
      NValueShift = QValueDenominator * this->PLevel;
    }


  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;
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
	      // identity 
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
			}
		    }
		}
	      // Psi_{+1}
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
			}
		    }
		}
	      // Psi_{-1}
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][2] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][2]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 2, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWLeft = OrthogonalBasisWLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWRight = OrthogonalBasisWRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductW = ScalarProductW[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // W
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][3] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][3]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 3, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
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
	      // identity 
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue * QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue * QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	      // Psi_{+1}
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	      // Psi_{-1}
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][2] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][2]; ++j)
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
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 2, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 2, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWLeft = OrthogonalBasisWLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWRight = OrthogonalBasisWRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductW = ScalarProductW[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // W
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][3] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][3]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductW(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisWRight(NeutralIndex4, NeutralIndex2);
				}
			      Tmp += TmpOrthogonalBasisWLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 3, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 3, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	    }
	}
    }

  // B^[1]  matrix evaluation
  if (this->BosonicVersion == false)
    {
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
		      RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				  
				}
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
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
			      
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 2, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= (i - 3); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 3, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 3); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
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
		      RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
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
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
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
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightPsiNumerical + WeightPsiNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 2, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= (i - 3); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisW1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 3, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 3); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisW2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightIdentityNumerical + WeightPsiNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 3, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
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
  else
    {
      // creation of the B matrices compatible with bosons

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
		      RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				  
				}
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
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
			      
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 2, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= (i - 3); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 3, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 3); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
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
      
      SparseRealMatrix V0Matrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
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
		      RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
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
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
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
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 2, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= (i - 3); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisW1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 3, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}
      
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 3); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		      RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		      RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisW2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
								    this->Get2RMatrixIndexV2(j, 3, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			    }
			}
		    }	      
		}
	    }
	}

      int TmpNbrBMatrices = 0;
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
  int TmpNbrBMatrices = 0;
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      if (BMatrices[i].GetNbrRow() > 0)
	this->RealBMatrices[TmpNbrBMatrices++] = BMatrices[i];
    }
  this->NbrBMatrices = TmpNbrBMatrices;

  delete[] BMatrices;

  delete[] ScalarProductIdentity;
  delete[] ScalarProductPsi;
  delete[] ScalarProductW;
  delete[] RationalScalarProductIdentity;
  delete[] RationalScalarProductPsi;
  delete[] RationalScalarProductW;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete[] MatrixPsi21[i];
      delete[] MatrixPsi12[i];
      delete[] MatrixPsi11[i];
      delete[] RationalMatrixPsi01[i];
      delete[] RationalMatrixPsi10[i];
      delete[] RationalMatrixPsi21[i];
      delete[] RationalMatrixPsi12[i];
      delete[] RationalMatrixPsi11[i];
      delete U1BosonBasis[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] MatrixPsi21;
  delete[] MatrixPsi12;
  delete[] MatrixPsi11;
  delete[] RationalMatrixPsi01;
  delete[] RationalMatrixPsi10;
  delete[] RationalMatrixPsi21;
  delete[] RationalMatrixPsi12;
  delete[] RationalMatrixPsi11;
  delete[] OrthogonalBasisIdentityLeft;
  delete[] OrthogonalBasisPsiLeft;
  delete[] OrthogonalBasisWLeft;
  delete[] OrthogonalBasisIdentityRight;
  delete[] OrthogonalBasisPsiRight;
  delete[] OrthogonalBasisWRight;
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSReadRezayi3Matrix::GetBondIndexRange(int pLevel, int qValue)
{
  if ((pLevel < 0) || (pLevel > this->PLevel))
    return 0;
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][0]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][0]))
    {
      int Tmp = this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][3]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][3]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][3][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][3]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][2]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][2]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][2][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][2]];
      return Tmp;
    }
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][1]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][1]))
    {
      int Tmp = this->NbrIndexPerPLevelCFTSectorQValue[pLevel][1][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][1]];
      if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][2]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][2]))
	Tmp += this->NbrIndexPerPLevelCFTSectorQValue[pLevel][2][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][2]];
      return Tmp;
    }
  if ((qValue >= this->NInitialValuePerPLevelCFTSector[pLevel][2]) && (qValue <= this->NLastValuePerPLevelCFTSector[pLevel][2]))
    return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][2][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][2]];
  return 0;
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int FQHEMPSReadRezayi3Matrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  if ((pLevel < 0) || (pLevel > this->PLevel) || (qValue < this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]) || 
      (qValue > this->NLastValuePerPLevelCFTSector[pLevel][cftSector]) || (cftSector > 2) ||  (cftSector < 0))
    return 0;
  if  (cftSector == 0)
    return (this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]
	    + this->NbrIndexPerPLevelCFTSectorQValue[pLevel][3][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][3]]);
  else
    return this->NbrIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]];
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int FQHEMPSReadRezayi3Matrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
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

int FQHEMPSReadRezayi3Matrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{
  if (cftSector != 0)
    return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]]  + localIndex);
  else
    {
      if (localIndex < this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]])
	return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]]  + localIndex);
      else
	return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][3][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][3]] + localIndex - this->NbrIndexPerPLevelCFTSectorQValue[pLevel][0][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][0]]);
    }
}


// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSReadRezayi3Matrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][0];
  if (this->NInitialValuePerPLevelCFTSector[pLevel][1] < minQ)
    minQ = this->NInitialValuePerPLevelCFTSector[pLevel][1];
  if (this->NLastValuePerPLevelCFTSector[pLevel][1] > maxQ)
    maxQ = this->NLastValuePerPLevelCFTSector[pLevel][1];  
  if (this->NInitialValuePerPLevelCFTSector[pLevel][2] < minQ)
    minQ = this->NInitialValuePerPLevelCFTSector[pLevel][2];
  if (this->NLastValuePerPLevelCFTSector[pLevel][2] > maxQ)
    maxQ = this->NLastValuePerPLevelCFTSector[pLevel][2];  
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSReadRezayi3Matrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSReadRezayi3Matrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      rowIndex = 3 * (this->PLevel + 1) - MinQ;
      columnIndex = rowIndex;
    }
  else
    {
      if (this->BosonicVersion == false)
	{
	  rowIndex = 3 * (this->PLevel + 2) - MinQ;
	  columnIndex = 3 * this->PLevel - MinQ;
	}
      else
	{
	  rowIndex = 3 * (this->PLevel + this->LaughlinIndex) - MinQ;
	  columnIndex = 3 * this->PLevel - MinQ;
	}
    }
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSReadRezayi3Matrix::GetQValueCFTSectorShift(int cftSector)
{
  switch (cftSector)
    {
    case 0:
      return 0;
    case 1:
      return 0;
    case 2:
      return 0;
    default:
      return 0;
    }
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSReadRezayi3Matrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
{
  if (this->UniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }
  int TmpMinQ = this->NbrNValue - 1;
  int TmpMaxQ = 0;    
  int NValueShift = this->PLevel;
  int QValue;
  int QValueDenominator;
  if (this->BosonicVersion == false)
    {
      QValue = 5;
      QValueDenominator = 3;
    }
  else
    {
      QValue = 2 + 3 * (this->LaughlinIndex - 1);
      QValueDenominator = 3;
    }
  if (this->BosonicVersion == false)
    {
      if ((cftSector == 0) || (cftSector == 3))
	{
	  for (int Q = 0; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - QValueDenominator);
		  TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - QValueDenominator);
		      TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - QValueDenominator);
			  TmpP += QPrime / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= QPrime / QValueDenominator - NValueShift;
		  QPrime += (QValue - QValueDenominator);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
		      QPrime += (QValue - QValueDenominator);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
			  QPrime += (QValue - QValueDenominator);
			}
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
      if (cftSector == 1) 
	{
	  for (int Q = 2; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - QValueDenominator);
		  TmpP += (QPrime) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - QValueDenominator);
		      TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - QValueDenominator);
			  TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
		  QPrime += (QValue - QValueDenominator);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
		      QPrime += (QValue - QValueDenominator);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime) / QValueDenominator - NValueShift;
			  QPrime += (QValue - QValueDenominator);
			}
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
      if (cftSector == 2) 
	{
	  for (int Q = 1; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue - QValueDenominator);
		  TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue - QValueDenominator);
		      TmpP += (QPrime) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue - QValueDenominator);
			  TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
		  QPrime += (QValue - QValueDenominator);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime) / QValueDenominator - NValueShift;
		      QPrime += (QValue - QValueDenominator);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
			  QPrime += (QValue - QValueDenominator);
			}
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
  else
    {
      if ((cftSector == 0) || (cftSector == 3))
	{
	  for (int Q = 0; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue);
		  TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += QPrime / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= QPrime / QValueDenominator - NValueShift;
		  QPrime += (QValue);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
			  QPrime += (QValue);
			}
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
      if (cftSector == 1) 
	{
	  for (int Q = 2; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue);
		  TmpP += (QPrime) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
		  QPrime += (QValue);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime) / QValueDenominator - NValueShift;
			  QPrime += (QValue);
			}
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
      if (cftSector == 2) 
	{
	  for (int Q = 1; Q < this->NbrNValue; Q += 3)
	    {
	      int QPrime = Q;
	      int TmpP = 0;
	      int TmpMaxP = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP)
		    TmpMaxP = TmpP;	    
		  QPrime -= (QValue);
		  TmpP += (QPrime - 2) / QValueDenominator - NValueShift;
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP)
			TmpMaxP = TmpP;	    
		      QPrime -= (QValue);
		      TmpP += (QPrime) / QValueDenominator - NValueShift;
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP)
			    TmpMaxP = TmpP;	    
			  QPrime -= (QValue);
			  TmpP += (QPrime - 4) / QValueDenominator - NValueShift;
			}
		    }
		}
	      QPrime = Q;
	      TmpP = 0;
	      int TmpMaxP2 = 0;
	      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		{
		  if (TmpP > TmpMaxP2)
		    TmpMaxP2 = TmpP;	    
		  TmpP -= (QPrime - 4) / QValueDenominator - NValueShift;
		  QPrime += (QValue);
		  if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= (QPrime) / QValueDenominator - NValueShift;
		      QPrime += (QValue);
		      if ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= (QPrime - 2) / QValueDenominator - NValueShift;
			  QPrime += (QValue);
			}
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

// get the matrix that into account the Jordan Wigner string on the torus geometry
//
// nbrFermions = number of fermions in the system
// return value = corresponding matrix

SparseRealMatrix FQHEMPSReadRezayi3Matrix::GetTorusStringMatrix(int nbrFermions)
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
  SparseRealMatrix StringMatrix (TmpDimension, TmpDimension, TmpNbrElementPerRow);

  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
 	{
	  int MinQValue;
	  int MaxQValue;
	  this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	  for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	    {
	      double Coefficient = 1.0;
	      //	      if (((CurrentQValue + this->GetQValueCFTSectorShift(CurrentCFTSector)) & 1) == 0)
	      if ((CurrentQValue & 1) == 0)
		{
		  Coefficient = -1.0;
		}
	      int NbrIndices = this->GetBondIndexRange(CurrentPLevel, CurrentQValue, CurrentCFTSector);
	      for (int i = 0; i < NbrIndices; ++i)
		{	      
		  int TmpIndex = this->GetBondIndexWithFixedChargePLevelCFTSector(i, CurrentPLevel, CurrentQValue, CurrentCFTSector);
		  StringMatrix.SetMatrixElement(TmpIndex, TmpIndex, Coefficient);
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

int* FQHEMPSReadRezayi3Matrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
{
  nbrIndices = 0;
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
 	{
 	  int MinQValue;
	  int MaxQValue;
	  this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	  for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	    {
	      if (((CurrentQValue + this->GetQValueCFTSectorShift(CurrentCFTSector)) % (2 + 3 * (this->LaughlinIndex - 1))) == topologicalSector)
		{
		  nbrIndices += this->GetBondIndexRange(CurrentPLevel, CurrentQValue, CurrentCFTSector);
		}
	    }
	}
    }
  int* TmpIndices =  new int [nbrIndices];
  nbrIndices = 0;
  for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
    {
      for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
 	{
 	  int MinQValue;
	  int MaxQValue;
	  this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	  for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
	    {
	      if (((CurrentQValue + this->GetQValueCFTSectorShift(CurrentCFTSector)) % (2 + 3 * (this->LaughlinIndex - 1))) == topologicalSector)
		{
		  int TmpBondIndexRange = this->GetBondIndexRange(CurrentPLevel, CurrentQValue, CurrentCFTSector);
		  for (int i = 0; i < TmpBondIndexRange; ++i)
		    {
		      TmpIndices[nbrIndices] = this->GetBondIndexWithFixedChargePLevelCFTSector(i, CurrentPLevel, CurrentQValue, CurrentCFTSector);
		      ++nbrIndices;
		    }
		}
	    }
	}
    }
  return TmpIndices;
}
  
