////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
//                           in its quasihole sector                          //
//                                                                            //
//                        last modification : 22/02/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3QuasiholeSectorMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

#include "GeneralTools/FilenameTools.h"



// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
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

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool useRational, 
										 bool trimChargeIndices, bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->BosonicVersion = bosonicVersion;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
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
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
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

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion, 
										 bool useRational, 
										 bool trimChargeIndices, bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->UseRationalFlag = useRational;
  this->BosonicVersion = bosonicVersion;
  this->UniformChargeIndexRange = !trimChargeIndices;
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

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa, 
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
  this->LoadMatrices(fileName);
}


// destructor
//

FQHEMPSReadRezayi3QuasiholeSectorMatrix::~FQHEMPSReadRezayi3QuasiholeSectorMatrix()
{
}
  
// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSReadRezayi3QuasiholeSectorMatrix::GetName()
{
  char* TmpName = new char[24];
  sprintf (TmpName, "readrezayi3_qh");
  return TmpName;
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSReadRezayi3QuasiholeSectorMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge (4l, 5l);
  cout << "central charge = " << CentralCharge << endl;
  LongRational CentralCharge12(CentralCharge);
  CentralCharge12 /= 12l;
  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  LongRational WeightPsi (2l, 3l);
  LongRational WeightSigma (1l, 15l);
  LongRational WeightPhi (7l, 5l);
  LongRational WeightEpsilon (2l, 5l);
  double WeightSigmaNumerical = WeightSigma.GetNumericalValue();
  double WeightPhiNumerical = WeightPhi.GetNumericalValue();
  double WeightEpsilonNumerical = WeightEpsilon.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductSigma = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPhi = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductEpsilon = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductSigma = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPhi = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductEpsilon = new LongRationalMatrix[this->PLevel + 1];
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
  RealMatrix* OrthogonalBasisSigmaLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisEpsilonLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisEpsilonRight = new RealMatrix[this->PLevel + 1];
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
  
  char* TmpScalarProductSigmaFileName = 0; 
  char* TmpScalarProductEpsilonFileName = 0;
  char* TmpScalarProductPhiFileName = 0;
  char* TmpMatrixElementSigmaEpsilonFileName = 0;
  char* TmpMatrixElementEpsilonSigmaFileName = 0;
  char* TmpMatrixElementSigmaSigmaFileName = 0;
  char* TmpMatrixElementPhiSigmaFileName = 0;
  char* TmpMatrixElementSigmaPhiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductEpsilonFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPhiFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaEpsilonFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementEpsilonSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementPhiSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaPhiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_readrezayi3_scalarproducts_sigma_level_%d.dat", cftDirectory, i);
	      sprintf (TmpScalarProductEpsilonFileName, "%s/cft_readrezayi3_scalarproducts_epsilon_level_%d.dat", cftDirectory, i);	  
	      if ((1 + i) <= this->PLevel)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_readrezayi3_scalarproducts_phi_level_%d.dat", cftDirectory, (i + 1));
	    }
	  else
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_readrezayi3_num_scalarproducts_sigma_level_%d.dat", cftDirectory, i);
	      sprintf (TmpScalarProductEpsilonFileName, "%s/cft_readrezayi3_num_scalarproducts_epsilon_level_%d.dat", cftDirectory, i);	  
	      if ((1 + i) <= this->PLevel)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_readrezayi3_num_scalarproducts_phi_level_%d.dat", cftDirectory, (i + 1));
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductSigmaFileName)))
	{
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
	    }
	  else
	    {
	      ScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
						     WeightSigma,
						     RationalScalarProductSigma,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      RationalScalarProductSigma[i] = Operation1.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical,
						     WeightSigmaNumerical,
						     ScalarProductSigma,  i- 1);
	      Operation1.ApplyOperation(architecture);
	      ScalarProductSigma[i] = Operation1.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		}
	    }
	}
      if ((cftDirectory != 0) && (IsFile(TmpScalarProductEpsilonFileName)))
	{
	  if (this->UseRationalFlag == true)
	    {
	      RationalScalarProductEpsilon[i].ReadMatrix(TmpScalarProductEpsilonFileName);
	    }
	  else
	    {
	      ScalarProductEpsilon[i].ReadMatrix(TmpScalarProductEpsilonFileName);
	    }
	}
      else
	{
	  if (this->UseRationalFlag == true)
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
						     WeightEpsilon,
						     RationalScalarProductEpsilon,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      RationalScalarProductEpsilon[i] = Operation2.GetRationalMatrixElements();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  RationalScalarProductEpsilon[i].WriteMatrix(TmpScalarProductEpsilonFileName);
		}
	    }
	  else
	    {
	      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12Numerical, 
						     WeightEpsilonNumerical,
						     ScalarProductEpsilon,  i - 1);
	      Operation2.ApplyOperation(architecture);
	      ScalarProductEpsilon[i] = Operation2.GetOverlapMatrix();
	      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		{
		  ScalarProductEpsilon[i].WriteMatrix(TmpScalarProductEpsilonFileName);
		}
	    }
	}
      if ((1 + i) <= this->PLevel)
	{
	  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPhiFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalScalarProductPhi[1 + i].ReadMatrix(TmpScalarProductPhiFileName);
		}
	      else
		{
		  ScalarProductPhi[1 + i].ReadMatrix(TmpScalarProductPhiFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
							 WeightPhi,
							 RationalScalarProductPhi + 1,  i - 1);
		  Operation1.ApplyOperation(architecture);
		  RationalScalarProductPhi[1 + i] = Operation1.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductPhi[1 + i].WriteMatrix(TmpScalarProductPhiFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12Numerical,
							 WeightPhiNumerical,
							 ScalarProductPhi + 1,  i- 1);
		  Operation1.ApplyOperation(architecture);
		  ScalarProductPhi[1 + i] = Operation1.GetOverlapMatrix();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      ScalarProductPhi[1 + i].WriteMatrix(TmpScalarProductPhiFileName);
		    }
		}
	    }
	}

      RealSymmetricMatrix TmpMatrix;
      if (this->UseRationalFlag == true)
 	{
 	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductSigma[i].GetNbrRow(), RationalScalarProductSigma[i].GetNbrColumn());
 	  for (int k = 0; k < RationalScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductSigma[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
      else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductSigma[i].GetNbrRow(), ScalarProductSigma[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductSigma[i].GetMatrixElement(k, l, Tmp);
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
      cout << "nbr of null vectors sigma sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisSigmaRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisSigmaLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisSigmaRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix();
	  OrthogonalBasisSigmaRight[i] = RealMatrix();
	}
      
      if ((1 + i) <= this->PLevel)
	{	  
	  if (this->UseRationalFlag == true)
	    {
	      LongRationalMatrix TmpRationalMatrix(RationalScalarProductPhi[1 + i].GetNbrRow(), RationalScalarProductPhi[1 + i].GetNbrColumn());
	      for (int k = 0; k < RationalScalarProductPhi[1 + i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductPhi[1 + i].GetNbrColumn(); ++l)
		  {
		    TmpRationalMatrix[l][k] = RationalScalarProductPhi[1 + i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      TmpMatrix = TmpRationalMatrix;
	    }
	  else
	    {
	      TmpMatrix = RealSymmetricMatrix (ScalarProductPhi[1 + i].GetNbrRow(), ScalarProductPhi[1 + i].GetNbrColumn());
	      for (int k = 0; k < ScalarProductPhi[1 + i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductPhi[1 + i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductPhi[1 + i].GetMatrixElement(k, l, Tmp);
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
	  cout << "nbr of null vectors phi sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisPhiLeft[1 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      OrthogonalBasisPhiRight[1 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisPhiLeft[1 + i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisPhiRight[1 + i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisPhiLeft[1 + i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisPhiRight[1 + i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisPhiLeft[1 + i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisPhiRight[1 + i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
		  }
	    }
	  else
	    {
	      OrthogonalBasisPhiLeft[1 + i] = RealMatrix();
	      OrthogonalBasisPhiRight[1 + i] = RealMatrix();
	    }
	}
      if (i < 1)
	{
	  OrthogonalBasisPhiLeft[i] = RealMatrix();
	  OrthogonalBasisPhiRight[i] = RealMatrix();
	}
      
      
      if (this->UseRationalFlag == true)
 	{
	  LongRationalMatrix TmpRationalMatrix(RationalScalarProductEpsilon[i].GetNbrRow(), RationalScalarProductEpsilon[i].GetNbrColumn());
	  for (int k = 0; k < RationalScalarProductEpsilon[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductEpsilon[i].GetNbrColumn(); ++l)
	      {
		TmpRationalMatrix[l][k] = RationalScalarProductEpsilon[i][l][k] * (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  TmpMatrix = TmpRationalMatrix;
 	}
       else
	{
	  TmpMatrix = RealSymmetricMatrix (ScalarProductEpsilon[i].GetNbrRow(), ScalarProductEpsilon[i].GetNbrColumn());
	  for (int k = 0; k < ScalarProductEpsilon[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductEpsilon[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductEpsilon[i].GetMatrixElement(k, l, Tmp);
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
      cout << "nbr of null vectors Epsilon sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;

      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisEpsilonLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisEpsilonRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisEpsilonLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisEpsilonRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisEpsilonLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisEpsilonRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisEpsilonLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisEpsilonRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisEpsilonLeft[i] = RealMatrix();
	  OrthogonalBasisEpsilonRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      if (this->UseRationalFlag == true)
 	{
 	  for (int k = 0; k < RationalScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		RationalScalarProductSigma[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  ScalarProductSigma[i] = RationalScalarProductSigma[i];
 	}
       else
	{
	  for (int k = 0; k < ScalarProductSigma[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductSigma[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductSigma[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductSigma[i].SetMatrixElement(k, l, Tmp);
	      }
	}
      if (this->UseRationalFlag == true)
 	{
	  for (int k = 0; k < RationalScalarProductEpsilon[i].GetNbrRow(); ++k)
	    for (int l = 0; l < RationalScalarProductEpsilon[i].GetNbrColumn(); ++l)
	      {
		RationalScalarProductEpsilon[i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
	      }
 	  ScalarProductEpsilon[i] = RationalScalarProductEpsilon[i];
 	}
       else
	{
	  for (int k = 0; k < ScalarProductEpsilon[i].GetNbrRow(); ++k)
	    for (int l = k; l < ScalarProductEpsilon[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		ScalarProductEpsilon[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		ScalarProductEpsilon[i].SetMatrixElement(k, l, Tmp);
	      }
	}
      if ((1 + i) <= this->PLevel)
	{	  
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < RationalScalarProductPhi[1 + i].GetNbrRow(); ++k)
		for (int l = 0; l < RationalScalarProductPhi[1 + i].GetNbrColumn(); ++l)
		  {
		    RationalScalarProductPhi[1 + i][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[i][l]);
		  }
	      ScalarProductPhi[1 + i] = RationalScalarProductPhi[1 + i];
	    }
	  else
	    {
	      for (int k = 0; k < ScalarProductPhi[1 + i].GetNbrRow(); ++k)
		for (int l = k; l < ScalarProductPhi[1 + i].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    ScalarProductPhi[1 + i].GetMatrixElement(k, l, Tmp);
		    Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[i][l]);
		    ScalarProductPhi[1 + i].SetMatrixElement(k, l, Tmp);
		  }
	    }
	}
      if (i < 1)
	ScalarProductPhi[i] = RealSymmetricMatrix();
    }

  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [4];
  this->NeutralSectorDimension[0] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[1] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[2] = new int [this->PLevel + 1];
  this->NeutralSectorDimension[3] = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisEpsilonLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[1][i] = OrthogonalBasisSigmaLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[2][i] = OrthogonalBasisSigmaLeft[i].GetNbrColumn();
      this->NeutralSectorDimension[3][i] = OrthogonalBasisPhiLeft[i].GetNbrColumn();
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
		  sprintf (TmpMatrixElementSigmaEpsilonFileName, "%s/cft_readrezayi3_matrixelement_sigmaepsilon_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementEpsilonSigmaFileName, "%s/cft_readrezayi3_matrixelement_epsilonsigma_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementSigmaSigmaFileName, "%s/cft_readrezayi3_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, i, j);
		  if ((1 + i) <= this->PLevel)
		    sprintf (TmpMatrixElementPhiSigmaFileName, "%s/cft_readrezayi3_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, (1 + i), j);
		  if ((1 + j) <= this->PLevel)
		    sprintf (TmpMatrixElementSigmaPhiFileName, "%s/cft_readrezayi3_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, i, (1 + j));
		}
	      else
		{
		  sprintf (TmpMatrixElementSigmaEpsilonFileName, "%s/cft_readrezayi3_num_matrixelement_sigmaepsilon_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementEpsilonSigmaFileName, "%s/cft_readrezayi3_num_matrixelement_epsilonsigma_level_%d_%d.dat", cftDirectory, i, j);
		  sprintf (TmpMatrixElementSigmaSigmaFileName, "%s/cft_readrezayi3_num_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, i, j);
		  if ((1 + i) <= this->PLevel)
		    sprintf (TmpMatrixElementPhiSigmaFileName, "%s/cft_readrezayi3_num_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, (1 + i), j);
		  if ((1 + j) <= this->PLevel)
		    sprintf (TmpMatrixElementSigmaPhiFileName, "%s/cft_readrezayi3_num_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, i, (1 + j));
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementEpsilonSigmaFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi01[i][j].ReadMatrix(TmpMatrixElementEpsilonSigmaFileName);
		}
	      else
		{
		  MatrixPsi01[i][j].ReadMatrix(TmpMatrixElementEpsilonSigmaFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightEpsilon, WeightSigma, Weight,
							 RationalMatrixPsi01,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi01[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpMatrixElementEpsilonSigmaFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightEpsilonNumerical, WeightSigmaNumerical, WeightNumerical,
							 MatrixPsi01,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi01[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi01[i][j].WriteMatrix(TmpMatrixElementEpsilonSigmaFileName);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaEpsilonFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi10[i][j].ReadMatrix(TmpMatrixElementSigmaEpsilonFileName);
		}
	      else
		{
		  MatrixPsi10[i][j].ReadMatrix(TmpMatrixElementSigmaEpsilonFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightSigma, WeightEpsilon, Weight,
							 RationalMatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi10[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpMatrixElementSigmaEpsilonFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightSigmaNumerical, WeightEpsilonNumerical, WeightNumerical,
							 MatrixPsi10,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi10[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi10[i][j].WriteMatrix(TmpMatrixElementSigmaEpsilonFileName);
		    }
		}
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaSigmaFileName)))
	    {
	      if (this->UseRationalFlag == true)
		{
		  RationalMatrixPsi11[i][j].ReadMatrix(TmpMatrixElementSigmaSigmaFileName);
		}
	      else
		{
		  MatrixPsi11[i][j].ReadMatrix(TmpMatrixElementSigmaSigmaFileName);
		}
	    }
	  else
	    {
	      if (this->UseRationalFlag == true)
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							 WeightSigma, WeightSigma, Weight,
							 RationalMatrixPsi11,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  RationalMatrixPsi11[i][j] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalMatrixPsi11[i][j].WriteMatrix(TmpMatrixElementSigmaSigmaFileName);
		    }
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							 WeightSigmaNumerical, WeightSigmaNumerical, WeightNumerical,
							 MatrixPsi11,  i - 1, j);
		  Operation2.ApplyOperation(architecture);
		  MatrixPsi11[i][j] = Operation2.GetMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      MatrixPsi11[i][j].WriteMatrix(TmpMatrixElementSigmaSigmaFileName);
		    }
		}
	    }
	  if ((1 + i) <= this->PLevel)
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementPhiSigmaFileName)))
		{
		  if (this->UseRationalFlag == true)
		    {
		      RationalMatrixPsi21[i][j].ReadMatrix(TmpMatrixElementPhiSigmaFileName);
		    }
		  else
		    {
		      MatrixPsi21[i][j].ReadMatrix(TmpMatrixElementPhiSigmaFileName);
		    }
		}
	      else
		{
		  if (this->UseRationalFlag == true)
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightPhi, WeightSigma, Weight,
							     RationalMatrixPsi21,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi21[i][j] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi21[i][j].WriteMatrix(TmpMatrixElementPhiSigmaFileName);
			}
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							     WeightPhiNumerical, WeightSigmaNumerical, WeightNumerical,
							     MatrixPsi21,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      MatrixPsi21[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  MatrixPsi21[i][j].WriteMatrix(TmpMatrixElementPhiSigmaFileName);
			}
		    }
		}
	    }
	  if ((1 + j) <= this->PLevel)
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaPhiFileName)))
		{
		  if (this->UseRationalFlag == true)
		    {
		      RationalMatrixPsi12[i][j].ReadMatrix(TmpMatrixElementSigmaPhiFileName);
		    }
		  else
		    {
		      MatrixPsi12[i][j].ReadMatrix(TmpMatrixElementSigmaPhiFileName);
		    }
		}
	      else
		{
		  if (this->UseRationalFlag == true)
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightSigma, WeightPhi, Weight,
							     RationalMatrixPsi12,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      RationalMatrixPsi12[i][j] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi12[i][j].WriteMatrix(TmpMatrixElementSigmaPhiFileName);
			}
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, j, CentralCharge12Numerical, 
							     WeightSigmaNumerical, WeightPhiNumerical, WeightNumerical,
							     MatrixPsi12,  i - 1, j);
		      Operation2.ApplyOperation(architecture);
		      MatrixPsi12[i][j] = Operation2.GetMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  MatrixPsi12[i][j].WriteMatrix(TmpMatrixElementSigmaPhiFileName);
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
	  for (int i = this->PLevel - 1; i >= 0; --i)
	    {
	      MatrixPsi21[1 + i][j] = MatrixPsi21[i][j];
	    }
	  MatrixPsi21[0][j] = RealMatrix();
	}
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int j = this->PLevel - 1; j >= 0; --j)
	    {
	      MatrixPsi12[i][1 + j] = MatrixPsi12[i][j];
	    }
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
	      if ((1 + j) <= this->PLevel)
		{	
		  for (int k = 0; k < RationalMatrixPsi11[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < RationalMatrixPsi11[i][j].GetNbrColumn(); ++l)
		      {
			RationalMatrixPsi12[i][j][l][k] *= (RationalMultiplicityFactor[i][k] * RationalMultiplicityFactor[j][l]);
		      }
		}
	      if ((1 + i) <= this->PLevel)
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
	      if ((1 + j) <= this->PLevel)
		{	
		  MatrixPsi12[i][1 + j] = RationalMatrixPsi12[i][j];
		}
	      if ((1 + i) <= this->PLevel)
		{	  
		  MatrixPsi21[1 + i][j] = RationalMatrixPsi21[i][j];
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
	      if ((1 + j) <= this->PLevel)
		{	
		  for (int k = 0; k < MatrixPsi01[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < MatrixPsi01[i][j].GetNbrColumn(); ++l)
		      {
			double Tmp;
			MatrixPsi12[i][1 + j].GetMatrixElement(k, l, Tmp);
			Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
			MatrixPsi12[i][1 + j].SetMatrixElement(k, l, Tmp);
		  }
		}
	      if ((1 + i) <= this->PLevel)
		{	  
		  for (int k = 0; k < MatrixPsi10[i][j].GetNbrRow(); ++k)
		    for (int l = 0; l < MatrixPsi10[i][j].GetNbrColumn(); ++l)
		      {
			double Tmp;
			MatrixPsi21[1 + i][j].GetMatrixElement(k, l, Tmp);
			Tmp *= (MultiplicityFactor[i][k] * MultiplicityFactor[j][l]);
			MatrixPsi21[1 + i][j].SetMatrixElement(k, l, Tmp);
		      }
		}
	    }
	  MatrixPsi01[i][j] *= sqrt(2.0 / 3.0);
	  MatrixPsi10[i][j] *= sqrt(2.0 / 3.0);
	  MatrixPsi11[i][j] *= 1.0 / sqrt(3.0);	  
	  if ((1 + j) <= this->PLevel)
	    {	
	      MatrixPsi12[i][1 + j] *= sqrt(7.0 / 2.0) / 3.0;
	    }
	  if ((1 + i) <= this->PLevel)
	    {	  
	      MatrixPsi21[1 + i][j] *= -sqrt(7.0 / 2.0) / 3.0;
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
      this->NbrNValue = QValueDenominator * (2 * this->PLevel) + 5;
      NValueShift = QValueDenominator * this->PLevel;
    }
  else
    {
      QValue = 2 + 3 * (this->LaughlinIndex - 1);
      QValueDenominator = 3;
      this->NbrNValue = QValueDenominator * (2 * this->PLevel) + QValue + 5;
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
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonRight = OrthogonalBasisEpsilonRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  RealSymmetricMatrix& TmpScalarProductEpsilon = ScalarProductEpsilon[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // Epsilon
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
			}
		    }
		}
	      // sigma1
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - QValueDenominator, p, ChargedIndex, NeutralIndex1)];
			}
		    }
		}
	      // sigma2
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][2] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][2]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
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
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // Phi
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][3] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][3]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
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
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonRight = OrthogonalBasisEpsilonRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  RealSymmetricMatrix& TmpScalarProductEpsilon = ScalarProductEpsilon[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // Epsilon
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductEpsilon(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisEpsilonRight(NeutralIndex4, NeutralIndex2);				  
				}
			      Tmp += TmpOrthogonalBasisEpsilonLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	      // sigma1
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductSigma(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigmaRight(NeutralIndex4, NeutralIndex2);				  
				}
				  Tmp += TmpOrthogonalBasisSigmaLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
									   + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
									   + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
				}
			    }
			  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - QValueDenominator, p, ChargedIndex, NeutralIndex1),
							this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
			}
		    }
		}
	      // sigma2
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][2] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][2]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductSigma(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigmaRight(NeutralIndex4, NeutralIndex2);				  
				}
				  Tmp += TmpOrthogonalBasisSigmaLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
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
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      // Phi
	      for (int j = this->NInitialValuePerPLevelCFTSector[i][3] + QValueDenominator; j <= this->NLastValuePerPLevelCFTSector[i][3]; ++j)
		{
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductPhi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhiRight(NeutralIndex4, NeutralIndex2);
				}
			      Tmp += TmpOrthogonalBasisPhiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    {
			      Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
								       + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								       + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			    }
			  else
			    {
			      if (this->TorusFlag)
				{
				  Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
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
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	  for (int p = 0; p <= (i - 1); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
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
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 1); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
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
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisEpsilon1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
											   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisEpsilon2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
											   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										       + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
					  BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightSigmaNumerical + ((double) (i + j))
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
	  for (int p = 0; p <= (i - 1); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
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
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 1); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue - QValueDenominator;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  if (this->CylinderFlag)
					    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
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
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
					{
					  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1)];
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	  for (int p = 0; p <= (i - 1); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
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
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 1); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
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
			      N2 = QValueDenominator * (j - i) + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisEpsilon1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 2, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisEpsilon2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
					    }
					  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				    }
				}
			      
			      N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][2]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][2]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
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
	  for (int p = 0; p <= (i - 1); ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
		      
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
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][3]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][3]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][2]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][2])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
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
	      RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= (j - 1); ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		      RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		      RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
		      
		      for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
			{	      
			  TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
			  for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			    {	      
			      TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			      int N2;
			      int N1;
			      N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			      N1 = N2 + QValue;
			      if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				  && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][3]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][3])))
				{ 
				  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
				    {
				      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
					{
					  double Tmp = 0.0;
					  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
					    {
					      double Tmp1 = 0.0;			      
					      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
						{
						  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
						}
					      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
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
  delete[] ScalarProductSigma;
  delete[] ScalarProductPhi;
  delete[] ScalarProductEpsilon;
  delete[] RationalScalarProductSigma;
  delete[] RationalScalarProductPhi;
  delete[] RationalScalarProductEpsilon;
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
  delete[] OrthogonalBasisSigmaLeft;
  delete[] OrthogonalBasisEpsilonLeft;
  delete[] OrthogonalBasisPhiLeft;
  delete[] OrthogonalBasisSigmaRight;
  delete[] OrthogonalBasisEpsilonRight;
  delete[] OrthogonalBasisPhiRight;
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSReadRezayi3QuasiholeSectorMatrix::GetQValueCFTSectorShift(int cftSector)
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

void FQHEMPSReadRezayi3QuasiholeSectorMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
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
		      TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
		      TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
      if (cftSector == 1) 
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
		  TmpP += (QPrime  - 2) / QValueDenominator - NValueShift;
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
			  TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
		  TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
      if (cftSector == 2) 
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
		  TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
			  TmpP += (QPrime) / QValueDenominator - NValueShift;
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
		  TmpP -= (QPrime) / QValueDenominator - NValueShift;
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
			  TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
		      TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
		      TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
      if (cftSector == 1) 
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
		  TmpP += (QPrime  - 2) / QValueDenominator - NValueShift;
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
			  TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
		  TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
      if (cftSector == 2) 
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
		  TmpP += (QPrime - 1) / QValueDenominator - NValueShift;
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
			  TmpP += (QPrime) / QValueDenominator - NValueShift;
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
		  TmpP -= (QPrime) / QValueDenominator - NValueShift;
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
			  TmpP -= (QPrime - 1) / QValueDenominator - NValueShift;
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
