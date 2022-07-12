////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                         in their quasihole sector                          //
//                                                                            //
//                        last modification : 11/02/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeSectorMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
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

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix()
{
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
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

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool useRational, 
										 bool trimChargeIndices, bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->BosonicVersion = bosonicVersion;
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
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
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  if (this->SelfDualFlag == true)
    {
      this->TransferMatrixDegeneracy /= 2;
      this->NbrCFTSectors = 1;
    }
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_qh_r_%d", this->RIndex);
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// cftDirectory = path to the directory where all the pure CFT matrices are stored
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

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion, 
										 bool useRational, bool trimChargeIndices,
										 bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->BosonicVersion = bosonicVersion;
  this->RIndex = rIndex;
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
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  if (this->SelfDualFlag == true)
    {
      this->TransferMatrixDegeneracy /= 2;
      this->NbrCFTSectors = 1;
    }
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_qh_r_%d", this->RIndex);
  this->CreateBMatrices(cftDirectory, architecture);
}

// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int pLevel, int nbrBMatrices, char* fileName, bool bosonicVersion, bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->BosonicVersion = bosonicVersion;
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
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightSigma", this->WeightSigma);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPhi", this->WeightPhi);
      ErrorFlag = StateDefinition.GetAsBoolean("SelfDual", this->SelfDualFlag);
      ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
      if (StateDefinition["PsiSquareMatrixElement"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleLongRational("PsiSquareMatrixElement", this->SquareMatrixElementNormalization);
	}
      else
	{
	  this->SquareMatrixElementNormalization = LongRational(1, 2);
	}
      this->MatrixElementNormalization = sqrt(fabs(this->SquareMatrixElementNormalization.GetNumericalValue()));
      if (StateDefinition["EMatrixDegeneracy"] != 0)
	{
	  ErrorFlag = StateDefinition.GetAsSingleInteger("EMatrixDegeneracy", this->TransferMatrixDegeneracy);	  
	}
      else
	{
	  switch (this->RIndex)
	    {
	    case 2:
	      this->TransferMatrixDegeneracy = 2;
	      break;
	    case 3:
	      this->TransferMatrixDegeneracy = 5;
	      break;
	    case 6:
	      this->TransferMatrixDegeneracy = 5;
	      break;
	    }
	}
      if (StateDefinition["Name"] != 0)
	{
	  this->BMatrixOutputName = new char[strlen(StateDefinition["Name"]) + 1]; 
	  strcpy(this->BMatrixOutputName, StateDefinition["Name"]);
	}
      else
	{
	  this->BMatrixOutputName = new char[256]; 
	  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
	}
      this->NbrCFTSectors = 2;
      if (this->SelfDualFlag == true)
	{
	  this->NbrCFTSectors = 1;
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
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction
// kappa = cylinder aspect ratio

FQHEMPSClustered2RQuasiholeSectorMatrix::FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa, 
										 bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion)
{
  this->RIndex = rIndex;
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
  this->LoadMatrices(fileName);
  this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SelfDualFlag = ((this->RIndex & 1) == 0);
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
  this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));
  this->MatrixElementNormalization = 1.0 / M_SQRT2;
  this->SquareMatrixElementNormalization = LongRational(1, 2);
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_qh_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2RQuasiholeSectorMatrix::~FQHEMPSClustered2RQuasiholeSectorMatrix()
{
}
  
// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2RQuasiholeSectorMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge12 (this->CentralCharge);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  double WeightSigmaNumerical = this->WeightSigma.GetNumericalValue();
  double WeightPhiNumerical = this->WeightPhi.GetNumericalValue();
  double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
  double WeightPrimaryFieldMatrixElementNumerical = this->WeightPrimaryFieldMatrixElement.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductSigma = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPhi = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductSigma = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPhi = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiRight = new RealMatrix[this->PLevel + 1];
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
	    {
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

  cout << "weight: " <<   this->WeightSigma << " " << this->WeightPhi << endl;
  char* TmpScalarProductSigmaFileName = 0; 
  char* TmpScalarProductPhiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPhiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      if (cftDirectory != 0)
	{
	  if (this->UseRationalFlag == true)
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      if (this->SelfDualFlag == false)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	  else
	    {
	      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	      if (this->SelfDualFlag == false)
		sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_num_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
	    }
	}
      this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductSigmaFileName, architecture, RationalScalarProductSigma, ScalarProductSigma, i, U1BosonBasis,
					   CentralCharge12, CentralCharge12Numerical, this->WeightSigma, WeightSigmaNumerical, "sigma",
					   OrthogonalBasisSigmaLeft, OrthogonalBasisSigmaRight, RationalMultiplicityFactor, MultiplicityFactor);
      if (this->SelfDualFlag == false)
	{
	  this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductPhiFileName, architecture, RationalScalarProductPhi, ScalarProductPhi, i, U1BosonBasis,
					       CentralCharge12, CentralCharge12Numerical, this->WeightPhi, WeightPhiNumerical, "phi",
					       OrthogonalBasisPhiLeft, OrthogonalBasisPhiRight, RationalMultiplicityFactor, MultiplicityFactor);
	}
      
      cout << "---------------------------------" << endl;
    }

  this->RescaleFullScalarProductMatrix(RationalScalarProductSigma, ScalarProductSigma, RationalMultiplicityFactor, MultiplicityFactor);
  if (this->SelfDualFlag == false)
    {
      this->RescaleFullScalarProductMatrix(RationalScalarProductPhi, ScalarProductPhi, RationalMultiplicityFactor, MultiplicityFactor);
    }

  this->U1BasisDimension = new int [this->PLevel + 1];	
  this->NeutralSectorDimension = new int* [this->NbrCFTSectors];
  for (int i = 0; i < this->NbrCFTSectors; ++i)
    this->NeutralSectorDimension[i] = new int [this->PLevel + 1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->NeutralSectorDimension[0][i] = OrthogonalBasisSigmaLeft[i].GetNbrColumn();
      if (this->SelfDualFlag == false)
	  this->NeutralSectorDimension[1][i] = OrthogonalBasisPhiLeft[i].GetNbrColumn();
      this->U1BasisDimension[i] = U1BosonBasis[i]->GetHilbertSpaceDimension();
    }
  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  double ExtraCylinderFactor = 1.0;
  if (this->BosonicVersion == false)
    {
      if ((this->RIndex & 1) == 0)
	{
	  QValue = 1 + (this->RIndex / 2);
	  this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 1;
	  NValueShift = 2 * this->PLevel - 1;
	  --this->NbrNValue;
	  ++NValueShift;
	  QValueDenominator = 1;
	}
      else
	{
	  QValue = 2 + this->RIndex;
	  this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 1;
	  NValueShift = 4 * this->PLevel - 2;
	  this->NbrNValue -= 2;
	  NValueShift += 2;
	  QValueDenominator = 2;
	  ExtraCylinderFactor = 4.0;
	}
    }
  else
    {
      if ((this->RIndex & 1) == 0)
	{
	  QValue = this->LaughlinIndex - 1 + (this->RIndex / 2);
	  this->NbrNValue = ((2 * this->PLevel) + QValue) + this->RIndex / 2 + 2;
	  NValueShift = 2 * this->PLevel - 1;
	  --this->NbrNValue;
	  ++NValueShift;
	  QValueDenominator = 1;
	}
      else
	{
	  QValue = ((2 * this->LaughlinIndex) - 2) + this->RIndex;
	  this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 2;
	  NValueShift = 4 * this->PLevel - 2;
	  this->NbrNValue -= 2;
	  NValueShift += 2;
	  QValueDenominator = 2;
	  ExtraCylinderFactor = 4.0;
	}
    }
  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (this->WeightPsi);
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      if (this->UseRationalFlag == true)
		{
		  if (this->SelfDualFlag == false)
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		      sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		  else
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		}
	      else
		{
		  if (this->SelfDualFlag == false)
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		      sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_num_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		  else
		    {
		      sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, i, j);
		    }
		}
	    }
	  if (this->SelfDualFlag == false)
	    {
	      this->ComputeFullMatrixElements(cftDirectory, TmpScalarProductSigmaFileName, architecture, 
					      RationalMatrixPsi01, MatrixPsi01, i, j, U1BosonBasis, 
					      CentralCharge12, CentralCharge12Numerical, 
					      this->WeightSigma, WeightSigmaNumerical, 
					      this->WeightPhi, WeightPhiNumerical,
					      this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	      this->ComputeFullMatrixElements(cftDirectory, TmpScalarProductPhiFileName, architecture, 
					      RationalMatrixPsi10, MatrixPsi10, i, j, U1BosonBasis, 
					      CentralCharge12, CentralCharge12Numerical, 
					      this->WeightPhi, WeightPhiNumerical, 
					      this->WeightSigma, WeightSigmaNumerical,
					      this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	    }
	  else
	    {
	      this->ComputeFullMatrixElements(cftDirectory, TmpScalarProductSigmaFileName, architecture, 
					      RationalMatrixPsi01, MatrixPsi01, i, j, U1BosonBasis, 
					      CentralCharge12, CentralCharge12Numerical, 
					      this->WeightSigma, WeightSigmaNumerical, 
					      this->WeightSigma, WeightSigmaNumerical,
					      this->WeightPrimaryFieldMatrixElement, WeightPrimaryFieldMatrixElementNumerical);
	    }
	}
    }


  this->RescaleFullMatrixElements(RationalMatrixPsi01, MatrixPsi01, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);
  if (this->SelfDualFlag == false)
    {
      this->RescaleFullMatrixElements(RationalMatrixPsi10, MatrixPsi10, RationalMultiplicityFactor, MultiplicityFactor, this->MatrixElementNormalization);
    }

  cout << "building B matrices" << endl;

  SparseRealMatrix* BMatrices = 0;
  SparseComplexMatrix* TmpComplexBMatrices = 0;
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
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
			{
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
				}
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
			{
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1)];
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
      BMatrices = new SparseRealMatrix[this->NbrBMatrices];
      BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
    }
  else
    {
      TmpComplexBMatrices = new SparseComplexMatrix[this->NbrBMatrices];
      TmpComplexBMatrices[0] = SparseComplexMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((this->RIndex & 1) == 0)
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 1; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				    }
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightSigmaNumerical +  ((double) i)
									      + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 1; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
				  if (this->TwistedTorusFlag == false)
				    {
				      if (this->CylinderFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				      else
					{
					  if (this->TorusFlag)
					    {
					      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
										       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					    }
					}
				      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
								    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
				    }
				  else
				    {
				      Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightPhiNumerical +  ((double) i)
										  + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										  + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				      TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
									      this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				    }
				}
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
									       + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				    }
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
										   + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightSigmaNumerical +  ((double) i)
									      + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		  if (this->SelfDualFlag == false)
		    {
		      for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
				  if (this->TwistedTorusFlag == false)
				    {
				      if (this->CylinderFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
										   + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
					}
				      else
					{
					  if (this->TorusFlag)
					    {
					      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPhiNumerical +  ((double) i)
										       + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
					    }
					}
				      BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1),
								    this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
				    }
				  else
				    {
				      Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightPhiNumerical +  ((double) i)
										  + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										  + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				      TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1),
									      this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				    }
				}
			    }
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
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
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
				  N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue - 1;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
				}	
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					    } 
					}
				    }
				}
			      else
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					    } 
					}
				    }
				}
			      
			      if ((this->RIndex & 1) == 0)
				{
				  N2 = (2 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 1;
				}
			      else
				{
				  N2 = (4 * (j - i) + 4 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
				}
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	}

      BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
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
				  N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue - 1;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
				}	
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
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
						      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      if (this->CylinderFlag)
						Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
											       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
											       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
					      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									    this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    } 
					}
				    }
				}
			      else
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
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
						      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      if (this->CylinderFlag)
						Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
											       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))
											       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift)) * QValueDenominator / (2.0 * ExtraCylinderFactor * QValue))));
					      BMatrices[1].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									    this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    } 
					}
				    }
				}
			      
			      if ((this->RIndex & 1) == 0)
				{
				  N2 = (2 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 1;
				}
			      else
				{
				  N2 = (4 * (j - i) + 4 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
				}
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
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
						      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      if (this->CylinderFlag)
						Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightPhiNumerical + ((double) (i + j))
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
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
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
				  N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue;
				}	
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					    } 
					}
				    }
				}
			      else
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
					    {
					      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1)];
					    } 
					}
				    }
				}
			      
			      if ((this->RIndex & 1) == 0)
				{
				  N2 = (2 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 4 + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
				    { 
				      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
					{
					  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
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
	}

      SparseRealMatrix V0Matrix (MatrixSize, MatrixSize, TmpNbrElementPerRow);
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  for (int p = 0; p <= i; ++p)
	    {
	      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	      RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	      RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	      for (int j = 0; j <= this->PLevel; ++j)
		{
		  for (int q = 0; q <= j; ++q)
		    {
		      BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		      BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		      RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		      RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
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
				  N2 = (2 * (j - i) + this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + NValueShift) / 2;
				  N1 = N2 + QValue;
				}	
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][1])))
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
						      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									    this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    } 
					}
				    }
				}
			      else
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][0]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
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
						      Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									    this->Get2RMatrixIndexV2(j, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					    } 
					}
				    }
				}
			      
			      if ((this->RIndex & 1) == 0)
				{
				  N2 = (2 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 4 + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      if (this->SelfDualFlag == false)
				{		  
				  if (((N1 >= this->NInitialValuePerPLevelCFTSector[i][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[i][1]))
				      && ((N2 >= this->NInitialValuePerPLevelCFTSector[j][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[j][0])))
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
						      Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
						    }
						  Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
						}
					      Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
					      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, N1, p, ChargedIndex1, NeutralIndex1),
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
	}

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
  delete[] ScalarProductSigma;
  delete[] ScalarProductPhi;
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
  delete[] OrthogonalBasisSigmaLeft;
  delete[] OrthogonalBasisPhiLeft;
  delete[] OrthogonalBasisSigmaRight;
  delete[] OrthogonalBasisPhiRight;
  delete[] RationalMultiplicityFactor;
  delete[] MultiplicityFactor;
}

// get the (k,r) exclude principle satisfied by the root configuration
// 
// pauliK = maximum number of particles in the pauliR consecutive orbitals
// pauliR = number of consecutive orbitals

void FQHEMPSClustered2RQuasiholeSectorMatrix::GetKRExclusionPrinciple(int& pauliK, int& pauliR)
{
  pauliK = 2;
  pauliR = 2 * (this->LaughlinIndex - 1) + this->RIndex;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSClustered2RQuasiholeSectorMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  if (((this->RIndex & 1) == 0) && (this->SelfDualFlag == true))
    {
      numerator = 1;
      denominator = (this->LaughlinIndex - 1) + (this->RIndex / 2);
    }
  else
    {
      numerator = 2;
      denominator = 2 * (this->LaughlinIndex - 1) + this->RIndex;
    }
}


// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSClustered2RQuasiholeSectorMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
      if (this->BosonicVersion == false)
	{
	  if ((this->RIndex & 1) == 0)
	    rowIndex = this->PLevel + (this->RIndex / 2) - MinQ;
	  else
	    rowIndex = 2 * this->PLevel + this->RIndex - 1 - MinQ;      
	  columnIndex = rowIndex;
	}
      else
	{
	  if ((this->RIndex & 1) == 0)
	    rowIndex = this->PLevel - MinQ;
	  else
	    rowIndex = 2 * this->PLevel - MinQ;      
	  columnIndex = rowIndex;
	}
    }
  else
    {
      if (this->BosonicVersion == false)
	{
	  if ((this->RIndex & 1) == 0)
	    {
	      rowIndex = this->PLevel + this->RIndex - MinQ;
	      columnIndex = this->PLevel + this->RIndex - 1 - MinQ ;
	    }
	  else
	    {
	      rowIndex = 2 * (this->PLevel + this->RIndex) - MinQ;
	      columnIndex = 2 * (this->PLevel + this->RIndex - 2) - MinQ;
	    }
	}
      else
	{
	  if ((this->RIndex & 1) == 0)
	    {
	      rowIndex = this->PLevel + this->LaughlinIndex - MinQ;
	      columnIndex = this->PLevel + this->LaughlinIndex - MinQ ;
	      cout << rowIndex << " " << columnIndex << " " <<  this->LaughlinIndex << endl;
	    }
	  else
	    {
	      rowIndex = 2 * (this->PLevel + this->LaughlinIndex) - MinQ;
	      columnIndex = 2 * this->PLevel - MinQ;
	    }
	}
    }
}


// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RQuasiholeSectorMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
{
  if (this->UniformChargeIndexRange == true)
    {
      minQ = 0;
      maxQ = this->NbrNValue - 1;
      return;
    }

  if (this->BosonicVersion == false)
    {
      if (this->SelfDualFlag == false)
	{
	  if ((this->RIndex & 1) == 0)
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
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
	      minQ = TmpMinQ;
	      maxQ = TmpMaxQ;
	    }
	  else
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
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
			      TmpP += (QPrime - 2) / 2 - NValueShift;
			    }
			}
		      QPrime = Q;
		      TmpP = 0;
		      int TmpMaxP2 = 0;
		      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= ((QPrime - 2) / 2) - NValueShift;
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
		  minQ = TmpMinQ;
		  maxQ = TmpMaxQ;
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
			  TmpP += (QPrime - 2 ) / 2 - NValueShift;
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
			      TmpP -= ((QPrime - 2) / 2) - NValueShift;
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
		  minQ = TmpMinQ;
		  maxQ = TmpMaxQ;
		}
	    }
	}
      else
	{
	  minQ = 0;
	  maxQ = this->NbrNValue - 1;
	  
	  if ((this->RIndex & 1) == 0)
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
	      int NValueShift = this->PLevel;
	      int QValue = 1 + (this->RIndex / 2);
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
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= QPrime - NValueShift - (this->RIndex / 2);
		      QPrime += (QValue - 1);
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	
		    }    
		}
	      minQ = TmpMinQ;
	      maxQ = TmpMaxQ;
	    }
	  else
	    {
	      minQ = 0;
	      maxQ = this->NbrNValue - 1;
	    }
	}
    }
  else
    {
      if (this->SelfDualFlag == false)
	{
	  if ((this->RIndex & 1) == 0)
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
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
	      minQ = TmpMinQ;
	      maxQ = TmpMaxQ;
	    }
	  else
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
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
			      TmpP += (QPrime - 2) / 2 - NValueShift;
			    }
			}
		      QPrime = Q;
		      TmpP = 0;
		      int TmpMaxP2 = 0;
		      while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
			{
			  if (TmpP > TmpMaxP2)
			    TmpMaxP2 = TmpP;	    
			  TmpP -= ((QPrime - 2) / 2) - NValueShift;
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
		  minQ = TmpMinQ;
		  maxQ = TmpMaxQ;
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
			  TmpP += (QPrime - 2 ) / 2 - NValueShift;
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
			      TmpP -= ((QPrime - 2) / 2) - NValueShift;
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
		  minQ = TmpMinQ;
		  maxQ = TmpMaxQ;
		}
	    }
	}
      else
	{
	  minQ = 0;
	  maxQ = this->NbrNValue - 1;
	  
	  if ((this->RIndex & 1) == 0)
	    {
	      int TmpMinQ = this->NbrNValue - 1;
	      int TmpMaxQ = 0;    
	      int NValueShift = this->PLevel;
	      int QValue = 1 + (this->RIndex / 2);
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
		    }
		  QPrime = Q;
		  TmpP = 0;
		  int TmpMaxP2 = 0;
		  while ((TmpP >= 0) && (QPrime < this->NbrNValue) && (QPrime >= 0))
		    {
		      if (TmpP > TmpMaxP2)
			TmpMaxP2 = TmpP;	    
		      TmpP -= QPrime - NValueShift - (this->RIndex / 2);
		      QPrime += (QValue);
		    }
		  if (((this->PLevel - TmpMaxP) >= pLevel) && ((this->PLevel - TmpMaxP2) >= pLevel))
		    {
		      if (Q < TmpMinQ)
			TmpMinQ = Q;
		      if (Q > TmpMaxQ)
			TmpMaxQ = Q;	
		    }    
		}
	      minQ = TmpMinQ;
	      maxQ = TmpMaxQ;
	    }
	  else
	    {
	      minQ = 0;
	      maxQ = this->NbrNValue - 1;
	    }
	}
    }
  cout << "range at " << pLevel << " : " << minQ << " " << maxQ << " (" << this->NbrNValue << ")" << endl;   
}


// get the matrix that into account the Jordan Wigner string on the torus geometry
//
// nbrFermions = number of fermions in the system
// return value = corresponding matrix

SparseRealMatrix FQHEMPSClustered2RQuasiholeSectorMatrix::GetTorusStringMatrix(int nbrFermions)
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
	      if ((CurrentQValue & 2) == 0)
		Coefficient = -1.0;
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

int* FQHEMPSClustered2RQuasiholeSectorMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
{
  nbrIndices = 0;
  if ((this->RIndex & 1) == 0)
    {
      int ModSector = (this->RIndex / 2) + (this->LaughlinIndex - 1);
      for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue;
	      int MaxQValue;
	      this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
		{
		  if ((CurrentQValue % ModSector) == topologicalSector)
		    {
		      nbrIndices += this->GetBondIndexRange(CurrentPLevel, CurrentQValue, CurrentCFTSector);
		    }
		}
	    }
	}
    }
  else
    {
      for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue;
	      int MaxQValue;
	      this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
		{
		  if ((CurrentQValue % (this->RIndex + 2 * (this->LaughlinIndex - 1))) == topologicalSector)
		    {
		      nbrIndices += this->GetBondIndexRange(CurrentPLevel, CurrentQValue, CurrentCFTSector);
		    }
		}
	    }
	}
    }
  int* TmpIndices =  new int [nbrIndices];
  nbrIndices = 0;
  if ((this->RIndex & 1) == 0)
    {
      int ModSector = (this->RIndex / 2) + (this->LaughlinIndex - 1);
      for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue;
	      int MaxQValue;
	      this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
		{
		  if ((CurrentQValue % ModSector) == topologicalSector)
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
    }
  else
    {
      for (int CurrentPLevel = 0; CurrentPLevel <= this->PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < this->NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue;
	      int MaxQValue;
	      this->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      for (int CurrentQValue = MinQValue; CurrentQValue <= MaxQValue; ++CurrentQValue)
		{
		  if ((CurrentQValue % (this->RIndex + 2 * (this->LaughlinIndex - 1))) == topologicalSector)
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
    }
  return TmpIndices;
}
  
// get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
// 
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true if we are dealing with fermions
// return value = minimum ky momentum 

int FQHEMPSClustered2RQuasiholeSectorMatrix::GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics)
{
  if (statistics == false)
    {
      if ((nbrParticles & 1) == 0)
	return (nbrParticles / 2);
      else
	return 0;
    }
  return 0;
}
