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
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
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

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix()
{
  this->UseRationalFlag = true;
  this->BosonicVersion = false;
  this->TwistedTorusFlag = false;
  this->NeutralSectorMaxPLevel = -1;
  this->CFTDirectory = 0;
  this->Architecture = 0;
}

// constructor
//
// pLevel = |P| level truncation
// centralCharge = value of the central charge
// outputName = name of the theory
// useRational = use arbitrary precision numbers for all the CFT calculations

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int pLevel, LongRational centralCharge, const char* outputName, bool useRational)
{
    this->PLevel = pLevel;
    this->NeutralSectorMaxPLevel = this->PLevel;
    this->CentralCharge = centralCharge;
    this->BMatrixOutputName = new char[512];
    sprintf(this->BMatrixOutputName, "%s", outputName);
    this->UseRationalFlag = useRational;

    this->TotalStartingIndexPerPLevel = 0;
    this->NbrIndicesPerPLevel = 0;
    this->NbrNValuesPerPLevel = 0;
    this->NInitialValuePerPLevel = 0;
    this->NLastValuePerPLevel = 0;
    this->RealBMatrices = 0;
    this->ComplexBMatrices = 0;
    this->QuasiholeBMatrices = 0;
    this->TwistedTorusFlag = false;
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

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool useRational,
						   bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RIndex = rIndex;
  this->BosonicVersion = bosonicVersion;
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  this->LaughlinIndex = laughlinIndex;
  this->UseRationalFlag = useRational;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevel;
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
  this->WeightPrimaryFieldMatrixElement = LongRational(this->RIndex, 4l);
  this->WeightIdentity = LongRational(0l, 1l);
  this->WeightPsi = LongRational(this->RIndex, 4l);
  this->CentralCharge = LongRational((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
  this->SquareMatrixElementNormalization = LongRational(1, 1);
  this->MatrixElementNormalization = 1.0;
  this->TransferMatrixDegeneracy = this->RIndex + 2;
  this->NbrCFTSectors = 2;
  this->BMatrixOutputName = new char[256]; 
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
  this->CreateBMatrices(0, architecture);
  this->CFTDirectory = 0;
  this->Architecture = architecture;
}

// constructor 
//
// rindex = r index (i.e. clustered (k=2,r) states) 
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
  
FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, const char* cftDirectory, bool bosonicVersion, bool useRational, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RIndex = rIndex;
  this->BosonicVersion = bosonicVersion;
  this->LaughlinIndex = laughlinIndex;
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  this->PLevel = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevel;
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
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
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


// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// bosonicVersion = use a version of the code that is compatible with bosonic wave functions
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction
// architecture = architecture to use for precalculation

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int pLevel, int nbrBMatrices, const char* fileName, bool bosonicVersion, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion,
						   AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->BosonicVersion = bosonicVersion;
  this->PhysicalIndices = new unsigned long[this->NbrBMatrices];
  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->PhysicalIndices[i] = (unsigned long) i;
    }
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
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
  this->PLevel = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevel;
  
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
	  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
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
  if (StateDefinition["CFTMatrixDirectory"] != 0)
    {
      this->CFTDirectory = new char [strlen(StateDefinition["CFTMatrixDirectory"]) + 1];
      strcpy (this->CFTDirectory, StateDefinition["CFTMatrixDirectory"]);
    }
  else
    {
      this->CFTDirectory = 0;
    }
  this->Architecture = architecture;
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
// torusFlag = true the torus geometry should be used instead of a genus-0 surface
// nbrFluxQuanta = number of flux quanta piercing the torus
// aspectRatio = aspect ratio of the torus(norm of tau)
// angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
// fluxInsertion = flux insertion along the tau direction

FQHEMPSClustered2RMatrix::FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, const char* fileName, 
						   bool trimChargeIndices, bool cylinderFlag, double kappa, 
						   bool torusFlag, int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion)
{
  this->RIndex = rIndex;
  this->LaughlinIndex = laughlinIndex;
  this->UniformChargeIndexRange = !trimChargeIndices;
  this->PLevel = pLevel;
  this->NeutralSectorMaxPLevel = this->PLevel;
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
  sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2RMatrix::~FQHEMPSClustered2RMatrix()
{
  delete[] this->BMatrixOutputName;
  if (this->CFTDirectory != 0)
    {
      delete[] this->CFTDirectory;
    }
}
  
// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void FQHEMPSClustered2RMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 2;
  denominator = 2  * (this->LaughlinIndex - 1) + this->RIndex;
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2RMatrix::CreateBMatrices (const char* cftDirectory, AbstractArchitecture* architecture)
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
  if (this->BosonicVersion == false)
    {
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
	  this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 1;
	  NValueShift = 4 * this->PLevel - 2;
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
	  QValueDenominator = 1;
	}
      else
	{
	  QValue = ((2 * this->LaughlinIndex) - 2) + this->RIndex;
	  this->NbrNValue = ((4 * this->PLevel) + QValue) + this->RIndex + 2;
	  NValueShift = 4 * this->PLevel - 2;
	  QValueDenominator = 2;
	  ExtraCylinderFactor = 4.0;
	}
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
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1)];
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
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				    }
				  
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightIdentityNumerical +  ((double) i)
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 1, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
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
			      if (this->TwistedTorusFlag == false)
				{
				  if (this->CylinderFlag)
				    {
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									       + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
										   + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightPsiNumerical +  ((double) i)
									      + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (4.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (4.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 1, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		}
	      else
		{
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][0] + 2; j <= this->NLastValuePerPLevelCFTSector[i][0]; ++j)
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
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
									       + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightIdentityNumerical +  ((double) i)
										   + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightIdentityNumerical +  ((double) i)
									      + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									      + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				  TmpComplexBMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, j - 2, p, ChargedIndex, NeutralIndex1),
									  this->Get2RMatrixIndexV2(i, 0, j, p, ChargedIndex, NeutralIndex2), Tmp2);
				}
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[i][1] + 2; j <= this->NLastValuePerPLevelCFTSector[i][1]; ++j)
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
				      Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
									       + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
									       + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
				    }			      
				  else
				    {
				      if (this->TorusFlag)
					{
					  Tmp *= exp(-this->Kappa * this->Kappa * (WeightPsiNumerical +  ((double) i)
										   + ((j - 2.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) * QValueDenominator / (16.0 * QValue))
										   + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) * QValueDenominator / (16.0 * QValue))));
					}
				    }
				  BMatrices[0].SetMatrixElement(this->Get2RMatrixIndexV2(i, 1, j - 2, p, ChargedIndex, NeutralIndex1),
								this->Get2RMatrixIndexV2(i, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
				}
			      else
				{
				  Complex Tmp2 = Tmp * exp(this->TauFactor * (WeightPsiNumerical +  ((double) i)
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
				  N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
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
				  N2 = (4 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
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
				  N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
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
				  N2 = (4 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue - 2;
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
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
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
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
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
      
      SparseRealMatrix V0Matrix (MatrixSize, MatrixSize, TmpNbrElementPerRow);
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
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 * this->RIndex + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
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
					  V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(i, 0, N1, p, ChargedIndex1, NeutralIndex1),
									this->Get2RMatrixIndexV2(j, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
					}
				      
				    }
				}
			      if ((this->RIndex & 1) == 0)
				{
				  N2 = (2 * (j - i) + 1 + NValueShift) / 2;
				  N1 = N2 + QValue;
				}
			      else
				{
				  N2 = (4 * (j - i) + 2 + NValueShift) / 2;
				  N1 = N2 + QValue;
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

// compute the scalar product matrices of the Virasoro descendant
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// return value = scalar product

LongRational FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									       LongRational& centralCharge12, LongRational& weight)
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
	  LongRational Tmp1 (centralCharge12);
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  LongRational Tmp2 (weight);
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
	}
    }
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
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight));
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
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

LongRational FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									       LongRational& centralCharge12, LongRational& weight,
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
	  LongRational Tmp1 (centralCharge12);
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  LongRational Tmp2 (weight);
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
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
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if (partition[TmpPosition - 1] < partition[TmpPosition])
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
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
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
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
	      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
	      return Tmp;
	    }
	}
    }
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
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight,
							   precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
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
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight,
						      precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}

// compute the scalar product matrices of the Virasoro descendant, using information from previous levels
// 
// partition = partition that desribes the product of Virasoro generators involved in the scalar product
// partitionLength = partition length
// position = position in partition starting from which all the indices are negative
// centralCharge12 = reference on the value of the central charge divided by 12
// weight = weight of the primary field that is considered
// precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
// precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
// basis = basis that related the partitions to their index
// temporaryOccupationNumber = local temporary to store the occupation numbers 
// return value = scalar product

double FQHEMPSClustered2RMatrix::ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
									 double& centralCharge12, double& weight,
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
    return 0.0;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0.0;
	}
      else
	{
	  double Tmp1 = centralCharge12;
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  double Tmp2 = weight;
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
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
      while ((TmpPosition < position) && (FlagSorted == true))
	{
	  TmpPLevel1 += partition[TmpPosition];
	  if (partition[TmpPosition - 1] < partition[TmpPosition])
	    FlagSorted = false;
	  ++TmpPosition;
	}
      if ((TmpPLevel1 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
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
	  if ((TmpPLevel2 <= precomputedScalarProductMaxPLevel) && (FlagSorted == true))
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
	      double Tmp;
	      precomputedScalarProduct[TmpPLevel1].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
	      return Tmp;
	    }
	}
    }
  double Tmp = 0l;
  if ((partition[position - 1] + partition[position]) == 0)
    {
      long TmpLength = 0l;
      long Store = partition[position - 1];
      for (int i = position + 1; i < partitionLength; ++i)
	TmpLength += partition[i];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 2] = partition[i];
      Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
	       + (2l * Store) * (weight - TmpLength)) * 
	      this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight,
							   precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
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
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, 
								 position - 1, centralCharge12, weight,
								 precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += this->ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight,
						      precomputedScalarProduct, precomputedScalarProductMaxPLevel, basis, temporaryOccupationNumber);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
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
  
LongRational FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
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
  
LongRational FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
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

// compute the matrix elements of any primary field in the Virasoro descendant basis, using double numbers instead of long rational
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
  
double FQHEMPSClustered2RMatrix::ComputeDescendantMatrixElement (long* partition, int partitionLength, 
								 int descendantPosition, int position,
								 double& centralCharge12, double& weight1, 
								 double& weight2, double& weight,
								 RealMatrix** precomputedDescendantMatrixElement, 
								 int precomputedDescendantMatrixElementMaxLeftPLevel, 
								 int precomputedDescendantMatrixElementMaxRightPLevel, 
								 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber)
{
  if (partitionLength == 0)
    {
      return 1.0;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0.0;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if (descendantPosition == partitionLength) 
    {
      double Tmp = 1.0;
      double TmpSum = weight1;
      TmpSum -= weight2;
      double Tmp2;
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
      double Tmp = 1.0;
      double TmpSum = weight1;
      TmpSum -= weight2;
      double Tmp2;
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
      Tmp *= 1.0 - (2.0 * (partitionLength & 1l));
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
		  double Tmp;
		  precomputedDescendantMatrixElement[TmpPLevel1][TmpPLevel2].GetMatrixElement(TmpIndex1, TmpIndex2, Tmp);
		  return Tmp;
		}
	    }
	}
    }
  if (descendantPosition < position)
    {
      double Tmp = 0.0;
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

  double Tmp1 = this->ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
							   centralCharge12, weight1, weight2, weight,
							   precomputedDescendantMatrixElement, precomputedDescendantMatrixElementMaxLeftPLevel,
							   precomputedDescendantMatrixElementMaxRightPLevel, basis, temporaryOccupationNumber);


  double Tmp2 = weight;
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

// compute the various arrays required to convert from quantum numbers and local indices to a global linearized index
//
// return value = dimension of the B matrix

int FQHEMPSClustered2RMatrix::ComputeLinearizedIndexArrays()
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
	      for (int k = 0; k <= i; ++k)
		{
		  this->StartingIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = TotalIndex;
		  int Tmp = this->U1BasisDimension[k] * this->NeutralSectorDimension[l][i - k];
		  this->NbrIndexPerPLevelCFTSectorQValueU1Sector[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]][k] = Tmp;
		  TotalIndex += Tmp;	      
		}
	      this->NbrIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]] = TotalIndex - this->StartingIndexPerPLevelCFTSectorQValue[i][l][j - this->NInitialValuePerPLevelCFTSector[i][l]];
	    }
	}
    }
  return TotalIndex;
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int FQHEMPSClustered2RMatrix::GetQValueCFTSectorShift(int cftSector)
{
  if ((this->RIndex & 1) == 1)
    return 0;
  if (cftSector == 0)
    return 0;
  return ((this->RIndex + 2 * (this->LaughlinIndex - 1)) / 2);
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int FQHEMPSClustered2RMatrix::GetBondIndexRange(int pLevel, int qValue)
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

int FQHEMPSClustered2RMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
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

int FQHEMPSClustered2RMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
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

int FQHEMPSClustered2RMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{
  return (this->StartingIndexPerPLevelCFTSectorQValue[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]]  + localIndex);
}


// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][0];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][0];
  if (this->NInitialValuePerPLevelCFTSector[pLevel][1] < minQ)
    minQ = this->NInitialValuePerPLevelCFTSector[pLevel][1];
  if (this->NLastValuePerPLevelCFTSector[pLevel][1] > maxQ)
    maxQ = this->NLastValuePerPLevelCFTSector[pLevel][1];  
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  minQ = this->NInitialValuePerPLevelCFTSector[pLevel][cftSector];
  maxQ = this->NLastValuePerPLevelCFTSector[pLevel][cftSector];
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSClustered2RMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  int MinQ;
  int MaxQ;
  this->GetChargeIndexRange(0, MinQ, MaxQ);
  if (padding == true)
    {
//       These lines should be activated in order for FQHEMPSClustered2RMatrix::GetSphereSiteDependentMatrices to work
//       You also need a root configuration such as 11001100
//       rowIndex = 2 * this->PLevel  + 1 - MinQ;
//       columnIndex = rowIndex;
//       return;
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
	      columnIndex = this->PLevel - MinQ;
	    }
	  else
	    {
	      rowIndex = 2 * (this->PLevel + this->RIndex) - MinQ;
	      columnIndex = 2 * this->PLevel - MinQ;
	    }
	}
      else
	{
	  if ((this->RIndex & 1) == 0)
	    {
	      rowIndex = this->PLevel + this->LaughlinIndex - MinQ;
	      columnIndex = this->PLevel - MinQ;
	    }
	  else
	    {
	      rowIndex = 2 * (this->PLevel + this->LaughlinIndex) - MinQ;
	      columnIndex = 2 * this->PLevel - MinQ;
	    }
	}
    }
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool FQHEMPSClustered2RMatrix::LoadHeader (ifstream& file)
{
    int HeaderSize = 0;
    ReadLittleEndian(file, HeaderSize);
    ReadLittleEndian(file, this->PLevel);
    ReadLittleEndian(file, this->LaughlinIndex);
    ReadLittleEndian(file, this->RIndex);
    ReadLittleEndian(file, this->NbrCFTSectors);
    ReadLittleEndian(file, this->NbrNValue);
    ReadLittleEndian(file, this->CylinderFlag);
    ReadLittleEndian(file, this->Kappa);
    ReadLittleEndian(file, this->UniformChargeIndexRange);

    this->NeutralSectorDimension = new int*[this->NbrCFTSectors];
    for (int x = 0; x < this->NbrCFTSectors; ++x)
        this->NeutralSectorDimension[x] = new int[this->PLevel + 1];

    this->U1BasisDimension = new int[this->PLevel + 1];
    this->NbrNValuesPerPLevelCFTSector = new int*[this->PLevel + 1];
    this->NInitialValuePerPLevelCFTSector = new int*[this->PLevel + 1];
    this->NLastValuePerPLevelCFTSector = new int*[this->PLevel + 1];
    this->StartingIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
    this->NbrIndexPerPLevelCFTSectorQValue = new int**[this->PLevel + 1];
    this->StartingIndexPerPLevelCFTSectorQValueU1Sector = new int***[this->PLevel + 1];
    this->NbrIndexPerPLevelCFTSectorQValueU1Sector = new int***[this->PLevel + 1];
    for (int p = 0; p <= this->PLevel; ++p)
    {
        ReadLittleEndian(file, this->U1BasisDimension[p]);

        this->NbrNValuesPerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->NInitialValuePerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->NLastValuePerPLevelCFTSector[p] = new int[this->NbrCFTSectors];
        this->StartingIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
        this->NbrIndexPerPLevelCFTSectorQValue[p] = new int*[this->NbrCFTSectors];
        this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p] = new int**[this->NbrCFTSectors];
        this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p] = new int**[this->NbrCFTSectors];
        for (int x = 0; x < this->NbrCFTSectors; ++x)
        {
            ReadLittleEndian(file, this->NeutralSectorDimension[x][p]);

            ReadLittleEndian(file, this->NbrNValuesPerPLevelCFTSector[p][x]);
            ReadLittleEndian(file, this->NInitialValuePerPLevelCFTSector[p][x]);
            ReadLittleEndian(file, this->NLastValuePerPLevelCFTSector[p][x]);

            this->StartingIndexPerPLevelCFTSectorQValue[p][x] = new int[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->NbrIndexPerPLevelCFTSectorQValue[p][x] = new int[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x] = new int*[this->NbrNValuesPerPLevelCFTSector[p][x]];
            this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x] = new int*[this->NbrNValuesPerPLevelCFTSector[p][x]];
            for (int n = 0; n < this->NbrNValuesPerPLevelCFTSector[p][x]; ++n)
            {
                ReadLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValue[p][x][n]);
                ReadLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValue[p][x][n]);
                this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n] = new int[p + 1];
                this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n] = new int[p + 1];
                for (int k = 0; k <= p; ++k)
                {
                    ReadLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                    ReadLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                }
            }
        }
    }

    return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool FQHEMPSClustered2RMatrix::SaveHeader (ofstream& file)
{
    int HeaderSize = sizeof(int) * 6 + sizeof(char) * 2 + sizeof(double);
    for (int p = 0; p <= this->PLevel; ++p)
    {
        HeaderSize += sizeof(int);
        for (int x = 0; x < this->NbrCFTSectors; ++x)
            HeaderSize += sizeof(int) * (4 + this->NbrNValuesPerPLevelCFTSector[p][x] * (2 + (p + 1) * 2));
    }

    WriteLittleEndian(file, HeaderSize);
    WriteLittleEndian(file, this->PLevel);
    WriteLittleEndian(file, this->LaughlinIndex);
    WriteLittleEndian(file, this->RIndex);
    WriteLittleEndian(file, this->NbrCFTSectors);
    WriteLittleEndian(file, this->NbrNValue);
    WriteLittleEndian(file, this->CylinderFlag);
    WriteLittleEndian(file, this->Kappa);
    WriteLittleEndian(file, this->UniformChargeIndexRange);

    for (int p = 0; p <= this->PLevel; ++p)
    {
        WriteLittleEndian(file, this->U1BasisDimension[p]);

        for (int x = 0; x < this->NbrCFTSectors; ++x)
        {
            WriteLittleEndian(file, this->NeutralSectorDimension[x][p]);

            WriteLittleEndian(file, this->NbrNValuesPerPLevelCFTSector[p][x]);
            WriteLittleEndian(file, this->NInitialValuePerPLevelCFTSector[p][x]);
            WriteLittleEndian(file, this->NLastValuePerPLevelCFTSector[p][x]);
            for (int n = 0; n < this->NbrNValuesPerPLevelCFTSector[p][x]; ++n)
            {
                WriteLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValue[p][x][n]);
                WriteLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValue[p][x][n]);
                for (int k = 0; k <= p; ++k)
                {
                    WriteLittleEndian(file, this->StartingIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                    WriteLittleEndian(file, this->NbrIndexPerPLevelCFTSectorQValueU1Sector[p][x][n][k]);
                }
            }
        }
    }

    return true;
}

// compute the charge index range at a given truncation level
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void FQHEMPSClustered2RMatrix::ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ)
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

// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int FQHEMPSClustered2RMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  if (this->TorusFlag == false)
    {
      nbrFluxQuanta += this->RIndex + 1;
      nbrFluxQuanta *= 2;
      return (nbrFluxQuanta / (this->RIndex + 2));
    }
  else
    {
      return ((nbrFluxQuanta * 2) / (this->RIndex + 2 * (this->LaughlinIndex - 1)));      
    }
}

// compute the scalar product matrix at a given level
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// scalarProductFileName = optional file name of the scalar porduct matrix if the CFT matrix storage is used
// architecture = architecture to use for precalculation
// pLevel = |P| truncation level 
// u1BosonBasis = basis that related the partitions to their index
// rationalScalarProduct = matrices where scalar product matrix elements are stored (rational version)
// scalarProduct = matrices where scalar product matrix elements are stored (double version)
// centralCharge12 = reference on the value of the central charge divided by 12
// centralCharge12Numerical = double accuracy version of centralCharge12
// weight = weight of the primary field that is considered
// weightNumerical = double accuracy version of weight
// sectorName = CFT sector name 
// orthogonalBasisLeft = left transformation matrices related to the complete orthogonal basis 
// orthogonalBasisRight = right transformation matrices related to the complete orthogonal basis 
// rationalMultiplicityFactor = array that contains the multiplicity factors
// multiplicityFactor = double accuracy version of rationalMultiplicityFactor

void FQHEMPSClustered2RMatrix::ComputeFullScalarProductMatrix(const char* cftDirectory, char* scalarProductFileName, AbstractArchitecture* architecture,
							      LongRationalMatrix* rationalScalarProduct, RealSymmetricMatrix* scalarProduct,
							      int pLevel, BosonOnDiskShort** u1BosonBasis, 
							      LongRational& centralCharge12, double centralCharge12Numerical, 
							      LongRational& weight, double weightNumerical,
							      const char* sectorName,
							      RealMatrix* orthogonalBasisLeft, RealMatrix* orthogonalBasisRight,
							      LongRational** rationalMultiplicityFactor, double** multiplicityFactor)
{
  if ((cftDirectory != 0) && (IsFile(scalarProductFileName)))
    {
      if (this->UseRationalFlag == true)
	{
	  rationalScalarProduct[pLevel].ReadMatrix(scalarProductFileName);
	}
      else
	{
	  scalarProduct[pLevel].ReadMatrix(scalarProductFileName);
	}
    }
  else
    {
      if (this->UseRationalFlag == true)
	{
	  FQHEMPSEvaluateCFTOperation Operation1(this, u1BosonBasis, pLevel, centralCharge12, 
						 weight,
						 rationalScalarProduct,  pLevel- 1);
	  Operation1.ApplyOperation(architecture);
	  rationalScalarProduct[pLevel] = Operation1.GetRationalMatrixElements();
	  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
	    {
	      rationalScalarProduct[pLevel].WriteMatrix(scalarProductFileName);
	    }
	}
      else
	{
	  FQHEMPSEvaluateCFTOperation Operation1(this, u1BosonBasis, pLevel, centralCharge12Numerical,
						 weightNumerical,
						 scalarProduct,  pLevel- 1);
	  Operation1.ApplyOperation(architecture);
	  scalarProduct[pLevel] = Operation1.GetOverlapMatrix();
	  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
	    {
	      scalarProduct[pLevel].WriteMatrix(scalarProductFileName);
	    }
	}
    }
  
  RealSymmetricMatrix TmpMatrix;
  if (this->UseRationalFlag == true)
    {
      LongRationalMatrix TmpRationalMatrix(rationalScalarProduct[pLevel].GetNbrRow(), rationalScalarProduct[pLevel].GetNbrColumn());
      for (int k = 0; k < rationalScalarProduct[pLevel].GetNbrRow(); ++k)
	for (int l = 0; l < rationalScalarProduct[pLevel].GetNbrColumn(); ++l)
	  {
	    TmpRationalMatrix[l][k] = rationalScalarProduct[pLevel][l][k] * (rationalMultiplicityFactor[pLevel][k] * rationalMultiplicityFactor[pLevel][l]);
	  }
      TmpMatrix = TmpRationalMatrix;
    }
  else
    {
      TmpMatrix = RealSymmetricMatrix (scalarProduct[pLevel].GetNbrRow(), scalarProduct[pLevel].GetNbrColumn());
      for (int k = 0; k < scalarProduct[pLevel].GetNbrRow(); ++k)
	for (int l = k; l < scalarProduct[pLevel].GetNbrColumn(); ++l)
	  {
	    double Tmp;
	    scalarProduct[pLevel].GetMatrixElement(k, l, Tmp);
	    Tmp *= (multiplicityFactor[pLevel][k] * multiplicityFactor[pLevel][l]);
	    TmpMatrix.SetMatrixElement(k, l, Tmp);
	  }
    }
  
  RealMatrix TmpBasis(u1BosonBasis[pLevel]->GetHilbertSpaceDimension(), u1BosonBasis[pLevel]->GetHilbertSpaceDimension());
  TmpBasis.SetToIdentity();
  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
  TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
  TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
  double Error = 0.0;
  for (int n = 0; n < u1BosonBasis[pLevel]->GetHilbertSpaceDimension(); ++n)
    if (fabs(TmpDiag(n, n)) > Error)
      Error = fabs(TmpDiag(n, n));
  Error *= 1e-14;
  if (Error < 1e-14)
    Error = 1e-14;
  int Count  = 0;
  for (int n = 0; n < u1BosonBasis[pLevel]->GetHilbertSpaceDimension(); ++n)
    {
      if (fabs(TmpDiag(n, n)) < Error)
	++Count;
    }

  int CountNegative = 0;
  for (int n = 0; n < u1BosonBasis[pLevel]->GetHilbertSpaceDimension(); ++n)
      if ((fabs(TmpDiag(n, n)) >= Error) && (TmpDiag(n, n) < 0))
          ++CountNegative;

  cout << "nbr of null vectors " << sectorName << " sector = " << Count << " (" << (u1BosonBasis[pLevel]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << " and " << CountNegative << " negatives" << endl;
  if (Count < u1BosonBasis[pLevel]->GetHilbertSpaceDimension())
    {
      orthogonalBasisLeft[pLevel] = RealMatrix (u1BosonBasis[pLevel]->GetHilbertSpaceDimension(), u1BosonBasis[pLevel]->GetHilbertSpaceDimension() - Count, true);
      orthogonalBasisRight[pLevel] = RealMatrix (u1BosonBasis[pLevel]->GetHilbertSpaceDimension(), u1BosonBasis[pLevel]->GetHilbertSpaceDimension() - Count, true);
      Count = 0;
      for (int n = 0; n < u1BosonBasis[pLevel]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  {
	    orthogonalBasisLeft[pLevel][Count].Copy(TmpBasis[n]);
	    orthogonalBasisRight[pLevel][Count].Copy(TmpBasis[n]);
	    if (TmpDiag(n, n) > 0)
	      {
		orthogonalBasisLeft[pLevel][Count] /=  sqrt(TmpDiag(n, n));
		orthogonalBasisRight[pLevel][Count] /=  sqrt(TmpDiag(n, n));
	      }
	    else
	      {
		orthogonalBasisLeft[pLevel][Count] /=  sqrt(-TmpDiag(n, n));
		orthogonalBasisRight[pLevel][Count] /=  -sqrt(-TmpDiag(n, n));
	      }
	    ++Count;
	  }
    }
  else
    {
      orthogonalBasisLeft[pLevel] = RealMatrix();
      orthogonalBasisRight[pLevel] = RealMatrix();
    }
}

// rescale the scalar product matrix at all levels
//
// rationalScalarProduct = matrices where scalar product matrix elements are stored (rational version)
// scalarProduct = matrices where scalar product matrix elements are stored (double version)
// rationalMultiplicityFactor = array that contains the multiplicity factors
// multiplicityFactor = double accuracy version of rationalMultiplicityFactor

void FQHEMPSClustered2RMatrix::RescaleFullScalarProductMatrix(LongRationalMatrix* rationalScalarProduct, RealSymmetricMatrix* scalarProduct,
							      LongRational** rationalMultiplicityFactor, double** multiplicityFactor)
{
  if (this->NeutralSectorMaxPLevel == -1)
    {
      this->NeutralSectorMaxPLevel = this->PLevel;
    }
  for (int i = 0; i <= this->NeutralSectorMaxPLevel; ++i)
    {
      if (this->UseRationalFlag == true)
 	{
    	  for (int k = 0; k < rationalScalarProduct[i].GetNbrRow(); ++k)
   	    for (int l = 0; l < rationalScalarProduct[i].GetNbrColumn(); ++l)
   	      {
   		rationalScalarProduct[i][l][k] *= (rationalMultiplicityFactor[i][k] * rationalMultiplicityFactor[i][l]);
   	      }
 	  scalarProduct[i] = rationalScalarProduct[i];
  	}
      else
	{
	  for (int k = 0; k < scalarProduct[i].GetNbrRow(); ++k)
	    for (int l = k; l < scalarProduct[i].GetNbrColumn(); ++l)
	      {
		double Tmp;
		scalarProduct[i].GetMatrixElement(k, l, Tmp);
		Tmp *= (multiplicityFactor[i][k] * multiplicityFactor[i][l]);
		scalarProduct[i].SetMatrixElement(k, l, Tmp);
	      }
	}
    }
}
  
// compute the matrix elements of a primary field at a given level
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// matrixElementsFileName = optional file name of the scalar porduct matrix if the CFT matrix storage is used
// architecture = architecture to use for precalculation
// pLevelLeft = |P| truncation level for the left state
// pLevelRight = |P| truncation level for the right state
// u1BosonBasis = basis that related the partitions to their index
// rationalMatrixelements = matrices where matrix elements are stored (rational version)
// matrixElements = matrices where matrix elements are stored (double version)
// centralCharge12 = reference on the value of the central charge divided by 12
// centralCharge12Numerical = double accuracy version of centralCharge12
// weightLeftState = weight of the left state primary field
// weightLeftStateNumerical = double accuracy version of weightLeftState
// weightRightState = weight of the left state primary field
// weightRightStateNumerical = double accuracy version of weightRightState
// weightPrimaryFieldMatrixElement = weight of primary field whose matrix elements are computed
// weightPrimaryFieldMatrixElementNumerical = double accuracy version of weightPrimaryFieldMatrixElement

void FQHEMPSClustered2RMatrix::ComputeFullMatrixElements(const char* cftDirectory, char* matrixElementsFileName, AbstractArchitecture* architecture,
							 LongRationalMatrix** rationalMatrixElements, RealMatrix** matrixElements,
							 int pLevelLeft, int pLevelRight, BosonOnDiskShort** u1BosonBasis, 
							 LongRational& centralCharge12, double centralCharge12Numerical, 
							 LongRational& weightLeftState, double weightLeftStateNumerical,
							 LongRational& weightRightState, double weightRightStateNumerical,
							 LongRational& weightPrimaryFieldMatrixElement, double weightPrimaryFieldMatrixElementNumerical)
{	
  if ((cftDirectory != 0) && (IsFile(matrixElementsFileName)))
    {
      if (this->UseRationalFlag == true)
	{
	  rationalMatrixElements[pLevelLeft][pLevelRight].ReadMatrix(matrixElementsFileName);
	}
      else
	{
	  matrixElements[pLevelLeft][pLevelRight].ReadMatrix(matrixElementsFileName);
	}
    }
  else
    {
      if (this->UseRationalFlag == true)
	{
	  FQHEMPSEvaluateCFTOperation Operation1(this, u1BosonBasis, pLevelLeft, pLevelRight, centralCharge12, 
						 weightLeftState, weightRightState, weightPrimaryFieldMatrixElement,
						 rationalMatrixElements, pLevelLeft - 1, pLevelRight);
	  Operation1.ApplyOperation(architecture);
	  rationalMatrixElements[pLevelLeft][pLevelRight] = Operation1.GetRationalMatrixElements();
	  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
	    {
	      rationalMatrixElements[pLevelLeft][pLevelRight].WriteMatrix(matrixElementsFileName);
	    }
	}
      else
	{
	  FQHEMPSEvaluateCFTOperation Operation1(this, u1BosonBasis, pLevelLeft, pLevelRight, centralCharge12Numerical, 
						 weightLeftStateNumerical, weightRightStateNumerical, weightPrimaryFieldMatrixElementNumerical,
						 matrixElements, pLevelLeft - 1, pLevelRight);
	  Operation1.ApplyOperation(architecture);
	  matrixElements[pLevelLeft][pLevelRight] = Operation1.GetMatrixElements();
	  if  ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
	    {
	      matrixElements[pLevelLeft][pLevelRight].WriteMatrix(matrixElementsFileName);
	    }
	}
    }
}

// rescale the matrix elements at all levels
//
// rationalMatrixElements = matrices where scalar product matrix elements are stored (rational version)
// matrixElements = matrices where scalar product matrix elements are stored (double version)
// rationalMultiplicityFactor = array that contains the multiplicity factors
// multiplicityFactor = double accuracy version of rationalMultiplicityFactor
// globalFactor = global rescaling factor

void FQHEMPSClustered2RMatrix::RescaleFullMatrixElements(LongRationalMatrix** rationalMatrixElements, RealMatrix** matrixElements,
							 LongRational** rationalMultiplicityFactor, double** multiplicityFactor,
							 double globalFactor)
{
  if (this->NeutralSectorMaxPLevel == -1)
    {
      this->NeutralSectorMaxPLevel = this->PLevel;
    }
  for (int j = 0; j <= this->NeutralSectorMaxPLevel; ++j)
    {
      for (int i = 0; i <= this->NeutralSectorMaxPLevel; ++i)
	{
	  if (this->UseRationalFlag == true)
	    {
	      for (int k = 0; k < rationalMatrixElements[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < rationalMatrixElements[i][j].GetNbrColumn(); ++l)
		  {
		    rationalMatrixElements[i][j][l][k] *= (rationalMultiplicityFactor[i][k] * rationalMultiplicityFactor[j][l]);
		  }
	      matrixElements[i][j] = rationalMatrixElements[i][j];
	    }
	  else
	    {
	      for (int k = 0; k < matrixElements[i][j].GetNbrRow(); ++k)
		for (int l = 0; l < matrixElements[i][j].GetNbrColumn(); ++l)
		  {
		    double Tmp;
		    matrixElements[i][j].GetMatrixElement(k, l, Tmp);
		    Tmp *= (multiplicityFactor[i][k] * multiplicityFactor[j][l]);
		    matrixElements[i][j].SetMatrixElement(k, l, Tmp);
		  }
	    }
	  matrixElements[i][j] *= globalFactor;
	}
    }
}

// compute the matrix elements for a primary field
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation
// fieldName = name of the primary field whose matrix elements are computed
// fieldWeight = weight of the primary field whose matrix elements are computed
// nbrSectors = number of conformal families in the basis states
// sectorNames = name of each primary field in the basis states
// weights = weight of each primary field in the basis states
// fusion = OPE structure coefficients, invoked only when fieldWeight != 0
// writeIntermediate = whether to output rational scalar products and matrix elements

void FQHEMPSClustered2RMatrix::ComputeMatrixElements(const char* cftDirectory, AbstractArchitecture* architecture,
						     const char* fieldName, LongRational fieldWeight, int nbrSectors, char** sectorNames, LongRational* weights, RealMatrix fusion, bool writeIntermediate)
{
    LongRational CentralCharge12(this->CentralCharge);
    cout << "central charge = " << CentralCharge12 << endl;
    CentralCharge12 /= 12l;
    double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();
    double fieldWeightNumerical = fieldWeight.GetNumericalValue();

    double* weightsNumerical = new double[nbrSectors];
    for (int s = 0; s < nbrSectors; ++s)
        weightsNumerical[s] = weights[s].GetNumericalValue();

    // descendents are truncated by the total conformal dimension, i.e. "shifted" level
    // all the intermediate matrices are indexed by the actual descendant level
    // only the filenames are indexed by the shifted descendant level
    int* StartingLevels = new int[nbrSectors];
    for (int s = 0; s < nbrSectors; ++s)
        StartingLevels[s] = (weightsNumerical[s] > 0) ? ((int) weightsNumerical[s]) : 0;

    BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
    LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevel + 1];
    double** MultiplicityFactor = new double*[this->PLevel + 1];
    unsigned long* TmpPartition = new unsigned long[this->PLevel + 2];
    for (int i = 0; i <= this->PLevel; ++i)
    {
        U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
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
    delete[] TmpPartition;

    RealSymmetricMatrix** ScalarProducts = new RealSymmetricMatrix*[nbrSectors];
    LongRationalMatrix** RationalScalarProducts = new LongRationalMatrix*[nbrSectors];
    RealMatrix** OrthogonalBasesLeft = new RealMatrix*[nbrSectors];
    RealMatrix** OrthogonalBasesRight = new RealMatrix*[nbrSectors];
    for (int s = 0; s < nbrSectors; ++s)
    {
        ScalarProducts[s] = new RealSymmetricMatrix[this->PLevel + 1];
        RationalScalarProducts[s] = new LongRationalMatrix[this->PLevel + 1];
        OrthogonalBasesLeft[s] = new RealMatrix[this->PLevel + 1];
        OrthogonalBasesRight[s] = new RealMatrix[this->PLevel + 1];
    }

    cout << "computing scalar products" << endl;
    for (int s = 0; s < nbrSectors; ++s)
        cout << sectorNames[s] << ", h = " << weights[s] << endl;

    char* TmpFileName = new char[1024];

    for (int i = 0; i <= this->PLevel; ++i)
      {
        cout << "Level = " <<  i << endl;
        for (int s = 0; s < nbrSectors; ++s)
	  {
            int o = StartingLevels[s];
            if (i + o > this->PLevel)
	      {
                RationalScalarProducts[s][i] = LongRationalMatrix();
                ScalarProducts[s][i] = RealSymmetricMatrix();
                OrthogonalBasesLeft[s][i] = RealMatrix();
                OrthogonalBasesRight[s][i] = RealMatrix();
                continue;
	      }
	    
            if (cftDirectory != 0)
	      {
                if (this->UseRationalFlag == true)
		  sprintf(TmpFileName, "%s/cft_%s_scalarproducts_%s_level_%d.dat", cftDirectory, this->BMatrixOutputName, sectorNames[s], i + o);
                else
		  sprintf(TmpFileName, "%s/cft_%s_num_scalarproducts_%s_level_%d.dat", cftDirectory, this->BMatrixOutputName, sectorNames[s], i + o);
	      }
            this->ComputeFullScalarProductMatrix(cftDirectory, TmpFileName, architecture, RationalScalarProducts[s], ScalarProducts[s], i, U1BosonBasis,
						 CentralCharge12, CentralCharge12Numerical, weights[s], weightsNumerical[s], sectorNames[s],
						 OrthogonalBasesLeft[s], OrthogonalBasesRight[s], RationalMultiplicityFactor, MultiplicityFactor);
	  }
        cout << "---------------------------------" << endl;
      }
    for (int s = 0; s < nbrSectors; ++s)
      this->RescaleFullScalarProductMatrix(RationalScalarProducts[s], ScalarProducts[s], RationalMultiplicityFactor, MultiplicityFactor);

    if (writeIntermediate && this->UseRationalFlag && architecture->CanWriteOnDisk())
      {
        for (int i = 0; i <= this->PLevel; ++i)
	  {
            for (int s = 0; s < nbrSectors; ++s)
	      {
                int o = StartingLevels[s];
                if (i + o > this->PLevel)
		  continue;
		
                sprintf(TmpFileName, "cft_%s_rescaledscalarproducts_%s_level_%d.dat", this->BMatrixOutputName, sectorNames[s], i + o);
                RationalScalarProducts[s][i].WriteMatrix(TmpFileName);
	      }
	  }
      }
    
    for (int i = 0; i <= this->PLevel; ++i)
      for (int s = 0; s < nbrSectors; ++s)
	OrthogonalBasesLeft[s][i].Transpose();
    
    if (fieldWeight == 0) // sandwiched field is identity
      {
        for (int s = 0; s < nbrSectors; ++s)
	  {
            int o = StartingLevels[s];
            for (int i = 0; i <= this->PLevel; ++i)
	      {
                sprintf(TmpFileName, "cft_%s_final_%s_identity_%s_level_%d_%d.dat", this->BMatrixOutputName, sectorNames[s], sectorNames[s], i, i);
                RealMatrix m = (i < o) ? RealMatrix() : (( OrthogonalBasesLeft[s][i - o] *  ((RealMatrix) ScalarProducts[s][i - o])) *  OrthogonalBasesRight[s][i - o]);
                m.WriteMatrix(TmpFileName);
	      }
	  }
      }
    else
      {
        if ((fusion.GetNbrRow() != nbrSectors) || (fusion.GetNbrColumn() != nbrSectors))
	  {
            cout << "fusion matrix size mismatch!" << endl;
            exit(1);
	  }
	
        int nbrChannels = 0;
        for (int l = 0; l < nbrSectors; ++l)
	  for (int r = 0; r < nbrSectors; ++r)
	    if (fusion[r][l] != 0)
	      ++nbrChannels;

        int* leftSectors = new int[nbrChannels];
        int* rightSectors = new int[nbrChannels];
        double* globalFactors = new double[nbrChannels];
	
        int channel = 0;
        for (int l = 0; l < nbrSectors; ++l)
	  {
            for (int r = 0; r < nbrSectors; ++r)
	      {
                if (fusion[r][l] != 0)
		  {
                    leftSectors[channel] = l;
                    rightSectors[channel] = r;
                    globalFactors[channel] = fusion[r][l];
                    ++channel;
                }
	      }
	  }

        cout << "computing matrix elements for " << fieldName << ", with h = " << fieldWeight << ", in channels:" << endl;
        for (int c = 0; c < nbrChannels; ++c)
            cout << " <" << sectorNames[leftSectors[c]] << "|" << fieldName << "|" << sectorNames[rightSectors[c]] << ">";
        cout << endl;

        RealMatrix*** Matrices = new RealMatrix**[nbrChannels];
        LongRationalMatrix*** RationalMatrices = new LongRationalMatrix**[nbrChannels];
        for (int c = 0; c < nbrChannels; ++c)
        {
            Matrices[c] = new RealMatrix*[this->PLevel + 1];
            RationalMatrices[c] = new LongRationalMatrix*[this->PLevel + 1];

            for (int i = 0; i <= this->PLevel; ++i)
            {
                Matrices[c][i] = new RealMatrix[this->PLevel + 1];
                RationalMatrices[c][i] = new LongRationalMatrix[this->PLevel + 1];
            }
        }

        for (int j = 0; j <= this->PLevel; ++j)
        {
            for (int i = 0; i <= this->PLevel; ++i)
            {
                cout << "Levels = " <<  i << " " << j << endl;
                for (int c = 0; c < nbrChannels; ++c)
                {
                    int l = StartingLevels[leftSectors[c]];
                    int r = StartingLevels[rightSectors[c]];
                    if ((i + l > this->PLevel) || (j + r > this->PLevel))
                    {
                        RationalMatrices[c][i][j] = LongRationalMatrix();
                        Matrices[c][i][j] = RealMatrix();
                        continue;
                    }

                    if (cftDirectory != 0)
                    {
                        if (this->UseRationalFlag == true)
                            sprintf(TmpFileName, "%s/cft_%s_matrixelement_%s_%s_%s_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, sectorNames[leftSectors[c]], fieldName, sectorNames[rightSectors[c]], i + l, j + r);
                        else
                            sprintf(TmpFileName, "%s/cft_%s_num_matrixelement_%s_%s_%s_level_%d_%d.dat", cftDirectory, this->BMatrixOutputName, sectorNames[leftSectors[c]], fieldName, sectorNames[rightSectors[c]], i + l, j + r);
                    }

                    this->ComputeFullMatrixElements(cftDirectory, TmpFileName, architecture,
                            RationalMatrices[c], Matrices[c], i, j, U1BosonBasis,
                            CentralCharge12, CentralCharge12Numerical,
                            weights[leftSectors[c]], weightsNumerical[leftSectors[c]],
                            weights[rightSectors[c]], weightsNumerical[rightSectors[c]],
                            fieldWeight, fieldWeightNumerical);
                }
            }
        }
        for (int c = 0; c < nbrChannels; ++c)
            this->RescaleFullMatrixElements(RationalMatrices[c], Matrices[c], RationalMultiplicityFactor, MultiplicityFactor, globalFactors[c]);

        if (writeIntermediate && this->UseRationalFlag && architecture->CanWriteOnDisk())
	  {
	    for (int j = 0; j <= this->PLevel; ++j)
	      {
                for (int i = 0; i <= this->PLevel; ++i)
		  {
                    for (int c = 0; c < nbrChannels; ++c)
                    {
		      int l = StartingLevels[leftSectors[c]];
		      int r = StartingLevels[rightSectors[c]];
		      if ((i + l > this->PLevel) || (j + r > this->PLevel))
			continue;
		      
		      sprintf(TmpFileName, "cft_%s_rescaledmatrixelement_%s_%s_%s_level_%d_%d.dat", this->BMatrixOutputName, sectorNames[leftSectors[c]], fieldName, sectorNames[rightSectors[c]], i + l, j + r);
		      RationalMatrices[c][i][j].WriteMatrix(TmpFileName); // OPE structure constants (globalFactors) are not included in the RationalMatrices
                    }
		  }
	      }
	  }
	
        for (int c = 0; c < nbrChannels; ++c)
	  {
            int l = StartingLevels[leftSectors[c]];
            int r = StartingLevels[rightSectors[c]];
	    
            for (int i = 0; i <= this->PLevel; ++i)
	      {
                for (int j = 0; j <= this->PLevel; ++j)
		  {
                    sprintf(TmpFileName, "cft_%s_final_%s_%s_%s_level_%d_%d.dat", this->BMatrixOutputName, sectorNames[leftSectors[c]], fieldName, sectorNames[rightSectors[c]], i, j);
                    RealMatrix m = ((i < l) || (j < r)) ? RealMatrix() : ((OrthogonalBasesLeft[leftSectors[c]][i - l] * Matrices[c][i - l][j - r]) * OrthogonalBasesRight[rightSectors[c]][j - r]);
                    m.WriteMatrix(TmpFileName);
		  }
	      }
	  }
	
        for (int c = 0; c < nbrChannels; ++c)
	  {
            for (int i = 0; i <= this->PLevel; ++i)
	      {
                delete[] RationalMatrices[c][i];
                delete[] Matrices[c][i];
	      }
            delete[] RationalMatrices[c];
            delete[] Matrices[c];
	  }
        delete[] RationalMatrices;
        delete[] Matrices;

        delete[] leftSectors;
        delete[] rightSectors;
        delete[] globalFactors;
    }
    
    delete[] TmpFileName;
    
    delete[] weightsNumerical;
    for (int i = 0; i <= this->PLevel; ++i)
      {
        delete U1BosonBasis[i];
        delete[] RationalMultiplicityFactor[i];
        delete[] MultiplicityFactor[i];
      }
    delete[] U1BosonBasis;
    delete[] RationalMultiplicityFactor;
    delete[] MultiplicityFactor;
    
    for (int s = 0; s < nbrSectors; ++s)
      {
        delete[] ScalarProducts[s];
        delete[] RationalScalarProducts[s];
        delete[] OrthogonalBasesLeft[s];
        delete[] OrthogonalBasesRight[s];
      }
    delete[] ScalarProducts;
    delete[] RationalScalarProducts;
    delete[] OrthogonalBasesLeft;
    delete[] OrthogonalBasesRight;

    delete[] StartingLevels;
}

// get the matrix that into account the Jordan Wigner string on the torus geometry
//
// nbrFermions = number of fermions in the system
// return value = corresponding matrix

SparseRealMatrix FQHEMPSClustered2RMatrix::GetTorusStringMatrix(int nbrFermions)
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

int* FQHEMPSClustered2RMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
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
	      if (((CurrentQValue + this->GetQValueCFTSectorShift(CurrentCFTSector)) % (this->RIndex + 2 * (this->LaughlinIndex - 1))) == topologicalSector)
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
	      if (((CurrentQValue + this->GetQValueCFTSectorShift(CurrentCFTSector)) % (this->RIndex + 2 * (this->LaughlinIndex - 1))) == topologicalSector)
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
  
// print a given state of the auxiliary space
//
// str = reference on the output stream
// index = index of the state
// return value = reference on the output stream

ostream& FQHEMPSClustered2RMatrix::PrintAuxiliarySpaceState(ostream& str, int index)
{
  int TmpPLevel;
  int TmpQ;
  int TmpSector;
  int TmpChargeSectorLevel;
  int TmpNeutralSectorLevel;
  int TmpChargeSectorIndex;
  int TmpNeutralSectorIndex;
  this->GetCFTSectorChargeAndPLevelFromMatrixIndex(index, TmpSector, TmpPLevel, TmpQ, TmpChargeSectorLevel, TmpNeutralSectorLevel, 
						   TmpChargeSectorIndex, TmpNeutralSectorIndex);
  str << "|" << index << ": x=" << TmpSector << " Q=" << TmpQ << " P=" << TmpPLevel 
      << " : P_mu=" << TmpChargeSectorLevel << " mu=" << TmpChargeSectorIndex 
      << " : P_l=" << TmpNeutralSectorLevel << " l=" << TmpNeutralSectorIndex << ">";
  return str;
}

// get the label associated to a given state of the auxiliary space
//
// index = auxiliary space index
// return value = string containing the label

char* FQHEMPSClustered2RMatrix::GetAuxiliarySpaceLabel(int index)
{
  int TmpPLevel;
  int TmpQ;
  int TmpCFTSector;
  int TmpChargeSectorLevel;
  int TmpNeutralSectorLevel;
  int TmpChargeSectorIndex;
  int TmpNeutralSectorIndex;
  this->GetCFTSectorChargeAndPLevelFromMatrixIndex(index, TmpCFTSector, TmpPLevel, TmpQ, TmpChargeSectorLevel, TmpNeutralSectorLevel, 
						   TmpChargeSectorIndex, TmpNeutralSectorIndex);
  char* TmpString = new char [256];
  sprintf (TmpString, "(x=%d, Q=%d, P=%d, P_mu=%d, mu=%d, P_l=%d, l=%d, i=%d)", TmpCFTSector, TmpQ, TmpPLevel, 
	   TmpChargeSectorLevel, TmpChargeSectorIndex, TmpNeutralSectorLevel, TmpNeutralSectorIndex, index);
  return TmpString;
}

// get the array where the site-dependent matrices for the geometry are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseRealMatrix** FQHEMPSClustered2RMatrix::GetSphereSiteDependentMatrices(int nbrFluxQuanta)
{
  // this code is not intended to be efficient but the code is closer to the PH Pfaffian code
  
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
  if (this->CFTDirectory != 0)
    {
      TmpScalarProductIdentityFileName = new char[512 + strlen(this->CFTDirectory)];
      TmpScalarProductPsiFileName = new char[512 + strlen(this->CFTDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
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

  QValue = this->LaughlinIndex - 1 + (this->RIndex / 2);
  cout << "this->LaughlinIndex = " << this->LaughlinIndex  << " " << QValue << endl;
// PH-Pfaffian values
//  this->NbrNValue = 4 * this->PLevel + 3 + this->RIndex;
//  NValueShift = 4 * this->PLevel + 2 + this->LaughlinIndex;
// Pfaffian values
  this->NbrNValue = 4 * this->PLevel + 5;
  NValueShift = 4 * this->PLevel + 2;
  QValueDenominator = 1;

  int MatrixSize = this->ComputeLinearizedIndexArrays();
  cout << "B matrix size = " << MatrixSize << endl;

  cout << "computing Psi matrix elements" << endl;
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
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
  for (int n = 0; n <= this->PLevel; ++n)
    {
      for (int p = 0; p <= this->PLevel; ++p)
	{
	  if (((p + n) <= this->PLevel)  && ((p + n) >= 0))
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
		  for (int j = this->NInitialValuePerPLevelCFTSector[p + n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[p + n][0]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(p + n, 0, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[p + n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[p + n][1]; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(p + n, 1, j - 1, p, ChargedIndex, NeutralIndex1)];
			    }
			}
		    }
		}
	    }
	}
    }

  BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  TmpB0Matrix = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);

  for (int n = 0; n <= this->PLevel; ++n)
    {
      for (int p = 0; p <= this->PLevel; ++p)
	{
	  if (((p + n) <= this->PLevel)  && ((p + n) >= 0))
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
		  for (int j = this->NInitialValuePerPLevelCFTSector[p + n][0] + 1; j <= this->NLastValuePerPLevelCFTSector[p + n][0]; ++j)
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
			      TmpB0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(p + n, 0, j - 1, p, ChargedIndex, NeutralIndex1),
							    this->Get2RMatrixIndexV2(p + n, 0, j, p, ChargedIndex, NeutralIndex2), Tmp);
			    }
			}
		    }
		  for (int j = this->NInitialValuePerPLevelCFTSector[p + n][1] + 1; j <= this->NLastValuePerPLevelCFTSector[p + n][1]; ++j)
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
			      TmpB0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(p + n, 1, j - 1, p, ChargedIndex, NeutralIndex1),
							    this->Get2RMatrixIndexV2(p + n, 1, j, p, ChargedIndex, NeutralIndex2), Tmp);
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

//   for (int i = 0; i < MatrixSize; ++i)
//     {
//       this->PrintAuxiliarySpaceState(cout, i) << endl;
//     }

  int NbrV0MatrixIndices = (2 * this->PLevel) + 2;
// PH-Pfaffian values
//   int V0MatrixIndexShift = this->PLevelShift + 1;
// Pfaffian values
  int V0MatrixIndexShift = this->PLevel;
  SparseRealMatrix* V0Matrices = new SparseRealMatrix[NbrV0MatrixIndices];

  for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
    {
      for (int i = 0; i < MatrixSize; ++i)
	TmpNbrElementPerRow[i] = 0;
      long TmpTotalNbrElements = 0l;
      for (int n = 0; n <= this->PLevel; ++n)
	{
	  for (int p = 0; p <= this->PLevel; ++p)
	    {
	      if (((p + n) <= this->PLevel) && ((p + n) >= 0))
		{
		  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
		  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
		  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
		  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
		  for (int m = 0; m <= this->PLevel; ++m)
		    {
		      for (int q = 0; q <= this->PLevel; ++q)
			{
			  if (((q + m) <= this->PLevel) && ((q + m) >= 0))
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
				      N2 = (2 * ((q + m) - (p + n))  + 2 + NValueShift) / 2;
				      N1 = N2 + QValue;
				      if (((N1 >= this->NInitialValuePerPLevelCFTSector[p + n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[p + n][0]))
					  && ((N2 >= this->NInitialValuePerPLevelCFTSector[q + m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[q + m][1]))
					  && ((V0MatrixIndex - V0MatrixIndexShift) == (m - n + 1)))
					{ 
					  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
					    {
					      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
						{
						  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(p + n, 0, N1, p, ChargedIndex1, NeutralIndex1)];
						  ++TmpTotalNbrElements;
						}
					    }
					  
					}
				      N2 = (2 * ((q + m) - (p + n)) + NValueShift) / 2;
				      N1 = N2 + QValue;
				      if (((N1 >= this->NInitialValuePerPLevelCFTSector[p + n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[p + n][1]))
					  && ((N2 >= this->NInitialValuePerPLevelCFTSector[q + m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[q + m][0]))
					  && ((V0MatrixIndex - V0MatrixIndexShift) == (m - n)))
					{ 
					  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
					    {
					      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
						{
						  ++TmpNbrElementPerRow[this->Get2RMatrixIndexV2(p + n, 1, N1, p, ChargedIndex1, NeutralIndex1)];
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
	  for (int n = 0; n <= this->PLevel; ++n)
	    {
	      for (int p = 0; p <= this->PLevel; ++p)
		{
		  if (((p + n) <= this->PLevel) && ((p + n) >= 0))
		    {
		      BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
		      BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[n];
		      RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[n];
		      RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[n];
		      for (int m = 0; m <= this->PLevel; ++m)
			{
			  for (int q = 0; q <= this->PLevel; ++q)
			    {
			      if (((q + m) <= this->PLevel) && ((q + m) >= 0))
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
					  N2 = (2 * ((q + m) - (p + n)) + 2 + NValueShift) / 2;
					  N1 = N2 + QValue;
					  if (((N1 >= this->NInitialValuePerPLevelCFTSector[p + n][0]) && (N1 <= this->NLastValuePerPLevelCFTSector[p + n][0]))
					      && ((N2 >= this->NInitialValuePerPLevelCFTSector[q + m][1]) && (N2 <= this->NLastValuePerPLevelCFTSector[q + m][1]))
					      && ((V0MatrixIndex - V0MatrixIndexShift) == (m - n + 1)))
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
						      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(p + n, 0, N1, p, ChargedIndex1, NeutralIndex1),
										this->Get2RMatrixIndexV2(q + m, 1, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
						    }
						  
						}
					    }
					  N2 = (2 * ((q + m) - (p + n)) + NValueShift) / 2;
					  N1 = N2 + QValue;
					  if (((N1 >= this->NInitialValuePerPLevelCFTSector[p + n][1]) && (N1 <= this->NLastValuePerPLevelCFTSector[p + n][1]))
					      && ((N2 >= this->NInitialValuePerPLevelCFTSector[q + m][0]) && (N2 <= this->NLastValuePerPLevelCFTSector[q + m][0]))
					      && ((V0MatrixIndex - V0MatrixIndexShift) == (m - n)))
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
						      V0Matrix.SetMatrixElement(this->Get2RMatrixIndexV2(p + n, 1, N1, p, ChargedIndex1, NeutralIndex1),
										this->Get2RMatrixIndexV2(q + m, 0, N2, q, ChargedIndex2, NeutralIndex2), Tmp);
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
//  	  cout << "V0[" << V0MatrixIndex << "] = " << endl;
//  	  V0Matrices[V0MatrixIndex].PrintNonZero(cout, TmpLabels, TmpLabels) << endl;
	}
    }

  cout.precision(14);
  double** TmpProjectorCoefficients = new double*[nbrFluxQuanta + 1];
  BinomialCoefficients TmpCoef (nbrFluxQuanta);
  for (int i = 0; i <= nbrFluxQuanta; ++i)
    {
      TmpProjectorCoefficients[i] = new double[NbrV0MatrixIndices];
      for (int V0MatrixIndex = 0; V0MatrixIndex < NbrV0MatrixIndices; ++V0MatrixIndex)
	{
	  TmpProjectorCoefficients[i][V0MatrixIndex] = 1.0 / (sqrt(4.0 * M_PI / ((double) (nbrFluxQuanta + 1)) * TmpCoef(nbrFluxQuanta, i)));
	}
    }


  SparseRealMatrix** TmpMatrices = new SparseRealMatrix*[nbrFluxQuanta + 1];
  SparseRealMatrix* TmpMatrices2 = new SparseRealMatrix[NbrV0MatrixIndices];
  double* TmpCoefficients = new double[NbrV0MatrixIndices];
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
      if (TmpNbrV0Matrices > 0)
	{
	  SparseRealMatrix TmpV0Matrix = SparseRealMatrixLinearCombination(TmpNbrV0Matrices, TmpCoefficients, TmpMatrices2);
	  for (int m = 1; m < this->NbrBMatrices; ++m)
	    {
	      TmpMatrices[i][m] = MemoryEfficientMultiply(TmpMatrices[i][m - 1], TmpV0Matrix);
	    }
//  	  cout << "B[1," << i << "] = " << endl;
//   	  TmpMatrices[i][1].PrintNonZero(cout, TmpLabels, TmpLabels) << endl;
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
  return TmpMatrices;
}

