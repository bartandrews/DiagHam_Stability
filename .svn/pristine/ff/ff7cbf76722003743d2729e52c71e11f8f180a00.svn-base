////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of qhe on disk main task                      //
//                                                                            //
//                        last modification : 08/10/2004                      //
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
#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Vector/RealVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicBlockLanczosAlgorithm.h"

#include "Options/Options.h"

#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// constructor
//  
// options = pointer to the options managers containing all running options
// space = pointer to the current Hilbert space
// lanczos = pointer to the Lanczos manager
// hamiltonian = pointer to the current Hamiltonian
// kyValue = total momentum value of the system along the y-axis
// shift = energy shift that is applied to the hamiltonian
// outputFileName = name of the file where results have to be stored
// firstRun = flag that indicates if it the first time the main task is used
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector
// kxValue = set the Kx value (-1 if the hamiltonian does not handle the Kx symmetry)
// explicitInitialVector = an optional pointer to an initial vector to be used in the Lanczos run, overriding command line arguments
// forceReal = assume that the hamiltonian is real even for kx>=0 (usually at the high symmetry points)

FQHEOnTorusMainTask::FQHEOnTorusMainTask(OptionManager* options, AbstractHilbertSpace* space, LanczosManager* lanczos, 
					 AbstractQHEHamiltonian* hamiltonian, int kyValue, double shift, char* outputFileName,
					 bool firstRun, char* eigenvectorFileName, int kxValue, Vector *explicitInitialVector, bool forceReal)
{
  this->OutputFileName = new char [strlen(outputFileName) + 1];
  strncpy(this->OutputFileName, outputFileName, strlen(outputFileName));
  this->OutputFileName[strlen(outputFileName)] = '\0';
  if (eigenvectorFileName == 0)
    {
      this->EigenvectorFileName = 0;
    }
  else
    {
      this->EigenvectorFileName = new char [strlen(eigenvectorFileName) + 1];
      strncpy(this->EigenvectorFileName, eigenvectorFileName, strlen(eigenvectorFileName));
      this->EigenvectorFileName[strlen(eigenvectorFileName)] = '\0';
    }
  this->Hamiltonian = hamiltonian;
  this->Space = space;
  if (this->Hamiltonian->GetHilbertSpaceDimension() !=  this->Space->GetHilbertSpaceDimension())
    {
      printf("Error: Hamiltonian does not match Space dimension\n");
      exit(1);
    }
  this->KyValue = kyValue;
  this->KxValue = -1;
  this->KyOnlyFlag = true;
  this->RealFlag = true;
  if (kxValue >= 0)
    {
      this->KxValue = kxValue;
      this->KyOnlyFlag = false;
      this->RealFlag = forceReal;      
    }
  //if (explicitInitialVector->GetVectorType()&Vector::ComplexDatas)
  //  this->RealFlag = false;
  this->MultiplicityFlag = false;
  this->AlgorithmManager = lanczos;
  this->EnergyShift = shift;
  this->ResumeFlag = ((BooleanOption*) (*options)["resume"])->GetBoolean();
  this->DiskFlag = ((BooleanOption*) (*options)["disk"])->GetBoolean();
  this->MaxNbrIterLanczos = ((SingleIntegerOption*) (*options)["iter-max"])->GetInteger();
  this->NbrIterLanczos = ((SingleIntegerOption*) (*options)["nbr-iter"])->GetInteger();
  this->NbrEigenvalue = ((SingleIntegerOption*) (*options)["nbr-eigen"])->GetInteger();
  if (this->NbrEigenvalue > this->Space->GetHilbertSpaceDimension())
    {
      this->NbrEigenvalue = this->Space->GetHilbertSpaceDimension();
    }
  this->FullDiagonalizationLimit = ((SingleIntegerOption*) (*options)["full-diag"])->GetInteger();
  this->BlockLanczosFlag = false;
  if ((*options)["block-lanczos"] != 0)
    {
      this->BlockLanczosFlag = ((BooleanOption*) (*options)["block-lanczos"])->GetBoolean();
    }
  this->SizeBlockLanczos = 1;
  if ((*options)["block-size"] != 0)
    {
      this->SizeBlockLanczos = ((SingleIntegerOption*) (*options)["block-size"])->GetInteger();
    }
  this->VectorMemory = ((SingleIntegerOption*) (*options)["nbr-vector"])->GetInteger();
  if ((*options)["save-precalculation"] != 0)
    {
      this->SavePrecalculationFileName = ((SingleStringOption*) (*options)["save-precalculation"])->GetString();
    }
    else
      SavePrecalculationFileName=0;
  this->FullReorthogonalizationFlag = ((BooleanOption*) (*options)["force-reorthogonalize"])->GetBoolean();
  this->EvaluateEigenvectors = ((BooleanOption*) (*options)["eigenstate"])->GetBoolean();
  this->EigenvectorConvergence = ((BooleanOption*) (*options)["eigenstate-convergence"])->GetBoolean();
  if ((*options)["show-itertime"] != 0)
    {
      this->ShowIterationTime = ((BooleanOption*) (*options)["show-itertime"])->GetBoolean();
    }
  else
    this->ShowIterationTime = false;
  if ((*options)["initial-vector"] != 0)
    {
      this->InitialVectorFileName = ((SingleStringOption*) (*options)["initial-vector"])->GetString();
    }
  else
    {
      this->InitialVectorFileName = 0;
    }
  ExplicitInitialVector=explicitInitialVector;
  if ((*options)["initial-blockvectors"] != 0)
    {
      this->InitialBlockVectorFileName = ((SingleStringOption*) (*options)["initial-blockvectors"])->GetString();
    }
  else
    {
      this->InitialBlockVectorFileName = 0;
    }
  if ((*options)["partial-lanczos"] != 0)
    {
      this->PartialLanczos = options->GetBoolean("partial-lanczos");
    }
  else
    {
      this->PartialLanczos = false;
    }
    
  if ((*options)["sr-save-interval"] != 0)
    {
      this->SpectralResponseSaveInterval = options->GetInteger("sr-save-interval");
    }
  else
    {
      this->SpectralResponseSaveInterval = 0;
    }
  if ((*options)["sr-omega-min"] != 0)
    {
      this-> SpectralResponseOmegaMin = options->GetDouble("sr-omega-min");
    }
  else
    {
      this->SpectralResponseOmegaMin = 0.0;
    }
  if ((*options)["sr-omega-max"] != 0)
    {
      this-> SpectralResponseOmegaMax = options->GetDouble("sr-omega-max");
    }
  else
    {
      this->SpectralResponseOmegaMax = 0.0;
    }
  if ((*options)["sr-epsilon"] != 0)
    {
      this-> SpectralResponseEpsilon = options->GetDouble("sr-epsilon");
    }
  else
    {
      this->SpectralResponseEpsilon = 0.0;
    }
  if ((*options)["sr-omega-interval"] != 0)
    {
      this-> SpectralResponseOmegaInterval = options->GetDouble("sr-omega-interval");
    }
  else
    {
      this->SpectralResponseOmegaInterval = 0.0;
    }
  if ((*options)["sr-spectral-resolution"] != 0)
    {
      this-> SpectralResponseSpectralResolution = options->GetDouble("sr-spectral-resolution");
    }
  else
    {
      this->SpectralResponseSpectralResolution = 0.0;
    }
    
  if ((*options)["use-lapack"] != 0)
    {
      this->LapackFlag = ((BooleanOption*) (*options)["use-lapack"])->GetBoolean();
    }
  else
    {
      this->LapackFlag = false;
    }
  if ((*options)["limit-time"] != 0)
    {
      this->MaximumAllowedTime = (((SingleIntegerOption*) (*options)["limit-time"])->GetInteger());
    }
  else
    {
      this->MaximumAllowedTime = 0;
    }
  if ((((*options)["use-hilbert"]) != 0) && (((SingleStringOption*) (*options)["use-hilbert"])->GetString() != 0))
    {
      this->ReducedHilbertSpaceDescription = ((SingleStringOption*) (*options)["use-hilbert"])->GetString();
      if (((*options)["export-hilberttransformation"] != 0) &&  (options->GetString("export-hilberttransformation") != 0))
	{
	  this->ReducedHilbertSpaceExportTransformation = options->GetString("export-hilberttransformation");
	}
      else
	{
	  this->ReducedHilbertSpaceExportTransformation = 0;
	}
    }
  else
    {
      this->ReducedHilbertSpaceDescription = 0;
      this->ReducedHilbertSpaceExportTransformation = 0;
    }
  if ((*options)["get-hvalue"] != 0)
    {
      this->ComputeEnergyFlag = ((BooleanOption*) (*options)["get-hvalue"])->GetBoolean();
    }
  else
    {
      this->ComputeEnergyFlag = false;
    }
  
  this->ShowHamiltonian = false;
  if (((*options)["show-hamiltonian"] != 0) && (((BooleanOption*) (*options)["show-hamiltonian"])->GetBoolean() == true))
    {
      this->ShowHamiltonian = true;
      if (this->ReducedHilbertSpaceDescription == 0)
	{
	  if (RealFlag)  
	    {
	      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
	      this->Hamiltonian->GetHamiltonian(HRep);
	      cout << HRep << endl;
	    }
	  else
	    {
	      ComplexMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension(), true);
	      this->Hamiltonian->GetHamiltonian(HRep);
	      cout << HRep << endl;
	    } 
	}
    }
  if (((*options)["friendlyshow-hamiltonian"] != 0) && (options->GetBoolean("friendlyshow-hamiltonian") == true))
    {
      this->ShowHamiltonian = true;
      if (this->ReducedHilbertSpaceDescription == 0)
	{
	  if (RealFlag)  
	    {	  
	      RealMatrix HRep (this->Space->GetHilbertSpaceDimension(), this->Space->GetHilbertSpaceDimension());
	      this->Hamiltonian->GetHamiltonian(HRep);
	      for (int i = 0; i < this->Space->GetHilbertSpaceDimension(); ++i)
		{
		  if (HRep[i].Norm() != 0.0)
		    {
		      cout << i << " : ";
		      this->Space->PrintState(cout, i) << endl;
		      for (int j = 0; j < this->Space->GetHilbertSpaceDimension(); ++j)
			{
			  if (HRep[i][j] != 0.0)
			    {
			      cout << "    " << j << " : ";
			      this->Space->PrintState(cout, j) << " : " << HRep[i][j] << endl;
			    }
			}
		    }
		}
	    }
	  else
	    {	  
	      ComplexMatrix HRep (this->Space->GetHilbertSpaceDimension(), this->Space->GetHilbertSpaceDimension());
	      this->Hamiltonian->GetHamiltonian(HRep);
	      for (int i = 0; i < this->Space->GetHilbertSpaceDimension(); ++i)
		{
		  if (HRep[i].Norm() != 0.0)
		    {
		      cout << i << " : ";
		      this->Space->PrintState(cout, i) << endl;
		      for (int j = 0; j < this->Space->GetHilbertSpaceDimension(); ++j)
			{
			  if (Norm(HRep[i][j]) != 0.0)
			    {
			      cout << "    " << j << " : ";
			      this->Space->PrintState(cout, j) << " : " << HRep[i][j] << endl;
			    }
			}
		    }
		}
	    }
	}
    }  
  if (((*options)["test-hermitian"] != 0) && (options->GetBoolean("test-hermitian") == true))
    {
      if (RealFlag)  
	{	 
	  RealMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
	  this->Hamiltonian->GetHamiltonian(HRep);
	  double Tmp1;
	  double Tmp2;
	  cout << "check hermiticity" << endl;
	  double AverageNorm = 0.0;
	  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    for (int j = i; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
	      {
		HRep.GetMatrixElement(i, j, Tmp1);
		AverageNorm += fabs(Tmp1);
	      }
	  AverageNorm /= 0.5 * ((double) this->Hamiltonian->GetHilbertSpaceDimension()) * ((double) (this->Hamiltonian->GetHilbertSpaceDimension() + 1));
	  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    for (int j = i; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
	      {
		HRep.GetMatrixElement(i, j, Tmp1);
		HRep.GetMatrixElement(j, i, Tmp2);
		if (fabs(Tmp1 - Tmp2) > (MACHINE_PRECISION * AverageNorm))
		  {
		    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << endl;
		  }
	      }
	  cout << "check done" << endl;
	}
      else
	{ 
	  ComplexMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
	  this->Hamiltonian->GetHamiltonian(HRep);
	  Complex Tmp1;
	  Complex Tmp2;
	  cout << "check hermiticity" << endl;
	  double AverageNorm = 0.0;
	  double Error = MACHINE_PRECISION;
	  if (((*options)["testhermitian-error"] != 0) && (options->GetDouble("testhermitian-error") != 0.0))
	    {
	      Error = options->GetDouble("testhermitian-error");
	    }
	  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    for (int j = i; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
	      {
		HRep.GetMatrixElement(i, j, Tmp1);
		AverageNorm += Norm(Tmp1);
	      }
	  AverageNorm /= 0.5 * ((double) this->Hamiltonian->GetHilbertSpaceDimension()) * ((double) (this->Hamiltonian->GetHilbertSpaceDimension() + 1));
	  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
	    for (int j = i; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
	      {
		HRep.GetMatrixElement(i, j, Tmp1);
		HRep.GetMatrixElement(j, i, Tmp2);
		if (Norm(Tmp1 - Conj(Tmp2)) > (Error * AverageNorm))
		  {
		    cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << " " << Norm(Tmp1 - Conj(Tmp2)) << " (should be lower than " << (Error * AverageNorm) << ")" << endl;
		  }
	      }
	  cout << "check done" << endl;
	}
    }
  if (((*options)["lanczos-precision"] != 0) && (((SingleDoubleOption*) (*options)["lanczos-precision"])->GetDouble() > 0))
    {
      this->LanczosPrecision = ((SingleDoubleOption*) (*options)["lanczos-precision"])->GetDouble();
    }
  else
    {
      this->LanczosPrecision = 0.0;
    }

  if (((*options)["fast-disk"] != 0) && (this->EvaluateEigenvectors == true))
    {
      this->FastDiskFlag = ((BooleanOption*) (*options)["fast-disk"])->GetBoolean();
      if ((*options)["resume-fastdisk"] != 0)
	{
	  this->ResumeFastDiskFlag = ((BooleanOption*) (*options)["resume-fastdisk"])->GetBoolean();
	}
    }
  else
    {
      this->FastDiskFlag = false;
      this->ResumeFastDiskFlag = false;
    }

  if ((((*options)["set-reorthogonalize"]) != 0) && (((SingleStringOption*) (*options)["set-reorthogonalize"])->GetString() != 0))
    {
      this->LanczosReorthogonalization = ((SingleStringOption*) (*options)["set-reorthogonalize"])->GetString();
    }
  else
    {
      this->LanczosReorthogonalization = 0;
    }
  this->FirstRun = firstRun;
}  
 
// destructor
//  

FQHEOnTorusMainTask::~FQHEOnTorusMainTask()
{
}

// execute the main task
// 
// return value = 0 if no error occurs, else return error code

int FQHEOnTorusMainTask::ExecuteMainTask()
{
  ofstream File;
  if (KyOnlyFlag)
    {
      if (this->FirstRun == true)
	{
	  File.open(this->OutputFileName, ios::binary | ios::out);
	  this->FirstRun = false;
	  File << "# Ky E";
	  if ((this->EvaluateEigenvectors == true) && (this->ComputeEnergyFlag == true))
	    File << " <H>";
	  File << endl;
	}
      else
	{
	  File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
	}
      File.precision(14);
      cout.precision(14);
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ky = " << this->KyValue << endl;
      cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
    }
  else
    {
      if (this->FirstRun == true)
	{
	  File.open(this->OutputFileName, ios::binary | ios::out);
	  this->FirstRun = false;
	  File << "# Kx Ky E";
	  if ((this->EvaluateEigenvectors == true) && (this->ComputeEnergyFlag == true))
	    File << " <H>";
	  File << endl;
	}
      else
	{
	  File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
	}
      File.precision(14);
      cout.precision(14);
      cout << "----------------------------------------------------------------" << endl;
      cout << " Kx = " << this->KxValue << " Ky = " << this->KyValue << endl;
      cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
    }
      
  if (RealFlag)
    {
      if (this->SavePrecalculationFileName != 0)
	{
	  this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
	}
      if (this->ReducedHilbertSpaceDescription != 0)
	{
	  this->DiagonalizeInHilbertSubspace(this->ReducedHilbertSpaceDescription, File);
	  cout << "----------------------------------------------------------------" << endl;
	  File.close(); 
	  return 0;
	}
  
      if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
	{
	  RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
	  this->Hamiltonian->GetHamiltonian(HRep);
	  if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
#ifdef __LAPACK__
	      if (this->LapackFlag == true)
		{
		  RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
		  if (this->EvaluateEigenvectors == false)
		    {
		      HRep.LapackDiagonalize(TmpDiag);
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File, TmpDiag[j] - this->EnergyShift);
			}
		    }
		  else
		    {
		      RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.LapackDiagonalize(TmpDiag, Q);
		      if (this->EvaluateEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  for (int j = 0; j < this->NbrEigenvalue; ++j)
			    {
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			      Q[j].WriteVector(TmpVectorName);
			      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
			    }
			  cout << endl;			  
			  delete[] TmpVectorName;
			}
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File, TmpDiag[j] - this->EnergyShift, false);
			  if (this->ComputeEnergyFlag == true)
			    {
			      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
			    }
			  File << endl;
			}
		    }
		}
	      else
		{
#endif
		  RealTriDiagonalSymmetricMatrix TmpTriDiag (this->Hamiltonian->GetHilbertSpaceDimension());
		  if (this->EvaluateEigenvectors == false)
		    {
		      HRep.Householder(TmpTriDiag, 1e-7);
		      TmpTriDiag.Diagonalize();
		      TmpTriDiag.SortMatrixUpOrder();
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift);
			}
		    }
		  else
		    {
		      RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.Householder(TmpTriDiag, 1e-7, Q);
		      TmpTriDiag.Diagonalize(Q);
		      TmpTriDiag.SortMatrixUpOrder(Q);
		      if (this->EvaluateEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  for (int j = 0; j < this->NbrEigenvalue; ++j)
			    {
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			      Q[j].WriteVector(TmpVectorName);
			      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
			    }	      
			  cout << endl;
			  delete[] TmpVectorName;
			}
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift, false);
			  if (this->ComputeEnergyFlag == true)
			    {
			      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
			    }
			  File << endl;
			}
		    }
#ifdef __LAPACK__
		}
#endif
	    }
	  else
	    {
	      this->WriteResult(File, HRep(0, 0)  - this->EnergyShift, false);
	      if (this->ComputeEnergyFlag == true)
		File << " " << (HRep(0, 0)  - this->EnergyShift) ;
	      File << endl;
	    }
	}
      else
	{
	  AbstractLanczosAlgorithm* Lanczos = AlgorithmManager->GetLanczosAlgorithm(this->Architecture, this->EvaluateEigenvectors, this->LapackFlag);
	  if (this->LanczosPrecision != 0.0)
	    Lanczos->SetEigenvaluePrecision(this->LanczosPrecision);
	  double GroundStateEnergy;
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(this->Hamiltonian);
	  if ((this->DiskFlag == true) && (this->ResumeFlag == true))
	    Lanczos->ResumeLanczosAlgorithm();
	  else
	    {
	      if (this->BlockLanczosFlag == false)
		{
		  if ((ExplicitInitialVector == 0) && (this->InitialVectorFileName == 0))
		    Lanczos->InitializeLanczosAlgorithm();
		  else
		    {
			  if (ExplicitInitialVector != 0)
				Lanczos->InitializeLanczosAlgorithm(*ExplicitInitialVector);
		      else
				{
					RealVector InitialVector;
					InitialVector.ReadVector(this->InitialVectorFileName);
					Lanczos->InitializeLanczosAlgorithm(InitialVector);
				}
		    }
		}
	      else
		{
		  if (this->InitialBlockVectorFileName == 0)
		    Lanczos->InitializeLanczosAlgorithm();
		  else
		    {
		      int TmpNbrInitialVectors;
		      ConfigurationParser InitialVectorDescription;
		      if (InitialVectorDescription.Parse(this->InitialBlockVectorFileName) == false)
			{
			  InitialVectorDescription.DumpErrors(cout) << endl;
			}
		      else
			{
			  char** VectorFileNames;
			  if (InitialVectorDescription.GetAsStringArray("InitialVectors", ' ', VectorFileNames, TmpNbrInitialVectors) == false)
			    {
			      cout << "Vectors are not defined or have a wrong value in " << this->InitialBlockVectorFileName << endl;
			    }
			  else
			    {
			      RealVector* InitialVectors = new RealVector[TmpNbrInitialVectors];
			      for (int i = 0; i < TmpNbrInitialVectors; ++i)
				InitialVectors[i].ReadVector(VectorFileNames[i]);
			      Lanczos->InitializeLanczosAlgorithm(InitialVectors, TmpNbrInitialVectors);		  
			      delete[] InitialVectors;
			    }
			}
		    }
		}
	    }
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  timeval TotalCurrentTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  int StartTimeSecond = TotalStartingTime.tv_sec;
	  if (this->ResumeFlag == false)
	    {
	      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	      CurrentNbrIterLanczos = NbrEigenvalue + 3;
	    }
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  gettimeofday (&(TotalCurrentTime), 0); 
	  int CurrentTimeSecond = TotalCurrentTime.tv_sec;
	  while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
											 ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
							   ((this->DiskFlag == false) && (((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
											  ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos))))))
	    {
	      if (this->BlockLanczosFlag == true)
		CurrentNbrIterLanczos += this->SizeBlockLanczos;
	      else
		++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << " ";
	      gettimeofday (&(TotalEndingTime), 0);
	      CurrentTimeSecond = TotalEndingTime.tv_sec;
	      if (this->ShowIterationTime == true)
		{
		  Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		    ((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
		  cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")";
		  TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
		  TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
		}
	      cout << endl;
	      if ((SpectralResponseSaveInterval > 0)&&(CurrentNbrIterLanczos%SpectralResponseSaveInterval == 0))
	        {
		  char* TmpName = new char [strlen(this->EigenvectorFileName) + 16];
                  sprintf (TmpName, "%s.omega_%g-%g_eps_%g.ni_%d.sr", this->EigenvectorFileName,SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, CurrentNbrIterLanczos);
	          ofstream File(TmpName, ios::out);
	          File.precision(14);
         	  Lanczos->SampleSpectralResponse(File, SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, SpectralResponseOmegaInterval, SpectralResponseSpectralResolution);
		  File.close();
		  delete [] TmpName;
		}
	    }
	  if ((Lanczos->TestConvergence() == true) && (CurrentNbrIterLanczos == 0))
	    {
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	    }
	  if (CurrentNbrIterLanczos >= this->MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      return 1;
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i < this->NbrEigenvalue; ++i)
	    {
	      cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
	      if  (this->ComputeEnergyFlag == false)
		{
		  this->WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift);
		}
	    }
	  cout << endl;
	  if ((this->EvaluateEigenvectors == true) && 
	      (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
					     ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
	       ((this->DiskFlag == false) && (((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
					      ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos))))))
	    {
	      RealVector* Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
	      if (Eigenvectors != 0)
		{
		  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		  if ((this->EigenvectorConvergence == true) && ((this->PartialLanczos == false) || (CurrentNbrIterLanczos <= this->NbrIterLanczos)))
		    {
		      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
		      Operation1.ApplyOperation(this->Architecture);
		      double Scalar = TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1];
		      Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));
		      while (Precision > 1e-7)
			{
			  ++CurrentNbrIterLanczos;
			  Lanczos->RunLanczosAlgorithm(1);
			  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
			  TmpMatrix.SortMatrixUpOrder();
			  Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
			  delete[] Eigenvectors;
			  Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
			  VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
			  Operation1.ApplyOperation(this->Architecture);
			  Scalar = TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1];
			  Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));		  
			  cout << (TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift) << " " << (Scalar - this->EnergyShift) << " " 
			       << Precision << " ";
			  if (this->ShowIterationTime == true)
			    {
			      gettimeofday (&(TotalEndingTime), 0);
			      Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
			      cout << "(" << Dt << " s)";
			      TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
			      TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
			    }
			  cout << endl;
			}
		    }
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 32];
		  for (int i = 0; i < this->NbrEigenvalue; ++i)
		    {
		      if (this->ComputeEnergyFlag == true)
			this->WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift);
		      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[i]), &TmpEigenvector);
		      Operation1.ApplyOperation(this->Architecture);
		      cout << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << " ";	
		      if (this->EvaluateEigenvectors == true)
			{
			  if ((this->PartialLanczos == false) || (CurrentNbrIterLanczos < this->NbrIterLanczos))
			    {	  
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, (Lanczos->EigenstateIndexShift() + i));
			    }
			  else
			    {
			      sprintf (TmpVectorName, "%s.%d.part.vec", this->EigenvectorFileName, (Lanczos->EigenstateIndexShift() + i));		  
			    }
			  Eigenvectors[i].WriteVector(TmpVectorName);
			}
		      if (this->ComputeEnergyFlag == true)
			{
			  File << " " << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << endl;
			}
		    }
		  cout << endl;
		  delete[] TmpVectorName;
		  delete[] Eigenvectors;
		}
	      else
		{
		  cout << "eigenvectors can't be computed" << endl;
		}
	    }

	    //print final spectral response	    
	    if (SpectralResponseSaveInterval != 0)
	    {
	      char* TmpName = new char [strlen(this->EigenvectorFileName) + 64];
	      sprintf (TmpName, "%s.omega_%g-%g_eps_%g.sr", this->EigenvectorFileName,SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon);
	      ofstream File(TmpName, ios::out);
	      File.precision(14);
	      Lanczos->SampleSpectralResponse(File, SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, SpectralResponseOmegaInterval, SpectralResponseSpectralResolution);
	      File.close();
	      delete [] TmpName;
	    }
	    
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      File.close();
    }
  else // have both kx and ky momentum -> complex problem
    {
      if (this->SavePrecalculationFileName != 0)
	{
	  this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
	}
      if (this->ReducedHilbertSpaceDescription != 0)
	{
	  this->ComplexDiagonalizeInHilbertSubspace(this->ReducedHilbertSpaceDescription, File);
	  cout << "----------------------------------------------------------------" << endl;
	  File.close(); 
	  return 0;
	}
      if (this->Hamiltonian->GetHilbertSpaceDimension() == 0) 
	return 0;
      if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
	{
	  HermitianMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
	  this->Hamiltonian->GetHamiltonian(HRep);
	  if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
#ifdef __LAPACK__
	      if (this->LapackFlag == true)
		{
		  RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
		  if (this->EvaluateEigenvectors == false)
		    {
		      HRep.LapackDiagonalize(TmpDiag);
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			this->WriteResult(File, TmpDiag[j] - this->EnergyShift);
		    }
		  else
		    {
		      ComplexMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.LapackDiagonalize(TmpDiag, Q);
		      if (this->EvaluateEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  ComplexVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  for (int j = 0; j < this->NbrEigenvalue; ++j)
			    {
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			      Q[j].WriteVector(TmpVectorName);
			      cout << (TmpDiag[j] - this->EnergyShift) << " " <<((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
			    }
			  cout << endl;			  
			  delete[] TmpVectorName;
			}
		      
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File,TmpDiag[j] - this->EnergyShift, false);
			  if (this->ComputeEnergyFlag == true)
			    {
			      ComplexVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
			    }
			  File << endl;
			}
		    }
		}
	      else
		{
#endif
		  RealDiagonalMatrix EVs(this->Hamiltonian->GetHilbertSpaceDimension());
		  if (this->EvaluateEigenvectors == false)
		    {		  
		      HRep.Diagonalize(EVs, /* error */ 1e-10 , /* maxIter */ 250);  
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			this->WriteResult(File, EVs[j] - this->EnergyShift);
		    }
		  else
		    {
		      ComplexMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.Diagonalize(EVs, Q, /* error */ 1e-10 , /* maxIter */ 250);
		  
		      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		      ComplexVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		      for (int j = 0; j < this->NbrEigenvalue; ++j)
			{
			  this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			  sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			  Q[j].WriteVector(TmpVectorName);
			  cout << ((Q[j]*TmpEigenvector) - this->EnergyShift) << " " << endl;
			}
		      cout << endl;
		      delete[] TmpVectorName;

		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  this->WriteResult(File, EVs[j] - this->EnergyShift, false);
			  if (this->ComputeEnergyFlag == true)
			    {
			      ComplexVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      File << " " << ((Q[j]*TmpEigenvector) - this->EnergyShift);
			    }
			  File << endl;
			}		  
		    }
#ifdef __LAPACK__
		}
#endif
	    }
	  else // Hilbert space of dimension one
	    {
	      this->WriteResult(File, HRep(0, 0)  - this->EnergyShift, false);
	      if (this->ComputeEnergyFlag == true)
		File << " " << (HRep(0, 0)  - this->EnergyShift) ;
	      File << endl;
	      if (this->EvaluateEigenvectors)
		{
	      
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  ComplexVector TmpEigenvector(1);
		  TmpEigenvector[0]=1.0;
		  sprintf (TmpVectorName, "%s.0.vec", this->EigenvectorFileName);
		  TmpEigenvector.WriteVector(TmpVectorName);
		  delete [] TmpVectorName;
		}
	    }
	}
      else
	{
	  AbstractLanczosAlgorithm* Lanczos = AlgorithmManager->GetLanczosAlgorithm(this->Architecture, this->EvaluateEigenvectors, this->LapackFlag);
	  if (this->LanczosPrecision != 0.0)
	    Lanczos->SetEigenvaluePrecision(this->LanczosPrecision);
	  double GroundStateEnergy;
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(this->Hamiltonian);
	  if ((this->DiskFlag == true) && (this->ResumeFlag == true))
	    Lanczos->ResumeLanczosAlgorithm();
	  else
	    {
	      if (this->InitialVectorFileName == 0)
		{
		  Lanczos->InitializeLanczosAlgorithm();
		}
	      else
		{	      
		  ComplexVector InitialVector;
		  InitialVector.ReadVector(this->InitialVectorFileName);
		  Lanczos->InitializeLanczosAlgorithm(InitialVector);
		}
	    }
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  timeval TotalCurrentTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  int StartTimeSecond = TotalStartingTime.tv_sec;
	  if (this->ResumeFlag == false)
	    {
	      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	      CurrentNbrIterLanczos = NbrEigenvalue + 3;
	    }
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  gettimeofday (&(TotalCurrentTime), 0); 
	  int CurrentTimeSecond = TotalCurrentTime.tv_sec;
	  while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
											 ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
							   ((this->DiskFlag == false) && (((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
											  ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos))))))
	    {
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << " ";
	      gettimeofday (&(TotalEndingTime), 0);
	      CurrentTimeSecond = TotalEndingTime.tv_sec;
	      if (this->ShowIterationTime == true)
		{
		  Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		    ((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
		  cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")";
		  TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
		  TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
		}
	      cout << endl;
	      if ((SpectralResponseSaveInterval > 0)&&(CurrentNbrIterLanczos%SpectralResponseSaveInterval == 0))
	        {
		  char* TmpName = new char [strlen(this->EigenvectorFileName) + 64];
                  sprintf (TmpName, "%s.omega_%g-%g_eps_%g.ni_%d.sr", this->EigenvectorFileName,SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, CurrentNbrIterLanczos);
	          ofstream File(TmpName, ios::out);
	          File.precision(14);
         	  Lanczos->SampleSpectralResponse(File, SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, SpectralResponseOmegaInterval, SpectralResponseSpectralResolution);
		  File.close();
		  delete [] TmpName;
		}
	    }
	    
	  if ((Lanczos->TestConvergence() == true) && (CurrentNbrIterLanczos == 0))
	    {
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	    }
	    
	  if (CurrentNbrIterLanczos >= this->MaxNbrIterLanczos)
	    {
	      cout << "too many Lanczos iterations" << endl;
	      File << "too many Lanczos iterations" << endl;
	      File.close();
	      return 1;
	    }
	    
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	       
	  for (int i = 0; i < this->NbrEigenvalue; ++i)
	    {
	      cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
	      if  (this->ComputeEnergyFlag == false)
		this->WriteResult(File, (TmpMatrix.DiagonalElement(i) - this->EnergyShift));
	    }
	    
	  cout << endl;
	  if ((this->EvaluateEigenvectors == true) && 
	      (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
					     ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
	       ((this->DiskFlag == false) && (((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
					      ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos))))))
	    {
	      ComplexVector* Eigenvectors = (ComplexVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
	      if (Eigenvectors != 0)
		{
		  ComplexVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		  if ((this->EigenvectorConvergence == true) && ((this->PartialLanczos == false) || (CurrentNbrIterLanczos <= this->NbrIterLanczos)))
		    {
		      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
		      Operation1.ApplyOperation(this->Architecture);
		      double Scalar = Norm(TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1]);
		      Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));
		      while (Precision > 1e-7)
			{
			  ++CurrentNbrIterLanczos;
			  Lanczos->RunLanczosAlgorithm(1);
			  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
			  TmpMatrix.SortMatrixUpOrder();
			  Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
			  delete[] Eigenvectors;
			  Eigenvectors = (ComplexVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
			  VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
			  Operation1.ApplyOperation(this->Architecture);
			  Scalar = Norm(TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1]);
			  Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));		  
			  cout << (TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift) << " " << (Scalar - this->EnergyShift) << " " 
			       << Precision << " ";
			  if (this->ShowIterationTime == true)
			    {
			      gettimeofday (&(TotalEndingTime), 0);
			      Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
			      cout << "(" << Dt << " s)";
			      TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
			      TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
			    }
			  cout << endl;
			}
		    }
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 32];
		  for (int i = 0; i < this->NbrEigenvalue; ++i)
		    {
		      if (this->ComputeEnergyFlag == true)
			File << (TmpMatrix.DiagonalElement(i) - this->EnergyShift);
		      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[i]), &TmpEigenvector);
		      Operation1.ApplyOperation(this->Architecture);
		      cout << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << " ";	
		      if (this->EvaluateEigenvectors == true)
			{
			  if ((this->PartialLanczos == false) || (CurrentNbrIterLanczos < this->NbrIterLanczos))
			    {	  
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, (Lanczos->EigenstateIndexShift() + i));
			    }
			  else
			    {
			      sprintf (TmpVectorName, "%s.%d.part.vec", this->EigenvectorFileName, (Lanczos->EigenstateIndexShift() + i));		  
			    }
			  Eigenvectors[i].WriteVector(TmpVectorName);
			}
		      if (this->ComputeEnergyFlag == true)
			{
			  File << " " << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << endl;
			}
		    }
		  cout << endl;
		  delete[] TmpVectorName;
		  delete[] Eigenvectors;
		}
	      else
		{
		  cout << "eigenvectors can't be computed" << endl;
		}
	    }
	    
	    //print final spectral repsonse	    
	    if (SpectralResponseSaveInterval != 0)
	    {
	      char* TmpName = new char [strlen(this->EigenvectorFileName) + 64];
	      sprintf (TmpName, "%s.omega_%g-%g_eps_%g.sr", this->EigenvectorFileName,SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon);
	      ofstream File(TmpName, ios::out);
	      File.precision(14);
	      Lanczos->SampleSpectralResponse(File, SpectralResponseOmegaMin, SpectralResponseOmegaMax, SpectralResponseEpsilon, SpectralResponseOmegaInterval, SpectralResponseSpectralResolution);
	      File.close();
	      delete [] TmpName;
	    }
	  
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
    }
  return 0;
}
  


// write a line of output to the results file
//
// file = stream to write to
// value = numerical value to be printed after columns for flux and momentum (if defined)
// terminate = indicate if line should be terminated with endl
// return value = stream to write to

ofstream& FQHEOnTorusMainTask::WriteResult(ofstream& file, double value, bool terminate)
{
  if (this->KyOnlyFlag)
    file << this->KyValue << " ";
  else
    file << this->KxValue << " " << this->KyValue << " ";
  file << value;
  if (this->MultiplicityFlag)
    file << " " << this->Multiplicity;
  if (terminate)
    file << endl;
  return file;
}

// do the Hamiltonian diagonalization in a given Hilbert subspace
//
// subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
// file = reference on the output file stream where eigenvalues have to be stored

void FQHEOnTorusMainTask::DiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file)
{
  ConfigurationParser ReducedBasis;
  if (ReducedBasis.Parse(subspaceDescription) == false)
    {
      ReducedBasis.DumpErrors(cout) << endl;
      return;
    }
  int TmpHilbertSpaceDimension;
  char** VectorFileNames;
  if (ReducedBasis.GetAsStringArray("Basis", ' ', VectorFileNames, TmpHilbertSpaceDimension) == false)
    {
      cout << "Vectors are not defined or have a wrong value in " << subspaceDescription << endl;
      return;
    }
  RealMatrix Basis (this->Space->GetHilbertSpaceDimension(), TmpHilbertSpaceDimension);
  char* DirectoryName = ReducedBasis["Directory"];
  char* TmpName;
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpName = VectorFileNames[i];
      if (DirectoryName != 0)
	{
	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	}
      cout << TmpName << endl;
      if (Basis[i].ReadVector(TmpName) == false)
	{
	  cout << "error while reading " << TmpName << endl;
	  if (DirectoryName != 0)
	    delete[] TmpName;
	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	    delete[] VectorFileNames[j];
	  delete[] VectorFileNames;
	  return;
	}
      if (DirectoryName != 0)
	delete[] TmpName;
    }
  RealSymmetricMatrix HRep (TmpHilbertSpaceDimension);
  RealVector* TmpVectors = new RealVector[TmpHilbertSpaceDimension];
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      RealVector TmpVector (Basis[0].GetVectorDimension(), true);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Basis[i]), &TmpVector);
      Operation1.ApplyOperation(this->Architecture);
      TmpVectors[i] = TmpVector;
    }
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      for (int j = i; j < TmpHilbertSpaceDimension; ++j)
	{
	  HRep(i ,j) = Basis[j] * TmpVectors[i];
	}
    }
  delete[] TmpVectors;
  if (this->ShowHamiltonian == true)
    cout << HRep << endl;
  bool DiagonalOnlyFlag = false;
  ReducedBasis.GetAsBoolean("DiagonalOnly", DiagonalOnlyFlag);
  if (DiagonalOnlyFlag  == true)
    {
      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	{
	  double Tmp = 0.0;
	  HRep.GetMatrixElement(j, j, Tmp);
	  this->WriteResult(file, (Tmp - this->EnergyShift), true);
	}      
    }
  else
    {
      if (TmpHilbertSpaceDimension > 1)
	{
#ifdef __LAPACK__
	  if (this->LapackFlag == true)
	    {
	      RealDiagonalMatrix TmpDiag (TmpHilbertSpaceDimension);
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.LapackDiagonalize(TmpDiag);
		}
	      else
		{
		  RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
		  for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		    TmpEigenvector(l, l) = 1.0;
		  HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
		  if (this->ReducedHilbertSpaceExportTransformation != 0)
		    {
		      TmpEigenvector.WriteMatrix(this->ReducedHilbertSpaceExportTransformation);
		    }
		  Basis.Multiply(TmpEigenvector);
		}
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  this->WriteResult(file, (TmpDiag[j] - this->EnergyShift), true);
		}
	    }
	  else
	    {
#endif
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (TmpHilbertSpaceDimension);
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.Householder(TmpTriDiag, 1e-7);
		  TmpTriDiag.Diagonalize();
		  TmpTriDiag.SortMatrixUpOrder();
		}
	      else
		{
		  RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
		  for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		    TmpEigenvector(l, l) = 1.0;
		  HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
		  TmpTriDiag.Diagonalize(TmpEigenvector);
		  TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
		  if (this->ReducedHilbertSpaceExportTransformation != 0)
		    {
		      TmpEigenvector.WriteMatrix(this->ReducedHilbertSpaceExportTransformation);
		    }
		  Basis.Multiply(TmpEigenvector);
		}
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  this->WriteResult(file, (TmpTriDiag.DiagonalElement(j) - this->EnergyShift), true);
		}
#ifdef __LAPACK__
	    }
#endif
	  if (this->EvaluateEigenvectors == true)
	    {
	      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
		  Basis[j].WriteVector(TmpVectorName);
		}
	      delete [] TmpVectorName;
	    }
	}
      else
	{
	  this->WriteResult(file, (HRep(0, 0) - this->EnergyShift), true);
	}
    }
  for (int j= 0; j < TmpHilbertSpaceDimension; ++j)
    delete[] VectorFileNames[j];
  delete[] VectorFileNames;
}

// do the Hamiltonian diagonalization in a given Hilbert subspace, when the hamiltonian is complex
//
// subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
// file = reference on the output file stream where eigenvalues have to be stored

void FQHEOnTorusMainTask::ComplexDiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file)
{
  ConfigurationParser ReducedBasis;
  if (ReducedBasis.Parse(subspaceDescription) == false)
    {
      ReducedBasis.DumpErrors(cout) << endl;
      return;
    }
  int TmpHilbertSpaceDimension;
  char** VectorFileNames;
  if (ReducedBasis.GetAsStringArray("Basis", ' ', VectorFileNames, TmpHilbertSpaceDimension) == false)
    {
      cout << "Vectors are not defined or have a wrong value in " << subspaceDescription << endl;
      return;
    }
  ComplexMatrix Basis (this->Space->GetHilbertSpaceDimension(), TmpHilbertSpaceDimension);
  char* DirectoryName = ReducedBasis["Directory"];
  char* TmpName;
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpName = VectorFileNames[i];
      if (DirectoryName != 0)
	{
	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	}
      cout << TmpName << endl;
      if (Basis[i].ReadVector(TmpName) == false)
	{
	  cout << "error while reading " << TmpName << endl;
	  if (DirectoryName != 0)
	    delete[] TmpName;
	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	    delete[] VectorFileNames[j];
	  delete[] VectorFileNames;
	  return;
	}
      if (DirectoryName != 0)
	delete[] TmpName;
    }
  HermitianMatrix HRep (TmpHilbertSpaceDimension);
  ComplexVector* TmpVectors = new ComplexVector[TmpHilbertSpaceDimension];
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      ComplexVector TmpVector (Basis[0].GetVectorDimension(), true);
      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Basis[i]), &TmpVector);
      Operation1.ApplyOperation(this->Architecture);
      TmpVectors[i] = TmpVector;
    }
  for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
    {
      for (int j = i; j < TmpHilbertSpaceDimension; ++j)
	{
	  Complex Tmp = Basis[j] * TmpVectors[i];
	  HRep.SetMatrixElement(i ,j, Tmp);
	}
    }
  delete[] TmpVectors;
  if (this->ShowHamiltonian == true)
    cout << HRep << endl;
  bool DiagonalOnlyFlag = false;
  ReducedBasis.GetAsBoolean("DiagonalOnly", DiagonalOnlyFlag);
  if (DiagonalOnlyFlag  == true)
    {
      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	{
	  double Tmp = 0.0;
	  HRep.GetMatrixElement(j, j, Tmp);
	  this->WriteResult(file, (Tmp - this->EnergyShift), true);
	}      
    }
  else
    {
      if (TmpHilbertSpaceDimension > 1)
	{
	  RealDiagonalMatrix TmpDiag (TmpHilbertSpaceDimension);
#ifdef __LAPACK__
	  if (this->LapackFlag == true)
	    {
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.LapackDiagonalize(TmpDiag);
		}
	      else
		{
		  ComplexMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
		  for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		    TmpEigenvector(l, l) = 1.0;
		  HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
	      if (this->ReducedHilbertSpaceExportTransformation != 0)
		{
		  TmpEigenvector.WriteMatrix(this->ReducedHilbertSpaceExportTransformation);
		}
	      Basis.Multiply(TmpEigenvector);
		}
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  this->WriteResult(file, (TmpDiag[j] - this->EnergyShift), true);
		}
	    }
	  else
	    {
#endif
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.Diagonalize(TmpDiag);
		}
	      else
		{
		  ComplexMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
		  for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		    TmpEigenvector(l, l) = 1.0;
		  HRep.Diagonalize(TmpDiag, TmpEigenvector);
		  if (this->ReducedHilbertSpaceExportTransformation != 0)
		    {
		      TmpEigenvector.WriteMatrix(this->ReducedHilbertSpaceExportTransformation);
		    }
		  Basis.Multiply(TmpEigenvector);
		}
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  this->WriteResult(file, (TmpDiag[j] - this->EnergyShift), true);
		}
#ifdef __LAPACK__
	    }
#endif
	  if (this->EvaluateEigenvectors == true)
	    {
	      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
	      for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
		{
		  sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
		  Basis[j].WriteVector(TmpVectorName);
		}
	      delete [] TmpVectorName;
	}
	}
      else
	{
	  this->WriteResult(file, (HRep(0, 0) - this->EnergyShift), true);
	}
    }
  for (int j= 0; j < TmpHilbertSpaceDimension; ++j)
    delete[] VectorFileNames[j];
  delete[] VectorFileNames;
}

// set a kx-value
//
// kxValue = kx value

void FQHEOnTorusMainTask::SetKxValue(int kxValue)
{
  this->KyOnlyFlag = false;
  this->RealFlag = false;
  this->KxValue = kxValue;
}

// set multiplicity of a given momentum sector
//
// multiplicity = sector multiplicity

void FQHEOnTorusMainTask::SetMultiplicity(int multiplicity) 
{
  this->MultiplicityFlag = true; 
  this->Multiplicity = multiplicity;
}

