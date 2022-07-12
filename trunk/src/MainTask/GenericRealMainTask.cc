////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                       class author: Gunnar Moeller                         //
//                                                                            //
//                      class of a generic real main task                     //
//                                                                            //
//                       last modification : 03/08/2009                       //
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

#include "MainTask/GenericRealMainTask.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/HamiltonianFullDiagonalizeOperation.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"

#include "Options/Options.h"

#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/GenericSignalHandler.h"

#include "Vector/ComplexVector.h"


#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



// constructor
//  
// options = pointer to the options managers containing all running options
// space = pointer to the current Hilbert space
// hamiltonian = pointer to the current Hamiltonian
// subspaceStr = string to prepend in output file for each eigenvalue in this subspace
// subspaceLegend = legend indicating contents of subspaceStr to include in output file
// shift = energy shift that is applied to the hamiltonian
// outputFileName = name of the file where results have to be stored
// firstRun = flag that indicates if it the first time the main task is used
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector

GenericRealMainTask::GenericRealMainTask(OptionManager* options, AbstractHilbertSpace* space, LanczosManager* lanczos, 
					 AbstractHamiltonian* hamiltonian, const char *subspaceStr, const char *subspaceLegend,
					 double shift, char* outputFileName, bool firstRun, char* eigenvectorFileName, bool fakeComplex)
{
  if (outputFileName != 0)
    {
      this->OutputFileName = new char [strlen(outputFileName) + 1];
      strcpy(this->OutputFileName, outputFileName);
      this->OutputFileName[strlen(outputFileName)] = '\0';
    }
  else
    {
      this->OutputFileName = 0;
    }
  if (eigenvectorFileName == 0)
    {
      this->EigenvectorFileName = 0;
    }
  else
    {
      this->EigenvectorFileName = new char [strlen(eigenvectorFileName) + 1];
      strcpy(this->EigenvectorFileName, eigenvectorFileName);
      this->EigenvectorFileName[strlen(eigenvectorFileName)] = '\0';
    }
  this->Hamiltonian = hamiltonian;
  this->Space = space;
  this->AlgorithmManager = lanczos;
  this->EnergyShift = shift;
  this->SubspaceStr = new char[strlen(subspaceStr)+1];
  strcpy(this->SubspaceStr,subspaceStr);
  this->SubspaceLegend = new char[strlen(subspaceLegend)+1];
  strcpy(this->SubspaceLegend,subspaceLegend);
  this->ResumeFlag = options->GetBoolean("resume");
  this->DiskFlag = options->GetBoolean("disk");
  this->MaxNbrIterLanczos = options->GetInteger("iter-max");
  this->NbrIterLanczos = options->GetInteger("nbr-iter");
  this->NbrEigenvalue = options->GetInteger("nbr-eigen");
  if (this->NbrEigenvalue > this->Space->GetHilbertSpaceDimension())
    {
      this->NbrEigenvalue = this->Space->GetHilbertSpaceDimension();
    }
  this->FullDiagonalizationLimit = options->GetInteger("full-diag");
  this->BlockLanczosFlag = false;
  if ((*options)["block-lanczos"] != 0)
    {
      this->BlockLanczosFlag = options->GetBoolean("block-lanczos");
    }
  this->SizeBlockLanczos = 1;
  if ((*options)["block-size"] != 0)
    {
      this->SizeBlockLanczos = options->GetInteger("block-size");
    }
  this->VectorMemory = options->GetInteger("nbr-vector");
  if ((*options)["save-precalculation"] != 0)
    {
      this->SavePrecalculationFileName = options->GetString("save-precalculation");
    }
  else
    this->SavePrecalculationFileName = 0;
  this->FullReorthogonalizationFlag = options->GetBoolean("force-reorthogonalize");
  this->EvaluateEigenvectors = options->GetBoolean("eigenstate");
  if ((*options)["all-eigenstates"] != 0)
    {
      this->EvaluateAllEigenvectors = options->GetBoolean("all-eigenstates");
    }
  else
    {
      this->EvaluateAllEigenvectors = false;
    }
  if ((*options)["first-eigenstate"] != 0)
    {
      this->FirstEigenstateIndex = options->GetInteger("first-eigenstate");
    }
  else
    {
      this->FirstEigenstateIndex = 0;
    }
  this->EigenvectorConvergence = options->GetBoolean("eigenstate-convergence");
  if ((*options)["show-itertime"] != 0)
    {
      this->ShowIterationTime = options->GetBoolean("show-itertime");
    }
  else
    this->ShowIterationTime = false;
  if ((*options)["initial-vector"] != 0)
    {
      this->InitialVectorFileName = options->GetString("initial-vector");
    }
  else
    {
      this->InitialVectorFileName = 0;
    }
  if ((*options)["initial-blockvectors"] != 0)
    {
      this->InitialBlockVectorFileName = options->GetString("initial-blockvectors");
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
  if ((*options)["use-lapack"] != 0)
    {
      this->LapackFlag = options->GetBoolean("use-lapack");
    }
  else
    {
      this->LapackFlag = false;
    }
  if ((*options)["use-scalapack"] != 0)
    {
      this->ScalapackFlag = options->GetBoolean("use-scalapack");
    }
  else
    {
      this->ScalapackFlag = false;
    }
  if ((*options)["limit-time"] != 0)
    {
      this->MaximumAllowedTime = (options->GetInteger("limit-time"));
    }
  else
    {
      this->MaximumAllowedTime = 0;
    }
  if ((((*options)["use-hilbert"]) != 0) && (options->GetString("use-hilbert") != 0))
    {
      this->ReducedHilbertSpaceDescription = options->GetString("use-hilbert");
    }
  else
    {
      this->ReducedHilbertSpaceDescription = 0;
    }
  if ((*options)["get-hvalue"] != 0)
    {
      this->ComputeEnergyFlag = options->GetBoolean("get-hvalue");
    }
  else
    {
      this->ComputeEnergyFlag = false;
    }
  this->ShowHamiltonian = false;
  if (((*options)["show-hamiltonian"] != 0) && (options->GetBoolean("show-hamiltonian") == true))
    {
      this->ShowHamiltonian = true;
      if (this->ReducedHilbertSpaceDescription == 0)
	{
	  RealMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
	  this->Hamiltonian->GetHamiltonian(HRep);
	  cout << HRep << endl;
	}
    }  
  this->FriendlyShowHamiltonian = false;
  if (((*options)["friendlyshow-hamiltonian"] != 0) && (options->GetBoolean("friendlyshow-hamiltonian") == true))
    {
      this->FriendlyShowHamiltonian = true;
      this->FriendlyShowHamiltonianError = 0.0;
      if ((*options)["friendlyshowhamiltonian-error"] != 0)
	{ 
	  this->FriendlyShowHamiltonianError  = options->GetDouble("friendlyshowhamiltonian-error");
	}
      if (this->ReducedHilbertSpaceDescription == 0)
	{
	  RealMatrix HRep (this->Space->GetHilbertSpaceDimension(), this->Space->GetHilbertSpaceDimension());
	  this->Hamiltonian->GetHamiltonian(HRep);
	  for (int i = 0; i < this->Space->GetHilbertSpaceDimension(); ++i)
	    {
	      if (HRep[i].Norm() > this->FriendlyShowHamiltonianError)
		{
		  cout << i << " : ";
		  this->Space->PrintState(cout, i) << endl;
		  for (int j = 0; j < this->Space->GetHilbertSpaceDimension(); ++j)
		    {
		      if (fabs(HRep[i][j]) > this->FriendlyShowHamiltonianError)
			{
			  cout << "    " << j << " : ";
			  this->Space->PrintState(cout, j) << " : " << HRep[i][j] << endl;
			}
		    }
		}
	    }
	}
    }  
  this->ExportBinaryHamiltonian = 0;
  if (((*options)["export-binhamiltonian"] != 0) && (options->GetString("export-binhamiltonian") != 0)) 
    {
      if (this->ReducedHilbertSpaceDescription == 0)
	{
	  RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
	  this->Hamiltonian->GetHamiltonian(HRep);
	  HRep.WriteMatrix(options->GetString("export-binhamiltonian"));	  
	}
      else
	{
	  this->ExportBinaryHamiltonian = new char [strlen(options->GetString("export-binhamiltonian")) + 1];
	  strcpy (this->ExportBinaryHamiltonian, options->GetString("export-binhamiltonian"));
	}
    }
  if (((*options)["export-hamiltonian"] != 0) && (options->GetString("export-hamiltonian") != 0))
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      HRep.SparseWriteAsciiMatrix(options->GetString("export-hamiltonian"));
    }  
  if (((*options)["test-hermitian"] != 0) && (options->GetBoolean("test-hermitian") == true))
    {
      RealMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      double Tmp1;
      double Tmp2;
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
	    AverageNorm += fabs(Tmp1);
	  }
      AverageNorm /= 0.5 * ((double) this->Hamiltonian->GetHilbertSpaceDimension()) * ((double) (this->Hamiltonian->GetHilbertSpaceDimension() + 1));
      for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
	for (int j = i; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
	  {
	    HRep.GetMatrixElement(i, j, Tmp1);
	    HRep.GetMatrixElement(j, i, Tmp2);
	    if (fabs(Tmp1 - Tmp2) > (Error * AverageNorm))
	      {
		cout << "error at " << i << " " << j << " : " << Tmp1 << " " << Tmp2 << endl;
	      }
	  }
      cout << "check done" << endl;
    }
  if (((*options)["compute-sparsity"] != 0) && (options->GetBoolean("compute-sparsity") == true))
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
      this->Hamiltonian->GetHamiltonian(HRep);
      long Tmp = HRep.ComputeNbrNonZeroMatrixElements();
      cout << "nbr of non zero matrix elements = " << Tmp << " (" 
	   << (100.0 * ((double) Tmp) / 
	       (((double) this->Hamiltonian->GetHilbertSpaceDimension()) * ((double) this->Hamiltonian->GetHilbertSpaceDimension()))) << "%)" << endl;
    }
  if (((*options)["lanczos-precision"] != 0) && (options->GetDouble("lanczos-precision") > 0))
    {
      this->LanczosPrecision = options->GetDouble("lanczos-precision");
    }
  else
    {
      this->LanczosPrecision = 0.0;
    }
  if (((*options)["fast-disk"] != 0) && (this->NbrEigenvalue == 1) && (this->EvaluateEigenvectors == true))
    {
      this->FastDiskFlag = options->GetBoolean("fast-disk");
      if ((*options)["resume-fastdisk"] != 0)
	{
	  this->ResumeFastDiskFlag = options->GetBoolean("resume-fastdisk");
	}
    }
  else
    {
      this->FastDiskFlag = false;
      this->ResumeFastDiskFlag = false;
    }
  this->PartialEigenstateFlag = 0;
  if (((*options)["partial-eigenstate"] != 0) && (this->EvaluateEigenvectors == true))
    {
      this->PartialEigenstateFlag = options->GetInteger("partial-eigenstate");
    }
  this->FirstRun = firstRun;
  this->FakeComplex=fakeComplex;
}
 
// destructor
//  

GenericRealMainTask::~GenericRealMainTask()
{
  if (this->OutputFileName != 0)
    delete[] this->OutputFileName;
  if (this->EigenvectorFileName != 0)
    delete[] this->EigenvectorFileName;
  delete [] this->SubspaceStr;
  delete [] this->SubspaceLegend;
}
  
// set architecture bound to the task
// 
// architecture = pointer to the architecture to use

void GenericRealMainTask::SetArchitecture(AbstractArchitecture* architecture)
{
  this->Architecture = architecture;
  if ((this->Architecture->GetArchitectureID() & AbstractArchitecture::WithCommunicator) != 0)
    if (this->OperationManagers.GetNbrElement() == 0)
      {
	this->OperationManagers += new ArchitectureBaseOperationManager((SimpleMPIArchitecture*) this->Architecture, this->Hamiltonian);
      }
}

// execute the main task
// 
// return value = 0 if no error occurs, else return error code

int GenericRealMainTask::ExecuteMainTask()
{
  ofstream File;
  if (this->OutputFileName != 0)
    {
      if (this->FirstRun == true)
	{
	  File.open(this->OutputFileName, ios::binary | ios::out);
	  this->FirstRun = false;
	  File << "# "<<SubspaceLegend;
	  File <<" E";
	  if ((this->EvaluateEigenvectors == true) && (this->ComputeEnergyFlag == true))
	    File << " <H>";
	  File << endl;
	}
      else
	{
	  File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
	}
      File.precision(14);
    }
  cout.precision(14);
  cout << "----------------------------------------------------------------" << endl;
  cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
  if (this->ReducedHilbertSpaceDescription != 0)
    {
      this->DiagonalizeInHilbertSubspace(this->ReducedHilbertSpaceDescription, File);
      cout << "----------------------------------------------------------------" << endl;
      if (this->OutputFileName != 0)
	{
	  File.close(); 
	}
      return 0;
    }
  if (this->SavePrecalculationFileName != 0)
    {
      this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
    }
  if (this->Hamiltonian->GetHilbertSpaceDimension() == 0)
    return 0;
  if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
    {
#ifdef __SCALAPACK__      
      if (this->ScalapackFlag == true)
	{
	  if ((this->Architecture->GetArchitectureID() & AbstractArchitecture::WithCommunicator) == 0)
	    {
	      cout << "error : SCALAPACK requires a MPI enable architecture" << endl;
	      return 1;
	    }	  
	  int TmpNbrEigenstates = this->NbrEigenvalue;
	  bool TmpEvaluateEigenvectors = this->EvaluateEigenvectors | this->EvaluateAllEigenvectors;
	  if (this->ComputeEnergyFlag == true)
	    {
	      TmpNbrEigenstates = 0;
	      TmpEvaluateEigenvectors = true;
	    }
	  if (this->EvaluateAllEigenvectors == true)
	    {
	      TmpNbrEigenstates = this->Hamiltonian->GetHilbertSpaceDimension();	      
	    }
	  HamiltonianFullDiagonalizeOperation Operation1 (this->Hamiltonian, false, TmpEvaluateEigenvectors, TmpNbrEigenstates);
	  Operation1.ApplyOperation(this->Architecture);
	  this->EigenvalueMatrix = Operation1.GetDiagonalizedHamiltonian();
	  if ((this->EvaluateEigenvectors == false) && (this->EvaluateAllEigenvectors == false))
	    {
	      if (this->OutputFileName != 0)
		{
		  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		    this->WriteResult(File, this->EigenvalueMatrix[j] - this->EnergyShift);
		}
	    }
	  else
	    {
	      RealMatrix Eigenstates;
	      Operation1.GetHamiltonianEigenstates(Eigenstates);
	      if (this->EvaluateAllEigenvectors == true)
		{
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  sprintf (TmpVectorName, "%s.eigenvec.mat", this->EigenvectorFileName);
		  Eigenstates.WriteMatrix(TmpVectorName);
		  delete[] TmpVectorName;
		}
	      else
		{
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		  int LastEigenstateIndex = this->FirstEigenstateIndex + this->NbrEigenvalue;
		  if (LastEigenstateIndex > this->Hamiltonian->GetHilbertSpaceDimension())
		    LastEigenstateIndex = this->Hamiltonian->GetHilbertSpaceDimension();
		  for (int j = this->FirstEigenstateIndex; j < LastEigenstateIndex; ++j)
		    {
		      this->Hamiltonian->LowLevelMultiply(Eigenstates[j], TmpEigenvector);
		      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
		      Eigenstates[j].WriteVector(TmpVectorName);
		      cout << ((TmpEigenvector * Eigenstates[j]) - this->EnergyShift) << " " << endl;		  
		    }
		  cout << endl;			  
		  delete[] TmpVectorName;
		}
	      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		{
		  if (this->OutputFileName != 0)
		    {
		      this->WriteResult(File,this->EigenvalueMatrix[j] - this->EnergyShift, false);
		    }
		  if (this->ComputeEnergyFlag == true)
		    {
		      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		      this->Hamiltonian->LowLevelMultiply(Eigenstates[j], TmpEigenvector);
		      if (this->OutputFileName != 0)
			{
			  File << " " << ((TmpEigenvector * Eigenstates[j]) - this->EnergyShift);
			}
		    }
		  if (this->OutputFileName != 0)
		    {
		      File << endl;
		    }
		}
	    }
	}
      else
	{
#endif
	  RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
	  this->Hamiltonian->GetHamiltonian(HRep);
	  if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
#ifdef __LAPACK__
	      if (this->LapackFlag == true)
		{
		  this->EigenvalueMatrix = RealDiagonalMatrix (this->Hamiltonian->GetHilbertSpaceDimension());
		  if ((this->EvaluateEigenvectors == false) && (this->EvaluateAllEigenvectors == false))
		    {
		      HRep.LapackDiagonalize(this->EigenvalueMatrix);
		      if (this->OutputFileName != 0)
			{
			  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			    this->WriteResult(File, this->EigenvalueMatrix[j] - this->EnergyShift);
			}
		    }
		  else
		    {
		      RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.LapackDiagonalize(this->EigenvalueMatrix, Q);
		      if (this->EvaluateAllEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  sprintf (TmpVectorName, "%s.eigenvec.mat", this->EigenvectorFileName);
			  Q.WriteMatrix(TmpVectorName);
			  delete[] TmpVectorName;
			}
		      if (this->EvaluateEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  int LastEigenstateIndex = this->FirstEigenstateIndex + this->NbrEigenvalue;
			  if (LastEigenstateIndex > this->Hamiltonian->GetHilbertSpaceDimension())
			    LastEigenstateIndex = this->Hamiltonian->GetHilbertSpaceDimension();
			  for (int j = this->FirstEigenstateIndex; j < LastEigenstateIndex; ++j)
			    {
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			      if (FakeComplex)
				{
				  ComplexVector TmpVector(Q[j],true);
				  TmpVector.WriteVector(TmpVectorName);
				}
			      else
				Q[j].WriteVector(TmpVectorName);
			      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
			    }
			  cout << endl;			  
			  delete[] TmpVectorName;
			}
		      
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  if (this->OutputFileName != 0)
			    {
			      this->WriteResult(File,this->EigenvalueMatrix[j] - this->EnergyShift, false);
			    }
			  if (this->ComputeEnergyFlag == true)
			    {
			      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      if (this->OutputFileName != 0)
				{
				  File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
				}
			    }
			  if (this->OutputFileName != 0)
			    {
			      File << endl;
			    }
			}
		    }
		}
	      else
		{
#endif
		  RealTriDiagonalSymmetricMatrix TmpTriDiag (this->Hamiltonian->GetHilbertSpaceDimension());
		  if ((this->EvaluateEigenvectors == false) && (this->EvaluateAllEigenvectors == false))
		    {
		      HRep.Householder(TmpTriDiag, 1e-7);
		      TmpTriDiag.Diagonalize();
		      TmpTriDiag.SortMatrixUpOrder();
		      if (this->OutputFileName != 0)
			{
			  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			    this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift);
			}
		    }
		  else
		    {
		      RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		      HRep.Householder(TmpTriDiag, 1e-7, Q);
		      TmpTriDiag.Diagonalize(Q);
		      TmpTriDiag.SortMatrixUpOrder(Q);
		      if (this->EvaluateAllEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  sprintf (TmpVectorName, "%s.eigenvec.mat", this->EigenvectorFileName);
			  Q.WriteMatrix(TmpVectorName);
			  delete[] TmpVectorName;
			}
		      if (this->EvaluateEigenvectors == true)
			{
			  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  int LastEigenstateIndex = this->FirstEigenstateIndex + this->NbrEigenvalue;
			  if (LastEigenstateIndex > this->Hamiltonian->GetHilbertSpaceDimension())
			    LastEigenstateIndex = this->Hamiltonian->GetHilbertSpaceDimension();
			  for (int j = this->FirstEigenstateIndex; j < LastEigenstateIndex; ++j)
			    {
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
			      if (FakeComplex)
				{
				  ComplexVector TmpVector(Q[j],true);
				  TmpVector.WriteVector(TmpVectorName);
				}
			      else
				Q[j].WriteVector(TmpVectorName);
			      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
			    }	      
			  cout << endl;
			  delete[] TmpVectorName;
			}
		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
			{
			  if (this->OutputFileName != 0)
			    {
			      this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift, false);
			    }
			  if (this->ComputeEnergyFlag == true)
			    {
			      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			      if (this->OutputFileName != 0)
				{
				  File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
				}
			    }
			  if (this->OutputFileName != 0)
			    {
			      File << endl;
			    }
			}
		    }
#ifdef __LAPACK__
		}
#endif
	    }
	  else
	    {
	      if (this->OutputFileName != 0)
		{
		  this->WriteResult(File, HRep(0, 0)  - this->EnergyShift, false);
		  if (this->ComputeEnergyFlag == true)
		    File << " " << (HRep(0, 0)  - this->EnergyShift) ;
		  File << endl;	
		}      
	      if (this->EvaluateAllEigenvectors == true)
		{
		  RealMatrix Q(1, 1);
		  Q.SetMatrixElement(0, 0, 1.0);
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  sprintf (TmpVectorName, "%s.eigenvec.mat", this->EigenvectorFileName);
		  Q.WriteMatrix(TmpVectorName);
		  delete[] TmpVectorName;
		}
	      if (this->EvaluateEigenvectors)
		{
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  RealVector TmpEigenvector(1);
		  TmpEigenvector[0] = 1.0;
		  sprintf (TmpVectorName, "%s.0.vec", this->EigenvectorFileName);
		  if (FakeComplex)
		    {
		      ComplexVector TmpVector(TmpEigenvector,true);
		      TmpVector.WriteVector(TmpVectorName);
		    }
		  else
		    TmpEigenvector.WriteVector(TmpVectorName);
		  delete [] TmpVectorName;
		}
	    }
#ifdef __SCALAPACK__      
	}
#endif
    }
  else
    {
      AbstractLanczosAlgorithm* Lanczos = AlgorithmManager->GetLanczosAlgorithm(this->Architecture, this->EvaluateEigenvectors);
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
	      if (this->InitialVectorFileName == 0)
		Lanczos->InitializeLanczosAlgorithm();
	      else
		{
		  if (FakeComplex)
		    {
		      ComplexVector TmpVector;
		      TmpVector.ReadVector(this->InitialVectorFileName);
		      long tmpL;
		      TmpVector.SetStandardPhase(tmpL);
		      RealVector InitialVector(TmpVector);
		      InitialVector.Normalize();
		      Lanczos->InitializeLanczosAlgorithm(InitialVector);
		    }
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
			  ComplexVector TmpVector;
			  long tmpL;
			  for (int i = 0; i < TmpNbrInitialVectors; ++i)
			    {
			      TmpVector.ReadVector(VectorFileNames[i]);
			      TmpVector.SetStandardPhase(tmpL);
			      InitialVectors[i]=RealVector(TmpVector);
			      InitialVectors[i].Normalize();
			    }
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
      GenericSignalHandler Usr1Handler(SIGUSR1);
      while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
										     ((this->MaximumAllowedTime > 0) 
										      && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
						       ((this->DiskFlag == false) && (((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
										      ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos))))))
	{
	  if (this->BlockLanczosFlag == true)
	    CurrentNbrIterLanczos += this->SizeBlockLanczos;
	  else
	    ++CurrentNbrIterLanczos;
	  Usr1Handler.StartToDeferSignal();
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

	  if (Usr1Handler.HavePendingSignal())
	    {
	      cout << "Terminating Lanczos iteration on user signal"<<endl;
	      if (this->OutputFileName != 0)
		{
		  File << "# Lanczos terminated at step "<<CurrentNbrIterLanczos<<" with precision "<<Precision<<endl;
		}
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      for (int i = 0; i < this->NbrEigenvalue; ++i)
		{
		  cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
		  if (this->OutputFileName != 0)
		    {
		      this->WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift, true);
		    }
		}
	      cout << endl;
	    }
	  
	  if ( ((Usr1Handler.HavePendingSignal()) && (this->EvaluateEigenvectors == true))
	       || ((this->PartialEigenstateFlag > 0) && ((CurrentNbrIterLanczos % (this->PartialEigenstateFlag * this->SizeBlockLanczos)) == 0)))
	    {
	      RealVector* Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
	      if (Eigenvectors != 0)
		{
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 32];
		  for (int i = 0; i < this->NbrEigenvalue; ++i)
		    {
		      sprintf (TmpVectorName, "%s.%d.part.%d.vec", this->EigenvectorFileName, (Lanczos->EigenstateIndexShift() + i), CurrentNbrIterLanczos);
		      if (FakeComplex)
			{
			  ComplexVector TmpVector(Eigenvectors[i],true);
			  TmpVector.WriteVector(TmpVectorName);
			}
		      else
			Eigenvectors[i].WriteVector(TmpVectorName);
		    }
		  delete[] TmpVectorName;
		  delete[] Eigenvectors;
		}
	      else
		{
		  cout << "eigenvectors can't be computed" << endl;
		}
	    }
	  Usr1Handler.ProcessDeferredSignal();
	}
      if ((Lanczos->TestConvergence() == true) && (CurrentNbrIterLanczos == 0))
	{
	  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	}
      if (CurrentNbrIterLanczos >= this->MaxNbrIterLanczos)
	{
	  cout << "too many Lanczos iterations" << endl;
	  if (this->OutputFileName != 0)
	    {
	      File << "too many Lanczos iterations" << endl;
	      File.close();
	    }
	  return 1;
	}
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      this->EigenvalueMatrix = RealDiagonalMatrix(this->NbrEigenvalue);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	{
	  cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
	  this->EigenvalueMatrix[i] = TmpMatrix.DiagonalElement(i);
	  if  (this->ComputeEnergyFlag == false)
	    {
	      if (this->OutputFileName != 0)
		{
		  WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift);
		}
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
		    {
		      if (this->OutputFileName != 0)
			{
			  WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift, false);
			}
		    }
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
		      if (FakeComplex)
			{
			  ComplexVector TmpVector(Eigenvectors[i],true);
			  TmpVector.WriteVector(TmpVectorName);
			}
		      else
			Eigenvectors[i].WriteVector(TmpVectorName);
		    }
		  if (this->OutputFileName != 0)
		    {
		      if (this->ComputeEnergyFlag == true)
			{
			  File << " " << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift);
			}
		      if (this->ComputeEnergyFlag == true)
			File << endl;
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
      gettimeofday (&(TotalEndingTime), 0);
      cout << "------------------------------------------------------------------" << endl << endl;;
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
      AlgorithmManager->FreeLanczosAlgorithm();
    }
  cout << "----------------------------------------------------------------" << endl;
  File.close(); 
  return 0;
}

// add optiongroup with options related to this module to the given OptionManager
//
void GenericRealMainTask::AddOptionGroup(OptionManager *optionManager)
{
  this->AlgorithmManager->AddOptionGroup(optionManager);
}


// write a line of output to the results file
//
// file = stream to write to
// value = numerical value to be printed after columns for flux and momentum (if defined)
// terminate = indicate if line should be terminated with endl
void GenericRealMainTask::WriteResult(ofstream& file, double value, bool terminate)
{
  if (SubspaceStr[0] != '\0')
    file << this->SubspaceStr <<" ";
  file << value;
  if (terminate)
    file << endl;
}

// do the Hamiltonian diagonalization in a given Hilbert subspace
//
// subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
// file = reference on the output file stream where eigenvalues have to be stored

void GenericRealMainTask::DiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file)
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
  int SpaceDimension = this->Space->GetHilbertSpaceDimension();
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
      if (Basis[i].GetVectorDimension() != SpaceDimension)
	{
	  cout << "Error: basis vector '" << TmpName << "' does not match the dimension of the Hilbert space."<<endl;
	  exit(1);
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
	  double Tmp = Basis[j] * TmpVectors[i];
	  HRep.SetMatrixElement(i ,j, Tmp);
	}
    }
  if (this->ExportBinaryHamiltonian != 0)
    {
      HRep.WriteMatrix(this->ExportBinaryHamiltonian);	        
    }
  if (this->ShowHamiltonian == true)
    {
      cout << HRep << endl;
    }
  delete[] TmpVectors;
  if (TmpHilbertSpaceDimension > 1)
    {
#ifdef __LAPACK__
      if (this->LapackFlag == true)
	{
	  this->EigenvalueMatrix = RealDiagonalMatrix(TmpHilbertSpaceDimension);
	  if (this->EvaluateEigenvectors == false)
	    {
	      HRep.LapackDiagonalize(this->EigenvalueMatrix);
	    }
	  else
	    {
	      RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
	      for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		TmpEigenvector(l, l) = 1.0;
	      HRep.LapackDiagonalize(this->EigenvalueMatrix, TmpEigenvector);
	      Basis.Multiply(TmpEigenvector);
	    }
	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	    {
	      WriteResult(file, (this->EigenvalueMatrix[j] - this->EnergyShift));
	    }
	}
      else
	{
#endif
	  this->EigenvalueMatrix = RealDiagonalMatrix(TmpHilbertSpaceDimension);
	  if (this->EvaluateEigenvectors == false)
	    {
	      HRep.Diagonalize(this->EigenvalueMatrix);
	    }
	  else
	    {
	      RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
	      for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
		TmpEigenvector(l, l) = 1.0;
	      HRep.Diagonalize(this->EigenvalueMatrix, TmpEigenvector);
	      Basis.Multiply(TmpEigenvector);
	    }
	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
	    {
	      WriteResult(file, (this->EigenvalueMatrix[j] - this->EnergyShift));
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
      WriteResult(file, (HRep(0, 0)  - this->EnergyShift));
    }
  for (int j= 0; j < TmpHilbertSpaceDimension; ++j)
    delete[] VectorFileNames[j];
  delete[] VectorFileNames;
}

