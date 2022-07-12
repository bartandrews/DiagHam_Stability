////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of main task to compute the E matrix properties         //
//                                                                            //
//                        last modification : 16/01/2013                      //
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

#include "MainTask/GenericNonSymmetricMainTask.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "HilbertSpace/AbstractHilbertSpace.h"

#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicComplexArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/PowerMethodAlgorithm.h"
#include "LanczosAlgorithm/ImplicitlyRestartedArnoldiAlgorithm.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <fstream>
#include <sys/time.h>


using std::ofstream;
using std::cout;
using std::endl;
using std::ios;



// constructor
//  
// options = pointer to the options managers containing all running options
// hamiltonian = pointer to the current Hamiltonian
// nbrEigenvalues = number of eigenvalues to compute
// computeEigenstates = true if the eigenstates have to be computed
// leftFlag = compute the left eigenstates if true 
// sortingError = numerical error when sorting real numbers
// energyShift = energy shift that has been applied to the Hamiltonian
// eigenstateFileName = if non zero, store the E matrix spectrum in this file
// eigenstateFileHeader = an optional header that can be added to EigenvectorFileName
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector 
//                       (eigenvectors are stored only if eigenvectorFileName is non zero)

GenericNonSymmetricMainTask::GenericNonSymmetricMainTask(OptionManager* options, AbstractHamiltonian* hamiltonian, 
					       int nbrEigenvalues, bool computeEigenstates, bool leftFlag, double sortingError,
					       double energyShift, char* eigenstateFileName, char* eigenstateFileHeader, char* eigenvectorFileName)
{
  this->Hamiltonian = hamiltonian;
  this->NbrEigenvalues = nbrEigenvalues;
  cout<<this->NbrEigenvalues<<endl;
  cout <<this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension()<<endl;
  if (this->NbrEigenvalues > this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension())
    {
      this->NbrEigenvalues = this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension();
    }
  this->Eigenvalues = new Complex[this->NbrEigenvalues];
  this->ComputeEigenstates = computeEigenstates;
  if (this->ComputeEigenstates == true)
    {
      this->Eigenstates = new ComplexVector[this->NbrEigenvalues];
      if (eigenvectorFileName != 0)
	{
	  this->EigenvectorFileName = new char [strlen(eigenvectorFileName) + 1];
	  strcpy (this->EigenvectorFileName, eigenvectorFileName);
	}
      else
	{
	  this->EigenvectorFileName = 0;
	}
    }
  else
    {
      this->Eigenstates = 0;
      this->EigenvectorFileName = 0;
    }
  if (eigenstateFileName != 0)
    {
      this->EigenstateFileName = new char [strlen(eigenstateFileName) + 1];
      strcpy (this->EigenstateFileName, eigenstateFileName);
    }
  else
    {
      this->EigenstateFileName = 0;      
    }
  if (eigenstateFileHeader != 0)
    {
      this->EigenstateFileHeader = new char [strlen(eigenstateFileName) + 1];
      strcpy (this->EigenstateFileHeader, eigenstateFileName);
    }
  else
    {
       this->EigenstateFileHeader = 0;      
   }
  this->LeftFlag = leftFlag;    
  this->SortingError = sortingError;

  this->SortEigenvalueRealPartFlag = false;
  if ((*options)["sort-real"] != 0)
    {
      this->SortEigenvalueRealPartFlag = options->GetBoolean("sort-real");
    }
  this->EigenvaluePrecision = MACHINE_PRECISION;
  if (((*options)["arnoldi-precision"] != 0) && (options->GetDouble("arnoldi-precision") > 0.0))
    {
      this->EigenvaluePrecision = options->GetDouble("arnoldi-precision");
    }
  this->ShowIterationTime = false;
  if ((*options)["show-itertime"] != 0)
    {
      this->ShowIterationTime = options->GetBoolean("show-itertime");
    }
  this->ResumeFlag = false;
  if ((*options)["resume"] != 0)
    {
      this->ResumeFlag = options->GetBoolean("resume");
    }
  this->DiskFlag = false;
  if ((*options)["disk"] != 0)
    {
      this->DiskFlag = options->GetBoolean("disk");
    }
  this->FullDiagonalizationLimit = 1000;
  if ((*options)["full-diag"] != 0)
    {
      this->FullDiagonalizationLimit = options->GetInteger("full-diag");
    }
  this->MaxNbrIterArnoldi = 3000;
  if ((*options)["iter-max"] != 0)
    {
      this->MaxNbrIterArnoldi = options->GetInteger("iter-max");
    }
  this->ArnoldiMemory = 0l;
  if ((*options)["arnoldi-memory"] != 0)
    {
      this->ArnoldiMemory = options->GetInteger("arnoldi-memory");
    }
  this->ImplicitlyRestartedFlag = false;
  if ((*options)["implicitly-restarted"] != 0)
    {
      this->ImplicitlyRestartedFlag = options->GetBoolean("implicitly-restarted");
    }

  this->EnergyShift = energyShift;
  this->PowerMethodFlag = false;
  if ((*options)["power-method"] != 0)
    {
      this->PowerMethodFlag = options->GetBoolean("power-method");
    }
}

// destructor
//  

GenericNonSymmetricMainTask::~GenericNonSymmetricMainTask()
{
  delete[] this->Eigenvalues;
  if (this->ComputeEigenstates == true)
    {
      delete[] this->Eigenstates;
    }
}

// set architecture bound to the task
// 
// architecture = pointer to the architecture to use

void GenericNonSymmetricMainTask::SetArchitecture(AbstractArchitecture* architecture)
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
  
int GenericNonSymmetricMainTask::ExecuteMainTask()
{
  double* TmpRealPart = new double [this->NbrEigenvalues];
  int* TmpIndices = new int[this->NbrEigenvalues];
  Complex* TmpEigenvalues = 0;
  ComplexVector* TmpEigenstates = 0;
  int LocalNbrEigenvalues = this->NbrEigenvalues;

  if (this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
    {	  
      RealMatrix HRepresentation (this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension(), 
				  this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension());
      this->Hamiltonian->GetHamiltonian(HRepresentation);

      ComplexDiagonalMatrix TmpDiag (HRepresentation.GetNbrRow(), true);  
      LocalNbrEigenvalues = this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension();
      TmpEigenvalues = new Complex [LocalNbrEigenvalues];

      if (this->ComputeEigenstates == true)
	{
	  ComplexMatrix TmpEigenstateMatrix (HRepresentation.GetNbrRow(), HRepresentation.GetNbrRow());  
	  HRepresentation.LapackDiagonalize(TmpDiag, TmpEigenstateMatrix, true);
	  TmpDiag.SortMatrixDownOrder(TmpEigenstateMatrix, true);
	  TmpEigenstates = new ComplexVector [this->NbrEigenvalues];
	  for (int i = 0; i < LocalNbrEigenvalues; ++i)
	    {
	      TmpEigenvalues[i] = TmpDiag[i];
	    }
	  for (int i = 0; i < this->NbrEigenvalues; ++i)
	    {
	      TmpEigenstates[i] = TmpEigenstateMatrix[i];
	    }	  
	}
      else
	{
	  HRepresentation.LapackDiagonalize(TmpDiag);
	  TmpDiag.SortMatrixDownOrder(true);
	  for (int i = 0; i < LocalNbrEigenvalues; ++i)
	    {
	      TmpEigenvalues[i] = TmpDiag[i];
	    }
	}
    }
  else
    {	 
      if (this->PowerMethodFlag == false)
	{
	  AbstractLanczosAlgorithm* Arnoldi = 0;
	  if (this->DiskFlag)
	    {
	      long TmpMemory = (ArnoldiMemory << 17) / this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension();
	      if (TmpMemory == 0)
		TmpMemory = 1;
	      Arnoldi = new BasicArnoldiAlgorithmWithDiskStorage (Architecture, this->NbrEigenvalues, this->MaxNbrIterArnoldi, true, false, 
								  this->ResumeFlag, TmpMemory, false);
	    }
	  else
	    {
	      if (this->ImplicitlyRestartedFlag == false)
		{
		  Arnoldi = new BasicArnoldiAlgorithm (Architecture, this->NbrEigenvalues, this->MaxNbrIterArnoldi, true, false, false, this->SortEigenvalueRealPartFlag);
		}
	      else
		{
		  int TmpMemory = (int) ((ArnoldiMemory << 17) / this->Hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension());
		  if (TmpMemory < (this->NbrEigenvalues + 3))
		    {
		      TmpMemory =  this->NbrEigenvalues + 3;
		    }
		  Arnoldi = new ImplicitlyRestartedArnoldiAlgorithm (Architecture, this->NbrEigenvalues, TmpMemory, this->NbrEigenvalues, this->MaxNbrIterArnoldi, true, false, false);
		}
	    }
	  Arnoldi->SetEigenvaluePrecision(this->EigenvaluePrecision);
	  Arnoldi->SetHamiltonian(this->Hamiltonian);
	  Arnoldi->InitializeLanczosAlgorithm();
	  Arnoldi->RunLanczosAlgorithm(this->NbrEigenvalues);
	  while (Arnoldi->TestConvergence() == false)
	    {
	      timeval TotalStartingTime;
	      timeval TotalEndingTime;
	      if (this->ShowIterationTime == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      Arnoldi->RunLanczosAlgorithm(1);
	      if (this->ShowIterationTime == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "iteration done in " << Dt << "s" << endl;
		}
	    }
	  if (this->ComputeEigenstates)
	    TmpEigenstates = (ComplexVector*) Arnoldi->GetEigenstates(this->NbrEigenvalues);  
	  Arnoldi->GetEigenvalues(TmpEigenvalues, this->NbrEigenvalues);
	}
      else
	{
	  this->NbrEigenvalues = 1;
	  
	  PowerMethodAlgorithm* PowerMethod = new PowerMethodAlgorithm (Architecture, this->MaxNbrIterArnoldi);
	  
	  PowerMethod->SetHamiltonian(this->Hamiltonian);
	  PowerMethod->InitializeLanczosAlgorithm();
	  PowerMethod->RunLanczosAlgorithm(this->NbrEigenvalues);
	  while (PowerMethod->TestConvergence() == false)
	    {
	      timeval TotalStartingTime;
	      timeval TotalEndingTime;
	      if (this->ShowIterationTime == true)
		{
		  gettimeofday (&(TotalStartingTime), 0);
		}
	      PowerMethod->RunLanczosAlgorithm(1);
	      if (this->ShowIterationTime == true)
		{
		  gettimeofday (&(TotalEndingTime), 0);
		  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
					((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
		  cout << "iteration done in " << Dt << "s" << endl;
		}
	    }
	  if (this->ComputeEigenstates)
	    TmpEigenstates = (ComplexVector*) PowerMethod->GetEigenstates(this->NbrEigenvalues);  
	  PowerMethod->GetEigenvalues(TmpEigenvalues, this->NbrEigenvalues);
	}
    }
  for (int i = 0; i < this->NbrEigenvalues; ++i)
    {
      TmpRealPart[i] = TmpEigenvalues[i].Re;
      TmpIndices[i] = i;
    }
  SortArrayDownOrdering<int>(TmpRealPart, TmpIndices, this->NbrEigenvalues);
  for (int i = 1; i < this->NbrEigenvalues; ++i)
    {
      if ((fabs(TmpEigenvalues[TmpIndices[i - 1]].Re - TmpEigenvalues[TmpIndices[i]].Re) < this->SortingError) && 
	  (TmpEigenvalues[TmpIndices[i - 1]].Im < TmpEigenvalues[TmpIndices[i]].Im))
	{
	  int Tmp = TmpIndices[i - 1];
	  TmpIndices[i - 1] = TmpIndices[i];
	  TmpIndices[i] = Tmp;
	}
    }
  for (int i = 0; i < this->NbrEigenvalues; ++i)
    {
      if (this->ComputeEigenstates)
	this->Eigenstates[i] = TmpEigenstates[TmpIndices[i]];
      this->Eigenvalues[i] = TmpEigenvalues[TmpIndices[i]] - this->EnergyShift;
    }
  if (this->ComputeEigenstates)
    {
      for (int i = 0; i < this->NbrEigenvalues; ++i)
	{
	  ComplexVector TestE (this->Eigenstates[i].GetVectorDimension());
	  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->Eigenstates[i]), &TestE);
	  Operation1.ApplyOperation(this->Architecture);
	  cout << ((this->Eigenstates[i] * TestE) - this->EnergyShift) << " " << this->Eigenvalues[i] << endl;
	}
      if (this->EigenvectorFileName != 0)
	{
	  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
	  for (int j = 0; j < this->NbrEigenvalues; ++j)
	    {
	      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
	      this->Eigenstates[j].WriteVector(TmpVectorName);
	    }
	}
    }
  if (this->EigenstateFileName != 0)
    {
      ofstream File;
      File.open(this->EigenstateFileName, ios::binary | ios::out);
      if (this->EigenstateFileHeader != 0)
	File << "# "<< this->EigenstateFileHeader << endl;
      File << "# E |E| arg(E)/pi" << endl;
      File.precision(14);
      for (int i = 0; i < LocalNbrEigenvalues; ++i)
	{
	  File << TmpEigenvalues[i] << " " << Norm(TmpEigenvalues[i]) << " " << (Arg(TmpEigenvalues[i]) / M_PI) << endl;
	}
      File.close();
    }
  delete[] TmpRealPart;
  delete[] TmpIndices;
  delete[] TmpEigenstates;
  delete[] TmpEigenvalues;
  return 0;
}
