////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                       class author: Gunnar Moeller                         //
//                                                                            //
//                    class of qhe on a lattice main task                     //
//                                                                            //
//                        last modification : 21/02/2008                      //
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

#ifdef HAVE_ARPACK
#include "arcomp.h"
#include "arrscomp.h"
#include "../examples/matprod/complex/cmatrixa.h"
#include "../examples/reverse/complex/rcompsol.h"
#endif


#include "MainTask/QHEOnLatticeMainTask.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Hamiltonian/AbstractQHEHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Options/Options.h"

#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

#ifdef HAVE_ARPACK
template<class T>
  void Test(T type)
  {
    
    // Defining a complex matrix.
    
    CompMatrixA<T> A(10); // n = 10*10.
    
    // Creating a complex eigenvalue problem and defining what we need:
    // the four eigenvectors of A with largest magnitude.
    
    ARrcCompStdEig<T> prob(A.ncols(), 4L);
    
    // Finding an Arnoldi basis.
    
    while (!prob.ArnoldiBasisFound()) {
      
      // Calling ARPACK FORTRAN code. Almost all work needed to
      // find an Arnoldi basis is performed by TakeStep.
      
      prob.TakeStep();
      
      if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {
	
	// Performing matrix-vector multiplication.
	// In regular mode, w = Av must be performed whenever
	// GetIdo is equal to 1 or -1. GetVector supplies a pointer
	// to the input vector, v, and PutVector a pointer to the
	// output vector, w.
	
	A.MultMv(prob.GetVector(), prob.PutVector());
	
      }
      
    }
    
    // Finding eigenvalues and eigenvectors.
    
    prob.FindEigenvectors();
    
    // Printing solution.
    
    Solution(prob);
    
  } // Test.
#endif


// constructor
//  
// options = pointer to the options managers containing all running options
// space = pointer to the current Hilbert space
// hamiltonian = pointer to the current Hamiltonian
// nbrFluxQuanta = number of flux quanta piercing the lattice (for output)
// shift = energy shift that is applied to the hamiltonian
// outputFileName = name of the file where results have to be stored
// firstRun = flag that indicates if it the first time the main task is used
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector
// ky = many-body momentum in y-direction

QHEOnLatticeMainTask::QHEOnLatticeMainTask(OptionManager* options, AbstractHilbertSpace* space, 
					   AbstractQHEHamiltonian* hamiltonian, int nbrFluxQuanta, double shift,
					   char* outputFileName, bool firstRun, char* eigenvectorFileName, int ky)
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
  this->NbrFluxQuanta = nbrFluxQuanta;
  this->Ky=ky;
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
  this->VectorMemory = ((SingleIntegerOption*) (*options)["nbr-vector"])->GetInteger();
  this->SavePrecalculationFileName = ((SingleStringOption*) (*options)["save-precalculation"])->GetString();
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
  if ((*options)["partial-lanczos"] != 0)
    {
      this->PartialLanczos = ((BooleanOption*) (*options)["partial-lanczos"])->GetBoolean();
    }
  else
    {
      this->PartialLanczos = false;
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
    }
  else
    {
      this->ReducedHilbertSpaceDescription = 0;
    }
  if ((*options)["get-hvalue"] != 0)
    {
      this->ComputeEnergyFlag = ((BooleanOption*) (*options)["get-hvalue"])->GetBoolean();
    }
  else
    {
      this->ComputeEnergyFlag = false;
    }
  if (((*options)["show-hamiltonian"] != 0) && (((BooleanOption*) (*options)["show-hamiltonian"])->GetBoolean() == true))
    {
      HermitianMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
      this->Hamiltonian->GetHamiltonian(HRep);
      cout << HRep << endl;
    }  
  if (((*options)["lanczos-precision"] != 0) && (((SingleDoubleOption*) (*options)["lanczos-precision"])->GetDouble() > 0))
    {
      this->LanczosPrecision = ((SingleDoubleOption*) (*options)["lanczos-precision"])->GetDouble();
    }
  else
    {
      this->LanczosPrecision = 0.0;
    }
  if (((*options)["fast-disk"] != 0) && (this->NbrEigenvalue == 1) && (this->EvaluateEigenvectors == true))
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
  this->FirstRun = firstRun;
#ifdef HAVE_ARPACK
  ArpackFlag = options->GetBoolean("use-arpack");
#endif
}
 
// destructor
//  

QHEOnLatticeMainTask::~QHEOnLatticeMainTask()
{
  delete[] this->OutputFileName;
  if (this->EigenvectorFileName != 0)
    delete[] this->EigenvectorFileName;
}
  
// set architecture bound to the task
// 
// architecture = pointer to the architecture to use

void QHEOnLatticeMainTask::SetArchitecture(AbstractArchitecture* architecture)
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

int QHEOnLatticeMainTask::ExecuteMainTask()
{
  ofstream File;
  if (this->FirstRun == true)
    {
      File.open(this->OutputFileName, ios::binary | ios::out);
      this->FirstRun = false;
      File << "# NbrFlux";
      if (Space->GetHilbertSpaceAdditionalSymmetry()&ParticleOnLattice::YTranslations)
	File <<" Ky";
      File <<" E";
      File << endl;
    }
  else
    {
      File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
    }
  File.precision(14);
  cout.precision(14);
  cout << "----------------------------------------------------------------" << endl;
  if (this->ReducedHilbertSpaceDescription != 0)
    {
      this->DiagonalizeInHilbertSubspace(this->ReducedHilbertSpaceDescription, File);
      cout << "----------------------------------------------------------------" << endl;
      File.close(); 
      return 0;
    }
  cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
  if (this->SavePrecalculationFileName != 0)
    {
      this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
    }
  if (this->Hamiltonian->GetHilbertSpaceDimension()==0) return 0;
  if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
    {
      HermitianMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	{
#ifdef __LAPACK__
	  if (this->LapackFlag == true)
	    {
	      cout << "Need to implement diagonalization with Lapack for complex matrices..." << endl;
// 	      RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
// 	      if ((this->EvaluateEigenvectors == false) && (this->ComputeLValueFlag == false))
// 		{
// 		  HRep.LapackDiagonalize(TmpDiag);
// 		  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
//	            this->WriteResult(File,TmpDiag[j] - this->EnergyShift);
// 		}
// 	      else
// 		{
// 		  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
// 		  HRep.LapackDiagonalize(TmpDiag, Q);
// 		  if (this->EvaluateEigenvectors == true)
// 		    {
// 		      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
// 		      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
// 		      for (int j = 0; j < this->NbrEigenvalue; ++j)
// 			{
// 			  this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
// 			  sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
// 			  Q[j].WriteVector(TmpVectorName);
// 			  cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
// 			}
// 		      cout << endl;			  
// 		      delete[] TmpVectorName;
// 		    }
// 		  if (this->ComputeLValueFlag == false)
// 		    for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
// 		      {
//      	        this->WriteResult(File,TmpDiag[j] - this->EnergyShift);
// 			if (this->ComputeEnergyFlag == true)
// 			  {
// 			    RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
// 			    this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
// 			    File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
// 			  }
// 			File << endl;
// 		      }
// 		  else
// 		    {
// 		      ParticleOnSphereSquareTotalMomentumOperator Oper((ParticleOnSphere*) Space, this->LzMax);
// 		      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
// 			{
//	                  this->WriteResult(File,TmpDiag[j] - this->EnergyShift);	      
// 			  if (this->ComputeEnergyFlag == true)
// 			    {
// 			      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
// 			      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
// 			      File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
// 			    }
// 			  double TmpMomentum = Oper.MatrixElement(Q[j], Q[j]).Re;
// 			  File << " "  << TmpMomentum << " " << (0.5 * (sqrt ((4.0 * TmpMomentum) + 1.0) - 1.0)) << endl;
// 			}
// 		    }
// 		}
	    }
	  else
	    {
#endif
	      // solution of eigenvalue problem of HRep with internal DiagHam routines
	      //cout << "HRep=" << endl << HRep << endl;
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
		  // cout.flush();
// 		  cout << endl << "eigenvector_0:"<<endl;
// 		  for (int i=0; i<this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
// 		    cout << Q[0][i].Re<<"+I*"<<Q[0][i].Im<<endl;		  
		  cout << endl;
		  // this->Hamiltonian->LowLevelMultiply(Q[0], TmpEigenvector);
// 		  cout.flush();
// 		  cout << endl << "H*eigenvector_0:"<<endl;
// 		  for (int i=0; i<this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
// 		    cout << TmpEigenvector[i].Re<<"+I*"<<TmpEigenvector[i].Im<<endl;		  
// 		  cout << endl;
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
#ifdef HAVE_ARPACK
      if (ArpackFlag)
	{
	  
	  // Solving a double precision problem with n = 100.
	  cout << "First Test" << endl;
	  cout.flush();
	  Test((double)0.0);

	  cout << "Second Test" << endl;
	  cout.flush();
	  // Creating a complex eigenvalue problem and defining what we need:
	  // the four eigenvalues of A with smallest real part.
	  int Dimension=this->Hamiltonian->GetHilbertSpaceDimension();
	  ComplexVector InputVector(Dimension);
	  ComplexVector OutputVector(Dimension);
	  
	  ARrcCompStdEig<double> prob(Dimension, 4L); //this->NbrEigenvalue, "SR");
	  // Finding an Arnoldi basis.

	  while (!prob.ArnoldiBasisFound()) {
	    
	    // Calling ARPACK FORTRAN code. Almost all work needed to
	    // find an Arnoldi basis is performed by TakeStep.
	    
	    prob.TakeStep();
	    
	    if ((prob.GetIdo() == 1)||(prob.GetIdo() == -1)) {
	      
	      // Performing matrix-vector multiplication.
	      // In regular mode, w = Av must be performed whenever
	      // GetIdo is equal to 1 or -1. GetVector supplies a pointer
	      // to the input vector, v, and PutVector a pointer to the
	      // output vector, w.

	      arcomplex<double> *ArpackInputVector = prob.GetVector();
	      arcomplex<double> *ArpackOutputVector = prob.PutVector();
	      cout << "ArpackInputVector="<<endl;
	      
	      for (int i=0; i<Dimension; ++i)
		{
		  cout << ArpackInputVector[i]<<"\t"<<ArpackOutputVector[i]<<endl;
		  InputVector.Re(i) = real(ArpackInputVector[i]);
		  InputVector.Im(i) = imag(ArpackInputVector[i]);
		}

	      cout << "InputVector="<<InputVector<<endl;
	      
	      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &InputVector, &OutputVector);
	      Operation1.ApplyOperation(this->Architecture);

	      for (int i=0; i<Dimension; ++i)
		{
		  ArpackOutputVector[i] = arcomplex<double>(OutputVector.Re(i),OutputVector.Im(i));
		}
	      
	      
	    }
	    
	  }
	  
	  // Finding eigenvalues and eigenvectors.
	  
	  prob.FindEigenvectors();
	  
	  // Printing solution.
	  
	  Solution(prob);
	}
      else
#endif	
	{  
	  AbstractLanczosAlgorithm* Lanczos;
	  if ((this->NbrEigenvalue == 1) && (this->FullReorthogonalizationFlag == false))
	    {
	      if (this->DiskFlag == false)
		if (this->EvaluateEigenvectors == true)
		  //Lanczos = new ComplexBasicLanczosAlgorithmWithGroundState(this->Architecture, this->MaxNbrIterLanczos);// replaced by more elaborate algorithm with fast-disk option -> still need to check that one
		  {
		    cout << "Using un-tested ComplexBasicLanczosAlgorithmWithGroundStateFastDisk!"<<endl;
		    Lanczos = new ComplexBasicLanczosAlgorithmWithGroundStateFastDisk(this->Architecture, this->MaxNbrIterLanczos , this->FastDiskFlag, this->ResumeFastDiskFlag);
		  }
		else
		  Lanczos = new ComplexBasicLanczosAlgorithm(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	      else
		if (this->EvaluateEigenvectors == true)
		  {
		    cout << "Complex Lanczos Algorithm with GroundState and Disk Storage not implemented!"<<endl;
		    exit(1);
		    //Lanczos = new ComplexBasicLanczosAlgorithmWithGroundStateDiskStorage(this->Architecture, this->NbrIterLanczos, this->MaxNbrIterLanczos);
		  }
		else
		  Lanczos = new ComplexBasicLanczosAlgorithmWithDiskStorage(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (this->DiskFlag == false)
		{
		  Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm (this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
		}
	      else
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage (this->Architecture, this->NbrEigenvalue, this->VectorMemory, this->MaxNbrIterLanczos);
	    }
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
							   ((this->DiskFlag == false) && ((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
							    ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)))))
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
	       ((this->DiskFlag == false) && ((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
		((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)))))
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
			      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, i);
			    }
			  else
			    {
			      sprintf (TmpVectorName, "%s.%d.part.vec", this->EigenvectorFileName, i);		  
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
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	  delete Lanczos;
	}
    }
  cout << "----------------------------------------------------------------" << endl;
  File.close(); 
  return 0;
}

// do the Hamiltonian diagonalization in a given Hilbert subspace
//
// subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
// file = reference on the output file stream where eigenvalues have to be stored

void QHEOnLatticeMainTask::DiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file)
{
//   ConfigurationParser ReducedBasis;
//   if (ReducedBasis.Parse(subspaceDescription) == false)
//     {
//       ReducedBasis.DumpErrors(cout) << endl;
//       return;
//     }
//   int TmpHilbertSpaceDimension;
//   char** VectorFileNames;
//   if (ReducedBasis.GetAsStringArray("Basis", ' ', VectorFileNames, TmpHilbertSpaceDimension) == false)
//     {
//       cout << "Vectors are not defined or have a wrong value in " << subspaceDescription << endl;
//       return;
//     }
//   RealMatrix Basis (this->Space->GetHilbertSpaceDimension(), TmpHilbertSpaceDimension);
//   char* DirectoryName = ReducedBasis["Directory"];
//   char* TmpName;
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     {
//       TmpName = VectorFileNames[i];
//       if (DirectoryName != 0)
// 	{
// 	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
// 	}
//       cout << TmpName << endl;
//       if (Basis[i].ReadVector(TmpName) == false)
// 	{
// 	  cout << "error while reading " << TmpName << endl;
// 	  if (DirectoryName != 0)
// 	    delete[] TmpName;
// 	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
// 	    delete[] VectorFileNames[j];
// 	  delete[] VectorFileNames;
// 	  return;
// 	}
//       if (DirectoryName != 0)
// 	delete[] TmpName;
//     }
//   RealSymmetricMatrix HRep (TmpHilbertSpaceDimension);
//   RealVector* TmpVectors = new RealVector[TmpHilbertSpaceDimension];
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     {
//       RealVector TmpVector (Basis[0].GetVectorDimension(), true);
//       VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Basis[i]), &TmpVector);
//       Operation1.ApplyOperation(this->Architecture);      
//       TmpVectors[i] = TmpVector;
//     }
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     {
//       for (int j = i; j < TmpHilbertSpaceDimension; ++j)
// 	{
// 	  HRep(i ,j) = Basis[j] * TmpVectors[i];
// 	}
//     }
//   delete[] TmpVectors;
//   if (TmpHilbertSpaceDimension > 1)
//     {
// #ifdef __LAPACK__
//       if (this->LapackFlag == true)
// 	{
// 	  RealDiagonalMatrix TmpDiag (TmpHilbertSpaceDimension);
// 	  if (this->EvaluateEigenvectors == false)
// 	    {
// 	      HRep.LapackDiagonalize(TmpDiag);
// 	    }
// 	  else
// 	    {
// 	      RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
// 	      for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
// 		TmpEigenvector(l, l) = 1.0;
// 	      HRep.LapackDiagonalize(TmpDiag, TmpEigenvector);
// 	      Basis.Multiply(TmpEigenvector);
// 	    }
// 	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
// 	    {
// 	      file << (this->LValue/ 2) << " " << (TmpDiag[j] - this->EnergyShift) << endl;
// 	    }
// 	}
//       else
// 	{
// #endif
// 	  RealTriDiagonalSymmetricMatrix TmpTriDiag (TmpHilbertSpaceDimension);
// 	  if (this->EvaluateEigenvectors == false)
// 	    {
// 	      HRep.Householder(TmpTriDiag, 1e-7);
// 	      TmpTriDiag.Diagonalize();
// 	      TmpTriDiag.SortMatrixUpOrder();
// 	    }
// 	  else
// 	    {
// 	      RealMatrix TmpEigenvector (TmpHilbertSpaceDimension, TmpHilbertSpaceDimension, true);	      
// 	      for (int l = 0; l < TmpHilbertSpaceDimension; ++l)
// 		TmpEigenvector(l, l) = 1.0;
// 	      HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
// 	      TmpTriDiag.Diagonalize(TmpEigenvector);
// 	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
// 	      Basis.Multiply(TmpEigenvector);
// 	    }
// 	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
// 	    {
// 	      file << (this->LValue/ 2) << " " << (TmpTriDiag.DiagonalElement(j) - this->EnergyShift) << endl;
// 	    }
// #ifdef __LAPACK__
// 	}
// #endif
//       if (this->EvaluateEigenvectors == true)
// 	{
// 	  for (int j = 0; j < TmpHilbertSpaceDimension; ++j)
// 	    {
// 	      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
// 	      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
// 	      Basis[j].WriteVector(TmpVectorName);
// 	    }
// 	}
//     }
//   else
//     {
//       file << (this->LValue/ 2) << " " << (HRep(0, 0)  - this->EnergyShift) << endl;
//     }
//   for (int j= 0; j < TmpHilbertSpaceDimension; ++j)
//     delete[] VectorFileNames[j];
//   delete[] VectorFileNames;
}


// add optiongroup with options related to this module to the given OptionManager
//
void QHEOnLatticeMainTask::AddOptionGroup(OptionManager *optionManager)
{  
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  (*optionManager) += LanczosGroup;
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
	"maximum Hilbert space dimension for which full diagonalization is applied", 300, true, 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup) += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "get-hvalue", "show energy expectation value for eigenstates", false);
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new  BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*LanczosGroup) += new  BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");
  (*LanczosGroup) += new  BooleanOption ('\n',"show-hamiltonian", "show Hamiltonian matrix, and exit");
#ifdef HAVE_ARPACK
  (*LanczosGroup) += new  BooleanOption ('\n',"use-arpack","use ARPACK routines for Lanczos algorithm");
#endif
}

// write a line of output to the results file
//
// file = stream to write to
// value = numerical value to be printed after columns for flux and momentum (if defined)
// terminate = indicate if line should be terminated with endl
void QHEOnLatticeMainTask::WriteResult(ofstream& file, double value, bool terminate)
{
  file << NbrFluxQuanta<<" ";
  if (this->Ky>=0) file << this->Ky << " ";
  file << value;
  // cout <<"wrote value "<<value<<" to file"<<endl;
  if (terminate)
    file << endl;
}
