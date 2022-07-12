////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of main task to compute the spectrum of                 //
//                      a non Hermitian complex matrix                        //
//                                                                            //
//                        last modification : 07/03/2016                      //
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


#ifndef GENERICNONHERMITIANMAINTASK_H
#define GENERICNONHERMITIANMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"
#include "Vector/ComplexVector.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include <iostream>


using std::ofstream;


class AbstractHamiltonian;
class OptionManager;


class GenericNonHermitianMainTask: public AbstractMainTask
{

 protected:

  // pointer to the current Hamiltonian
  AbstractHamiltonian* Hamiltonian;

  // number of eigenvalues to compute
  int NbrEigenvalues;
  // array where the eigenvalues are stored
  Complex* Eigenvalues;
  // true if the eigenstates have to be computed
  bool ComputeEigenstates;
  // array where the eigenstates are stored
  ComplexVector* Eigenstates;
  // Compute the left eigenstates if true 
  bool LeftFlag;
    

  // string describing present subspace to prepend in output
  char *SubSpaceStr;
  // legend indicating contents of subspaceStr to include in output file
  char *SubSpaceLegend;

  // numerical error when sorting real numbers
  double SortingError;

  // use the power method instead of Arnoldi
  bool PowerMethodFlag;

  // energy shift that has been applied to the Hamiltonian
  double EnergyShift;

  // show time spent for each Lanczos iteration
  bool ShowIterationTime;
  // enable Arnoldi disk resume capabilities
  bool DiskFlag;
  // resume from disk datas
  bool ResumeFlag;
  // maximum Hilbert space dimension for which full diagonalization is applied
  int FullDiagonalizationLimit;
  // maximum number of Arnoldi iteration
  int MaxNbrIterArnoldi;
  // amount of memory (in Mbytes) that can be used for the Arnoldi algorithm with disk storage
  long ArnoldiMemory;
  // use the implicitly restarted Arnoldi algorithm
  bool ImplicitlyRestartedFlag;
  // sort the eigenvalues only with respect to their real part
  bool SortEigenvalueRealPartFlag;
  // eigenvalue precision
  double EigenvaluePrecision;
  // flag that indicates if it the first time the main task is used
  bool FirstRun;

  // if non zero, store the E matrix spectrum in this file
  char* EigenstateFileName;
  // prefix to add to the name of each file that will contain an eigenvector
  char* EigenvectorFileName;

 public:

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
  GenericNonHermitianMainTask(OptionManager* options, AbstractHamiltonian* hamiltonian, 
			      int nbrEigenvalues, bool computeEigenstates, bool leftFlag, double sortingError,  const char* subSpaceStr, const char* subSpaceLegend = 0, double energyShift = 0.0, bool firstRun=true,
			 char* eigenstateFileName = 0, char* eigenvectorFileName = 0);
  
  // destructor
  //  
  ~GenericNonHermitianMainTask();
  
  // set architecture binded to the task
  // 
  // architecture = pointer to the architecture to use
  virtual void SetArchitecture(AbstractArchitecture* architecture);

  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  virtual int ExecuteMainTask();

  // get the E matrix eigenvalues
  //
  // return value = pointer to the eigenvalue array
  virtual Complex* GetEigenvalues();

  // get the E matrix eigenstates
  //
  // return value = pointer to the eigenstate array
  virtual ComplexVector* GetEigenstates();

 protected:

  // write a line of output to the results file
  //
  // file = stream to write to
  // value = numerical value to be printed after columns for flux and momentum (if defined)
  // terminate = indicate if line should be terminated with endl
  virtual void WriteResult(ofstream& file, Complex value, bool terminate=true);
  
};

// get the E matrix eigenvalues
//
// return value = pointer to the eigenvalue array

inline Complex* GenericNonHermitianMainTask::GetEigenvalues()
{
  return this->Eigenvalues;
}

// get the E matrix eigenstates
//
// return value = pointer to the eigenstate array

inline ComplexVector* GenericNonHermitianMainTask::GetEigenstates()
{
  return this->Eigenstates;
}

#endif
