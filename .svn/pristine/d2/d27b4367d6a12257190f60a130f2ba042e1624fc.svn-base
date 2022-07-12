////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
// class of main task to compute the spectrum of a non symmetric real matrix  //
//                                                                            //
//                        last modification : 17/02/2016                      //
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


#ifndef GENERICNONSYMMETRICMAINTASK_H
#define GENERICNONSYMMETRICMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"
#include "Vector/ComplexVector.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include <iostream>


using std::ofstream;


class AbstractHamiltonian;
class OptionManager;


class GenericNonSymmetricMainTask: public AbstractMainTask
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

  // if non zero, store the E matrix spectrum in this file
  char* EigenstateFileName;
  // prefix to add to the name of each file that will contain an eigenvector
  char* EigenvectorFileName;
  // an optional header that can be added to EigenvectorFileName
  char* EigenstateFileHeader;

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
  GenericNonSymmetricMainTask(OptionManager* options, AbstractHamiltonian* hamiltonian, 
			 int nbrEigenvalues, bool computeEigenstates, bool leftFlag, double sortingError, double energyShift = 0.0,
			 char* eigenstateFileName = 0, char* eigenstateFileHeader = 0, char* eigenvectorFileName = 0);
  
  // destructor
  //  
  ~GenericNonSymmetricMainTask();
  
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
  
};

// get the E matrix eigenvalues
//
// return value = pointer to the eigenvalue array

inline Complex* GenericNonSymmetricMainTask::GetEigenvalues()
{
  return this->Eigenvalues;
}

// get the E matrix eigenstates
//
// return value = pointer to the eigenstate array

inline ComplexVector* GenericNonSymmetricMainTask::GetEigenstates()
{
  return this->Eigenstates;
}

#endif
