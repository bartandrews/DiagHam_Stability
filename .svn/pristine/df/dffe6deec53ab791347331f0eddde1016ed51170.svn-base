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


#ifndef QHEONDISKMAINTASK_H
#define QHEONDISKMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"


class AbstractQHEHamiltonian;
class AbstractHilbertSpace;
class OptionManager;


class QHEOnDiskMainTask: public AbstractMainTask
{

 protected:

  // name of the file where results have to be stored
  char* OutputFileName;

  // twice the total momentum value of the system
  int LValue;

  // energy shift applied to the hamiltonian
  double EnergyShift;

  // pointer to the current Hamiltonian
  AbstractQHEHamiltonian* Hamiltonian;
  // pointer to the current Hilbert space
  AbstractHilbertSpace* Space;

  // name of the file where hamiltonian precalculations have to be saved (null if no precalculation has to be saved)
  char* SavePrecalculationFileName;
  // maximum Hilbert space dimension for which full diagonalization is applied
  int FullDiagonalizationLimit;
  // enable Lanczos disk resume capabilities
  bool DiskFlag;
  // resume from disk datas
  bool ResumeFlag;
  // enable block Lanczos algorithm
  bool BlockLanczosFlag;
  // size of the blocks used in the block Lanczos algorithm
  int SizeBlockLanczos;
  // number of eigenvalues to evaluate 
  int NbrEigenvalue;
  // number of lanczos iteration (for the current run)
  int NbrIterLanczos;
  // maximum time allowed for Lanczos iterations (in seconds)
  int MaximumAllowedTime;
  // maximum number of Lanczos iteration
  int MaxNbrIterLanczos;
  // maximum number of vector in RAM during Lanczos iteration
  int VectorMemory;
  // force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1
  bool FullReorthogonalizationFlag;
  // evaluate eigenstates
  bool EvaluateEigenvectors;
  // prefix to add to the name of each file that will contain an eigenvector
  char* EigenvectorFileName;
  // evaluate Lanczos convergence from eigenstate convergence
  bool EigenvectorConvergence;
  // show time spent for each Lanczos iteration
  bool ShowIterationTime;
  // name of the file that contains the vector to use as initial vector for the Lanczos algorithm (null if a random vector has to be picked)
  char* InitialVectorFileName;
  // name of the file that describes the set of vectors to use as initial set of vectors for the block Lanczos algorithm (null if a random vectors have to be picked)  
  char* InitialBlockVectorFileName;
  // allow to only run a given number of Lanczos iterations
  bool PartialLanczos;
  // use LAPACK libraries instead of DiagHam libraries
  bool LapackFlag;
  // name of the file that contains the vector files used to describe the Hilbert subspace
  char* ReducedHilbertSpaceDescription;
  // name of the file where the transformation basis from the Hilbert subspace to the diagonal basis is stored
  char* ReducedHilbertSpaceExportTransformation;
  // show the hamiltonian
  bool ShowHamiltonian;
  // compute hamiltonian mean value for each eigenvalue
  bool ComputeEnergyFlag;
  // define Lanczos precision for eigenvalues (0 if automatically defined by the program)
  double LanczosPrecision;
  // use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm to get the ground state
  bool FastDiskFlag;
  // indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeFastDiskFlag;
  // flag that indicates if it the first time the main task is used
  bool FirstRun;
  // name of the file that contains a optional set of vectors to which eigenstates have to be orthogonal
  char* LanczosReorthogonalization;

 public:

  // default constructor
  //
  QHEOnDiskMainTask ();

  // constructor
  //  
  // options = pointer to the options managers containing all running options
  // space = pointer to the current Hilbert space
  // hamiltonian = pointer to the current Hamiltonian
  // lValue = twice the total momentum value of the system
  // shift = energy shift that is applied to the hamiltonian
  // outputFileName = name of the file where results have to be stored
  // firstRun = flag that indicates if it the first time the main task is used
  // eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector
  QHEOnDiskMainTask(OptionManager* options, AbstractHilbertSpace* space, 
		      AbstractQHEHamiltonian* hamiltonian, int lValue, double shift, char* outputFileName,
		      bool firstRun = true, char* eigenvectorFileName = 0);
  
  // destructor
  //  
  ~QHEOnDiskMainTask();
  
  // set architecture binded to the task
  // 
  // architecture = pointer to the architecture to use
  void SetArchitecture(AbstractArchitecture* architecture);

  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  int ExecuteMainTask();

};

#endif
