////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                         class author: Gunnar Moeller                       //
//                                                                            //
//                      class of qhe on lattice main task                     //
//                                                                            //
//                        last modification : 13/02/2008                      //
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


#ifndef GENERICCOMPLEXMAINTASK_H
#define GENERICCOMPLEXMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"
#include "MathTools/Complex.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include <iostream>


using std::ofstream;


class AbstractHamiltonian;
class AbstractHilbertSpace;
class OptionManager;


class GenericComplexMainTask: public AbstractMainTask
{

 protected:

  // name of the file where results have to be stored
  char* OutputFileName;

  // string describing present subspace to prepend in output
  char *SubspaceStr;
  // legend indicating contents of subspaceStr to include in output file
  char *SubspaceLegend;
  // energy shift applied to the hamiltonian
  double EnergyShift;

  // pointer to the current Hamiltonian
  AbstractHamiltonian* Hamiltonian;
  // pointer to the current Hilbert space
  AbstractHilbertSpace* Space;
  // pointer to Lanczos manager
  LanczosManager *AlgorithmManager;
  

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
  // evaluate all eigenstates
  bool EvaluateAllEigenvectors;
  // index of the first eigenstate to compute
  int FirstEigenstateIndex;
  // prefix to add to the name of each file that will contain an eigenvector
  char* EigenvectorFileName;
  // evaluate expectation value of energy for each trial state
  bool ComputeEnergyFlag;
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
  // use SCALAPACK libraries instead of DiagHam and/or Lapack libraries
  bool ScalapackFlag;
  // name of the file that contains the vector files used to describe the Hilbert subspace
  char* ReducedHilbertSpaceDescription;  
  // show the hamiltonian
  bool ShowHamiltonian;
  // show the hamiltonian in a friendly way, showing only non-zero matrix elements
  bool FriendlyShowHamiltonian;
  // when showing the hamiltonian in a friendly way, error below which a matrix is considered to be equal to zero
  double FriendlyShowHamiltonianError;  
  // define Lanczos precision for eigenvalues (0 if automatically defined by the program)
  double LanczosPrecision;
  // use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm to get the ground state
  bool FastDiskFlag;
  // indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeFastDiskFlag;
  // flag that indicates if it the first time the main task is used
  bool FirstRun;
  // flag that indicate if eigenstates have to be computed and saved at a given frequency
  int PartialEigenstateFlag;
  //  if non-zero, store the Hamiltonian as a binary file
  char* ExportBinaryHamiltonian;

 public:

  // constructor
  //  
  // options = pointer to the options managers containing all running options
  // space = pointer to the current Hilbert space
  // hamiltonian = pointer to the current Hamiltonian
  // lanczos = pointer to the Lanczos algorithm manager
  // subspaceStr = string to prepend in output file for each eigenvalue in this subspace
  // subspaceLegend = legend indicating contents of subspaceStr to include in output file
  // shift = energy shift that is applied to the hamiltonian
  // outputFileName = name of the file where results have to be stored
  // firstRun = flag that indicates if it the first time the main task is used
  // eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector
  GenericComplexMainTask(OptionManager* options, AbstractHilbertSpace* space, LanczosManager* lanczos, 
			 AbstractHamiltonian* hamiltonian, const char* subspaceStr, const char* subspaceLegend,
			 double shift, char* outputFileName, bool firstRun=true, char* eigenvectorFileName=0);
  
  // destructor
  //  
  ~GenericComplexMainTask();
  
  // set architecture binded to the task
  // 
  // architecture = pointer to the architecture to use
  void SetArchitecture(AbstractArchitecture* architecture);

  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  int ExecuteMainTask();

  // add optiongroup with options related to this module to the given OptionManager
  //
  void AddOptionGroup(OptionManager *optionManager);

 protected:

  // write a line of output to the results file
  //
  // file = stream to write to
  // value = numerical value to be printed after columns for flux and momentum (if defined)
  // terminate = indicate if line should be terminated with endl
  virtual void WriteResult(ofstream& file, double value, bool terminate=true);

  // do the Hamiltonian diagonalization in a given Hilbert subspace
  //
  // subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
  // file = reference on the output file stream where eigenvalues have to be stored
  virtual void DiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file);
    
};

#endif
