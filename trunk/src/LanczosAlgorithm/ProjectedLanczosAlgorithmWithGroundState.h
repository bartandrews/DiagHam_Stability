////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of basic Lanczos algorithm with real vectors           //
//                         and ground state evaluation                        //
//                      (without any re-orthogonalization)                    //
//                                                                            //
//                        last modification : 17/09/2002                      //
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


#ifndef PROJECTEDLANCZOSALGORITHMWITHGROUNDSTATE_H
#define PROJECTEDLANCZOSALGORITHMWITHGROUNDSTATE_H


#include "config.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "GeneralTools/GarbageFlag.h"


class ProjectedLanczosAlgorithmWithGroundState : public AbstractLanczosAlgorithm
{

 protected:

  RealVector V1;
  RealVector V2;
  RealVector V3;

  // vector that contains the initial state
  RealVector InitialState;

  //vector that contains the ground state (if GroundStateFlag is set to true)
  RealVector GroundState;

  // flag to indicate if ground state has already been computed
  bool GroundStateFlag;

  // flag to indicate the use of disk storage to increase speed of ground state calculation
  bool DiskFlag;
  // flag to indicate that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  bool ResumeDiskFlag;

  int Index;

  // number of wanted eigenvalues
  int NbrEigenvalue;
  // value of the last wanted eigenvalue at previous Lanczos iteration
  double PreviousLastWantedEigenvalue;

  // number of vectors to which eigenstates have to be othogonal
  int OrthogonalizationSetSize;
  // vectors to which eigenstates have to be othogonal
  RealVector* OrthogonalizationSet;
  // names of the files corresponding to vectors to which eigenstates have to be othogonal
  char** OrthogonalizationSetFileNames;

  // projector Hamiltonian(s)
  AbstractHamiltonian** Projectors;
  // number thereof
  int NbrProjectors;

  // flag for restarting
  bool RestartProjection;
    

  // Index for number of iterations for internal projector operations
  int InternalIndex;

  // value of the last estimate of the projector groundstate at the previous internal Lanczos iteration
  double PreviousProjectorGroundstate;
  double ProjectorGroundstate;
  RealTriDiagonalSymmetricMatrix InternalTridiagonalizedMatrix;
  RealTriDiagonalSymmetricMatrix InternalDiagonalizedMatrix;
  // precision needed for eigenvalues of projectors
  double ProjectorPrecision;

  // status flags for each of the lanczos vectors, and size thereof
  int *MainLanczosVectorFlags;
  int MainIterMax;
  int *ProjectorLanczosVectorFlags;
  int ProjectorIterMax;

  enum VectorFlags
    {
      VectorIndexMask = 0xffff,
      SavedOnDisk = 0x10000,
      SavedInMemory = 0x20000,
      SavedOnSwap1 = 0x40000,
      SavedOnSwap2 = 0x80000
    };

  // additional vectors for storage
  RealVector *LanczosVectorStorage;

  // number of storage vectors
  int NbrStorageVectors;
  
  // status flags for LanczosVectorStorage
  int *VectorStorageFlags;

  // count number of stored vectors of each sort
  int NbrStoredMain;
  int NbrStoredInternal;

  // Vectors stored in swap files
  int Swap1Index;
  int Swap2Index;
  
  enum StorageFlagDescriptions
    {
      Empty = 0x0,
      StorageIndexMask = 0xffff,
      MainLanczosVector = 0x10000,
      ProjectorLanczosVector = 0x20000
    };

  // temporary array for filenames
  char *TmpOutputName;

  // garbage flag
  GarbageFlag Flag;
  

 public:

  // default constructor
  //
  // projectors = operators to use for projection after each application of the Hamiltonian
  // nbrProjectors = number of separate projectors
  // architecture = architecture to use for matrix operations
  // maxIter = an approximation of maximal number of iteration
  // nbrStorageVectors = number of vectors that can be held in memory (in addition to minimum of 3)
  // projectorIterMax = max number of iterations before restarting internal projector lanczos
  // diskFlag = use disk storage to increase speed of ground state calculation
  // resumeDiskFlag = indicates that the Lanczos algorithm has to be resumed from an unfinished one (loading initial Lanczos algorithm state from disk)
  // projectorPrecision = precision required for projector operations
  // restartProjection = flag indicating whether projection should be restarted in precision not reached
  ProjectedLanczosAlgorithmWithGroundState(AbstractHamiltonian** projectors, int nbrProjectors, AbstractArchitecture* architecture, int maxIter = 0, int nbrStorageVectors = 0, int projectorIterMax = 100, bool diskFlag = false, bool resumeDiskFlag = false, double projectorPrecision = 1e-13, bool restartProjection = false);

  // copy constructor
  //
  // algorithm = algorithm from which new one will be created
  ProjectedLanczosAlgorithmWithGroundState(const ProjectedLanczosAlgorithmWithGroundState& algorithm);

  // destructor
  //
  ~ProjectedLanczosAlgorithmWithGroundState();

  // initialize Lanczos algorithm with a random vector
  //
  void InitializeLanczosAlgorithm();
  
  // initialize Lanczos algorithm with a given vector
  //
  // vector = reference to the vector used as first step vector
  void InitializeLanczosAlgorithm(const Vector& vector);

  // force orthogonalization with respect to a set of vectors
  //
  // fileName = name of the file describing the set of vectors
  // return value = true if no error occured
  bool ForceOrthogonalization(char* fileName);

  // get the n first eigenstates (limited to the ground state fro this class, return NULL if nbrEigenstates > 1)
  //
  // nbrEigenstates = number of needed eigenstates
  // return value = array containing the eigenstates
  Vector* GetEigenstates(int nbrEigenstates);

  // get ground state (by re-running Lanczos algorithm)
  //
  // return value = reference on ground state
  Vector& GetGroundState();

  // run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrIter = number of iteration to do 
  void RunLanczosAlgorithm (int nbrIter);
  
  // test if convergence has been reached
  //
  // return value = true if convergence has been reached
  bool TestConvergence ();
  
 protected:
  
  // read current Lanczos state from disk
  //
  // return value = true if no error occurs
  bool ReadState();

  // write current Lanczos state on disk
  //
  // return value = true if no error occurs
  bool WriteState();

  // orthogonalize a vector with respect to a set of external vectors
  //
  // inputVector = reference on the vector whose component on the external set has to be removed 
  void ExternalOrthogonalization(RealVector& inputVector);


  // project the vector stored in V1 to the groundstate of the given Projector number
  // nbrProjector = number of projector to use for this run
  // return = true if Lanczos was run, or false if no projection needed to be done
  bool ProjectVector(int nbrProjector);
  

  // run current internal Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
  //
  // nbrProjector = number of projector to use for this run
  // nbrIter = number of iteration to do 
  // return = true, if initial vector satisfied the condition that projector was zero
  bool RunProjectorLanczosAlgorithm (int nbrProjector, int nbrIter);


  // initialize Lanczos algorithm starting from the vector currently held in V1
  //
  void InitializeProjectorLanczosAlgorithm();

  // clear storage associated with internal Lanczos algorithm
  //
  void ClearProjectorLanczosAlgorithm();

  // diagonalize internal tridiagonalized matrix and find ground state energy
  //
  void ProjectorDiagonalize ();

  // test if convergence has been reached in internal lanczos
  //
  // return value = true if convergence has been reached
  bool TestProjectorConvergence ();

  // get last vector produced by projector lanczos; result is stored in internal vector V1
  //
  // return value = reference on produced vector
  Vector& GetProjectorGroundState();


  // save vector, either to internal memory, or disk, or both
  // vec = vector to be saved
  // index = vector index in Lanczos routine
  // mainLanczos = flag indicating whether vector is part of main lanczos algorithm
  // keepOriginal = flag indicating whether original vector needs to be kept in place
  // return = true on success
  bool SaveVector(RealVector &vec, int index, bool mainLanczos = true, bool keepOriginal = false); // test version
  bool SaveVectorTest(RealVector &vec, int index, bool mainLanczos = true, bool keepOriginal = false);
  // reread vector
  // vec = vector to be retrieved
  // index = vector index in Lanczos routine
  // mainLanczos = flag indicating whether vector is part of main lanczos algorithm
  // keepCopy = flag indicating whether the saved vector still needs to be kept in (live) memory after reloading
  void ReadVector(RealVector &vec, int index, bool mainLanczos = true, bool keepCopy = true);  // test version
  void ReadVectorTest(RealVector &vec, int index, bool mainLanczos = true, bool keepCopy = true);
};

#endif
