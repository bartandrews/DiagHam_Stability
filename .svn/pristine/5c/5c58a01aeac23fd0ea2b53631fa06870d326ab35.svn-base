////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                         in their quasihole sector                          //
//                                                                            //
//                        last modification : 11/02/2013                      //
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


#ifndef FQHEMPSCLUSTERED2RQUASIHOLESECTORMATRIX_H
#define FQHEMPSCLUSTERED2RQUASIHOLESECTORMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;


class FQHEMPSClustered2RQuasiholeSectorMatrix : public FQHEMPSClustered2RMatrix
{

 protected:

  // conformal weights of the sigma field
  LongRational WeightSigma;
  // conformal weights of the phi field (if r is even, then the sigma field has to be identical to the phi field)
  LongRational WeightPhi;

  // flag that indicates that the field is self conjugate
  bool SelfDualFlag;

  // matrix element of <phi|Psi(1)|sigma>
//  double MatrixElementNormalization;
  // conformal weight of the primary field used in the matrix element
//  LongRational WeightPrimaryFieldMatrixElement;

 public:
  
  // default constructor 
  //
  FQHEMPSClustered2RQuasiholeSectorMatrix();

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
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool bosonicVersion = false, bool useRational = true, 
					  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion = false, bool useRational = true, 
					  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RQuasiholeSectorMatrix(int pLevel, int nbrBMatrices, char* fileName, bool bosonicVersion = false, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0);

  // destructor
  //
  ~FQHEMPSClustered2RQuasiholeSectorMatrix();
  
  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture);

  // get the (k,r) exclude principle satisfied by the root configuration
  // 
  // pauliK = maximum number of particles in the pauliR consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual void GetKRExclusionPrinciple(int& pauliK, int& pauliR); 

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding);

  // compute the charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ);

  // get the matrix that into account the Jordan Wigner string on the torus geometry
  //
  // nbrFermions = number of fermions in the system
  // return value = corresponding matrix
  virtual SparseRealMatrix GetTorusStringMatrix(int nbrFermions);

  // get the auxiliary space indices that are related to a given topological scetor
  //
  // topologicalSector = index of the topological sector to select
  // nbrIndices = reference on the integer that will be set to the number of indices
  // return value = array that contains the auxiliary space indices related to the selected topological sector
  virtual int* GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices);
  
  // get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
  // 
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // statistics = true if we are dealing with fermions
  // return value = minimum ky momentum 
  virtual int GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics);
  
 protected:

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  // nbrCharges = total number of charge indices
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
  // globalIndexShift = index of the first state at the considered level
  // return value = linearized index
  virtual int Get2RMatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, 
			       int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);
  
  // compute the linearized index of a block of the matrices with fixed charge and p-level values
  //
  // chargedPartitionIndex = index of the partition in the charge sector
  // nbrCharges = total number of partitions at the current level
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
  // globalIndexShift = index of the first state at the considered level 
  // return value = linearized index
  virtual int Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
				      int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);

};

// compute the linearized index of the B matrix for the (k=2,r) clustered states
//
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector
// nbrCharges = total number of charge indices
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// globalIndexShift = index of the first state at the considered level
// return value = linearized index

inline int FQHEMPSClustered2RQuasiholeSectorMatrix::Get2RMatrixIndex(int charge, int chargedPartitionIndex, 
								     int nbrCharges, int chargeSectorDimension, 
								     int fieldIndex, int neutralPartitionIndex, 
								     int nbrIdentityDescendant, int globalIndexShift)
{
  return (((((nbrIdentityDescendant * fieldIndex) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}

// compute the linearized index of a block of the matrices with fixed charge and p-level values
//
// chargedPartitionIndex = index of the partition in the charge sector
// nbrCharges = total number of partitions at the current level
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// globalIndexShift = index of the first state at the considered level 
// return value = linearized index

inline int FQHEMPSClustered2RQuasiholeSectorMatrix::Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
									    int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift)
{
  return (globalIndexShift + (chargedPartitionIndex + chargeSectorDimension * ((fieldIndex * nbrIdentityDescendant) + neutralPartitionIndex)));
}


#endif
