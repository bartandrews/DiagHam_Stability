////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of MPS matrix for the Laughlin state                //
//                                                                            //
//                        last modification : 30/10/2012                      //
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


#ifndef FQHEMPSLAUGHLINMATRIX_H
#define FQHEMPSLAUGHLINMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "MathTools/FactorialCoefficient.h"


class FQHEMPSLaughlinMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // use a version of the code that is compatible with bosonic wave functions
  bool BosonicVersion;

  // index of the Laughlin (i.e. 1/nu) 
  int LaughlinIndex;

  // |P| level truncation
  int PLevel;
  
  // number of charge indices
  int NbrNValue;
  // global shift of charge indices to make it non-negative
  int NValueGlobalShift;
  // first linearized index for each truncation level
  int* TotalStartingIndexPerPLevel;
  // number of linearized index per truncation level
  int* NbrIndicesPerPLevel;
  // number of N (i.e. charge) values per  truncation level
  int* NbrNValuesPerPLevel;
  // initial N (i.e. charge) value per  truncation level
  int* NInitialValuePerPLevel;
  // last N (i.e. charge) value per  truncation level
  int* NLastValuePerPLevel;
  // indicate that the charge index range does not depend on the truncation level
  bool UniformChargeIndexRange;

  // true if B_0 has to be normalized on the cylinder geometry
  bool CylinderFlag;
  // cylinder aspect ratio
  double Kappa;

  // inverse of the B[0] matrix, used to compute the site dependent MPS from the site independent MPS
  SparseRealMatrix InverseBMatrixZero;

  // number of orbitals covered by the site-dependent matrices
  int SiteDependentMatrixNbrOrbitals;
  // site dependent MPS
  SparseRealMatrix** SiteDependentMatrices;
  // signed orbital indices associated to each site-dependent matrix
  int* SiteDependentMatrixOrbitalIndices;
  //  number of site-dependent matrices for each orbital
  int* NbrSiteDependentMatrices;
  // physical indices for each site-dependent matrix
  unsigned long** SiteDependentPhysicalIndices;
  
  // true if B matrix have to be evaluated on the twisted torus geometry (requiring complex B matrices)
  bool TwistedTorusFlag;
  // number of flux quanta piercing the torus
  int TorusNbrFluxQuanta;
  // angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  double TorusAngle;
  // aspect ratio of the torus(norm of tau)
  double TorusAspectRatio;
  // factor in the exponential due to tau
  Complex TauFactor;
  // flux insertion along the tau direction
  double TorusFluxInsertion;

 public:
  
  // default constructor 
  //
  FQHEMPSLaughlinMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool bosonicVersion = false, 
			bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0);

  // constructor for the torus geometry
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // trimChargeIndices = trim the charge indices
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool bosonicVersion, bool trimChargeIndices, 
			int nbrFluxQuanta, double aspectRatio, double angle, double fluxInsertion);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSLaughlinMatrix();
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName ();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // get the MPS truncation level
  //
  // return value = truncation level
  virtual int GetTruncationLevel();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // create the B matrices for the laughlin state
  //
  virtual void CreateBMatrices ();

  // create the B matrices for the laughlin state
  //
  virtual void AlternateCreateBMatrices ();

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = tuncation level of the block left indices
  // q1 = charge index of the block left indices
  // pLevel1 = tuncation level of the block right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

  // get the range for the bond index when fixing the tuncation level and the charge index
  //
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = range for the bond index with fixed tuncation level and charge index
  virtual int GetBondIndexRange(int pLevel, int qValue);

  // get the bond index for a fixed truncation level and the charge index 
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue);

  // get the charge index range
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int& minQ, int& maxQ);

  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of particles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding);

  // get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
  // 
  // nbrParticles = number of particles
  // nbrFluxQuanta = number of flux quanta
  // statistics = true if we are dealing with fermions
  // return value = minimum ky momentum 
  virtual int GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics);
  
  // compute P, N from the linearized index of the B matrix for the Laughlin states
  //
  // index = linearized index
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  void GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex);

  // compute the level and the charge index of a given matrix index
  //
  // index = matrix index
  // pLevel = reference on the level
  // qValue = reference on the charge index
  virtual void GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // get the matrix that takes into account the Jordan Wigner string on the torus geometry
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
  
  // print a given state of the auxiliary space
  //
  // str = reference on the output stream
  // index = index of the state
  // return value = reference on the output stream
  virtual ostream& PrintAuxiliarySpaceState(ostream& str, int index);

  // get a given physical indiex
  //
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  virtual void GetPhysicalIndex(int index, unsigned long* configuration);

  // get the array where the site-dependent matrices are stored
  //
  // nbrFluxQuanta = number of flux quanta in the finite size system
  // return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)
  virtual SparseRealMatrix** GetSiteDependentMatrices(int nbrFluxQuanta);

  // compute the site-dependent matrices
  //
  // initialOrbitalIndex = index of the first orbital
  // lastOrbitalIndex = index of the last orbital
  virtual void ComputeSiteDependentMatrices(int initialOrbitalIndex, int lastOrbitalIndex);

  // get the site-dependent matrices (real version) computed through ComputeSiteDependentMatrices
  //
  // siteDependentMatrices = reference on the site-dependent matrices
  // nbrSiteDependentMatrices = reference on the array providing the number of site-dependent matrices per orbital
  // siteDependentMatrixOrbitalIndices = reference on the array providing the orbital indices 
  // siteDependentPhysicalIndices = reference on the array providing the physical indices associated to each site-dependent matrix
  // return value = number of orbitals covered by the site-dependent matrices
  virtual int GetSiteDependentMatrices(SparseRealMatrix**& siteDependentMatrices, int*& nbrSiteDependentMatrices, int*& siteDependentMatrixOrbitalIndices, unsigned long**& siteDependentPhysicalIndices);
  
 protected:

  // load the specific informations from the file header
  // 
  // file = reference on the input file stream
  // return value = true if no error occurred  
  virtual bool LoadHeader (ifstream& file);

  // save the specific informations to the file header 
  // 
  // file = reference on the output file stream
  // return value = true if no error occurred  
  virtual bool SaveHeader (ofstream& file);

  // compute the linearized index of the B matrix for the Laughlin states
  //
  // pLevel = truncation level
  // partitionIndex = index of the partition at the given pLevel
  // charge = charge index (with the global shift by Pmax)
  // return value = linearized index
  int GetMatrixIndex(int pLevel, int partitionIndex, int charge);

  // create the matrix element of the B matrix U(1) part
  //
  // chargeNumerator = numerator of the charge (in sqrt(q) unit)
  // chargeDenominator = denominator of the charge (in sqrt(q) unit)
  // partition1 = U(1) partition associated to the left state
  // p1Level = length of partition1
  // partition2 = U(1) partition associated to the left state
  // p1Level = length of partition2
  // coef = reference on a temporary factorial coefficient
  // return value = matrix element
  double CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, 
				       unsigned long* partition1, unsigned long* partition2, 
				       int p1Level, int p2Level, FactorialCoefficient& coef);

  // compute the charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeChargeIndexRange(int pLevel, int& minQ, int& maxQ);

  // compute the global charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ);

};

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSLaughlinMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return this->LaughlinIndex;
}

// get the maximum occupation per orbital
//
// return value = aximum occupation per orbital

inline int FQHEMPSLaughlinMatrix::GetMaximumOccupation()
{
  return (this->NbrBMatrices - 1);
}

// compute the linearized index of the B matrix for the Laughlin states
//
// pLevel = truncation level
// partitionIndex = index of the partition at the given pLevel
// charge = charge index (already with NValueGlobalShift)
// return value = linearized index
inline int FQHEMPSLaughlinMatrix::GetMatrixIndex(int pLevel, int partitionIndex, int charge)
{
    return this->TotalStartingIndexPerPLevel[pLevel] + this->NbrNValuesPerPLevel[pLevel] * partitionIndex + charge - this->NInitialValuePerPLevel[pLevel];
}

// compute P, N from the linearized index of the B matrix for the Laughlin states
//
// index = linearized index
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector

inline void FQHEMPSLaughlinMatrix::GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex)
{
   int i = 0;
   bool found = false;
   int U1dim;  
   int TmpSum = this->NbrIndicesPerPLevel[0];
   
   while ((i <= this->PLevel) && (found == false))
     {
        if (index < TmpSum)
           {
             found = true;
             chargedPartitionIndex = i;
             U1dim = this->NbrIndicesPerPLevel[i]/this->NbrNValue;
             charge = (this->NbrIndicesPerPLevel[i] - (TmpSum - index)) % this->NbrNValue;
           }
        else
          {
             i++;
             TmpSum += this->NbrIndicesPerPLevel[i];      
          }
     }
     
}

// compute the  level and the charge index of a given matrix index
//
// index = matrix index
// pLevel = reference on the level
// qValue = reference on the charge index

inline void FQHEMPSLaughlinMatrix::GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue)
{
  pLevel = 0;
  while ((pLevel <= this->PLevel) && (this->TotalStartingIndexPerPLevel[pLevel] <= index))
    ++pLevel;
  --pLevel;
  qValue = this->NInitialValuePerPLevel[pLevel] + (index - this->TotalStartingIndexPerPLevel[pLevel]) % this->NbrNValuesPerPLevel[pLevel];
}

// get the MPS truncation level
//
// return value = truncation level

inline int FQHEMPSLaughlinMatrix::GetTruncationLevel()
{
  return this->PLevel;
}

#endif
