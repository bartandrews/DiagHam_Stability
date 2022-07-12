////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract MPS matrix for the FQHE                 //
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



#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/Endian.h"
#include "MathTools/BinomialCoefficients.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

AbstractFQHEMPSMatrix::AbstractFQHEMPSMatrix()
{
  this->NbrBMatrices = 0;
  this->RealBMatrices = 0;
  this->ComplexBMatrices = 0;
  this->QuasiholeBMatrices = 0;
  this->PhysicalIndices = new unsigned long[2];
  this->PhysicalIndices[0] = 0x0ul;
  this->PhysicalIndices[1] = 0x1ul;
  this->TorusFlag = false;
}

// destructor
//

AbstractFQHEMPSMatrix::~AbstractFQHEMPSMatrix()
{
  if (this->RealBMatrices != 0)
    delete[] this->RealBMatrices;
  if (this->ComplexBMatrices != 0)
    delete[] this->ComplexBMatrices;
  if (this->QuasiholeBMatrices != 0)
    delete[] this->QuasiholeBMatrices;
}
  
// save the matrices 
// 
// fileName = name of the file where the matrices have to be stored
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::SaveMatrices (const char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrBMatrices);
  this->SaveHeader(File);
  if (this->RealBMatrices != 0)
    {
      int ComplexFlag  = 0;
      WriteLittleEndian(File, ComplexFlag);
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (this->RealBMatrices[i].WriteMatrix(File) == false)
	    {
	      File.close();  
	      cout << "error while storing matrix " << i << ", can't create " << fileName << endl;
	      return false;
	    }
	}
    }
  File.close();  
  return true;
}

// load the matrices 
// 
// fileName = name of the file where the matrices are stored
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::LoadMatrices (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrBMatrices);
  this->LoadHeader(File);
  int ComplexFlag  = 0;
  ReadLittleEndian(File, ComplexFlag);
  if (ComplexFlag == 0)
    {
      this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
      for (int i = 0; i < this->NbrBMatrices; ++i)
	{
	  if (this->RealBMatrices[i].ReadMatrix(File) == false)
	    {
	      File.close();  
	      cout << "error while reading matrix " << i << ", can't read " << fileName << endl;
	      return false;
	    }
	}
    }
  File.close();  
  return true;
}

// get the name describing the B matrices 
// 
// return value = name 

char* AbstractFQHEMPSMatrix::GetName()
{
  char* TmpName = new char[16];
  sprintf(TmpName, "dummy");
  return TmpName;
}

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

void AbstractFQHEMPSMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  numerator = 0;
  denominator = 0;
}

// get the (k,r) exclude principle satisfied by the root configuration
// 
// pauliK = maximum number of particles in the pauliR consecutive orbitals
// pauliR = number of consecutive orbitals

void AbstractFQHEMPSMatrix::GetKRExclusionPrinciple(int& pauliK, int& pauliR)
{
  this->GetFillingFactor(pauliK, pauliR);
}

// get the minimum ky momentum (i.e. within the reduced Brillouin zone) on the torus compatible with the current state
// 
// nbrParticles = number of particles
// nbrFluxQuanta = number of flux quanta
// statistics = true if we are dealing with fermions
// return value = minimum ky momentum 

int AbstractFQHEMPSMatrix::GetTorusMinimumKyMomentum(int nbrParticles, int nbrFluxQuanta, bool statistics)
{
  return 0;
}
  
// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = tuncation level of the block left indices
// q1 = charge index of the block left indices
// pLevel1 = tuncation level of the block right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

SparseRealMatrix AbstractFQHEMPSMatrix::ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2)
{
  int BlockNbrRow = this->GetBondIndexRange(pLevel1, q1);
  int BlockNbrColumn = this->GetBondIndexRange(pLevel2, q2);
  SparseRealMatrix TmpMatrix(BlockNbrRow, BlockNbrColumn);
  double Tmp = 0.0;
  for (int i = 0; i < BlockNbrRow; ++i)
    {
      for (int j = 0; j < BlockNbrColumn; ++j)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargeAndPLevel(i, pLevel1, q1),
				  this->GetBondIndexWithFixedChargeAndPLevel(j, pLevel2, q2), Tmp);
	  if (Tmp != 0.0)
	    {
	      TmpMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
    }
  return TmpMatrix;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = tuncation level of the block left indices
// q1 = charge index of the block left indices
// pLevel1 = tuncation level of the block right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

RealMatrix AbstractFQHEMPSMatrix::ExtractBlock(RealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2)
{
  int BlockNbrRow = this->GetBondIndexRange(pLevel1, q1);
  int BlockNbrColumn = this->GetBondIndexRange(pLevel2, q2);
  RealMatrix TmpMatrix(BlockNbrRow, BlockNbrColumn, true);
  double Tmp = 0.0;
  for (int i = 0; i < BlockNbrRow; ++i)
    {
      for (int j = 0; j < BlockNbrColumn; ++j)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargeAndPLevel(i, pLevel1, q1),
				  this->GetBondIndexWithFixedChargeAndPLevel(j, pLevel2, q2), Tmp);
	  if (Tmp != 0.0)
	    {
	      TmpMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
    }
  return TmpMatrix;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = truncation level of the block left indices
// cftSector1 = CFT sector of the blck left indices
// q1 = charge index of the block left indices
// pLevel1 = truncation level of the block right indices
// cftSector2 = CFT sector of the blck right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

SparseRealMatrix AbstractFQHEMPSMatrix::ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int cftSector1, int q1, int pLevel2, int cftSector2, int q2)
{
  int BlockNbrRow = this->GetBondIndexRange(pLevel1, q1, cftSector1);
  int BlockNbrColumn = this->GetBondIndexRange(pLevel2, q2, cftSector2);
  SparseRealMatrix TmpMatrix(BlockNbrRow, BlockNbrColumn);
  double Tmp = 0.0;
  for (int i = 0; i < BlockNbrRow; ++i)
    {
      for (int j = 0; j < BlockNbrColumn; ++j)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargePLevelCFTSector(i, pLevel1, q1, cftSector1),
				  this->GetBondIndexWithFixedChargePLevelCFTSector(j, pLevel2, q2, cftSector2), Tmp);
	  if (Tmp != 0.0)
	    {
	      TmpMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
    }
  return TmpMatrix;
}

// extract a block with fixed quantum numbers of a given matrix written the MPS basis
//
// matrix = reference on the matrix
// pLevel1 = truncation level of the block left indices
// cftSector1 = CFT sector of the blck left indices
// q1 = charge index of the block left indices
// pLevel1 = truncation level of the block right indices
// cftSector2 = CFT sector of the blck right indices
// q2 = charge index of the block left indices
// return value = block corresponding to the quantum numbers

RealMatrix AbstractFQHEMPSMatrix::ExtractBlock(RealMatrix& matrix, int pLevel1, int cftSector1, int q1, int pLevel2, int cftSector2, int q2)
{
  int BlockNbrRow = this->GetBondIndexRange(pLevel1, q1, cftSector1);
  int BlockNbrColumn = this->GetBondIndexRange(pLevel2, q2, cftSector2);
  RealMatrix TmpMatrix(BlockNbrRow, BlockNbrColumn, true);
  double Tmp = 0.0;
  for (int i = 0; i < BlockNbrRow; ++i)
    {
      for (int j = 0; j < BlockNbrColumn; ++j)
	{
	  matrix.GetMatrixElement(this->GetBondIndexWithFixedChargePLevelCFTSector(i, pLevel1, q1, cftSector1),
				  this->GetBondIndexWithFixedChargePLevelCFTSector(j, pLevel2, q2, cftSector2), Tmp);
	  if (Tmp != 0.0)
	    {
	      TmpMatrix.SetMatrixElement(i, j, Tmp);
	    }
	}
    }
  return TmpMatrix;
}

// get the charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::GetChargeIndexRange (int pLevel, int& minQ, int& maxQ)
{
  minQ = 1;
  maxQ = 0;
  return;
}

// get the charge index range at a given truncation level and in a given CFT sector
// 
// pLevel = tuncation level
// cftSector = CFT sector
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ)
{
  return this->GetChargeIndexRange(pLevel, minQ, maxQ);
}

// compute the global charge index range at a given truncation level
// 
// pLevel = tuncation level
// minQ = reference on the lowest charge index
// maxQ = reference on the lowest charge index

void AbstractFQHEMPSMatrix::ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ)
{
  minQ = 1;
  maxQ = 0;
  return;
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::LoadHeader (ifstream& file)
{
  int HeaderSize = 0;
  ReadLittleEndian(file, HeaderSize);
  file.seekg (HeaderSize, ios::cur);
  return true;
}

// save the specific informations to the file header 
// 
// file = reference on the output file stream
// return value = true if no error occurred  

bool AbstractFQHEMPSMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = 0;
  WriteLittleEndian(file, HeaderSize);
  return true;
}

// get the edge matrix for localized quasiholes, with normal ordering
//
// nbrQuasiholes = number of quasiholes
// quasiholePositions = quasihole positions (for cylinder, positions have to be expressed in perimeter units)
// return value = pointer to the edge matrix

SparseComplexMatrix* AbstractFQHEMPSMatrix::GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions)
{
  return 0;
}
  
// get the matrix that into account the Jordan Wigner string on the torus geometry
//
// nbrFermions = number of fermions in the system
// return value = corresponding matrix

SparseRealMatrix AbstractFQHEMPSMatrix::GetTorusStringMatrix(int nbrFermions)
{
  int TmpDimension = 0;
  if (this->RealBMatrices != 0)
    {
      TmpDimension = this->RealBMatrices[0].GetNbrColumn();
    }
  else
    {
      TmpDimension = this->ComplexBMatrices[0].GetNbrColumn();
    }
  int* TmpNbrElementPerRow =  new int [TmpDimension];
  for (int i = 0; i < TmpDimension; ++i)
    {
      TmpNbrElementPerRow[i] = 1;
    }
  SparseRealMatrix StringMatrix (TmpDimension, this->RealBMatrices[0].GetNbrColumn(), TmpNbrElementPerRow);
  for (int i = 0; i < StringMatrix.GetNbrColumn(); ++i)
    {
      StringMatrix.SetMatrixElement(i, i, 1.0);
    } 
  return StringMatrix;
}

// get the range for the bond index when fixing the tuncation level and the charge index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = range for the bond index with fixed tuncation level and charge index

int AbstractFQHEMPSMatrix::GetBondIndexRange(int pLevel, int qValue)
{
  return 0;
}

// get the range for the bond index when fixing the tuncation level, charge and CFT sector index
//
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = range for the bond index with fixed tuncation level, charge and CFT sector index

int AbstractFQHEMPSMatrix::GetBondIndexRange(int pLevel, int qValue, int cftSector)
{
  return this->GetBondIndexRange(pLevel, qValue);
}

// get the bond index for a fixed truncation level and the charge index 
//
// localIndex = bond index in the pLevel and qValue restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// return value = bond index in the full bond index range

int AbstractFQHEMPSMatrix::GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue)
{
  return 0;
}

// get the bond index for a fixed truncation level, charge and CFT sector index
//
// localIndex = bond index in the pLevel and qValue and cftSector restricted range
// pLevel = tuncation level of the block
// qValue = charge index of the block
// cftSector = CFT sector index of the block
// return value = bond index in the full bond index range

int AbstractFQHEMPSMatrix::GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector)
{
  return this->GetBondIndexWithFixedChargeAndPLevel(localIndex, pLevel, qValue);
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void AbstractFQHEMPSMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
  rowIndex = -1;
  columnIndex = -1;
}

// get the number of particles that fit the root configuration once the number of flux quanta is fixed
// 
// nbrFluxQuanta = number of flux quanta
// padding = assume that the state has the extra padding
// return value = number of partciles

int AbstractFQHEMPSMatrix::GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding)
{
  int Numerator;
  int Denominator;
  this->GetFillingFactor(Numerator, Denominator);
  if (this->TorusFlag == false)
    {
      int NbrParticles = ((nbrFluxQuanta + 1) * Numerator);
      if ((NbrParticles % Denominator) == 0)
	return (NbrParticles / Denominator);
      else
	return ((NbrParticles / Denominator) + 1);
    }
  else
    {
      return ((nbrFluxQuanta * Numerator) / Denominator);
    }
}

// get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
//
// cftSector = index of the CFT sector
// return value = Q sector shift

int AbstractFQHEMPSMatrix::GetQValueCFTSectorShift(int cftSector)
{
  return 0;
}

// get the auxiliary space indices that are related to a given topological scetor
//
// topologicalSector = index of the topological sector to select
// nbrIndices = reference on the integer that will be set to the number of indices
// return value = array that contains the auxiliary space indices related to the selected topological sector

int* AbstractFQHEMPSMatrix::GetTopologicalSectorIndices(int topologicalSector, int& nbrIndices)
{
  nbrIndices = 0;
  if (this->RealBMatrices != 0)
    {
      nbrIndices = this->RealBMatrices[0].GetNbrColumn();
    }
  else
    {
      nbrIndices = this->ComplexBMatrices[0].GetNbrColumn();
    }
  int* TmpIndices =  new int [nbrIndices];
  for (int i = 0; i < nbrIndices; ++i)
    {
      TmpIndices[i] = i;
    }
  return TmpIndices;
}
  
// print a given state of the auxiliary space
//
// str = reference on the output stream
// index = index of the state
// return value = reference on the output stream

ostream& AbstractFQHEMPSMatrix::PrintAuxiliarySpaceState(ostream& str, int index)
{
  str << "|" << index << ">";
  return str;
}

// get the parent MPS matrices if the current MPS matrices have ones
//
// return value = pointer to the parent MPS matrices

AbstractFQHEMPSMatrix* AbstractFQHEMPSMatrix::GetParentMPSMatrices()
{
  return 0;
}

// get the array that gives the index of each entry of the current B matrix within its parent B matrix
//
// return value = index mapping array

int* AbstractFQHEMPSMatrix::GetIndexMappingArray()
{
  return 0;
}

// get a given physical indiex
//
// index = index to retrieve
//  configuration = array where the description of the physical index will be stored

void AbstractFQHEMPSMatrix::GetPhysicalIndex(int index, unsigned long* configuration)
{  
  for (int i = 0; i < this->GetNbrOrbitals(); ++i)
    {
      configuration[i] = (this->PhysicalIndices[index] >> i) & 0x1ul;
    }
}

// get a given physical index at a given orbital for the site-dependent MPS
//
// orbitalIndex = orbital index 
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void AbstractFQHEMPSMatrix::GetSiteDependentPhysicalIndex(int orbitalIndex, int index, unsigned long* configuration)
{
  this->GetPhysicalIndex(index, configuration);
}


// print a given physical index
//
// str = reference on the output stream
// index = integer associated to the  physical index 
// return value = reference on the output stream

ostream& AbstractFQHEMPSMatrix::PrintPhysicalIndex(ostream& str, int index)
{
  unsigned long* TmpConfiguration = new unsigned long[this->GetNbrOrbitals()];
  this->GetPhysicalIndex(index, TmpConfiguration);
  str << TmpConfiguration[0];
  for (int i = 1; i < this->GetNbrOrbitals(); ++i)
    {
      str << " " << TmpConfiguration[i];
    }
  delete[] TmpConfiguration;
  return str;
}

// get the array where the site-dependent matrices are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseRealMatrix** AbstractFQHEMPSMatrix::GetSiteDependentMatrices(int nbrFluxQuanta)
{
  cout << "site dependent real MPS is not available for " << this->GetName() << endl;
  return 0;
}

// get the array where the site-dependent matrices are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseComplexMatrix** AbstractFQHEMPSMatrix::GetSiteDependentComplexMatrices(int nbrFluxQuanta)
{
  cout << "site dependent complex MPS is not available for " << this->GetName() << endl;
  return 0;
}

// get the array where the site-dependent matrices for the geometry are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseRealMatrix** AbstractFQHEMPSMatrix::GetSphereSiteDependentMatrices(int nbrFluxQuanta)
{
  BinomialCoefficients TmpCoef (nbrFluxQuanta);
  SparseRealMatrix** TmpMatrices = new SparseRealMatrix*[nbrFluxQuanta + 1];
  SparseRealMatrix* TmpSiteIndependentMatrices = this->GetMatrices();
  for (int i = 0; i <= nbrFluxQuanta; ++i)
    {
      TmpMatrices[i] = new  SparseRealMatrix[this->NbrBMatrices];
      double TmpFactorialCoefficent = 1.0;
      TmpMatrices[i][0].Copy(TmpSiteIndependentMatrices[0]);
      for (int j = 1; j < this->NbrBMatrices; ++j)
	{
	  TmpMatrices[i][j].Copy(TmpSiteIndependentMatrices[j]);
	  TmpFactorialCoefficent /= sqrt((double) j);
	  TmpMatrices[i][j] *= (pow(TmpCoef.GetNumericalCoefficient(nbrFluxQuanta, i), -0.5 * ((double) j))
				* TmpFactorialCoefficent);
	}
    }
  return TmpMatrices;
}

// get the array where the site-dependent matrices for the geometry are stored
//
// nbrFluxQuanta = number of flux quanta in the finite size system
// return value = pointer to the array of matrices (first entry being the orbital index, the second being the occupation number)

SparseComplexMatrix** AbstractFQHEMPSMatrix::GetSphereSiteDependentComplexMatrices(int nbrFluxQuanta)
{
  BinomialCoefficients TmpCoef (nbrFluxQuanta);
  SparseComplexMatrix** TmpMatrices = new SparseComplexMatrix*[nbrFluxQuanta + 1];
  SparseComplexMatrix* TmpSiteIndependentMatrices = this->GetComplexMatrices();
  for (int i = 0; i <= nbrFluxQuanta; ++i)
    {
      TmpMatrices[i] = new  SparseComplexMatrix[this->NbrBMatrices];
      double TmpFactorialCoefficent = 1.0;
      TmpMatrices[i][0].Copy(TmpSiteIndependentMatrices[0]);
      for (int j = 1; i < this->NbrBMatrices; ++j)
	{
	  TmpMatrices[i][j].Copy(TmpSiteIndependentMatrices[j]);
	  TmpFactorialCoefficent /= sqrt((double) j);
	  TmpMatrices[i][j] *= (pow(TmpCoef.GetNumericalCoefficient(nbrFluxQuanta, i), 0.5 * ((double) j))
				* TmpFactorialCoefficent);
	}
    }
  return TmpMatrices;
}


// get the label associated to a given state of the auxiliary space
//
// index = auxiliary space index
// return value = string containing the label

char* AbstractFQHEMPSMatrix::GetAuxiliarySpaceLabel(int index)
{
  char* TmpString = new char [128];
  int TmpQ = 0;
  int TmpPLevel = 0;
  int TmpCFTSector = 0;
  if (this->GetNbrCFTSectors() == 1)
    {
      this->GetChargeAndPLevelFromMatrixIndex(index, TmpPLevel, TmpQ);
      sprintf (TmpString, "(Q=%d, P=%d, i=%d)", TmpQ, TmpPLevel, index);
    }
  else
    {
      this->GetCFTSectorChargeAndPLevelFromMatrixIndex(index, TmpCFTSector, TmpPLevel, TmpQ);
      sprintf (TmpString, "(x=%d, Q=%d, P=%d, i=%d)", TmpCFTSector, TmpQ, TmpPLevel, index);
    }
  return TmpString;
}

// compute the site-dependent matrices
//
// initialOrbitalIndex = index of the first orbital
// lastOrbitalIndex = index of the last orbital

void AbstractFQHEMPSMatrix::ComputeSiteDependentMatrices(int initialOrbitalIndex, int lastOrbitalIndex)
{
  cout << "ComputeSiteDependentMatrices is not available for " << this->GetName() << endl;
}


// get the site-dependent matrices (real version) computed through ComputeSiteDependentMatrices
//
// siteDependentMatrices = reference on the site-dependent matrices
// nbrSiteDependentMatrices = reference on the array providing the number of site-dependent matrices per orbital
// siteDependentMatrixOrbitalIndices = reference on the array providing the orbital indices 
// siteDependentPhysicalIndices = reference on the array providing the physical indices associated to each site-dependent matrix
// return value = number of orbitals covered by the site-dependent matrices

int AbstractFQHEMPSMatrix::GetSiteDependentMatrices(SparseRealMatrix**& siteDependentMatrices, int*& nbrSiteDependentMatrices, int*& siteDependentMatrixOrbitalIndices, unsigned long**& siteDependentPhysicalIndices)
{
  cout << "GetSiteDependentMatrices is not available for " << this->GetName() << endl;
  return 0;
}
