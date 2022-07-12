////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract MPO matrix for the FQHE                //
//                                                                            //
//                        last modification : 29/07/2016                      //
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


#ifndef ABSTRACTFQHEMPOMATRIX_H
#define ABSTRACTFQHEMPOMATRIX_H


#include "config.h"
#include "MathTools/Complex.h" 

#include <fstream>

class SparseRealMatrix;
class SparseComplexMatrix;


using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

class AbstractFQHEMPOMatrix
{

 protected:

  // total number of O matrices
  int NbrOMatrices;

  // array that describes all the physical indices for the upper bond
  unsigned long* UpperPhysicalIndices;
  // array that describes all the physical indices for the lower bond
  unsigned long* LowerPhysicalIndices;

  // array where the O matrices are stored (for real matrices)
  SparseRealMatrix* RealOMatrices;

  // array where the B matrices are stored (for complex matrices)
  SparseComplexMatrix* ComplexOMatrices;

  // MPO bond dimension
  int BondDimension;

 public:
  
  // default constructor 
  //
  AbstractFQHEMPOMatrix();

  // destructor
  //
  ~AbstractFQHEMPOMatrix();
  
  // save the matrices 
  // 
  // fileName = name of the file where the matrices have to be stored
  // return value = true if no error occurred  
  virtual bool SaveMatrices (const char* fileName);

  // load the matrices 
  // 
  // fileName = name of the file where the matrices are stored
  // return value = true if no error occurred  
  virtual bool LoadMatrices (const char* fileName);

  // get the number of B matrices
  //
  // return value = number of B matrices
  virtual int GetNbrMatrices();

  // get the number of orbitals that associated to a set of B matrices
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the maximum occupation per orbital
  //
  // return value = aximum occupation per orbital
  virtual int GetMaximumOccupation();

  // get the MPO bond dimension
  //
  // return value = MPO bond dimension
  virtual int GetBondDimension();

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseRealMatrix* GetMatrices();

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseComplexMatrix* GetComplexMatrices();

  // get the name describing the O matrices 
  // 
  // return value = name 
  virtual char* GetName();

  // get the boundary indices of the MPO representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex);

  // get a given physical index for the upper bond
  //
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  void GetUpperPhysicalIndex(int index, unsigned long* configuration);

  // get a given physical index for the lower bond
  //
  // index = index to retrieve
  // configuration = array where the description of the physical index will be stored
  void GetLowerPhysicalIndex(int index, unsigned long* configuration);

  // get the array of physical indices for the upper bond
  //
  // return value  = array of physical indices
  virtual unsigned long* GetUpperPhysicalIndices();

  // get the array of physical indices for the lower bond
  //
  // return value  = array of physical indices
  virtual unsigned long* GetLowerPhysicalIndices();

  // print a given physical index for the upper bond
  //
  // str = reference on the output stream
  // index = integer associated to the  physical index 
  // return value = reference on the output stream
  virtual ostream& PrintUpperPhysicalIndex(ostream& str, int index);

  // print a given physical index for the lower bond
  //
  // str = reference on the output stream
  // index = integer associated to the  physical index 
  // return value = reference on the output stream
  virtual ostream& PrintLowerPhysicalIndex(ostream& str, int index);

  // print a given state of the auxiliary space
  //
  // str = reference on the output stream
  // index = index of the state
  // return value = reference on the output stream
  virtual ostream& PrintAuxiliarySpaceState(ostream& str, int index);

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

};

// get the number of B matrices
//
// return value = number of B matrices

inline int AbstractFQHEMPOMatrix::GetNbrMatrices()
{
  return this->NbrOMatrices;
}

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseRealMatrix* AbstractFQHEMPOMatrix::GetMatrices()
{
  return this->RealOMatrices;
}

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseComplexMatrix* AbstractFQHEMPOMatrix::GetComplexMatrices()
{
  return this->ComplexOMatrices;
}

// get the array of physical indices for the upper bond
//
// return value  = array of physical indices

inline unsigned long* AbstractFQHEMPOMatrix::GetUpperPhysicalIndices()
{
  return this->UpperPhysicalIndices;
}

// get the array of physical indices for the lower bond
//
// return value  = array of physical indices

inline unsigned long* AbstractFQHEMPOMatrix::GetLowerPhysicalIndices()
{
  return this->LowerPhysicalIndices;
}

// get the MPO bond dimension
//
// return value = MPO bond dimension

inline int AbstractFQHEMPOMatrix::GetBondDimension()
{
  return this->BondDimension;
}

// get the number of orbitals that associated to a set of B matrices
//
// return value = number of orbitals

inline int AbstractFQHEMPOMatrix::GetNbrOrbitals()
{
  return 1;
}

// get the maximum occupation per orbital
//
// return value = aximum occupation per orbital

inline int AbstractFQHEMPOMatrix::GetMaximumOccupation()
{
  return 1;
}

#endif
