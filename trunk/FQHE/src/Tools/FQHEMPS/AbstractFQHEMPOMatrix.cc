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



#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPOMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

AbstractFQHEMPOMatrix::AbstractFQHEMPOMatrix()
{
  this->NbrOMatrices = 0;
  this->RealOMatrices = 0;
  this->ComplexOMatrices = 0;
  this->UpperPhysicalIndices = 0;
  this->LowerPhysicalIndices = 0;
}

// destructor
//

AbstractFQHEMPOMatrix::~AbstractFQHEMPOMatrix()
{
  if (this->RealOMatrices != 0)
    delete[] this->RealOMatrices;
  if (this->ComplexOMatrices != 0)
    delete[] this->ComplexOMatrices;
  if (this->UpperPhysicalIndices != 0)
    delete[] this->UpperPhysicalIndices;
  if (this->LowerPhysicalIndices != 0)
    delete[] this->LowerPhysicalIndices;
}
  
// save the matrices 
// 
// fileName = name of the file where the matrices have to be stored
// return value = true if no error occurred  

bool AbstractFQHEMPOMatrix::SaveMatrices (const char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrOMatrices);
  for (int i = 0; i < this->NbrOMatrices; ++i)
    {
      WriteLittleEndian(File, this->UpperPhysicalIndices[i]);
      WriteLittleEndian(File, this->LowerPhysicalIndices[i]);
    }
  this->SaveHeader(File);
  if (this->RealOMatrices != 0)
    {
      int ComplexFlag  = 0;
      WriteLittleEndian(File, ComplexFlag);
      for (int i = 0; i < this->NbrOMatrices; ++i)
	{
	  if (this->RealOMatrices[i].WriteMatrix(File) == false)
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

bool AbstractFQHEMPOMatrix::LoadMatrices (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrOMatrices);
  for (int i = 0; i < this->NbrOMatrices; ++i)
    {
      ReadLittleEndian(File, this->UpperPhysicalIndices[i]);
      ReadLittleEndian(File, this->LowerPhysicalIndices[i]);
    }

  this->LoadHeader(File);
  int ComplexFlag  = 0;
  ReadLittleEndian(File, ComplexFlag);
  if (ComplexFlag == 0)
    {
      this->RealOMatrices = new SparseRealMatrix [this->NbrOMatrices];
      for (int i = 0; i < this->NbrOMatrices; ++i)
	{
	  if (this->RealOMatrices[i].ReadMatrix(File) == false)
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

char* AbstractFQHEMPOMatrix::GetName()
{
  char* TmpName = new char[16];
  sprintf(TmpName, "dummy_mpo");
  return TmpName;
}

// load the specific informations from the file header
// 
// file = reference on the input file stream
// return value = true if no error occurred  

bool AbstractFQHEMPOMatrix::LoadHeader (ifstream& file)
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

bool AbstractFQHEMPOMatrix::SaveHeader (ofstream& file)
{
  int HeaderSize = 0;
  WriteLittleEndian(file, HeaderSize);
  return true;
}

// get the boundary indices of the MPO representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index

void AbstractFQHEMPOMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex)
{
  rowIndex = -1;
  columnIndex = -1;
}
  
// print a given state of the auxiliary space
//
// str = reference on the output stream
// index = index of the state
// return value = reference on the output stream

ostream& AbstractFQHEMPOMatrix::PrintAuxiliarySpaceState(ostream& str, int index)
{
  str << "|" << index << ">";
  return str;
}

// get a given physical index for the upper bond
//
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void AbstractFQHEMPOMatrix::GetUpperPhysicalIndex(int index, unsigned long* configuration)
{  
  for (int i = 0; i < this->GetNbrOrbitals(); ++i)
    {
      configuration[i] = (this->UpperPhysicalIndices[index] >> i) & 0x1ul;
    }
}

// get a given physical index for the lower bond
//
// index = index to retrieve
// configuration = array where the description of the physical index will be stored

void AbstractFQHEMPOMatrix::GetLowerPhysicalIndex(int index, unsigned long* configuration)
{  
  for (int i = 0; i < this->GetNbrOrbitals(); ++i)
    {
      configuration[i] = (this->LowerPhysicalIndices[index] >> i) & 0x1ul;
    }
}

// print a given physical index for the upper bond
//
// str = reference on the output stream
// index = integer associated to the  physical index 
// return value = reference on the output stream

ostream& AbstractFQHEMPOMatrix::PrintUpperPhysicalIndex(ostream& str, int index)
{
  unsigned long* TmpConfiguration = new unsigned long[this->GetNbrOrbitals()];
  this->GetUpperPhysicalIndex(index, TmpConfiguration);
  str << TmpConfiguration[0];
  for (int i = 1; i < this->GetNbrOrbitals(); ++i)
    {
      str << " " << TmpConfiguration[i];
    }
  delete[] TmpConfiguration;
}

// print a given physical index for the lower bond
//
// str = reference on the output stream
// index = integer associated to the  physical index 
// return value = reference on the output stream

ostream& AbstractFQHEMPOMatrix::PrintLowerPhysicalIndex(ostream& str, int index)
{
  unsigned long* TmpConfiguration = new unsigned long[this->GetNbrOrbitals()];
  this->GetLowerPhysicalIndex(index, TmpConfiguration);
  str << TmpConfiguration[0];
  for (int i = 1; i < this->GetNbrOrbitals(); ++i)
    {
      str << " " << TmpConfiguration[i];
    }
  delete[] TmpConfiguration;
}

