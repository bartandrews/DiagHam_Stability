////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class for n dimensional complex vector                    //
//             where only part of the elements are actually stored            //
//                                                                            //
//                        last modification : 25/12/2007                      //
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


#include "Vector/PartialComplexVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <fstream>

#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;


// default constructor
//

PartialComplexVector::PartialComplexVector()
{
  this->VectorType = Vector::ComplexDatas | Vector::PartialData;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->Components = 0;
  this->VectorId = 0;
  this->IndexShift = 0;
  this->RealDimension = 0;
}

// constructor for an empty complex vector (all coordinates set to zero)
//
// size = effective vector dimension (i.e. number of stored elements)
// realSize = size of the full vector
// indexShift = index of the first element which is actually stored in the current partial vector
// zeroFlag = true if all coordinates have to be set to zero

PartialComplexVector::PartialComplexVector(int size, int realSize, int indexShift, bool zeroFlag)
{
  this->VectorType = Vector::ComplexDatas | Vector::PartialData;
  this->Dimension = size;
  this->IndexShift = indexShift;
  this->RealDimension = realSize;
  this->TrueDimension = this->Dimension;
  this->Components = new Complex [this->Dimension + 1]; 
  this->Flag.Initialize();
  this->VectorId = 0;
  if (zeroFlag == true)
    for (int i = 0; i < this->Dimension; i++)
      {
	this->Components[i] = 0.0;
      }
}

// constructor from an array of doubles
//
// real = array of doubles corresponding to real part
// imaginary = array of doubles corresponding to imaginary part
// size = effective vector dimension (i.e. number of stored elements)
// realSize = size of the full vector
// indexShift = index of the first element which is actually stored in the current partial vector
 
PartialComplexVector::PartialComplexVector(double* real, double* imaginary, int size, int realSize, int indexShift)
{
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
  this->VectorId = 0;
  this->IndexShift = indexShift;
  this->RealDimension = realSize;
  this->Components = new Complex [this->Dimension + 1];
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re = real[i];
      this->Components[i].Im = imaginary[i];
    }
  // delete [] real;
  // delete [] imaginary;
}

// copy constructor
//
// vector = vector to copy
// DuplicateFlag = true if datas have to be duplicated

PartialComplexVector::PartialComplexVector(const PartialComplexVector& vector, bool duplicateFlag)
{
  this->VectorType = Vector::ComplexDatas | Vector::PartialData;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
  if (duplicateFlag == false)
    {
      this->Components = vector.Components;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Flag.Initialize();
	  this->Components = new Complex [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    {
	      this->Components[i] = vector.Components[i];
	    }
	}
      else
	{
	  this->Components = 0;
	}
    }
}

#ifdef __MPI__

// constructor from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

PartialComplexVector::PartialComplexVector(MPI::Intracomm& communicator, int id, bool broadcast)
{
  this->VectorType = Vector::ComplexDatas | Vector::PartialData;
  int TmpArray[5];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 5, MPI::INT, id);
  else
    communicator.Recv(TmpArray, 5, MPI::INT, id, 1);   
  this->Dimension = TmpArray[0];
  this->VectorId = TmpArray[1];
  this->IndexShift = TmpArray[3];
  this->RealDimension = TmpArray[4];
  this->Components = new Complex [this->Dimension + 1];
  if (TmpArray[2] == 1)
    for (int i = 0; i <= this->Dimension; ++i) 
      {
	this->Components[i] = 0.0;
      }
  else
    {
      if (TmpArray[2] == 2)
	{
	  if (broadcast == true)
	    {
	      communicator.Bcast(this->Components, 2*this->Dimension, MPI::DOUBLE, id);      
	    }
	  else
	    {
	      communicator.Recv(this->Components, 2*this->Dimension, MPI::DOUBLE, id, 1);   
	    }
	}
      else
	{
	  if (TmpArray[2] == 3) // using Scatterv
	    communicator.Scatterv(/*sendbuf*/ NULL, /* *sendcounts */NULL, /* *displs */ NULL, /*sendtype*/ 0, this->Components, 2*this->Dimension, MPI::DOUBLE, id);
	  else
	    {
	      cout << "Unknown mode for creating scattered data"<<endl;
	      exit(1);
	    }
	}
    }
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
}

#endif


// destructor
//

PartialComplexVector::~PartialComplexVector ()
{
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

PartialComplexVector& PartialComplexVector::operator = (const PartialComplexVector& vector)
{
  if ((this->Dimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Components;
    }
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.Dimension;
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

PartialComplexVector& PartialComplexVector::Copy (PartialComplexVector& vector)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  this->Localize();
  vector.Localize();
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i].Re = vector.Components[i].Re;
      this->Components[i].Im = vector.Components[i].Im;
    }
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

PartialComplexVector& PartialComplexVector::Copy (PartialComplexVector& vector, double coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  this->Localize();
  vector.Localize();
  this->IndexShift = vector.IndexShift;
  this->RealDimension = vector.RealDimension;
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i].Re = vector.Components[i].Re * coefficient;
      this->Components[i].Im = vector.Components[i].Im * coefficient;
    }
  this->Delocalize();
  vector.Delocalize();
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* PartialComplexVector::EmptyClone(bool zeroFlag)
{
  return new PartialComplexVector(this->Dimension, this->Dimension, 0, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors

Vector* PartialComplexVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  PartialComplexVector* TmpVectors = new PartialComplexVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = PartialComplexVector(this->Dimension, this->Dimension, 0, zeroFlag);
  return TmpVectors;
}


// put select vector components to zero
//
// start = start index
// nbrComponent = number of components to set to zero
// return value = reference on current vector

Vector& PartialComplexVector::ClearVectorSegment (long start, long nbrComponent)
{
  start -= this->IndexShift;
  nbrComponent += start;
  if (nbrComponent >= ((long) this->Dimension))
    nbrComponent = (long) this->Dimension;     
  for (;start < nbrComponent; ++ start)
    this->Components[start] = 0.0;  
  return *this;
}

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool PartialComplexVector::WriteVector (char* fileName)
{
  this->Localize();
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Dimension);
  WriteLittleEndian(File, this->RealDimension);
  WriteLittleEndian(File, this->IndexShift);
  WriteBlockLittleEndian(File, &(this->Components[0].Re), 2*this->Dimension);
//   for (int i = 0; i < this->Dimension; ++i)
//     {
//       WriteLittleEndian(File, this->Components[i].Re);
//       WriteLittleEndian(File, this->Components[i].Im);
//     }
  File.close();
  this->Delocalize();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool PartialComplexVector::ReadVector (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);
  this->Resize(TmpDimension);
  ReadLittleEndian(File, this->RealDimension);
  ReadLittleEndian(File, this->IndexShift);
  ReadBlockLittleEndian(File, &(this->Components[0].Re), 2*this->Dimension);
//   for (int i = 0; i < this->Dimension; ++i)
//     {
//       ReadLittleEndian(File, this->Components[i].Re);
//       ReadLittleEndian(File, this->Components[i].Im);
//     }
  File.close();
  return true;
}

#ifdef __MPI__

// send a vector to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& PartialComplexVector::SendVector(MPI::Intracomm& communicator, int id)
{
  communicator.Send(&this->VectorType, 1, MPI::INT, id, 1);
  communicator.Send(&this->Dimension, 1, MPI::INT, id, 1); 
  int Acknowledge = 0;
  communicator.Recv(&Acknowledge, 1, MPI::INT, id, 1);
  if (Acknowledge != 0)
    return *this;
  communicator.Send(&this->IndexShift, 1, MPI::INT, id, 1);
  communicator.Send(&this->RealDimension, 1, MPI::INT, id, 1); 
  communicator.Send(this->Components, 2*this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// broadcast a vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// return value = reference on the current vector

Vector& PartialComplexVector::BroadcastVector(MPI::Intracomm& communicator,  int id)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  communicator.Bcast(&this->IndexShift, 1, MPI::INT, id);
  communicator.Bcast(&this->RealDimension, 1, MPI::INT, id); 
  if (this->VectorType != TmpVectorType)
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    return *this;
  if (TmpDimension != this->Dimension)
    {
      this->Resize(TmpDimension);      
    }
  communicator.Bcast(this->Components, 2*this->Dimension, MPI::DOUBLE, id); 
  return *this;
}

// broadcast part of vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
// nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
// return value = reference on the current vector

Vector& PartialComplexVector::BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  communicator.Bcast(&firstComponent, 1, MPI::INT, id);
  communicator.Bcast(&nbrComponent, 1, MPI::INT, id);
  communicator.Bcast(&this->IndexShift, 1, MPI::INT, id);
  communicator.Bcast(&this->RealDimension, 1, MPI::INT, id); 
  if (this->VectorType != TmpVectorType)
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    return *this;
  if (TmpDimension != this->Dimension)
    {
      this->Resize(TmpDimension);      
    }
  communicator.Bcast(this->Components + firstComponent, 2*nbrComponent, MPI::DOUBLE, id);
  return *this;
}

// receive a vector from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current vector

Vector& PartialComplexVector::ReceiveVector(MPI::Intracomm& communicator, int id)
{
  int TmpVectorType = 0;
  int TmpDimension = 0;
  int TmpRealDimension = 0;
  communicator.Recv(&TmpVectorType, 1, MPI::INT, id, 1);
  communicator.Recv(&TmpDimension, 1, MPI::INT, id, 1); 
  if (TmpVectorType != this->VectorType)
    {
      TmpDimension = 1;
      communicator.Send(&TmpDimension, 1, MPI::INT, id, 1);
      return *this;
    }
  else
    {
      if (TmpDimension != this->Dimension)
	{
	  this->Resize(TmpDimension);      
	}
      TmpDimension = 0;
      communicator.Send(&TmpDimension, 1, MPI::INT, id, 1);
    }
  communicator.Recv(&this->IndexShift, 1, MPI::INT, id, 1);
  communicator.Recv(&this->RealDimension, 1, MPI::INT, id, 1); 
  communicator.Recv(this->Components, 2*this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// add current vector to the current vector of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& PartialComplexVector::SumVector(MPI::Intracomm& communicator, int id)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  if ((this->VectorType != TmpVectorType) || (this->Dimension != TmpDimension))
    {
      Acknowledge = 1;
    }
  if (id != communicator.Get_rank())
    communicator.Send(&Acknowledge, 1, MPI::INT, id, 1);      
  else
    {
      int NbrMPINodes = communicator.Get_size();
      bool Flag = false;
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    communicator.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	    if (Acknowledge == 1)
	      Flag = true;
	  }
      if (Flag == true)
	Acknowledge = 1;
    }
  communicator.Bcast(&Acknowledge, 1, MPI::INT, id);
  if (Acknowledge != 0)
    {
      return *this;
    }
  Complex* TmpComponents = 0;
  if (id == communicator.Get_rank())
    {
      TmpComponents = new Complex [this->Dimension];
    }
  communicator.Bcast(&this->IndexShift, 1, MPI::INT, id);
  communicator.Bcast(&this->RealDimension, 1, MPI::INT, id); 
  communicator.Reduce(this->Components, TmpComponents, 2*this->Dimension, MPI::DOUBLE, MPI::SUM, id); 
  if (id == communicator.Get_rank())
    {
      for (int i = 0; i < this->Dimension; ++i)
	this->Components[i] = TmpComponents[i];
      delete[] TmpComponents;
    }
  return *this;
}

// reassemble vector from a scattered one
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& PartialComplexVector::ReassembleVector(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      int NbrMPINodes = communicator.Get_size();
      int TmpArray[2];
      for (int i = 0; i < NbrMPINodes; ++i)
	if (id != i)
	  {
	    TmpArray[0] = 0;
	    TmpArray[1] = 0;
	    communicator.Recv(TmpArray, 2, MPI::INT, i, 1); 
	    TmpArray[0] -= this->IndexShift;
	    communicator.Recv(this->Components + TmpArray[0], 2 * TmpArray[1], MPI::DOUBLE, id, 1);   	    
	  }      
    }
  else
    {
      int TmpArray[2];
      TmpArray[0] = this->IndexShift;
      TmpArray[1] = this->Dimension;
      communicator.Send(TmpArray, 2, MPI::INT, id, 1);
      communicator.Send(this->Components, 2 * this->Dimension, MPI::DOUBLE, id, 1);  
    }
  return *this;
}

// create a new vector on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* PartialComplexVector::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[5];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 2;
      TmpArray[3] = this->IndexShift;
      TmpArray[4] = this->RealDimension;
      communicator.Bcast(TmpArray, 5, MPI::INT, id);      
      communicator.Bcast(this->Components, 2*this->Dimension, MPI::DOUBLE, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new PartialComplexVector(communicator, id);
    }
  return 0;
}

// create a new vector on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* PartialComplexVector::BroadcastEmptyClone(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[5];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 0;
      TmpArray[3] = this->IndexShift;
      TmpArray[4] = this->RealDimension;
      if (zeroFlag == true)
	{
	  TmpArray[2] = 1;
	}
      communicator.Bcast(TmpArray, 5, MPI::INT, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new PartialComplexVector(communicator, id);
    }
  return 0;
}

#endif
