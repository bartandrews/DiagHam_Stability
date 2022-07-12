////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for n dimensional real vector                     //
//            whose memory allocation is done only on one process             //
//                                                                            //
//                        last modification : 15/06/2004                      //
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


#ifdef USE_CLUSTER_ARCHITECTURE

#include "Vector/DelocalizedRealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "GeneralTools/ListIterator.h"

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

DelocalizedRealVector::DelocalizedRealVector()
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->Dimension = 0;
  this->TrueDimension = 0;
  this->Components = 0;
  this->VectorId = 0;
  this->LocalizationId = 0;
  this->LocalId = 0;
  this->Architecture = 0;
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
}

// constructor for an empty real vector (all coordinates set to zero)
//
// size = Vector Dimension 
// architecture = pointer to the cluster architecture in use
// vectorId = id of the vector
// localizationId = id of the process where the vector is localized
// zeroFlag = true if all coordinates have to be set to zero

DelocalizedRealVector::DelocalizedRealVector(int size, AbstractClusterArchitecture* architecture, int vectorId, 
					     int localizationId, bool zeroFlag)
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
  this->Architecture = architecture;
  this->LocalizationId = localizationId;
  this->LocalId = this->Architecture->GetProcessId();
  this->VectorId = vectorId;
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
  if (this->LocalizationId == this->LocalId)
    {
      this->Components = new double [this->Dimension + 1]; 
      if (zeroFlag == true)
	for (int i = 0; i < this->Dimension; i++)
	  {
	    this->Components[i] = 0.0;
	  }
    }
  else
    {
      this->Components = 0;
    }
}

// copy constructor
//
// vector = vector to copy
// DuplicateFlag = true if datas have to be duplicated

DelocalizedRealVector::DelocalizedRealVector(const DelocalizedRealVector& vector, bool duplicateFlag)
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->Architecture = vector.Architecture;
  this->DummyElementPosition = vector.DummyElementPosition;
  this->DummyElement = vector.DummyElement;
  this->FlushDummyElement();
  if (duplicateFlag == false)
    {
      this->Components = vector.Components;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Components = new double [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    this->Components[i] = vector.Components[i];
	}
      else
	{
	  this->Components = 0;
	}
      this->Flag.Initialize();
    }
}

// copy constructor
//
// vector = vector to copy
// architecture = pointer to the cluster architecture in use
// DuplicateFlag = true if datas have to be duplicated

DelocalizedRealVector::DelocalizedRealVector(const RealVector& vector, AbstractClusterArchitecture* architecture, bool duplicateFlag)
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->VectorId = vector.VectorId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->Architecture = architecture;
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
  if (duplicateFlag == false)
    {
      this->Components = vector.Components;
      this->Flag = vector.Flag;
    }
  else
    {
      if (vector.Dimension > 0)
	{
	  this->Components = new double [this->TrueDimension + 1]; 
	  for (int i = 0; i < this->Dimension; i++)
	    this->Components[i] = vector.Components[i];
	}
      else
	{
	  this->Components = 0;
	}
      this->Flag.Initialize();
    }
}

// copy constructor from a complex vector (keep only real part and datas are duplicated)
//
// vector = vector to copy
// architecture = pointer to the cluster architecture in use

DelocalizedRealVector::DelocalizedRealVector(const ComplexVector& vector, AbstractClusterArchitecture* architecture)
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->Architecture = architecture;
  this->LocalId = this->Architecture->GetProcessId();
  this->LocalizationId = this->LocalId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
  this->VectorId = 0;
  if (this->Dimension > 0)
    {
      this->Components = new double[this->Dimension + 1];
      for (int i = 0; i < this->Dimension; ++i)
	{
	  this->Components[i] = vector.RealComponents[i];
	}
    }
  else
    this->Components = 0;
  this->Flag.Initialize();
}

// copy constructor from a vector (duplicate datas if necessary)
//
// vector = vector to copy
// architecture = pointer to the cluster architecture in use

DelocalizedRealVector::DelocalizedRealVector(const Vector& vector, AbstractClusterArchitecture* architecture)
{
  this->VectorType = Vector::RealDatas | Vector::NonLocalDatas;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->Architecture = architecture;
  this->LocalId = this->Architecture->GetProcessId();
  switch (vector.VectorType)
    {
    case (Vector::RealDatas):
      {
	this->VectorId = vector.VectorId;
	this->LocalizationId = this->LocalId;
	this->Components = ((RealVector&) vector).Components;
	this->Flag = ((RealVector&) vector).Flag;
	this->DummyElementPosition = new int;
	(*(this->DummyElementPosition)) = -1;
      }
      break;
    case (Vector::RealDatas | Vector::NonLocalDatas):
      {
 	this->VectorId = vector.VectorId;
	this->LocalizationId = ((DelocalizedRealVector&) vector).LocalId;
	this->DummyElementPosition = ((DelocalizedRealVector&) vector).DummyElementPosition;
	this->DummyElement = ((DelocalizedRealVector&) vector).DummyElement;
	this->FlushDummyElement();
	this->Components = ((DelocalizedRealVector&) vector).Components;
	this->Flag = ((DelocalizedRealVector&) vector).Flag;
      }
      break;
    case (Vector::ComplexDatas):
      {
	this->DummyElementPosition = new int;
	(*(this->DummyElementPosition)) = -1;
	this->VectorId = 0;
	if (this->Dimension > 0)
	  {
	    this->Components = new double[this->Dimension + 1];
	    for (int i = 0; i < this->Dimension; ++i)
	      {
		this->Components[i] = ((ComplexVector&) vector).RealComponents[i];
	      }
	  }
	else
	  this->Components = 0;
      }
      break;
    }
  if (vector.VectorType == Vector::RealDatas)
    {
      this->VectorType = Vector::RealDatas;
    }
  else
    if (vector.VectorType == Vector::ComplexDatas)
      {
	if (this->Dimension > 0)
	  {
	    this->Components = new double[this->Dimension + 1];
	    for (int i = 0; i < this->Dimension; ++i)
	      {
		this->Components[i] = ((ComplexVector&) vector).RealComponents[i];
	      }
	  }
	else
	  this->Components = 0;
	this->Flag.Initialize();
      }
    else
      {
	this->Components = 0;
	this->Flag.Initialize();
      }
}

// destructor
//

DelocalizedRealVector::~DelocalizedRealVector ()
{
  this->FlushDummyElement();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Dimension != 0)
	delete[] this->Components;
      delete this->DummyElementPosition;
    }
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

DelocalizedRealVector& DelocalizedRealVector::operator = (const DelocalizedRealVector& vector)
{
  this->FlushDummyElement();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Dimension != 0)
	delete[] this->Components;
      delete this->DummyElementPosition;
    }
  this->DummyElementPosition = vector.DummyElementPosition;
  this->DummyElement = vector.DummyElement;
  this->FlushDummyElement();
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->LocalizationId = vector.LocalizationId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.Dimension;
  return *this;
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

DelocalizedRealVector& DelocalizedRealVector::operator = (const RealVector& vector)
{
  this->FlushDummyElement();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Dimension != 0)
	delete[] this->Components;
      delete this->DummyElementPosition;
    }
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
  this->Flag = vector.Flag;
  this->VectorId = vector.VectorId;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->LocalizationId = this->LocalId;
  this->TrueDimension = vector.Dimension;
  return *this;
}

// assignement from a complex vector (keep only real part and datas are duplicated)
//
// vector = vector to assign
// return value = reference on current vector

DelocalizedRealVector& DelocalizedRealVector::operator = (const ComplexVector& vector)
{
  this->FlushDummyElement();
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Dimension != 0)
	delete[] this->Components;
      delete this->DummyElementPosition;
    }
  this->LocalizationId = this->LocalId;
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  this->Components = new double[this->Dimension + 1];
  this->VectorId = 0;
  this->DummyElementPosition = new int;
  (*(this->DummyElementPosition)) = -1;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i] = vector.RealComponents[i];
    }
  this->Flag.Initialize();
  return *this;
}

// assignement from a vector (duplicate datas if necessary)
//
// vector = vector to assign
// return value = reference on current vector

DelocalizedRealVector& DelocalizedRealVector::operator = (const Vector& vector)
{
  switch (vector.VectorType)
    {
    case (Vector::RealDatas):
      return ((*this) = (RealVector&) vector);
      break;
    case (Vector::RealDatas | Vector::NonLocalDatas):
      return ((*this) = (DelocalizedRealVector&) vector);
      break;
    case (Vector::ComplexDatas):
      return ((*this) = (ComplexVector&) vector);
      break;
    }
  return *this;
}

// Resize vector
//
// dimension = new dimension

void DelocalizedRealVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  this->FlushDummyElement();
  if (this->LocalizationId == this->LocalId)
    {
      double* TmpVector = new double [dimension + 1];
      for (int i = 0; i < this->Dimension; i++)
	TmpVector[i] = this->Components[i];
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->Components;
	}
      this->DummyElementPosition = new int;
      (*(this->DummyElementPosition)) = -1;
      this->Dimension = dimension;
      this->TrueDimension = dimension;
      this->Components = TmpVector;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
    }
  else
    {
      this->Localize();
      double* TmpVector = new double [dimension + 1];
      for (int i = 0; i < this->Dimension; i++)
	TmpVector[i] = this->Components[i];
      delete[] this->Components;
      this->Dimension = dimension;
      this->TrueDimension = dimension;
      this->Components = TmpVector;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
      this->Delocalize(true);
    }
  return;
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void DelocalizedRealVector::ResizeAndClean (int dimension)
{
  if ((dimension <= this->TrueDimension) || (this->LocalizationId != this->LocalId))
    {
      this->Dimension = dimension;
      return;
    }
  this->FlushDummyElement();
  if (this->LocalizationId == this->LocalId)
    {
      double* TmpVector = new double [dimension + 1];
      for (int i = 0; i < this->Dimension; i++)
	TmpVector[i] = this->Components[i];
      for (int i = this->Dimension; i < dimension; i++)
	TmpVector[i] = 0.0;  
      if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
	{
	  delete[] this->Components;
	}
      this->Dimension = dimension;
      this->TrueDimension = dimension;
      this->Components = TmpVector;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
    }
  else
    {
      this->Localize();
      double* TmpVector = new double [dimension + 1];
      for (int i = 0; i < this->Dimension; i++)
	TmpVector[i] = this->Components[i];
      for (int i = this->Dimension; i < dimension; i++)
	TmpVector[i] = 0.0;  
      delete[] this->Components;
      this->Dimension = dimension;
      this->TrueDimension = dimension;
      this->Components = TmpVector;
      this->Flag = GarbageFlag();
      this->Flag.Initialize();
      this->Delocalize(true);
    }
  return;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* DelocalizedRealVector::EmptyClone(bool zeroFlag)
{
  this->FlushDummyElement();
  return new DelocalizedRealVector(this->Dimension, this->Architecture, this->VectorId, this->LocalizationId, zeroFlag);
}

// localize the current vector to the current process
// 

void DelocalizedRealVector::Localize()
{
  if (this->LocalizationId != this->LocalId)
    {
      this->FlushDummyElement();
      int TmpLocalizationId = this->LocalizationId;
      RealVector TmpVector (this->Architecture->GetRealVector(this->VectorId));
      TmpVector.SetVectorId(this->VectorId);
      *this = TmpVector;
      this->LocalizationId = TmpLocalizationId;
    }
}

// delocalize the current vector from the current process
// 
// transfertFlag = indicates if the current vector datas have to sent to the vector real location

void DelocalizedRealVector::Delocalize(bool transfertFlag)
{
  if (this->LocalizationId != this->LocalId)
    {
      this->FlushDummyElement();
      if (transfertFlag == true)
	{
	  this->Architecture->SetRealVector(*this, this->VectorId);
	}
      delete[] this->Components;
    }
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool DelocalizedRealVector::ReadVector (char* fileName)
{
  if (this->LocalizationId == this->LocalId)
    {
      ifstream File;
      File.open(fileName, ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "Cannot open the file: " << fileName << endl;
	  return false;
	}
      int TmpDimension;
      File.read ((char*) &(TmpDimension), sizeof(int));
      this->Resize(TmpDimension);
      for (int i = 0; i < this->Dimension; ++i)
	File.read ((char*) (&(this->Components[i])), sizeof(double));
      File.close();
      return true;
    }
  else
    {
      RealVector TmpRealVector;
      TmpRealVector.ReadVector(fileName);
      this->Architecture->SetRealVector(TmpRealVector, this->VectorId);
      this->Resize(TmpRealVector.Dimension);
      return true;
    }
}

// input file stream overload
//
// file = reference on input file stream
// vector = reference on vector to save
// return value = reference on output file stream

ifstream& operator >> (ifstream& file, DelocalizedRealVector& vector)
{
  if (vector.LocalizationId == vector.LocalId)
    {
      file.read ((char*) &(vector.Dimension), sizeof(int));
      if (vector.Dimension > 0)
	{
	  vector.Resize(vector.Dimension);
	  file.read ((char*) (vector.Components), sizeof(double) * vector.Dimension);
	}
      else
	{
	  vector.Dimension = 0;
	  vector.Resize(vector.Dimension);      
	}
    }
  else
    {
      RealVector TmpRealVector;
      file >> TmpRealVector;
      vector.Architecture->SetRealVector(TmpRealVector, vector.VectorId);
      vector.Resize(TmpRealVector.GetVectorDimension());
   }
  return file;
}

#endif
