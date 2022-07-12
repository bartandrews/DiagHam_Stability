////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                   class for n dimensional complex vector                   //
//                                                                            //
//                        last modification : 17/01/2001                      //
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


#include "Vector/ComplexVector.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "GeneralTools/Endian.h"


#include <fstream>

using std::cout;
using std::ofstream;
using std::ifstream;
using std::ios;
using std::endl;

  
// default constructor
//

ComplexVector::ComplexVector() 
{
  this->Components = 0;
  this->Dimension = 0;
  this->VectorId = 0;
  this->TrueDimension = this->Dimension;
  this->VectorType = Vector::ComplexDatas;
}

// constructor for an empty real vector
//
// size = Vector Dimension 
// zeroFlag = true if all coordinates have to be set to zero

ComplexVector::ComplexVector(int size, bool zeroFlag) 
{
  this->Dimension = size;
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
  this->Components = new Complex [this->Dimension];
  this->VectorType = Vector::ComplexDatas;
  this->VectorId = 0;
  if (zeroFlag == true)
    for (int i = 0; i < this->Dimension; ++i)
      this->Components[i] = 0.0;

}

// constructor from arrays of doubles
//
// real = array of doubles corresponding to real part
// imaginary = array of doubles corresponding to imaginary part
// size = Vector Dimension 

ComplexVector::ComplexVector(double* real, double* imaginary, int size) 
{
  this->Flag.Initialize();
  this->Dimension = size;
  this->VectorType = Vector::ComplexDatas;
  this->VectorId = 0;
  this->TrueDimension = this->Dimension;
  this->Components = new Complex [this->Dimension];
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re = real[i];
      this->Components[i].Im = imaginary[i];
    }
  // delete [] real;
  // delete [] imaginary;
}

// constructor from Complex array
//
// components = array of Complex values
// size = Vector Dimension 
ComplexVector::ComplexVector(Complex *components, int size)
{
  this->Flag.Initialize();
  this->Dimension = size;
  this->VectorType = Vector::ComplexDatas;
  this->VectorId = 0;
  this->TrueDimension = this->Dimension;
  this->Components = components;
}

// copy constructor
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

ComplexVector::ComplexVector(const ComplexVector& vector, bool duplicateFlag) 
{
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->VectorType = Vector::ComplexDatas;
  this->VectorId = vector.VectorId;
  if (vector.Dimension == 0)
    {
      this->Flag.Initialize();
      this->Components = 0;
    }
  else
    if (duplicateFlag == false)
      {
	this->Flag = vector.Flag;
	this->Components = vector.Components;
      }
    else
      {
	this->Flag.Initialize();
	this->Components = new Complex [this->TrueDimension];
	for (int i = 0; i < this->Dimension; i++)
	  this->Components[i] = vector.Components[i];
      }
}

// copy constructor from a real vector
//
// vector = vector to copy
// duplicateFlag = true if datas have to be duplicated

ComplexVector::ComplexVector(const RealVector& vector, bool duplicateFlag) 
{
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->VectorType = Vector::ComplexDatas;
  this->VectorId = 0;
  if (vector.Dimension == 0)
    {
      this->Flag.Initialize();
      this->Components = 0;
    }
  else
    {
      this->Flag.Initialize();
      this->Components = new Complex [this->Dimension];
      for (int i = 0; i < this->Dimension; ++i)
	this->Components[i] = vector.Components[i];
    }
}


#ifdef __MPI__

// constructor from informations sent using MPI
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts or sends the vector
// broadcast = true if the vector is broadcasted

ComplexVector::ComplexVector(MPI::Intracomm& communicator, int id, bool broadcast)
{
  this->VectorType = Vector::ComplexDatas;
  int TmpArray[3];
  if (broadcast == true)
    communicator.Bcast(TmpArray, 3, MPI::INT, id);      
  else
    communicator.Recv(TmpArray, 3, MPI::INT, id, 1);   
  this->Dimension = TmpArray[0];
  this->VectorId = TmpArray[1];
  this->Components = new Complex [this->Dimension + 1];
  if (TmpArray[2] == 1)
    for (int i = 0; i <= this->Dimension; ++i) 
      this->Components[i] = 0.0;
  else
    if (TmpArray[2] == 2)
      {
	if (broadcast == true)
	  communicator.Bcast(this->Components, 2 * this->Dimension, MPI::DOUBLE, id);      
	else
	  communicator.Recv(this->Components, 2 * this->Dimension, MPI::DOUBLE, id, 1);   
      }
  this->TrueDimension = this->Dimension;
  this->Flag.Initialize();
}

#endif

// destructor
//

ComplexVector::~ComplexVector () 
{
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
}

// assignement
//
// vector = vector to assign
// return value = reference on current vector

ComplexVector& ComplexVector::operator = (const ComplexVector& vector) 
{
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag = vector.Flag;
  this->Components = vector.Components;
  this->Dimension = vector.Dimension;
  this->TrueDimension = vector.TrueDimension;
  this->VectorId = vector.VectorId;
  return *this;
}

// Resize vector
//
// dimension = new dimension

void ComplexVector::Resize (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      return;
    }
  Complex* TmpVector = new Complex [dimension];
  int i = 0;
  for (; i < this->Dimension; ++i)
    TmpVector[i] = this->Components[i];
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
}

// Resize vector and set to zero all components that have been added
//
// dimension = new dimension

void ComplexVector::ResizeAndClean (int dimension)
{
  if (dimension <= this->TrueDimension)
    {
      this->Dimension = dimension;
      for (int i = this->Dimension; i < dimension; ++i)
	this->Components[i] = 0.0;  
      return;
    }
  Complex* TmpVector = new Complex [dimension];
  int i = 0;
  for (; i < this->Dimension; ++i)
    TmpVector[i] = this->Components[i];
  for (; i < dimension; i++)
    TmpVector[i] = 0.0;  
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->Dimension = dimension;
  this->TrueDimension = dimension;
  this->Components = TmpVector;
}

// assignement from a real vector
//
// vector = vector to assign
// return value = reference on current vector

ComplexVector& ComplexVector::operator = (const RealVector& vector) 
{
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Dimension = vector.Dimension;
  this->TrueDimension = this->Dimension;
  if (vector.Dimension == 0)
    {
      this->Components = 0;
    }
  this->Components = new Complex [this->Dimension];
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i] = vector.Components[i];
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

ComplexVector& ComplexVector::Copy (ComplexVector& vector, double coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  Complex *Ci=this->Components;
  Complex *Cj=vector.Components;
  for (int i = 0; i < this->Dimension; i++, Ci++, Cj++)
    {
      Ci->Re = Cj->Re * coefficient;
      Ci->Im = Cj->Im * coefficient;
    }
  return *this;
}

// copy a vector into another
//
// vector = vector to copy
// coefficient = optional coefficient which multiply source to copy
// return value = reference on current vector

ComplexVector& ComplexVector::Copy (ComplexVector& vector, const Complex& coefficient)
{
  if (this->Dimension != vector.Dimension)
    this->Resize(vector.Dimension);
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i].Re = vector.Components[i].Re * coefficient.Re - vector.Components[i].Im * coefficient.Im;
      this->Components[i].Im = vector.Components[i].Re * coefficient.Im + vector.Components[i].Im * coefficient.Re;
    }
  return *this;
}

// create a new vector with same size and same type but non-initialized components
//
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* ComplexVector::EmptyClone(bool zeroFlag)
{
  return new ComplexVector(this->Dimension, zeroFlag);
}

// create an array of new vectors with same size and same type but non-initialized components
//
// nbrVectors = number of vectors to sreate
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to the array of new vectors

Vector* ComplexVector::EmptyCloneArray(int nbrVectors, bool zeroFlag)
{
  ComplexVector* TmpVectors = new ComplexVector [nbrVectors];
  for (int i = 0; i < nbrVectors; ++i)
    TmpVectors[i] = ComplexVector(this->Dimension, zeroFlag);
  return TmpVectors;
}

// put all vector components to zero
//
// return value = reference on current vector

Vector& ComplexVector::ClearVector ()
{
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] = 0.0;  
  return *this;
}

// change sign of a vector
//
// return value = reference on current vector

ComplexVector& ComplexVector::operator - () 
{
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i] *= -1;        
  return *this;
}

// return a new vector with opposite sign form a given source vector
//
// V1 = source vector
// return value = new vector

ComplexVector operator - (const ComplexVector& V1) 
{
  if (V1.Dimension > 0)
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i] = -V1.Components[i];
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector
// return value = result of scalar product

Complex operator * (const ComplexVector& V1, const ComplexVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  if (min == 0)
    return Complex();
  Complex tmp (V1.Components[0].Re * V2.Components[0].Re + V1.Components[0].Im * V2.Components[0].Im, 
	       V1.Components[0].Re * V2.Components[0].Im - V1.Components[0].Im * V2.Components[0].Re);
  for (int i = 1; i < min; ++i)
    {
      tmp.Re += V1.Components[i].Re * V2.Components[i].Re + V1.Components[i].Im * V2.Components[i].Im;
      tmp.Im += V1.Components[i].Re * V2.Components[i].Im - V1.Components[i].Im * V2.Components[i].Re;
    }
  return tmp;
}

// scalar product between two vectors
//
// V1 = first vector
// V2 = second vector (real vector)
// return value = result of scalar product

Complex operator * (const ComplexVector& V1, const RealVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  Complex tmp (V1.Components[0].Re * V2.Components[0], - V1.Components[0].Im * V2.Components[0]);
  for (int i = 1; i < min; ++i)
    {
      tmp.Re += V1.Components[i].Re * V2.Components[i];
      tmp.Im -= V1.Components[i].Im * V2.Components[i];
    }
  return tmp;
}

// scalar product between two vectors
//
// V1 = first vector (real vector)
// V2 = second vector
// return value = result of scalar product

Complex operator * (const RealVector& V1, const ComplexVector& V2) 
{
  int min = V1.Dimension;
  if (min > V2.Dimension)
    min = V2.Dimension;
  Complex tmp (V1.Components[0] * V2.Components[0].Re, V2.Components[0].Im * V1.Components[0]);
  for (int i = 1; i < min; i++)
    {
      tmp.Re += V2.Components[i].Re * V1.Components[i];
      tmp.Im += V2.Components[i].Im * V1.Components[i];
    }
  return tmp;
}

// do part of the scalar product between two vectors in a given range of indices
//
// vRight = right vector of the scalar product
// firstComponent = index of the first component to consider
// nbrComponent = number of components to consider
// step = increment between to consecutive indices to consider
// return value = result of the partial scalar product

Complex ComplexVector::PartialScalarProduct (const ComplexVector& vRight, int firstComponent, int nbrComponent, int step)
{
  Complex tmp (this->Components[firstComponent].Re * vRight.Components[firstComponent].Re + 
	       this->Components[firstComponent].Im * vRight.Components[firstComponent].Im, 
	       this->Components[firstComponent].Re * vRight.Components[firstComponent].Im - 
	       this->Components[firstComponent].Im * vRight.Components[firstComponent].Re);
  nbrComponent *= step;
  nbrComponent += firstComponent;
  firstComponent += step;
  for (; firstComponent < nbrComponent; firstComponent += step)
    {
      tmp.Re += this->Components[firstComponent].Re * vRight.Components[firstComponent].Re + 
	this->Components[firstComponent].Im * vRight.Components[firstComponent].Im;
      tmp.Im += this->Components[firstComponent].Re * vRight.Components[firstComponent].Im - 
	this->Components[firstComponent].Im * vRight.Components[firstComponent].Re;
    }
  return tmp;
}

// sum two vectors
//
// V1 = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::operator += (ComplexVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i] += V1.Components[i];
  return *this;
}

// sum two vectors
//
// V1 = real vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::operator += (RealVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i].Re += V1.Components[i];
  return *this;
}

// sum two vectors
//
// vector = vector to add
// return value = reference on current vector

Vector& ComplexVector::operator += (Vector& vector)
{
  if (vector.VectorType == Vector::RealDatas)
    return (*this += ((RealVector&) vector));
  if (vector.VectorType == Vector::ComplexDatas)
    return (*this += ((ComplexVector&) vector));
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x, const ComplexVector& V) 
{
  if ((this->Dimension == 0) || (this->Dimension != V.Dimension))
    {
      cout << "Wrong vector dimension!"<<endl;
      return *this;
    }
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].AddMultiply(V.Components[i], x);
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x, const ComplexVector& V, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->Components[i].AddMultiply(V.Components[i], x);
    }
  return *this;
}

// add a linear combination to a given vector
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x, const ComplexVector& V) 
{
  if ((this->Dimension == 0) || (this->Dimension != V.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re += x.Re * V.Components[i].Re - x.Im * V.Components[i].Im;
      this->Components[i].Im += x.Re * V.Components[i].Im + x.Im * V.Components[i].Re;
    }
  return *this;
}

// add a linear combination to a given vector, for a given range of indices
//
// x = multiplicative coefficient
// V = vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x, const ComplexVector& V, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > V.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->Components[i].Re += x.Re * V.Components[i].Re - x.Im * V.Components[i].Im;
      this->Components[i].Im += x.Re * V.Components[i].Im + x.Im * V.Components[i].Re;
    }
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
						    const ComplexVector& v2)
{
  if ((this->Dimension == 0) || (this->Dimension != v1.Dimension) || (this->Dimension != v2.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re += x1 * v1.Components[i].Re + x2 * v2.Components[i].Re;
      this->Components[i].Im += x1 * v1.Components[i].Im + x2 * v2.Components[i].Im;
    }
  return *this;
}

// add a linear combination of two vectors to a given vector, for a given range of indices
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (double x1, const ComplexVector& v1, double x2, 
						    const ComplexVector& v2, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->Components[i].Re += x1 * v1.Components[i].Re + x2 * v2.Components[i].Re;
      this->Components[i].Im += x1 * v1.Components[i].Im + x2 * v2.Components[i].Im;
    }
  return *this;
}

// add a linear combination of two vectors to a given vector
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
						    const ComplexVector& v2)
{
  if ((this->Dimension == 0) || (this->Dimension != v1.Dimension) || (this->Dimension != v2.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re += x1.Re * v1.Components[i].Re  - x1.Im * v1.Components[i].Im 
	+ x2.Re * v2.Components[i].Re - x2.Im * v2.Components[i].Im;
      this->Components[i].Im += x1.Im * v1.Components[i].Re  + x1.Re * v1.Components[i].Im 
	+ x2.Im * v2.Components[i].Re + x2.Re * v2.Components[i].Im;
    }
  return *this;
}

// add a linear combination of two vectors to a given vector, for a given range of indices
//
// x1 = multiplicative coefficient of first vector
// v1 = first vector to add
// x2 = multiplicative coefficient of first vector
// v2 = first vector to add
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on current vector

ComplexVector& ComplexVector::AddLinearCombination (const Complex& x1, const ComplexVector& v1, const Complex& x2, 
						    const ComplexVector& v2, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if ((LastComponent > this->Dimension) || (LastComponent > v2.Dimension) || 
      (LastComponent > v1.Dimension))
    return *this;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      this->Components[i].Re += x1.Re * v1.Components[i].Re  - x1.Im * v1.Components[i].Im 
	+ x2.Re * v2.Components[i].Re - x2.Im * v2.Components[i].Im;
      this->Components[i].Im += x1.Im * v1.Components[i].Re  + x1.Re * v1.Components[i].Im 
	+ x2.Im * v2.Components[i].Re + x2.Re * v2.Components[i].Im;
    }
  return *this;
}

// substract two vectors
//
// V1 = first vector
// return value = reference on current vector

ComplexVector& ComplexVector::operator -= (const ComplexVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i] -= V1.Components[i];
    }
  return *this;
}

// substract two vectors
//
// V1 = first real vector
// return value = reference on current vector

ComplexVector& ComplexVector::operator -= (const RealVector& V1) 
{
  if ((this->Dimension == 0) || (this->Dimension != V1.Dimension))
    return *this;
  for (int i = 0; i < this->Dimension; i++)
    {
      this->Components[i].Re -= V1.Components[i];
    }
  return *this;
}

// sum two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

ComplexVector operator + (const ComplexVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i].Re + V2.Components[i].Re;
	  TmpComponents[i].Im = V1.Components[i].Im + V2.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// sum two vectors with left one real
//
// V1 = first vector (real)
// V2 = second vector
// return value = resulting vector

ComplexVector operator + (const RealVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i] + V2.Components[i].Re;
	  TmpComponents[i].Im = V2.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// sum two vectors with right one real
//
// V1 = first vector
// V2 = second vector (real)
// return value = resulting vector

ComplexVector operator + (const ComplexVector& V1, const RealVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V2.Components[i] + V1.Components[i].Re;
	  TmpComponents[i].Im = V1.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors
//
// V1 = first vector
// V2 = second vector
// return value = resulting vector

ComplexVector operator - (const ComplexVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i].Re - V2.Components[i].Re;
	  TmpComponents[i].Im = V1.Components[i].Im - V2.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors with left one real
//
// V1 = first vector (real)
// V2 = second vector
// return value = resulting vector

ComplexVector operator - (const RealVector& V1, const ComplexVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i] - V2.Components[i].Re;
	  TmpComponents[i].Im = -V2.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// substract two vectors with rightt one real
//
// V1 = first vector 
// V2 = second vector (real)
// return value = resulting vector

ComplexVector operator - (const ComplexVector& V1, const RealVector& V2) 
{
  if ((V1.Dimension != 0) && (V2.Dimension == V1.Dimension))
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i].Re - V2.Components[i];
	  TmpComponents[i].Im = V1.Components[i].Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the right hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

ComplexVector operator * (const ComplexVector& V1, double d) 
{
  if (V1.Dimension != 0)
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].AddMultiply(V1.Components[i], d);
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a complex number on the right hand side
//
// V1 = vector to multiply
// d = complex to use
// return value = resulting vector

ComplexVector operator * (const ComplexVector& V1, const Complex& d) 
{
  if (V1.Dimension != 0)
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i].Re* d.Re - V1.Components[i].Im * d.Im;
	  TmpComponents[i].Im = V1.Components[i].Im * d.Re + V1.Components[i].Re* d.Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the left hand side
//
// V1 = vector to multiply
// d = real to use
// return value = resulting vector

ComplexVector operator * (double d, const ComplexVector& V1) 
{
  if (V1.Dimension != 0)
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].AddMultiply(V1.Components[i], d);
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a complex number on the left hand side
//
// V1 = vector to multiply
// d = complex to use
// return value = resulting vector

ComplexVector operator * (const Complex& d, const ComplexVector& V1) 
{
  if (V1.Dimension != 0)
    {
      Complex* TmpComponents = new Complex [V1.Dimension];
      for (int i = 0; i < V1.Dimension; ++i)
	{
	  TmpComponents[i].Re = V1.Components[i].Re* d.Re - V1.Components[i].Im * d.Im;
	  TmpComponents[i].Im = V1.Components[i].Im * d.Re + V1.Components[i].Re* d.Im;
	}
      return ComplexVector(TmpComponents, V1.Dimension);
    }
  else
    return ComplexVector();
}

// multiply a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (double d) 
{
  for (int i = 0; i < this->Dimension; ++i)
    this->Components[i] *= d;  
  return *this;
}

// divide a vector with a real number on the right hand side
//
// d = real to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator /= (double d) 
{
  d = 1.0 / d;
  for (int i = 0; i < this->Dimension; i++)
    this->Components[i] *= d;  
  return *this;
}

// multiply a vector with a complex number on the right hand side
//
// d = complex to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const Complex& d) 
{
  double tmp;
  for (int i = 0; i < this->Dimension; ++i)
    {
      tmp = d.Re * this->Components[i].Re - d.Im * this->Components[i].Im;
      this->Components[i].Im *= d.Re;
      this->Components[i].Im += this->Components[i].Re * d.Im;  
      this->Components[i].Re = tmp;  
    }
  return *this;
}

// divide a vector with a complex number on the right hand side
//
// d = complex to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator /= (const Complex& d) 
{
  double fac = 1.0 / (d.Re * d.Re + d.Im * d.Im);
  double dRe = d.Re * fac;
  double dIm = - d.Im * fac;  
  double tmp;
  for (int i = 0; i < this->Dimension; i += 2)
    {
      tmp = dRe * this->Components[i].Re - dIm * this->Components[i].Im;
      this->Components[i].Im *= dRe;
      this->Components[i].Im += this->Components[i].Re * dIm;  
      this->Components[i].Re = tmp;  
    }
  return *this;
}

// left multiply a vector with a real matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const RealMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  Complex* tmp = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      tmp[i].Re = 0.0;
      tmp[i].Im = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp[i].Re += M.Columns[j].Components[i] * this->Components[j].Re;
	  tmp[i].Im += M.Columns[j].Components[i] * this->Components[j].Im;
	}      
    }
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->Components = tmp;
  return *this;
}

// left multiply a vector with a complex matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const ComplexMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrColumn))
    return *this;
  Complex* tmp = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      tmp[i].Re = 0.0;
      tmp[i].Im = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp[i].Re += (M.Columns[j].Components[i].Re * this->Components[j].Re - 
			M.Columns[j].Components[i].Im * this->Components[j].Im);
	  tmp[i].Im += (M.Columns[j].Components[i].Re * this->Components[j].Im + 
			M.Columns[j].Components[i].Im * this->Components[j].Re);
	}      
    }
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->Components = tmp;
  return *this;
}

// left multiply a vector with an hermtian conjugated  complex matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator &= (const ComplexMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this; 
  Complex* tmp = new Complex [M.NbrRow];  
  for (int i = 0; i < M.NbrColumn; i++)
    {
      tmp[i].Re = 0.0;
      tmp[i].Im = 0.0;
      for (int j = 0; j < M.NbrColumn; ++j)
	{
	  tmp[i].Re += (M.Columns[j].Components[i].Re * this->Components[j].Re + 
			M.Columns[j].Components[i].Im * this->Components[j].Im);
	  tmp[i].Im += (M.Columns[j].Components[i].Re * this->Components[j].Im - 
			M.Columns[j].Components[i].Im * this->Components[j].Re);
	}      
    }
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Dimension = M.NbrRow;
  this->Components = tmp;
  return *this;
}

// left multiply a vector with an hermitian matrix (using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const HermitianMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  Complex* tmp = new Complex [M.NbrRow];
  for (int i = 0; i < M.NbrRow; ++i)
    {
      tmp[i].Re = M.DiagonalElements[i] * this->Components[i].Re;
      tmp[i].Im = M.DiagonalElements[i] * this->Components[i].Im;
      int j = 0;
      int pos = i - 1;
      for (; j < i; ++j)
	{
	  tmp[i].Re += M.RealOffDiagonalElements[pos] * this->Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[pos] * this->Components[j].Im;
	  tmp[i].Im += M.RealOffDiagonalElements[pos] * this->Components[j].Im - 
	    M.ImaginaryOffDiagonalElements[pos] * this->Components[j].Re;
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  tmp[i].Re += M.RealOffDiagonalElements[pos] * this->Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[pos] * this->Components[j].Im;
	  tmp[i].Im += M.RealOffDiagonalElements[pos] * this->Components[j].Im + 
	    M.ImaginaryOffDiagonalElements[pos] * this->Components[j].Re;
	  ++pos;
	}
    }
  if ((this->Components != 0) && (this->Flag.Used() == true))
    if (this->Flag.Shared() == false)
      {
	delete[] this->Components;
      }
  this->Flag.Initialize();
  this->TrueDimension = M.NbrRow;
  this->Components = tmp;
  return *this;
}

// left multiply a vector with a complex tridiagonal hermitian matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const ComplexTriDiagonalHermitianMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  int ReducedDim = this->Dimension - 1;
  Complex Tmp1 (this->Components[0].Re, this->Components[0].Im);
  Complex Tmp2;
  this->Components[0].Re *= M.DiagonalElements[0];
  this->Components[0].Re += this->Components[1].Re * M.RealUpperDiagonalElements[0]
    - this->Components[1].Im * M.ImaginaryUpperDiagonalElements[0];
  this->Components[0].Im *= M.DiagonalElements[0];
  this->Components[0].Im += this->Components[1].Im * M.RealUpperDiagonalElements[0]
    + this->Components[1].Re * M.ImaginaryUpperDiagonalElements[0];
  for (int i = 1; i < ReducedDim; i++)
    {
      Tmp2 = Complex (this->Components[i].Re, this->Components[i].Im);
      this->Components[i].Re *= M.DiagonalElements[i];
      this->Components[i].Re += this->Components[i + 1].Re * M.RealUpperDiagonalElements[i] 
	- this->Components[i + 1].Im * M.ImaginaryUpperDiagonalElements[i] + Tmp1.Re * M.RealUpperDiagonalElements[i - 1] 
	+ Tmp1.Im * M.ImaginaryUpperDiagonalElements[i - 1];
      this->Components[i].Im *= M.DiagonalElements[i];
      this->Components[i].Im += this->Components[i + 1].Im * M.RealUpperDiagonalElements[i]
	+ this->Components[i + 1].Re * M.ImaginaryUpperDiagonalElements[i] + Tmp1.Im * M.RealUpperDiagonalElements[i - 1] 
	- Tmp1.Re * M.ImaginaryUpperDiagonalElements[i - 1];
      Tmp1 = Tmp2;
    }
  this->Components[this->Dimension - 1].Re *= M.DiagonalElements[this->Dimension - 1];
  this->Components[this->Dimension - 1].Re += Tmp1.Re * M.RealUpperDiagonalElements[this->Dimension - 2] 
    + Tmp1.Im * M.ImaginaryUpperDiagonalElements[this->Dimension - 2];
  this->Components[this->Dimension - 1].Im *= M.DiagonalElements[this->Dimension - 1];  
  this->Components[this->Dimension - 1].Im += Tmp1.Im * M.RealUpperDiagonalElements[this->Dimension - 2] 
    - Tmp1.Re * M.ImaginaryUpperDiagonalElements[this->Dimension - 2];
  return *this;
}

// left multiply a vector with a real tridiagonal symmetric matrix (without using temporary vector)
//
// M = matrix to use
// return value = reference on current vector

ComplexVector& ComplexVector::operator *= (const RealTriDiagonalSymmetricMatrix&  M)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  int ReducedDim = this->Dimension - 1;
  Complex Tmp1 (this->Components[0].Re, this->Components[1].Im);
  Complex Tmp2;
  this->Components[0].Re *= M.DiagonalElements[0];
  this->Components[0].Re += this->Components[1].Re * M.UpperDiagonalElements[0];
  this->Components[0].Im *= M.DiagonalElements[0];
  this->Components[0].Im += this->Components[1].Im * M.UpperDiagonalElements[0];
  for (int i = 1; i < ReducedDim; i++)
    {
      Tmp2 = Complex (this->Components[i].Re, this->Components[i].Im);
      this->Components[i].Re *= M.DiagonalElements[i];
      this->Components[i].Re += this->Components[i + 1].Re * M.UpperDiagonalElements[i] 
	+ Tmp1.Re * M.UpperDiagonalElements[i - 1]; 
      this->Components[i].Im *= M.DiagonalElements[i];
      this->Components[i].Im += this->Components[i + 1].Im * M.UpperDiagonalElements[i]
	+ Tmp1.Im * M.UpperDiagonalElements[i - 1];
      Tmp1 = Tmp2;
    }
  this->Components[this->Dimension - 1].Re *= M.DiagonalElements[this->Dimension - 1];
  this->Components[this->Dimension - 1].Re += Tmp1.Re * M.UpperDiagonalElements[this->Dimension - 2];
  this->Components[this->Dimension - 1].Im *= M.DiagonalElements[this->Dimension - 1];  
  this->Components[this->Dimension - 1].Im += Tmp1.Im * M.UpperDiagonalElements[this->Dimension - 2];
  return *this;
}

// left multiply a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const HermitianMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re = M.DiagonalElements[i] * V.Components[i].Re;
      this->Components[i].Im = M.DiagonalElements[i] * V.Components[i].Im;
      int pos = i - 1;
      int j = 0;
      for (; j < i; ++j)
	{
	  this->Components[i].Re += M.RealOffDiagonalElements[pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[pos ] * V.Components[j].Im;
	  this->Components[i].Im += M.RealOffDiagonalElements[pos] * V.Components[j].Im 
	    - M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Re;
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  this->Components[i].Re += M.RealOffDiagonalElements[pos] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Im;
	  this->Components[i].Im += M.RealOffDiagonalElements[pos] * V.Components[j].Im + 
	    M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Re;
	  ++pos;
	}
    }
  return *this;
}

// do a partial left multication of a vector with an hermitian matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, 
					int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) || (this->Dimension != M.NbrRow))
    return *this;
  int Last = sourceStart + sourceNbrComponent;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  int j;
  double x;
  double y;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; ++i)
    {
      x = 0.0;
      y = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im 
	    + M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i].Re = x;
      this->Components[i].Im = y;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im 
	    - M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.Components[i].Re;
      y += M.DiagonalElements[i] * V.Components[i].Im;
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos2] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos2] * V.Components[j].Im + 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.Components[j].Re;
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i].Re = x;
      this->Components[i].Im = y;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i].Re = x;
      this->Components[i].Im = y;
   }
  return *this;
}

// left multiply a vector with an hermitian matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const HermitianMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (V.Dimension != M.NbrColumn) || (this->Dimension != M.NbrRow))
    return *this;
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re += M.DiagonalElements[i] * V.Components[i].Re;
      this->Components[i].Im += M.DiagonalElements[i] * V.Components[i].Im;
      int pos = i - 1;
      int j = 0;
      for (; j < i; ++j)
	{
	  this->Components[i].Re += M.RealOffDiagonalElements[pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[pos ] * V.Components[j].Im;
	  this->Components[i].Im += M.RealOffDiagonalElements[pos] * V.Components[j].Im 
	    - M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Re;
	  pos += this->Dimension - j - 2 + M.Increment;
	}
      ++pos;
      ++j;
      for (; j < this->Dimension; ++j)
	{
	  this->Components[i].Re += M.RealOffDiagonalElements[pos] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Im;
	  this->Components[i].Im += M.RealOffDiagonalElements[pos] * V.Components[j].Im + 
	    M.ImaginaryOffDiagonalElements[pos] * V.Components[j].Re;
	  ++pos;
	}
    }
  return *this;
}

// do a partial left multication of a vector with an hermitian matrix and add result to the current vector
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const HermitianMatrix&  M, ComplexVector& V, int sourceStart, 
					   int sourceNbrComponent)
{
  if ((this->Dimension == 0) || (V.Dimension < (sourceNbrComponent + sourceStart)) || (this->Dimension != M.NbrRow))
    return *this;
  int Last = sourceStart + sourceNbrComponent;
  int Inc1 =  this->Dimension - sourceNbrComponent + M.Increment - 2;
  int Inc2 =  this->Dimension + M.Increment - 2;
  int j;
  double x;
  double y;
  int i = 0;
  int Pos = sourceStart - 1;
  for (; i < sourceStart; ++i)
    {
      x = 0.0;
      y = 0.0;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im 
	    + M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  ++Pos;
	}
      Pos += Inc1 - i;
      this->Components[i].Re += x;
      this->Components[i].Im += y;
    }
  Inc1 = this->Dimension - Last + M.Increment;
  int Pos2 = Pos;
  int Pos3 = Pos;
  ++Pos2;
  for (; i < Last; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < i; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im 
	    - M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  Pos += Inc2 - j;
	}
      x += M.DiagonalElements[i] * V.Components[i].Re;
      y += M.DiagonalElements[i] * V.Components[i].Im;
      ++j;
      for (; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos2] * V.Components[j].Re - 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos2] * V.Components[j].Im + 
	    M.ImaginaryOffDiagonalElements[Pos2] * V.Components[j].Re;
	  ++Pos2;
	}
      Pos2 += Inc1; 
      ++Pos3;
      this->Components[i].Re += x;
      this->Components[i].Im += y;
    }
  for (; i < this->Dimension; ++i)
    {
      x = 0.0;
      y = 0.0;
      Pos = Pos3;
      for (j = sourceStart; j < Last; ++j)
	{
	  x += M.RealOffDiagonalElements[Pos] * V.Components[j].Re + 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Im;
	  y += M.RealOffDiagonalElements[Pos] * V.Components[j].Im - 
	    M.ImaginaryOffDiagonalElements[Pos] * V.Components[j].Re;
	  Pos += Inc2 - j;
	}
      ++Pos3;
      this->Components[i].Re += x;
       this->Components[i].Im += y;
   }
  return *this;
}

// left multiply a vector with a real matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const RealMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  if (V.Dimension != M.NbrColumn)
    V.Resize(M.NbrRow);
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re = 0.0;
      this->Components[i].Im = 0.0;
      for (int j = 0; j < V.Dimension; ++j)
	{
	  this->Components[i].Re += M.Columns[j].Components[i] * V.Components[j].Re;
	  this->Components[i].Im += M.Columns[j].Components[i] * V.Components[j].Im;
	}
    }
  return *this;
}

// left multiply a vector with a complex matrix and use to store result in current vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const ComplexMatrix&  M, ComplexVector& V)
{
  if ((this->Dimension == 0) || (this->Dimension != M.NbrRow))
    return *this;
  if (V.Dimension != M.NbrColumn)
    V.Resize(M.NbrColumn);
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re = 0.0;
      this->Components[i].Im = 0.0;
      for (int j = 0; j < V.Dimension; ++j)
	{
	  this->Components[i].Re += (M.Columns[j].Components[i].Re * V.Components[j].Re -  
				      M.Columns[j].Components[i].Im * V.Components[j].Im);
	  this->Components[i].Im += (M.Columns[j].Components[i].Re * V.Components[j].Im + 
					   M.Columns[j].Components[i].Im * V.Components[j].Re);
	}
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V);
      break;
    case (Matrix::ComplexElements):
      return this->Multiply((ComplexMatrix&) M, V);
      break;
    case (Matrix::RealElements):
      return this->Multiply((RealMatrix&) M, V);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add result to current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					int destStart, int destStep)
{
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					   int destStart, int destStep)
{
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceNbrComponent = number of component to take into account in the source vector
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceNbrComponent)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V, sourceStart, sourceNbrComponent);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and use to store result in current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::Multiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					int sourceNbrComponent, int destStart, int destStep)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->Multiply((HermitianMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
  return *this;
}

// left multiply a vector with a matrix and add current 
// vector (without creating temporary vector)
//
// M = matrix to use
// V = vector to multiply
// sourceStart = source vector first coordinate to modify
// sourceStep = step to add to go to the following source vector coordinate
// sourceNbrComponent = number of component to take into account in the source vector
// destStart = destination vector first coordinate to modify
// destStep = step to add to go to the following destination vector coordinate
// return value = reference on current vector

ComplexVector& ComplexVector::AddMultiply (const Matrix&  M, ComplexVector& V, int sourceStart, int sourceStep, 
					   int sourceNbrComponent, int destStart, int destStep)
{
  switch (M.MatrixType)
    {
    case (Matrix::ComplexElements | Matrix::Hermitian):
      return this->AddMultiply((HermitianMatrix&) M, V, sourceStart, sourceStep, sourceNbrComponent, destStart, destStep);
      break;
    default:
      return *this;
    }
  return *this;
}

// get vector norm
//
// return value = vector norm

double ComplexVector::Norm() 
{
  if (this->Dimension == 0)
    return 0.0;
  double tmp = this->Components[0].Re * this->Components[0].Re + this->Components[0].Im * this->Components[0].Im;
  for (int i = 1; i  < this->Dimension; ++i)
    tmp += this->Components[i].Re * this->Components[i].Re + this->Components[i].Im * this->Components[i].Im;
  return sqrt(tmp);
}
  
// get square of vector norm
//
// return value = square of vector norm

double ComplexVector::SqrNorm () 
{
  if (this->Dimension == 0)
    return 0.0;
  double tmp = this->Components[0].Re * this->Components[0].Re + this->Components[0].Im * this->Components[0].Im;
  for (int i = 2; i  < this->Dimension; ++i)
    tmp += this->Components[i].Re * this->Components[i].Re + this->Components[i].Im * this->Components[i].Im;
  return tmp;
}
  
// normalize vector
//
// return value = reference on current vector

ComplexVector& ComplexVector::Normalize()
{
  if (this->Dimension == 0)
    return *this;
  double tmp = this->Components[0].Re * this->Components[0].Re + this->Components[0].Im * this->Components[0].Im;
  for (int i = 2; i  < this->Dimension; ++i)
    tmp += this->Components[i].Re * this->Components[i].Re + this->Components[i].Im * this->Components[i].Im;
  tmp = 1.0 / sqrt(tmp);
  for (int i = 0; i < this->Dimension; ++i)
    {
      this->Components[i].Re *= tmp;
      this->Components[i].Im *= tmp;
    }
  return *this;
}
  
// Extract a subvector from a given vector
//
// firstCoordinate = Coordinate where extraction has to begin
// lastCoordinate = Coordinate where extraction has to stop (extract also include this last coordinate)
// Step = distance to the next coordinate (1 means to take the following)
// return value = return corresponding subvector

ComplexVector ComplexVector::Extract(int firstCoordinate, int lastCoordinate, int step) 
{
  if (this->Dimension == 0)
    return ComplexVector();
  int TmpDimension = (lastCoordinate - firstCoordinate) / step;
  int TrueLast = (firstCoordinate + TmpDimension * step);
  TmpDimension++;
  Complex* TmpComponents = new Complex [TmpDimension];
  int j = 0;
  for (int i = firstCoordinate; i <= TrueLast; i += step)
    {
      TmpComponents[j] = this->Components[i];
      ++j;
    }
  return ComplexVector(TmpComponents, TmpDimension);  
}
  
// Merge a subvector into a given vector
//
// V = vector to merge
// firstCoordinate = Coordinate where merge has to begin
// step = distance to the next coordinate in the destination vector (1 means to take the following)
// return value = reference to the current Vector

ComplexVector& ComplexVector::Merge(const ComplexVector& V, int firstCoordinate, int step) 
{
  if ((this->Dimension == 0) || (V.Dimension == 0))
    return *this;
  int Max = firstCoordinate + (V.Dimension * step);
  if (Max > this->Dimension)
    return *this;
  int  j = 0;
  for (int i = firstCoordinate; i < Max; i += step)
    {
      this->Components[i] = V.Components[j];
      ++j;
    }
  return *this;
}
  
// Merge a real subvector into a complex given vector
//
// V = real vector to merge
// firstCoordinate = Coordinate where merge has to begin
// step = distance to the next coordinate in the destination vector (1 means to take the following)
// return value = reference to the current Vector

ComplexVector& ComplexVector::Merge(const RealVector& V, int firstCoordinate, int step) 
{
  if ((this->Dimension == 0) || (V.Dimension == 0))
    return *this;
  int Max = firstCoordinate + (V.Dimension * step);
  if (Max > this->Dimension)
    return *this;
  int  j = 0;
  for (int i = firstCoordinate; i < Max; i += step)
    {
      this->Components[i] = V.Components[j];
      ++j;
    }
  return *this;
}

// reverse elements of the current vector (i.e. exchanging i <-> N - i)
//
// return value = reference to the current Vector

ComplexVector& ComplexVector::ReverseVector()
{
  if (this->Dimension != 0)
    {
      this->Localize();
      int Max = this->Dimension - 1;
      Complex Tmp;
      for (int i = 0; i < this->Dimension/2; ++i)
	{
	  Tmp = this->Components[i];
	  this->Components[i] = this->Components[Max - i];
	  this->Components[Max - i] = Tmp;
	}
      this->Delocalize(true);
    }
  return *this;
}

  
// Output Stream overload
//

ostream& operator << (ostream& Str, const ComplexVector& P)
{
  for (int i = 0; i < P.Dimension; ++i)
    {
      Str << P.Components[i].Re;
      if (P.Components[i].Im < 0.0)
	Str << P.Components[i].Im << "i    ";
      else
	if (P.Components[i].Im != 0.0)
	  Str << "+" << P.Components[i].Im << "i    ";
	else
	  Str << "    ";
      Str << endl;
    }
  return Str;
}

#ifdef USE_OUTPUT
// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// v = vector to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, const ComplexVector& v)
{
  Str << "{";
  int i = 0;
  for (; i < (v.Dimension - 1); ++i)
    {
      Str << v.Components[i].Re;
      if (v.Components[i].Im < 0.0)
	Str << v.Components[i].Im << "I    ";
      else
	if (v.Components[i].Im != 0.0)
	  Str << "+" << v.Components[i].Im << "I,";
	else
	  Str << ",";
    }
  Str << v.Components[i++].Re;
  if (v.Components[i].Im < 0.0)
    Str << v.Components[i].Im << "I}";
  else
    if (v.Components[i].Im != 0.0)
      Str << "+" << v.Components[i].Im << "I}";
    else
      Str << "}";
  return Str;
}
#endif

// write vector in a file 
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool ComplexVector::WriteVector (const char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->Dimension);
  for (int i = 0; i < this->Dimension; ++i)
    {
      WriteLittleEndian(File, this->Components[i].Re);
      WriteLittleEndian(File, this->Components[i].Im);
    }
  File.close();
  return true;
}

// write vector in a file in ascii mode
//
// fileName = name of the file where the vector has to be stored
// return value = true if no error occurs

bool ComplexVector::WriteAsciiVector (const char* fileName)
{
  ofstream File;
  File.precision(14);
  File.open(fileName, ios::binary | ios::out);
  int ReducedDimension = this->Dimension - 1;
  for (int i = 0; i < ReducedDimension; ++i)
    File << this->Components[i].Re << " " << this->Components[i].Im << " ";
  File << this->Components[ReducedDimension].Re << " " << this->Components[ReducedDimension].Im << endl;  
  File.close();
  return true;
}

// read vector from a file 
//
// fileName = name of the file where the vector has to be read
// return value = true if no error occurs

bool ComplexVector::ReadVector (const char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "Cannot open the file: " << fileName << endl;
      return false;
    }

  unsigned ZeroPos, MaxPos;
  File.seekg (0, ios::beg);
  ZeroPos = File.tellg();
  File.seekg (0, ios::end);
  MaxPos = File.tellg ();
  
  unsigned Length = MaxPos-ZeroPos-sizeof(int);  
  File.seekg (0, ios::beg);
  int TmpDimension;
  ReadLittleEndian(File, TmpDimension);

  if (Length/sizeof(double)<(unsigned)TmpDimension)
    {      
      cout << "Error reading complex vector "<<fileName<<": estimated length "<<Length/(2*sizeof(double))<<" vs dimension "<<TmpDimension<<endl;
      if ((unsigned)TmpDimension==Length/sizeof(double))
	cout << "This could be a real vector!"<<endl;
      exit(1);
    }

  this->Resize(TmpDimension);
  for (int i = 0; i < this->Dimension; ++i)
    {
      ReadLittleEndian(File, this->Components[i].Re);
      ReadLittleEndian(File, this->Components[i].Im);
    }
  File.close();
  return true;
}

#ifdef __MPI__

// send a vector to a given MPI process
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& ComplexVector::SendVector(MPI::Intracomm& communicator, int id)
{
  communicator.Send(&this->VectorType, 1, MPI::INT, id, 1);
  communicator.Send(&this->Dimension, 1, MPI::INT, id, 1); 
  int Acknowledge = 0;
  communicator.Recv(&Acknowledge, 1, MPI::INT, id, 1);
  if (Acknowledge != 0)
    return *this;
  communicator.Send(this->Components, 2 * this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// broadcast a vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// return value = reference on the current vector

Vector& ComplexVector::BroadcastVector(MPI::Intracomm& communicator,  int id)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
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
  communicator.Bcast(this->Components, 2 * this->Dimension, MPI::DOUBLE, id); 
  return *this;
}

// broadcast part of vector to all MPI processes associated to the same communicator
// 
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// firstComponent = index of the first component (useless if the method is not called by the MPI process which broadcasts the vector)
// nbrComponent = number of component (useless if the method is not called by the MPI process which broadcasts the vector)
// return value = reference on the current vector

Vector& ComplexVector::BroadcastPartialVector(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  int TmpVectorType = this->VectorType;
  int TmpDimension = this->Dimension;
  int Acknowledge = 0;
  communicator.Bcast(&TmpVectorType, 1, MPI::INT, id);
  communicator.Bcast(&TmpDimension, 1, MPI::INT, id);
  communicator.Bcast(&firstComponent, 1, MPI::INT, id);
  communicator.Bcast(&nbrComponent, 1, MPI::INT, id);
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
  communicator.Bcast(this->Components + firstComponent, 2 * nbrComponent, MPI::DOUBLE, id); 
  return *this;
}

// receive a vector from a MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the source MPI process
// return value = reference on the current vector

Vector& ComplexVector::ReceiveVector(MPI::Intracomm& communicator, int id)
{
  int TmpVectorType = 0;
  int TmpDimension = 0;
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
  communicator.Recv(this->Components, 2 * this->Dimension, MPI::DOUBLE, id, 1); 
  return *this;
}

// add current vector to the current vector of a given MPI process
// 
// communicator = reference on the communicator to use 
// id = id of the destination MPI process
// return value = reference on the current vector

Vector& ComplexVector::SumVector(MPI::Intracomm& communicator, int id)
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
  communicator.Reduce(this->Components, TmpComponents, 2 * this->Dimension, MPI::DOUBLE, MPI::SUM, id); 
  if (id == communicator.Get_rank())
    {
      for (int i = 0; i < this->Dimension; ++i)
	this->Components[i] = TmpComponents[i];
      delete [] TmpComponents;
    }
  return *this;
}

// create a new vector on each MPI node which is an exact clone of the broadcasted one
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* ComplexVector::BroadcastClone(MPI::Intracomm& communicator, int id)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 2;
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
      communicator.Bcast(this->Components, 2 * this->Dimension, MPI::DOUBLE, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new ComplexVector(communicator, id);
    }
  return 0;
}

// create a new vector on given MPI node which is an exact clone of the sent one but with only part of the data
// 
// communicator = reference on the communicator to use
// id = id of the destination MPI process
// firstComponent = index of the first component 
// nbrComponent = number of component to send
// return value = reference on the current vector

Vector& ComplexVector::SendPartialClone(MPI::Intracomm& communicator, int id, int firstComponent, int nbrComponent)
{
  communicator.Send(&this->VectorType, 1, MPI::INT, id, 1);
  int TmpArray[5];
  TmpArray[0] = nbrComponent;
  TmpArray[1] = this->VectorId;
  TmpArray[2] = 2;
  TmpArray[3] = firstComponent;
  TmpArray[4] = this->Dimension;
  communicator.Send(TmpArray, 5, MPI::INT, id, 1); 
  communicator.Send(this->Components + firstComponent, 2 * nbrComponent, MPI::DOUBLE, id, 1); 
  return *this;
}

// create a new vector on each MPI node with same size and same type but non-initialized components
//
// communicator = reference on the communicator to use 
// id = id of the MPI process which broadcasts the vector
// zeroFlag = true if all coordinates have to be set to zero
// return value = pointer to new vector 

Vector* ComplexVector::BroadcastEmptyCloneVector(MPI::Intracomm& communicator, int id, bool zeroFlag)
{
  if (id == communicator.Get_rank())
    {
      communicator.Bcast(&this->VectorType, 1, MPI::INT, id);
      int TmpArray[3];
      TmpArray[0] = this->Dimension;
      TmpArray[1] = this->VectorId;
      TmpArray[2] = 0;
      if (zeroFlag == true)
	{
	  TmpArray[2] = 1;
	}
      communicator.Bcast(TmpArray, 3, MPI::INT, id);      
    }
  else
    {
      int Type = 0;
      communicator.Bcast(&Type, 1, MPI::INT, id);  
      return new ComplexVector(communicator, id);
    }
  return 0;
}

#endif
