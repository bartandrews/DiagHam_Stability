////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract hamiltonian                      //
//                                                                            //
//                        last modification : 28/02/2001                      //
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


#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Complex.h"
#include "BitmapPicture/AbstractBitmapPicture.h"
#include "BitmapPicture/TgaFormat.h" 
#include "Color/PicRGB.h"


#include <stdlib.h>


using std::cout;
using std::endl;


// destructor
//

AbstractHamiltonian::~AbstractHamiltonian() 
{
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix
HermitianMatrix& AbstractHamiltonian::GetHamiltonian (HermitianMatrix& M)
{
  ComplexVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = Complex(1.0, 0.0);
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = Complex(0.0, 0.0);	
    }
  return M;
}
  
// store real part of Hamiltonian into a real symmetric matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding real symmetric matrix 

RealSymmetricMatrix& AbstractHamiltonian::GetHamiltonian (RealSymmetricMatrix& M)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	M.SetMatrixElement(i, j, TmpV2[j]);
      TmpV1[i] = 0.0;	
    }
  return M;
}
  
// store Hamiltonian into a matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding matrix 

Matrix& AbstractHamiltonian::GetHamiltonian (Matrix& M)
{
  switch (M.GetMatrixType())
    {
      case (Matrix::RealElements | Matrix::Symmetric):
	return this->GetHamiltonian((RealSymmetricMatrix&) M);
      break;
      case (Matrix::ComplexElements | Matrix::Hermitian):
	return this->GetHamiltonian((HermitianMatrix&) M);
      break;
    default:
      return M;
      break;
    }
  return M;
}
  
// return matrix representation of current Hamiltonian
//
// return value = reference to representation

Matrix* AbstractHamiltonian::GetHamiltonian ()
{
  HermitianMatrix* TmpH = new HermitianMatrix(this->GetHilbertSpaceDimension());
  this->GetHamiltonian(*TmpH);
  return TmpH;
}
  
// store Hamiltonian into a picture (drawing non zero element with a color scale)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color BlackColor (0.0, 0.0, 0.0);
  Color WhiteColor (1.0, 1.0, 1.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, WhiteColor);
	    }
	  else
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// store Hamiltonian into a picture (drawing non zero element in black)
//
// error = absolute minimum value to be considered as non zero element
// return value = pointer to the picture associated to the matrix

AbstractBitmapPicture* AbstractHamiltonian::GetHamiltonianColorPicture (double error)
{
  RealVector TmpV1 (this->GetHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetHilbertSpaceDimension(), true);
  Color RedColor (1.0, 0.0, 0.0);
  Color GreenColor (0.0, 1.0, 0.0);
  Color BlueColor (0.0, 0.0, 1.0);
  Color BlackColor (0.0, 0.0, 0.0);
  TgaFormat* TmpPicture = new TgaFormat (this->GetHilbertSpaceDimension(), this->GetHilbertSpaceDimension());
  double Max = 0.0;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) > error)
	    {
	      if (fabs(TmpV2[j]) > Max)
		{
		  Max = fabs(TmpV2[j]);
		}
	    }
	}
      TmpV1[i] = 0.0;	
    }
  Max = 4.0 / Max;
  for (int i = 0; i < this->GetHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiply(TmpV1, TmpV2);
      for (int j = 0; j < this->GetHilbertSpaceDimension(); j++)
	{
	  if (fabs(TmpV2[j]) < error)
	    {
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, BlackColor);
	    }
	  else
	    {
	      double Fac = (fabs(TmpV2[j]) * Max);
	      Color TmpColor;
	      if (Fac >= 4.0)
		{
		  TmpColor = RedColor;
		}
	      else
		if (Fac >= 3.0)
		  {
		    TmpColor = GreenColor * (4.0 - Fac) + RedColor;
		  }
		else
		  if (Fac >= 2.0)
		    {
		      TmpColor = GreenColor + RedColor * (Fac - 2.0);
		    }
		  else
		    if (Fac >= 1.0)
		      {
			TmpColor = BlueColor * (2.0 - Fac) + GreenColor;
		      }
		    else
		      {
			TmpColor = GreenColor * (Fac - 1.0) + BlueColor;
		      }
	      ((AbstractBitmapPicture*) TmpPicture)->SetPixel(i, j, TmpColor);
	    }
	}
      TmpV1[i] = 0.0;	
    }
  return TmpPicture;
}
  
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractHamiltonian::LeftInteractionOperators()
{
  return List<Matrix*>();
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractHamiltonian::RightInteractionOperators()
{
  return List<Matrix*>();
}

// multiply a vector by the current hamiltonian and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination)
{
  if (vSource.GetVectorType() != vDestination.GetVectorType())
    {
      cout << "error" << endl;
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::Multiply(Vector& vSource, Vector& vDestination, 
				      int firstComponent, int nbrComponent)
{
  if (vSource.GetVectorType() != vDestination.GetVectorType())
    return vDestination;
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination)
{
  if (vSource.GetVectorType() != vDestination.GetVectorType())
    {
      cout << "error" << endl;
      return vDestination;
    }
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination);
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

Vector& AbstractHamiltonian::AddMultiply(Vector& vSource, Vector& vDestination, 
					 int firstComponent, int nbrComponent)
{
  if (vSource.GetVectorType() != vDestination.GetVectorType())
    return vDestination;
  if (vSource.GetVectorType() == Vector::RealDatas)
    {
      return this->LowLevelAddMultiply((RealVector&) vSource, (RealVector&) vDestination, firstComponent, nbrComponent);
    }
  else
    {
      return this->LowLevelAddMultiply((ComplexVector&) vSource, (ComplexVector&) vDestination, firstComponent, nbrComponent);
    }
  return vDestination;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on real vector used as first vector (will contain last produced vector at the end)
// return value = reference on real tridiagonal symmetric  matrix

RealTriDiagonalSymmetricMatrix& AbstractHamiltonian::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
							      ComplexVector& V1) 
{
  ComplexVector V2(V1.GetVectorDimension());
  ComplexVector V3(V1.GetVectorDimension());
  int Index = 0;
  if (M.GetNbrRow() != dimension)
    M.Resize(dimension, dimension);
  V1 /= V1.Norm();
  this->LowLevelMultiply(V1, V2);
  M.DiagonalElement(Index) = (V1 * V2).Re;
  V2.AddLinearCombination(-M.DiagonalElement(Index), V1);
  V2 /= V2.Norm();
  for (int i = 2; i < dimension; i++)
    {
      this->LowLevelMultiply(V2, V3);
      M.UpperDiagonalElement(Index) = (V1 * V3).Re;
      M.DiagonalElement(Index + 1) = (V2 * V3).Re;
      V3.AddLinearCombination(-M.DiagonalElement(Index + 1), V2);
      V3.AddLinearCombination(-M.UpperDiagonalElement(Index), V1);
      V3 /= V3.Norm();
      ComplexVector TmpV (V1);
      V1 = V2;
      V2 = V3;
      V3 = TmpV;
      Index++;
    }  
  M.UpperDiagonalElement(Index) = this->MatrixElement(V1, V2).Re;
  M.DiagonalElement(Index + 1) = this->MatrixElement(V2, V2).Re;
  return M;
}

// Tridiagonalize hamiltonian using Lanczos algorithm with base re-orthogonalization
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on real vector used as first vector (will contain last produced vector at the end)
// return value = reference on real tridiagonal symmetric  matrix

RealTriDiagonalSymmetricMatrix& AbstractHamiltonian::FullReorthogonalizedLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
										  RealVector& V1) 
{
  RealVector* Vectors = new RealVector [dimension];
  Vectors[0] = V1;
  Vectors[1] = RealVector (V1.GetVectorDimension());
  int Index = 0;
  if (M.GetNbrRow() != dimension)
    M.Resize(dimension, dimension);
  Vectors[0] /= Vectors[0].Norm();
  this->LowLevelMultiply(Vectors[0], Vectors[1]);
  M.DiagonalElement(Index) = (Vectors[0] * Vectors[1]);

  Vectors[1].AddLinearCombination(-M.DiagonalElement(Index), Vectors[0]);
  Vectors[1] /= Vectors[1].Norm();
  for (int i = 2; i < dimension; i++)
    {
      Vectors[i] = RealVector (V1.GetVectorDimension());
      this->LowLevelMultiply(Vectors[i - 1], Vectors[i]);
      M.UpperDiagonalElement(Index) = (Vectors[i - 2] * Vectors[i]);
      M.DiagonalElement(Index + 1) = (Vectors[i - 1] * Vectors[i]);
      Vectors[i].AddLinearCombination(-M.DiagonalElement(Index + 1), Vectors[i - 1], 
				      -M.UpperDiagonalElement(Index), Vectors[i - 2]);
      for (int j = 0; j < i; j++)
	Vectors[i].AddLinearCombination(- (Vectors[j] * Vectors[i]), Vectors[j]);	    
      double VectorNorm = Vectors[i].Norm();
      while (VectorNorm < 1e-5)
	{
	  cout << "subspace !!! " << i << endl;
	  double tmp = 0;
	  for (int j = 0; j < V1.GetVectorDimension(); j++)
	    {
	      Vectors[i][j] = rand ();
	      tmp += Vectors[i][j] * Vectors[i][j];
	    }
	  tmp = sqrt(tmp);
	  Vectors[i] /= tmp;
	  RealVector TmpVector(V1.GetVectorDimension());
	  this->LowLevelMultiply(Vectors[i], TmpVector);
	  Vectors[i] = TmpVector;
	  for (int j = 0; j < i; j++)
	    {
	      Vectors[i].AddLinearCombination(- (Vectors[j] * Vectors[i]), Vectors[j]);
	    }
	  VectorNorm = Vectors[i].Norm();
	}
      Vectors[i] /= VectorNorm;
//      cout << (Vectors[0] * Vectors[i]) << " | ";
      Index++;
    }  
  M.UpperDiagonalElement(Index) = this->MatrixElement(Vectors[V1.GetVectorDimension() - 2], 
						      Vectors[V1.GetVectorDimension() - 1]).Re;
  M.DiagonalElement(Index + 1) = this->MatrixElement(Vectors[V1.GetVectorDimension() - 1], 
						     Vectors[V1.GetVectorDimension() - 1]).Re;
  delete[] Vectors;
  return M;
}

// Tridiagonalize hamiltonian using Lanczos algorithm with partial base re-orthogonalization
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric complex tridiagonal matrix where result has to be stored
// V1 = reference on real vector used as first vector (will contain last produced vector at the end)
// step = number of iterations before re-orthogonalizing whole base
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& AbstractHamiltonian::ReorthogonalizedLanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
									      RealVector& V1, int step) 
{
  RealVector* Vectors = new RealVector [dimension];
  Vectors[0] = V1;
  Vectors[1] = RealVector (V1.GetVectorDimension());
  int Index = 0;
  if (M.GetNbrRow() != dimension)
    M.Resize(dimension, dimension);
  Vectors[0] /= Vectors[0].Norm();
  this->LowLevelMultiply(Vectors[0], Vectors[1]);
  M.DiagonalElement(Index) = (Vectors[0] * Vectors[1]);

  Vectors[1].AddLinearCombination(-M.DiagonalElement(Index), Vectors[0]);
  Vectors[1] /= Vectors[1].Norm();
  int Step = 0;
  for (int i = 2; i < dimension; i++)
    {
      Vectors[i] = RealVector (V1.GetVectorDimension());
      this->LowLevelMultiply(Vectors[i - 1], Vectors[i]);
      M.UpperDiagonalElement(Index) = (Vectors[i - 2] * Vectors[i]);
      M.DiagonalElement(Index + 1) = (Vectors[i - 1] * Vectors[i]);
      Vectors[i].AddLinearCombination(-M.DiagonalElement(Index + 1), Vectors[i - 1], 
				      -M.UpperDiagonalElement(Index), Vectors[i - 2]);
      if (++Step == step)
	{
	  for (int j = 0; j < (i - 1); j++)
	    Vectors[i - 1].AddLinearCombination(- (Vectors[j] * Vectors[i - 1]), Vectors[j]);	    
	  for (int j = 0; j < i; j++)
	    Vectors[i].AddLinearCombination(- (Vectors[j] * Vectors[i]), Vectors[j]);	    
	  Step = 0;
	}
      double VectorNorm = Vectors[i].Norm();
      while (VectorNorm < 1e-5)
	{
	  cout << "subspace !!!" << endl;
	  double tmp = 0;
	  for (int j = 0; j < V1.GetVectorDimension(); j++)
	    {
	      Vectors[i][j] = rand ();
	      tmp += Vectors[i][j] * Vectors[i][j];
	    }
	  tmp = sqrt(tmp);
	  Vectors[i] /= tmp;
	  RealVector TmpVector(V1.GetVectorDimension());
	  this->LowLevelMultiply(Vectors[i], TmpVector);
	  Vectors[i] = TmpVector;
	  for (int j = 0; j < i; j++)
	    {
	      Vectors[i].AddLinearCombination(- (Vectors[j] * Vectors[i]), Vectors[j]);
	    }
	  VectorNorm = Vectors[i].Norm();
	}
      Vectors[i] /= VectorNorm;
//      cout << (Vectors[0] * Vectors[i]) << " | ";
      Index++;
    }  
  M.UpperDiagonalElement(Index) = this->MatrixElement(Vectors[V1.GetVectorDimension() - 2], 
						      Vectors[V1.GetVectorDimension() - 1]).Re;
  M.DiagonalElement(Index + 1) = this->MatrixElement(Vectors[V1.GetVectorDimension() - 1], 
						     Vectors[V1.GetVectorDimension() - 1]).Re;
  delete[] Vectors;
  return M;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on real vector used as first vector (will contain last produced vector at the end)
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& AbstractHamiltonian::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
							      RealVector& V1) 
{
  RealVector V2(V1.GetVectorDimension());
  RealVector V3(V1.GetVectorDimension());
  int Index = 0;
  if (M.GetNbrRow() != dimension)
    M.Resize(dimension, dimension);
  V1 /= V1.Norm();
  this->LowLevelMultiply(V1, V2);
  M.DiagonalElement(Index) = (V1 * V2);
  V2.AddLinearCombination(-M.DiagonalElement(Index), V1);
  double tmpNorm;
  V2 /= V2.Norm(); 
  for (int i = 2; i < dimension; i++)
    {
      this->LowLevelMultiply(V2, V3);
      M.UpperDiagonalElement(Index) = (V1 * V3);
      M.DiagonalElement(Index + 1) = (V2 * V3);
      V3.AddLinearCombination(-M.DiagonalElement(Index + 1), V2,-M.UpperDiagonalElement(Index), V1);
      tmpNorm = V3.Norm();
      V3 /= tmpNorm;
      RealVector TmpV (V1);
      V1 = V2;
      V2 = V3;
      V3 = TmpV;
      Index++;
    }  
  M.UpperDiagonalElement(Index) = this->MatrixElement(V1, V2).Re;
  M.DiagonalElement(Index + 1) = this->MatrixElement(V2, V2).Re;
  return M;
}

// Tridiagonalize an hermitian matrix using Lanczos algorithm without re-orthogonalizing base at each step
//
// dimension = maximum iteration number
// M = reference on real tridiagonal symmetric matrix where result has to be stored
// V1 = reference on real vector used as first vector (will contain last produced vector at the end)
// V2 = reference on real vector used as second vector (will contain next to last produced vector at the end)
// useV2 = true if V2 already contains result of second Lanczos iteration, in that case M is supposed to give
//         results of previous Lanczos iteration
// return value = reference on real tridiagonal symmetric matrix

RealTriDiagonalSymmetricMatrix& AbstractHamiltonian::Lanczos (int dimension, RealTriDiagonalSymmetricMatrix& M, 
							      RealVector& V1, RealVector& V2, bool useV2) 
{
  RealVector V3(V1.GetVectorDimension());
  int Index;
  int Start;
  double tmpNorm;
  if (useV2 == false)
    {
      Start = 2;
      Index = 0;
      if (M.GetNbrRow() != dimension)
	M.Resize(dimension, dimension);
      V1 /= V1.Norm();
      this->LowLevelMultiply(V1, V2);
      M.DiagonalElement(Index) = (V1 * V2);
      V2.AddLinearCombination(-M.DiagonalElement(Index), V1);
      V2 /= V2.Norm(); 
    }
  else 
    {
      Start = M.GetNbrRow();
      Index = M.GetNbrRow() - 1;
      if (M.GetNbrRow() != dimension)
	M.Resize(dimension, dimension);
    }
  for (int i = Start; i < dimension; i++)
    {
      this->LowLevelMultiply(V2, V3);
      M.UpperDiagonalElement(Index) = (V1 * V3);
      M.DiagonalElement(Index + 1) = (V2 * V3);
      V3.AddLinearCombination(-M.DiagonalElement(Index + 1), V2,-M.UpperDiagonalElement(Index), V1);
      tmpNorm = V3.Norm();
      V3 /= tmpNorm;
      RealVector TmpV (V1);
      V1 = V2;
      V2 = V3;
      V3 = TmpV;
      Index++;
    }  
  M.UpperDiagonalElement(Index) = this->MatrixElement(V1, V2).Re;
  M.DiagonalElement(Index + 1) = this->MatrixElement(V2, V2).Re;
  return M;
}

// multiply a vector by the current hamiltonian using threads
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// nbrProcess = number of process to run
// return value = reference on vector where result has been stored

#ifdef __SMP__

RealVector& AbstractHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int nbrProcess)
{
  return this->LowLevelMultiply(vSource, vDestination);
}

// multiply a vector by the current hamiltonian using threads
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// nbrProcess = number of process to run
// return value = reference on vector where result has been stored

ComplexVector& AbstractHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int nbrProcess)
{
  return this->LowLevelMultiply(vSource, vDestination);
}

#endif
