////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of fermions on cylinder that allow to use MPS with operator     //
//                                                                            //
//                        last modification : 15/10/2012                      //
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
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>
#include <stdlib.h>
#include <fstream>

//#define MPSWRAPPER_MULT
#define MPSWRAPPER_CONJ
#define MPSWRAPPER_SMP

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constuctor
//

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// referenceState = array that describes the root configuration
// rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// columnIndex = column index of the MPS element that has to be evaluated
// bMatrices = array that gives the B matrices 
// memory = amount of memory granted for precalculations

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,  
							  int rowIndex, int columnIndex, SparseRealMatrix* bMatrices, AbstractArchitecture* architecture,
							  unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalLz = 0;
  this->StateDescription = 0x0ul;
  int TmpIndex = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->StateDescription |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  this->TotalLz += i;
	  this->StateLzMax = i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
  this->LargeHilbertSpaceDimension = 1l;
  this->HilbertSpaceDimension = 1;
  this->Flag.Initialize();
  this->StateDescription = 0l;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->MPSRowIndex = rowIndex;
  this->MPSColumnIndex = columnIndex;
  this->Architecture = architecture;

  int NbrBMatrices = 2;
  this->BMatrices = new SparseRealMatrix[NbrBMatrices];
  this->ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      this->BMatrices[i] = bMatrices[i];
      this->ConjugateBMatrices[i] = bMatrices[i].Transpose();
    }

  this->MaxTmpMatrixElements = (((long) this->BMatrices[0].GetNbrRow()) * 
				((long) this->BMatrices[0].GetNbrRow() / 1l));
  cout << "Requested memory for sparse matrix multiplications = " << ((this->MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  this->TmpMatrixElements = new double [this->MaxTmpMatrixElements];
  this->TmpColumnIndices = new int [this->MaxTmpMatrixElements];
  this->TmpElements = new double [this->BMatrices[0].GetNbrRow()];
  SparseRealMatrix TmpMatrixNorm (this->ComputeBMatrixNormalization());
  double Tmp;
  TmpMatrixNorm.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  this->StateNormalization = Tmp;
  cout<<"Cylinder normalization "<<this->StateNormalization<<endl;
}

// constructor in presence of quasiholes
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// referenceState = array that describes the root configuration
// rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// columnIndex = column index of the MPS element that has to be evaluated
// bMatrices = array that gives the B matrices 
// quasiholeBMatrices = array that contains the quasihole B matrices
// nbrQuasiholes = number of quasihole B matrices
// architecture = pointer to the archiecture
// memory = amount of memory granted for precalculations

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,  
							  int rowIndex, int columnIndex, SparseRealMatrix* bMatrices, 
							  SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholes,
							  AbstractArchitecture* architecture, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalLz = 0;
  this->StateDescription = 0x0ul;
  int TmpIndex = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->StateDescription |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  this->TotalLz += i;
	  this->StateLzMax = i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
  this->LargeHilbertSpaceDimension = 1l;
  this->HilbertSpaceDimension = 1;
  this->Flag.Initialize();
  this->StateDescription = 0l;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->MPSRowIndex = rowIndex;
  this->MPSColumnIndex = columnIndex;
  this->Architecture = architecture;

  int NbrBMatrices = 2;
  this->BMatrices = new SparseRealMatrix[NbrBMatrices];
  this->ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      this->BMatrices[i] = bMatrices[i];
      this->ConjugateBMatrices[i] = bMatrices[i].Transpose();
    }

  this->MaxTmpMatrixElements = (((long) this->BMatrices[0].GetNbrRow()) * 
				((long) this->BMatrices[0].GetNbrRow() / 1l));
  cout << "Requested memory for sparse matrix multiplications = " << ((this->MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  this->TmpMatrixElements = new double [this->MaxTmpMatrixElements];
  this->TmpColumnIndices = new int [this->MaxTmpMatrixElements];
  this->TmpElements = new double [this->BMatrices[0].GetNbrRow()];
  this->QuasiholeBMatrixContribution.Copy(quasiholeBMatrices[0]);
  for (int i = 1; i < nbrQuasiholes; ++i)
    {
      this->QuasiholeBMatrixContribution.Multiply(quasiholeBMatrices[i]);
    }
  this->ConjugateQuasiholeBMatrixContribution = this->QuasiholeBMatrixContribution.HermitianTranspose();
  SparseRealMatrix TmpMatrixNorm (this->ComputeBMatrixNormalization());
  SparseComplexMatrix TmpMatrixNorm2 = Conjugate(this->ConjugateQuasiholeBMatrixContribution, TmpMatrixNorm, this->QuasiholeBMatrixContribution);
  Complex Tmp;
  TmpMatrixNorm2.PrintNonZero(cout) << endl;
  TmpMatrixNorm2.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  this->ComplexStateNormalization = Tmp;
  cout << "Cylinder normalization " << this->ComplexStateNormalization << endl;
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper(const FermionOnCylinderMPSWrapper& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->MPSRowIndex = fermions.MPSRowIndex;
  this->MPSColumnIndex = fermions.MPSColumnIndex;
  this->BMatrices = fermions.BMatrices;
  this->ConjugateBMatrices = fermions.ConjugateBMatrices;
  this->StateNormalization = fermions.StateNormalization;
  this->MaxTmpMatrixElements = fermions.MaxTmpMatrixElements;
  this->TmpMatrixElements = new double [this->MaxTmpMatrixElements];
  this->TmpColumnIndices = new int [this->MaxTmpMatrixElements];
  this->TmpElements = new double [this->BMatrices[0].GetNbrRow()];
  this->Architecture = fermions.Architecture;
}

// destructor
//

FermionOnCylinderMPSWrapper::~FermionOnCylinderMPSWrapper ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnCylinderMPSWrapper& FermionOnCylinderMPSWrapper::operator = (const FermionOnCylinderMPSWrapper& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->BMatrices;
      delete[] this->ConjugateBMatrices;
    }
  delete[] this->TmpMatrixElements;
  delete[] this->TmpColumnIndices;
  delete[] this->TmpElements;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->MPSRowIndex = fermions.MPSRowIndex;
  this->MPSColumnIndex = fermions.MPSColumnIndex;
  this->StateNormalization = fermions.StateNormalization;
  this->MaxTmpMatrixElements = fermions.MaxTmpMatrixElements;
  this->TmpMatrixElements = new double [this->MaxTmpMatrixElements];
  this->TmpColumnIndices = new int [this->MaxTmpMatrixElements];
  this->TmpElements = new double [this->BMatrices[0].GetNbrRow()];
  this->Architecture = fermions.Architecture;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnCylinderMPSWrapper::Clone()
{
  return new FermionOnCylinderMPSWrapper(*this);
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnCylinderMPSWrapper::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  if ((m1 == m2) || (n1 == n2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  SparseRealMatrix TmpMatrix (this->BMatrices[0].GetNbrRow(), this->BMatrices[0].GetNbrRow());
  TmpMatrix.SetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, 1.0);

  double Sign = 1.0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      if (i == m1)
	{
	  if ((m1 == n1) || (m1 == n2))
	    {
	      SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[1],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[1]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif

	      Sign = 1.0;
	      if (m1 > i)
		Sign *= -1.0;
	      if (m2 > i)
		Sign *= -1.0;
	      if (n1 > i)
		Sign *= -1.0;
	      if (n2 > i)
		Sign *= -1.0; 
	      TmpMatrix2 *= Sign;
	      TmpMatrix = TmpMatrix2;
	    }
	  else
	    {
	      SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[0],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[0]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif

	      Sign = 1.0;
	      if (m1 > i)
		Sign *= -1.0;
	      if (m2 > i)
		Sign *= -1.0;
	      TmpMatrix2 *= Sign;
	      TmpMatrix = TmpMatrix2;
	    }
	}
      else
	{
	  if (i == m2)
	    {
	      if ((m2 == n1) || (m2 == n2))
		{
		  SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
		  if (this->Architecture == 0)
		    {
#endif
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
		      TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[1],  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		    }
		  else
		    {
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		      TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[1]),  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		    }
#endif

		  Sign = 1.0;
		  if (m1 > i)
		    Sign *= -1.0;
		  if (m2 > i)
		    Sign *= -1.0;
		  if (n1 > i)
		    Sign *= -1.0;
		  if (n2 > i)
		    Sign *= -1.0; 
		  TmpMatrix2 *= Sign;
		  TmpMatrix = TmpMatrix2;
		}
	      else
		{
		  SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
		  if (this->Architecture == 0)
		    {
#endif
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
		      TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[0],  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		    }
		  else
		    {
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		      
#ifdef MPSWRAPPER_CONJ	      
		      TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[0]),  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		    }
#endif
		  Sign = 1.0;
		  if (m1 > i)
		    Sign *= -1.0;
		  if (m2 > i)
		    Sign *= -1.0;
		  TmpMatrix2 *= Sign;
		  TmpMatrix = TmpMatrix2;
		}
	    }
	  else
	    {
	      if ((i == n1) || (i == n2))
		{
		  SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
		  if (this->Architecture == 0)
		    {
#endif
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif		      
#ifdef MPSWRAPPER_CONJ
		      TmpMatrix2 = Conjugate(this->ConjugateBMatrices[0], TmpMatrix, this->BMatrices[1],  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		    }
		  else
		    {
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[0]), &TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		      TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[0]), &TmpMatrix, &(this->BMatrices[1]),  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		    }
#endif

		  Sign = 1.0;
		  if (n1 > i)
		    Sign *= -1.0;
		  if (n2 > i)
		    Sign *= -1.0; 
		  TmpMatrix2 *= Sign;
		  TmpMatrix = TmpMatrix2;
		}
	      else
		{
		  SparseRealMatrix TmpMatrix2;
		  SparseRealMatrix TmpMatrix3;
#ifdef MPSWRAPPER_SMP
		  if (this->Architecture == 0)
		    {
#endif
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix3 = Multiply(TmpMatrix1, this->BMatrices[1],  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
		      TmpMatrix2 = Conjugate(this->ConjugateBMatrices[0], TmpMatrix, this->BMatrices[0],  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		      TmpMatrix3 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[1],  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		    }
		  else
		    {		      
#ifdef MPSWRAPPER_MULT
		      SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[0]), &TmpMatrix, 
							     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix3 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					    this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		      TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[0]), &TmpMatrix, &(this->BMatrices[0]),  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		      TmpMatrix3 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[1]),  
					     this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		    }
#endif

		  Sign = 1.0;
		  if (m1 > i)
		    Sign *= -1.0;
		  if (m2 > i)
		    Sign *= -1.0;
		  if (n1 > i)
		    Sign *= -1.0;
		  if (n2 > i)
		    Sign *= -1.0; 
		  TmpMatrix3 *= Sign;
		  TmpMatrix = TmpMatrix2 + TmpMatrix3;
		}
	    }
	}
    }
  double Tmp = 0.0;
  TmpMatrix.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  coefficient = -Tmp / this->StateNormalization;
  return 0;

}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnCylinderMPSWrapper::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  cout<<"Unimplemented"<<endl;
  return this->HilbertSpaceDimension;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnCylinderMPSWrapper::ProdA (int index, int* n, int nbrIndices)
{
  cout<<"ProdA not implemented"<<endl;
  return 0.0;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnCylinderMPSWrapper::AA (int index, int n1, int n2)
{
  cout<<"AA not implemented"<<endl;
  return 0.0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnCylinderMPSWrapper::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  cout<<"ProdAd not implemented"<<endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnCylinderMPSWrapper::AdAd (int m1, int m2, double& coefficient)
{
  cout<<"AdAd not implemented"<<endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnCylinderMPSWrapper::AdA (int index, int m)
{
   SparseRealMatrix TmpMatrix (this->BMatrices[0].GetNbrRow(), this->BMatrices[0].GetNbrRow());
   TmpMatrix.SetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, 1.0);

   for (int i = 0; i <= this->LzMax; ++i)
     {
       if (i == m)
	 {
	   SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
						  this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   SparseRealMatrix TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
						  this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   TmpMatrix = TmpMatrix2;
	 }
       else
	 {
	   SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrix, 
						  this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   SparseRealMatrix TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
						  this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
				 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   SparseRealMatrix TmpMatrix3 = Multiply(TmpMatrix1, this->BMatrices[1],  
						  this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	   TmpMatrix = TmpMatrix2 + TmpMatrix3;
	 }      
     }
  double Tmp = 0.0;
  TmpMatrix.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  return Tmp / this->StateNormalization;
}


// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnCylinderMPSWrapper::AdA (int index, int m, int n, double& coefficient)
{
  SparseRealMatrix TmpMatrix (this->AdACore(m ,n));
  double Tmp = 0.0;
  TmpMatrix.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  coefficient = Tmp / this->StateNormalization;
  return 0;
}

// core part of the  a^+_m a_n operator calculation
//
// m = index of the creation operator
// n = index of the annihilation operator
// return value = matrix that results from the a^+_m a_n operator calculation

SparseRealMatrix FermionOnCylinderMPSWrapper::AdACore (int m, int n)
{
  SparseRealMatrix TmpMatrix (this->BMatrices[0].GetNbrRow(), this->BMatrices[0].GetNbrRow());
  TmpMatrix.SetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, 1.0);
  double Sign = 1.0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      if (i == m)
	{
	  if (m == n)
	    {
	      SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[1],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[1]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif

	      Sign = 1.0;
	      if (m > i)
		Sign *= -1.0;
	      if (n > i)
		Sign *= -1.0;
	      TmpMatrix2 *= Sign;
	      TmpMatrix = TmpMatrix2;
	    }
	  else
	    {
	      SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[0],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{		  
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif		 
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[0]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif

	      Sign = 1.0;
	      if (m > i)
		Sign *= -1.0;
	      TmpMatrix2 *= Sign;
	      TmpMatrix = TmpMatrix2;
	    }
	}
      else
	{
	  if (i == n)
	    {
	      SparseRealMatrix TmpMatrix2;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[1],  
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[0], TmpMatrix, this->BMatrices[1],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{	      	      
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[0]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif	      
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[0]), &TmpMatrix, &(this->BMatrices[1]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif
	      
	      Sign = 1.0;
	      if (n > i)
		Sign *= -1.0;
	      TmpMatrix2 *= Sign;
	      TmpMatrix = TmpMatrix2;
	    }
	  else
	    {
	      SparseRealMatrix TmpMatrix2;
	      SparseRealMatrix TmpMatrix3;
#ifdef MPSWRAPPER_SMP
	      if (this->Architecture == 0)
		{
#endif
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
				       this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrix, 
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix3 = Multiply(TmpMatrix1, this->BMatrices[1],  
					this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
		  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[0], TmpMatrix, this->BMatrices[0],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
		  TmpMatrix3 = Conjugate(this->ConjugateBMatrices[1], TmpMatrix, this->BMatrices[1],  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
		}
	      else
		{
#ifdef MPSWRAPPER_MULT
		  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[0]), &TmpMatrix, 
							 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrix, 
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix3 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
					this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif	      	      
#ifdef MPSWRAPPER_CONJ	      
		  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[0]), &TmpMatrix, &(this->BMatrices[0]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
		  TmpMatrix3 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrix, &(this->BMatrices[1]),  
					 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
		}
#endif
	      Sign = 1.0;
	      if (m > i)
		Sign *= -1.0;
	      if (n > i)
		Sign *= -1.0;
	      TmpMatrix3 *= Sign;
	      TmpMatrix = TmpMatrix2 + TmpMatrix3;
	    }
	}
    }
  return TmpMatrix;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

long FermionOnCylinderMPSWrapper::AdA (long index, int m, int n, Complex& coefficient)
{
  SparseRealMatrix TmpMatrix (this->AdACore(m ,n));
  SparseComplexMatrix TmpMatrix2 = Conjugate(this->ConjugateQuasiholeBMatrixContribution, TmpMatrix, this->QuasiholeBMatrixContribution);
  Complex Tmp;
  TmpMatrix2.GetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, Tmp);
  coefficient = Tmp / this->ComplexStateNormalization;
  return 0;
}
  
// compute the B matrix contribution to the state normalization
//
// return value = matrix that provides the B matrix contribution

SparseRealMatrix FermionOnCylinderMPSWrapper::ComputeBMatrixNormalization()
{
  SparseRealMatrix TmpMatrixNorm (this->BMatrices[0].GetNbrRow(), this->BMatrices[0].GetNbrRow());
  TmpMatrixNorm.SetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, 1.0);

  for (int i = 0; i <= this->LzMax; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
#ifdef MPSWRAPPER_SMP
      if (this->Architecture == 0)
	{
#endif
#ifdef MPSWRAPPER_MULT
	  SparseRealMatrix TmpMatrix1 = Multiply(this->ConjugateBMatrices[0], TmpMatrixNorm, 
						 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	  TmpMatrix2 = Multiply(TmpMatrix1, this->BMatrices[0],  
						 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	  TmpMatrix1 = Multiply(this->ConjugateBMatrices[1], TmpMatrixNorm, 
				this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	  TmpMatrix3 = Multiply(TmpMatrix1, this->BMatrices[1],  
						 this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_CONJ
	  TmpMatrix2 = Conjugate(this->ConjugateBMatrices[0], TmpMatrixNorm, this->BMatrices[0], 
						this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
	  TmpMatrix3 = Conjugate(this->ConjugateBMatrices[1], TmpMatrixNorm, this->BMatrices[1], 
					      this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements); 
#endif
#ifdef MPSWRAPPER_SMP
	}
      else
	{
#ifdef MPSWRAPPER_MULT
	  SparseRealMatrix TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[0]), &TmpMatrixNorm, 
						 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
	  TmpMatrix2 = Multiply(&TmpMatrix1, &(this->BMatrices[0]),  
						 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
	  TmpMatrix1 = Multiply(&(this->ConjugateBMatrices[1]), &TmpMatrixNorm, 
				this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
	  TmpMatrix3 = Multiply(&TmpMatrix1, &(this->BMatrices[1]),  
						 this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
#ifdef MPSWRAPPER_CONJ
	  TmpMatrix2 = Conjugate(&(this->ConjugateBMatrices[0]), &TmpMatrixNorm, &(this->BMatrices[0]), 
						  this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
	  TmpMatrix3 = Conjugate(&(this->ConjugateBMatrices[1]), &TmpMatrixNorm, &(this->BMatrices[1]), 
						  this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture); 
#endif
	}
#endif  
       TmpMatrixNorm = TmpMatrix2 + TmpMatrix3;
     }
  return TmpMatrixNorm;
}
