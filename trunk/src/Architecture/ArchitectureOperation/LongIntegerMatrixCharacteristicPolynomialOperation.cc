////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of operation to evaluate the characteristic polynomial       //
//                           of a long integer matrix                         //
//                                                                            //
//                        last modification : 13/01/2022                      //
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
#include "Architecture/ArchitectureOperation/LongIntegerMatrixCharacteristicPolynomialOperation.h"
#include "Matrix/SparseRealMatrix.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"


// constructor 
//
// sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated

LongIntegerMatrixCharacteristicPolynomialOperation::LongIntegerMatrixCharacteristicPolynomialOperation (LongIntegerMatrix* sourceMatrix)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceMatrix->GetNbrRow();
  this->SourceMatrix = sourceMatrix;
  this->NbrNonZeroMatrixElements = 0;
  this->NonZeroMatrixElementPositions = 0;
#ifdef __GMP__
  this->CharacteristicPolynomial = new mpz_t[sourceMatrix->GetNbrRow() + 1];
#else
  this->CharacteristicPolynomial = new LONGLONG[sourceMatrix->GetNbrRow() + 1];
#endif
  this->TemporaryMatrix1 = LongIntegerMatrix(sourceMatrix->NbrRow, sourceMatrix->NbrColumn);
  this->TemporaryMatrix1.Copy(*(this->SourceMatrix));
  this->TemporaryMatrix2 = LongIntegerMatrix(sourceMatrix->NbrRow, sourceMatrix->NbrColumn, true);
  this->OperationType = AbstractArchitectureOperation::LongIntegerMatrixMultiply;
}

// constructor providing the sparse structure of the input matrix
//
// sourceMatrix = pointer to the matrix for which the characteric polynomial should be evaluated
// nbrNonZeroMatrixElements = number of non-zero matrix element per row
// nonZeroMatrixElementPositions = positions of the non-zero matrix element per row

LongIntegerMatrixCharacteristicPolynomialOperation::LongIntegerMatrixCharacteristicPolynomialOperation (LongIntegerMatrix* sourceMatrix, int* nbrNonZeroMatrixElements, int** nonZeroMatrixElementPositions)
{
  this->FirstComponent = 0;
  this->NbrComponent = sourceMatrix->GetNbrRow();
  this->SourceMatrix = sourceMatrix;
  this->NbrNonZeroMatrixElements = nbrNonZeroMatrixElements;
  this->NonZeroMatrixElementPositions = nonZeroMatrixElementPositions;
#ifdef __GMP__
  this->CharacteristicPolynomial = new mpz_t[sourceMatrix->GetNbrRow() + 1];
#else
  this->CharacteristicPolynomial = new LONGLONG[sourceMatrix->GetNbrRow() + 1];
#endif    
  this->TemporaryMatrix1 = LongIntegerMatrix(sourceMatrix->NbrRow, sourceMatrix->NbrColumn);
  this->TemporaryMatrix1.Copy(*(this->SourceMatrix));
  this->TemporaryMatrix2 = LongIntegerMatrix(sourceMatrix->NbrRow, sourceMatrix->NbrColumn, true);
  this->OperationType = AbstractArchitectureOperation::SparseMatrixMatrixMultiply;
}

// copy constructor 
//
// operation = reference on operation to copy

LongIntegerMatrixCharacteristicPolynomialOperation::LongIntegerMatrixCharacteristicPolynomialOperation(const LongIntegerMatrixCharacteristicPolynomialOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->SourceMatrix = operation.SourceMatrix;
  this->NbrNonZeroMatrixElements = operation.NbrNonZeroMatrixElements;
  this->NonZeroMatrixElementPositions = operation.NonZeroMatrixElementPositions;
  this->CharacteristicPolynomial = operation.CharacteristicPolynomial;
  this->TemporaryMatrix1 = operation.TemporaryMatrix1;
  this->TemporaryMatrix2 = operation.TemporaryMatrix2;
  this->OperationType = AbstractArchitectureOperation::SparseMatrixMatrixMultiply;
}
  
// destructor
//

LongIntegerMatrixCharacteristicPolynomialOperation::~LongIntegerMatrixCharacteristicPolynomialOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void LongIntegerMatrixCharacteristicPolynomialOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* LongIntegerMatrixCharacteristicPolynomialOperation::Clone()
{
  return new LongIntegerMatrixCharacteristicPolynomialOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool LongIntegerMatrixCharacteristicPolynomialOperation::RawApplyOperation()
{
  int LastComponent = this->FirstComponent + this->NbrComponent;
  for (int j = this->FirstComponent; j < LastComponent; ++j)
    {
      LongIntegerVector& TmpInputVector = this->TemporaryMatrix1.Columns[j];
      LongIntegerVector& TmpOutputVector = this->TemporaryMatrix2.Columns[j];
      
 #ifdef __GMP__
     for (int i = 0; i < this->SourceMatrix->NbrRow; ++i)
	{
	  mpz_set_ui(TmpOutputVector[i], 0ul);
	  int* TmpNonZeroMatrixElementPositions = this->NonZeroMatrixElementPositions[i];
	  for (int l = 0; l < this->NbrNonZeroMatrixElements[i]; ++l)
	    {
	      mpz_addmul(TmpOutputVector[i], this->SourceMatrix->Columns[TmpNonZeroMatrixElementPositions[l]][i], TmpInputVector[TmpNonZeroMatrixElementPositions[l]]);
	    }
	}	  
#else
     for (int i = 0; i < this->SourceMatrix->NbrRow; ++i)
	{
	  TmpOutputVector[i] = (LONGLONG) 0ul;
	  int* TmpNonZeroMatrixElementPositions = this->NonZeroMatrixElementPositions[i];
	  for (int l = 0; l < this->NbrNonZeroMatrixElements[i]; ++l)
	    {
	      TmpOutputVector[i] += this->SourceMatrix->Columns[TmpNonZeroMatrixElementPositions[l]][i] * TmpInputVector[TmpNonZeroMatrixElementPositions[l]];
	    }
	}	       
#endif     
    }
  LongIntegerMatrix TmpMatrix = this->TemporaryMatrix2;
  this->TemporaryMatrix2 = this->TemporaryMatrix1;
  this->TemporaryMatrix1 = TmpMatrix;
  return true;
}

// apply operation for mono processor architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool LongIntegerMatrixCharacteristicPolynomialOperation::ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture)
{
  this->FirstComponent = 0;
  this->NbrComponent = this->SourceMatrix->GetNbrRow();

#ifdef __GMP__
  mpz_t TmpTrace;
#else
  LONGLONG TmpTrace;
#endif
  this->AlgorithmInitialization(&TmpTrace);

  for (int k = this->SourceMatrix->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->SourceMatrix->NbrRow; ++i)
	{
#ifdef __GMP__
	  mpz_add (this->TemporaryMatrix1.Columns[i][i], this->TemporaryMatrix1.Columns[i][i], this->CharacteristicPolynomial[k + 1]);
#else
	  this->TemporaryMatrix1.Columns[i][i] +=  this->CharacteristicPolynomial[k + 1];
#endif
	}      
      this->RawApplyOperation();
      this->TemporaryMatrix1.Trace(TmpTrace);
#ifdef __GMP__
      mpz_divexact_ui(TmpTrace, TmpTrace, (unsigned long) (this->SourceMatrix->NbrRow - k));
      mpz_neg(TmpTrace, TmpTrace);
      mpz_set(this->CharacteristicPolynomial[k], TmpTrace);
#else
      TmpTrace /= (LONGLONG) (this->SourceMatrix->NbrRow - k);
      TmpTrace *= (LONGLONG) -1l;
      this->CharacteristicPolynomial[k] = TmpTrace;
#endif
    }
#ifdef __GMP__
  mpz_clear(TmpTrace);
#endif
  return true;
}
  
// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool LongIntegerMatrixCharacteristicPolynomialOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->SourceMatrix->GetNbrRow() / architecture->GetNbrThreads();
  int TmpFirstComponent = 0;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  LongIntegerMatrixCharacteristicPolynomialOperation** TmpOperations = new LongIntegerMatrixCharacteristicPolynomialOperation* [architecture->GetNbrThreads()];
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      TmpOperations[i] = (LongIntegerMatrixCharacteristicPolynomialOperation*) this->Clone();
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpFirstComponent += Step;
    }
  TmpOperations[ReducedNbrThreads] = (LongIntegerMatrixCharacteristicPolynomialOperation*) this->Clone();
  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->SourceMatrix->GetNbrRow() - TmpFirstComponent);  
  architecture->SetThreadOperation(TmpOperations[ReducedNbrThreads], ReducedNbrThreads);

#ifdef __GMP__
  mpz_t TmpTrace;
#else
  LONGLONG TmpTrace;
#endif
  this->AlgorithmInitialization(&TmpTrace);

  for (int k = this->SourceMatrix->NbrRow - 2; k >= 0; --k)
    {
      for (int i = 0; i < this->SourceMatrix->NbrRow; ++i)
	{
#ifdef __GMP__
	  mpz_add (this->TemporaryMatrix1.Columns[i][i], this->TemporaryMatrix1.Columns[i][i], this->CharacteristicPolynomial[k + 1]);
#else
	  this->TemporaryMatrix1.Columns[i][i] +=  this->CharacteristicPolynomial[k + 1];
#endif
	}      
      architecture->SendJobs();
      
      LongIntegerMatrix TmpMatrix = this->TemporaryMatrix2;
      this->TemporaryMatrix2 = this->TemporaryMatrix1;
      this->TemporaryMatrix1 = TmpMatrix;

      
      bool ZeroFlag = true;
      for (int i = 0; (i < this->TemporaryMatrix1.NbrColumn) && (ZeroFlag == true); ++i)
	{
	  ZeroFlag = this->TemporaryMatrix1.Columns[i].IsNullVector();
	}
      if (ZeroFlag == true)
	{
	  while (k >= 0)
	    {
#ifdef __GMP__
	      mpz_set_ui(this->CharacteristicPolynomial[k], 0ul);
#else
	      this->CharacteristicPolynomial[k] = (LONGLONG) 0l;
#endif
	      --k;
	    }	  
	}
      else
	{
#ifdef __GMP__
	  this->TemporaryMatrix1.Trace(TmpTrace);
	  mpz_divexact_ui(TmpTrace, TmpTrace, (unsigned long) (this->SourceMatrix->NbrRow - k));
	  mpz_neg(TmpTrace, TmpTrace);
	  mpz_set(this->CharacteristicPolynomial[k], TmpTrace);      
#else
	  TmpTrace /= (LONGLONG) (this->SourceMatrix->NbrRow - k);
	  TmpTrace *= (LONGLONG) -1l;
	  this->CharacteristicPolynomial[k] = TmpTrace;
#endif
	}
    }
#ifdef __GMP__
  mpz_clear(TmpTrace);
#endif

  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

// common code for the algorithm initialization
//
// trace = pointer to the trace
#ifdef __GMP__
void LongIntegerMatrixCharacteristicPolynomialOperation::AlgorithmInitialization(mpz_t* trace)
#else
void LongIntegerMatrixCharacteristicPolynomialOperation::AlgorithmInitialization(LONGLONG* trace)
#endif  
{
#ifdef __GMP__
  for (int i = 0; i <= this->SourceMatrix->NbrRow; ++i)
    {
      mpz_init(this->CharacteristicPolynomial[i]);
    }
  mpz_init(*trace);
  mpz_set_ui(this->CharacteristicPolynomial[this->SourceMatrix->NbrRow], 1ul);
#else
  for (int i = 0; i <= this->SourceMatrix->NbrRow; ++i)
    {
      this->CharacteristicPolynomial[i] = (LONGLONG) 0l;
    }
  this->CharacteristicPolynomial[this->SourceMatrix->NbrRow] = (LONGLONG) 1l;
#endif
  
  this->SourceMatrix->Trace(*trace);
#ifdef __GMP__
  mpz_neg(*trace, *trace);
  mpz_set(this->CharacteristicPolynomial[this->SourceMatrix->NbrRow - 1], (*trace));
#else
  (*trace) *= (LONGLONG) -1l;
  this->CharacteristicPolynomial[this->SourceMatrix->NbrRow - 1] = (*trace);
#endif
}
