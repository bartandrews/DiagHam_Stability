////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                          class author: Yang-Le Wu                          //
//                                                                            //
//      class of fermions on sphere that allow to use MPS with operator       //
//                                                                            //
//                        last modification : 17/05/2013                      //
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
#include "HilbertSpace/FermionOnSphereMPSAlternativeWrapper.h"
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

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constuctor
//

FermionOnSphereMPSAlternativeWrapper::FermionOnSphereMPSAlternativeWrapper()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// lzMax = maximum Lz value reached by a fermion
// rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// columnIndex = column index of the MPS element that has to be evaluated
// bMatrices = array that gives the B matrices 
// quasiholeBMatrices = array that contains the quasihole B matrices
// nbrQuasiholeBMatrices = number of quasihole B matrices
// memory = amount of memory granted for precalculations
FermionOnSphereMPSAlternativeWrapper::FermionOnSphereMPSAlternativeWrapper(int nbrFermions, int lzMax,
        int rowIndex, int columnIndex, bool unnormalizedFlag,
        SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholeBMatrices,
        AbstractArchitecture* architecture, unsigned long memory)
{
    this->NbrFermions = nbrFermions;
    this->IncNbrFermions = this->NbrFermions + 1;
    this->LzMax = lzMax;
    this->NbrLzValue = this->LzMax + 1;

    this->LargeHilbertSpaceDimension = 1l;
    this->HilbertSpaceDimension = 1;
    this->Flag.Initialize();
    this->Architecture = architecture;

    this->MPSRowIndex = rowIndex;
    this->MPSColumnIndex = columnIndex;
    this->UnnormalizedFlag = unnormalizedFlag;

    this->Operators = 0;
    this->InitializeOperators();

    int NbrBMatrices = 2;
    this->BMatrices = new SparseRealMatrix[NbrBMatrices];
    this->ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
    for (int i = 0; i < NbrBMatrices; ++i)
    {
        this->BMatrices[i] = bMatrices[i];
        this->ConjugateBMatrices[i] = bMatrices[i].Transpose();
    }

    this->NbrQuasiholeBMatrices = nbrQuasiholeBMatrices;
    if (this->NbrQuasiholeBMatrices > 0)
    {
        this->QuasiholeBMatrices = new SparseComplexMatrix[this->NbrQuasiholeBMatrices];
        this->ConjugateQuasiholeBMatrices = new SparseComplexMatrix[this->NbrQuasiholeBMatrices];
        for (int i = 0; i < this->NbrQuasiholeBMatrices; ++i)
        {
            this->QuasiholeBMatrices[i] = quasiholeBMatrices[i];
            this->ConjugateQuasiholeBMatrices[i] = quasiholeBMatrices[i].HermitianTranspose();
        }
    }
    else
    {
        this->QuasiholeBMatrices = 0;
        this->ConjugateQuasiholeBMatrices = 0;
    }

    this->MaxTmpMatrixElements = ((long) this->BMatrices[0].GetNbrRow()) * ((long) this->BMatrices[0].GetNbrRow());
    cout << "Requested memory for sparse matrix multiplications = ";
    cout << ((this->MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
    this->TmpMatrixElements = new double[this->MaxTmpMatrixElements];
    this->TmpColumnIndices = new int[this->MaxTmpMatrixElements];
    this->TmpElements = new double[this->BMatrices[0].GetNbrRow()];

    this->ComputeNormalization();
    cout << "MPS norm = " << this->StateNormalization << endl;
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereMPSAlternativeWrapper::FermionOnSphereMPSAlternativeWrapper(const FermionOnSphereMPSAlternativeWrapper& fermions)
{
    this->NbrFermions = fermions.NbrFermions;
    this->IncNbrFermions = fermions.IncNbrFermions;
    this->LzMax = fermions.LzMax;
    this->NbrLzValue = fermions.NbrLzValue;
    this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
    this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
    this->Flag = fermions.Flag;
    this->MPSRowIndex = fermions.MPSRowIndex;
    this->MPSColumnIndex = fermions.MPSColumnIndex;
    this->UnnormalizedFlag = fermions.UnnormalizedFlag;
    this->Operators = fermions.Operators;
    this->BMatrices = fermions.BMatrices;
    this->ConjugateBMatrices = fermions.ConjugateBMatrices;
    this->NbrQuasiholeBMatrices = fermions.NbrQuasiholeBMatrices;
    this->QuasiholeBMatrices = fermions.QuasiholeBMatrices;
    this->ConjugateQuasiholeBMatrices = fermions.ConjugateQuasiholeBMatrices;
    this->StateNormalization = fermions.StateNormalization;
    this->MaxTmpMatrixElements = fermions.MaxTmpMatrixElements;
    this->TmpMatrixElements = new double [this->MaxTmpMatrixElements];
    this->TmpColumnIndices = new int [this->MaxTmpMatrixElements];
    this->TmpElements = new double [this->BMatrices[0].GetNbrRow()];
    this->Architecture = fermions.Architecture;
}

// destructor
//

FermionOnSphereMPSAlternativeWrapper::~FermionOnSphereMPSAlternativeWrapper ()
{
    if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
        delete[] this->Operators;
        delete[] this->BMatrices;
        delete[] this->ConjugateBMatrices;
        delete[] this->QuasiholeBMatrices;
        delete[] this->ConjugateQuasiholeBMatrices;
    }
    delete[] this->TmpMatrixElements;
    delete[] this->TmpColumnIndices;
    delete[] this->TmpElements;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereMPSAlternativeWrapper& FermionOnSphereMPSAlternativeWrapper::operator = (const FermionOnSphereMPSAlternativeWrapper& fermions)
{
    if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
        delete[] this->Operators;
        delete[] this->BMatrices;
        delete[] this->ConjugateBMatrices;
        delete[] this->QuasiholeBMatrices;
        delete[] this->ConjugateQuasiholeBMatrices;
    }
    delete[] this->TmpMatrixElements;
    delete[] this->TmpColumnIndices;
    delete[] this->TmpElements;

    this->NbrFermions = fermions.NbrFermions;
    this->IncNbrFermions = fermions.IncNbrFermions;
    this->LzMax = fermions.LzMax;
    this->NbrLzValue = fermions.NbrLzValue;
    this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
    this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
    this->Flag = fermions.Flag;
    this->MPSRowIndex = fermions.MPSRowIndex;
    this->MPSColumnIndex = fermions.MPSColumnIndex;
    this->UnnormalizedFlag = fermions.UnnormalizedFlag;
    this->Operators = fermions.Operators;
    this->BMatrices = fermions.BMatrices;
    this->ConjugateBMatrices = fermions.ConjugateBMatrices;
    this->NbrQuasiholeBMatrices = fermions.NbrQuasiholeBMatrices;
    this->QuasiholeBMatrices = fermions.QuasiholeBMatrices;
    this->ConjugateQuasiholeBMatrices = fermions.ConjugateQuasiholeBMatrices;
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

AbstractHilbertSpace* FermionOnSphereMPSAlternativeWrapper::Clone()
{
    return new FermionOnSphereMPSAlternativeWrapper(*this);
}

// initialize .Operators, the 2x2 matrices for those operators in OperatorName
//

void FermionOnSphereMPSAlternativeWrapper::InitializeOperators()
{
    RealMatrix* op;
    this->Operators = new RealMatrix[NbrOperatorName];

    // row marks final state, column marks initial state
    // also, remember the ComplexMatrix[column][row] convention
    // so, [initial][final] = amplitude

    op = this->Operators + IdentityOperator;
    op->ResizeAndClean(2, 2);
    (*op)[0][0] = 1.0;
    (*op)[1][1] = 1.0;

    op = this->Operators + CreationOperator;
    op->ResizeAndClean(2, 2);
    (*op)[0][1] = 1.0;

    op = this->Operators + AnnihilationOperator;
    op->ResizeAndClean(2, 2);
    (*op)[1][0] = 1.0;

    op = this->Operators + NumberOperator;
    op->ResizeAndClean(2, 2);
    (*op)[1][1] = 1.0;

    op = this->Operators + JordanWignerOperator;
    op->ResizeAndClean(2, 2);
    (*op)[0][0] = 1.0;
    (*op)[1][1] = -1.0;
}


// compute the matrix element of the product of a chain of operators between two MPSs 
//
// ops = array of operators, of length NbrLzValue
// return value = < MPS | chain of operators | MPS >

Complex FermionOnSphereMPSAlternativeWrapper::ComputeSandwich(OperatorName ops[])
{
    // matrix configurations:
    //
    //    row ----B---- column
    // column ---B^T--- row
    //
    //           m1
    //            |
    //        operator
    //            |
    //           m2
    //
    //  column ---B^T-----+----------+
    //             |      |          |
    //         operator   | Sandwich |
    //             |      |          |
    //    row  ----B------+----------+

    int NbrRow = this->BMatrices[0].GetNbrRow();
    SparseRealMatrix Sandwich(NbrRow, NbrRow);
    Sandwich.SetMatrixElement(this->MPSColumnIndex, this->MPSColumnIndex, 1.0);

    // square norm of the unnormalized one-body wave-function, 2 * pi * 2^j * j!
    double SqrNormAnnulus = 2 * 3.1415926535897932385;
    for (int j = 1; j <= this->LzMax + 1; ++j)
        SqrNormAnnulus *= 2 * j;

    for (int j = this->LzMax; j >= 0; --j)
    {
        OperatorName op = ops[j];
        SqrNormAnnulus /= 2 * (j + 1);

        SparseRealMatrix M(NbrRow, NbrRow);
        for (int m1 = 0; m1 < 2; ++m1) // final (bra)
            for (int m2 = 0; m2 < 2; ++m2) // initial (ket)
                if (this->Operators[op][m2][m1] != 0.0)
                {
                    double coeff = this->Operators[ops[j]][m2][m1];
                    if (this->UnnormalizedFlag)
                    {
                        if ((op == CreationOperator) || (op == AnnihilationOperator))
                            coeff *= SqrNormAnnulus;
                        else if (((op == IdentityOperator) || (op == JordanWignerOperator)) && (m1 == 1) && (m2 == 1))
                            coeff *= SqrNormAnnulus;
                        else if (op == NumberOperator)
                            coeff *= SqrNormAnnulus * SqrNormAnnulus;
                    }

                    SparseRealMatrix Tmp = Conjugate(this->BMatrices + m2, &Sandwich, this->ConjugateBMatrices + m1,
                            this->TmpMatrixElements, this->TmpColumnIndices, this->MaxTmpMatrixElements, this->Architecture);
                    Tmp *= coeff;
                    M = M + Tmp;
                }
        Sandwich = M;
    }

    if (this->NbrQuasiholeBMatrices == 0)
    {
        double Result = 0;
        Sandwich.GetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, Result);
        return Complex(Result, 0);
    }
    else
    {
        Complex Result;
        SparseComplexMatrix M(Sandwich);
        M = Conjugate(this->QuasiholeBMatrices[0], M, this->ConjugateQuasiholeBMatrices[0]); // 1x1 matrix
        M.GetMatrixElement(this->MPSRowIndex, this->MPSRowIndex, Result);
        return Result;
    }
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

int FermionOnSphereMPSAlternativeWrapper::AdAdAA(int index, int m1, int m2, int n1, int n2, double& coefficient)
{
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

int FermionOnSphereMPSAlternativeWrapper::ProdAdProdA(int index, int* m, int* n, int nbrIndices, double& coefficient)
{
    return 0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereMPSAlternativeWrapper::AdA(int index, int m, int n, Complex& coefficient)
{
    if ((m < 0) || (m > this->LzMax) || (n < 0) || (n > this->LzMax))
    {
        coefficient = 0;
        return 0;
    }

    OperatorName* ops = new OperatorName[this->NbrLzValue];
    if (m == n)
    {
        for (int j = 0; j < this->NbrLzValue; ++j)
            ops[j] = (j == m)?NumberOperator:IdentityOperator;
    }
    else
    {
        for (int j = 0; j < this->NbrLzValue; ++j)
        {
            if (j == m)
                ops[j] = CreationOperator;
            else if (j == n)
                ops[j] = AnnihilationOperator;
            else if ((j < m) != (j < n))
                ops[j] = JordanWignerOperator;
            else
                ops[j] = IdentityOperator;
        }
    }
    coefficient = this->ComputeSandwich(ops) / this->StateNormalization; // no extra sign
    delete[] ops;
    return 0;
}

// compute the normalization < MPS | MPS > and store it in .StateNormalization
//
void FermionOnSphereMPSAlternativeWrapper::ComputeNormalization()
{
    OperatorName *ops = new OperatorName[this->NbrLzValue];
    for (int i = 0; i < this->NbrLzValue; ++i)
        ops[i] = IdentityOperator;
    this->StateNormalization = this->ComputeSandwich(ops);
    delete[] ops;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& FermionOnSphereMPSAlternativeWrapper::PrintState(ostream& Str, int state)
{
    return Str;
}
