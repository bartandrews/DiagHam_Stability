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


#ifndef FERMIONONSPHEREMPSALTERNATIVEWRAPPER_H
#define FERMIONONSPHEREMPSALTERNATIVEWRAPPER_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/SparseComplexMatrix.h"


#include <iostream>

class FermionOnSphereMPSAlternativeWrapper :  public ParticleOnSphere
{
protected:

    // number of fermions
    int NbrFermions;
    // number of fermions plus 1
    int IncNbrFermions;
    // maximum Lz value reached by a fermion
    int LzMax;
    // number of Lz values in a state
    int NbrLzValue;

    // row index of the MPS element that has to be evaluated
    int MPSRowIndex;
    // column index of the MPS element that has to be evaluated
    int MPSColumnIndex;

    // true if the input B matrices are constructed for unnormalized annulus orbitals
    bool UnnormalizedFlag;
    // array of 2x2 matrices for those operators in enum OperatorName
    RealMatrix* Operators;

    // array where the B matrice are stored
    SparseRealMatrix* BMatrices;
    SparseRealMatrix* ConjugateBMatrices;

    // number of quasihole B matrices
    int NbrQuasiholeBMatrices;
    // array where the (optional) quasihole matrices is stored
    SparseComplexMatrix* QuasiholeBMatrices;
    SparseComplexMatrix* ConjugateQuasiholeBMatrices;

    // state normalization. use only the real part in the absence of quasiholes
    Complex StateNormalization;

    // temporary arrays needed for the sparse matrix multiplications
    double* TmpMatrixElements ;
    int* TmpColumnIndices;
    double* TmpElements;
    // maximum number of matrix elements that can ba stored in TmpMatrixElements;
    long MaxTmpMatrixElements;

    // pointer to the architecture
    AbstractArchitecture* Architecture;

public:

    // default constuctor
    //
    FermionOnSphereMPSAlternativeWrapper();

    // basic constructor
    // 
    // nbrFermions = number of fermions
    // lzMax = maximum Lz value reached by a fermion
    // rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
    // columnIndex = column index of the MPS element that has to be evaluated
    // unnormalizedFlag = true if the input B matrices are constructed for unnormalized annulus orbitals
    // bMatrices = array that gives the B matrices 
    // quasiholeBMatrices = array that contains the quasihole B matrices
    // nbrQuasiholeBMatrices = number of quasihole B matrices
    // memory = amount of memory granted for precalculations
    FermionOnSphereMPSAlternativeWrapper(int nbrFermions, int lzMax,
            int rowIndex, int columnIndex, bool unnormalizedFlag,
            SparseRealMatrix* bMatrices, SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholeBMatrices,
            AbstractArchitecture* architecture, unsigned long memory = 10000000);

    // copy constructor (without duplicating datas)
    //
    // fermions = reference on the hilbert space to copy to copy
    FermionOnSphereMPSAlternativeWrapper(const FermionOnSphereMPSAlternativeWrapper& fermions);

    // destructor
    //
    virtual ~FermionOnSphereMPSAlternativeWrapper();

    // assignment (without duplicating datas)
    //
    // fermions = reference on the hilbert space to copy to copy
    // return value = reference on current hilbert space
    FermionOnSphereMPSAlternativeWrapper& operator = (const FermionOnSphereMPSAlternativeWrapper& fermions);

    // clone Hilbert space (without duplicating datas)
    //
    // return value = pointer to cloned Hilbert space
    AbstractHilbertSpace* Clone();

    // get the particle statistic 
    //
    // return value = particle statistic
    virtual int GetParticleStatistic();

    // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
    //
    // index = index of the state on which the operator has to be applied
    // m1 = first index for creation operator
    // m2 = second index for creation operator
    // n1 = first index for annihilation operator
    // n2 = second index for annihilation operator
    // coefficient = reference on the double where the multiplicative factor has to be stored
    // return value = index of the destination state 
    virtual int AdAdAA(int index, int m1, int m2, int n1, int n2, double& coefficient);

    // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
    //
    // index = index of the state on which the operator has to be applied
    // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
    // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
    // nbrIndices = number of creation (or annihilation) operators
    // coefficient = reference on the double where the multiplicative factor has to be stored
    // return value = index of the destination state 
    virtual int ProdAdProdA(int index, int* m, int* n, int nbrIndices, double& coefficient);

    // apply a^+_m a_n operator to a given state 
    //
    // index = index of the state on which the operator has to be applied
    // m = index of the creation operator
    // n = index of the annihilation operator
    // coefficient = reference on the double where the multiplicative factor has to be stored
    // return value = index of the destination state 
    virtual int AdA(int index, int m, int n, Complex& coefficient);
  
    // print a given State
    //
    // Str = reference on current output stream 
    // state = ID of the state to print
    // return value = reference on current output stream 
    virtual ostream& PrintState (ostream& Str, int state);

protected:

    // names of operators stored in .Operators
    enum OperatorName
    {
        IdentityOperator,
        CreationOperator,
        AnnihilationOperator,
        NumberOperator,
        JordanWignerOperator,
        NbrOperatorName
    };

    // initialize .Operators, the 2x2 matrices for those operators in OperatorName
    //
    void InitializeOperators();

    // compute the matrix element of the product of a chain of operators between the two MPSs 
    //
    // ops = array of operators, of length NbrLzValue
    // return value = < MPS | chain of operators | MPS > (not properly normalized)
    Complex ComputeSandwich(OperatorName ops[]);

    // compute the normalization < MPS | MPS > and store it in .StateNormalization
    //
    void ComputeNormalization();
};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereMPSAlternativeWrapper::GetParticleStatistic()
{
    return ParticleOnSphere::FermionicStatistic;
}

#endif

