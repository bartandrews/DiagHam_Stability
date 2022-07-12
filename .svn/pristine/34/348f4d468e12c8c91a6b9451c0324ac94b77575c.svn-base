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


#ifndef FERMIONONCYLINDERMPSWRAPPER_H
#define FERMIONONCYLINDERMPSWRAPPER_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "Matrix/SparseRealMatrix.h"


#include <iostream>

class FermionOnCylinderMPSWrapper :  public FermionOnSphereMPSWrapper
{

 protected:

  // number of matrix products that are precomputed

  int NbrPrecalculatedMatrixProducts;

 public:

  // default constuctor
  //
  FermionOnCylinderMPSWrapper();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  // lzMax = maximum Lz value reached by a fermion
  // referenceState = array that describes the root configuration
  // rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // columnIndex = column index of the MPS element that has to be evaluated
  // bMatrices = array that gives the B matrices 
  // architecture = pointer to the archiecture
  // memory = amount of memory granted for precalculations
  FermionOnCylinderMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,  
			     int rowIndex, int columnIndex, SparseRealMatrix* bMatrices, AbstractArchitecture* architecture, 
			     unsigned long memory = 10000000);

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
  FermionOnCylinderMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,  
			       int rowIndex, int columnIndex, SparseRealMatrix* bMatrices, 
			       SparseComplexMatrix* quasiholeBMatrices, int nbrQuasiholes,
			       AbstractArchitecture* architecture, unsigned long memory = 10000000);


  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnCylinderMPSWrapper(const FermionOnCylinderMPSWrapper& fermions);

  // destructor
  //
  virtual ~FermionOnCylinderMPSWrapper ();

  // assignment (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnCylinderMPSWrapper& operator = (const FermionOnCylinderMPSWrapper& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AAA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdA (int index, int m, int n, double& coefficient);

  // apply a^+_m a_n operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual long AdA (long index, int m, int n, Complex& coefficient);
    

 protected:
  
  // core part of the  a^+_m a_n operator calculation
  //
  // m = index of the creation operator
  // n = index of the annihilation operator
  // return value = matrix that results from the a^+_m a_n operator calculation
  virtual SparseRealMatrix AdACore (int m, int n);
    
  // compute the B matrix contribution to the state normalization
  //
  // return value = matrix that provides the B matrix contribution
  virtual SparseRealMatrix ComputeBMatrixNormalization();

};



#endif


