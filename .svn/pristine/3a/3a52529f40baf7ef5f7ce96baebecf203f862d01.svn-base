////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of abstract spin chain                       //
//                                                                            //
//                        last modification : 18/04/2001                      //
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


#ifndef ABSTRACTSPINCHAIN_H
#define ABSTRACTSPINCHAIN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"


class Matrix;
class AbstractArchitecture;
class RealSymmetricMatrix;
class HermitianMatrix;
class ComplexMatrix;
class SparseRealMatrix;
class SparseComplexMatrix;

using std::cout;
using std::endl;

class AbstractSpinChain : public AbstractHilbertSpace
{

 protected:

  int ChainLength;
  // table of indices mapping each spin to the mirror-symmetry conjugate point
  int* MirrorTransformationTable;
  
 public:

  // virtual destructor
  //
  virtual ~AbstractSpinChain ();
  
  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);

  //return value length of the spin chain
  //
  virtual int GetSpinChainLength() const {return this->ChainLength;}

  // get the value of the spin (i.e. S) at a given site
  // 
  // site = site index
  // return value = twice the spin
  virtual int GetLocalSpin(int site);

  // get the value of the spin (i.e. S) at a given site for a give state
  // 
  // site = site index
  // state = state index in Chain Description
  // return value = twice the spin
  virtual int GetLocalSpin(int site, int state);

  // return value of spin projection on (Oz) for a given state
  //
  // Str = reference on current output stream 
  // return value = spin projection on (Oz)
  virtual int TotalSz (int state);

  // return matrix representation of Sx
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Sxi (int i, Matrix& M);

  // return matrix representation of i * Sy
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Syi (int i, Matrix& M);

  // return matrix representation of Sz
  //
  // i = operator position
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Szi (int i, Matrix& M);

  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient);

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient);

  // return index of resulting state from application of Sz_i operator on a given state
  //
  // i = position of Sz operator
  // state = index of the state to be applied on Sz_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int Szi (int i, int state, double& coefficient);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  virtual double SziSzj (int i, int j, int state);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSzj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state 
  virtual int SmiSpj (int i, int j, int state, double& coefficient);
 
  // return index of resulting state from application of S+_i S-_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // state = index of the state to be applied on S+_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmj (int i, int j, int state, double& coefficient);

  // return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // k = position of Sz operator
  // state = index of the state to be applied on S+_i S-_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SpiSmjSzk (int i, int j, int k, int state, double& coefficient);

  // return index of resulting state from application of S-_i S+_j Sz_k operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // k = position of Sz operator
  // state = index of the state to be applied on S-_i S+_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpiSzk (int i, int j, int k, int state, double& coefficient);

  // return index of resulting state from application of S+_i operator on a given state
  //
  // i = position of S+ operator
  // state = index of the state to be applied on S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int Spi (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of S-_i operator on a given state
  //
  // i = position of S- operator
  // state = index of the state to be applied on S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int Smi (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of leftmost S- operator
  // j1 = position of leftmost S+ operator
  // i2 = position of rightmost S- operator
  // j2 = position of rightmost S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
  //
  // i1 = position of first Sz operator
  // j1 = position of second Sz operator
  // i2 = position of S- operator
  // j2 = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  virtual int SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
   // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of S-_i S+_j S-_k S+_l S-_m S+_n operator on a given state
  //
  // i = position of the first S- operator
  // j = position of the first S+ operator
  // k = position of the second S- operator
  // l = position of the second S+ operator
  // m = position of the third S- operator
  // n = position of the third S+ operator
  // state = index of the state to be applied on S-_i S+_j S-_k S+_l S-_m S+_n operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SmiSpjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
				  int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of S-_i S+_j Sz_k Sz_l S-_m S+_n operator on a given state
  //
  // i = position of the first S- operator
  // j = position of the first S+ operator
  // k = position of the first Sz operator
  // l = position of the second Sz operator
  // m = position of the second S- operator
  // n = position of the second S+ operator
  // state = index of the state to be applied on the S-_i S+_j Sz_k Sz_l S-_m S+_n operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SmiSpjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
				  int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of S-_i S+_j Sz_k Sz_l S-_m S+_n operator on a given state
  //
  // i = position of the first Sz operator
  // j = position of the second Sz operator
  // k = position of the first S- operator
  // l = position of the first S+ operator
  // m = position of the second S- operator
  // n = position of the second S+ operator
  // state = index of the state to be applied on the Sz_i Sz_j S-_k S+_l S-_m S+_n operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SziSzjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
				  int& nbrTranslationX, int& nbrTranslationY);
   
  // return index of resulting state from application of S-z_i Sz_j Sz_k Sz_l S-_m S+_n operator on a given state
  //
  // i = position of the first Sz operator
  // j = position of the second Sz operator
  // k = position of the third Sz operator
  // l = position of the fourth Sz operator
  // m = position of the first S- operator
  // n = position of the first S+ operator
  // state = index of the state to be applied on the Sz_i Sz_j Sz_k Sz_l S-_m S+_n operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SziSzjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
				  int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
  //
  // i = position of S+ operator
  // j = position of S- operator
  // k = position of Sz operator
  // state = index of the state to be applied on S+_i S-_j Sz_k operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int SpiSmjSzk (int i, int j, int k, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // return index of resulting state from application of two-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state);
  
  // return index of resulting state from application of a 3 sites permutation operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // return value = index of resulting state
  virtual int Pijk (int i, int j, int k, int state);

  // return index of resulting state from application of a 3 sites permutation inverse operator on a given state
  //
  // i = first position
  // j = second position (j > i)
  // k = third position  (k > j)
  // state = index of the state to be applied on P_ijk operator
  // return value = index of resulting state
  virtual int Pminusijk (int i, int j, int k, int state);

  // return index of resulting state from application of a 3 sites permutation operator on a given state
  //
  // i = first position
  // j = second position
  // k = third position 
  // state = index of the state to be applied on P_ijk operator
  // coefficient = reference on the numerical coefficient
  // nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
  // return value = index of resulting state
  virtual int Pijk (int i, int j, int k, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of a 3 sites permutation inverse operator on a given state
  //
  // i = first position
  // j = second position 
  // k = third position 
  // state = index of the state to be applied on P_ijk operator
  // coefficient = reference on the numerical coefficient
  // nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
  // return value = index of resulting state
  virtual int Pminusijk (int i, int j, int k, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of four-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pijkl (int i, int j, int k, int l, int state);
  
  // return index of resulting state from application of Pij operator (permutation of two spins) on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on Pij operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  // return value = index of resulting state
  virtual int Pij (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);
  
  // return index of resulting state from application of four-site exchange operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = index of resulting state
  virtual int Pijkl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long state){cout <<"using non defined function FindStateIndex in AbstractSpinChain" <<endl; return -1;};

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. Sz is not conserved.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // shift = position of the A part leftmost site within the full system
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
	
  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // sites = list of sites that define the A subsystem 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
  // 
  // sites = list of sites that define the A subsystem 
  // nbrSites = number of sites that are part of the A subsystem 
  // szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

  // convert a state defined in the real space basis into a state in the (Kx,Ky) basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space);

  // convert a state defined in the (Kx,Ky) basis into a state in the real space basis
  //
  // state = reference on the state to convert
  // space = pointer to the Hilbert space where state is defined
  // return value = state in the (Kx,Ky) basis
  virtual ComplexVector ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space);

  // return the Bosonic Occupation of a given state in the basis
  //
  // index = index of the state in the basis
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetBosonicOccupation (unsigned int index, int * finalState);
  
  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  virtual inline unsigned long EncodeSiteState(int physicalState, int sitePosition);
  
  // apply the mirror symmetry to a state
  //
  // stateIndex = index of the current state in the Hlibert space
  // coefficient = coefficient associated with the transformation
  // nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
  // nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
  virtual inline int ApplyMirrorSymmetry(int stateIndex, double& coefficient,
                                             int& nbrTranslationX, int& nbrTranslationY);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long initialIndex = 0l, long nbrComponents = 0l);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long initialIndex = 0l, long nbrComponents = 0l);

  // create a state from its MPS description
  //
  // bMatrices = array that gives the B matrices 
  // state = reference to vector that will contain the state description
  // mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
  // mPSColumnIndex = column index of the MPS element that has to be evaluated
  // initialIndex = initial index to compute
  // nbrComponents = number of components to compute
  virtual void CreateStateFromMPSDescription (ComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
					      long initialIndex = 0l, long nbrComponents = 0l);

  // get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
  //
  // return value = pointer to normalization factors
  virtual double* GetBasisNormalization();
 
};

// return index of resulting state from application of S+_i S-_j operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// state = index of the state to be applied on S+_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

inline int AbstractSpinChain::SpiSmj (int i, int j, int state, double& coefficient)
{
  return this->SmiSpj(j, i, state, coefficient);
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

inline unsigned long AbstractSpinChain::EncodeSiteState(int physicalState, int sitePosition)
{
  cout << "warning, using undefined function AbstractSpinChain::EncodeSiteState(int state, int sitePosition)" << endl;
  return - 1;
}

// apply the mirror symmetry to a state
//
// stateIndex = index of the current state in the Hlibert space
// coefficient = coefficient associated with the transformation
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
inline int AbstractSpinChain::ApplyMirrorSymmetry(int stateIndex, double& coefficient,
                                             int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy function AbstractSpinChain::ApplyMirrorSymmetry" << endl;
  coefficient = 1.0;
  return this->HilbertSpaceDimension;
}


#endif


