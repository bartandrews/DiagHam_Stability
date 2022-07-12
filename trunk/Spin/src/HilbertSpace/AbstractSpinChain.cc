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


#include "HilbertSpace/AbstractSpinChain.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>

using std::cout;
using std::endl;


// virtual destructor
//

AbstractSpinChain::~AbstractSpinChain () 
{
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool AbstractSpinChain::WriteHilbertSpace (char* fileName)
{
  cout << "Warning: using dummy method AbstractSpinChain::WriteHilbertSpace" << endl;
  return false;
}

// return value of spin projection on (Oz) for a given state
//
// Str = reference on current output stream 
// return value = spin projection on (Oz)

int AbstractSpinChain::TotalSz (int state)
{
  cout << "warning : TotalSz is not implemented" << endl;
  return 0;
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

int AbstractSpinChain::GetLocalSpin(int site)
{
  return -1;
}

// get the value of the spin (i.e. S) at a given site for a give state
// 
// site = site index
// state = state index in Chain Description
// return value = twice the spin

int AbstractSpinChain::GetLocalSpin(int site, int state)
{
  return this->GetLocalSpin(site);
}



// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Sxi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Syi (int i, Matrix& M)
{
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& AbstractSpinChain::Szi (int i, Matrix& M)
{
  return M;
}


// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::Spi (int i, int state, double& coefficient)
{
  cout << "warning : Spi is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::Smi (int i, int state, double& coefficient)
{
  cout << "warning : Smi is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::Szi (int i, int state, double& coefficient)
{
  cout << "warning : Szi is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double AbstractSpinChain::SziSzj (int i, int j, int state)
{
  cout << "warning : SziSzj is not implemented" << endl;
  return 0.0;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SpiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning : SpiSpj is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SmiSmj (int i, int j, int state, double& coefficient)
{
  cout << "warning : SmiSmj is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SpiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning : SpiSpj is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SmiSzj (int i, int j, int state, double& coefficient)
{
  cout << "warning : SmiSzj is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int AbstractSpinChain::SmiSpj (int i, int j, int state, double& coefficient)
{
  cout << "warning : SmiSpj is not implemented" << endl;
  coefficient = 0.0;
  return this->HilbertSpaceDimension;
}
 


// return index of resulting state from application of S+_i S-_j Sz_k operator on a given state
//
// i = position of S+ operator
// j = position of S- operator
// k = position of Sz operator
// state = index of the state to be applied on S+_i S-_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SpiSmjSzk (int i, int j, int k, int state, double& coefficient)
{
  cout << "warning : SpiSmjSzk is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j Sz_k operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// k = position of Sz operator
// state = index of the state to be applied on S-_i S+_j Sz_k operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SmiSpiSzk (int i, int j, int k, int state, double& coefficient)
{
  cout << "warning : SmiSpjSzk is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::Spi (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : Spi with 2d translations is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::Smi (int i, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : Smi with 2d translations is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : SmiSpj with 2d translations is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : SpiSpj with 2d translations is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : SmiSmj with 2d translations is not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i1 S+_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of leftmost S- operator
// j1 = position of leftmost S+ operator
// i2 = position of rightmost S- operator
// j2 = position of rightmost S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SmiSpjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  cout << "Warning: using dummy method AbstractSpinChain::SmiSpjSmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i1 Sz_j1 S-_i2 S+_j2 operator on a given state
//
// i1 = position of first Sz operator
// j1 = position of second Sz operator
// i2 = position of S- operator
// j2 = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SziSzjSmiSpj (int i1, int j1, int i2, int j2, int state, double& coefficient)
{
  cout << "Warning: using dummy method AbstractSpinChain::SziSzjSmiSpj" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int AbstractSpinChain::SmiSpjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SmiSpjSmkSpl" << endl;
  return this->HilbertSpaceDimension;
}
  
// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state
int AbstractSpinChain::SziSzjSmkSpl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SziSzjSmkSpl" << endl;
  return this->HilbertSpaceDimension;
}
  

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

int AbstractSpinChain::SmiSpjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
					   int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SmiSpjSmkSplSmmSpn" << endl;
  return this->HilbertSpaceDimension;
}

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

int AbstractSpinChain::SmiSpjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
					   int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SmiSpjSzkSzlSmmSpn" << endl;
  return this->HilbertSpaceDimension;
}

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

int AbstractSpinChain::SziSzjSmkSplSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
					   int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SziSzjSmkSplSmmSpn" << endl;
  return this->HilbertSpaceDimension;
}
  
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

int AbstractSpinChain::SziSzjSzkSzlSmmSpn (int i, int j, int k, int l, int m, int n, int state, double& coefficient, 
					   int& nbrTranslationX, int& nbrTranslationY)
 {
  cout << "Warning: using dummy method AbstractSpinChain::SziSzjSzkSzlSmmSpn" << endl;
  return this->HilbertSpaceDimension;
}
 
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

int AbstractSpinChain::SpiSmjSzk (int i, int j, int k, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "Warning: using dummy method AbstractSpinChain::SpiSmjSzk" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of two-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int AbstractSpinChain::Pij (int i, int j, int state)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pij" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of a 3 sites permutation operator on a given state
//
// i = first position
// j = second position (j > i)
// k = third position  (k > j)
// state = index of the state to be applied on P_ijk operator
// return value = index of resulting state

int AbstractSpinChain::Pijk (int i, int j, int k, int state)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pijk" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of a 3 sites permutation inverse operator on a given state
//
// i = first position
// j = second position (j > i)
// k = third position  (k > j)
// state = index of the state to be applied on P_ijk operator
// return value = index of resulting state

int AbstractSpinChain::Pminusijk (int i, int j, int k, int state)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pminusijk" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of a 3 sites permutation operator on a given state
//
// i = first position
// j = second position
// k = third position 
// state = index of the state to be applied on P_ijk operator
// coefficient = reference on the numerical coefficient
// nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
// return value = index of resulting state

int AbstractSpinChain::Pijk (int i, int j, int k, int state, double& coefficient, int& nbrTranslation)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pijk" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of a 3 sites permutation inverse operator on a given state
//
// i = first position
// j = second position 
// k = third position 
// state = index of the state to be applied on P_ijk operator
// coefficient = reference on the numerical coefficient
// nbrTranslations = reference on the number of translations to apply to the resulting state to obtain the canonical state
// return value = index of resulting state

int AbstractSpinChain::Pminusijk (int i, int j, int k, int state, double& coefficient, int& nbrTranslation)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pminusijk" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of four-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int AbstractSpinChain::Pijkl (int i, int j, int k, int l, int state)
{
  cout << "Warning: using dummy method AbstractSpinChain::Pijkl" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Pij operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on Pij operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslationX = reference on the number of translations in the x direction to obtain the canonical form of the resulting state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of resulting state

int AbstractSpinChain::Pij (int i, int j, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : using dummy method AbstractSpinChain::Pij" << endl;
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of four-site exchange operator on a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = index of resulting state

int AbstractSpinChain::Pijkl (int i, int j, int k, int l, int state, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  cout << "warning : using dummy method AbstractSpinChain::Pijkl" << endl;
  return this->HilbertSpaceDimension;
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialDensityMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialDensityMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialEntanglementMatrix(nbrSites, szSector, 0, groundState, architecture);
}	

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. Sz is not conserved.
// 
// nbrSites = number of sites that are part of the A subsystem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return ComplexMatrix();
}


// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return this->EvaluatePartialEntanglementMatrix(nbrSites, szSector, 0, groundState, architecture);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{
  return RealSymmetricMatrix();
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsystem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix AbstractSpinChain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  return HermitianMatrix();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, RealVector& groundState, AbstractArchitecture* architecture)
{
  cout << "warning, EvaluatePartialEntanglementMatrix is not defined" << endl;
  return RealMatrix();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// shift = position of the A part leftmost site within the full system
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, int shift, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  cout << "warning, EvaluatePartialEntanglementMatrix is not defined" << endl;
  return ComplexMatrix();
}
	
// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
// 
// sites = list of sites that define the A subsystem 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  cout << "warning, EvaluatePartialEntanglementMatrix is not defined" << endl;
  return RealMatrix();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. 
// 
// sites = list of sites that define the A subsystem 
// nbrSites = number of sites that are part of the A subsystem 
// szSector = Sz sector in which the density matrix has to be evaluated (disregarded here)
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsystem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix AbstractSpinChain::EvaluatePartialEntanglementMatrix (int* sites, int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  cout << "warning, EvaluatePartialEntanglementMatrix is not defined" << endl;
  return ComplexMatrix();
}

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector AbstractSpinChain::ConvertToKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  cout << "warning, ConvertToKxKyBasis is not defined" << endl;
  return ComplexVector();
}

// convert a state defined in the (Kx,Ky) basis into a state in the real space basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector AbstractSpinChain::ConvertFromKxKyBasis(ComplexVector& state, AbstractSpinChain* space)
{
  cout << "warning, ConvertFromKxKyBasis is not defined" << endl;
  return ComplexVector();
}
  
// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

void AbstractSpinChain::GetBosonicOccupation (unsigned int index, int * finalState)
{
  cout << "warning, using undefined function AbstractSpinChain::GetBosonicOccupation" << endl;
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void AbstractSpinChain::CreateStateFromMPSDescription (SparseRealMatrix* bMatrices, RealVector& state, int mPSRowIndex, int mPSColumnIndex, 
						       long initialIndex, long nbrComponents)
{
  cout << "warning, using undefined function AbstractSpinChain::CreateStateFromMPSDescription" << endl;
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void AbstractSpinChain::CreateStateFromMPSDescription (SparseComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
						       long initialIndex, long nbrComponents)
{
  cout << "warning, using undefined function AbstractSpinChain::CreateStateFromMPSDescription" << endl;
}

// create a state from its MPS description
//
// bMatrices = array that gives the B matrices 
// state = reference to vector that will contain the state description
// mPSRowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// mPSColumnIndex = column index of the MPS element that has to be evaluated
// initialIndex = initial index to compute
// nbrComponents = number of components to compute

void AbstractSpinChain::CreateStateFromMPSDescription (ComplexMatrix* bMatrices, ComplexVector& state, int mPSRowIndex, int mPSColumnIndex, 
						       long initialIndex, long nbrComponents)
{
  cout << "warning, using undefined function AbstractSpinChain::CreateStateFromMPSDescription" << endl;
}

// get the normalization factor in front of each basis state (i.e. 1/sqrt(orbit size))
//
// return value = pointer to normalization factors

double* AbstractSpinChain::GetBasisNormalization()
{
  cout << "warning, using dummy method AbstractSpinChain::GetBasisNormalization" << endl;
  return 0;
}

