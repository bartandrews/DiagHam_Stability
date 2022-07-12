////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere with SU(3) spin                //
//                                                                            //
//                        last modification : 20/01/2008                      //
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
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"

#include <iostream>


using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnSphereWithSU3Spin::~ParticleOnSphereWithSU3Spin ()
{
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSU3Spin::AdA (int index, int m)
{
  return (this->Ad1A1(index, m) + this->Ad2A2(index, m) + this->Ad3A3(index, m));
}

// apply sum_s a^+_m_s a_m_s operator to a given state (sum over all spin states)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double ParticleOnSphereWithSU3Spin::AdA (long index, int m)
{
  return (this->Ad1A1(index, m) + this->Ad2A2(index, m) + this->Ad3A3(index, m));
}

// apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_1 a_m_1

double ParticleOnSphereWithSU3Spin::Ad1A1 (int index, int m)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A1 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_2 a_m_2

double ParticleOnSphereWithSU3Spin::Ad2A2 (int index, int m)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A2 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_3 a_m_3

double ParticleOnSphereWithSU3Spin::Ad3A3 (int index, int m)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A3 not implemented" << endl;
  return 0.0;
}

// apply a^+_m_1 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A1 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A2 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A3 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A1 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A2 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A3 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A1 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A2 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A3 (int index, int m, int n, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A1 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A2 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1A3 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A1 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A2 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2A3 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A1 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A2 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3A3 (int index, int m, int n, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3A3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0, 1 or 2 
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(3) index for the first annihilation operator
// sigma2 = SU(3) index for the second annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{  
  cout << "ParticleOnSphereWithSU3Spin::AsigmaAsigma not implemented" << endl;
  return 0.0;
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0, 1 or 2 
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(3) index for the first annihilation operator
// sigma2 = SU(3) index for the second annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::AsigmaAsigma not implemented" << endl;
  return 0.0;
}

// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A1 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A1 not implemented" << endl;
  return 0.0;
}

// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A2 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A2 not implemented" << endl;
  return 0.0;
}

// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A3 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A3 not implemented" << endl;
  return 0.0;
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A2A2 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A2A2 not implemented" << endl;
  return 0.0;
}


// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A2A3 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A2A3 not implemented" << endl;
  return 0.0;
}


// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A3A3 (int index, int n1, int n2)
{  
  cout << "ParticleOnSphereWithSU3Spin::A3A3 not implemented" << endl;
  return 0.0;
}


// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A1 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A1 not implemented" << endl;
  return 0.0;
}


// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A2 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A2 not implemented" << endl;
  return 0.0;
}


// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A1A3 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A1A3 not implemented" << endl;
  return 0.0;
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A2A2 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A2A2 not implemented" << endl;
  return 0.0;
}

// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A2A3 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A2A3 not implemented" << endl;
  return 0.0;
}

// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value =  multiplicative factor 

double ParticleOnSphereWithSU3Spin::A3A3 (int index, int n1, int n2, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::A3A3 not implemented" << endl;
  return 0.0;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0, 1 or 2
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::AdsigmaAdsigma not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is 0, 1 or 2
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::AdsigmaAdsigma not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad1 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad2 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad3 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2Ad2 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2Ad2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2Ad3 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3Ad3 (int m1, int m2, double& coefficient)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad1 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad1 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad1Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad1Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2Ad2 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad2Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad2Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int ParticleOnSphereWithSU3Spin::Ad3Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation)
{  
  cout << "ParticleOnSphereWithSU3Spin::Ad3Ad3 not implemented" << endl;
  return this->HilbertSpaceDimension;
}

// evaluate wave function in real space using a given basis
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSU3Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis)
{
  return this->EvaluateWaveFunction(state, position, basis, 0, this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSU3Spin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
								 AbstractFunctionBasis& basis, int nextCoordinates)
{
  return this->EvaluateWaveFunctionWithTimeCoherence(state, position, basis, nextCoordinates, 0, 
						     this->HilbertSpaceDimension);
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location
Complex ParticleOnSphereWithSU3Spin::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
						int firstComponent, int nbrComponent)
{
  return Complex(0.0, 0.0);
}


// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex ParticleOnSphereWithSU3Spin::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, 
									 AbstractFunctionBasis& basis, 
									 int nextCoordinates, int firstComponent, 
									 int nbrComponent)
{
  return Complex(0.0, 0.0);
}
                                                                                                                        
// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void ParticleOnSphereWithSU3Spin::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
                                    
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// nbrN1Sector = number of type 1 particles  that belong to the subsytem 
// nbrN2Sector = number of type 1 particles  that belong to the subsytem 
// nbrN3Sector = number of type 1 particles  that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnSphereWithSU3Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, int nbrN1Sector, int nbrN2Sector, int nbrN3Sector, RealVector& groundState, AbstractArchitecture* architecture)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnSphereWithSU3Spin::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int lzSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  HermitianMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}
