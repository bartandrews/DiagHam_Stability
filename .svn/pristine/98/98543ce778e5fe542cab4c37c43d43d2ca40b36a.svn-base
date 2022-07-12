////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                       Copyright (C) 2006 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//            class for elementary factor in expansion of CF Orbitals         //
//                                                                            //
//                        last modification : 17/04/2006                      //
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

#ifndef DERIVATIVE_PRODUCT
#define DERIVATIVE_PRODUCT

#include "config.h"
#include "GeneralTools/List.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "DerivativeProductFactor.h"
#include "JainCFOnSphereOrbitals.h"

class SumDerivativeProduct;

class DerivativeProduct
{
 private:  
  friend class JainCFOnSphereOrbitals;
  friend class SumDerivativeProduct;

 protected:
  void fastGetValues(Complex *result);
  void hardwire();
  
 public:
  DerivativeProduct();
  DerivativeProduct( const DerivativeProductFactor &toWrap);
  DerivativeProduct( const DerivativeProduct &toCopy);
  DerivativeProduct(DerivativeProduct &Reference, List<DerivativeProductFactor> PriorFactors,
		    List<DerivativeProductFactor> LaterFactors);
  
  ~DerivativeProduct();

  void Simplify();

  bool isNonZero();

  void getValues(Complex* result);
  Complex getValue(int particle);

  void TestHighestPowers();
  
  int NumberOfDerivativeProductFactors();
  
  SumDerivativeProduct Derivative( int DeriveU, int DeriveV=0);

  // assignment operator
  DerivativeProduct& operator = (const DerivativeProduct& Assign);
  
  DerivativeProduct& operator*= (const DerivativeProduct &toMultiply);
  DerivativeProduct& operator*= (DerivativeProductFactor &toMultiply);

  bool operator ^ ( DerivativeProduct &other);
  bool operator < ( DerivativeProduct &other);
  bool operator > ( DerivativeProduct &other);

  // Output Stream overload
  //
  // str = reference on output stream
  // D = DerivativeProductFactor
  // return value = referenceint GetNumSites(){return this->NSites;} on output stream
  friend ostream& operator << (ostream& str, DerivativeProduct& D);

  
 protected:
  JainCFOnSphereOrbitals *CFOrbitals;
  List<DerivativeProductFactor> ProductFactors;
  double PreFactor;
  GarbageFlag Flag;
  Complex **FastProductFactors;
  int NFactors;
  int NbrParticles;
};

// fast evaluation routine:
// same as getValues, but uses faster access
// the return value is stored in the array "result" that has to be already reserved

inline void DerivativeProduct::fastGetValues(Complex *result)
{
  Complex *CPtr=result, *CPtr2;
  for (int i=0; i<this->NbrParticles; ++i)
    *(CPtr++) = this->PreFactor;
  for (int f=0;f<NFactors; ++f)
    {      
      CPtr=result;
      CPtr2=FastProductFactors[f];
      for (int i=0; i<this->NbrParticles; ++i)
	*(CPtr++) *= *(CPtr2++);
    }
}

#endif //DERIVATIVE_PRODUCT
