////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2006 Gunnar Moeller                     //
//                                                                            //
//                                                                            //
//             class for elementary factor in expansion of CF Orbitals        //
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

#ifndef DERIVATIVE_PRODUCT_FACTOR
#define DERIVATIVE_PRODUCT_FACTOR

#include "config.h"
#include "MathTools/Complex.h"
#include "JainCFOnSphereOrbitals.h"

class DerivativeProduct;

class DerivativeProductFactor
{
 private:  
  friend class JainCFOnSphereOrbitals;
  friend class DerivativeProduct;
  friend class SumDerivativeProduct;
  
 public:
  DerivativeProductFactor();
  DerivativeProductFactor(JainCFOnSphereOrbitals *CFOs, int UDerivatives, int VDerivatives=0, int Power=1, double preFactor=1.0);
  DerivativeProductFactor(const DerivativeProductFactor &toCopy);
  
  DerivativeProductFactor& operator = ( const DerivativeProductFactor &toCopy);

  
  ~DerivativeProductFactor(){};
  
  Complex* getValues();
  Complex getValue(int particle);
  void TestHighestPowers();
  DerivativeProduct Derivative( int DeriveU, int DeriveV=0);

  bool isScalar();
  bool isZero();
  bool isNonZero();
  
  bool Multiply(DerivativeProductFactor &other);

  // defines an order for DerivativeProductFactor's
  bool operator < (const DerivativeProductFactor &other);
  bool operator > (const DerivativeProductFactor &other);

  // tests whether same derivatives are taken, so that result may be simplified.
  bool operator ^ (const DerivativeProductFactor &other);

  // tests whether same derivatives are taken, and the same total power is present 
  bool operator == (const DerivativeProductFactor &other);
  // or whether this is false
  bool operator != (const DerivativeProductFactor &other);

  // Output Stream overload
  //
  // str = reference on output stream
  // D = DerivativeProductFactor
  // return value = referenceint GetNumSites(){return this->NSites;} on output stream
  friend ostream& operator << (ostream& str, DerivativeProductFactor& D);

 protected:
  JainCFOnSphereOrbitals *CFOrbitals;
  int UDerivatives;
  int VDerivatives;
  int Power;
  double PreFactor;
};


#endif //DERIVATIVE_PRODUCT_FACTOR
