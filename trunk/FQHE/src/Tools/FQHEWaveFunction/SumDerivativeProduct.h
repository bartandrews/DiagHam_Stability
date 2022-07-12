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

#ifndef SUM_DERIVATIVE_PRODUCT
#define SUM_DERIVATIVE_PRODUCT

#include "config.h"
#include "GeneralTools/List.h"
#include "MathTools/Complex.h"
#include "DerivativeProductFactor.h"
#include "DerivativeProduct.h"
#include "JainCFOnSphereOrbitals.h"

class SumDerivativeProduct
{
 private:  
  friend class JainCFOnSphereOrbitals;
  
 public:
  SumDerivativeProduct();
  SumDerivativeProduct(JainCFOnSphereOrbitals *CFOrbitals);
  SumDerivativeProduct( const DerivativeProduct &toWrap);
  SumDerivativeProduct( const DerivativeProductFactor &toWrap);
  SumDerivativeProduct( const SumDerivativeProduct &toCopy);
  
  ~SumDerivativeProduct();

  void Simplify();

  void TestHighestPowers();

  int NumberOfDerivativeProductFactors();

  Complex getValue(int particle);
  void getValues(Complex *result);
  void CommentValues(int particle);

  void fastGetValues(Complex *result);
  void hardwire();
  
  SumDerivativeProduct Derivative( int DeriveU, int DeriveV=0);

  SumDerivativeProduct& operator = (const SumDerivativeProduct& Assign);
  
  SumDerivativeProduct& operator*= (const DerivativeProduct &toMultiply);
  SumDerivativeProduct& operator*= (const DerivativeProductFactor &toMultiply);

  SumDerivativeProduct& operator+= (DerivativeProduct &toAdd);
  SumDerivativeProduct& operator+= (const SumDerivativeProduct &toAdd);
  SumDerivativeProduct& operator+= (const DerivativeProductFactor &toAdd);
  
  friend ostream& operator << (ostream& str, SumDerivativeProduct& S);
  
 protected:
  JainCFOnSphereOrbitals *CFOrbitals;
  List<DerivativeProduct> Summands;  
  GarbageFlag Flag;
  Complex* TmpSum;
  DerivativeProduct **FastSummands;
  int NSummands;
};

// use this routine for a faster calculation,
// calling fastGetValues before hardwire() will cause a segfault!
//
inline void SumDerivativeProduct::fastGetValues(Complex *result)
{
  Complex *CPtr, *CPtr2;
  FastSummands[0]->fastGetValues(result);  
  for (int s=1;s<NSummands; ++s)
    {
      FastSummands[s]->fastGetValues(TmpSum);
      CPtr=result;
      CPtr2=TmpSum;
      for (int i=0; i<CFOrbitals->GetNbrParticles(); ++i)
	*(CPtr++) += *(CPtr2++);
    }
}


#endif //SUM_DERIVATIVE_PRODUCT
