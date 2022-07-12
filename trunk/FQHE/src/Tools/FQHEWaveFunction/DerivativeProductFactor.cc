////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                     Copyright (C) 2006 Gunnar Moeller                      //
//                                                                            //
//                                                                            //
//              class for elementary factor in expansion of CF Orbitals       //
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

#include "DerivativeProductFactor.h"
#include "DerivativeProduct.h"

#include <iostream>
using std::cout;
using std::endl;

DerivativeProductFactor::DerivativeProductFactor(JainCFOnSphereOrbitals *CFOs,
						 int UDerivatives, int VDerivatives,
						 int Power, double preFactor)
{
  this->UDerivatives=UDerivatives;
  this->VDerivatives=VDerivatives;
  this->Power=Power;
  this->CFOrbitals=CFOs;
  this->PreFactor=preFactor;
}

DerivativeProductFactor::DerivativeProductFactor( const DerivativeProductFactor &toCopy)
{
  this->UDerivatives=toCopy.UDerivatives;
  this->VDerivatives=toCopy.VDerivatives;
  this->Power=toCopy.Power;
  this->CFOrbitals=toCopy.CFOrbitals;
  this->PreFactor=toCopy.PreFactor;
}

DerivativeProductFactor& DerivativeProductFactor::operator = ( const DerivativeProductFactor &toCopy)
{
  this->UDerivatives=toCopy.UDerivatives;
  this->VDerivatives=toCopy.VDerivatives;
  this->Power=toCopy.Power;
  this->CFOrbitals=toCopy.CFOrbitals;
  this->PreFactor=toCopy.PreFactor;
  return *this;
}

DerivativeProductFactor::DerivativeProductFactor()
{
  this->UDerivatives=0;
  this->VDerivatives=0;
  this->Power=1;
  this->CFOrbitals=NULL;
  this->PreFactor=1.0;
}

Complex* DerivativeProductFactor::getValues()
{  
  return (this->CFOrbitals->DerivativeFactors.GetVector(this->UDerivatives, this->VDerivatives, this->Power-1));
}

Complex DerivativeProductFactor::getValue(int particle)
{
//   cout << "Value: UD=" <<this->UDerivatives<<" VD="<<this->VDerivatives << " Power: "<< this->Power<<
//     " particle " << particle << " result: "<<
//     this->CFOrbitals->DerivativeFactors.Get(this->UDerivatives, this->VDerivatives, this->Power-1, particle) ;
//   if (this->Power>1) cout << "Power 1 was: "<< this->CFOrbitals->DerivativeFactors.Get(this->UDerivatives, this->VDerivatives, 0, particle) << endl;
//   else cout << endl;
  return this->CFOrbitals->DerivativeFactors.Get(this->UDerivatives, this->VDerivatives, this->Power-1, particle);
}

void DerivativeProductFactor::TestHighestPowers()
{
  if (CFOrbitals!=NULL)
    if(CFOrbitals->MaxDerivativePower[this->UDerivatives][this->VDerivatives] < this->Power)
      CFOrbitals->MaxDerivativePower[this->UDerivatives][this->VDerivatives] = this->Power;
}

bool DerivativeProductFactor::isScalar()
{
  if ( (this->UDerivatives==0) && (this->VDerivatives==0) ) return true;
  else return false;
}

bool DerivativeProductFactor::isZero()
{
  if (this->PreFactor==0.0) return true;
  else return false;
}

bool DerivativeProductFactor::isNonZero()
{
  if (this->PreFactor!=0.0) return true;
  else return false;
}

bool DerivativeProductFactor::Multiply(DerivativeProductFactor &other)
{
  if (*this ^ other)
    {
      this->PreFactor *= other.PreFactor;
      this->Power+=other.Power;
      return true;
    }
  else if (this->isScalar())
    {
      double storePrefactor=this->PreFactor;
      *this=other;
      this->PreFactor*=storePrefactor;
      return true;
    }
  else if (other.isScalar())
    {
      this->PreFactor*=other.PreFactor;
      return true;
    }
  else return false;
}


// apply derivatives, prefactors are not included in consideration, yet and have to be corrected!

DerivativeProduct DerivativeProductFactor::Derivative( int DeriveU, int DeriveV)
{
  if (this->isScalar()) // just a number -> derivative is zero!
    return DerivativeProduct(DerivativeProductFactor(this->CFOrbitals, 0,0,1,0.0));  
  else if (this->Power==1)
    {
      DerivativeProductFactor result(this->CFOrbitals, this->UDerivatives+DeriveU,
				     this->VDerivatives+DeriveV, this->Power, this->PreFactor);
      return DerivativeProduct(result);
    }
  else if (this->Power>1)
    {
      DerivativeProduct tmp(DerivativeProductFactor(this->CFOrbitals, this->UDerivatives,
						     this->VDerivatives, this->Power-1, this->PreFactor));
      tmp*=DerivativeProductFactor(this->CFOrbitals, this->UDerivatives+DeriveU,
				      this->VDerivatives+DeriveV, 1, this->Power);
      return tmp;
    }
  // negative Powers are not allowed, thus only power 0 may occur, which yields result: zero!
  else return DerivativeProduct(DerivativeProductFactor(this->CFOrbitals, 0,0,1,0.0));
  
}

bool DerivativeProductFactor::operator < (const DerivativeProductFactor &other)
{
  int order1 = this->VDerivatives | (this->UDerivatives << 10) | (this->Power<<20);
  int order2 = other.VDerivatives | (other.UDerivatives << 10) | (other.Power<<20);
  return (order1 < order2);
}

bool DerivativeProductFactor::operator > (const DerivativeProductFactor &other)
{
  int order1 = this->VDerivatives | (this->UDerivatives << 10) | (this->Power<<20);
  int order2 = other.VDerivatives | (other.UDerivatives << 10) | (other.Power<<20);
  return (order1 > order2);
}
  
bool DerivativeProductFactor::operator ^ (const DerivativeProductFactor &other)
{
  if ((this->UDerivatives == other.UDerivatives) && (this->VDerivatives == other.VDerivatives)
      && (this->CFOrbitals == other.CFOrbitals ))
    return true;
  else 
    return false;
}

bool DerivativeProductFactor::operator == (const DerivativeProductFactor &other)
{
  if ((this->UDerivatives == other.UDerivatives) && (this->VDerivatives == other.VDerivatives)
      && (this->Power == other.Power) && (this->CFOrbitals == other.CFOrbitals ))
    return true;
  else 
    return false;
}

bool DerivativeProductFactor::operator != (const DerivativeProductFactor &other)
{
  if ((this->UDerivatives == other.UDerivatives) && (this->VDerivatives == other.VDerivatives)
      && (this->Power == other.Power) && (this->CFOrbitals == other.CFOrbitals ))
    return false;
  else 
    return true;
}

/*
ostream& operator << (ostream& str, DerivativeProductFactor& D)
{
  if (D.PreFactor != 1.0) str << D.PreFactor << "*";
  if (D.Power>1)
    str<<"(du^"<<D.UDerivatives<<" dv^"<<D.VDerivatives<<"X)^"<<D.Power<<" ";
  else
    str <<"du^"<<D.UDerivatives<<" dv^"<<D.VDerivatives<<"X ";
  return str;
}
*/

ostream& operator << (ostream& str, DerivativeProductFactor& D)
{
  if (D.PreFactor != 1.0) str << D.PreFactor << "*";
  if (D.Power>1)
    {
      if (D.VDerivatives > 0)
	{
	  if (D.UDerivatives+D.VDerivatives>1) str<<"(G_("<<D.UDerivatives<<","<<D.VDerivatives-1<<"))^"<<D.Power;
	  else str<<"G^"<<D.Power;
	}
      else
	{
	   if (D.UDerivatives>1) str<<"(F_"<<D.UDerivatives-1<<")^"<<D.Power;
	   else str<<"F^"<<D.Power;
	}
    }
  else
    {
      if (D.VDerivatives > 0)
	{
	  if (D.UDerivatives+D.VDerivatives>1)
	    str<<"G_("<<D.UDerivatives<<","<<D.VDerivatives-1<<")";
	  else str<<"G";
	}
      else if (D.UDerivatives > 0)
	{
	  if (D.UDerivatives>1)  str<<"F_"<<D.UDerivatives-1;
	  else str << "F";
	}
      else str << "1.0 ";
    }
  return str;
}
