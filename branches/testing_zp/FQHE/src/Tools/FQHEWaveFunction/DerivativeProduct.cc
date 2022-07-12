////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2006 Gunnar Moeller                       //
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

#include "DerivativeProduct.h"
#include "DerivativeProductFactor.h"
#include "SumDerivativeProduct.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/List.h"

#include <iostream>
using std::cout;
using std::endl;

DerivativeProduct::DerivativeProduct()
{
  this->CFOrbitals=NULL;
  this->PreFactor=1.0;
  this->Flag.Initialize();
  this->FastProductFactors=NULL;
  this->NFactors=0;
  this->NbrParticles=0;
}

DerivativeProduct::DerivativeProduct( const DerivativeProductFactor &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->PreFactor=1.0;
  this->ProductFactors+=toWrap;
  this->Flag.Initialize();
  this->FastProductFactors=NULL;
  this->NFactors=0;
  this->NbrParticles=this->CFOrbitals->GetNbrParticles();
}

DerivativeProduct::DerivativeProduct(const DerivativeProduct &toCopy)
{
  this->CFOrbitals=toCopy.CFOrbitals;
  this->PreFactor=toCopy.PreFactor;
  this->ProductFactors=toCopy.ProductFactors;
  this->Flag = toCopy.Flag;
  this->FastProductFactors=toCopy.FastProductFactors;
  this->NFactors=toCopy.NFactors;
  this->NbrParticles=toCopy.NbrParticles;
}

DerivativeProduct::DerivativeProduct(DerivativeProduct &Reference, List<DerivativeProductFactor> PriorFactors,
		    List<DerivativeProductFactor> LaterFactors)
{
  this->CFOrbitals=Reference.CFOrbitals;
  this->PreFactor=Reference.PreFactor;
  this->ProductFactors=PriorFactors;
  this->ProductFactors.Link(LaterFactors);
  this->Flag = Reference.Flag;
  this->FastProductFactors=Reference.FastProductFactors;
  this->NFactors=Reference.NFactors;
  this->NbrParticles=this->CFOrbitals->GetNbrParticles();
  //cout << "Nbr Elements after Link: " << this->ProductFactors.GetNbrElement();
}

DerivativeProduct::~DerivativeProduct()
{
  if ((this->FastProductFactors != NULL) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->FastProductFactors;
    }
}

bool DerivativeProduct::isNonZero()
{
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    if(Factor->isZero()) return false;
  return true;
}

Complex DerivativeProduct::getValue(int particle)
{
  Complex result = this->PreFactor;
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    result *=Factor->getValue(particle);
  return result;
}

void DerivativeProduct::getValues(Complex *result)
{
  Complex *TmpProductFactors;
  for (int i=0; i<CFOrbitals->GetNbrParticles(); ++i)
    result[i] = this->PreFactor;
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    {      
      TmpProductFactors = Factor->getValues();
      for (int i=0; i<CFOrbitals->GetNbrParticles(); ++i)
	{
	  result[i] *= TmpProductFactors[i];
	  //	  cout << TmpProductFactors[i] << " =?= "<<  Factor->getValue(i)<<endl;
	}
    }
}

void DerivativeProduct::hardwire()
{
  if (this->NFactors==0)
    {
      this->Simplify(); // get a single prefactor in this DerivativeProduct
      this->NFactors = this->ProductFactors.GetNbrElement();
      this->FastProductFactors = new Complex*[NFactors];
      DerivativeProductFactor *ProductFactor;
      int i=0;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (ProductFactor=LI())!=NULL; )
	FastProductFactors[i++]=ProductFactor->getValues();
    }
}

void DerivativeProduct::TestHighestPowers()
{
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    Factor->TestHighestPowers();
}

int DerivativeProduct::NumberOfDerivativeProductFactors()
{
  return this->ProductFactors.GetNbrElement();
} 

void DerivativeProduct::Simplify()
{
  DerivativeProductFactor *Factor;
  int pos=0;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;++pos)
    {
      this->PreFactor*=Factor->PreFactor;
      Factor->PreFactor=1.0;
      if (Factor->isScalar())
	{
	  //cout << "found scalar factor to simplify!" << endl;
	  this->ProductFactors.DeleteElement(pos);
	}
    }
  if (this->PreFactor==0.0)
    this->ProductFactors.DeleteList();
}


SumDerivativeProduct DerivativeProduct::Derivative(int DeriveU, int DeriveV)
{
  SumDerivativeProduct result(this->CFOrbitals);
  List<DerivativeProductFactor> PriorFactors;
  List<DerivativeProductFactor> LaterFactors;
  int elements=this->ProductFactors.GetNbrElement();
  //cout <<endl<< "Derivative of " << *this << endl;
  for (int pos=0; pos < elements; ++pos)
    {
      DerivativeProduct derivative = this->ProductFactors[pos].Derivative(DeriveU,DeriveV);
      //cout << "Element " << pos << "= " << this->ProductFactors[pos]<< endl;
      //cout << "its der. " << derivative << endl;      
      if (derivative.isNonZero())
	{
	  PriorFactors.DeleteList();
	  if (pos>0) PriorFactors=Extract(this->ProductFactors, 0, pos);
	  LaterFactors.DeleteList();
	  if (pos+1 < elements) LaterFactors=Extract(this->ProductFactors, pos+1);      
	  DerivativeProduct tmp(*this,PriorFactors,LaterFactors);
	  //cout << "other terms: " << tmp << endl;
	  tmp*=derivative;
	  //cout << "Full term: " << tmp << endl;
	  result +=tmp;
	}
    }
  return result;
}

// assignment operator
DerivativeProduct& DerivativeProduct::operator = (const DerivativeProduct& Assign)
{
  this->CFOrbitals=Assign.CFOrbitals;
  this->PreFactor=Assign.PreFactor;
  this->ProductFactors=Assign.ProductFactors;
  this->Flag = Assign.Flag;
  this->FastProductFactors=Assign.FastProductFactors;
  this->NFactors=Assign.NFactors;
  this->NbrParticles=Assign.NbrParticles;  
  return *this;
}


/*DerivativeProduct& DerivativeProduct::operator*= (DerivativeProductFactor &toMultiply)
{
  DerivativeProductFactor *Factor;
  bool toBeInserted=true;
  int pos=0;
  if (toMultiply.isScalar()) // simple number to be multiplied...
    {
      this->PreFactor*=toMultiply.PreFactor;
      return *this;
    }
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL; ++pos)
    {
      if ( (*Factor) < toMultiply )
	{
	  if (Factor->isScalar())
	    Factor->Multiply(toMultiply);
	  else
	    this->ProductFactors.Insert(toMultiply,pos);
	  toBeInserted=false;
	  break;
	}
      else if ( Factor->Multiply(toMultiply) )
	{
	  toBeInserted=false;
	  break;
	}
    }
  if (toBeInserted) this->ProductFactors+=toMultiply;
  return *this;
}
*/

DerivativeProduct& DerivativeProduct::operator*= (DerivativeProductFactor &toMultiply)
{
  DerivativeProductFactor *Factor;
  bool toBeInserted=true;
  if (this->ProductFactors.GetNbrElement()>0)
    {
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  if ( Factor->Multiply(toMultiply) )
	    {
	      toBeInserted=false;
	      break;
	    }
	}
    }
  if (toBeInserted) this->ProductFactors+=toMultiply;
  this->Simplify();
  return *this;
}

DerivativeProduct& DerivativeProduct::operator*= (const DerivativeProduct &toMultiply)
{
  DerivativeProductFactor *Factor;
  if (this->ProductFactors.GetNbrElement()>0)
    {
      this->PreFactor*=toMultiply.PreFactor;
      for (ListIterator<DerivativeProductFactor> LI(toMultiply.ProductFactors); (Factor=LI())!=NULL; )
	(*this)*=(*Factor);
    }
  else
    {
      double storePreFactor=this->PreFactor;
      *this = toMultiply;
      this->PreFactor*=storePreFactor;
    }
  return *this;
}


bool DerivativeProduct::operator ^ (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() != other.ProductFactors.GetNbrElement())
    return false;
  else
    {
      this->ProductFactors.UpOrder();
      other.ProductFactors.UpOrder();
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  if ( (*Factor) !=  (*OtherFactor) ) return false;
	}
      return true;
    }
}

/*  // Alternative formulation of the comparation operator
bool DerivativeProduct::operator ^ (DerivativeProduct &other)
{
  this->Simplify();
  other.Simplify();
  if (this->ProductFactors.GetNbrElement() != other.ProductFactors.GetNbrElement())
    return false;
  else
  {
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      bool isPresent;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  isPresent=false;
	  for (ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors); (OtherFactor=OtherLI())!=NULL;)
	    if ( (*Factor) == (*OtherFactor) )
	      {
		isPresent=true;
		break;
	      }
	  if (isPresent==false) return false;
	}
      return true;
    }
}
*/

bool DerivativeProduct::operator < (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() < other.ProductFactors.GetNbrElement())
    return true;
  else
    {     
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  //cout << "Comparing factors " << *Factor << " and " << *OtherFactor << endl;
	  if ( (*Factor) <  (*OtherFactor) ) return true;
	}
      return false;
    }
}

bool DerivativeProduct::operator > (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() > other.ProductFactors.GetNbrElement())
    return true;
  else
    {     
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  //cout << "Comparing factors " << *Factor << " and " << *OtherFactor << endl;
	  if ( (*Factor) >  (*OtherFactor) ) return true;
	}
      return false;
    }
}

ostream& operator << (ostream& str, DerivativeProduct& D)
{
  DerivativeProductFactor *Factor;
  if (D.ProductFactors.GetNbrElement()>=1)
    {
      if (D.PreFactor!=1.0) str << D.PreFactor << "*" ;
      ListIterator<DerivativeProductFactor> LI(D.ProductFactors);
      Factor=LI();
      str << (*Factor);
      for (; (Factor=LI())!=NULL;)
	str <<"*"<< (*Factor);    
    }
  else str << "1.0";
  str << " ";
  return str;
}
