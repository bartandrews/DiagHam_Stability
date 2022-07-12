////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                     Copyright (C) 2006 Gunnar Moeller                      //
//                                                                            //
//                                                                            //
//           class for elementary factor in expansion of CF Orbitals          //
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

#include "SumDerivativeProduct.h"
#include "GeneralTools/ListIterator.h"

#include <iostream>
using std::cout;
using std::endl;


SumDerivativeProduct::SumDerivativeProduct()
{
  this->CFOrbitals=NULL;
  this->TmpSum=NULL;
  this->Flag.Initialize();
  this->FastSummands=NULL;
  this->NSummands=0;
}

SumDerivativeProduct::SumDerivativeProduct(JainCFOnSphereOrbitals *CFOrbitals)
{
  this->CFOrbitals=CFOrbitals;
  this->Flag.Initialize();
  this->TmpSum=NULL;
  this->FastSummands=NULL;
  this->NSummands=0;
}

SumDerivativeProduct::SumDerivativeProduct(const DerivativeProductFactor &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->Summands+=DerivativeProduct(toWrap);
  this->Flag.Initialize();
  if (this->CFOrbitals!=NULL) this->TmpSum = new Complex[this->CFOrbitals->GetNbrParticles()];
  else this->TmpSum = NULL;
  this->FastSummands=NULL;
  this->NSummands=0;
}

SumDerivativeProduct::SumDerivativeProduct(const DerivativeProduct &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->Summands+=toWrap;
  this->Flag.Initialize();
  if (this->CFOrbitals!=NULL) this->TmpSum = new Complex[this->CFOrbitals->GetNbrParticles()];
  else this->TmpSum = NULL;
  this->FastSummands=NULL;
  this->NSummands=0;
}

SumDerivativeProduct::SumDerivativeProduct(const SumDerivativeProduct &toCopy)
{
  this->CFOrbitals=toCopy.CFOrbitals;
  this->Summands=toCopy.Summands;
  this->TmpSum=toCopy.TmpSum;
  this->Flag = toCopy.Flag;
  this->FastSummands=toCopy.FastSummands;
  this->NSummands=toCopy.NSummands;
}

SumDerivativeProduct::~SumDerivativeProduct()
{
  if ((this->TmpSum != NULL) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->TmpSum;     
    }
  if ((this->FastSummands != NULL) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->FastSummands;     
    }
}

Complex SumDerivativeProduct::getValue(int particle)
{
  Complex Sum(0.0);
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    Sum+=Product->getValue(particle);
  return Sum;
}

void SumDerivativeProduct::getValues(Complex *result)
{
  Complex *CPtr, *CPtr2;
  DerivativeProduct *Product;
  ListIterator<DerivativeProduct> LI(this->Summands);
  Product=LI();
  Product->getValues(result);  
  while((Product=LI())!=NULL)
    {
      Product->getValues(TmpSum);
      CPtr=result;
      CPtr2=TmpSum;
      for (int i=0; i<CFOrbitals->GetNbrParticles(); ++i)
	//result[i] += TmpSum[i];
	*(CPtr++) += *(CPtr2++);
    }
}

void SumDerivativeProduct::CommentValues(int particle)
{
  DerivativeProduct *Product;
  ListIterator<DerivativeProduct> LI(this->Summands);
  while((Product=LI())!=NULL)
    {
      Product->getValues(TmpSum);
      cout << "P="<<*Product<<"= " <<TmpSum[particle]<<endl;
    }
}

void SumDerivativeProduct::hardwire()
{
  if (this->NSummands==0)
    {
      this->NSummands = this->Summands.GetNbrElement();
      this->FastSummands = new DerivativeProduct*[NSummands];
      DerivativeProduct *Product;
      int i=0;
      for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
	{
	  FastSummands[i++]=Product;
	  Product->hardwire();
	}
    }
}

void SumDerivativeProduct::TestHighestPowers()
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    Product->TestHighestPowers();
}

int SumDerivativeProduct::NumberOfDerivativeProductFactors()
{
  int sum=0;
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    sum+=Product->NumberOfDerivativeProductFactors();
  return sum;
} 

SumDerivativeProduct SumDerivativeProduct::Derivative(int DeriveU, int DeriveV)
{
  SumDerivativeProduct result(this->CFOrbitals);
  DerivativeProduct *Product;
  //  cout << "Derivative of SumDerivativeProduct " << *this << endl;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    result += Product->Derivative(DeriveU, DeriveV);
  return result;
}

SumDerivativeProduct& SumDerivativeProduct::operator = (const SumDerivativeProduct& Assign)
{
  this->CFOrbitals=Assign.CFOrbitals;
  this->Summands=Assign.Summands;
  this->Flag=Assign.Flag;
  this->TmpSum=Assign.TmpSum;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator*=(const DerivativeProduct &toMultiply)
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    (*Product)*=toMultiply;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator*= (const DerivativeProductFactor &toMultiply)
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    (*Product)*=toMultiply;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= ( DerivativeProduct &toAdd)
{
  DerivativeProduct *Summand;
  bool toBeInserted=true;
  int pos=0;
  //cout << "Adding " << toAdd << " to Sum " << *this << endl;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Summand=LI())!=NULL; ++pos)
    {
      if ( *Summand ^ toAdd ) // may be added together
	{
	  toBeInserted=false;
	  toAdd.Simplify();
	  Summand->Simplify();
	  Summand->PreFactor+=toAdd.PreFactor;
	}
    }
  if (toBeInserted) this->Summands+=toAdd;
  //cout << "Result: " << *this<< endl;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= (const SumDerivativeProduct &toAdd)
{
  DerivativeProduct *Summand;
  for (ListIterator<DerivativeProduct> LI(toAdd.Summands); (Summand=LI())!=NULL;)
    {
      (*this)+=(*Summand);
    }
  this->Summands.UpOrder();  // clean up things a little...
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= (const DerivativeProductFactor &toAdd)
{
  (*this)+=DerivativeProduct(toAdd);
  return *this;
}

ostream& operator << (ostream& str, SumDerivativeProduct& S)
{
  DerivativeProduct *Summand;
  if (S.Summands.GetNbrElement()>0)
    {
      ListIterator<DerivativeProduct> LI(S.Summands);
      Summand=LI();
      str << (*Summand);
      for (; (Summand=LI())!=NULL;)
	str <<"+" << (*Summand);    
    }
  else str << "0.0";
  str << " ";
  return str;
}
