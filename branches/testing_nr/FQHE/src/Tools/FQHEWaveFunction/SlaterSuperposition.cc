////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2008 Gunnar Moeller                   //
//                                                                            //
//                                                                            //
//           class implementing a single slater determinant form              //
//                                                                            //
//                        last modification : 16/01/2008                      //
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

#include "SlaterSuperposition.h"
#include "GeneralTools/ListIterator.h"
#include <iostream>
#include <cmath>

using std::cout;
using std::endl;
using std::sqrt;

// default constructor
SlaterSuperposition::SlaterSuperposition()
{
  this->Flag.Initialize();
  this->NbrElements=0;
  this->LTotal=0;
  this->LzTotal=0;
}

// constructor from list of elements
// lTotal = twice the total angular momentum
// lzTotal = twice the total Lz component
SlaterSuperposition::SlaterSuperposition(int lTotal, int lzTotal, int nbrElements, SlaterComponent *elements)
{
  this->Flag.Initialize();
  this->NbrElements=0;
  for (int i=0; i<nbrElements; ++i)
    {
      (*this)+=elements[i];
    }
  this->LTotal=lTotal;
  this->LzTotal=lzTotal;
}

// constructor for highest angular momentum eigenstate
SlaterSuperposition::SlaterSuperposition(int nbrParticles, int lzMax)
{
  this->Flag.Initialize();
  this->NbrElements=1;
  this->Elements+=SlaterComponent(nbrParticles, lzMax);
  this->LzTotal=Elements[0].GetTotalLz();
  this->LTotal=LzTotal;
}
  


// constructor beginning with a single term
SlaterSuperposition::SlaterSuperposition(SlaterComponent &element)
{
  this->Flag.Initialize();
  this->NbrElements=1;
  Elements+=element;
  this->LzTotal=Elements[0].GetTotalLz();
  this->LTotal=LzTotal;
}


// copy constructor 
SlaterSuperposition::SlaterSuperposition(const SlaterSuperposition &toCopy)
{
  this->Flag=toCopy.Flag;
  this->NbrElements=toCopy.NbrElements;
  this->Elements=toCopy.Elements;
  this->LzTotal=toCopy.LzTotal;
  this->LTotal=toCopy.LTotal;
}


// destructor
SlaterSuperposition::~SlaterSuperposition()
{  
}


// assignment operator
SlaterSuperposition& SlaterSuperposition::operator = (const SlaterSuperposition& toCopy)
{
  this->Flag=toCopy.Flag;
  this->NbrElements=toCopy.NbrElements;
  this->Elements=toCopy.Elements;
  this->LzTotal=toCopy.LzTotal;
  this->LTotal=toCopy.LTotal;
  return *this;
}

  
// addition of an element
SlaterSuperposition& SlaterSuperposition::operator += (const SlaterComponent& toAdd)
{
  SlaterComponent *component;
  bool added=false;
  for (ListIterator<SlaterComponent> LI(this->Elements); (component=LI())!=0;)
    {
      if ((component!=0)&&(*component==toAdd))
	{
	  *component+=toAdd;
	  added=true;
	}
    }
  if (!added)
    {
      Elements+=toAdd;
      ++NbrElements;
    }
  return *this;
}

// addition of a superposition
SlaterSuperposition& SlaterSuperposition::operator += (const SlaterSuperposition& toAdd)
{
  SlaterComponent *component;
  cout << "before: " << *this<<endl;
  cout << "adding " << toAdd << endl;
  for (ListIterator<SlaterComponent> LI(toAdd.Elements); (component=LI())!=0;)
    (*this)+=*component;
  cout << "after: " << *this<<endl;
  return *this;
}

SlaterSuperposition& SlaterSuperposition::operator *=(double factor)
{
  SlaterComponent *component;
  for (ListIterator<SlaterComponent> LI(this->Elements); (component=LI())!=0;)
    *component *= factor;
  return *this;
}

// accessor methods
void SlaterSuperposition::SetLTotal(int lTotal)
{
  this->LTotal=lTotal;
}

void SlaterSuperposition::SetLzTotal(int lzTotal)
{
  this->LzTotal=lzTotal;
}


SlaterSuperposition SlaterSuperposition::ApplyLMinus()
{
  if (this->LzTotal < -this->LTotal)
    return SlaterSuperposition();
  SlaterSuperposition result;
  SlaterComponent *component;
  for (ListIterator<SlaterComponent> LI(this->Elements); (component=LI())!=0;)
    result+= component->ApplyLMinus();
  result.SetLTotal(this->LTotal);
  result.SetLzTotal(this->LzTotal-2);
  int tmp = this->LzTotal-2;
  result*=1.0/sqrt(0.25*LTotal*(LTotal+2)-0.25*tmp*(tmp-2));
  return result;
}


ostream& operator << (ostream & str, const SlaterSuperposition & s)
{
  SlaterComponent *component;
  ListIterator<SlaterComponent> LI(s.Elements);
  str << "State[L="<<s.LTotal/2.0<<", M="<<s.LzTotal<<"]= ";
  if ((component=LI())!=0)
    str << *component;
  else str << "0";
  for (; (component=LI())!=0;)
    str << " + "<<*component;
  return str;
}
