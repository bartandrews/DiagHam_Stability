////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              general class OrderedList                            //
//                                                                            //
//                        last modification : 07/06/2000                      //
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


#ifndef ORDEREDLIST_H
#define ORDEREDLIST_H


#include "config.h"
#include "GeneralTools/ListElement.h"

#include <iostream>
using std::cout;
using std::endl;

// flag switching debug output
//#define TEST_ORDEREDLIST

template <class ClassName>
class OrderedList
{

public:
  
  ListElement<ClassName> *FirstElement;
  ListElement<ClassName> *CurrentElement;
  int NbrElement;
  bool EliminateDuplicates;
  
  //constructors
  OrderedList(bool eliminateDuplicates=true);
  OrderedList(const ClassName &Element, bool eliminateDuplicates=true);
  OrderedList(const OrderedList &L);
  
  //destructors
  ~OrderedList();	// delete the list
  
  // delete all elements
  void DeleteList();

  //assignment
  OrderedList<ClassName>& operator = (const OrderedList<ClassName>& L);
  
  //add an element to the list
  OrderedList<ClassName>& operator += (const ClassName& Element);
  OrderedList<ClassName>& operator += (const ClassName *Element);

  // add an element to the list, and get its position
  OrderedList<ClassName>& Insert (const ClassName& Element, int &Pos, ClassName* &Duplicate);
  
  //merge two lists
  OrderedList<ClassName>& Merge (const OrderedList<ClassName>& L, int Pos);
  
  //extract a sublist
  OrderedList<ClassName> Extract(const OrderedList<ClassName>& L, int Pos);
  OrderedList<ClassName> Extract(const OrderedList<ClassName>& L, int Pos, int Nbr);
  
  //return an element of the list
  ClassName& operator [] (const int index);
  
  //delete an element of the list
  void DeleteElement (int Pos);
  
  //return the number of element
  int GetNbrElement () {return NbrElement;}

    
};



//constructors

template <class ClassName>
OrderedList<ClassName>::OrderedList(bool eliminateDuplicates)
{
  this->FirstElement = 0;
  this->CurrentElement = 0;
  this->NbrElement = 0;
  this->EliminateDuplicates=eliminateDuplicates;
}

template <class ClassName>
OrderedList<ClassName>::OrderedList(const ClassName &Element, bool eliminateDuplicates)
{
  this->FirstElement = new ListElement<ClassName>(Element);
  this->CurrentElement = this->FirstElement;
  this->NbrElement = 1;
  this->EliminateDuplicates=eliminateDuplicates;
}

template <class ClassName>
OrderedList<ClassName>::OrderedList(const OrderedList<ClassName> &L)
{
  if (L.NbrElement == 0)
    {
      this->FirstElement = 0;
      this->CurrentElement = 0;
      this->NbrElement = 0;
      this->EliminateDuplicates=L.EliminateDuplicates;
   }
  else
    {
      ListElement<ClassName> *Tmp = L.FirstElement;
      this->FirstElement = new ListElement<ClassName>(*Tmp);
      this->CurrentElement = this->FirstElement;      
      while (Tmp->NextPointer != 0)
	{
	  Tmp = Tmp->NextPointer;
	  this->CurrentElement->NextPointer = new ListElement<ClassName>(Tmp->Element, this->CurrentElement);
	  this->CurrentElement = this->CurrentElement->NextPointer;
	}
      this->NbrElement = L.NbrElement;
      this->EliminateDuplicates=L.EliminateDuplicates;
    }
}

//destructors

template <class ClassName>
OrderedList<ClassName>::~OrderedList()
{
  ListElement<ClassName> *ElementPointer = this->FirstElement;
  ListElement<ClassName> *OldElement;
  while (ElementPointer != 0)
    {
      OldElement = ElementPointer;
      ElementPointer = ElementPointer->NextPointer;
      delete OldElement;
    }
}

// delete all elements

template <class ClassName>
void OrderedList<ClassName>::DeleteList()
{
  ListElement<ClassName> *ElementPointer = this->FirstElement;
  ListElement<ClassName> *OldElement;
  while (ElementPointer != 0)
    {
      OldElement = ElementPointer;
      ElementPointer = ElementPointer->NextPointer;
      delete OldElement;
    }
  this->FirstElement = 0;
  this->CurrentElement = 0;
  this->NbrElement = 0;
}

//assignement

template <class ClassName>
OrderedList<ClassName>& OrderedList<ClassName>::operator = (const OrderedList<ClassName>& L)
{
  if (this->NbrElement != 0)
    {
      ListElement<ClassName> *ElementPointer = this->FirstElement;
      ListElement<ClassName> *OldElement;
      while (ElementPointer != 0)
   	{
	  OldElement = ElementPointer;
	  ElementPointer = ElementPointer->NextPointer;
	  delete OldElement;
   	}
    }
  if (L.NbrElement == 0)
    {
      this->FirstElement = 0;
      this->CurrentElement = 0;
      this->NbrElement = 0;
    }
  else
    {
      ListElement<ClassName> *Tmp = L.FirstElement;
      this->FirstElement = new ListElement<ClassName>(*Tmp);
      this->CurrentElement = this->FirstElement;
      while (Tmp->NextPointer != 0)
	{
	  Tmp = Tmp->NextPointer;
	  this->CurrentElement->NextPointer = new ListElement<ClassName>(Tmp->Element, this->CurrentElement);
	  this->CurrentElement = this->CurrentElement->NextPointer;
	}
      this->NbrElement = L.NbrElement;
    }
  return *this;
}

//add an element to the list

template <class ClassName>
inline OrderedList<ClassName>& OrderedList<ClassName>::operator += (const ClassName& Element)
{
  int pos;
  ClassName *duplicate;
  this->Insert(Element, pos, duplicate);
  return *this;
}

template <class ClassName>
inline OrderedList<ClassName>& OrderedList<ClassName>::operator += (const ClassName* Element)
{
  int pos;
  ClassName *duplicate;
  this->Insert(*Element, pos, duplicate);
  return *this;
}



// insert an element in the list
// Pos = returns position where the element was inserted
//
template <class ClassName>
OrderedList<ClassName>& OrderedList<ClassName>::Insert (const ClassName& Element, int &Pos, ClassName* &Duplicate)
{
  if (this->NbrElement !=0)
    {
      int i = 0;
      ListElement<ClassName> *Tmp = this->FirstElement;
      while (((Tmp->Element<Element)||((!this->EliminateDuplicates)&&(Tmp->Element==Element))) && (Tmp->NextPointer != 0))
	{
	  ++i;
#ifdef TEST_ORDEREDLIST
 	  cout << "1 - Found "<<Tmp->Element<< " < " ;
	  cout << (ClassName)Element << " -> proceed to next "<< endl;
#endif
	  Tmp = Tmp->NextPointer;
	}
#ifdef TEST_ORDEREDLIST
      cout<< "2 - Positioning new element "<<Element<<" near " << Tmp->Element << endl;
      if (Tmp->Element<Element)
	{
	  cout << "3 - " << Tmp->Element<<" < " << Element << endl;
	}
      else
	{
	  cout << "4 - " << Tmp->Element<<" >= " << Element << endl;
	}
#endif
      if ((this->EliminateDuplicates)&&(Tmp->Element==Element))
	{
	  Pos=i;
	  Duplicate = &(Tmp->Element);
	  return *this;
	}
      if (Tmp->NextPointer == 0)
	{
	  if (Tmp->Element<Element) // insert after Tmp
	    {
#ifdef TEST_ORDEREDLIST
	      cout << "5a - Next=0 && " << Tmp->Element<<" < " << Element << endl;
	      cout<< "5b - Insert new element "<<Element<<" after existing "<<Tmp->Element<<endl;
#endif
	      Tmp->NextPointer = new ListElement<ClassName>(Element, Tmp); 	      
	      Duplicate = NULL;
	      Pos = i+1;
	      
	    }
	  else // insert before tmp
	    {
#ifdef TEST_ORDEREDLIST
	      cout << "6a - Next=0 && " << Tmp->Element<<" >= " << endl; //Element << endl;
	      cout << "6b - Insert new element "<<Element<<" before existing "<<Tmp->Element<<endl;
#endif
	      Tmp->PreviousPointer = new ListElement<ClassName>(Element, Tmp, Tmp->PreviousPointer);
	      if (Tmp->PreviousPointer->PreviousPointer==0)
		{
#ifdef TEST_ORDEREDLIST
 		  cout << "7 - New first element!"<<endl;
#endif
		  this->FirstElement = Tmp->PreviousPointer;
		}
	      Duplicate = NULL;
	      Pos = i;
	    }
	}
      else
	{
	  if ((Tmp->Element<Element)||(Tmp->Element==Element)) // insert after Tmp
	    {
#ifdef TEST_ORDEREDLIST
	      cout<< "8 - NonZero NextPointer, insert new element "<<Element<<" after existing "<<Tmp->Element<<endl;
#endif
	      Tmp->NextPointer = new ListElement<ClassName>(Element, Tmp->NextPointer, Tmp->PreviousPointer);
	      Duplicate = NULL;
	      Pos = i+1;
	    }
	  else // insert before tmp
	    {
#ifdef TEST_ORDEREDLIST
 	      cout << "9a - NonZero NextPointer, " << Tmp->Element<<" >= " << Element << endl;
	      cout<< "9b - NonZero NextPointer, insert new element "<<Element<<" before existing "<<Tmp->Element<<endl;
#endif
	      Tmp->PreviousPointer = new ListElement<ClassName>(Element, Tmp, Tmp->PreviousPointer);
	      if (Tmp->PreviousPointer->PreviousPointer==0)
		{
#ifdef TEST_ORDEREDLIST
 		  cout << "10 - New first element!"<<endl;
#endif
		  this->FirstElement = Tmp->PreviousPointer;
		}
	      Duplicate = NULL;
	      Pos = i;
	    }
	}
      this->NbrElement++;
    }
  else
    {
      this->CurrentElement = new ListElement<ClassName>(Element);
      this->FirstElement = this->CurrentElement;
      Pos=0;
      Duplicate = NULL;
      this->NbrElement++;
    }
  return *this;
}


//merge two lists

template <class ClassName>
OrderedList<ClassName>& OrderedList<ClassName>::Merge (const OrderedList<ClassName>& L, int Pos)
{
  if (Pos < this->NbrElement)
    {
      ListElement<ClassName> *Tmp = L.FirstElement;
      ListElement<ClassName> *Tmp2 = this->FirstElement;
      int i = 0;
      while ((i++ < (Pos-1)) && (Tmp2 !=0))
	Tmp2 = Tmp2->NextPointer;
      if (Tmp2 !=0)
	{
	  ListElement<ClassName> *Tmp3 = Tmp2->NextPointer;
	  if ((Pos == 0) && (Tmp != 0))
	    {
	      Tmp3 = this->FirstElement;
	      this->FirstElement = new ListElement<ClassName>(Tmp->Element);
	      Tmp2 = this->FirstElement;
	      Tmp = Tmp->NextPointer;
	      this->NbrElement++;
	    }
	  while (Tmp != 0)
	    {
	      Tmp2->NextPointer = new ListElement<ClassName>(Tmp->Element, Tmp2);
	      Tmp2 = Tmp2->NextPointer;
	      Tmp = Tmp->NextPointer;
	      this->NbrElement++;
	    }
	  Tmp2->NextPointer = Tmp3;
	  if (Tmp3 != 0)
	    Tmp3->PreviousPointer = Tmp2;
	  else
	    this->CurrentElement = Tmp2;
	}
    }
  return *this;
}



//extract a sublist

template <class ClassName>
OrderedList<ClassName> Extract(const OrderedList<ClassName>& L, int Pos)
{
  OrderedList<ClassName> TmpList;
  ListElement<ClassName> *Tmp = L.FirstElement;
  int i =0;
  while ((i++ < Pos) && (Tmp != 0))
    Tmp = Tmp->NextPointer;
  while (Tmp != 0)
    {
      TmpList += Tmp->Element;
      Tmp = Tmp->NextPointer;
    }
  return TmpList;
}

template <class ClassName>
OrderedList<ClassName> Extract(const OrderedList<ClassName>& L, int Pos, int Nbr)
{
  OrderedList<ClassName> TmpList;
  ListElement<ClassName> *Tmp = L.FirstElement;
  int i =0;
  while ((i++ < Pos) && (Tmp != 0))
    Tmp = Tmp->NextPointer;
  i =0;
  while ((i++ < Nbr) && (Tmp != 0))
    {
      TmpList += Tmp->Element;
      Tmp = Tmp->NextPointer;
    }
  return TmpList;
}

//return an element of the list

template <class ClassName>
ClassName& OrderedList<ClassName>::operator [] (const int index)
{
  ListElement<ClassName> *ElementPointer = this->FirstElement;
  int i = 0;
  while (i++<index)
    ElementPointer = ElementPointer->NextPointer;
  return ElementPointer->Element;
}

//delete an element of the list

template <class ClassName>
void OrderedList<ClassName>::DeleteElement (int Pos)
{
  if (Pos == 0)
    {
      if (this->FirstElement != 0)
	{
	  if (this->NbrElement == 1)
	    {
	      delete this->FirstElement;
	      this->FirstElement = 0;
	      this->CurrentElement = 0;
	    }
	  else
	    {
	      ListElement<ClassName> *ElementPointer = this->FirstElement->NextPointer;
	      delete this->FirstElement;
	      this->FirstElement = ElementPointer;
	    }
	}
    }
  else
    {
      ListElement<ClassName> *ElementPointer = this->FirstElement;
      int i = 0;
      while (i++<Pos)
	ElementPointer = ElementPointer->NextPointer;
      if (ElementPointer == this->CurrentElement)
      	this->CurrentElement = ElementPointer->PreviousPointer;
      delete ElementPointer;
    }
  this->NbrElement--;
}


#endif

