////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              general class List                            //
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


#ifndef LIST_H
#define LIST_H


#include "config.h"
#include "GeneralTools/ListElement.h"


template <class ClassName>
class List
{

public:
  
  ListElement<ClassName> *FirstElement;
  ListElement<ClassName> *CurrentElement;
  int NbrElement;
  
  //constructors
  List();
  List(const ClassName &Element);
  List(const List &L);
  
  //destructors
  ~List();	// delete the list
  
  // delete all elements
  void DeleteList();

  //assignement
  List<ClassName>& operator = (const List<ClassName>& L);
  
  //add an element to the list
  List<ClassName>& operator += (const ClassName& Element);
  List<ClassName>& operator += (ClassName* Element);
  
  //insert an element in the list
  List<ClassName>& Insert (const ClassName& Element, int Pos);
  List<ClassName>& Insert (ClassName* Element, int Pos);
  
  //link two lists
  List<ClassName>& Link (const List<ClassName>& L);
  
  //merge two lists
  List<ClassName>& Merge (const List<ClassName>& L, int Pos);
  
  //extract a sublist
  List<ClassName> Extract(const List<ClassName>& L, int Pos);
  List<ClassName> Extract(const List<ClassName>& L, int Pos, int Nbr);
  
  //return an element of the list
  ClassName& operator [] (const int index);
  
  //delete an element of the list
  void DeleteElement (int Pos);
  
  //return the number of element
  int GetNbrElement () {return NbrElement;}
  
  //sort list
  List<ClassName>& UpOrder ();
  List<ClassName>& DownOrder ();
  
};


//constructors

template <class ClassName>
List<ClassName>::List()
{
  this->FirstElement = 0;
  this->CurrentElement = 0;
  this->NbrElement = 0;
}

template <class ClassName>
List<ClassName>::List(const ClassName &Element)
{
  this->FirstElement = new ListElement<ClassName>(Element);
  this->CurrentElement = this->FirstElement;
  this->NbrElement = 1;
}

template <class ClassName>
List<ClassName>::List(const List<ClassName> &L)
{
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
}

//destructors

template <class ClassName>
List<ClassName>::~List()
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
void List<ClassName>::DeleteList()
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
List<ClassName>& List<ClassName>::operator = (const List<ClassName>& L)
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
inline List<ClassName>& List<ClassName>::operator += (const ClassName& Element)
{
  this->NbrElement++;
  if (this->CurrentElement !=0)
    this->CurrentElement = new ListElement<ClassName>(Element, this->CurrentElement);
  else
    {
      this->CurrentElement = new ListElement<ClassName>(Element);
      this->FirstElement = this->CurrentElement;
    }
  return *this;
}

template <class ClassName>
inline List<ClassName>& List<ClassName>::operator += (ClassName* Element)
{
  this->NbrElement++;
  if (this->CurrentElement !=0)
    this->CurrentElement = new ListElement<ClassName>(*Element, this->CurrentElement);
  else
    {
      this->CurrentElement = new ListElement<ClassName>(*Element);
      this->CurrentElement = this->FirstElement;
   }
  return *this;
}

//insert an element in the list

template <class ClassName>
List<ClassName>& List<ClassName>::Insert (const ClassName& Element, int Pos)
{
  if ((this->NbrElement !=0) && (Pos !=0))
    {
      int i = 0;
      ListElement<ClassName> *Tmp = this->FirstElement;
      while ((i++ < Pos) && (Tmp->NextPointer != 0))
	Tmp = Tmp->NextPointer;
      if (Tmp->NextPointer == 0)
	this->CurrentElement = new ListElement<ClassName>(Element, this->CurrentElement);
      else
      	Tmp->NextPointer = new ListElement<ClassName>(Element, Tmp->NextPointer, Tmp->PreviousPointer);
      this->NbrElement++;
    }
  else
    {
      if (this->NbrElement == 0)
	{
	  this->CurrentElement = new ListElement<ClassName>(Element);
	  this->FirstElement = this->CurrentElement;
	}
      else
	{
	  this->FirstElement->PreviousPointer =  new ListElement<ClassName>(Element, this->FirstElement, 0);
	  this->FirstElement = this->FirstElement->PreviousPointer;
	}
      this->NbrElement++;
    }
  return *this;
}

//link two lists

template <class ClassName>
List<ClassName>& List<ClassName>::Link (const List<ClassName>& L)
{
  ListElement<ClassName> *Tmp = L.FirstElement;
  if (Tmp == 0)
    return *this;
  if (this->CurrentElement != 0)
    {
      while (Tmp != 0)
	{
	  this->CurrentElement->NextPointer = new ListElement<ClassName>(Tmp->Element, this->CurrentElement);
	  this->CurrentElement = this->CurrentElement->NextPointer;
	  this->NbrElement++;
	  Tmp = Tmp->NextPointer;
	}
      return *this;
    }
  else
    {
      this->CurrentElement = new ListElement<ClassName>(Tmp->Element, this->CurrentElement);
      this->FirstElement = this->CurrentElement;
      this->NbrElement++;
      Tmp = Tmp->NextPointer;
      while (Tmp != 0)
	{
	  this->CurrentElement->NextPointer = new ListElement<ClassName>(Tmp->Element, this->CurrentElement);
	  this->CurrentElement = this->CurrentElement->NextPointer;
	  this->NbrElement++;
	  Tmp = Tmp->NextPointer;
	}
      return *this;
    }
}

//merge two lists

template <class ClassName>
List<ClassName>& List<ClassName>::Merge (const List<ClassName>& L, int Pos)
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
List<ClassName> Extract(const List<ClassName>& L, int Pos)
{
  List<ClassName> TmpList;
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
List<ClassName> Extract(const List<ClassName>& L, int Pos, int Nbr)
{
  List<ClassName> TmpList;
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
ClassName& List<ClassName>::operator [] (const int index)
{
  ListElement<ClassName> *ElementPointer = this->FirstElement;
  int i = 0;
  while (i++<index)
    ElementPointer = ElementPointer->NextPointer;
  return ElementPointer->Element;
}

//delete an element of the list

template <class ClassName>
void List<ClassName>::DeleteElement (int Pos)
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

//sort list

template <class ClassName>
List<ClassName>& List<ClassName>::UpOrder ()
{
  if (this->NbrElement>1)
    for (int i = 1; i < this->NbrElement; i++)
      {
	ListElement<ClassName> *E1 = this->FirstElement;
	ListElement<ClassName> *E2 = this->FirstElement->NextPointer;
      	for (int j = 1; j <= (this->NbrElement-i); j++)
	  {
	    if ((E1->Element) > (E2->Element))
	      {
            	ClassName TmpE = E1->Element;
		E1->Element = E2->Element;
		E2->Element = TmpE;
	      }
            E1 = E2;
            E2 = E2->NextPointer;
	  }
      }
  return *this;
}

template <class ClassName>
List<ClassName>& List<ClassName>::DownOrder ()
{
  if (this->NbrElement>1)
    for (int i = 1; i < this->NbrElement; i++)
      {
	ListElement<ClassName> *E1 = this->FirstElement;
	ListElement<ClassName> *E2 = this->FirstElement->NextPointer;
      	for (int j = 1; j <= (this->NbrElement-i); j++)
	  {
	    if ((E1->Element) < (E2->Element))
	      {
            	ClassName TmpE = E1->Element;
		E1->Element = E2->Element;
		E2->Element = TmpE;
	      }
            E1 = E2;
            E2 = E2->NextPointer;
	  }
      }
  return *this;
}

#endif

