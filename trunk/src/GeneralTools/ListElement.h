////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   List Elements for general class List                     //
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


#ifndef LISTELEMENT_H
#define LISTELEMENT_H

#include "config.h"

/* template <class ClassName> */
/* ostream& operator << (ostream & Str, const ClassName &Element) */
/* { */
/*   Str << Element; */
/*   return Str;  */
/* } */


template <class ClassName>
class ListElement
{

public:

  ClassName Element;
  ListElement* NextPointer;
  ListElement* PreviousPointer;
  
  //constructors
  ListElement();
  ListElement(const ClassName &Element);
  ListElement(const ClassName &Element, ListElement *PP);
  ListElement(const ClassName &Element, ListElement *PP, ListElement *NP);
  ListElement(ListElement &LElement);

  //destructors
  ~ListElement();
  
  void Exchange(ListElement *LElement);

  //template<class U> friend ostream& operator << (ostream & Str, const typename ListElement<U>::Element &Element);

};



//constructors

template <class ClassName>
ListElement<ClassName>::ListElement()
{
  this->Element = ClassName();
  this->NextPointer = 0;
  this->PreviousPointer = 0;
}

template <class ClassName>
ListElement<ClassName>::ListElement(const ClassName &Element)
{
  this->Element = Element;
  this->NextPointer = 0;
  this->PreviousPointer = 0;
}

template <class ClassName>
ListElement<ClassName>::ListElement(const ClassName &Element, ListElement *PP)
{
  this->Element = Element;
  this->NextPointer = 0;
  this->PreviousPointer = PP;
  if (PP != 0)
    this->PreviousPointer->NextPointer = this;
}

template <class ClassName>
ListElement<ClassName>::ListElement(const ClassName &Element, ListElement *NP, ListElement *PP)
{
  this->Element = Element;
  this->NextPointer = NP;
  this->PreviousPointer = PP;
  if (PP != 0)
    this->PreviousPointer->NextPointer = this;
  if (NP != 0)
    this->NextPointer->PreviousPointer = this;
}

template <class ClassName>
ListElement<ClassName>::ListElement(ListElement<ClassName> &LElement)
{
  this->Element = LElement.Element;
  this->NextPointer = LElement.NextPointer;
  this->PreviousPointer = LElement.PreviousPointer;
}

//destructors

template <class ClassName>
ListElement<ClassName>::~ListElement()
{
  if (this->NextPointer != 0)
    this->NextPointer->PreviousPointer = this->PreviousPointer;
  if (this->PreviousPointer != 0)
    this->PreviousPointer->NextPointer = this->NextPointer;
}

template <class ClassName>
void ListElement<ClassName>::Exchange(ListElement *LElement)
{
  ListElement* TmpNP = this->NextPointer;
  ListElement* TmpPP = this->PreviousPointer;
  if (TmpNP!=0)
    TmpNP->PreviousPointer = LElement;
  if (TmpPP!=0)
    TmpPP->NextPointer = LElement;
  this->NextPointer = LElement->NextPointer;
  this->PreviousPointer = LElement->PreviousPointer;
  if (LElement->NextPointer!=0)
    LElement->NextPointer->PreviousPointer = this;
  if (LElement->NextPointer!=0)
    LElement->PreviousPointer->NextPointer = this;
  LElement->NextPointer = TmpNP;
  LElement->PreviousPointer = TmpPP;
}


/* template <class ClassName> */
/* ostream& operator << (ostream & Str, const typename ListElement<ClassName>::Element &Element) */
/* { */
/*   Str << Element; */
/*   return Str;  */
/* } */

  


#endif
