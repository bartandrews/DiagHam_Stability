////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                           RayTrace version  0.10                           //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class List Iterator                             //
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


#ifndef LISTITERATOR_H
#define LISTITERATOR_H


#include "config.h"
#include "GeneralTools/List.h"


template <class ClassName>
class ListIterator : List<ClassName>
{
public:
  
  //constructors
  ListIterator();
  ListIterator(const List<ClassName>& L);
  
  //destructors
  ~ListIterator() {}
  
  //define list
  void DefineList(const List<ClassName>& L);
  
  //take an element and go to the next one
  ClassName* operator () ();
  
  //go to the next one element
  void Next();

  //go to the previous one element
  void Previous();
  
private:
  
  ListElement<ClassName> *CurrentElement;

};


//constructors

template <class ClassName>
ListIterator<ClassName>::ListIterator()
{
  this->CurrentElement = 0;
}

template <class ClassName>
ListIterator<ClassName>::ListIterator(const List<ClassName>& L)
{
  this->CurrentElement = L.FirstElement;
}

//define list

template <class ClassName>
inline void ListIterator<ClassName>::DefineList(const List<ClassName>& L)
{
  this->CurrentElement = L.FirstElement;
}

//go to the next one element

template <class ClassName>
inline void ListIterator<ClassName>::Next()
{
  if (this->CurrentElement != 0)
    this->CurrentElement = this->CurrentElement->NextPointer;
}

//go to the previous one element

template <class ClassName>
inline void ListIterator<ClassName>::Previous()
{
  if (this->CurrentElement != 0)
    this->CurrentElement = this->CurrentElement->PreviousPointer;
}

template <class ClassName>
inline ClassName* ListIterator<ClassName>::operator () ()
{
  ClassName *TmpE = &(this->CurrentElement->Element);
  if (this->CurrentElement != 0)
    this->CurrentElement = this->CurrentElement->NextPointer;
  return TmpE;
}

#endif

