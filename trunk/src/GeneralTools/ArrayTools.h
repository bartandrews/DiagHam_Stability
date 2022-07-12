////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      various generic tools about arrays                    //
//                                                                            //
//                        last modification : 15/09/2003                      //
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

#ifndef ARRAYTOOLS_H
#define ARRAYTOOLS_H

#include "config.h"
#include "GeneralTools/List.h"
#include "MathTools/Complex.h"
#include "MathTools/LongRational.h"
#include "MathTools/FactorialCoefficient.h"


// less operator (A < B) for two elements made of two smaller elements
//
// elementA1 = first part of A
// elementA2 = second part of A
// elementB1 = first part of B
// elementB2 = second part of B
// return = true if A < B

template <class ClassName>
inline bool DoubleElementLessThan (ClassName& elementA1, ClassName& elementA2, ClassName& elementB1, ClassName& elementB2)
{
  return ((elementA1 < elementB1) || ((elementA1 == elementB1) && (elementA2 < elementB2)));
}

// less operator (A < B) for two elements made of three smaller elements
//
// elementA1 = first part of A
// elementA2 = second part of A
// elementA3 = third part of A
// elementB1 = first part of B
// elementB2 = second part of B
// elementB3 = third part of B
// return = true if A < B

template <class ClassName>
inline bool TripleElementLessThan (ClassName& elementA1, ClassName& elementA2, ClassName& elementA3, 
				   ClassName& elementB1, ClassName& elementB2, ClassName& elementB3)
{
  return ((elementA1 < elementB1) || ((elementA1 == elementB1) && 
				      ((elementA2 < elementB2) || ((elementA2 == elementB2) && (elementA3 < elementB3)))));
}

// less operator (A < B) for two elements made of four smaller elements
//
// elementA1 = first part of A
// elementA2 = second part of A
// elementA3 = third part of A
// elementA4 = fourth part of A
// elementB1 = first part of B
// elementB2 = second part of B
// elementB3 = third part of B
// elementB4 = fourth part of B
// return = true if A < B

template <class ClassName>
inline bool QuadElementLessThan (ClassName& elementA1, ClassName& elementA2, ClassName& elementA3, ClassName& elementA4, 
				 ClassName& elementB1, ClassName& elementB2, ClassName& elementB3, ClassName& elementB4)
{
  return ((elementA1 < elementB1) || ((elementA1 == elementB1) && 
				      ((elementA2 < elementB2) || ((elementA2 == elementB2) && ((elementA3 < elementB3) ||
												((elementA3 == elementB3) && (elementA4 < elementB4)))))));
}

// up ordering array sort using quick sort
//
// array = pointer to the array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	  }
	--j;
	ClassName Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	SortArrayUpOrdering(array, i);
	SortArrayUpOrdering(&(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// down ordering array sort using quick sort
//
// array = pointer to the array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0l:
      return;
    case 1l:
      return;
    case 2l:
      {
	if (array[0l] < array[1l])
	  {
	    ClassName TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (array[0l] < array[1l])
	  {
	    TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	  }
	if (array[1l] < array[2l])
	  {
	    TmpElement = array[1l];
	    array[1l] = array[2l];
	    array[2l] = TmpElement;
	  }	
	if (array[0l] < array[1l])
	  {
	    TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	  }	
	return;
      }
      break;
    default:
      {
	long j = nbrValue - 1l;
	long i = nbrValue >> 1;
	ClassName TmpElement;
	if (array[0l] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0l];
	    array[0l] = TmpElement;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	  }
	if (array[0l] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0l];
	    array[0l] = TmpElement;
	  }
	--j;
	ClassName Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	i = 0l;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2l] = array[i];
	array[i] = Pivot;
	SortArrayDownOrdering(array, i);
	SortArrayDownOrdering(&(array[i + 1l]), nbrValue - i - 1l);	
      }
    }
  return;
}

// sort in down ordering an array of doubles and apply the same sort on another array (using quick sort)
//
// doubleArray = pointer to the array of doubles
// array = pointer to the second array to sort in the same way as doubleArray
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(double* doubleArray, ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (doubleArray[0] < doubleArray[1])
	  {
	    double TmpElement = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement;
	    ClassName TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	double TmpElement;
	ClassName TmpElement2;
	if (doubleArray[0] < doubleArray[1])
	  {
	    TmpElement = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	if (doubleArray[1] < doubleArray[2])
	  {
	    TmpElement = doubleArray[1];
	    doubleArray[1] = doubleArray[2];
	    doubleArray[2] = TmpElement;
	    TmpElement2 = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement2;
	  }	
	if (doubleArray[0] < doubleArray[1])
	  {
	    TmpElement = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	double TmpElement;
	ClassName TmpElement2;
	if (doubleArray[0] <  doubleArray[i])
	  {
	    TmpElement = doubleArray[i];
	    doubleArray[i] = doubleArray[0];
	    doubleArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	if (doubleArray[i] <  doubleArray[j])
	  {
	    TmpElement = doubleArray[i];
	    doubleArray[i] = doubleArray[j];
	    doubleArray[j] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement2;
	  }
	if (doubleArray[0] <  doubleArray[i])
	  {
	    TmpElement = doubleArray[i];
	    doubleArray[i] = doubleArray[0];
	    doubleArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	--j;
	double Pivot = doubleArray[i];
	ClassName Pivot2 = array[i];
	doubleArray[i] = doubleArray[j];
	doubleArray[j] = Pivot;
	array[i] = array[j];
	array[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (doubleArray[++i] > Pivot);
	    while (doubleArray[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = doubleArray[i];
		doubleArray[i] = doubleArray[j];
		doubleArray[j] = TmpElement;	    
		TmpElement2 = array[i];
		array[i] = array[j];
		array[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	doubleArray[nbrValue - 2] = doubleArray[i];
	doubleArray[i] = Pivot;
	array[nbrValue - 2] = array[i];
	array[i] = Pivot2;
	SortArrayDownOrdering(doubleArray, array, i);
	SortArrayDownOrdering(&(doubleArray[i + 1]), &(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in down ordering an arrayand apply the same sort on another array of integer (using quick sort)
//
// array = pointer to the array
// integerArray = pointer to the second array of integers to sort in the same way as array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(ClassName* array, int* integerArray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0l:
      return;
    case 1l:
      return;
    case 2l:
      {
	if (array[0l] < array[1l])
	  {
	    ClassName TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	    int TmpElement2 = integerArray[0l];
	    integerArray[0l] = integerArray[1l];
	    integerArray[1l] = TmpElement2;
	  }
	return;
      }
      break;
    case 3l:
      {
	ClassName TmpElement;
	int TmpElement2;
	if (array[0l] < array[1l])
	  {
	    TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	    TmpElement2 = integerArray[0l];
	    integerArray[0l] = integerArray[1l];
	    integerArray[1l] = TmpElement2;
	  }
	if (array[1l] < array[2l])
	  {
	    TmpElement = array[1l];
	    array[1l] = array[2l];
	    array[2l] = TmpElement;
	    TmpElement2 = integerArray[1l];
	    integerArray[1l] = integerArray[2l];
	    integerArray[2l] = TmpElement2;
	  }	
	if (array[0l] < array[1l])
	  {
	    TmpElement = array[0l];
	    array[0l] = array[1l];
	    array[1l] = TmpElement;
	    TmpElement2 = integerArray[0l];
	    integerArray[0l] = integerArray[1l];
	    integerArray[1l] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	long j = nbrValue - 1l;
	long i = nbrValue >> 1;
	ClassName TmpElement;
	int TmpElement2;
	if (array[0l] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0l];
	    array[0l] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[0l];
	    integerArray[0l] = TmpElement2;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[j];
	    integerArray[j] = TmpElement2;
	  }
	if (array[0l] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0l];
	    array[0l] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[0l];
	    integerArray[0l] = TmpElement2;
	  }
	--j;
	ClassName Pivot = array[i];
	int Pivot2 = integerArray[i];
	array[i] = array[j];
	array[j] = Pivot;
	integerArray[i] = integerArray[j];
	integerArray[j] = Pivot2;	    
	i = 0l;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = integerArray[i];
		integerArray[i] = integerArray[j];
		integerArray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2l] = array[i];
	array[i] = Pivot;
	integerArray[nbrValue - 2l] = integerArray[i];
	integerArray[i] = Pivot2;
	SortArrayDownOrdering(array, integerArray, i);
	SortArrayDownOrdering(&(array[i + 1l]), &(integerArray[i + 1l]), nbrValue - i - 1l);	
      }
    }
  return;
}

// sort in down ordering an arrayand apply the same sort on two other arrays of integer (using quick sort)
//
// array = pointer to the array
// integerArray1 = pointer to the second array of integers to sort in the same way as array
// integerArray2 = pointer to the third array of integers to sort in the same way as array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(ClassName* array, int* integerArray1, int* integerArray2, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] < array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    int TmpElement2 = integerArray1[0];
	    integerArray1[0] = integerArray1[1];
	    integerArray1[1] = TmpElement2;
	    TmpElement2 = integerArray2[0];
	    integerArray2[0] = integerArray2[1];
	    integerArray2[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	int TmpElement2;
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = integerArray1[0];
	    integerArray1[0] = integerArray1[1];
	    integerArray1[1] = TmpElement2;
	    TmpElement2 = integerArray2[0];
	    integerArray2[0] = integerArray2[1];
	    integerArray2[1] = TmpElement2;
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = integerArray1[1];
	    integerArray1[1] = integerArray1[2];
	    integerArray1[2] = TmpElement2;
	    TmpElement2 = integerArray2[1];
	    integerArray2[1] = integerArray2[2];
	    integerArray2[2] = TmpElement2;
	  }	
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = integerArray1[0];
	    integerArray1[0] = integerArray1[1];
	    integerArray1[1] = TmpElement2;
	    TmpElement2 = integerArray2[0];
	    integerArray2[0] = integerArray2[1];
	    integerArray2[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	int TmpElement2;
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = integerArray1[i];
	    integerArray1[i] = integerArray1[0];
	    integerArray1[0] = TmpElement2;
	    TmpElement2 = integerArray2[i];
	    integerArray2[i] = integerArray2[0];
	    integerArray2[0] = TmpElement2;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = integerArray1[i];
	    integerArray1[i] = integerArray1[j];
	    integerArray1[j] = TmpElement2;
	    TmpElement2 = integerArray2[i];
	    integerArray2[i] = integerArray2[j];
	    integerArray2[j] = TmpElement2;
	  }
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = integerArray1[i];
	    integerArray1[i] = integerArray1[0];
	    integerArray1[0] = TmpElement2;
	    TmpElement2 = integerArray2[i];
	    integerArray2[i] = integerArray2[0];
	    integerArray2[0] = TmpElement2;
	  }
	--j;
	ClassName Pivot = array[i];
	int Pivot2 = integerArray1[i];
	int Pivot3 = integerArray2[i];
	array[i] = array[j];
	array[j] = Pivot;
	integerArray1[i] = integerArray1[j];
	integerArray1[j] = Pivot2;	    
	integerArray2[i] = integerArray2[j];
	integerArray2[j] = Pivot3;	    
	i = 0;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = integerArray1[i];
		integerArray1[i] = integerArray1[j];
		integerArray1[j] = TmpElement2;	    
		TmpElement2 = integerArray2[i];
		integerArray2[i] = integerArray2[j];
		integerArray2[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	integerArray1[nbrValue - 2] = integerArray1[i];
	integerArray1[i] = Pivot2;
	integerArray2[nbrValue - 2] = integerArray2[i];
	integerArray2[i] = Pivot3;
	SortArrayDownOrdering(array, integerArray1, integerArray2, i);
	SortArrayDownOrdering(&(array[i + 1]), &(integerArray1[i + 1]), &(integerArray2[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in down ordering an arrayand apply the same sort on another array of integer (using quick sort)
//
// array = pointer to the array
// integerArray = pointer to the second array of integers to sort in the same way as array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(ClassName* array, long* integerArray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] < array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    long TmpElement2 = integerArray[0];
	    integerArray[0] = integerArray[1];
	    integerArray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	long TmpElement2;
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = integerArray[0];
	    integerArray[0] = integerArray[1];
	    integerArray[1] = TmpElement2;
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = integerArray[1];
	    integerArray[1] = integerArray[2];
	    integerArray[2] = TmpElement2;
	  }	
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = integerArray[0];
	    integerArray[0] = integerArray[1];
	    integerArray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	long j = nbrValue - 1;
	long i = nbrValue >> 1;
	ClassName TmpElement;
	long TmpElement2;
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[0];
	    integerArray[0] = TmpElement2;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[j];
	    integerArray[j] = TmpElement2;
	  }
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = integerArray[i];
	    integerArray[i] = integerArray[0];
	    integerArray[0] = TmpElement2;
	  }
	--j;
	ClassName Pivot = array[i];
	long Pivot2 = integerArray[i];
	array[i] = array[j];
	array[j] = Pivot;
	integerArray[i] = integerArray[j];
	integerArray[j] = Pivot2;	    
	i = 0;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = integerArray[i];
		integerArray[i] = integerArray[j];
		integerArray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	integerArray[nbrValue - 2] = integerArray[i];
	integerArray[i] = Pivot2;
	SortArrayDownOrdering(array, integerArray, i);
	SortArrayDownOrdering(&(array[i + 1]), &(integerArray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in down ordering an array and apply the same sort on another array of double (using quick sort)
//
// array = pointer to the array
// integerArray = pointer to the second array of integers to sort in the same way as array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(ClassName* array, double* doubleArray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] < array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    double TmpElement2 = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	double TmpElement2;
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement2;
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = doubleArray[1];
	    doubleArray[1] = doubleArray[2];
	    doubleArray[2] = TmpElement2;
	  }	
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = doubleArray[0];
	    doubleArray[0] = doubleArray[1];
	    doubleArray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	double TmpElement2;
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = doubleArray[i];
	    doubleArray[i] = doubleArray[0];
	    doubleArray[0] = TmpElement2;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = doubleArray[i];
	    doubleArray[i] = doubleArray[j];
	    doubleArray[j] = TmpElement2;
	  }
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = doubleArray[i];
	    doubleArray[i] = doubleArray[0];
	    doubleArray[0] = TmpElement2;
	  }
	--j;
	ClassName Pivot = array[i];
	double Pivot2 = doubleArray[i];
	array[i] = array[j];
	array[j] = Pivot;
	doubleArray[i] = doubleArray[j];
	doubleArray[j] = Pivot2;	    
	i = 0;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = doubleArray[i];
		doubleArray[i] = doubleArray[j];
		doubleArray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	doubleArray[nbrValue - 2] = doubleArray[i];
	doubleArray[i] = Pivot2;
	SortArrayDownOrdering(array, doubleArray, i);
	SortArrayDownOrdering(&(array[i + 1]), &(doubleArray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in down ordering an array of complex numbers and apply the same sort on another array (using quick sort)
//
// complexArray = pointer to the array of complex numbers
// array = pointer to the second array to sort in the same way as complexArray
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayDownOrdering(Complex* complexArray, ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (complexArray[0] < complexArray[1])
	  {
	    Complex TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    ClassName TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	Complex TmpElement;
	ClassName TmpElement2;
	if (complexArray[0] < complexArray[1])
	  {
	    TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	if (complexArray[1] < complexArray[2])
	  {
	    TmpElement = complexArray[1];
	    complexArray[1] = complexArray[2];
	    complexArray[2] = TmpElement;
	    TmpElement2 = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement2;
	  }	
	if (complexArray[0] < complexArray[1])
	  {
	    TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	Complex TmpElement;
	ClassName TmpElement2;
	if (complexArray[0] <  complexArray[i])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[0];
	    complexArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	if (complexArray[i] <  complexArray[j])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[j];
	    complexArray[j] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement2;
	  }
	if (complexArray[0] <  complexArray[i])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[0];
	    complexArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	--j;
	Complex Pivot = complexArray[i];
	ClassName Pivot2 = array[i];
	complexArray[i] = complexArray[j];
	complexArray[j] = Pivot;
	array[i] = array[j];
	array[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (complexArray[++i] > Pivot);
	    while (complexArray[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = complexArray[i];
		complexArray[i] = complexArray[j];
		complexArray[j] = TmpElement;	    
		TmpElement2 = array[i];
		array[i] = array[j];
		array[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	complexArray[nbrValue - 2] = complexArray[i];
	complexArray[i] = Pivot;
	array[nbrValue - 2] = array[i];
	array[i] = Pivot2;
	SortArrayDownOrdering(complexArray, array, i);
	SortArrayDownOrdering(&(complexArray[i + 1]), &(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in up ordering an array of complex numbers and apply the same sort on another array (using quick sort)
//
// complexArray = pointer to the array of complex numbers
// array = pointer to the second array to sort in the same way as complexArray
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(Complex* complexArray, ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (complexArray[0] > complexArray[1])
	  {
	    Complex TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    ClassName TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	Complex TmpElement;
	ClassName TmpElement2;
	if (complexArray[0] > complexArray[1])
	  {
	    TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	if (complexArray[1] > complexArray[2])
	  {
	    TmpElement = complexArray[1];
	    complexArray[1] = complexArray[2];
	    complexArray[2] = TmpElement;
	    TmpElement2 = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement2;
	  }	
	if (complexArray[0] > complexArray[1])
	  {
	    TmpElement = complexArray[0];
	    complexArray[0] = complexArray[1];
	    complexArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	Complex TmpElement;
	ClassName TmpElement2;
	if (complexArray[0] >  complexArray[i])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[0];
	    complexArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	if (complexArray[i] >  complexArray[j])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[j];
	    complexArray[j] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement2;
	  }
	if (complexArray[0] >  complexArray[i])
	  {
	    TmpElement = complexArray[i];
	    complexArray[i] = complexArray[0];
	    complexArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	--j;
	Complex Pivot = complexArray[i];
	ClassName Pivot2 = array[i];
	complexArray[i] = complexArray[j];
	complexArray[j] = Pivot;
	array[i] = array[j];
	array[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (complexArray[++i] < Pivot);
	    while (complexArray[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = complexArray[i];
		complexArray[i] = complexArray[j];
		complexArray[j] = TmpElement;	    
		TmpElement2 = array[i];
		array[i] = array[j];
		array[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	complexArray[nbrValue - 2] = complexArray[i];
	complexArray[i] = Pivot;
	array[nbrValue - 2] = array[i];
	array[i] = Pivot2;
	SortArrayDownOrdering(complexArray, array, i);
	SortArrayDownOrdering(&(complexArray[i + 1]), &(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// sort in up ordering an array of long int numbers and apply the same sort on another array (using quick sort)
//
// longArray = pointer to the array of long int numbers
// array = pointer to the second array to sort in the same way as longArray
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(unsigned long* longArray, ClassName* array, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (longArray[0] > longArray[1])
	  {
	    unsigned long TmpElement = longArray[0];
	    longArray[0] = longArray[1];
	    longArray[1] = TmpElement;
	    ClassName TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	unsigned long TmpElement;
	ClassName TmpElement2;
	if (longArray[0] > longArray[1])
	  {
	    TmpElement = longArray[0];
	    longArray[0] = longArray[1];
	    longArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }
	if (longArray[1] > longArray[2])
	  {
	    TmpElement = longArray[1];
	    longArray[1] = longArray[2];
	    longArray[2] = TmpElement;
	    TmpElement2 = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement2;
	  }	
	if (longArray[0] > longArray[1])
	  {
	    TmpElement = longArray[0];
	    longArray[0] = longArray[1];
	    longArray[1] = TmpElement;
	    TmpElement2 = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	unsigned long TmpElement;
	ClassName TmpElement2;
	if (longArray[0] >  longArray[i])
	  {
	    TmpElement = longArray[i];
	    longArray[i] = longArray[0];
	    longArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	if (longArray[i] >  longArray[j])
	  {
	    TmpElement = longArray[i];
	    longArray[i] = longArray[j];
	    longArray[j] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement2;
	  }
	if (longArray[0] >  longArray[i])
	  {
	    TmpElement = longArray[i];
	    longArray[i] = longArray[0];
	    longArray[0] = TmpElement;
	    TmpElement2 = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement2;
	  }
	--j;
	unsigned long Pivot = longArray[i];
	ClassName Pivot2 = array[i];
	longArray[i] = longArray[j];
	longArray[j] = Pivot;
	array[i] = array[j];
	array[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (longArray[++i] < Pivot);
	    while (longArray[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = longArray[i];
		longArray[i] = longArray[j];
		longArray[j] = TmpElement;	    
		TmpElement2 = array[i];
		array[i] = array[j];
		array[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	longArray[nbrValue - 2] = longArray[i];
	longArray[i] = Pivot;
	array[nbrValue - 2] = array[i];
	array[i] = Pivot2;
	SortArrayUpOrdering(longArray, array, i);
	SortArrayUpOrdering(&(longArray[i + 1]), &(array[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array sort using quick sort, and sort row of an 2d-array in the way
//
// array = pointer to the array
// coarray = pointer to the 2d-array to sort in the way as array
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(ClassName* array, ClassName** coarray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    ClassName* TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	ClassName* TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = coarray[1];
	    coarray[1] = coarray[2];
	    coarray[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	ClassName* TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[j];
	    coarray[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	--j;
	ClassName Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	ClassName* Pivot2 = coarray[i];
	coarray[i] = coarray[j];
	coarray[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = coarray[i];
		coarray[i] = coarray[j];
		coarray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray[nbrValue - 2] = coarray[i];
	coarray[i] = Pivot2;
	SortArrayUpOrdering(array, coarray, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array of integers using quick sort, and sort element of a second array in the same manner
//
// array = array of integers
// coarray = second array to sort
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(int* array, ClassName* coarray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    int TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    ClassName TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	int TmpElement;
	ClassName TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = coarray[1];
	    coarray[1] = coarray[2];
	    coarray[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	int TmpElement;
	ClassName TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[j];
	    coarray[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	--j;
	int Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	ClassName Pivot2 = coarray[i];
	coarray[i] = coarray[j];
	coarray[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = coarray[i];
		coarray[i] = coarray[j];
		coarray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray[nbrValue - 2] = coarray[i];
	coarray[i] = Pivot2;
	SortArrayUpOrdering(array, coarray, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array of ULONGLONG using quick sort, and sort element of a second array in the same manner
//
// array = array of ULONGLONG
// coarray = second array to sort
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(ULONGLONG* array, ClassName* coarray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    ULONGLONG TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    ClassName TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	ULONGLONG TmpElement;
	ClassName TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = coarray[1];
	    coarray[1] = coarray[2];
	    coarray[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ULONGLONG TmpElement;
	ClassName TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[j];
	    coarray[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	--j;
	ULONGLONG Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	ClassName Pivot2 = coarray[i];
	coarray[i] = coarray[j];
	coarray[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = coarray[i];
		coarray[i] = coarray[j];
		coarray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray[nbrValue - 2] = coarray[i];
	coarray[i] = Pivot2;
	SortArrayUpOrdering(array, coarray, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array of doubles using quick sort, and sort element of a second array in the same manner
//
// array = array of integers
// coarray = second array to sort
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(double* array, ClassName* coarray, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    double TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    ClassName TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	double TmpElement;
	ClassName TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = coarray[1];
	    coarray[1] = coarray[2];
	    coarray[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray[0];
	    coarray[0] = coarray[1];
	    coarray[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	double TmpElement;
	ClassName TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[j];
	    coarray[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray[i];
	    coarray[i] = coarray[0];
	    coarray[0] = TmpElement2;
	  }
	--j;
	double Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	ClassName Pivot2 = coarray[i];
	coarray[i] = coarray[j];
	coarray[j] = Pivot2;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = coarray[i];
		coarray[i] = coarray[j];
		coarray[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray[nbrValue - 2] = coarray[i];
	coarray[i] = Pivot2;
	SortArrayUpOrdering(array, coarray, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array of doubles using quick sort, and sort element of two secondary arrays in the same manner
//
// array = array of integers
// coarray1 = first secondary array to sort
// coarray2 = second secondary array to sort
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(double* array, ClassName* coarray1, ClassName* coarray2, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    double TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    ClassName TmpElement2 = coarray1[0];
	    coarray1[0] = coarray1[1]; 
	    coarray1[1] = TmpElement2;
	    TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	double TmpElement;
	ClassName TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray1[0];
	    coarray1[0] = coarray1[1];
	    coarray1[1] = TmpElement2;
	    TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement2 = coarray1[1];
	    coarray1[1] = coarray1[2];
	    coarray1[2] = TmpElement2;
	    TmpElement2 = coarray2[1];
	    coarray2[1] = coarray2[2];
	    coarray2[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement2 = coarray1[0];
	    coarray1[0] = coarray1[1];
	    coarray1[1] = TmpElement2;
	    TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	double TmpElement;
	ClassName TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray1[i];
	    coarray1[i] = coarray1[0];
	    coarray1[0] = TmpElement2;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[0];
	    coarray2[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement2 = coarray1[i];
	    coarray1[i] = coarray1[j];
	    coarray1[j] = TmpElement2;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[j];
	    coarray2[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement2 = coarray1[i];
	    coarray1[i] = coarray1[0];
	    coarray1[0] = TmpElement2;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[0];
	    coarray2[0] = TmpElement2;
	  }
	--j;
	double Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	ClassName Pivot2 = coarray1[i];
	coarray1[i] = coarray1[j];
	coarray1[j] = Pivot2;
	ClassName Pivot3 = coarray2[i];
	coarray2[i] = coarray2[j];
	coarray2[j] = Pivot3;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement2 = coarray1[i];
		coarray1[i] = coarray1[j];
		coarray1[j] = TmpElement2;	    
		TmpElement2 = coarray2[i];
		coarray2[i] = coarray2[j];
		coarray2[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray1[nbrValue - 2] = coarray1[i];
	coarray1[i] = Pivot2;
	coarray2[nbrValue - 2] = coarray2[i];
	coarray2[i] = Pivot3;
	SortArrayUpOrdering(array, coarray1, coarray2, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray1[i + 1]), &(coarray2[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// up ordering array of doubles using quick sort, and sort element of two secondary arrays in the same manner
//
// array = array of integers
// coarray1 = first secondary array to sort
// coarray2 = second secondary array to sort
// nbrValue = nbr of value in the array

template <class ClassName>
void SortArrayUpOrdering(int* array, int* coarray1, ClassName* coarray2, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    int TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement = coarray1[0];
	    coarray1[0] = coarray1[1]; 
	    coarray1[1] = TmpElement;
	    ClassName TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }
	return;
      }
      break;
    case 3:
      {
	int TmpElement;
	ClassName TmpElement2;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement = coarray1[0];
	    coarray1[0] = coarray1[1];
	    coarray1[1] = TmpElement;
	    TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    TmpElement = coarray1[1];
	    coarray1[1] = coarray1[2];
	    coarray1[2] = TmpElement;
	    TmpElement2 = coarray2[1];
	    coarray2[1] = coarray2[2];
	    coarray2[2] = TmpElement2;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    TmpElement = coarray1[0];
	    coarray1[0] = coarray1[1];
	    coarray1[1] = TmpElement;
	    TmpElement2 = coarray2[0];
	    coarray2[0] = coarray2[1];
	    coarray2[1] = TmpElement2;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	int TmpElement;
	ClassName TmpElement2;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement = coarray1[i];
	    coarray1[i] = coarray1[0];
	    coarray1[0] = TmpElement;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[0];
	    coarray2[0] = TmpElement2;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    TmpElement = coarray1[i];
	    coarray1[i] = coarray1[j];
	    coarray1[j] = TmpElement;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[j];
	    coarray2[j] = TmpElement2;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    TmpElement = coarray1[i];
	    coarray1[i] = coarray1[0];
	    coarray1[0] = TmpElement;
	    TmpElement2 = coarray2[i];
	    coarray2[i] = coarray2[0];
	    coarray2[0] = TmpElement2;
	  }
	--j;
	int Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	int Pivot2 = coarray1[i];
	coarray1[i] = coarray1[j];
	coarray1[j] = Pivot2;
	ClassName Pivot3 = coarray2[i];
	coarray2[i] = coarray2[j];
	coarray2[j] = Pivot3;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		TmpElement = coarray1[i];
		coarray1[i] = coarray1[j];
		coarray1[j] = TmpElement;	    
		TmpElement2 = coarray2[i];
		coarray2[i] = coarray2[j];
		coarray2[j] = TmpElement2;	    
	      }
	    else
	      break;
	  }	
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	coarray1[nbrValue - 2] = coarray1[i];
	coarray1[i] = Pivot2;
	coarray2[nbrValue - 2] = coarray2[i];
	coarray2[i] = Pivot3;
	SortArrayUpOrdering(array, coarray1, coarray2, i);
	SortArrayUpOrdering(&(array[i + 1]), &(coarray1[i + 1]), &(coarray2[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}

// merge a list of arrays in a smart way i.e. using as less as possible temporary storage 
//
// arrayList = list of arrays to merge (list will be deleted after use and memory associate to each array free)
// arraySizeList = list of array sizes (list will be deleted after use)
// return value = array corresponding to the merge

template <class ClassName>
ClassName* SmartMergeArrayListIntoArray(List<ClassName*>& arrayList, List<long>& arraySizeList)
{
  List<ClassName*> TmpArrayList;
  List<long> TmpArraySize;
  int NbrArray = arrayList.GetNbrElement();
  if (NbrArray == 0)
    return 0;
  if (NbrArray == 1)
    {
      ClassName* TmpArray = arrayList[0];
      arrayList.DeleteList();
      arraySizeList.DeleteList();
      return TmpArray;
    }
  int Pos = 0;
  while (NbrArray > 1)
    { 
      ClassName* TmpArray1 = arrayList[Pos];
      ClassName* TmpArray2 = arrayList[Pos + 1];
      long TmpArraySize1 = arraySizeList[Pos];
      long TmpArraySize2 = arraySizeList[Pos + 1];      
      TmpArraySize += (TmpArraySize1 + TmpArraySize2);     
      ClassName* TmpArrayMerge = new ClassName [TmpArraySize1 + TmpArraySize2];
      long i = 0;
      for (; i < TmpArraySize1; ++i)
	TmpArrayMerge[i] = TmpArray1[i];
      for (TmpArraySize1 = 0; TmpArraySize1 < TmpArraySize2; ++TmpArraySize1)
	TmpArrayMerge[i++] = TmpArray2[TmpArraySize1];  
      delete[] TmpArray1;
      delete[] TmpArray2;
      TmpArrayList += TmpArrayMerge;
      Pos += 2;
      NbrArray -= 2;
    }
  if (NbrArray == 1)
    {
      TmpArrayList += arrayList[Pos];
      TmpArraySize += arraySizeList[Pos];
    }
  arrayList.DeleteList();
  arraySizeList.DeleteList();
  ClassName* TmpArray = SmartMergeArrayListIntoArray(TmpArrayList, TmpArraySize);
  return TmpArray;
}

// down ordering array sort using quick sort, counting the number of permutations and apply the same sort on another array
//
// array = pointer to the array
// ulArray = pointer to the unsigned long second array to sort in the same way as array
// nbrValue = nbr of value in the array
// nbrPermutation = reference on the number of permutation that have to be done

template <class ClassName>
void SortArrayDownOrdering(ClassName* array,unsigned long* ulArray, long nbrValue, int& nbrPermutation)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] < array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    unsigned long TmpElement2=ulArray[0];
	    ulArray[0] = ulArray[1];
	    ulArray[1] = TmpElement2;
	    nbrPermutation++;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	unsigned long TmpElement2;
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[0];
	    ulArray[0] = ulArray[1];
	    ulArray[1] = TmpElement2;
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[1];
	    ulArray[1] = ulArray[2];
	    ulArray[2] = TmpElement2;
	  }
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[0];
	    ulArray[0] = ulArray[1];
	    ulArray[1] = TmpElement2;
	  }
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	unsigned long TmpElement2;
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[i];
	    ulArray[i] = ulArray[0];
	    ulArray[0] = TmpElement2;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[i];
	    ulArray[i] = ulArray[j];
	    ulArray[j] = TmpElement2;
	  }
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	    TmpElement2=ulArray[i];
	    ulArray[i] = ulArray[0];
	    ulArray[0] = TmpElement2;
	  }
	--j;
	ClassName Pivot;
	unsigned long Pivot2;
	if((i==j)||(array[i]==array[j]))
	  {
	    Pivot = array[i];
	    Pivot2=ulArray[i];
	  }
	else
	  {
	    Pivot = array[i];
	    array[i] = array[j];
	    array[j] = Pivot;
	    nbrPermutation++;
	    Pivot2=ulArray[i];
	    ulArray[i] = ulArray[j];
	    ulArray[j] = Pivot2;
	  }
	i = 0;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;
		nbrPermutation++;
		TmpElement2=ulArray[i];
		ulArray[i] = ulArray[j];
		ulArray[j] = TmpElement2;
	      }
	    else
	      break;
	  }
	if(((nbrValue - 2)==i)||(array[nbrValue - 2]==array[i]))
	  ;
	else
	  {
	    array[nbrValue - 2] = array[i];
	    array[i] = Pivot;
	    nbrPermutation++;
	    ulArray[nbrValue - 2] = ulArray[i];
	    ulArray[i] = Pivot2;
	  }
	SortArrayDownOrdering(array,ulArray, i,nbrPermutation);
	SortArrayDownOrdering(&(array[i + 1]),&(ulArray[i + 1]), nbrValue - i - 1,nbrPermutation);
      }
    }
  return;
}

// down ordering array sort using quick sort, counting the number of permutations
//
// array = pointer to the array
// nbrValue = nbr of value in the array
// nbrPermutation = reference on the number of permutation that have to be done (it has to be set to zero at the first call)

template <class ClassName>
void SortArrayDownOrderingPermutation(ClassName* array, long nbrValue, int& nbrPermutation)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] < array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[0] < array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	  }
	--j;
	ClassName Pivot;
	if((i==j)||(array[i]==array[j]))
	  {
	    Pivot = array[i];
	  }
	else
	  {
	    Pivot = array[i];
	    array[i] = array[j];
	    array[j] = Pivot;
	    nbrPermutation++;
	  }
	i = 0;
	while (true)
	  {
	    while (array[++i] > Pivot);
	    while (array[--j] < Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;
		nbrPermutation++;
	      }
	    else
	      break;
	  }
	if(((nbrValue - 2)==i)||(array[nbrValue - 2]==array[i]))
	  ;
	else
	  {
	    array[nbrValue - 2] = array[i];
	    array[i] = Pivot;
	    nbrPermutation++;
	  }
	SortArrayDownOrderingPermutation(array, i,nbrPermutation);
	SortArrayDownOrderingPermutation(&(array[i + 1]), nbrValue - i - 1,nbrPermutation);
      }
    }
  return;
}

// down ordering array sort using bubble sort, counting the number of permutations of neighboring elements
//
// array = pointer to the array
// nbrValue = nbr of value in the array
// nbrPermutation = reference on the number of permutation that have to be done (it has to be set to zero at the first call)

template <class ClassName>
void SortArrayDownOrderingPermutationBubbleSort(ClassName* array, long nbrValue, int& nbrPermutation)
{
  long ReducedNbrValues = nbrValue - 1;
  ClassName TmpElement;
  for (long i = 0; i < ReducedNbrValues; ++i)
    {
      for (long j = ReducedNbrValues; j > i; --j)
	{
	  if (array[j] > array[j - 1])
	    {
	      TmpElement = array[j - 1];
	      array[j - 1] = array[j];
	      array[j] = TmpElement;
	      nbrPermutation++;
	    }
	}
    }
}

// up ordering array sort using quick sort, counting the number of permutations
//
// array = pointer to the array
// nbrValue = nbr of value in the array
// nbrPermutation = reference on the number of permutation that have to be done (it has to be set to zero at the first call)

template <class ClassName>
void SortArrayUpOrderingPermutation(ClassName* array, long nbrValue, int& nbrPermutation)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (array[0] > array[1])
	  {
	    ClassName TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[1] > array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	    nbrPermutation++;
	  }	
	if (array[0] > array[1])
	  {
	    TmpElement = array[0];
	    array[0] = array[1];
	    array[1] = TmpElement;
	    nbrPermutation++;
	  }	
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[i] >  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	    nbrPermutation++;
	  }
	if (array[0] >  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	    nbrPermutation++;
	  }
	--j;
	ClassName Pivot = array[i];
	array[i] = array[j];
	array[j] = Pivot;
	i = 0;
	while (true)
	  {
	    while (array[++i] < Pivot);
	    while (array[--j] > Pivot);
	    if (i < j)
	      {
		TmpElement = array[i];
		array[i] = array[j];
		array[j] = TmpElement;	    
		nbrPermutation++;
	      }
	    else
	      break;
	  }	
	if(((nbrValue - 2)==i)||(array[nbrValue - 2]==array[i]))
	  ;
	else
	  {
	    array[nbrValue - 2] = array[i];
	    array[i] = Pivot;
	    nbrPermutation++;
	  }
	SortArrayUpOrderingPermutation(array, i, nbrPermutation);
	SortArrayUpOrderingPermutation(&(array[i + 1]), nbrValue - i - 1, nbrPermutation);	
      }
    }
  return;
}


// find an element in a sorted array (from smallest to largest)
//
// element = element to find
// array = sorted array where to search 
// nbrValue = number of values in array
// return value = element position (negative if not in the array)

template <class ClassName>
int SearchInArray(const ClassName& element, ClassName* array, int nbrValue)
{

  int StartIndex = 0;
  int EndIndex = nbrValue - 1;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {       
      MidIndex = (StartIndex + EndIndex) >> 1;       
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
   
  if (array[StartIndex] == element)
    return StartIndex;
  if (array[EndIndex] == element)
    return EndIndex;
  return -1;
}

// find an element in a sorted array (from largest to smallest)
//
// element = element to find
// array = sorted array where to search 
// nbrValue = number of values in array
// return value = element position (nbrValue if not in the array)

template <class ClassName>
int SearchInArrayDownOrdering(const ClassName& element, ClassName* array, int nbrValue)
{

  int StartIndex = 0;
  int EndIndex = nbrValue - 1;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {       
      MidIndex = (StartIndex + EndIndex) >> 1;       
      if(array[MidIndex] > element)
	StartIndex = MidIndex;
      else
	EndIndex = MidIndex;
    }
   
  if (array[StartIndex] == element)
    return StartIndex;
  if (array[EndIndex] == element)
    return EndIndex;
  return nbrValue;
}

// find an element in a unsorted array
//
// element = element to find
// array = sorted array where to search 
// nbrValue = number of values in array
// return value = element position (negative if not in the array)

template <class ClassName>
int SearchInUnsortedArray(const ClassName& element, ClassName* array, int nbrValue)
{
  for (int i = 0; i < nbrValue; ++i)
    if (array[i] == element)
      return i;
  return -1;
}

// find an element in an array, if found add c to the corresponding weight in weight array, if not found insert the element and set the weight to c
//
// element = element to find
// array = array where to search 
// weight = weight array
// nbrValue = number of values in array
// c = weight to add
// return value = number of inserted value

template <class ClassName>
int SearchInArrayAndSetWeight(ClassName element, ClassName*& array, long*& weigth, unsigned long nbrValue, long c)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      weigth[0] = c;
      return 1;
    }

  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {
       
      MidIndex = (StartIndex + EndIndex) >> 1;
       
      if(array[MidIndex] == element)
        {
	  weigth[MidIndex] += c;
	  return 0;
        }
       
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
   
  if(array[StartIndex] == element)
    {
      weigth[StartIndex] += c;
      return 0;
    }
   
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  long TmpWeigth;
  long TmpWeigth1;
   
  if(array[StartIndex]<element)
    StartIndex++;
   
  TmpElement = array[StartIndex];
  TmpWeigth = weigth[StartIndex];
  array[StartIndex] = element;
  weigth[StartIndex] = c;
  for (unsigned int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      TmpWeigth1 = weigth[i];
      array[i] = TmpElement;
      weigth[i] = TmpWeigth;
      TmpElement = TmpElement1;
      TmpWeigth = TmpWeigth1;
    }
  return 1;
}


// find an element in an array, if found add c to the corresponding weight in weight array, if not found insert the element and set the weight to c
//
// element = element to find
// array = array where to search 
// weight = weight array
// nbrValue = number of values in array
// c = weight to add
// return value = number of inserted value

template <class ClassName>
int SearchInArrayAndDefinedWeight(ClassName element, ClassName*& array, long*& weigth, unsigned long nbrValue, long c)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      weigth[0] = c;
      return 1;
    }

  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {
       
      MidIndex = (StartIndex + EndIndex) >> 1;
       
      if(array[MidIndex] == element)
        {
	  return 0;
        }
       
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
   
  if(array[StartIndex] == element)
    {
      return 0;
    }
   
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  long TmpWeigth;
  long TmpWeigth1;
   
  if(array[StartIndex]<element)
    StartIndex++;
   
  TmpElement = array[StartIndex];
  TmpWeigth = weigth[StartIndex];
  array[StartIndex] = element;
  weigth[StartIndex] = c;
  for (unsigned int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      TmpWeigth1 = weigth[i];
      array[i] = TmpElement;
      weigth[i] = TmpWeigth;
      TmpElement = TmpElement1;
      TmpWeigth = TmpWeigth1;
    }
  return 1;
}

// find an element in an array, if found add c to the corresponding weight in weight array, if not found insert the element and set the weight to c
//
// element = element to find
// array = array where to search 
// weight = weight array
// nbrValue = number of values in array
// c = weight to add
// return value = number of inserted value

template <class ClassName>
int SearchInArrayAndSetWeight(ClassName element, ClassName*& array, double*& weigth, unsigned long nbrValue, double c)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      weigth[0] = c;
      return 1;
    }

  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {
       
      MidIndex = (StartIndex + EndIndex) >> 1;
       
      if(array[MidIndex] == element)
        {
	  weigth[MidIndex] += c;
	  return 0;
        }
       
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
   
  if(array[StartIndex] == element)
    {
      weigth[StartIndex] += c;
      return 0;
    }
   
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  double TmpWeigth;
  double TmpWeigth1;
   
  if(array[StartIndex]<element)
    StartIndex++;
   
  TmpElement = array[StartIndex];
  TmpWeigth = weigth[StartIndex];
  array[StartIndex] = element;
  weigth[StartIndex] = c;
  for (unsigned int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      TmpWeigth1 = weigth[i];
      array[i] = TmpElement;
      weigth[i] = TmpWeigth;
      TmpElement = TmpElement1;
      TmpWeigth = TmpWeigth1;
    }
  return 1;
}

// find an element in an array, if found add c to the corresponding weight in weight array, if not found insert the element and set the weight to c
//
// element = element to find
// array = array where to search 
// weight = weight array
// nbrValue = number of values in array
// c = weight to add
// return value = number of inserted value

template <class ClassName>
int SearchInArrayAndSetWeight(ClassName element, ClassName*& array, LongRational*& weigth, unsigned long nbrValue, LongRational c)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      weigth[0] = c;
      return 1;
    }
  
  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
  
  while((EndIndex - StartIndex) > 1)
    {
      MidIndex = (StartIndex + EndIndex) >> 1;
      
      if(array[MidIndex] == element)
	{
	  weigth[MidIndex] += c;
	  return 0;
	}
      
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
  
  if(array[StartIndex] == element)
    {
      weigth[StartIndex] += c;
      return 0;
    }
  
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  LongRational TmpWeigth;
  LongRational TmpWeigth1;
  
  if(array[StartIndex]<element)
    StartIndex++;
  
  TmpElement = array[StartIndex];
  TmpWeigth = weigth[StartIndex];
  array[StartIndex] = element;
  weigth[StartIndex] = c;
  for (unsigned int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      TmpWeigth1 = weigth[i];
      array[i] = TmpElement;
      weigth[i] = TmpWeigth;
      TmpElement = TmpElement1;
      TmpWeigth = TmpWeigth1;
    }
  return 1;
}

// find an element in an array, if found add c to the corresponding weight in weight array, if not found insert the element and set the weight to c
//
// element = element to find
// array = array where to search 
// weight = weight array
// nbrValue = number of values in array
// c = weight to add
// return value = number of inserted value

template <class ClassName>
int SearchInArrayAndSetWeight(ClassName element, ClassName*& array, Complex *& weigth, unsigned long nbrValue, Complex c)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      weigth[0] = c;
      return 1;
    }

  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
   
  while((EndIndex - StartIndex) > 1)
    {
       
      MidIndex = (StartIndex + EndIndex) >> 1;
       
      if(array[MidIndex] == element)
        {
	  weigth[MidIndex] += c;
	  return 0;
        }
       
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
   
  if(array[StartIndex] == element)
    {
      weigth[StartIndex] += c;
      return 0;
    }
   
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  Complex TmpWeigth;
  Complex TmpWeigth1;
   
  if(array[StartIndex]<element)
    StartIndex++;
   
  TmpElement = array[StartIndex];
  TmpWeigth = weigth[StartIndex];
  array[StartIndex] = element;
  weigth[StartIndex] = c;
  for (unsigned int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      TmpWeigth1 = weigth[i];
      array[i] = TmpElement;
      weigth[i] = TmpWeigth;
      TmpElement = TmpElement1;
      TmpWeigth = TmpWeigth1;
    }
  return 1;
}


// find an element in an array, if not found insert the element such that the array is still sorted
//
// element = element to find
// array = array where to search 
// nbrValue = number of values in array
// return value = number of inserted value

template <class ClassName>
int SearchInSortedArrayAndInsert(ClassName element, ClassName*& array, unsigned long nbrValue)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      return 1;
    }
  
  unsigned long StartIndex = 0l;
  unsigned long EndIndex = nbrValue;
  unsigned long MidIndex;
  
  while((EndIndex - StartIndex) > 1)
    {
      MidIndex = (StartIndex + EndIndex) >> 1;
      
      if(array[MidIndex] == element)
	{
	  return 0;
	}
      
      if(array[MidIndex] > element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
  
  if(array[StartIndex] == element)
    {
      return 0;
    }
  
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  
  if(array[StartIndex] < element)
    StartIndex++;
  
  TmpElement = array[StartIndex];
  array[StartIndex] = element;

  for (unsigned long i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      array[i] = TmpElement;
      TmpElement = TmpElement1;
    }
  return 1;
}


// find an element in an array, if not found insert the element such that the array is still sorted with down ordering
//
// element = element to find
// array = array where to search 
// nbrValue = number of values in array
// return value = negative if no value was inserted, otherwise return the position where the value was inserted

template <class ClassName>
int SearchInDownSortedArrayAndInsert(ClassName element, ClassName*& array, unsigned long nbrValue)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      return 1;
    }
  
  unsigned long StartIndex = 0l;
  unsigned long EndIndex = nbrValue;
  unsigned long MidIndex;
  
  while((EndIndex - StartIndex) > 1)
    {
      MidIndex = (StartIndex + EndIndex) >> 1;
      
      if(array[MidIndex] == element)
	{
	  return -1;
	}
      
      if (array[MidIndex] < element)
	EndIndex = MidIndex;
      else
	StartIndex = MidIndex;
    }
  
  if ((array[StartIndex] == element) || (array[EndIndex] == element))
    {
      return -1;
    }
  
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  
  if (array[StartIndex] > element)
    StartIndex++;
  
  TmpElement = array[StartIndex];
  array[StartIndex] = element;
  MidIndex = StartIndex;

  for (unsigned long i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      array[i] = TmpElement;
      TmpElement = TmpElement1;
    }
  return MidIndex;
}

// find an element in an array, if not found insert the element such that the array is still sorted with down ordering
//
// element = element to find
// array = array where to search 
// nbrValue = number of values in array
// return value = negative if no value was inserted, otherwise return the position where the value was inserted

template <class ClassName>
int SearchInDownSortedArrayAndInsert(ClassName element, ClassName*& array, int nbrValue)
{
  if (nbrValue == 0)
    {
      array[0] = element;
      return 1;
    }
  
  int StartIndex = 0;
  int EndIndex = nbrValue;
  int MidIndex;
  
  while((EndIndex - StartIndex) > 1)
    {
      MidIndex = (StartIndex + EndIndex) >> 1;
      
      if(array[MidIndex] == element)
	{
	  return -1;
	}
      
      MidIndex = (StartIndex + EndIndex) >> 1;       
      if(array[MidIndex] > element)
	StartIndex = MidIndex;
      else
	EndIndex = MidIndex;
    }
  
  if(array[StartIndex] == element)
    {
      return -1;
    }
  
  ++nbrValue;
  ClassName TmpElement;
  ClassName TmpElement1;
  
  if (array[StartIndex] > element)
    StartIndex++;
  
  TmpElement = array[StartIndex];
  array[StartIndex] = element;
  MidIndex = StartIndex;

  for (int i = StartIndex + 1; i < nbrValue; ++i)
    {
      TmpElement1 = array[i];
      array[i] = TmpElement;
      TmpElement = TmpElement1;
    }
  return MidIndex;
}


// compute the product of occupation number factorial of an array
//
// array = array where to search 
// nbrValue = number of values in array
// return value = product of occupation number factorial

template <class ClassName>
unsigned long MultiplicitiesFactorial(ClassName*& array, unsigned long nbrValue)
{
  ClassName Temp;
  int p;
  int Pos = 0;
  FactorialCoefficient Coefficient;
  Coefficient.SetToOne();
	
  while(Pos < nbrValue)
    {
      Temp = array[Pos];
      Pos++;
      p = 1;
      while((Pos < nbrValue)&&(array[Pos] == Temp))
	{
	  Pos++;
	  p++;
	}
      Coefficient.FactorialMultiply(p);
    }
  return Coefficient.GetIntegerValue();
}

// down ordering array sort using quick sort, each element being made of two smaller elements
//
// array1 = pointer to the array of  the element first part 
// array2 = pointer to the array of  the element second part 
// nbrValue = nbr of value in the array

template <class ClassName>
void SortDoubleElementArrayDownOrdering(ClassName* array1, ClassName* array2, long nbrValue)
{
  switch (nbrValue)
    {
    case 0l:
      return;
    case 1l:
      return;
    case 2l:
      {
	if (DoubleElementLessThan(array1[0l], array2[0l], array1[1l], array2[1l]) == true)
	  {
	    ClassName TmpElement = array1[0l];
	    array1[0l] = array1[1l];
	    array1[1l] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[1l];
	    array2[1l] = TmpElement;
	  }
	return;
      }
      break;
    case 3l:
      {
	ClassName TmpElement;
	if (DoubleElementLessThan(array1[0l], array2[0l], array1[1l], array2[1l]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[1l];
	    array1[1l] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[1l];
	    array2[1l] = TmpElement;
	  }
	if (DoubleElementLessThan(array1[1l], array2[1l], array1[2l], array2[2l]) == true)
	  {
	    TmpElement = array1[1l];
	    array1[1l] = array1[2l];
	    array1[2l] = TmpElement;
	    TmpElement = array2[1l];
	    array2[1l] = array2[2l];
	    array2[2l] = TmpElement;
	  }	
	if (DoubleElementLessThan(array1[0l], array2[0l], array1[1l], array2[1l]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[1l];
	    array1[1l] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[1l];
	    array2[1l] = TmpElement;
	  }
	return;
      }
      break;
    default:
      {
	long j = nbrValue - 1l;
	long i = nbrValue >> 1;
	ClassName TmpElement;
	if (DoubleElementLessThan(array1[0l], array2[0l], array1[i], array2[i]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[i];
	    array2[i] = TmpElement;
	  }
	if (DoubleElementLessThan(array1[i], array2[i], array1[j], array2[j]) == true)
	  {
	    TmpElement = array1[i];
	    array1[i] = array1[j];
	    array1[j] = TmpElement;
	    TmpElement = array2[i];
	    array2[i] = array2[j];
	    array2[j] = TmpElement;
	  }
	if (DoubleElementLessThan(array1[0l], array2[0l], array1[i], array2[i]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[i];
	    array2[i] = TmpElement;
	  }
	--j;
	ClassName Pivot1 = array1[i];
	array1[i] = array1[j];
	array1[j] = Pivot1;
	ClassName Pivot2 = array2[i];
	array2[i] = array2[j];
	array2[j] = Pivot2;
	i = 0l;
	while (true)
	  {
	    ++i; 
 	    --j;
	    while (DoubleElementLessThan(Pivot1, Pivot2, array1[i], array2[i]) == true)
	      ++i;
	    while (DoubleElementLessThan(array1[j], array2[j], Pivot1, Pivot2) == true)
	      --j;
	    if (i < j)
	      {
		TmpElement = array1[i];
		array1[i] = array1[j];
		array1[j] = TmpElement;	    
		TmpElement = array2[i];
		array2[i] = array2[j];
		array2[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array1[nbrValue - 2l] = array1[i];
	array1[i] = Pivot1;
	array2[nbrValue - 2l] = array2[i];
	array2[i] = Pivot2;
	SortDoubleElementArrayDownOrdering(array1, array2, i);
	SortDoubleElementArrayDownOrdering(&(array1[i + 1]), &(array2[i + 1]), nbrValue - i - 1l);	
      }
    }
  return;
}

// down ordering array sort using quick sort, each element being made of three smaller elements
//
// array1 = pointer to the array of  the element first part 
// array2 = pointer to the array of  the element second part 
// array3 = pointer to the array of  the element third part 
// nbrValue = nbr of value in the array

template <class ClassName>
void SortTripleElementArrayDownOrdering(ClassName* array1, ClassName* array2, ClassName* array3, long nbrValue)
{
  switch (nbrValue)
    {
    case 0l:
      return;
    case 1l:
      return;
    case 2l:
      {
	if (TripleElementLessThan(array1[0], array2[0], array3[0], 
				  array1[1], array2[1], array3[1]) == true)
	  {
	    ClassName TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	  }
	return;
      }
      break;
    case 3l:
      {
	ClassName TmpElement;
	if (TripleElementLessThan(array1[0], array2[0], array3[0], 
				  array1[1], array2[1], array3[1]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	  }
	if (TripleElementLessThan(array1[1], array2[1], array3[1], 
				  array1[2], array2[2], array3[2]) == true)
	  {
	    TmpElement = array1[1];
	    array1[1] = array1[2];
	    array1[2] = TmpElement;
	    TmpElement = array2[1];
	    array2[1] = array2[2];
	    array2[2] = TmpElement;
	    TmpElement = array3[1];
	    array3[1] = array3[2];
	    array3[2] = TmpElement;
	  }	
	if (TripleElementLessThan(array1[0], array2[0], array3[0], 
				  array1[1], array2[1], array3[1]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	  }
	return;
      }
      break;
    default:
      {
	long j = nbrValue - 1l;
	long i = nbrValue >> 1;
	ClassName TmpElement;
	if (TripleElementLessThan(array1[0l], array2[0l], array3[0l], 
				  array1[i], array2[i], array3[i]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[i];
	    array2[i] = TmpElement;
	    TmpElement = array3[0l];
	    array3[0l] = array3[i];
	    array3[i] = TmpElement;
	  }
	if (TripleElementLessThan(array1[i], array2[i], array3[i], 
				  array1[j], array2[j], array3[j]) == true)
	  {
	    TmpElement = array1[i];
	    array1[i] = array1[j];
	    array1[j] = TmpElement;
	    TmpElement = array2[i];
	    array2[i] = array2[j];
	    array2[j] = TmpElement;
	    TmpElement = array3[i];
	    array3[i] = array3[j];
	    array3[j] = TmpElement;
	  }
	if (TripleElementLessThan(array1[0l], array2[0l], array3[0l], 
				  array1[i], array2[i], array3[i]) == true)
	  {
	    TmpElement = array1[0l];
	    array1[0l] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0l];
	    array2[0l] = array2[i];
	    array2[i] = TmpElement;
	    TmpElement = array3[0l];
	    array3[0l] = array3[i];
	    array3[i] = TmpElement;
	  }
	--j;
	ClassName Pivot1 = array1[i];
	array1[i] = array1[j];
	array1[j] = Pivot1;
	ClassName Pivot2 = array2[i];
	array2[i] = array2[j];
	array2[j] = Pivot2;
	ClassName Pivot3 = array3[i];
	array3[i] = array3[j];
	array3[j] = Pivot3;
	i = 0l;
	while (true)
	  {
	    ++i; 
 	    --j;
	    while (TripleElementLessThan(Pivot1, Pivot2, Pivot3, 
					 array1[i], array2[i], array3[i]) == true)
	      ++i;
	    while (TripleElementLessThan(array1[j], array2[j], array3[j], 
					 Pivot1, Pivot2, Pivot3) == true)
	      --j;
	    if (i < j)
	      {
		TmpElement = array1[i];
		array1[i] = array1[j];
		array1[j] = TmpElement;	    
		TmpElement = array2[i];
		array2[i] = array2[j];
		array2[j] = TmpElement;	    
		TmpElement = array3[i];
		array3[i] = array3[j];
		array3[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array1[nbrValue - 2l] = array1[i];
	array1[i] = Pivot1;
	array2[nbrValue - 2l] = array2[i];
	array2[i] = Pivot2;
	array3[nbrValue - 2l] = array3[i];
	array3[i] = Pivot3;
	SortTripleElementArrayDownOrdering(array1, array2, array3, i);
	SortTripleElementArrayDownOrdering(&(array1[i + 1]), &(array2[i + 1]), &(array3[i + 1]), nbrValue - i - 1l);	
      }
    }
  return;
}

// down ordering array sort using quick sort, each element being made of four smaller elements
//
// array1 = pointer to the array of  the element first part 
// array2 = pointer to the array of  the element second part 
// array3 = pointer to the array of  the element third  part 
// array4 = pointer to the array of  the element fourth part 
// nbrValue = nbr of value in the array

template <class ClassName>
void SortQuadElementArrayDownOrdering(ClassName* array1, ClassName* array2, ClassName* array3, ClassName* array4, long nbrValue)
{
  switch (nbrValue)
    {
    case 0:
      return;
    case 1:
      return;
    case 2:
      {
	if (QuadElementLessThan(array1[0], array2[0], array3[0], array4[0], 
				array1[1], array2[1], array3[1], array4[1]) == true)
	  {
	    ClassName TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	    TmpElement = array4[0];
	    array4[0] = array4[1];
	    array4[1] = TmpElement;
	  }
	return;
      }
      break;
    case 3:
      {
	ClassName TmpElement;
	if (QuadElementLessThan(array1[0], array2[0], array3[0], array4[0], 
				array1[1], array2[1], array3[1], array4[1]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	    TmpElement = array4[0];
	    array4[0] = array4[1];
	    array4[1] = TmpElement;
	  }
	if (QuadElementLessThan(array1[1], array2[1], array3[1], array4[1], 
				array1[2], array2[2], array3[2], array4[2]) == true)
	  {
	    TmpElement = array1[1];
	    array1[1] = array1[2];
	    array1[2] = TmpElement;
	    TmpElement = array2[1];
	    array2[1] = array2[2];
	    array2[2] = TmpElement;
	    TmpElement = array3[1];
	    array3[1] = array3[2];
	    array3[2] = TmpElement;
	    TmpElement = array4[1];
	    array4[1] = array4[2];
	    array4[2] = TmpElement;
	  }	
	if (QuadElementLessThan(array1[0], array2[0], array3[0], array4[0], 
				array1[1], array2[1], array3[1], array4[1]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[1];
	    array1[1] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[1];
	    array2[1] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[1];
	    array3[1] = TmpElement;
	    TmpElement = array4[0];
	    array4[0] = array4[1];
	    array4[1] = TmpElement;
	  }
	return;
      }
      break;
    default:
      {
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	if (QuadElementLessThan(array1[0], array2[0], array3[0], array4[0], 
				array1[i], array2[i], array3[i], array4[i]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[i];
	    array2[i] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[i];
	    array3[i] = TmpElement;
	    TmpElement = array4[0];
	    array4[0] = array4[i];
	    array4[i] = TmpElement;
	  }
	if (QuadElementLessThan(array1[i], array2[i], array3[i], array4[i], 
				array1[j], array2[j], array3[j], array4[j]) == true)
	  {
	    TmpElement = array1[i];
	    array1[i] = array1[j];
	    array1[j] = TmpElement;
	    TmpElement = array2[i];
	    array2[i] = array2[j];
	    array2[j] = TmpElement;
	    TmpElement = array3[i];
	    array3[i] = array3[j];
	    array3[j] = TmpElement;
	    TmpElement = array4[i];
	    array4[i] = array4[j];
	    array4[j] = TmpElement;
	  }
	if (QuadElementLessThan(array1[0], array2[0], array3[0], array4[0], 
				array1[i], array2[i], array3[i], array4[i]) == true)
	  {
	    TmpElement = array1[0];
	    array1[0] = array1[i];
	    array1[i] = TmpElement;
	    TmpElement = array2[0];
	    array2[0] = array2[i];
	    array2[i] = TmpElement;
	    TmpElement = array3[0];
	    array3[0] = array3[i];
	    array3[i] = TmpElement;
	    TmpElement = array4[0];
	    array4[0] = array4[i];
	    array4[i] = TmpElement;
	  }
	--j;
	ClassName Pivot1 = array1[i];
	array1[i] = array1[j];
	array1[j] = Pivot1;
	ClassName Pivot2 = array2[i];
	array2[i] = array2[j];
	array2[j] = Pivot2;
	ClassName Pivot3 = array3[i];
	array3[i] = array3[j];
	array3[j] = Pivot3;
	ClassName Pivot4 = array4[i];
	array4[i] = array4[j];
	array4[j] = Pivot4;
	i = 0;
	while (true)
	  {
	    ++i; 
 	    --j;
	    while (QuadElementLessThan(Pivot1, Pivot2, Pivot3, Pivot4,
				       array1[i], array2[i], array3[i], array4[i]) == true)
	      ++i;
	    while (QuadElementLessThan(array1[j], array2[j], array3[j], array4[j], 
				       Pivot1, Pivot2, Pivot3, Pivot4) == true)
	      --j;
	    if (i < j)
	      {
		TmpElement = array1[i];
		array1[i] = array1[j];
		array1[j] = TmpElement;	    
		TmpElement = array2[i];
		array2[i] = array2[j];
		array2[j] = TmpElement;	    
		TmpElement = array3[i];
		array3[i] = array3[j];
		array3[j] = TmpElement;	    
		TmpElement = array4[i];
		array4[i] = array4[j];
		array4[j] = TmpElement;	    
	      }
	    else
	      break;
	  }	
	array1[nbrValue - 2] = array1[i];
	array1[i] = Pivot1;
	array2[nbrValue - 2] = array2[i];
	array2[i] = Pivot2;
	array3[nbrValue - 2] = array3[i];
	array3[i] = Pivot3;
	array4[nbrValue - 2] = array4[i];
	array4[i] = Pivot4;
	SortQuadElementArrayDownOrdering(array1, array2, array3, array4, i);
	SortQuadElementArrayDownOrdering(&(array1[i + 1]), &(array2[i + 1]), &(array3[i + 1]), &(array4[i + 1]), nbrValue - i - 1);	
      }
    }
  return;
}


// Sort an array from the smallest element to the largest, removing all the duplicate elements and resizing the array if needed
//
// array = reference to the array pointer
// nbrValues = reference on the number of elements (it will be modified to indicate the number of distinct elements

template <class ClassName>
void SortArrayUpOrderingAndRemoveDuplicates(ClassName*& array,  long& nbrValues)
{
  SortArrayUpOrdering(array, nbrValues);
  long NbrDistinctValues = 0;
  long TmpIndex = 0l;
  while (TmpIndex < nbrValues)
    {
      long TmpIndex2 = TmpIndex + 1l;
      while ((TmpIndex2 < nbrValues) && (array[TmpIndex] == array[TmpIndex2]))
	{
	  ++TmpIndex2;
	}
      ++NbrDistinctValues;
      TmpIndex = TmpIndex2;
    }
  ClassName* TmpArray = new ClassName[NbrDistinctValues];
  NbrDistinctValues = 0;
  TmpIndex = 0l;
  while (TmpIndex < nbrValues)
    {
      long TmpIndex2 = TmpIndex + 1l;
      while ((TmpIndex2 < nbrValues) && (array[TmpIndex] == array[TmpIndex2]))
	{
	  ++TmpIndex2;
	}
      TmpArray[NbrDistinctValues] = array[TmpIndex];
      ++NbrDistinctValues;
      TmpIndex = TmpIndex2;
    }

  delete[] array;  
  array = TmpArray;
  nbrValues = NbrDistinctValues;
  return;
}

#endif
