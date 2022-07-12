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
	  }
	if (array[1] < array[2])
	  {
	    TmpElement = array[1];
	    array[1] = array[2];
	    array[2] = TmpElement;
	  }	
	if (array[0] < array[1])
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
	if (array[0] <  array[i])
	  {
	    TmpElement = array[i];
	    array[i] = array[0];
	    array[0] = TmpElement;
	  }
	if (array[i] <  array[j])
	  {
	    TmpElement = array[i];
	    array[i] = array[j];
	    array[j] = TmpElement;
	  }
	if (array[0] <  array[i])
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
	array[nbrValue - 2] = array[i];
	array[i] = Pivot;
	SortArrayDownOrdering(array, i);
	SortArrayDownOrdering(&(array[i + 1]), nbrValue - i - 1);	
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
	    int TmpElement2 = integerArray[0];
	    integerArray[0] = integerArray[1];
	    integerArray[1] = TmpElement2;
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
	int j = nbrValue - 1;
	int i = nbrValue >> 1;
	ClassName TmpElement;
	int TmpElement2;
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
	int Pivot2 = integerArray[i];
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

#endif
