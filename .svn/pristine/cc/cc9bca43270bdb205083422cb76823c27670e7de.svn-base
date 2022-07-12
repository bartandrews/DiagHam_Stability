

////////////////////////////////////////////////////////////////////////////////



#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include <algorithm>
#include <set>

#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/OrderedList.h"
#include "GeneralTools/List.h"
#include "GeneralTools/SmallIntegerArray.h"
#include "GeneralTools/ArrayTools.h"

// evaluate all the different way to split the permutations of a given in two sub groups (not compulsorily the same number of elements in the two subgroups)
//
// nbrPermutations = number of different permutations
// nbrElement = number of element in the permutations
// nbrElementPerColor = number of element in the first group
// permutations1 = array where will be stored the permutations of the first group
// permutations2 = array where will be stored the permutations of the second group

inline void EvaluatePermutationsOfSubGroups(unsigned long nbrPermutations, int nbrElement, int nbrElementPerColor, unsigned long * permutations1, unsigned long * permutations2)
{
  unsigned long NbrElementOtherColor = nbrElement - nbrElementPerColor;
  unsigned long MinValue = (0x1ul << nbrElementPerColor) - 0x1ul;
  unsigned long MaxValue = MinValue << (NbrElementOtherColor);
  unsigned long* TmpArrayPerm = new unsigned long [nbrElement];
  nbrPermutations = 0;
  for (; MinValue <= MaxValue; ++MinValue)
    {
      int Count = 0;
      int Pos = 0;
      while ((Pos < nbrElement) && (Count <= nbrElementPerColor))
	{
	  if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	    ++Count;
	  ++Pos;
	}
      if (Count == nbrElementPerColor)
	{
	  int Pos1 = 0;
	  int Pos2 = nbrElementPerColor;
	  for (Pos = 0; Pos < nbrElement; ++Pos)
	    {
	      if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
		{
		  TmpArrayPerm[Pos1] = (unsigned long) Pos;
		  ++Pos1;
		}
	      else
		{
		  TmpArrayPerm[Pos2] = (unsigned long) Pos;
		  ++Pos2;
		}
	    }
	  unsigned long TmpPerm2 = 0ul;
	  unsigned long TmpPerm3 = 0ul;
	  for (int i = 0; i < nbrElementPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i * 5);
	    }
	  for (unsigned i = 0; i < NbrElementOtherColor; ++i)
	    {
	      TmpPerm3 |= TmpArrayPerm[i + nbrElementPerColor]<< (i*5);
	    }
	  permutations1[nbrPermutations] = TmpPerm2;
	  permutations2[nbrPermutations] = TmpPerm3;	      
	  ++nbrPermutations;
	}
    }
  delete[] TmpArrayPerm;
  return;
}

// evaluate all the different way to split the permutations of a given in two sub groups but exclude ones which are symmetric with higher index (not compulsorily the same number of elements in the two subgroups)
//
// nbrPermutations = number of different permutations
// nbrElement = number of element in the permutations
// nbrElementPerColor = number of element in the first group
// permutations1 = array where will be stored the permutations of the first group
// permutations2 = array where will be stored the permutations of the second group

inline void EvaluatePermutationsOfSubGroupsSymmetric(unsigned long nbrPermutations, int nbrElement, int nbrElementPerColor, unsigned long * permutations1, unsigned long * permutations2)
{
  unsigned long MinValue = (0x1ul << nbrElementPerColor) - 0x1ul;
  unsigned long MaxValue = MinValue << (nbrElement - nbrElementPerColor);
  unsigned long* TmpArrayPerm = new unsigned long [nbrElement];
  unsigned long Mask = (0x1ul << nbrElement) - 1;
  nbrPermutations = 0;
  for (; MinValue <= MaxValue; ++MinValue)
    {
      int Count = 0;
      int Pos = 0;
      while ((Pos < nbrElement) && (Count <= nbrElementPerColor))
	{
	  if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
	    ++Count;
	  ++Pos;
	}
      if (Count == nbrElementPerColor && (MinValue < ((~MinValue) & Mask  )) )
	{
	  int Pos1 = 0;
	  int Pos2 = nbrElementPerColor;
	  for (Pos = 0; Pos < nbrElement; ++Pos)
	    {
	      if (((MinValue >> Pos) & 0x1ul) != 0x0ul)
		{
		  TmpArrayPerm[Pos1] = (unsigned long) Pos;
		  ++Pos1;
		}
	      else
		{
		  TmpArrayPerm[Pos2] = (unsigned long) Pos;
		  ++Pos2;
		}
	    }
	  unsigned long TmpPerm2 = 0ul;
	  unsigned long TmpPerm3 = 0ul;
	  for (int i = 0; i < nbrElementPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i * 5);
	      TmpPerm3 |= TmpArrayPerm[i + nbrElementPerColor] << (i *5);
	    }
	  permutations1[nbrPermutations] = TmpPerm2;
	  permutations2[nbrPermutations] = TmpPerm3;	      
	  ++nbrPermutations;
	}
    }
  delete[] TmpArrayPerm;
  return;
}

// evaluate all the different permutations of a given monomial and keep track of the signs
//
// nbrParticles = number of different permutations
// monomial = array storing the monomial
// nbrPermutations = reference to the number of permutations
// monomialPermutatiosn = reference to array where permutations will be stored
// permutationSigns = reference to array where permutation signs will be stored

inline void EvaluateMonomialPermutations(int nbrParticles, unsigned long* monomial, int& nbrPermutations, unsigned long**& monomialPermutations, double*& permutationSigns)
{
  FactorialCoefficient *FactCoef = new FactorialCoefficient();
  FactCoef->FactorialMultiply(nbrParticles);
  monomialPermutations = new unsigned long*[FactCoef->GetIntegerValue()];  
  permutationSigns = new double[FactCoef->GetIntegerValue()];
  delete FactCoef;
  unsigned long* TmpPermutation = new unsigned long[nbrParticles];
  unsigned long* TmpPermutation2 = new unsigned long[nbrParticles];
  unsigned long* TmpPermutation3 = new unsigned long[nbrParticles];
  memcpy(TmpPermutation, monomial, sizeof(unsigned long) * nbrParticles);
  int NumPermutations;
  
  nbrPermutations = 0;
  do 
    {
      monomialPermutations[nbrPermutations] = new unsigned long[nbrParticles];      
      memcpy(monomialPermutations[nbrPermutations], TmpPermutation, sizeof(unsigned long) * nbrParticles);
      
      memcpy(TmpPermutation2, monomial, sizeof(unsigned long) * nbrParticles);
      memcpy(TmpPermutation3, monomialPermutations[nbrPermutations], sizeof(unsigned long) * nbrParticles);
      NumPermutations = 0;
      SortArrayDownOrdering(TmpPermutation3,  TmpPermutation2, nbrParticles, NumPermutations);
      if ( (NumPermutations & 0x1) == 1 )
	{
	  permutationSigns[nbrPermutations] = -1.0;
	}
      else
	{
	  permutationSigns[nbrPermutations] = 1.0;
	}
      nbrPermutations++;
    }
  while (std::prev_permutation(TmpPermutation, TmpPermutation + nbrParticles));
  delete [] TmpPermutation;
  delete [] TmpPermutation2;
  delete [] TmpPermutation3;
}

//Generate all possible partitions of a set into two, one of size sizeOfSubset.
//Return value is the number of possible partitions, i.e. the size of SubSetFermionOccupationResults

inline int AllGivenSizeSubsets(unsigned int * setToPartitionOccupationBasis, const int nbrElements, const int maxElement, int sizeOfSubset, OrderedList< SmallIntegerArray > & SubsetOccupationBasisResults, List< SmallIntegerArray > & SubsetOccupationBasisComplements) {
  //initialise SubsetOccupationBasisComplements to setToPartitionOccupationBasis
  cout << "AllGivenSizeSubsets: Initialising storage\n";

  int NbrElements = nbrElements;
  int MaxElement = maxElement;

  unsigned int * setToPartitionAndComplement = new unsigned int [2*MaxElement];
  //SmallIntegerArray NewPartitionAndComplement(2*MaxElement, NbrElements); 

  for(int i=0; i<MaxElement; i++) {
    setToPartitionAndComplement[i] = 0;
  }
  for(int i=0; i<MaxElement; i++) {
    setToPartitionAndComplement[MaxElement+i] = setToPartitionOccupationBasis[i];
  }

  SmallIntegerArray * ZerothPartitionAndComplement = new SmallIntegerArray(2*MaxElement, NbrElements, setToPartitionAndComplement);

  delete [] setToPartitionAndComplement;

  //OrderedList< SmallIntegerArray > * CurrentSizeSubsetOccupations = new OrderedList< SmallIntegerArray >(true); //eliminateDuplicates = true
  //*CurrentSizeSubsetOccupations += (*ZerothPartitionAndComplement);
  //OrderedList< SmallIntegerArray > * NextSizeSubsetOccupations = new OrderedList< SmallIntegerArray>(true); //eliminateDuplicates = true
  std::set< SmallIntegerArray > * CurrentSizeSubsetOccupations = new std::set< SmallIntegerArray >;
  CurrentSizeSubsetOccupations->insert(*ZerothPartitionAndComplement);

  cout << "AllGivenSizeSubsets: Beginning iteration\n";
  if( sizeOfSubset > 0 && sizeOfSubset < NbrElements )
    {
      
      for( int tempSubsetSize = 0; tempSubsetSize < sizeOfSubset; tempSubsetSize++ )
	{
	  
	  //clear the holders for the next occupations
	  //NextSizeSubsetOccupations->DeleteList();
	  
	  // OrderedList< SmallIntegerArray > * NextSizeSubsetOccupations = new OrderedList< SmallIntegerArray >(true); //eliminateDuplicates = true
	  std::set< SmallIntegerArray > * NextSizeSubsetOccupations = new std::set< SmallIntegerArray >; 
	  int pass = 0;
	  for( std::set< SmallIntegerArray >::iterator tempSubsetIndex = CurrentSizeSubsetOccupations->begin(); tempSubsetIndex != CurrentSizeSubsetOccupations->end(); tempSubsetIndex++ )
	    {
	      //cout << "Calculating partitions from smaller partition " << pass << " of " << CurrentSizeSubsetOccupations->size() << "\n\t";
	      int complementElement;
	      //int insertionPosition;
	      for( int subsetElement = 0; subsetElement < MaxElement; subsetElement++) 
		{
		  complementElement = subsetElement + MaxElement;
		  SmallIntegerArray NewPartitionAndComplement = *tempSubsetIndex;
		  if( NewPartitionAndComplement.GetElement( complementElement ) > 0 ) {
		    //cout << "Element to add " << subsetElement << "\n";
		    //SmallIntegerArray NewPartitionAndComplement = *tempSubsetIndex;
		    //cout << "Copied\n";
		    NewPartitionAndComplement.SetElement( complementElement, NewPartitionAndComplement.GetElement(complementElement) - 1 );
		    //cout << "complement element decremented\n";
		    NewPartitionAndComplement.SetElement( subsetElement, NewPartitionAndComplement.GetElement(subsetElement) + 1 );
		    //cout << "subset element incremented\n";
		    NextSizeSubsetOccupations->insert(NewPartitionAndComplement);
		    
		  }
		}
	      pass++;
	    }
	  *CurrentSizeSubsetOccupations = *NextSizeSubsetOccupations;
	  delete NextSizeSubsetOccupations;
	  //cout << "nextsize set to currentsize\n";
	}
      // delete NextSizeSubsetOccupations;
    
      //Output results into SubsetOccupationBasisResults and SubsetOccupationBasisComplements. Might want to consider using two lists for the results
      //cout << "Saving results to SubsetOccupationBasisResults and SubsetOccupationBasisComplements\n";
      int insertionPosition;
      unsigned int * subsetHolder = new unsigned int [MaxElement];
      unsigned int * complementHolder = new unsigned int [MaxElement];
      SmallIntegerArray subsetcomplement;
      int iteration = 0;
      for(std::set< SmallIntegerArray >::iterator i = CurrentSizeSubsetOccupations->begin(); i != CurrentSizeSubsetOccupations->end(); i++) {
	//cout << "Saving partition " << iteration << " of " << CurrentSizeSubsetOccupations->size() << "\n";
	subsetcomplement = *i;
	for(int j=0; j<MaxElement; j++) {
	  //cout << "\tSaving element " << j << " of " << MaxElement << " value " << subsetcomplement.GetElement(j) << "\n";
	  subsetHolder[j] = subsetcomplement.GetElement(j);
	  //cout << "\tsubset element saved\n";
	  complementHolder[j] = subsetcomplement.GetElement(MaxElement+j);
	  //cout << "\tcomplement element saved\n";
	    }
	SmallIntegerArray subset(MaxElement, sizeOfSubset, subsetHolder);
	SmallIntegerArray complement(MaxElement, NbrElements - sizeOfSubset, complementHolder);
	SmallIntegerArray * duplicate;
	SubsetOccupationBasisResults.Insert(subset, insertionPosition, duplicate);
	//cout << "Partition subset " << iteration << " saved\n";
	SubsetOccupationBasisComplements.Insert(complement, insertionPosition);
	//cout << "Partition complement " << iteration << " saved\n";
	iteration++;
      }
	  
      delete ZerothPartitionAndComplement;
      delete CurrentSizeSubsetOccupations;
      delete [] subsetHolder;
      delete [] complementHolder;
      
      return SubsetOccupationBasisResults.GetNbrElement();
    }
  else 
    {
      if( sizeOfSubset == 0 ) 
	{
	  int insertionPosition;
	  unsigned int * subsetHolder = new unsigned int [MaxElement];
	  unsigned int * complementHolder = new unsigned int [MaxElement];
	  SmallIntegerArray subsetcomplement;
	  int iteration = 0;
	  for(std::set< SmallIntegerArray >::iterator i = CurrentSizeSubsetOccupations->begin(); i != CurrentSizeSubsetOccupations->end(); i++) {
	    //cout << "Saving partition " << iteration << " of " << CurrentSizeSubsetOccupations->size() << "\n";
	    subsetcomplement = *i;
	    for(int j=0; j<MaxElement; j++) {
	      //cout << "\tSaving element " << j << " of " << MaxElement << " value " << subsetcomplement.GetElement(j) << "\n";
	      subsetHolder[j] = subsetcomplement.GetElement(j);
	      //cout << "\tsubset element saved\n";
	      complementHolder[j] = subsetcomplement.GetElement(MaxElement+j);
	      //cout << "\tcomplement element saved\n";
	    }
	    SmallIntegerArray subset(MaxElement, 1, subsetHolder); //set highest integer to store to one to avoid SmallIntegerArray throwing an exception
	    SmallIntegerArray complement(MaxElement, NbrElements - sizeOfSubset, complementHolder);
	    SmallIntegerArray * duplicate;
	    SubsetOccupationBasisResults.Insert(subset, insertionPosition, duplicate);
	    //cout << "Partition subset " << iteration << " saved\n";
	    SubsetOccupationBasisComplements.Insert(complement, insertionPosition);
	    //cout << "Partition complement " << iteration << " saved\n";
	    iteration++;
	  }
	  delete ZerothPartitionAndComplement;
	  delete CurrentSizeSubsetOccupations;
	  delete [] subsetHolder;
	  delete [] complementHolder;
	  
	  return 1; //one trivial subset found

	}
      else if( sizeOfSubset == NbrElements )
	{
	  int insertionPosition;
	  unsigned int * subsetHolder = new unsigned int [MaxElement];
	  unsigned int * complementHolder = new unsigned int [MaxElement];
	  SmallIntegerArray subsetcomplement;
	  int iteration = 0;
	  for(std::set< SmallIntegerArray >::iterator i = CurrentSizeSubsetOccupations->begin(); i != CurrentSizeSubsetOccupations->end(); i++) {
	    //cout << "Saving partition " << iteration << " of " << CurrentSizeSubsetOccupations->size() << "\n";
	    subsetcomplement = *i;
	    for(int j=0; j<MaxElement; j++) {
	      //cout << "\tSaving element " << j << " of " << MaxElement << " value " << subsetcomplement.GetElement(j) << "\n";
	      subsetHolder[j] = subsetcomplement.GetElement(MaxElement+j); //swap partition and complement
	      //cout << "\tsubset element saved\n";
	      complementHolder[j] = subsetcomplement.GetElement(j); //swap partition and complement
	      //cout << "\tcomplement element saved\n";
	    }
	    SmallIntegerArray subset(MaxElement, sizeOfSubset, subsetHolder);
	    SmallIntegerArray complement(MaxElement, 1, complementHolder); //highestInteger set to one to avoid SmallIntegerArray throwing an exception
	    SmallIntegerArray * duplicate;
	    SubsetOccupationBasisResults.Insert(subset, insertionPosition, duplicate);
	    //cout << "Partition subset " << iteration << " saved\n";
	    SubsetOccupationBasisComplements.Insert(complement, insertionPosition);
	    //cout << "Partition complement " << iteration << " saved\n";
	    iteration++;
	  }
	  delete ZerothPartitionAndComplement;
	  delete CurrentSizeSubsetOccupations;
	  delete [] subsetHolder;
	  delete [] complementHolder;

	  return 1; //one trivial subset found
	}
      else
	{
	  delete ZerothPartitionAndComplement;
	  delete CurrentSizeSubsetOccupations;
	  cout << "sizeOfSubset should be in range {0, ..., NbrElements  = " << NbrElements << "}\n";
	  return -1;
	}
    }
}

#endif
