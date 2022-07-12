#include "GroupedPermutations.h"

#include <iostream>
using std::cout;
using std::endl;


// default constructor
// nbrGroups = number of groups
// elementsPerGroup = number of elements per group
// orderedGroups = flag indicating whether order of groups matters
GroupedPermutations::GroupedPermutations(int nbrGroups, unsigned elementsPerGroup, bool orderedGroups)
{
  this->NbrGroups=nbrGroups;
  this->ElementsPerGroup=elementsPerGroup;
  this->NbrElements = nbrGroups*elementsPerGroup;
  this->OrderedGroups = orderedGroups;
  this->NbrBitsForElements = getHighestBit(NbrElements);
  this->NbrBitsPerGroup = getHighestBit(NbrGroups);  
  this->PermutationList = OrderedList<PermutationElement>(/* eliminateDuplicates */ true);
  this->MyArray = new unsigned[NbrElements];
  this->MapOfGroups = new unsigned[NbrGroups];
  this->InverseMapOfGroups = new unsigned[NbrGroups];
  this->CountOfGroups = new unsigned[NbrGroups];
  for (int i=0; i<NbrGroups; ++i)
    {
      this->MapOfGroups[i]=i;
      this->InverseMapOfGroups[i]=i;
    }

  this->CentralRecursion(this->GetInitialString(),SmallIntegerArray(this->NbrGroups), 1);

  this->NbrPermutations = PermutationList.GetNbrElement();
  this->Permutations = new SmallIntegerArray[NbrPermutations];
  this->Multiplicities = new unsigned long[NbrPermutations];

  ListElement<PermutationElement> *ElementPointer = PermutationList.FirstElement;
  int i = 0;
  while (i<NbrPermutations)
    {
      this->Permutations[i] = ElementPointer->Element.Value;
      this->Multiplicities[i] = ElementPointer->Element.Multiplicity;      
      ElementPointer = ElementPointer->NextPointer;
      ++i;
    }
  this->PermutationList.DeleteList();
  
}


// destructor
GroupedPermutations::~GroupedPermutations()
{
  delete [] MyArray;
  delete [] MapOfGroups;
  delete [] CountOfGroups;
}


// central recursive function that generates all different permutations
void GroupedPermutations::CentralRecursion(SmallIntegerArray remainingElements, SmallIntegerArray permutation, unsigned long multiplicity)
{
//   for (int i=0; i<permutation.GetNbrElements(); ++i) cout << " ";
//   cout << "CentralRecursion ("<<remainingElements<<", "<<permutation<<", "<<multiplicity<<")"<<endl;
  int NbrRemainingElements = remainingElements.GetNbrElements();
  if (NbrRemainingElements>1)
    {      
      for (int i=0; i<NbrRemainingElements; ++i)
	{
	  for (int j=0; j<i; ++j)
	    MyArray[j]=remainingElements.GetElement(j);
	  for (int j=i+1; j<NbrRemainingElements; ++j)
	    MyArray[j-1]=remainingElements.GetElement(j);	  
	  CentralRecursion(SmallIntegerArray(NbrRemainingElements-1, this->NbrGroups, MyArray),
			   SmallIntegerArray(permutation, remainingElements.GetElement(i)), 1);
	}
    }
  else // NbrRemainingElements==1
    {
      SmallIntegerArray FinalPermutation(permutation, remainingElements.GetElement(0));
      PermutationElement PE(this->GetPermutationString(FinalPermutation),1);
      int Pos;
      PermutationElement *Duplicate;
      //cout << "Inserting " << PE.Value << endl;
      this->PermutationList.Insert(PE, Pos, Duplicate);
      //cout << "Found duplicate at pos "<<Pos<<": " << Duplicate << endl;
      if (Duplicate!=NULL)
	{
	  Duplicate->Multiplicity+=1;
	}
//       cout << "Values in List are now:"<<endl;
//       for (int i=0; i<PermutationList.GetNbrElement(); ++i)
// 	{
// 	  cout << PermutationList[i];
// 	  if (PermutationList[i] < PermutationList[PermutationList.GetNbrElement()-1])
// 	    cout << " < " <<PermutationList[PermutationList.GetNbrElement()-1] << endl;
// 	  else
// 	    cout << " >= " <<PermutationList[PermutationList.GetNbrElement()-1] << endl;
// 	}
//       cout << "End List"<<endl;
    }
  
}

// translate internal form of permutations to a canonic expression
SmallIntegerArray GroupedPermutations::GetPermutationString(SmallIntegerArray &permutation)
{
  for (int i=0; i<NbrGroups; ++i)
    CountOfGroups[i]=0;
//  cout << "translating "<<permutation;
  if (!this->OrderedGroups)
    {
//       cout << "getting map for "<<permutation<<endl;
      int PresentGroup=0;
      int PresentElement=1;
      unsigned PresentElementValue;
      this->MapOfGroups[PresentGroup++]=permutation.GetElement(0);
      while (PresentGroup < NbrGroups)
	{
	  PresentElementValue=permutation.GetElement(PresentElement++);
	  bool GroupKnown=false;
	  for (int i=0; (i<PresentGroup)&&(!GroupKnown); ++i)
	    if (PresentElementValue==this->MapOfGroups[i])
	      GroupKnown=true;
	  if (!GroupKnown)
	    {
	      this->MapOfGroups[PresentGroup++]=PresentElementValue;
	    }
	}
      for (int i=0; i<NbrGroups; ++i)
	{
	  this->InverseMapOfGroups[MapOfGroups[i]]=i;
// 	  cout << "Map["<<i<<"]="<<MapOfGroups[i]<<endl;
	}
    }
  for (int i=0; i<NbrElements; ++i)
    {            
      unsigned PresentElementValue=permutation.GetElement(i);
      unsigned PresentGroup=MapOfGroups[PresentElementValue];
      MyArray[PresentGroup*ElementsPerGroup+CountOfGroups[PresentGroup]]=i;
      ++CountOfGroups[PresentGroup];
    }
  // SmallIntegerArray tmp(this->NbrElements, this->NbrElements, MyArray);
  // cout << " -> result: "<< tmp <<endl;
  return SmallIntegerArray(this->NbrElements, this->NbrElements, MyArray);     
}



// get an initial string without permutations
SmallIntegerArray GroupedPermutations::GetInitialString()
{  
  int count=0;
  for (int g=0; g<NbrGroups; ++g)
    for (int e=0; e<ElementsPerGroup; ++e)
      MyArray[count++]=g;
  return SmallIntegerArray(this->NbrElements, this->NbrGroups, MyArray);
}
