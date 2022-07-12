////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2008 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//  allows the calculation of permutations of groups and their multiplicity   //
//                                                                            //
//                        last modification : 16/04/2008                      //
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


#include "KCombinations.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/SmallIntegerArray.h"
#include "GeneralTools/OrderedList.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <cstdlib>

using std::cout;
using std::endl;


// default constructor
// nbrChoices = number of elements to choose
// nbrElements = total number of elements
// orderedGroups = flag indicating whether order of groups matters
KCombinations::KCombinations(int nbrChoices, int nbrElements)
{
  this->NbrElements=nbrElements;
  this->NbrChoices=nbrChoices;
  FactorialCoefficient Coeff;
  Coeff.SetToOne();
  Coeff.FactorialMultiply(nbrElements);
  Coeff.FactorialDivide(nbrElements-nbrChoices);
  Coeff.FactorialDivide(nbrChoices);
  this->NbrCombinations=Coeff.GetIntegerValue();
  this->Combinations = new int*[NbrCombinations];
  int TmpNbr = this->GenerateCombinations(nbrChoices-1, 0, 0);
  if (TmpNbr != NbrCombinations)
    {
      cout << "Discrepancy of count in KCombinations"<<endl;
      exit(1);
    }
}

// destructor
KCombinations::~KCombinations()
{
  if (NbrCombinations!=0)
    {
      for (int i=0; i<NbrCombinations; ++i)
	delete [] Combinations[i];
      delete [] Combinations;
    }
}

// central recursive function that generates all different permutations
// currentPos = current position to be chosen (start from zero)
// minVal = minimum value for current position
// pos = place where to start storing things
// return = total number of quantum numbers
int KCombinations::GenerateCombinations(int currentPos, int minVal, int pos)
{
  if (currentPos<NbrChoices-1)
    {
      for (int val=minVal; val<NbrElements-(NbrChoices-1-currentPos); ++val)
	{
	  int TmpPos = pos;
	  pos = GenerateCombinations(currentPos+1, val+1, TmpPos);
	  for (int i=TmpPos; i<pos; ++i)
	    Combinations[i][currentPos]=val;
	}
    }
  else
    {
      for (int val=minVal; val<NbrElements; ++val)
	{
	  Combinations[pos]=new int[NbrChoices];
	  Combinations[pos++][NbrChoices-1]=val;
	}
    }
  return pos;
}



// output stream overload
ostream& operator << (ostream & Str, KCombinations& a)
{
  for (int i=0; i<a.NbrCombinations; ++i)
    {
      Str << i<<": "<<a.Combinations[i][0];
      for (int k=1; k<a.NbrChoices; ++k)
	Str <<" "<<a.Combinations[i][k];
      Str << endl;
    }
  return Str;
}
