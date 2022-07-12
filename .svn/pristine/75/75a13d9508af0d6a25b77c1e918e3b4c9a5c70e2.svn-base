////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-20016 Nicolas Regnault                 //
//                         class author: Antoine Sterdyniak                   //
//                                                                            //
//                                                                            //
//                 class of local spins around a single PEPS tensor           //
//                                                                            //
//                        last modification : 09/05/2016                      //
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


#include "PEPSLocalPhysicalAndVirtualSpin.h"
#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;



PEPSLocalPhysicalAndVirtualSpin::PEPSLocalPhysicalAndVirtualSpin (int numberSpin, int physicalSpinValue, int nbrVirtualSpinRepresentation, int * tableOfVirtualSpinRepresentations, int sz)
{
  this->ChainLength =  numberSpin;
  this->Sz = sz;
  this->PhysicalSpinValue  = physicalSpinValue;
  this->NbrVirtualSpinRepresentation =  nbrVirtualSpinRepresentation;
  this->TableOfVirtualSpinRepresentations = tableOfVirtualSpinRepresentations;
  this->NbrSpinStatesForVirtualSpins = 0;
  for(int i =0; i <  this->NbrVirtualSpinRepresentation;i++)
    {
      this->NbrSpinStatesForVirtualSpins += (this->TableOfVirtualSpinRepresentations[i]+1);
    }
  this->TableOfSzValuesForVirtualSpin = new int [this->NbrSpinStatesForVirtualSpins];
  this->TableOfSzValuesForPhysicalSpin = new int [this->PhysicalSpinValue+1];
  this->TableOfSValuesForVirtualSpin = new int [this->NbrSpinStatesForVirtualSpins];

  int Index = 0;
  for(int i =0; i < this->PhysicalSpinValue+1 ;i++)  
    { 
       this->TableOfSzValuesForPhysicalSpin[i] = 2 * i -  this->PhysicalSpinValue;
    }
  for(int i =0; i < this->NbrVirtualSpinRepresentation;i++)
    {
      for(int p=0; p <this->TableOfVirtualSpinRepresentations[i]+1;p++)
	{
	  this->TableOfSzValuesForVirtualSpin[Index] = 2*p - this->TableOfVirtualSpinRepresentations[i];
	  this->TableOfSValuesForVirtualSpin[Index] =  this->TableOfVirtualSpinRepresentations[i];
	  Index++;
	}
    }

  PowerOfNbrSpinStatesForVirtualSpins = new int[this->ChainLength];
  PowerOfNbrSpinStatesForVirtualSpins[0] = 1;
  for(int i = 1 ; i < this->ChainLength ; i++)
    {
      PowerOfNbrSpinStatesForVirtualSpins[i] = PowerOfNbrSpinStatesForVirtualSpins[i-1] * NbrSpinStatesForVirtualSpins;
    }
  
  
  int TmpHilbertSpaceDimension = 0;
  for(int i =0 ; i <  this->PhysicalSpinValue+1; i++)
    {
      TmpHilbertSpaceDimension+=this->EvaluateHilbertSpaceDimension(this->ChainLength -1,this->Sz- (2*i - this->PhysicalSpinValue));
    }
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->StateDescription = new unsigned long[TmpHilbertSpaceDimension];
  
  cout<<"HilbertSpaceDimension = "  <<this->HilbertSpaceDimension<<endl;
  TmpHilbertSpaceDimension = 0;
  for(int i =0 ; i <   this->PhysicalSpinValue+1; i++)
    {
      TmpHilbertSpaceDimension=this->GeneratesStates(this->ChainLength -1,this->Sz- (2*i - this->PhysicalSpinValue), i * PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength - 1], TmpHilbertSpaceDimension);
    }  
  if (TmpHilbertSpaceDimension !=   this->HilbertSpaceDimension ) 
    {
      cout <<"Error while generating HilbertSpace "<<TmpHilbertSpaceDimension << " "<<this->HilbertSpaceDimension<<endl;
    }
}

AbstractHilbertSpace* PEPSLocalPhysicalAndVirtualSpin::Clone()
{
  return new PEPSLocalPhysicalAndVirtualSpin(this->ChainLength ,this->PhysicalSpinValue, this->NbrVirtualSpinRepresentation, this->TableOfVirtualSpinRepresentations,this->Sz);
}

PEPSLocalPhysicalAndVirtualSpin::~PEPSLocalPhysicalAndVirtualSpin()
{
  delete [] this->StateDescription;
  delete [] this->TableOfVirtualSpinRepresentations;
  delete [] this->TableOfSzValuesForVirtualSpin;
  delete [] this->TableOfSValuesForVirtualSpin;
  delete [] this->PowerOfNbrSpinStatesForVirtualSpins;
  delete [] this->TableOfSzValuesForPhysicalSpin;
}



int PEPSLocalPhysicalAndVirtualSpin::EvaluateHilbertSpaceDimension(int nbrVirtualSpinToBeAttributed, int szToBeRealized)
{
//  cout <<"nbrVirtualSpinToBeAttributed = "<< nbrVirtualSpinToBeAttributed<<" " <<szToBeRealized<<endl;
  int Count = 0; 
  if (nbrVirtualSpinToBeAttributed == 1 )
    {
      for(int i =0; i <this->NbrSpinStatesForVirtualSpins;i++)
	{
	  if (this->TableOfSzValuesForVirtualSpin[i] ==  szToBeRealized )
	    Count++;
	}
      return Count;
    }
  
  for(int i =0; i <this->NbrSpinStatesForVirtualSpins;i++)
    {
      Count+=this->EvaluateHilbertSpaceDimension(nbrVirtualSpinToBeAttributed-1, szToBeRealized-this->TableOfSzValuesForVirtualSpin[i]);
    } 
  return Count;
}


int PEPSLocalPhysicalAndVirtualSpin::GeneratesStates(int nbrVirtualSpinToBeAttributed, int szToBeRealized, int beginningOfStateRepresentation, int pos)
{
  if (nbrVirtualSpinToBeAttributed == 1 )
    {
      for(int i =0; i <this->NbrSpinStatesForVirtualSpins;i++)
	{
	  if (this->TableOfSzValuesForVirtualSpin[i] ==  szToBeRealized )
	    {
	      this->StateDescription[pos] =  beginningOfStateRepresentation + i;
	      pos++;
	    }
	}
      return pos;
    }
  
  for(int i =0; i <this->NbrSpinStatesForVirtualSpins;i++)
    {
      pos = this->GeneratesStates(nbrVirtualSpinToBeAttributed-1, szToBeRealized-this->TableOfSzValuesForVirtualSpin[i],beginningOfStateRepresentation+i * PowerOfNbrSpinStatesForVirtualSpins[nbrVirtualSpinToBeAttributed-1],pos );
    }
  
  return pos;  
}


void PEPSLocalPhysicalAndVirtualSpin::PrintConversionTable()
{
  for(int i =0; i <this->NbrSpinStatesForVirtualSpins;i++)
    {
      cout <<i << " "<<this->TableOfSValuesForVirtualSpin[i]<<" "<<this->TableOfSzValuesForVirtualSpin[i]<<endl;
    }
  cout <<"Physical Table"<<endl;
  for(int i =0; i < this->PhysicalSpinValue+1;i++)
    {
      cout <<i << " "<<this->TableOfSzValuesForPhysicalSpin[i]<<endl;
    }
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& PEPSLocalPhysicalAndVirtualSpin::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  
  unsigned long TmpState = this->StateDescription[state];
  
  int ArrayRepresentation[this->ChainLength];
  this->ConvertFromIntegerToArray(TmpState ,ArrayRepresentation);
  
  for(int i = 0 ; i < this->ChainLength; i++)
    {
      Str << ArrayRepresentation[i]<<" ";
    }
  
  return Str;
}


// find state index
//
// state = state description
// return value = corresponding index

int PEPSLocalPhysicalAndVirtualSpin::FindStateIndex(unsigned long state)
{
  int PosMin = 0;
  int PosMax = this->HilbertSpaceDimension-1;
  int PosMid = ( PosMin + PosMax)>>1;

  while ( ( PosMax - PosMin ) > 1 )
    {
      if (this->StateDescription[PosMid] < state ) 
	{
	  PosMin = PosMid;
	}
      else
	{
	  if (this->StateDescription[PosMid] == state )  
	    return PosMid;
	  PosMax = PosMid;
	}
      PosMid = ( PosMin + PosMax) >> 1;
    }
  if (this->StateDescription[PosMin] ==  state )
    return PosMin;
  if (this->StateDescription[PosMax] ==  state )
    return PosMax;
  return  this->HilbertSpaceDimension;   
}

// get the value of the spin (i.e. S) at a given site
// 
// site = site index
// return value = twice the spin

int PEPSLocalPhysicalAndVirtualSpin::GetLocalSpin(int site, int state)
{
  if (site ==0)
    return this->PhysicalSpinValue;
  int LocalIndexSite =  (this->StateDescription[state]/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength - 1 - site])%this->PowerOfNbrSpinStatesForVirtualSpins[1];
//  cout << site<< " "<<state << " " <<this->TableOfSValuesForVirtualSpin[LocalIndexSite]<<endl;
  return  this->TableOfSValuesForVirtualSpin[LocalIndexSite];
}




// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue
double PEPSLocalPhysicalAndVirtualSpin::SziSzj (int i, int j, int state)
{
  int TmpState = this->StateDescription[state];
  int LocalIndexSiteI =  (TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1 - i])%this->PowerOfNbrSpinStatesForVirtualSpins[1];
  int LocalIndexSiteJ  =  (TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength - 1 - j])%this->PowerOfNbrSpinStatesForVirtualSpins[1]; 
//  cout <<i<< " "<<j<< " "<<state<< " "<<this->TableOfSzValuesForVirtualSpin[LocalIndexSiteI] *  this->TableOfSzValuesForVirtualSpin[LocalIndexSiteJ]*0.25<<endl
  
  if (i==0 ) 
    {
      LocalIndexSiteI =  TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1];
      if (j==0 ) 
	{
	  LocalIndexSiteJ =  TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1];
	  return this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteI] *   this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteJ]*0.25; 
	}
      return  this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteI] *  this->TableOfSzValuesForVirtualSpin[LocalIndexSiteJ]*0.25; 
    }
    
  if (j==0)
    {
      LocalIndexSiteJ =  TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1];
      return  this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteJ] *  this->TableOfSzValuesForVirtualSpin[LocalIndexSiteI]*0.25; 
    }
  
  return  this->TableOfSzValuesForVirtualSpin[LocalIndexSiteI] *  this->TableOfSzValuesForVirtualSpin[LocalIndexSiteJ]*0.25; 
}


// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

int PEPSLocalPhysicalAndVirtualSpin::SmiSpj (int i, int j, int state, double& coefficient)
{
  int TmpState = this->StateDescription[state];
  int LocalIndexSiteI =  (TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1 - i])%this->PowerOfNbrSpinStatesForVirtualSpins[1];
  int LocalIndexSiteJ  =  (TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength -1 - j])%this->PowerOfNbrSpinStatesForVirtualSpins[1];   
  
  int CoefI = 0;
  int CoefJ = 0;
  
  if  (i==0 )
    {
      LocalIndexSiteI =  TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1];
      CoefI =  (this->PhysicalSpinValue * (this->PhysicalSpinValue + 2) -  (  this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteI] * ( this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteI] - 2 ) ));
    }
  else
    CoefI =   (this->TableOfSValuesForVirtualSpin[LocalIndexSiteI] * (this->TableOfSValuesForVirtualSpin[LocalIndexSiteI] + 2) -  ( this->TableOfSzValuesForVirtualSpin[LocalIndexSiteI] * (this->TableOfSzValuesForVirtualSpin[LocalIndexSiteI] - 2 ) ));

 if  (j==0 )
   {
     LocalIndexSiteJ =  TmpState/this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1];
     CoefJ =  (this->PhysicalSpinValue * (this->PhysicalSpinValue + 2) -  (  this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteJ] * ( this->TableOfSzValuesForPhysicalSpin[LocalIndexSiteJ] + 2 ) ));

   }
 else
   CoefJ =   (this->TableOfSValuesForVirtualSpin[LocalIndexSiteJ] * (this->TableOfSValuesForVirtualSpin[LocalIndexSiteJ] + 2) -  ( this->TableOfSzValuesForVirtualSpin[LocalIndexSiteJ] * (this->TableOfSzValuesForVirtualSpin[LocalIndexSiteJ] + 2 ) ));
 
 if ( (CoefI  == 0 ) || (CoefJ  == 0 ))
   return this->HilbertSpaceDimension; 
 coefficient = sqrt((double)  CoefJ * CoefI)*0.25;
 TmpState+=this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength -1 - j] - this->PowerOfNbrSpinStatesForVirtualSpins[this->ChainLength-1 - i];
 return this->FindStateIndex(TmpState);
}



// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

void PEPSLocalPhysicalAndVirtualSpin::ApplyRotation (ComplexVector & initialState, ComplexVector &  FinalState)
{
  int ArrayRepresentation[this->ChainLength];
  
  for( int i = 0 ; i < this->HilbertSpaceDimension; i++)
    {
      this->ConvertFromIntegerToArray(this->StateDescription[i], ArrayRepresentation);

      int Tmp = ArrayRepresentation[1];
      for(int j = 1 ; j < this->ChainLength - 1; j++)
	{
	  ArrayRepresentation[j] = ArrayRepresentation[j+1];
	}
      ArrayRepresentation[this->ChainLength-1] = Tmp;
      
      FinalState[this->FindStateIndex(this->ConvertFromArrayToInteger(ArrayRepresentation))] =  initialState[i];
    }
}


// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

void PEPSLocalPhysicalAndVirtualSpin::ApplyHReflexion (ComplexVector & initialState, ComplexVector &  FinalState)
{
  int ArrayRepresentation[this->ChainLength];
  for(int i = 0 ; i < this->HilbertSpaceDimension; i++)
    {
      this->ConvertFromIntegerToArray(this->StateDescription[i], ArrayRepresentation);
      
      for (int t = 0 ; t < (this->ChainLength-1)/2 ; t++)
	{
	  int Tmp = ArrayRepresentation[t+2];
	  ArrayRepresentation[t+2] = ArrayRepresentation[this->ChainLength - 1 - t];
	  ArrayRepresentation[this->ChainLength - 1 - t] = Tmp;
	}
      
      FinalState[this->FindStateIndex(this->ConvertFromArrayToInteger(ArrayRepresentation))] =  initialState[i];
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state 

void PEPSLocalPhysicalAndVirtualSpin::ApplyDXReflexion (ComplexVector & initialState, ComplexVector &  FinalState)
{
  if ( ( this->ChainLength -1 )% 2 != 0)
    {
      cout << "Diagonal Reflexion is not a symmetry of the tensor"<<endl;
      return ;
    }
  
  int ArrayRepresentation[this->ChainLength];
  
  for( int i = 0 ; i < this->HilbertSpaceDimension; i++)
    {
      this->ConvertFromIntegerToArray(this->StateDescription[i], ArrayRepresentation);
      
      for (int t = 1 ; t <= (this->ChainLength-1)/2 ; t++)
	{
	  int Tmp = ArrayRepresentation[t];
	  ArrayRepresentation[t] = ArrayRepresentation[this->ChainLength - t];
	  ArrayRepresentation[this->ChainLength -  t] = Tmp;
	}

/*      if ( i ==0)
	{
	  for(int p = 0 ; p <this->ChainLength ; p++)
	    cout <<ArrayRepresentation[p]<<" ";
	  cout <<endl;
	}*/
      FinalState[this->FindStateIndex(this->ConvertFromArrayToInteger(ArrayRepresentation))] =  initialState[i];

    }
}
