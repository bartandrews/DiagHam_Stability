////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                         to particles on a torus with                       //
//                                                                            //
//                        last modification : 27/06/2003                      //
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


#include "config.h"
#include "Hamiltonian/AbstractQHEOnTorusHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

AbstractQHEOnTorusHamiltonian::AbstractQHEOnTorusHamiltonian()
{
  this->FilterInteractionFlag = false;
  this->FilterInteractionIndices = 0;
  this->NbrFilterInteractionIndices = 0;
}

// destructor
//

AbstractQHEOnTorusHamiltonian::~AbstractQHEOnTorusHamiltonian()
{
  if (this->FilterInteractionFlag == true)
    {
      delete[] this->FilterInteractionIndices ;
    }
}

// get all the indices that should appear in the annihilation/creation operators
//

void AbstractQHEOnTorusHamiltonian::GetIndices()
{
  this->NbrSectorSums = this->NbrLzValue;
  this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
  for (int i = 0; i < this->NbrSectorSums; ++i)
    this->NbrSectorIndicesPerSum[i] = 0;      
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      if (this->FilterInteractionFlag == true)
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      {
		if (SearchInArray<int>(m1 * this->NbrLzValue + m2, this->FilterInteractionIndices, this->NbrFilterInteractionIndices) >= 0)
		  ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	      }
	}
      else
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	}
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      if (this->FilterInteractionFlag == true)
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      {
		if (SearchInArray<int>(m1 * this->NbrLzValue + m2, this->FilterInteractionIndices, this->NbrFilterInteractionIndices) >= 0)
		  {
		    int TmpSum = (m1 + m2) % this->NbrLzValue;
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
	}
      else
	{
	  for (int m1 = 0; m1 < this->LzMax; ++m1)
	    for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	}
    }
  else
    {
      if (this->FilterInteractionFlag == true)
	{
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      {
		if (SearchInArray<int>(m1 * this->NbrLzValue + m2, this->FilterInteractionIndices, this->NbrFilterInteractionIndices) >= 0)
		  ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	      }
	}
      else
	{
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      ++this->NbrSectorIndicesPerSum[(m1 + m2) % this->NbrLzValue];
	}
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      if (this->FilterInteractionFlag == true)
	{
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      {
		if (SearchInArray<int>(m1 * this->NbrLzValue + m2, this->FilterInteractionIndices, this->NbrFilterInteractionIndices) >= 0)
		  {
		    int TmpSum = (m1 + m2) % this->NbrLzValue;
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
	}
      else
	{
	  for (int m1 = 0; m1 <= this->LzMax; ++m1)
	    for (int m2 = m1; m2 <= this->LzMax; ++m2)
	      {
		int TmpSum = (m1 + m2) % this->NbrLzValue;
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = m1;
		this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = m2;
		++this->NbrSectorIndicesPerSum[TmpSum];    
	      }
	}
    }
}

// find all the indices that have to be kept when filtering the hamiltionian
//
// filterInteractionFile = name of the file that describe which terms in the interaction should be kept

void AbstractQHEOnTorusHamiltonian::FindFilteredIndices(char* filterInteractionFile)
{
  MultiColumnASCIIFile FilterRules;
  if (FilterRules.Parse(filterInteractionFile) == false)
    {
      FilterRules.DumpErrors(cout);
      return;
    }
  int TmpNbrNBody = FilterRules.GetNbrColumns();
  if (TmpNbrNBody <= 1)
    {
      cout << "error, " << filterInteractionFile << " should have at least two columns" << endl;
    }
  return;
}
