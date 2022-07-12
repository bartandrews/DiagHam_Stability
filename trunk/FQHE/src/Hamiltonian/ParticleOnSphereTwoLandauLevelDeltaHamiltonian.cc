////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          Copyright (C) 2001-2004 Niall Moran and Nicolas Regnault          //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                   two Landau levels and delta interaction                  //
//                                                                            //
//                        last modification : 14/03/2011                      //
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


#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "MathTools/ThreeJSymbol.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>



using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state in the lower Landau level
// landauLevelIndexDifference = difference of indices between the lower and higher Landau levels
// pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
//                   first index refered to the spin sector (sorted as down-dowm, up-up, up-down, down-up)
// cyclotronEnergy = cyclotron energy in e^2/(epsilon l_b) unit
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereTwoLandauLevelDeltaHamiltonian::ParticleOnSphereTwoLandauLevelDeltaHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax,
											       double** pseudoPotential, double* cyclotronEnergy,
											       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName, bool showIntFactorsFlag)
{
  this->Particles = particles;
  this->LzMax = lzmax; //maximum Lz which is on the SLL.
  this->NbrLzValue = this->LzMax + 1; //Number of values in Lz range (0 to LzMaz).
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->OneBodyTermFlag = false;
  this->Architecture = architecture;
  this->TotalCyclotronEnergy = cyclotronEnergy;
  this->LzMaxDown = this->LzMax - 1; // Largest Lz value on LLL
  this->LzMaxUp = this->LzMax; // Largest Lz value on SLL  
  this->ShowIntFactorsFlag = showIntFactorsFlag;
  this->LzFermionDownShift = 0;
  this->LzFermionUpShift = 0;
  
  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      this->LzFermionDownShift = 1;      
    }
  
  if ( pseudoPotential != NULL )
    {
      this->UsePseudoPotentials = true;
      // {1,1,1,1}, {1,1,0,0}, {1,1,1,0}, {0,0,1,1}, {0,0,0,0}, {0,0,1,0}, {1,0,1,1}, {1,0,0,0}, {1,0,1,0}, {1,0,0,1}
      this->PseudoPotentials = new double* [10];
      this->PseudoPotentialMins = new int [10];
      this->NbrPseudoPotentialCoeffs = new int[10];
      int TwoQ = this->LzMax - 2;
      this->NbrPseudoPotentialCoeffs[0] = TwoQ + 3; this->PseudoPotentialMins[0] = 0;
      this->NbrPseudoPotentialCoeffs[1] = TwoQ + 1; this->PseudoPotentialMins[1] = 0;
      this->NbrPseudoPotentialCoeffs[2] = TwoQ + 1; this->PseudoPotentialMins[2] = 1;
      this->NbrPseudoPotentialCoeffs[3] = TwoQ + 1; this->PseudoPotentialMins[3] = 0; 
      this->NbrPseudoPotentialCoeffs[4] = TwoQ + 1; this->PseudoPotentialMins[4] = 0;
      this->NbrPseudoPotentialCoeffs[5] = TwoQ; this->PseudoPotentialMins[5] = 1;
      this->NbrPseudoPotentialCoeffs[6] = TwoQ + 1; this->PseudoPotentialMins[6] = 1;
      this->NbrPseudoPotentialCoeffs[7] = TwoQ; this->PseudoPotentialMins[7] = 1;
      this->NbrPseudoPotentialCoeffs[8] = TwoQ + 1; this->PseudoPotentialMins[8] = 1;
      this->NbrPseudoPotentialCoeffs[9] = TwoQ + 1; this->PseudoPotentialMins[9] = 1;
      
      for ( int i = 0 ; i < 10; i++ ) 
	{
	  this->PseudoPotentials[i] = new double[this->NbrPseudoPotentialCoeffs[i]]; 
	  cout << "Pseudo type "<<i<<endl;
	  for ( int j = 0; j < this->NbrPseudoPotentialCoeffs[i] ; j++ )
	    {
	      this->PseudoPotentials[i][j] = pseudoPotential[i][this->NbrPseudoPotentialCoeffs[i] - 1 - j];
	      cout << this->PseudoPotentials[i][j] << " ";
	    }
          cout<<endl;	  
	}            
    }
  else
    this->UsePseudoPotentials = false;
      
  // Calculation interaction factors.
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex); 
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;    
  this->L2Hamiltonian = 0;
  this->S2Hamiltonian = 0;  
      
    
  if (precalculationFileName == 0)
    {
      if (memory > 0)
        {
          cout << "done" << endl;
          long TmpMemory = this->FastMultiplicationMemory(memory);
          cout << "done" << endl;
          if (TmpMemory < 1024)
            cout  << "fast = " <<  TmpMemory << "b ";
          else
            if (TmpMemory < (1 << 20))
              cout  << "fast = " << (TmpMemory >> 10) << "kb ";
            else
	      if (TmpMemory < (1 << 30))
		cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
          else
            {
              cout  << "fast = " << (TmpMemory >> 30) << ".";
              TmpMemory -= ((TmpMemory >> 30) << 30);
              TmpMemory *= 100l;
              TmpMemory >>= 30;
              if (TmpMemory < 10l)
                cout << "0";
              cout  << TmpMemory << " Gb ";
            }
          if (this->DiskStorageFlag == false)
            {
              this->EnableFastMultiplication();
            }
          else
            {
              char* TmpFileName = this->Architecture->GetTemporaryFileName();
              this->EnableFastMultiplicationWithDiskStorage(TmpFileName);             
              delete[] TmpFileName;
            }
        }
    }
  else
    this->LoadPrecalculation(precalculationFileName);
 
}

// destructor
//

ParticleOnSphereTwoLandauLevelDeltaHamiltonian::~ParticleOnSphereTwoLandauLevelDeltaHamiltonian() 
{
  if ( this->UsePseudoPotentials )
    {
      for (int j = 0; j < 10; ++j)
        delete[] this->PseudoPotentials[j];
      delete[] this->PseudoPotentials;
      delete[] this->NbrPseudoPotentialCoeffs;
      delete[] this->PseudoPotentialMins;
    }
}

// evaluate all interaction factors
//   

void ParticleOnSphereTwoLandauLevelDeltaHamiltonian::EvaluateInteractionFactors()
{
  //set these to 0 so old method not used.
  this->NbrIntraSectorSums = 0;
  this->NbrInterSectorSums = 0;
  this->M1IntraValue = 0;
  this->M2IntraValue = 0;
  this->M1InterValue = 0;
  this->M2InterValue = 0;
    
  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {
      // Set the number of possible sums for each sector. 
      this->NbrUpUpSectorSums = 2 * this->LzMaxUp + 1; // Ranges from 0 to 2*LzMax for given sector.
      this->NbrDownDownSectorSums = 2 * this->LzMaxDown + 1 - 2; // the -2 is because a sum of 0 or 1 is not possible in the LLL. 
      this->NbrUpDownSectorSums = this->LzMaxUp + this->LzMaxDown; // goes from 1 to LzMaxUp + LzMaxDown
      
      //Allocate space for sums and set counts to zero.
      this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
      this->NbrUpDownSectorIndicesPerSum = new int [this->NbrUpDownSectorSums];
      this->NbrDownDownSectorIndicesPerSum = new int [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	this->NbrUpUpSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	this->NbrUpDownSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	this->NbrDownDownSectorIndicesPerSum[i] = 0;

      // Count number of combinations that sum to each.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  this->NbrUpUpSectorIndicesPerSum[m1 + m2]++;   
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrUpDownSectorIndicesPerSum[m1 + m2 - 1]++;
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2]++;

      // Allocate sapce for indices and reset counters to 0 so can be used as indices.
      this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
	{
	  this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
	  this->NbrUpUpSectorIndicesPerSum[i] = 0;
	}
      this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
      if (this->NbrUpDownSectorIndicesPerSum[i] > 0)
      {
	  this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
	  this->NbrUpDownSectorIndicesPerSum[i] = 0;
      }
      this->DownDownSectorIndicesPerSum = new int* [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
	{
	  this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
	  this->NbrDownDownSectorIndicesPerSum[i] = 0;
	}

      // set the indices.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  {
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1] = m1;
	    this->UpUpSectorIndicesPerSum[(m1 + m2)][1 + (this->NbrUpUpSectorIndicesPerSum[(m1 + m2)] << 1)] = m2;
	    ++this->NbrUpUpSectorIndicesPerSum[(m1 + m2)];
	  }
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1] = m1;
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][1 + (this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1)] = m2;
	    ++this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)];
	  }
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->DownDownSectorIndicesPerSum[m1 + m2 - 2][this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2] << 1] = m1;
	    this->DownDownSectorIndicesPerSum[m1 + m2 - 2][1 + (this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2] << 1)] = m2;
	    ++this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 2];
	  }	      
    } 
  else if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      // Set the number of possible sums for each sector. 
      this->NbrUpUpSectorSums     = 2 * this->LzMaxUp - 1; // Ranges from 1 to 2*LzMaxUp-1 for given sector.
      this->NbrDownDownSectorSums = 2 * this->LzMaxDown - 1 - 2; // (Ranges from 3 to LzMaxDown + LzMaxDown - 1)  the -2 is because a sum of 0 or 1 is not possible in the LLL. 
      this->NbrUpDownSectorSums   = this->LzMaxUp + this->LzMaxDown ; // goes from 1 to LzMaxUp + LzMaxDown
      
      //Allocate space for sums and set counts to zero.
      this->NbrUpUpSectorIndicesPerSum = new int[this->NbrUpUpSectorSums];
      this->NbrUpDownSectorIndicesPerSum = new int [this->NbrUpDownSectorSums];
      this->NbrDownDownSectorIndicesPerSum = new int [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	this->NbrUpUpSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	this->NbrUpDownSectorIndicesPerSum[i] = 0;
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	this->NbrDownDownSectorIndicesPerSum[i] = 0;

      // Count number of combinations that sum to each.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  if ( m1 != m2 ) this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1]++;   
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  this->NbrUpDownSectorIndicesPerSum[m1 + m2 - 1]++;
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  if ( m1 != m2 ) this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3]++;

/*      cout << "Indices UpUp "<<endl;	
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	cout << i << " " << this->NbrUpUpSectorIndicesPerSum[i] << endl;

      cout << "Indices UpDown "<<endl;	
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
	cout << i << " " << this->NbrUpDownSectorIndicesPerSum[i] << endl;

      cout << "Indices DownDown "<<endl;	
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	cout << i << " " << this->NbrDownDownSectorIndicesPerSum[i] << endl;
*/
      

      // Allocate sapce for indices and reset counters to 0 so can be used as indices.
      this->UpUpSectorIndicesPerSum = new int* [this->NbrUpUpSectorSums];
      for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
	  {
	    this->UpUpSectorIndicesPerSum[i] = new int[2 * this->NbrUpUpSectorIndicesPerSum[i]];      
	    this->NbrUpUpSectorIndicesPerSum[i] = 0;
	  } 
	else 
	  {
	    this->UpUpSectorIndicesPerSum[i] = 0;
	  }
      this->UpDownSectorIndicesPerSum = new int* [this->NbrUpDownSectorSums];
      for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
        if (this->NbrUpDownSectorIndicesPerSum[i] > 0)
          {
	    this->UpDownSectorIndicesPerSum[i] = new int[2 * this->NbrUpDownSectorIndicesPerSum[i]];      
	    this->NbrUpDownSectorIndicesPerSum[i] = 0;
          }
	else 
	  {
	    this->UpDownSectorIndicesPerSum[i] = 0;
	  }
      this->DownDownSectorIndicesPerSum = new int* [this->NbrDownDownSectorSums];
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
	if (this->NbrDownDownSectorIndicesPerSum[i] > 0)
	  {
	    this->DownDownSectorIndicesPerSum[i] = new int[2 * this->NbrDownDownSectorIndicesPerSum[i]];      
	    this->NbrDownDownSectorIndicesPerSum[i] = 0;
	  }
	else 
	  {
	    this->DownDownSectorIndicesPerSum[i] = 0;
	  }  
	

      // set the indices.
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 0; m2 <= this->LzMaxUp; ++m2)
	  {
	    if ( m1 != m2 ) 
	      {
		this->UpUpSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1] << 1] = m1;
		this->UpUpSectorIndicesPerSum[m1 + m2 - 1][(this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1] << 1) + 1] = m2;
		++this->NbrUpUpSectorIndicesPerSum[m1 + m2 - 1];
	      }
	  }
	  
      for (int m1 = 0; m1 <= this->LzMaxUp; ++m1)
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1] = m1;
	    this->UpDownSectorIndicesPerSum[m1 + m2 - 1][(this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)] << 1) + 1] = m2 - this->LzFermionDownShift;
	    ++this->NbrUpDownSectorIndicesPerSum[(m1 + m2 - 1)];
	  }
	  
      for (int m1 = 1; m1 <= this->LzMaxDown; ++m1)
	{
	for (int m2 = 1; m2 <= this->LzMaxDown; ++m2)
	  {
	    if ( m1 != m2 ) 
	      {
		this->DownDownSectorIndicesPerSum[m1 + m2 - 3][this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3] << 1] = m1 - this->LzFermionDownShift;
		this->DownDownSectorIndicesPerSum[m1 + m2 - 3][(this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3] << 1) + 1] = m2 - this->LzFermionDownShift;
		++this->NbrDownDownSectorIndicesPerSum[m1 + m2 - 3];
	      }
	  }
	}	
	
    }
    
  // Create interaction factor arrays and initialise to 0.
  
  // the three that end in UpUp
  this->InteractionFactorsUpUpUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsUpDownUpUp = new double* [this->NbrUpUpSectorSums];
  this->InteractionFactorsDownDownUpUp = new double* [this->NbrUpUpSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpUp[i] = 0;
      this->InteractionFactorsUpDownUpUp[i] = 0;
      this->InteractionFactorsDownDownUpUp[i] = 0;
    }
  
  // finally the three that end in Up Down
  this->InteractionFactorsUpUpUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsUpDownUpDown = new double* [this->NbrUpDownSectorSums];
  this->InteractionFactorsDownDownUpDown = new double* [this->NbrUpDownSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrUpDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpUpDown[i] = 0;
      this->InteractionFactorsUpDownUpDown[i] = 0;
      this->InteractionFactorsDownDownUpDown[i] = 0;
    }
  
  //now the trhee that end in Down Down
  this->InteractionFactorsDownDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpDownDownDown = new double* [this->NbrDownDownSectorSums];
  this->InteractionFactorsUpUpDownDown = new double* [this->NbrDownDownSectorSums];
  //this is prob unnecessary but no harm.
  for (int i = 0; i < this->NbrDownDownSectorSums; ++i)
    {
      this->InteractionFactorsUpUpDownDown[i] = 0;
      this->InteractionFactorsUpDownDownDown[i] = 0;
      this->InteractionFactorsDownDownDownDown[i] = 0;
    }
        
  

  long TotalNbrInteractionFactors = 0;
  if ( this->UsePseudoPotentials )  // if pseudopotentials were specified these are used to set the interaction coefficients. 
    {

      ClebschGordanCoefficients ClebschDownDown (this->LzMaxDown - 1, this->LzMaxDown - 1);
      ClebschGordanCoefficients ClebschUpUp (this->LzMaxUp, this->LzMaxUp);
      ClebschGordanCoefficients ClebschUpDown (this->LzMaxUp, this->LzMaxDown - 1);
      ClebschGordanCoefficients ClebschDownUp (this->LzMaxDown - 1, this->LzMaxUp);

      int J;
      
      // int m4;
      double ClebschCoef;

      double FermionicSign = -1.0;
      

      int Sign = 0; 
      if (this->LzMax & 1)
	Sign = 1;
      double TmpCoefficient = 0.0;

            
      //Arrays of pseudopotentials from 0-9 where each corresponds to following intereaction.
      //"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
      //"PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
      //"PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown", "PseudopotentialsUpDownDownUp";
	
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
	{
	  //upup-upup term  
	  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	    {
	      if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
		{
		  this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]];
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[0] ; J < this->NbrPseudoPotentialCoeffs[0] + this->PseudoPotentialMins[0] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.CarefulGetCoefficient(m1,m2,J*2);			     	
				  //cout << "ClebschCeof: (" << m1 << ", " << m2 << ", " << (J * 1) << ") = " <<  ClebschUpUp.CarefulGetCoefficient(m1,m2,J*2) << endl;
				  TmpCoefficient = ClebschCoef * ClebschUpUp.CarefulGetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpUpUp[i][Index] += this->PseudoPotentials[0][J-this->PseudoPotentialMins[0]] * TmpCoefficient;
				}
			    }
			  TotalNbrInteractionFactors++;
			  ++Index;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Down Down" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpUpDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
			  this->InteractionFactorsUpUpDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[1] ; J < this->NbrPseudoPotentialCoeffs[1] + this->PseudoPotentialMins[1] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpUp.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpDownDown[i][Index] += this->PseudoPotentials[1][J-this->PseudoPotentialMins[1]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsUpUpDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/ 2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0 ,1.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Up Down" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpUpUpDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i+1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
			  this->InteractionFactorsUpUpUpDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[2] ; J < this->NbrPseudoPotentialCoeffs[2] + this->PseudoPotentialMins[2] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpUp.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpUpDown[i][Index] += 2.0 * this->PseudoPotentials[2][J-this->PseudoPotentialMins[2]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsUpUpUpDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0,1.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
      
      
		    
	  //cout << "Down Down Up Up" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownUpUp[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsDownDownUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[3] ; J < this->NbrPseudoPotentialCoeffs[3] + this->PseudoPotentialMins[3] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownUpUp[i][Index] += this->PseudoPotentials[3][J-this->PseudoPotentialMins[3]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsDownDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,0.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Down Down Down Down" << endl;  
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[4] ; J < this->NbrPseudoPotentialCoeffs[4] + this->PseudoPotentialMins[4] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownDownDown[i][Index] += this->PseudoPotentials[4][J-this->PseudoPotentialMins[4]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsDownDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Down Down Up Down" << endl;
	  for (int i = 1; i < this->NbrUpDownSectorSums-1; ++i) // go through the possible sums of Lz values with one particle on LLL and one on FLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownUpDown[i-1] = new double[this->NbrDownDownSectorIndicesPerSum[i-1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-1]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i-1][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  int m4 = (this->DownDownSectorIndicesPerSum[i-1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsDownDownUpDown[i-1][Index] = 0;
			  for ( J = this->PseudoPotentialMins[5] ; J < this->NbrPseudoPotentialCoeffs[5] + this->PseudoPotentialMins[5] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownUpDown[i-1][Index] += 2.0 * this->PseudoPotentials[5][J-this->PseudoPotentialMins[5]] * TmpCoefficient;
				} 				 
			    }
			  //this->InteractionFactorsDownDownUpDown[i-1][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1/2.0 ,0.0,(double)m2/2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpDown[i-1][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Up Down Up Up" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownUpUp[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+1][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+1][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsUpDownUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[6] ; J < this->NbrPseudoPotentialCoeffs[6] + this->PseudoPotentialMins[6] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpDownUpUp[i][Index] += 2.0 * this->PseudoPotentials[6][J-this->PseudoPotentialMins[6]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsUpDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,1.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownUpUp[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Up Down Down Down" << endl;
	  //now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift ;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+1]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsUpDownDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[7] ; J < this->NbrPseudoPotentialCoeffs[7] + this->PseudoPotentialMins[7] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpDownDownDown[i][Index] += 2.0 * this->PseudoPotentials[7][J-this->PseudoPotentialMins[7]] * TmpCoefficient;
				}
			    }
			  //this->InteractionFactorsUpDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,0.0 ,(double)m1 /2.0 ,0.0 ,(double)m2 /2.0,1.0,(double)m3/2.0,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Up Down Up Down" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownUpDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
			  
			  this->InteractionFactorsUpDownUpDown[i][Index] = 0;
        /*
			  for ( J = this->PseudoPotentialMins[8] ; J < this->NbrPseudoPotentialCoeffs[8] + this->PseudoPotentialMins[8] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpDownUpDown[i][Index] += this->PseudoPotentials[8][J-this->PseudoPotentialMins[8]] * TmpCoefficient;
				}
			    }
        */


        // Now antisymmetrize so there are 4 contributions (two are the same)
                          
        double tmp1234 = 0.0;
        for ( J = this->PseudoPotentialMins[8] ; J < this->NbrPseudoPotentialCoeffs[8] + this->PseudoPotentialMins[8] ; J++ ) 
          {
            if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
        {
          ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);                  
          TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
          tmp1234 += this->PseudoPotentials[8][J-this->PseudoPotentialMins[8]] * TmpCoefficient;
        }
          }

        double tmp1243 = 0.0;
        for ( J = this->PseudoPotentialMins[9] ; J < this->NbrPseudoPotentialCoeffs[9] + this->PseudoPotentialMins[9] ; J++ ) 
          {
            if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
        {
          ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);                  
          TmpCoefficient = ClebschCoef * ClebschDownUp.GetCoefficient(m4, m3, J*2);
          tmp1243 += this->PseudoPotentials[9][J-this->PseudoPotentialMins[9]] * TmpCoefficient;
        }
          }

        double tmp2134 = 0.0;
        for ( J = this->PseudoPotentialMins[9] ; J < this->NbrPseudoPotentialCoeffs[9] + this->PseudoPotentialMins[9] ; J++ ) 
          {
            if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
        {
          ClebschCoef = ClebschDownUp.GetCoefficient(m2,m1,J*2);                  
          TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
          tmp2134 += this->PseudoPotentials[9][J-this->PseudoPotentialMins[9]] * TmpCoefficient;
        }
          }

        double tmp2143 = 0.0;
        for ( J = this->PseudoPotentialMins[8] ; J < this->NbrPseudoPotentialCoeffs[8] + this->PseudoPotentialMins[8] ; J++ ) 
          {
            if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
        {
          ClebschCoef = ClebschDownUp.GetCoefficient(m2,m1,J*2);                  
          TmpCoefficient = ClebschCoef * ClebschDownUp.GetCoefficient(m4, m3, J*2);
          tmp2143 += this->PseudoPotentials[8][J-this->PseudoPotentialMins[8]] * TmpCoefficient;
        }
          }
                            
        this->InteractionFactorsUpDownUpDown[i][Index] = tmp1234 + tmp2134 + tmp1243 + tmp2143;


			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	}
      else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
	{

	  //upup-upup term  
	  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	    {
	      if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
		{
		  this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]];
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
                          
			  this->InteractionFactorsUpUpUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[0] ; J < this->NbrPseudoPotentialCoeffs[0] + this->PseudoPotentialMins[0] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.CarefulGetCoefficient(m1,m2,J*2);			     	
				  //cout << "ClebschCeof: (" << m1 << ", " << m2 << ", " << (J * 1) << ") = " <<  ClebschUpUp.CarefulGetCoefficient(m1,m2,J*2) << endl;
				  TmpCoefficient = ClebschCoef * ClebschUpUp.CarefulGetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpUpUp[i][Index] += FermionicSign * this->PseudoPotentials[0][J-this->PseudoPotentialMins[0]] * TmpCoefficient;
				}
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpUpUpUp[i][Index] = FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ+2,m2,TwoQ+2,m3,TwoQ+2,m4);

			  TotalNbrInteractionFactors++;
			  ++Index;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Down Down" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpUpDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
                          
			  this->InteractionFactorsUpUpDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[1] ; J < this->NbrPseudoPotentialCoeffs[1] + this->PseudoPotentialMins[1] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpUp.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpDownDown[i][Index] += FermionicSign * this->PseudoPotentials[1][J-this->PseudoPotentialMins[1]] * TmpCoefficient;
				}
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpUpDownDown[i][Index] = FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ+2,m2,TwoQ,m3,TwoQ,m4);

			  //this->InteractionFactorsUpUpDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/ 2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0 ,1.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Up Down" << endl; Extra factor of 2 in the matrix element because of antisymmetrization!
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpUpUpDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  
			  this->InteractionFactorsUpUpUpDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[2] ; J < this->NbrPseudoPotentialCoeffs[2] + this->PseudoPotentialMins[2] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpUp.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpUpUpDown[i][Index] += 2.0 * FermionicSign * this->PseudoPotentials[2][J-this->PseudoPotentialMins[2]] * TmpCoefficient;
				}
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpUpUpDown[i][Index] = 2.0 * FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ+2,m2,TwoQ+2,m3,TwoQ,m4);

			  //this->InteractionFactorsUpUpUpDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0,1.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
      
     		    
	  //cout << "Down Down Up Up" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownUpUp[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i][j2 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsDownDownUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[3] ; J < this->NbrPseudoPotentialCoeffs[3] + this->PseudoPotentialMins[3] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownUpUp[i][Index] += FermionicSign * this->PseudoPotentials[3][J-this->PseudoPotentialMins[3]] * TmpCoefficient;
				}
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsDownDownUpUp[i][Index] = FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ,m1,TwoQ,m2,TwoQ+2,m3,TwoQ+2,m4);

			  //this->InteractionFactorsDownDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,0.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	  // Down Down Down Down  
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i][j2 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1 ;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsDownDownDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[4] ; J < this->NbrPseudoPotentialCoeffs[4] + this->PseudoPotentialMins[4] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownDownDown[i][Index] += FermionicSign * this->PseudoPotentials[4][J-this->PseudoPotentialMins[4]] * TmpCoefficient;
				}
			    }
			  
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsDownDownDownDown[i][Index] = FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ,m1,TwoQ,m2,TwoQ,m3,TwoQ,m4);

			  //this->InteractionFactorsDownDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	    //cout << "Down Down Up Down" << endl;Extra factor of 2 in the matrix element because of antisymmetrization!
	    for (int i = 2; i < this->NbrDownDownSectorSums + 2; ++i) // go through the possible sums of Lz values with one particle on LLL and one on FLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsDownDownUpDown[i-2] = new double[this->NbrDownDownSectorIndicesPerSum[i-2] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-2]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i-2][j2 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i-2][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsDownDownUpDown[i-2][Index] = 0;
			  for ( J = this->PseudoPotentialMins[5] ; J < this->NbrPseudoPotentialCoeffs[5] + this->PseudoPotentialMins[5] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsDownDownUpDown[i-2][Index] += 2.0 * FermionicSign * this->PseudoPotentials[5][J-this->PseudoPotentialMins[5]] * TmpCoefficient;
				} 				 
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsDownDownUpDown[i-2][Index] = 2.0 * FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ,m1,TwoQ,m2,TwoQ+2,m3,TwoQ,m4);
 
			  //this->InteractionFactorsDownDownUpDown[i-1][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1/2.0 ,0.0,(double)m2/2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpDown[i-1][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	    //cout << "Up Down Up Up" << endl;Extra factor of 2 in the matrix element because of antisymmetrization!
	    for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownUpUp[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsUpDownUpUp[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[6] ; J < this->NbrPseudoPotentialCoeffs[6] + this->PseudoPotentialMins[6] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpUp.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpDownUpUp[i][Index] += 2.0 * FermionicSign * this->PseudoPotentials[6][J-this->PseudoPotentialMins[6]] * TmpCoefficient;
				}
			    }
                          
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpDownUpUp[i][Index] = 2.0 * FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ,m2,TwoQ+2,m3,TwoQ+2,m4);

			  //this->InteractionFactorsUpDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,1.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownUpUp[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	    //cout << "Up Down Down Down" << endl; Extra factor of 2 in the matrix element because of antisymmetrization!
	    //now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
	    for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i+2][(j2 << 1) + 1] + this->LzFermionDownShift )<< 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsUpDownDownDown[i][Index] = 0;
			  for ( J = this->PseudoPotentialMins[7] ; J < this->NbrPseudoPotentialCoeffs[7] + this->PseudoPotentialMins[7] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  this->InteractionFactorsUpDownDownDown[i][Index] += 2.0 * FermionicSign * this->PseudoPotentials[7][J-this->PseudoPotentialMins[7]] * TmpCoefficient;
				}
			    }
			  
                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpDownDownDown[i][Index] = 2.0 * FermionicSign * this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ,m2,TwoQ,m3,TwoQ,m4);
 
			  //this->InteractionFactorsUpDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,0.0 ,(double)m1 /2.0 ,0.0 ,(double)m2 /2.0,1.0,(double)m3/2.0,0.0,(double)m4/2.0);
			  //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	    
	    //cout << "Up Down Up Down" << endl; This term is tricky; it contains 4 different terms (1,0,1,0), (1,0,0,1), (0,1,0,1), (0,1,1,0), so it will require both Pseudopotentials[8] and Pseudopotentials[9] to evaluate!
	    for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  this->InteractionFactorsUpDownUpDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  
			  this->InteractionFactorsUpDownUpDown[i][Index] = 0;

			  // Now antisymmetrize so there are 4 contributions (two are the same)
                          
                          double tmp1234 = 0.0;
			  for ( J = this->PseudoPotentialMins[8] ; J < this->NbrPseudoPotentialCoeffs[8] + this->PseudoPotentialMins[8] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  tmp1234 += FermionicSign * this->PseudoPotentials[8][J-this->PseudoPotentialMins[8]] * TmpCoefficient;
				}
			    }

                          double tmp1243 = 0.0;
			  for ( J = this->PseudoPotentialMins[9] ; J < this->NbrPseudoPotentialCoeffs[9] + this->PseudoPotentialMins[9] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschUpDown.GetCoefficient(m1,m2,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownUp.GetCoefficient(m4, m3, J*2);
				  tmp1243 += FermionicSign * this->PseudoPotentials[9][J-this->PseudoPotentialMins[9]] * TmpCoefficient;
				}
			    }

                          double tmp2134 = 0.0;
			  for ( J = this->PseudoPotentialMins[9] ; J < this->NbrPseudoPotentialCoeffs[9] + this->PseudoPotentialMins[9] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownUp.GetCoefficient(m2,m1,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschUpDown.GetCoefficient(m3, m4, J*2);
				  tmp2134 += FermionicSign * this->PseudoPotentials[9][J-this->PseudoPotentialMins[9]] * TmpCoefficient;
				}
			    }

                          double tmp2143 = 0.0;
			  for ( J = this->PseudoPotentialMins[8] ; J < this->NbrPseudoPotentialCoeffs[8] + this->PseudoPotentialMins[8] ; J++ ) 
			    {
			      if ( abs(m1 + m2) <= (J*2) && abs(m3 + m4) <= (J*2) ) 
				{
				  ClebschCoef = ClebschDownUp.GetCoefficient(m2,m1,J*2);			     			  
				  TmpCoefficient = ClebschCoef * ClebschDownUp.GetCoefficient(m4, m3, J*2);
				  tmp2143 += FermionicSign * this->PseudoPotentials[8][J-this->PseudoPotentialMins[8]] * TmpCoefficient;
				}
			    }
                            
			  this->InteractionFactorsUpDownUpDown[i][Index] = tmp1234 - tmp2134 - tmp1243 + tmp2143;
                          

                          //int TwoQ = this->LzMax - 2;
                          //this->InteractionFactorsUpDownUpDown[i][Index] = FermionicSign * (this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ,m2,TwoQ+2,m3,TwoQ,m4) - this->CalculateLaplacianDeltaInteractionFactor(TwoQ,m2,TwoQ+2,m1,TwoQ+2,m3,TwoQ,m4) - this->CalculateLaplacianDeltaInteractionFactor(TwoQ+2,m1,TwoQ,m2,TwoQ,m4,TwoQ+2,m3) + this->CalculateLaplacianDeltaInteractionFactor(TwoQ,m2,TwoQ+2,m1,TwoQ,m4,TwoQ+2,m3));
  
			  ++Index;
			  TotalNbrInteractionFactors++;
			}
		    }
		}
	    }
	}	      	
    }
  else // otherwise the interaction factors for the delta interaction are set explicitly. 
    {
      cout << "Delta interaction not tested... provide pseudopotentials if you require delta interaction." << endl;
/*
      //when not using pseudopotentials will use expression that was worked out for the delta interaction.
      
      //will start with the interactions confined completely to the LLL (Down).
      //cout << "Down Down Down Down" << endl;
      for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsDownDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrDownDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsDownDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Up Down Down Down" << endl;
	//now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
	for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpDownDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+1]; ++j2)
		    {
		      int m3 = (this->UpDownSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpDownSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsUpDownDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,0.0 ,(double)m1 /2.0 ,0.0 ,(double)m2 /2.0,1.0,(double)m3/2.0,0.0,(double)m4/2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Up Up Down Down" << endl;
	for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpUpDownDown[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
		    {
		      int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		      
		      this->InteractionFactorsUpUpDownDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/ 2.0 ,0.0,(double)m1/2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0 ,1.0,(double)m4/2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Down Down Up Down" << endl;
	for (int i = 1; i < this->NbrUpDownSectorSums-1; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsDownDownUpDown[i-1] = new double[this->NbrDownDownSectorIndicesPerSum[i-1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-1]; ++j2)
		    {
		      int m3 = (this->DownDownSectorIndicesPerSum[i-1][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      int m4 = (this->DownDownSectorIndicesPerSum[i-1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsDownDownUpDown[i-1][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1/2.0 ,0.0,(double)m2/2.0 ,0.0,(double)m3/2.0 ,0.0,(double)m4/2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpDown[i-1][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Up Down Up Down" << endl;
	for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpDownUpDown[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsUpDownUpDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownUpDown[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}

	//cout << "Up Up Up Down" << endl;
	for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpUpUpDown[i] = new double[this->NbrUpUpSectorIndicesPerSum[i+1] * this->NbrUpDownSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		  for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j2)
		    {
		      int m3 = (this->UpUpSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpUpSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxUp;
		      
		      this->InteractionFactorsUpUpUpDown[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,0.0,(double)m2 /2.0,1.0,(double)m3 /2.0,1.0,(double)m4/2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Down Down Up Up" << endl;
	for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsDownDownUpUp[i] = new double[this->NbrDownDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+2]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
		{
		  int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		  for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsDownDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1) /2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,0.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	
	//cout << "Up Down Up Up" << endl;
	for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpDownUpUp[i] = new double[this->NbrUpDownSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i+1]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j1)
		{
		  int m1 = (this->UpUpSectorIndicesPerSum[i+1][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpUpSectorIndicesPerSum[i+1][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		  for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1 + this->LzFermionDownShift;
		      
		      this->InteractionFactorsUpDownUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,1.0,(double)m3 /2.0,0.0,(double)m4 /2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1 << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpDownUpUp[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
	  
	//cout << "Up Up Up Up" << endl;    
	for (int i = 0; i < this->NbrUpUpSectorSums; ++i) // go through the possible sums of Lz values on LLL
	{
	  if (this->NbrUpUpSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
	    {
	      this->InteractionFactorsUpUpUpUp[i] = new double[this->NbrUpUpSectorIndicesPerSum[i] * this->NbrUpUpSectorIndicesPerSum[i]]; //for all m1, m2, m3, m4 such that m1 + m2 = m3 + m4 = current_sum
	      int Index = 0;
	      for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		{
		  int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		  int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		  for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
		    {
		      int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
		      int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp ;
		      
		      this->InteractionFactorsUpUpUpUp[i][Index] = this->CalculateDeltaInteractionFactor((double)(this->LzMaxDown-1)/2.0 ,1.0,(double)m1 /2.0 ,1.0,(double)m2 /2.0,1.0,(double)m3 /2.0,1.0,(double)m4 / 2.0);
		      //cout << this->LzMaxDown-1.0 << ": " << m1  << ", " << m2 << ", " << m3 << ", " << m4 << ": " << this->InteractionFactorsUpUpUpUp[i][Index] << endl;
		      ++Index;
		      TotalNbrInteractionFactors++;
		    }
		}
	    }
	}
*/	    
    }
    
  if ( this->ShowIntFactorsFlag )
    {
      if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
        {
	  //upup-upup term  
	  cout << "UpUpUpUp Terms" << endl ;
	  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	    {	   
	      if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpUpUp[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpUp[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Down Down" << endl;
	  cout << "UpUpDownDown Terms" << endl ;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpDownDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
		    
	  cout << "UpUpUpDown Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpUpDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
      
      
		    
	  cout << "DownDownUpUp Terms" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	    
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1;
			  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownUpUp[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "DownDownDownDown Terms" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxDown - 1;
			  int m4 = (this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownDownDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "DownDownUpDown Terms" << endl;
	  for (int i = 1; i < this->NbrUpDownSectorSums-1; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 ;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-1]; ++j2)
			{
			  int m3 = (this->DownDownSectorIndicesPerSum[i-1][j2 << 1] << 1) - this->LzMaxDown - 1;
			  int m4 = (this->DownDownSectorIndicesPerSum[i-1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownUpDown[i-1][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownUpDown[i-1][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "UpDownUpUp Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+1]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+1][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+1][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsUpUpUpDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	cout << "UpDownDownDown Terms" << endl;
	//now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
	for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->DownDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxDown - 1;
		      int m2 = (this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+1]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i+1][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i+1][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsUpDownDownDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "UpDownUpDown Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum. 
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxDown - 1 ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxDown - 1;		      
			  if (this->InteractionFactorsUpDownUpDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpDownUpDown[i][Index] << endl;	
			  ++Index;
			}
		    }
		}
	    }
	}
      else if ( Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic )
        {

	   //upup-upup term  
	  cout << "UpUpUpUp Terms" << endl ;
	  for (int i = 0; i < this->NbrUpUpSectorSums; ++i)
	    {	   
	      if (this->NbrUpUpSectorIndicesPerSum[i] > 0)
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpUpUp[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpUp[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  //cout << "Up Up Down Down" << endl;
	  cout << "UpUpDownDown Terms" << endl ;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1]  + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1]  + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i+2][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpDownDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
		    
	  cout << "UpUpUpDown Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
		      for (int j2 = 0; j2 < this->NbrUpUpSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpUpSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = (this->UpUpSectorIndicesPerSum[i][(j2 << 1) + 1] << 1) - this->LzMaxUp;
			  if (this->InteractionFactorsUpUpUpDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
      
	  cout << "DownDownUpUp Terms" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	    
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i+2]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i+2][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i+2][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i][j2 << 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownUpUp[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownUpUp[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }

    
	  cout << "DownDownDownDown Terms" << endl;
	  for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{	      
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1] + this->LzFermionDownShift )  << 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift )<< 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i][j2 << 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownDownDown[i][Index]  != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "DownDownUpDown Terms" << endl;
	  for (int i = 2; i < this->NbrUpDownSectorSums-2; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
		      for (int j2 = 0; j2 < this->NbrDownDownSectorIndicesPerSum[i-2]; ++j2)
			{
			  int m3 = ((this->DownDownSectorIndicesPerSum[i-2][j2 << 1]  + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  int m4 = ((this->DownDownSectorIndicesPerSum[i-2][(j2 << 1) + 1]  + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsDownDownUpDown[i-2][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsDownDownUpDown[i-2][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "UpDownUpUp Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpUpSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpUpSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = (this->UpUpSectorIndicesPerSum[i][(j1 << 1) + 1] << 1) - this->LzMaxUp ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1]+ this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
			  if (this->InteractionFactorsUpUpUpDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpUpUpDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	cout << "UpDownDownDown Terms" << endl;
	//now we set the interaction terms where a single operator acts on the first LL and the rest on the LLL 
	for (int i = 0; i < this->NbrDownDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrDownDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum.
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrDownDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = ((this->DownDownSectorIndicesPerSum[i][j1 << 1] + this->LzFermionDownShift )<< 1) - this->LzMaxDown - 1;
		      int m2 = ((this->DownDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1 ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i+2]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i+2][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i+2][(j2 << 1) + 1] + this->LzFermionDownShift) << 1) - this->LzMaxDown - 1;
			  if (this->InteractionFactorsUpDownDownDown[i][Index] != 0)
				  cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpDownDownDown[i][Index] << endl;
			  ++Index;
			}
		    }
		}
	    }
	    
	  cout << "UpDownUpDown Terms" << endl;
	  for (int i = 0; i < this->NbrUpDownSectorSums; ++i) // go through the possible sums of Lz values on LLL
	    {
	      if (this->NbrUpDownSectorIndicesPerSum[i] > 0) // if there are m1 and m2 values that give this sum. 
		{
		  int Index = 0;
		  for (int j1 = 0; j1 < this->NbrUpDownSectorIndicesPerSum[i]; ++j1)
		    {
		      int m1 = (this->UpDownSectorIndicesPerSum[i][j1 << 1] << 1) - this->LzMaxUp ;
		      int m2 = ((this->UpDownSectorIndicesPerSum[i][(j1 << 1) + 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
		      for (int j2 = 0; j2 < this->NbrUpDownSectorIndicesPerSum[i]; ++j2)
			{
			  int m3 = (this->UpDownSectorIndicesPerSum[i][j2 << 1] << 1) - this->LzMaxUp;
			  int m4 = ((this->UpDownSectorIndicesPerSum[i][(j2 << 1) + 1] + this->LzFermionDownShift ) << 1) - this->LzMaxDown - 1 ;
			  if (this->InteractionFactorsUpDownUpDown[i][Index] != 0)			      
			  	cout << "<" << (double)m1/2.0 << ", " << (double)m2/2.0 << "| V | " << (double)m3/2.0 << ", " << (double)m4/2.0 << "> = " << this->InteractionFactorsUpDownUpDown[i][Index] << endl;	
			  ++Index;
			}
		    }
		}
	    }
	
	}
    }

  this->NbrOneBodyInteractionFactorsUpUp = this->LzMaxUp + 1;
  this->NbrOneBodyInteractionFactorsUpDown = 0;
  this->NbrOneBodyInteractionFactorsDownUp = 0;
  this->NbrOneBodyInteractionFactorsDownDown = this->LzMaxDown;
  this->OneBodyInteractionFactorsUpUp = new double [this->NbrOneBodyInteractionFactorsUpUp];
  this->OneBodyInteractionFactorsUpDown = 0;
  this->OneBodyInteractionFactorsDownUp = 0;
  this->OneBodyInteractionFactorsDownDown = new double [this->NbrOneBodyInteractionFactorsDownDown]; 
  this->OneBodyMValuesUpUp = new int[this->NbrOneBodyInteractionFactorsUpUp];
  this->OneBodyMValuesUpDown = 0;
  this->OneBodyMValuesDownUp = 0;
  this->OneBodyMValuesDownDown = new int[this->NbrOneBodyInteractionFactorsDownDown];
  for (int i = 0; i <= this->LzMaxUp; ++i)
    {
      this->OneBodyMValuesUpUp[i] = i;
      this->OneBodyInteractionFactorsUpUp[i] = this->TotalCyclotronEnergy[1];
    }

  for (int i = 1; i <= this->LzMaxDown; ++i)
    {
      this->OneBodyMValuesDownDown[i-1] = i - this->LzFermionDownShift;
      this->OneBodyInteractionFactorsDownDown[i-1] = this->TotalCyclotronEnergy[0];
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// evaluate a particular delta interaction factor  <Q,l1,m1; Q,l2,m2|V|Q,l3,m3; Q,l4,m4> in front of c_{l1,m1}^+ c_{l2,m2}^+ c_{l4,m4} c_{l3,m3}
//
// <Q,l1,m1; Q,l2,m2|V_C|Q,l3,m3; Q,l4,m4> = Sum[-L(L+1)/Q (-1)^(2Q-j1-j2-M)SL[-Q,l1,-j1,0,L,-M,Q,l3,j3]SL[-Q,l2,-j2,0,L,M,Q,l4,j4],{L,Max[Abs[l1-l3],Abs[l2-l4]],Min[l1+l3,l2+l4]}]
//
// Q = max angular momentum on LLL.
// l1 = LL of m1, 0 or 1 supported.
// m1 = angular momentum of first operator.
// l2 = LL of m2, 0 or 1 supported.
// m2 = angular momentum of second operator.
// l3 = LL of m3, 0 or 1 supported.
// m3 = angular momentum of third operator.
// l4 = LL of m4, 0 or 1 supported.	
// m4 = angular momentum of fourth operator.
//
// return value = interaction factor  <Q,l1,m1; Q,l2,m2|V|Q,l3,m3; Q,l4,m4>

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateLaplacianDeltaInteractionFactor(int Twol1,int Twom1,int Twol2,int Twom2,int Twol3,int Twom3,int Twol4,int Twom4)
{
  int TwoQ = (this->LzMax-2);
  int TwoM = Twom3 - Twom1;
  double tmpfactor = 1/(0.5*TwoQ);
  int TwoLmin = abs(Twol1-Twol3);
  if (abs(Twol2-Twol4) > TwoLmin)
    TwoLmin = abs(Twol2-Twol4); 
  int TwoLmax = Twol1+Twol3; 
  if (TwoLmax > (Twol2+Twol4)) 
    TwoLmax = Twol2+Twol4;

  double tmpsum = 0.0;
  int TwoL = TwoLmin;
  while (TwoL <= TwoLmax)
   {
     tmpsum += tmpfactor * (-0.5*TwoL*(0.5*TwoL+1)) * pow(-1.0, TwoQ - 0.5*Twom1 - 0.5*Twom2 - 0.5*TwoM) * CalculateSLFactor(-TwoQ,Twol1,-Twom1,0,TwoL,-TwoM,TwoQ,Twol3,Twom3) * CalculateSLFactor(-TwoQ,Twol2,-Twom2,0,TwoL,TwoM,TwoQ,Twol4,Twom4);
     TwoL=TwoL+2;
   }

 //cout<<"Interaction factor " << Twol1<<" "<<Twom1<< "; "<<Twol2<<" "<<Twom2<< "; "<<Twol3<<" "<<Twom3<< "; "<<Twol4<<" "<<Twom4<< " =  "<<tmpsum << endl; 
 return tmpsum; 
}

// evaluate a particular Coulomb interaction factor  <Q,l1,m1; Q,l2,m2|V_C|Q,l3,m3; Q,l4,m4> in front of c_{l1,m1}^+ c_{l2,m2}^+ c_{l4,m4} c_{l3,m3}
//
// <Q,l1,m1; Q,l2,m2|V_C|Q,l3,m3; Q,l4,m4> = Sum[4Pi/Sqrt[Q] 1/(2L+1)(-1)^(2Q-j1-j2-M)SL[-Q,l1,-j1,0,L,-M,Q,l3,j3]SL[-Q,l2,-j2,0,L,M,Q,l4,j4],{L,Max[Abs[l1-l3],Abs[l2-l4]],Min[l1+l3,l2+l4]}]
//
// Q = max angular momentum on LLL.
// l1 = LL of m1, 0 or 1 supported.
// m1 = angular momentum of first operator.
// l2 = LL of m2, 0 or 1 supported.
// m2 = angular momentum of second operator.
// l3 = LL of m3, 0 or 1 supported.
// m3 = angular momentum of third operator.
// l4 = LL of m4, 0 or 1 supported.	
// m4 = angular momentum of fourth operator.
//
// return value = interaction factor  <Q,l1,m1; Q,l2,m2|V_C|Q,l3,m3; Q,l4,m4>

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateCoulombInteractionFactor(int Twol1,int Twom1,int Twol2,int Twom2,int Twol3,int Twom3,int Twol4,int Twom4)
{
  int TwoQ = (this->LzMax-2);
  int TwoM = Twom3 - Twom1;
  double tmpfactor = 4.0 * M_PI/sqrt(0.5*TwoQ);
  int TwoLmin = abs(Twol1-Twol3);
  if (abs(Twol2-Twol4) > TwoLmin)
    TwoLmin = abs(Twol2-Twol4); 
  int TwoLmax = Twol1+Twol3; 
  if (TwoLmax > (Twol2+Twol4)) 
    TwoLmax = Twol2+Twol4;

  double tmpsum = 0.0;
  int TwoL = TwoLmin;
  while (TwoL <= TwoLmax)
   {
     tmpsum += tmpfactor * (1.0/(TwoL+1)) * pow(-1.0, TwoQ - 0.5*Twom1 - 0.5*Twom2 - 0.5*TwoM) * CalculateSLFactor(-TwoQ,Twol1,-Twom1,0,TwoL,-TwoM,TwoQ,Twol3,Twom3) * CalculateSLFactor(-TwoQ,Twol2,-Twom2,0,TwoL,TwoM,TwoQ,Twol4,Twom4);
     TwoL=TwoL+2;
   }

 //cout<<"Interaction factor " << Twol1<<" "<<Twom1<< "; "<<Twol2<<" "<<Twom2<< "; "<<Twol3<<" "<<Twom3<< "; "<<Twol4<<" "<<Twom4<< " =  "<<tmpsum << endl; 
 return tmpsum; 
}

// evaluate a particular SL interaction factor needed by CalculateCoulombInteractionFactor
//
// SL(Q1,l1,m1;Q2,l2,m2;Q3=-Q1-Q2,l3,m3=-m1-m3) = (-1)^(l1+l2+l3)Sqrt[(2l1+1)(2l2+1)(2l3+1)/(4Pi)]ThreeJSymbol[{l1,-alpha1},{l2,-alpha2},{l3,-alpha3}]ThreeJSymbol[{l1,Q1},{l2,Q2},{l3,Q3}]
//
// Q = max angular momentum on LLL.
// l1 = LL of m1, 0 or 1 supported.
// m1 = angular momentum of first operator.
// l2 = LL of m2, 0 or 1 supported.
// m2 = angular momentum of second operator.
// l3 = LL of m3, 0 or 1 supported.
// m3 = angular momentum of third operator.
// l4 = LL of m4, 0 or 1 supported.	
// m4 = angular momentum of fourth operator.
//
// return value = interaction factor with delta interaction

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateSLFactor(int TwoQ1,int Twol1,int Twom1,int TwoQ2,int Twol2,int Twom2,int TwoQ3,int Twol3,int Twom3)
{
   if ((TwoQ1+TwoQ2+TwoQ3) != 0)
      return 0.0;
   if ((Twom1+Twom2+Twom3) != 0)
      return 0.0;
   if ( (Twol3>(Twol1+Twol2)) or (Twol3<abs(Twol1-Twol2)) or (abs(Twom1)>Twol1) or (abs(Twom2)>Twol2) or (abs(Twom3)>Twol3) or (abs(TwoQ1)>Twol1) or (abs(TwoQ2)>Twol2) or (abs(TwoQ3)>Twol3))
	return 0.0;
   
   double tmpnorm = pow(-1.0, 0.5*Twol1+0.5*Twol2+0.5*Twol3) * sqrt((Twol1+1)*(Twol2+1)*(Twol3+1)/(4.0*M_PI));
  
   ThreeJSymbol ThreeJ(Twol1, Twol2);

   double tmp3j1 = ThreeJ.GetCoefficient (-Twom1, -Twom2, Twol3);
   double tmp3j2 = ThreeJ.GetCoefficient (TwoQ1, TwoQ2, Twol3);

   return (tmpnorm * tmp3j1 * tmp3j2);
}

/*

//ZP: Untested! Provide pseudopotential if you need delta interaction...
//
// evaluate a particular interaction factor <l1,m1,l2,m2|\delta|l3,m3,l4.m4>
//
// Q = max angular momentum on LLL.
// l1 = LL of m1, 0 or 1 supported.
// m1 = angular momentum of first operator.
// l2 = LL of m2, 0 or 1 supported.
// m2 = angular momentum of second operator.
// l3 = LL of m3, 0 or 1 supported.
// m3 = angular momentum of third operator.
// l4 = LL of m4, 0 or 1 supported.	
// m4 = angular momentum of fourth operator.
//
// return value = interaction factor with delta interaction

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateDeltaInteractionFactor(double Q,double l1,double m1, double l2, double m2, double l3, double m3, double l4, double m4)
{
   double NormCoeff = 1.0;
   double Value = 0.0;
   
   //work out normalisation.
   NormCoeff *= ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateNormalization(Q, Q+l1, m1)*pow(-1.0,(double)l1);
   NormCoeff *= ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateNormalization(Q, Q+l2, m2)*pow(-1.0,(double)l2);
   NormCoeff *= ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateNormalization(Q, Q+l3, m3)*pow(-1.0,(double)l3);
   NormCoeff *= ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateNormalization(Q, Q+l4, m4)*pow(-1.0,(double)l4);
    
   NormCoeff *= 4.0*M_PI;
   
   for ( int S1 = 0 ; S1 <= l1 ; S1++ ) 
    {
      double S1Coeff = 0.0;
      if ( (l1 + Q - m1 - S1) >= 0 && (Q + m1 + S1) >= 0 ) 
	{
	  FactorialCoefficient* S1Factorial = new FactorialCoefficient(); 
	  S1Factorial->FactorialMultiply((long)(l1+2*Q));
	  S1Factorial->FactorialDivide((long)(l1 + Q - m1 - S1));
	  S1Factorial->FactorialDivide((long)(Q + m1 + S1));
	  S1Coeff = pow(-1.0,(double)S1)*S1Factorial->GetNumericalValue();
	}
	  
      for ( int S2 = 0 ; S2 <= l2 ; S2++ )
        {
	  double S2Coeff = 0.0;
	  if ( (l2 + Q - m2 - S2) >= 0 && (Q + m2 + S2) >= 0 ) 
	    {
	      FactorialCoefficient* S2Factorial = new FactorialCoefficient(); 
	      S2Factorial->FactorialMultiply((long)(l2+2*Q));
	      S2Factorial->FactorialDivide((long)(l2 + Q - m2 - S2));
	      S2Factorial->FactorialDivide((long)(Q + m2 + S2));
	      S2Coeff = pow(-1.0,(double)S2)*S2Factorial->GetNumericalValue();
	    }
          for ( int S3 = 0 ; S3 <= l3 ; S3++ )
            {
	      double S3Coeff = 0.0;
              if ( (l3 + Q - m3 - S3) >= 0 && (Q + m3 + S3) >= 0 ) 
		{
		  FactorialCoefficient* S3Factorial = new FactorialCoefficient(); 
		  S3Factorial->FactorialMultiply((long)(l3+2*Q));
		  S3Factorial->FactorialDivide((long)(l3 + Q - m3 - S3));
		  S3Factorial->FactorialDivide((long)(Q + m3 + S3));
		  S3Coeff = pow(-1.0,(double)S3)*S3Factorial->GetNumericalValue();
		}
		  
              for ( int S4 = 0 ; S4 <= l4 ; S4++ ) 
                {
		    double S4Coeff = 0.0;
                    if ( (l4 + Q - m4 - S4) >= 0 && (Q + m4 + S4) >= 0 ) 
		      {
			FactorialCoefficient* S4Factorial = new FactorialCoefficient(); 
			S4Factorial->FactorialMultiply((long)(l4+2*Q));
			S4Factorial->FactorialDivide((long)(l4 + Q - m4 - S4));
			S4Factorial->FactorialDivide((long)(Q + m4 + S4));
			S4Coeff = pow(-1.0,(double)S4)*S4Factorial->GetNumericalValue();
		      }
                    Value += NormCoeff * S1Coeff * S2Coeff * S3Coeff * S4Coeff * ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateBetaFunction((long)(2*Q - m1 - m2 + 1 + l1 + l2 + l3 + l4 - S1 - S2 - S3 - S4),(long)(2*Q + m1 + m2 + 1 + S1 + S2 + S3 + S4));
                }  
            }
        }
    }
  return Value;
}

// evaluate normalisation for Q, l, m
//
// Q = max angular momentum on LLL.
// l = Q + LL 
// m = angular momentum
//
// return value = normalisation

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateNormalization(double Q,double l,double m)
{
  double Value = (2.0*(double)l + 1.0)/(4.0*M_PI);
  FactorialCoefficient* MyFactorial = new FactorialCoefficient(); 
  MyFactorial->FactorialMultiply((long)(l-m));
  MyFactorial->FactorialMultiply((long)(l+m));
  MyFactorial->FactorialDivide((long)(l-Q));
  MyFactorial->FactorialDivide((long)(l+Q));
  return sqrt(Value*MyFactorial->GetNumericalValue());
}

// evaluate the beta function B(x,y) = (x-1)!(y-1)!/(x+y-1)!
//
// x = first arg
// y = second arg
//
// return value = value of beta function with args x,y

double ParticleOnSphereTwoLandauLevelDeltaHamiltonian::CalculateBetaFunction(long x, long y)
{
  FactorialCoefficient* MyFactorial = new FactorialCoefficient(); 
  MyFactorial->FactorialMultiply(x-1);
  MyFactorial->FactorialMultiply(y-1);
  MyFactorial->FactorialDivide(x+y-1);
  return MyFactorial->GetNumericalValue();
}
*/