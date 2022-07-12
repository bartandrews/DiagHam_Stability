////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
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


#ifndef PARTICLEONSPHERETWOLANDAULEVELDELTAHAMILTONIAN_H
#define PARTICLEONSPHERETWOLANDAULEVELDELTAHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinFullHamiltonian.h"

#include <iostream>


using std::ostream;


class ParticleOnSphereTwoLandauLevelDeltaHamiltonian : public AbstractQHEOnSphereWithSpinFullHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  // first index refered to the spin sector (sorted as up-up, down-down, up-down)
  double** PseudoPotentials;
  
  // array to store the nubmer of pseudopotential coefficients expected for each combination of uuuu, uudd....
  int *NbrPseudoPotentialCoeffs;
  
  // array to store the lowest value of total angular momentum for each case. this can be either 0 or 1
  int *PseudoPotentialMins;
  
  // bool value which signifies whether or not pseudo potentials should be used.
  bool UsePseudoPotentials;

  double* TotalCyclotronEnergy;
  
  int LandauLevelIndexDifference;
  
  // twice the maximum momentum a particle can reach in the upper Landau level
  int LzMaxUp;
  
  // twice the maximum momentum a particle can reach in the lower Landau level
  int LzMaxDown;
    
  // Flag whcih indicates whether or not to print the values of the interaction factors for debugging purposes
  bool ShowIntFactorsFlag;
  
  // shift to compensate for shift on down level of fermion 2LL class.
  int LzFermionDownShift;
  
  // shift to compensate for shift on up level of fermion 2LL class.
  int LzFermionUpShift;
  
  
 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state in the lower Landau level
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // cyclotronEnergy = cyclotron energy in e^2/(epsilon l_b) unit
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereTwoLandauLevelDeltaHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
                                                                                     double** pseudoPotential, double *cyclotronEnergy,
                                                                                     AbstractArchitecture* architecture,long memory, bool onDiskCacheFlag, char* precalculationFileName, bool showIntFactorsFlag);
										     
  // destructor
  //
  ~ParticleOnSphereTwoLandauLevelDeltaHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

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

  double CalculateLaplacianDeltaInteractionFactor(int Twol1,int Twom1,int Twol2,int Twom2,int Twol3,int Twom3,int Twol4,int Twom4);

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

  double CalculateCoulombInteractionFactor(int Twol1,int Twom1,int Twol2,int Twom2,int Twol3,int Twom3,int Twol4,int Twom4);

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

  double CalculateSLFactor(int TwoQ1,int Twol1,int Twom1,int TwoQ2,int Twol2,int Twom2,int TwoQ3,int Twol3,int Twom3);

  /*  
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
  static double CalculateDeltaInteractionFactor(double Q,double l1,double m1, double l2, double m2, double l3, double m3, double l4, double m4);
  
  // evaluate normalisation for Q, l, m
  //
  // Q = max angular momentum on LLL.
  // l = LL
  // m = angular momentum
  //
  // return value = normalisation
  static double CalculateNormalization(double Q,double l,double m);
  
  // evaluate the beta function B(x,y) = (x-1)!(y-1)!/(x+y-1)!
  //
  // x = first arg
  // y = second arg
  //
  // return value = value of beta function with args x,y
  static double CalculateBetaFunction(long x, long y);

  */
 
 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();

  // core part of the AddMultiply method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added  
  void EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination);
  
  // core part of the PartialFastMultiplicationMemory method involving two-body term
  // 
  // particles = pointer to the Hilbert space
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // memory = reference on the amount of memory required for precalculations
  void EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory);
  
  // core part of the FastMultiplication method involving the two-body interaction
  // 
  // particles = pointer to the Hilbert space
  // index = index of the component on which the Hamiltonian has to act on
  // indexArray = array where indices connected to the index-th component through the Hamiltonian
  // coefficientArray = array of the numerical coefficients related to the indexArray
  // position = reference on the current position in arrays indexArray and coefficientArray
  void EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position);
};


// core part of the AddMultiply method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added

inline void ParticleOnSphereTwoLandauLevelDeltaHamiltonian::EvaluateMNTwoBodyAddMultiplyComponent(ParticleOnSphereWithSpin* particles, int index, RealVector& vSource, RealVector& vDestination)
{
  int Dim = particles->GetHilbertSpaceDimension();
  double Coefficient;
  double Coefficient3;;
  double* TmpInteractionFactor;
  int Index;

  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {  
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownUpUp
		  if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpUp[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
			    ++TmpInteractionFactor;
			  }	  
		    }
				
			  // now UpDownUpUp
		  if ( j >= 1 && j < (this->NbrUpDownSectorSums + 1) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsUpDownUpUp[j-1] + (i * this->NbrUpDownSectorIndicesPerSum[j-1]);	
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j-1] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j-1][k << 1], this->UpDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
			    ++TmpInteractionFactor;
			  }	  
		    }
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpUpDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+1]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+1] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+1][k << 1], this->UpUpSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
		      ++TmpInteractionFactor;
		    }
		      
		    // now DownDownUpDown
		    if ( j >= 1 && j < (this->NbrDownDownSectorSums + 1) ) 
		      {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpDown[j-1] + (i * this->NbrDownDownSectorIndicesPerSum[j-1]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-1] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-1][k << 1], this->DownDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
			      }
			    ++TmpInteractionFactor;
			  }	  
		      }
				  
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
	
	    
	    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpDownDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{				      
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			 
			}
		      ++TmpInteractionFactor;
		    }	  

				    
		  // now UpDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpDownDownDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j+1]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+1] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+1][k << 1], this->UpDownSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
			}
		      ++TmpInteractionFactor;
		    }	  
		}
	    }
	}

    } 
  else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {  
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			  {
			     vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
                             //cout<<"i= "<<index<<"| "; particles->PrintState(cout,index); cout<<"a+a+aa "<< this->UpUpSectorIndicesPerSum[j][k << 1]<<" "<<this->UpUpSectorIndicesPerSum[j][(k << 1) + 1]<<" "<< this->UpUpSectorIndicesPerSum[j][i << 1] << " "<< this->UpUpSectorIndicesPerSum[j][(i << 1) + 1] << " "; 
                             //cout<<" --> Index="<<Index<<"| "; particles->PrintState(cout,Index);
                             //cout<<" Coeff= "<<Coefficient3*Coefficient/vSource[index]<<endl;  
                          }
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownUpUp
		  if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpUp[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
			    ++TmpInteractionFactor;
			  }	  
		    }
				
		  // now UpDownUpUp		  
		  TmpInteractionFactor = this->InteractionFactorsUpDownUpUp[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
		      ++TmpInteractionFactor;
		    }	  		    
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpUpDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
		      ++TmpInteractionFactor;
		    }
		      
		    // now DownDownUpDown
		    if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		      {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpDown[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			    
			      }
			    ++TmpInteractionFactor;
			  }	  
		      }
				  
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;			
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
		    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  Coefficient3 *= vSource[index];
				
		  // first UpUpDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpDownDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{				      
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;
                          //cout<<"i= "<<index<<"| "; particles->PrintState(cout,index); cout<<"a+a+aa "<< this->UpUpSectorIndicesPerSum[j][k << 1]<<" "<<this->UpUpSectorIndicesPerSum[j][(k << 1) + 1]<<" "<< this->UpUpSectorIndicesPerSum[j][i << 1] << " "<< this->UpUpSectorIndicesPerSum[j][(i << 1) + 1] << " "; 
                          //cout<<" --> Index="<<Index<<"| "; particles->PrintState(cout,Index);
                          //cout<<" Coeff= "<<Coefficient3*Coefficient/vSource[index]<<endl;  			  
			}
		      ++TmpInteractionFactor;
		    }	  

				    
		  // now UpDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpDownDownDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+2][k << 1], this->UpDownSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  vDestination[Index] += Coefficient * (*TmpInteractionFactor) * Coefficient3;		      
			}
		      ++TmpInteractionFactor;
		    }	  
		}
	    }
	}      
    }
}


// core part of the PartialFastMultiplicationMemory method involving two-body term
// 
// particles = pointer to the Hilbert space
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// memory = reference on the amount of memory required for precalculations

inline void ParticleOnSphereTwoLandauLevelDeltaHamiltonian::EvaluateMNTwoBodyFastMultiplicationMemoryComponent(ParticleOnSphereWithSpin* particles, int firstComponent, int lastComponent, long& memory)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;;
  int Dim = particles->GetHilbertSpaceDimension();

  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {  
      for (int idx = firstComponent; idx < lastComponent; ++idx)
	{
	  // Annihilation operators acting on first LL (UpUp)
	  for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAu(idx, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {			
		      // first UpUpUpUp
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
		      // now DownDownUpUp
		      if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
			{
			    for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			      {
				Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
				if (Index < Dim)
				  {
				    ++memory;
				    ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
				  }
			      }	  
			}
				    
		      // now UpDownUpUp
		      if ( j >= 1 && j < (this->NbrUpDownSectorSums + 1) ) 
			{
			    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j-1] ; k++ ) 
			      {
				Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j-1][k << 1], this->UpDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
				if (Index < Dim)
				  {
				    ++memory;
				    ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
				  }
			      }	  
			}
		    }
		}
	    }
		
	  // Annihilation operators acting on both levels (UpDown)
	  for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAd(idx, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {
				    
		      // first UpUpUpDown
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+1] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+1][k << 1], this->UpUpSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
			// now DownDownUpDown
			if ( j >= 1 && j < (this->NbrDownDownSectorSums + 1) ) 
			  {
			    for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-1] ; k++ ) 
			      {
				Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-1][k << 1], this->DownDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
				if (Index < Dim)
				{
				    ++memory;
				    ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
				  }
			      }	  
			  }
				      
			// now UpDownUpDown
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			  }	  
		    }
		}
	    }
	    
		
		
	  // Annihilation operators acting on LLL (DownDown)
	  for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AdAd(idx, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {
		      // first UpUpDownDown
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
		      // now DownDownDownDown
		      for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  

					
		      // now UpDownDownDown	
		      for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+1] ; k++ ) 
			{
			  Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+1][k << 1], this->UpDownSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  
		    }
		}
	    }
	}
    }
  else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      for (int idx = firstComponent; idx < lastComponent; ++idx)
	{
	  // Annihilation operators acting on first LL (UpUp)
	  for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAu(idx, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {			
		      // first UpUpUpUp
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
		      // now DownDownUpUp
		      if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
			{
			    for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			      {
				Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
				if (Index < Dim)
				  {
				    ++memory;
				    ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
				  }
			      }	  
			}
				    
		      // now UpDownUpUp		      
		      for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  			
		    }
		}
	    }
		
	  // Annihilation operators acting on both levels (UpDown)
	  for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AuAd(idx, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {
				    
		      // first UpUpUpDown
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
			// now DownDownUpDown
			if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
			  {
			    for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			      {
				Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
				if (Index < Dim)
				{
				    ++memory;
				    ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
				  }
			      }	  
			  }
				      
			// now UpDownUpDown
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			    {
				++memory;
				++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			  }	  
		    }
		}
	    }
	    
		
		
	  // Annihilation operators acting on LLL (DownDown)
	  for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	    {
	      for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
		{
		  Coefficient3 = particles->AdAd(idx, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
		  if ( Coefficient3 != 0.0 )
		    {
		      // first UpUpDownDown
		      for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
			{
			  Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}
			  
		      // now DownDownDownDown
		      for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
			{
			  Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  

					
		      // now UpDownDownDown	
		      for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+2] ; k++ ) 
			{
			  Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+2][k << 1], this->UpDownSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
			  if (Index < Dim)
			    {
			      ++memory;
			      ++this->NbrInteractionPerComponent[idx - this->PrecalculationShift];
			    }
			}	  
		    }
		}
	    }
	}            
    }
	   
  if ((this->OneBodyInteractionFactorsUpDown != 0) || (this->OneBodyInteractionFactorsUpUp != 0))
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  ++memory;
	  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];	  
	}
    }
  if (this->OneBodyInteractionFactorsUpDown != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsUpDown; ++j)
	    {
	      Index = particles->AduAd(i, this->OneBodyMValuesUpDown[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactorsDownUp != 0)
    {
      for (int i = firstComponent; i < lastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactorsDownUp; ++j)
	    {
	      Index = particles->AddAu(i, this->OneBodyMValuesDownUp[j], Coefficient);
	      if (Index < Dim)
		{
		  ++memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
}


// core part of the FastMultiplication method involving the two-body interaction
// 
// particles = pointer to the Hilbert space
// index = index of the component on which the Hamiltonian has to act on
// indexArray = array where indices connected to the index-th component through the Hamiltonian
// coefficientArray = array of the numerical coefficients related to the indexArray
// position = reference on the current position in arrays indexArray and coefficientArray

inline void ParticleOnSphereTwoLandauLevelDeltaHamiltonian::EvaluateMNTwoBodyFastMultiplicationComponent(ParticleOnSphereWithSpin* particles, int index, int* indexArray, double* coefficientArray, long& position)
{
  int Index;
  double Coefficient = 0.0;
  double Coefficient3 = 0.0;
  double* TmpInteractionFactor;
  int Dim = particles->GetHilbertSpaceDimension();
 
  
  if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic ) 
    {
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{              
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;		     
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownUpUp
		  if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpUp[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
				++position;
			      }
			    ++TmpInteractionFactor;
			  }	  
		    }
				
			  // now UpDownUpUp
		  if ( j >= 1 && j < (this->NbrUpDownSectorSums + 1) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsUpDownUpUp[j-1] + (i * this->NbrUpDownSectorIndicesPerSum[j-1]);	
			for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j-1] ; k++ ) 
			  {
			    Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j-1][k << 1], this->UpDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  			    
				++position;
			      }
			    ++TmpInteractionFactor;
			  }	  
		    }
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  // first UpUpUpDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+1]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+1] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+1][k << 1], this->UpUpSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		    // now DownDownUpDown
		    if ( j >= 1 && j < (this->NbrDownDownSectorSums + 1) ) 
		      {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpDown[j-1] + (i * this->NbrDownDownSectorIndicesPerSum[j-1]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-1] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-1][k << 1], this->DownDownSectorIndicesPerSum[j-1][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
				++position;
			      }
			    ++TmpInteractionFactor;
			  }	  
		      }
				  
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    indexArray[position] = Index;
			    coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			    ++position;
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
	
	    
	    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		// first UpUpDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpDownDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  

				    
		  // now UpDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpDownDownDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j+1]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+1] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+1][k << 1], this->UpDownSectorIndicesPerSum[j+1][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  
		}
	    }
	}
    }
  else if ( this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic ) 
    {
      // Annihilation operators acting on first LL (UpUp)
      for (int j = 0; j < this->NbrUpUpSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpUpSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAu(index, this->UpUpSectorIndicesPerSum[j][i << 1], this->UpUpSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{              
		  // first UpUpUpUp
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpUp[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;		     
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownUpUp
		  if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		    {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpUp[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  
				++position;
			      }
			    ++TmpInteractionFactor;
			  }	  
		    }
				
		  // now UpDownUpUp		  
		  TmpInteractionFactor = this->InteractionFactorsUpDownUpUp[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  			    
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  		    
		}
	    }
	}
	    
      // Annihilation operators acting on both levels (UpDown)
      for (int j = 0; j < this->NbrUpDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrUpDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AuAd(index, this->UpDownSectorIndicesPerSum[j][i << 1], this->UpDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		  // first UpUpUpDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpUpDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j][k << 1], this->UpUpSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		    // now DownDownUpDown
		    if ( j >= 2 && j < (this->NbrDownDownSectorSums + 2) ) 
		      {
			TmpInteractionFactor = this->InteractionFactorsDownDownUpDown[j-2] + (i * this->NbrDownDownSectorIndicesPerSum[j-2]);	
			for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j-2] ; k++ ) 
			  {
			    Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j-2][k << 1], this->DownDownSectorIndicesPerSum[j-2][(k << 1) + 1], Coefficient);
			    if (Index < Dim)
			      {
				indexArray[position] = Index;
				coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;
				++position;
			      }
			    ++TmpInteractionFactor;
			  }	  
		      }
				  
		    // now UpDownUpDown
		    TmpInteractionFactor = this->InteractionFactorsUpDownUpDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j]);	
		    for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j] ; k++ ) 
		      {
			Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j][k << 1], this->UpDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
			if (Index < Dim)
			  {
			    indexArray[position] = Index;
			    coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			    ++position;
			  }
			++TmpInteractionFactor;
		      }	  
		  }
	    }
	}
	
	    
	    
      // Annihilation operators acting on LLL (DownDown)
      for (int j = 0; j < this->NbrDownDownSectorSums; ++j) 
	{
	  for ( int i = 0 ; i < this->NbrDownDownSectorIndicesPerSum[j] ; i++ ) 
	    {
	      Coefficient3 = particles->AdAd(index, this->DownDownSectorIndicesPerSum[j][i << 1], this->DownDownSectorIndicesPerSum[j][(i << 1) + 1]);
	      if ( Coefficient3 != 0.0 )
		{
		// first UpUpDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpUpDownDown[j] + (i * this->NbrUpUpSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpUpSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdu(this->UpUpSectorIndicesPerSum[j+2][k << 1], this->UpUpSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }
		      
		  // now DownDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsDownDownDownDown[j] + (i * this->NbrDownDownSectorIndicesPerSum[j]);	
		  for ( int k = 0 ; k < this->NbrDownDownSectorIndicesPerSum[j] ; k++ ) 
		    {
		      Index = particles->AddAdd(this->DownDownSectorIndicesPerSum[j][k << 1], this->DownDownSectorIndicesPerSum[j][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  

				    
		  // now UpDownDownDown
		  TmpInteractionFactor = this->InteractionFactorsUpDownDownDown[j] + (i * this->NbrUpDownSectorIndicesPerSum[j+2]);	
		  for ( int k = 0 ; k < this->NbrUpDownSectorIndicesPerSum[j+2] ; k++ ) 
		    {
		      Index = particles->AduAdd(this->UpDownSectorIndicesPerSum[j+2][k << 1], this->UpDownSectorIndicesPerSum[j+2][(k << 1) + 1], Coefficient);
		      if (Index < Dim)
			{
			  indexArray[position] = Index;
			  coefficientArray[position] = Coefficient * (*TmpInteractionFactor) * Coefficient3;			  		      
			  ++position;
			}
		      ++TmpInteractionFactor;
		    }	  
		}
	    }
	}      
    }
}


#endif
