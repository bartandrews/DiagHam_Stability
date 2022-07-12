////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//               forcing vanishing properties with multiple groups            //
//                        of n-body hard core interaction                     //
//                                                                            //
//                        last modification : 04/06/2011                      //
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


#ifndef PARTICLEONSPHEREMULTIPLEGROUPSNBODYHARDCOREHAMILTONIAN_H
#define PARTICLEONSPHEREMULTIPLEGROUPSNBODYHARDCOREHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"

#include <iostream>


using std::ostream;


class ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian : public ParticleOnSphereNBodyHardCoreHamiltonian
{

 protected:

  // number of groups that have a given vanishing property
  int NbrGroups;
  // number of particle that interact simultaneously through the hard core interaction for the each group
  int* NbrNBodys;
  
 // bool * GroupsNBodyFlags;
  
 public:

  // default constructor
  //
  ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrGroups = number of groups that have a given vanishing property
  // nbrNBodys = number of particle that interact simultaneously through the hard core interaction for the each group
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrGroups, int* nbrNBodys, double l2Factor, 
							AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
							char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();
  
  double EvaluateInteractionCoefficient(int m1, int m2, int m3,int m4, int n1,int n2, int n3, int n4, double ** TmpInteractionCoeffients);
  
  inline double EvaluateSymmetrizedInteractionCoefficient(int m1, int m2, int m3,int m4, int n1,int n2, int n3,int n4, double ** TmpInteractionCoeffients);


};


// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term, symmetrized on the n indices
//
// m1 = first creation operator index
// m2 = second creation operator index
// m3 = third creation operator index
// n1 = first annihilation operator index
// n2 = second annihilation operator index
// n3 = thrid annihilation operator index
// return value = numerical coefficient

inline double ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::EvaluateSymmetrizedInteractionCoefficient(int m1, int m2, int m3,int m4, int n1,int n2, int n3,int n4,double ** TmpInteractionCoeffients)
{  
  double UltimateCoef = 0.0;
		    
		    if ( n1 !=  n2)
		    {
		      if ( n2 !=  n3)
			{
			  if ( n3 !=  n4)
			  {
			    
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n4,  n2,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n2,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n3,  n2,  TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n4,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n4,  n1,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n4,  n3,  n1,  TmpInteractionCoeffients);
			  
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n4,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n1,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n1,  n4,  n2,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n2,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n1,  n2,  TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n3,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n1,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n3,  n2,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n3,  n1,  n2,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n1,  n2,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n1,  n3,  n2,  TmpInteractionCoeffients);
			  
			  }
			  else //  n3 ==  n4
			  {
			    UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n4,  n2,  TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n4,  n1,  TmpInteractionCoeffients);
			  
			  
 //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n4,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n1,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n1,  n4,  n2,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n2,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n1,  n2,  TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  }
			}
			else
			{
			  if ( n3 !=  n4) //( n2 ==  n3) 
			  {
			    
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n2,  n3,  TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n4,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n4,  n1,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n4,  n3,  n1,  TmpInteractionCoeffients);
			  
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n3,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n1,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n1,  n2,  n3,  TmpInteractionCoeffients);
			    
			  }
			  else  //( n3 ==  n4) && ( n2 ==  n3)
			  {
			    UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n1,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n2, n3,  n4,  n1,  TmpInteractionCoeffients);
			  
			  }
			}
		    }
			  else //( n1 ==  n2)
			    {
			      if ( n2 !=  n3)
			      {
				if ( n3 !=  n4)
				{
				  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n4,  n2,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n2,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n3,  n2,  TmpInteractionCoeffients);
			  
			  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n4,  n1,  TmpInteractionCoeffients);
			  
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n1,  n2,  TmpInteractionCoeffients);
			  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 			  
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n3,  n1,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n2,  n1,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n3,  n1,  n2,  TmpInteractionCoeffients);
			      
				}
				else //( n1 ==  n2) && ( n3 ==  n4)
				{
				  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n2,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n3,  n4,  n2,  TmpInteractionCoeffients);
			  
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n1,  n4,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n2,  n4,  n1,  TmpInteractionCoeffients);
			  
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n3, n4,  n1,  n2,  TmpInteractionCoeffients);
			  
				  
				}
			      }
			      else //( n1 ==  n2 ==  n3)
				{
				  if ( n3 !=  n4)
				  {
				    UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
			  
			   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n4,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n4,  n2,  n3,  TmpInteractionCoeffients);
			  
			  UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n4, n1,  n2,  n3,  TmpInteractionCoeffients);
				    
				    
				  }
				else
				{
				   UltimateCoef += this->EvaluateInteractionCoefficient(m1,m2,m3,m4, n1, n2,  n3,  n4,  TmpInteractionCoeffients);
				}
				}
			    }
				return UltimateCoef;
}


#endif
