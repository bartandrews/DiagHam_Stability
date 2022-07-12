////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class of particle on a torus with spin                  //
//                                                                            //
//                        last modification : 10/09/2002                      //
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


#ifndef PARTICLEONTORUSWITHSPIN_H
#define PARTICLEONTORUSWITHSPIN_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"

#include <iostream>


using std::ostream;


class ParticleOnTorusWithSpin :  public AbstractHilbertSpace
{

 public:

  enum 
    {
      BosonicStatistic = 0x1,
      FermionicStatistic = 0x2,
    };
  
  // virtual destructor
  //
  virtual ~ParticleOnTorusWithSpin ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // apply sum_m au^+_m au_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int SumAudAu (int index, double& coefficient) = 0;

  // apply sum_m ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int SumAddAd (int index, double& coefficient) = 0;

  // apply au^+_m au_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAu (int index, int m, double& coefficient) = 0;

  // apply ad^+_m ad_m operator to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for density operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, double& coefficient) = 0;

  // apply au^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply ad^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply ad^+_m1 ad^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply au^+_m1 au^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAudAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply au^+_m1 au^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAudAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply ad^+_m1 au^+_m2 au_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAudAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply ad^+_m1 ad^+_m2 ad_n1 au_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply au^+_m1 ad^+_m2 ad_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // apply au^+_m1 ad^+_m2 au_n1 ad_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AudAddAuAd (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

};

#endif


