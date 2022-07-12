////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of internal interaction for                    //
//                 particles on torus with coulombian interaction             //
//                                                                            //
//                        last modification : 08/11/2002                      //
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


#ifndef PARTICLEONTORUSCOULOMBINTERNALINTERACTION_H
#define PARTICLEONTORUSCOULOMBINTERNALINTERACTION_H


#include "config.h"
#include "Interaction/InternalInteraction/AbstractInternalInteraction.h"
#include "GeneralTools/List.h"


class Matrix;


class ParticleOnTorusCoulombInternalInteraction: public AbstractInternalInteraction
{

 private:

  // flag to indicate if particles are fermions
  bool FermionFlag;

  // precision on interaction coefficient
  double Precision;

  // ratio between the width in the x direction and the width in the y direction
  double Ratio;
  // ratio between the width in the y direction and the width in the x direction
  double InvRatio;

 public:

  // constructor
  //
  // fermion = flag to indicate if particles are fermions
  // nbrSite = total number of sites
  // ratio = ratio between the width in the x direction and the width in the y direction
  // precision = precision on interaction coefficient 
  ParticleOnTorusCoulombInternalInteraction(bool fermion, int nbrSite, double ratio, double precision = MACHINE_PRECISION);

  // destructor
  //
  ~ParticleOnTorusCoulombInternalInteraction();

  // evaluate and add to hamitonian internal interaction terms
  //
  // hamiltonian = reference on hamiltonian matrix representation
  // operators = list of operators used to construct interaction terms
  // return value = reference on hamiltonian matrix representation
  Matrix& AddInteraction (Matrix& hamiltonian, List<Matrix*>& operators);

 private:

  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  double EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

};

#endif
