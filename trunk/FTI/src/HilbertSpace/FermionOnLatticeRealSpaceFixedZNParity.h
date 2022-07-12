////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                 where only the particle number is fixed modulo N           //
//                                                                            //
//                        last modification : 22/12/2017                      //
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


#ifndef FERMIONONLATTICEREALSPACEFIXEDZNPARITY_H
#define FERMIONONLATTICEREALSPACEFIXEDZNPARITY_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class FermionOnLatticeRealSpaceFixedZNParity : public FermionOnLatticeRealSpace
{


  friend class FermionOnLatticeRealSpaceAnd2DTranslation;
  
 protected:
  
  // N of the particle number modulo N
  int ZNValue;
  
  // "parity" of the particle number (between 0 and ZNValue - 1)
  int ParitySector;
  
  // total number of sites
  int NbrSite;
  
 public:

  // default constructor
  // 
  FermionOnLatticeRealSpaceFixedZNParity ();

  // basic constructor
  // 
  // nbrSite = number of sites
  // zNValue = N of the particle number modulo N
  // paritySector = "parity" of the particle number (between 0 and ZNValue - 1)
  FermionOnLatticeRealSpaceFixedZNParity (int nbrSite, int zNValue, int paritySector);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnLatticeRealSpaceFixedZNParity(const FermionOnLatticeRealSpaceFixedZNParity& fermions);

  // destructor
  //
  ~FermionOnLatticeRealSpaceFixedZNParity ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnLatticeRealSpaceFixedZNParity& operator = (const FermionOnLatticeRealSpaceFixedZNParity& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();


 protected:

};


#endif


