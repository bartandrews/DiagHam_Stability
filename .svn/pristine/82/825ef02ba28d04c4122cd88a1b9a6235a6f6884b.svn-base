////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of manager for particles on sphere                 //
//                                                                            //
//                        last modification : 05/10/2008                      //
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


#ifndef PARTICLEONSPHEREMANAGER_H
#define PARTICLEONSPHEREMANAGER_H


#include "config.h"


class OptionManager;
class ParticleOnSphere;


class ParticleOnSphereManager
{

 protected:

  // pointer to the option manager
  OptionManager* Options;

  // use fermionic statistics
  bool FermionFlag;
  // use bosonic statistics
  bool BosonFlag;
  
  // indicate which SU(K) internal degree of freedom as to be used
  int SUKIndex;

  // true if some option values can/will be retrieved from a file name
  bool FilenameFlag;

 public:

  // default constructor
  //
  // fermionFlag = use fermionic statistics if true
  // bosonFlag = use bosonic statistics if true, if both are fermionFlag and bosonFlag are equal, add an option to handle statistic selection
  // sUKIndex = indicate which SU(K) internal degree of freedom as to be used
  // filenameFlag = true if some option values can/will be retrieved from a file name
  ParticleOnSphereManager(bool fermionFlag = true, bool bosonFlag = false, int sUKIndex = 1, bool filenameFlag = false);

  // destructor
  //
  ~ParticleOnSphereManager();
 
  // add an option group containing all options related to the Hilbert space construction 
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  // get the Hilbert space defined by the running options and a given Total Lz value
  //
  // totalLz = twice the system total Lz value
  // return value = pointer to the Lanczos algorithm
  ParticleOnSphere* GetHilbertSpace(int totalLz);

 
 protected:

  // get the Hilbert space defined by the running options and a given Total Lz value and for the U(1) case
  //
  // totalLz = twice the system total Lz value
  // return value = pointer to the Lanczos algorithm
  ParticleOnSphere* GetHilbertSpaceU1(int totalLz);

  // get the Hilbert space defined by the running options and a given Total Lz value and for the SU(2) case
  //
  // totalLz = twice the system total Lz value
  // return value = pointer to the Lanczos algorithm
  ParticleOnSphere* GetHilbertSpaceSU2(int totalLz);

  // get the Hilbert space defined by the running options and a given Total Lz value and for the SU(3) case
  //
  // totalLz = twice the system total Lz value
  // return value = pointer to the Lanczos algorithm
  ParticleOnSphere* GetHilbertSpaceSU3(int totalLz);

  // get the Hilbert space defined by the running options and a given Total Lz value and for the SU(4) case
  //
  // totalLz = twice the system total Lz value
  // return value = pointer to the Lanczos algorithm
  ParticleOnSphere* GetHilbertSpaceSU4(int totalLz);


};

#endif
