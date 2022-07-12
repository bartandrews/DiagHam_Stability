////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of Lanczos Manager                         //
//                                                                            //
//                        last modification : 01/05/2008                      //
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


#ifndef LANCZOSMANAGER_H
#define LANCZOSMANAGER_H


#include "config.h"


class AbstractLanczosAlgorithm;
class AbstractArchitecture;
class OptionManager;


class LanczosManager
{

 protected:

  // pointer to the best avalaible Lancozs algorithm in agrement with the option constraints
  AbstractLanczosAlgorithm* LanczosAlgorithm;

  // pointer to the option manager
  OptionManager* Options;

  // use complex Lanczos algorithms if true
  bool ComplexFlag;

 public:

  // default constructor
  //
  // complexFlag = use complex Lanczos algorithms if true
  LanczosManager(bool complexFlag = false);

  // destructor
  //
  ~LanczosManager();
 
  // add an option group containing all options related to the architecture 
  //
  // manager = pointer to the option manager
  void AddOptionGroup(OptionManager* manager);

  // get the best avalaible Lanczos algorithm in agrement with the option constraints
  //
  // architecture = pointer to the architecture to use within the Lanczos algorithm
  // forceEigenstateComputation = if true, force computation of eigenstates
  // useLapack = true if some Lanczos algorithms should rely on the Lapack library
  // return value = pointer to the Lanczos algorithm
  AbstractLanczosAlgorithm* GetLanczosAlgorithm(AbstractArchitecture* architecture, bool forceEigenstateComputation = false, bool useLapack = false);

  // delete last created Lanczos object
  //
  // return = true if object deleted
  bool FreeLanczosAlgorithm();

  // set Lanczos to complex algorithms
  //
  void SetComplexAlgorithms();

  // set Lanczos to real algorithms
  //
  void SetRealAlgorithms();

};


#endif

