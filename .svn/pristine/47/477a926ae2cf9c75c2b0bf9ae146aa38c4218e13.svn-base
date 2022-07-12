////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of abstract hamiltonian precalculation operation        //
//                                                                            //
//                        last modification : 13/11/2003                      //
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


#ifndef ABSTRACTPRECALCULATIONOPERATION_H
#define ABSTRACTPRECALCULATIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractPrecalculationOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;
  // index of the first component (if long numbers are required)
  long LargeFirstComponent;
  // number of component  (if long numbers are required)
  long LargeNbrComponent;

 public:
  
  // destructor
  //
  virtual ~AbstractPrecalculationOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  virtual int GetHilbertSpaceDimension () = 0;

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // mpiNodeNbr = node on which architecture is run
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // mpiNodeNbr = provide the additional MPI node ID
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr);

  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

#endif
