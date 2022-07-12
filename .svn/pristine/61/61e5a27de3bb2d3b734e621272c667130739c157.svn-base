////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of FTI band structure calculation operation             //
//                                                                            //
//                        last modification : 26/09/2012                      //
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


#ifndef FTICOMPUTEBANDSTRUCTUREOPERATION_H
#define FTICOMPUTEBANDSTRUCTUREOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"


class FTIComputeBandStructureOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the tight binding model
  AbstractTightBindingModel* TightBindingModel;

 public:
  
  // constructor 
  //
  // tightBindingModel = pointer to the tight binding model
  FTIComputeBandStructureOperation(AbstractTightBindingModel* tightBindingModel);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FTIComputeBandStructureOperation(const FTIComputeBandStructureOperation& operation);
  
  // destructor
  //
  ~FTIComputeBandStructureOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();

  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  virtual int GetHilbertSpaceDimension();

  
 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for simple MPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension  

int FTIComputeBandStructureOperation::GetHilbertSpaceDimension ()
{
  return this->TightBindingModel->GetNbrStatePerBand();
}



#endif
