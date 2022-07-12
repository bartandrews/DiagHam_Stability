#ifndef _COMPLEXMPOPeratorDefinedByFiles_H
#define _COMPLEXMPOPeratorDefinedByFiles_H

#include "MPSObjects/ComplexMPOperatorOBC.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class ComplexMPOPeratorDefinedByFiles : public ComplexMPOperatorOBC
{
 protected:

  
 public:

  ComplexMPOPeratorDefinedByFiles();
  ComplexMPOPeratorDefinedByFiles(MultiColumnASCIIFile & tensorElementsFile, MultiColumnASCIIFile & boundaryVectorsFile, AbstractArchitecture * architecture  = 0); 
  ~ComplexMPOPeratorDefinedByFiles();
  
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
  virtual void InitializeBoundaryVectors(MultiColumnASCIIFile & boundaryVectorsFile);
 
};

#endif
