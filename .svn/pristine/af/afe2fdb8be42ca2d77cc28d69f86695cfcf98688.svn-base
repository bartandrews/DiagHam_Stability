#ifndef _REALMPOPeratorDefinedByFiles_H
#define _REALMPOPeratorDefinedByFiles_H

#include "MPSObjects/RealMPOperatorOBC.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "Architecture/AbstractArchitecture.h"

class RealMPOPeratorDefinedByFiles : public RealMPOperatorOBC
{
 protected:

  
 public:

  RealMPOPeratorDefinedByFiles();
  RealMPOPeratorDefinedByFiles( MultiColumnASCIIFile & tensorElementsFile, MultiColumnASCIIFile & boundaryVectorsFile,AbstractArchitecture * architecture = 0); 
  ~RealMPOPeratorDefinedByFiles();
  
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
  virtual void InitializeBoundaryVectors(MultiColumnASCIIFile & boundaryVectorsFile);

};

#endif
