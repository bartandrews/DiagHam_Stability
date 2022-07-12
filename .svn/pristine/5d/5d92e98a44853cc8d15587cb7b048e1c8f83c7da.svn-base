#ifndef _COMPLETELYPOSITIVEMAP_H
#define _COMPLETELYPOSITIVEMAP_H 

#include "Architecture/AbstractArchitecture.h"
#include "Matrix/ComplexMatrix.h"


class CompletelyPositiveMap 
{

  unsigned int KrausRank;
  
  ComplexMatrix * DirectOperators;
  ComplexMatrix * DaggersOperators;
  
 public:  
  CompletelyPositiveMap (unsigned int krausRank, ComplexMatrix * directOperators, ComplexMatrix * daggersOperators); 


  void ApplyTheMap(ComplexMatrix & sourceMatrices, ComplexMatrix & destinationMatrices);
  inline unsigned int GetDimension() const { return this->DirectOperators[0].GetNbrRow();}
  
};

#endif
