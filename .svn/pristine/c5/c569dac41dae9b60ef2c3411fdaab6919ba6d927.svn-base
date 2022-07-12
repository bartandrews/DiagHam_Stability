#include "CompletelyPositiveMap.h"

CompletelyPositiveMap::CompletelyPositiveMap (unsigned int krausRank,  ComplexMatrix * directOperators, ComplexMatrix * daggersOperators) : KrausRank (krausRank) ,  DirectOperators (directOperators), DaggersOperators(daggersOperators)
{
  for(int i = 0 ; i < KrausRank; i++)
    this->DaggersOperators[i].HermitianTranspose();
}

void CompletelyPositiveMap::ApplyTheMap(ComplexMatrix & sourceMatrices, ComplexMatrix & destinationMatrices)
{
  destinationMatrices = this->DaggersOperators[0] * sourceMatrices * this->DirectOperators[0];
  for (int i = 1 ; i < this->KrausRank; i++)
    {
      destinationMatrices = this->DaggersOperators[i] * sourceMatrices * this->DirectOperators[i]; 
    }
  return; 
}
