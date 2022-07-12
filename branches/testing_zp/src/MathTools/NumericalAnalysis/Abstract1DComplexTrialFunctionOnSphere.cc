
#include "Abstract1DComplexTrialFunctionOnSphere.h"

// virtual destructor
//
Abstract1DComplexTrialFunctionOnSphere::~Abstract1DComplexTrialFunctionOnSphere()
{
}

// get function properties, and possible extensions of interface 
// 
unsigned Abstract1DComplexTrialFunctionOnSphere::GetProperties()
{
  return Abstract1DComplexFunction::Basic | Abstract1DComplexFunction::OnSphere | Abstract1DComplexFunction::Trial;
}
