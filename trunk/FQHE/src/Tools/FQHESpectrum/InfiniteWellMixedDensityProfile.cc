#include "InfiniteWellMixedDensityProfile.h"

#include <cmath>

using std::sqrt;
using std::sin;

// the standard constructor -> well with 1 magnetic length width
InfiniteWellMixedDensityProfile::InfiniteWellMixedDensityProfile()
{
  this->Width=1.0;
  this->SqNorm=2.0;
  this->Band=1;
}

// constructor
// width = width of the well
// band = band index (allows for excited states to be obtained)
// 
InfiniteWellMixedDensityProfile::InfiniteWellMixedDensityProfile(double width, int band)
{
  this->Width=width;
  this->SqNorm=2.0/width;
  this->Band=band+1;
}
  
// virtual destructor
InfiniteWellMixedDensityProfile::~InfiniteWellMixedDensityProfile()
{
}

// get minimum and maximum value of density profile where the probability density is larger than precision
// min = minimum value of z offset
// max = maximum value of z offset
// precision = requested precision
void InfiniteWellMixedDensityProfile::GetSupport(double &min, double &max, double precision)
{
  min=-this->Width/2.0;
  max=this->Width/2.0;
}

// evaluate the density for a given offset
// z = offset of distribution
double InfiniteWellMixedDensityProfile::GetValue(double z)
{
  z+=Width/2.0;
  double cz1=sin(Band*M_PI*z/this->Width);
  double cz2=sin((Band+1)*M_PI*z/this->Width);
  return this->SqNorm*cz1*cz2;
}
