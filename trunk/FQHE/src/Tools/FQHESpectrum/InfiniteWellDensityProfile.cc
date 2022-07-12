#include "InfiniteWellDensityProfile.h"

#include <cmath>

using std::sqrt;
using std::sin;

// the standard constructor -> well with 1 magnetic length width
InfiniteWellDensityProfile::InfiniteWellDensityProfile()
{
  this->Width=1.0;
  this->SqNorm=2.0;
  this->Band=1;
}

// constructor
// width = width of the well
// band = band index (allows for excited states to be obtained)
// 
InfiniteWellDensityProfile::InfiniteWellDensityProfile(double width, int band)
{
  this->Width=width;
  this->SqNorm=2.0/width;
  this->Band=band+1;
}
  
// virtual destructor
InfiniteWellDensityProfile::~InfiniteWellDensityProfile()
{
}

// get minimum and maximum value of density profile where the probability density is larger than precision
// min = minimum value of z offset
// max = maximum value of z offset
// precision = requested precision
void InfiniteWellDensityProfile::GetSupport(double &min, double &max, double precision)
{
  min=-this->Width/2.0;
  max=this->Width/2.0;
}

// evaluate the density for a given offset
// z = offset of distribution
double InfiniteWellDensityProfile::GetValue(double z)
{
  z+=Width/2.0;
  double cz=sin(Band*M_PI*z/this->Width);
  return this->SqNorm*cz*cz;
}
