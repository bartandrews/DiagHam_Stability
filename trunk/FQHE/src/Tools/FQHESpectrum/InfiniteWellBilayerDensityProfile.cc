#include "InfiniteWellBilayerDensityProfile.h"

#include <cmath>

using std::sqrt;
using std::sin;

// the standard constructor -> well with 1 magnetic length width
InfiniteWellBilayerDensityProfile::InfiniteWellBilayerDensityProfile()
{
  this->Width=1.0;
  this->SqNorm=1.0;
  this->Symmetry=1.0;
}

// constructor
// width = width of the well
// band = band index (>0 : right well, <0 : left well)
// 
InfiniteWellBilayerDensityProfile::InfiniteWellBilayerDensityProfile(double width, int band)
{
  this->Width=width;
  this->SqNorm=1.0/width;
  this->Symmetry=(band>=0?1.0:-1.0);
}
  
// virtual destructor
InfiniteWellBilayerDensityProfile::~InfiniteWellBilayerDensityProfile()
{
}

// get minimum and maximum value of density profile where the probability density is larger than precision
// min = minimum value of z offset
// max = maximum value of z offset
// precision = requested precision
void InfiniteWellBilayerDensityProfile::GetSupport(double &min, double &max, double precision)
{
  min=-this->Width/2.0;
  max=this->Width/2.0;
}

// evaluate the density for a given offset
// z = offset of distribution
double InfiniteWellBilayerDensityProfile::GetValue(double z)
{
  z+=Width/2.0;
  double tmp=M_PI*z/this->Width;
  double cz=sin(tmp)+Symmetry*sin(2.0*tmp);
  return this->SqNorm*cz*cz;
}
