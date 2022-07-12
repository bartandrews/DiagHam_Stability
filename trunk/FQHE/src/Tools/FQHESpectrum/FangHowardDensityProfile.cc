#include "config.h"
#include "FangHowardDensityProfile.h"
#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

using std::exp;
using std::log;


// the standard constructor -> well with 1 magnetic length width
FangHowardDensityProfile::FangHowardDensityProfile()
{
  this->Width=1.0;
  this->SqNorm=13.5;
  this->ExpIndex=-3.0;
  this->Shift=-1.0;
  this->Flip=false;
}

// constructor
// width = width of the well
// shift = flag indicating whether to shift the potential such that it's centre of mass is at zero
// flip = flip direction of z-axis
//
FangHowardDensityProfile::FangHowardDensityProfile(double width, bool shift, bool flip)
{
  this->Width=width;
  this->SqNorm=13.5/(width*width*width);
  this->ExpIndex=-3.0/width;
  this->Flip=flip;
  if (shift)
    this->Shift=-width;
  else
    this->Shift=0.0;
  if (flip)
    {
      this->ExpIndex *= -1.0;
      this->Shift *= -1.0;
    }
}
  
// virtual destructor
FangHowardDensityProfile::~FangHowardDensityProfile()
{
}

// get minimum and maximum value of density profile where the probability density is larger than precision
// min = minimum value of z offset
// max = maximum value of z offset
// precision = requested precision
void FangHowardDensityProfile::GetSupport(double &min, double &max, double precision)
{
  if (this->Flip)
    {
      max=-Shift;
      min=this->Width/3.0*log(precision/SqNorm);
      min-=2.0*this->Width/3.0*log(max)-Shift;
      // cout << "Test: "<<this->GetValue(max)<< endl;
    }
  else
    {
      min=-Shift;
      max=-this->Width/3.0*log(precision/SqNorm);
      max+=2.0*this->Width/3.0*log(max)-Shift;
      // cout << "Test: "<<this->GetValue(max)<< endl;
    }
}

// evaluate the density for a given offset
// z = offset of distribution
double FangHowardDensityProfile::GetValue(double z)
{
  z-=Shift;
  return this->SqNorm*z*z*exp(this->ExpIndex*z);
}

