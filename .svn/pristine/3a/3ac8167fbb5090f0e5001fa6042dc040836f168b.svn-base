#include "config.h"
#include "TabulatedDensityProfile.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <cmath>

#include <iostream>
using std::cout;
using std::endl;

using std::exp;
using std::log;


// the standard constructor
TabulatedDensityProfile::TabulatedDensityProfile()
{
  this->NbrDataPoints=0;
}


// the norm of the (square) wavefunction
double SqNorm;


// constructor
// dataFile = file name of text-file to be opened
// scale = scale-factor in the z-direction
// squareValues = option to square the values of individual points (wavefunction provided)
// splineInterpolate = flag indicating whether a spline interpolation should be used  
// checkNorm = flag indicating whether the norm should be recalculated
TabulatedDensityProfile::TabulatedDensityProfile(char *dataFile, double scale, bool squareValues, bool checkNorm)
{
  MultiColumnASCIIFile DataParser('\t');

  if (!(DataParser.Parse(dataFile)))
    {
      cout << "Problem parsing density profile "<<dataFile<<endl;
      exit(1);
    }

  if (DataParser.GetNbrColumns()!=2)
    {
      cout << "error: TabulatedDensityProfile requires two columns [z, rho(z)]"<<endl;
      exit(1);
    }

  this->NbrDataPoints = DataParser.GetNbrLines();
  ZValues = DataParser.GetAsDoubleArray(0);
  Density = DataParser.GetAsDoubleArray(1);

  if ((scale!=1.0)&&(scale!=0.0))
    for (int i=0; i<NbrDataPoints; ++i)      
      {
	ZValues[i]*=scale;
	Density[i]/=scale;  
      }

  this->Zmin=ZValues[0];
  this->Zmax=ZValues[NbrDataPoints-1];

  if (Zmin > Zmax)
    {
      cout << "TabulatedDensityProfile requires a density profile with increasing z-coordinate in the first column"<<endl;
      exit(1);
    }

  
  if (squareValues)
    for (int i=0; i<NbrDataPoints; ++i) Density[i]*=Density[i];  

#ifdef HAVE_GSL
  this->SplineInterpolate=true;

  // alternative choices
  const gsl_interp_type *SplineType = gsl_interp_cspline;
  //const gsl_interp_type *t = gsl_interp_akima;  

  this->Spline = gsl_spline_alloc (SplineType, NbrDataPoints);
  this->SplineAccelerator = gsl_interp_accel_alloc ();
  gsl_spline_init (Spline, ZValues, Density, NbrDataPoints);

#else
  this->SplineInterpolate=false;
#endif
}
  
// virtual destructor
TabulatedDensityProfile::~TabulatedDensityProfile()
{
  if (this->NbrDataPoints!=0)
    {
      delete [] this->ZValues;
      delete [] this->Density;
    }

#ifdef HAVE_GSL
  if (this->SplineInterpolate)
    {
      gsl_spline_free (this->Spline);
      gsl_interp_accel_free (this->SplineAccelerator);
    }
#endif
}

// get minimum and maximum value of density profile where the probability density is larger than precision
// min = minimum value of z offset
// max = maximum value of z offset
// precision = requested precision
void TabulatedDensityProfile::GetSupport(double &min, double &max, double precision)
{
  min=this->Zmin;
  max=this->Zmax;

  int PosMin=0;
  while ((fabs(this->Density[PosMin])<precision)&&(PosMin<NbrDataPoints-1)
	 && (!((Density[PosMin]==0.0)&&(Density[PosMin+1]!=0.0))) )
    ++PosMin;

  int PosMax=NbrDataPoints-1;
  while ((fabs(this->Density[PosMax])<precision)&&(PosMax>1)
	 && (!((Density[PosMax]==0.0)&&(Density[PosMax-1]!=0.0))) )
    --PosMax;
  if (PosMin>=PosMax)
    {
      cout << "Attention: Non-sensical bounds in TabulatedDensityProfile::GetSupport!"<<endl;
    }
  min = ZValues[PosMin];
  max = ZValues[PosMax];
}

// evaluate the density for a given offset
// z = offset of distribution
double TabulatedDensityProfile::GetValue(double z)
{
#ifdef HAVE_GSL
  return gsl_spline_eval (this->Spline, z, this->SplineAccelerator);
#else
  cout << "Attention: Native implementation of splines not provided, please link to GSL!"<<endl;
  return 0.0;
#endif
}

