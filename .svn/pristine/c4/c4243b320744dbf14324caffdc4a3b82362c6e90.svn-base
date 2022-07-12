#include "PseudoPotentials.h"

#include "Vector/RealVector.h"

#include "MathTools/Complex.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "PseudoPotentialGenerator.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include <iostream>

using std::cout;
using std::endl;

// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double* EvaluatePseudopotentials(int nbrFlux, int landauLevel, double layerSeparation, bool quiet)
{
  cout.precision(14);
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* Pseudopotentials = new double [MaxMomentum + 1];
  ClebschGordanCoefficients MainCoefficients(MaxMomentum, MaxMomentum);
  ClebschGordanCoefficients* Coefficients = new ClebschGordanCoefficients[MaxMomentum + 1];
  // new formfactors for finite thickness/separation:
  double *FormFactors = new double[MaxMomentum + 1];
  double dd = layerSeparation*layerSeparation;
  double base = ( sqrt(2.0*nbrFlux + dd) - layerSeparation ) / ( sqrt(2.0*nbrFlux + dd) + layerSeparation );
  FormFactors[0]=sqrt(base);
  for (int l = 1; l <= MaxMomentum; ++l)
    FormFactors[l] = FormFactors[l-1]*base;
  for (int l = 0; l <= MaxMomentum; ++l)
    Coefficients[l] = ClebschGordanCoefficients(MaxMomentum, l << 1);  
  for (int l = 0; l <= MaxMomentum; ++l)
    {
      double TmpPseudopotentials = 0.0;
      for (int m1 = -MaxMomentum; m1 <= MaxMomentum; m1 +=2)
	for (int m2 = -MaxMomentum; m2 <= MaxMomentum; m2 +=2)
	  {
	    double TmpCoef = 0.0;
	    double TmpCoef2;
	    int Min = abs(m1 - m2) >> 1;
	    double Sign = 1.0;
	    for (int j = Min; j <= MaxMomentum; ++j)
	      {
		TmpCoef2 = Coefficients[j].GetCoefficient(m1, m2 - m1, MaxMomentum) * Coefficients[j].GetCoefficient(nbrFlux, 0, MaxMomentum);
		TmpCoef += FormFactors[j] * Sign * TmpCoef2 * TmpCoef2;
		Sign *= -1.0;
	      }
	    TmpPseudopotentials += (MainCoefficients.GetCoefficient(m1, -m1, l << 1) * 
				    MainCoefficients.GetCoefficient(m2, -m2, l << 1) * TmpCoef);
	  }
      Pseudopotentials[MaxMomentum - l] = TmpPseudopotentials / sqrt (0.5 * nbrFlux);
      if (quiet == false)
	cout << "V[" << (MaxMomentum - l) << "] = " << Pseudopotentials[MaxMomentum - l] << endl;
    }
  delete[] Coefficients;
  delete[] FormFactors;
  return Pseudopotentials;
}

// evalute one body potentials for two impurities located at the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPolePotential = potential of the impurity located at the north pole
// southPolePotential = potential of the impurity located at the south pole
// return value = array that conatins the pseudopotentials

double* EvaluateOneBodyPotentials(int nbrFlux, int landauLevel, double northPolePotential, double southPolePotential)
{
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* OneBodyPotentials = new double [MaxMomentum + 1];
  for (int i = 0; i <= MaxMomentum; ++i)
   OneBodyPotentials[i] = 0.0;
  ParticleOnSphereGenericLLFunctionBasis Basis (nbrFlux, landauLevel);
  RealVector Value(2, true);
  Complex TmpValue;
  Value[0] = M_PI;
  Basis.GetFunctionValue(Value, TmpValue, landauLevel);
  OneBodyPotentials[landauLevel] = southPolePotential * SqrNorm(TmpValue);
  Value[0] = 0.0;
  Basis.GetFunctionValue(Value, TmpValue, MaxMomentum - landauLevel);
  OneBodyPotentials[MaxMomentum - landauLevel] = northPolePotential * SqrNorm(TmpValue);  
  return OneBodyPotentials;
}


// evaluate pseudopotentials coefficients of the monomials r^n in units of 1/R
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle of the corresponding LLL problem)
// exponentN = exponent of the monomial
// onlyOdd = boolean indidicating whether it's sufficient to reproduce only the odd pseudopotentials V_(2m+1)
// return value = array that conatins the coefficients V_m(r^n)
// where m runs over 0,...,nbrFlux, or if option onlyOdd given, from 0 to nbrFlux/2 with entries V_(2m+1)(r^n)
//
double* GetMonomialPseudopotentials(int nbrFlux, int exponentN, bool onlyOdd, bool verbose)
{
  FactorialCoefficient TmpCoeff;
  int sizeRst = (nbrFlux+1);
  if (onlyOdd) // fit only the odd pseudopotentials?
    {
      sizeRst /= 2;
      if (nbrFlux&1==0) --sizeRst;
    }    
  int spacing = (onlyOdd ?2:1);
  int M= (onlyOdd?1:0);
  double *rst = new double[sizeRst];  
  if ( exponentN%2 == 0)  // even exponents:
    {      
      TmpCoeff.Power2Multiply(exponentN);
      TmpCoeff.FactorialDivide(nbrFlux+exponentN/2 +1);
      TmpCoeff.FactorialMultiply(nbrFlux+ 1);
      TmpCoeff.FactorialDivide(nbrFlux+exponentN/2 +1);
      TmpCoeff.FactorialMultiply(nbrFlux+ 1);      
      double Prefactor =  TmpCoeff.GetNumericalValue();
      if (verbose) cout << "Pseudopotential Coefficients of <r^" << exponentN<<">"<<endl;
      for (int i=0; i<sizeRst; ++i, M+=spacing)
	{
	  int J=nbrFlux-M;	  
	  TmpCoeff.SetToOne();
	  TmpCoeff.FactorialMultiply(nbrFlux+exponentN/2-J);
	  TmpCoeff.FactorialDivide(nbrFlux-J);
	  TmpCoeff.FactorialMultiply(nbrFlux+exponentN/2+J+1);
	  TmpCoeff.FactorialDivide(nbrFlux+J+1);	  
	  rst[i]= Prefactor * TmpCoeff.GetNumericalValue() / sqrt (0.5 * nbrFlux);
	  if (verbose) cout << "V_"<<M<<"="<<rst[i]<<endl;
	}      
    }
  else // odd exponents:
    {
      int m=(exponentN+1)/2; // N = 2m-1
      if (verbose) cout << "Pseudopotential Coefficients of <r^" << exponentN<<">"<<endl;
      for (int i=0; i<sizeRst; ++i, M+=spacing)
	{
	  int J=nbrFlux-M;	  
	  TmpCoeff.SetToOne();
	  TmpCoeff.Power2Multiply(2*m+1);
	  
	  TmpCoeff.FactorialMultiply(2*nbrFlux+2*m+2*J+2);
	  TmpCoeff.FactorialDivide(nbrFlux+J+1);
	  TmpCoeff.FactorialDivide(nbrFlux+m+J+1);
	  
	  TmpCoeff.FactorialMultiply(2*nbrFlux+2*m-2*J);
	  TmpCoeff.FactorialDivide(nbrFlux-J);
	  TmpCoeff.FactorialDivide(nbrFlux+m-J);
	  
	  TmpCoeff.FactorialMultiply(nbrFlux+1);
	  TmpCoeff.FactorialMultiply(nbrFlux+m+1);
	  TmpCoeff.FactorialDivide(2*nbrFlux+2*m+2);

	  TmpCoeff.FactorialMultiply(nbrFlux+1);
	  TmpCoeff.FactorialMultiply(nbrFlux+m+1);
	  TmpCoeff.FactorialDivide(2*nbrFlux+2*m+2);
	  
	  rst[i]= TmpCoeff.GetNumericalValue() / sqrt (0.5 * nbrFlux);
	  if (verbose) cout << "V_"<<M<<"="<<rst[i]<<endl;
	}
    }
  return rst;
}



#ifdef HAVE_GSL  

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

namespace DiagPseudoPotentials
{

  struct paramZ2
  {
    double z1;
    double d;
    AbstractZDensityProfile *Rho2;
    PseudoPotentialGenerator *VmGen;
    gsl_spline *Vm;
    gsl_interp_accel *VmAcc;
    int relativeM;
  };

  struct paramZ1
  {
    double Z2min;
    double Z2max;
    double epsAbs;
    double epsRel;
    paramZ2 *innerVar;
    AbstractZDensityProfile *Rho1;
    gsl_function *innerIntegral;
    size_t sizeW;
    gsl_integration_workspace *w;
  };

  // inner integration of density profile
  double IntegrandZ2 (double z2, void * params)
  {
    paramZ2 *parameters=(paramZ2*)params;
    double result = gsl_spline_eval (parameters->Vm, fabs(parameters->z1-z2-parameters->d), parameters->VmAcc);
    result *=parameters->Rho2->GetValue(z2);    
    return result;
  }


  // inner integration of density profile without interpolation
  double IntegrandZ2NoInterpolation (double z2, void * params)
  {
    paramZ2 *parameters=(paramZ2*)params;
    double result = parameters->VmGen->Evaluate(parameters->relativeM, fabs(parameters->z1-z2-parameters->d));
    result *=parameters->Rho2->GetValue(z2);    
    return result;
  }

  
  // outer integration of density profile
  double IntegrandZ1 (double z1, void * params)
  {
    paramZ1 *parameters=(paramZ1*)params;
    parameters->innerVar->z1=z1;
    double result, error;
    // call inner integral over z2, maybe change algorithm by choosing other KEY values
    gsl_integration_qag (parameters->innerIntegral, parameters->Z2min, parameters->Z2max,
			 parameters->epsAbs, parameters->epsRel, parameters->sizeW,
			 /* KEY */ GSL_INTEG_GAUSS41, parameters->w, &result, &error);        
    return result*parameters->Rho1->GetValue(z1);
  }
  
}

#endif

// evalute pseudopotentials for coulomb interaction in a given Landau level with a given density profile
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// zDensity = density distribution of the layer 
// points = number of points where exact pseudopotentials are calculated
// multiplier = number of integration intervals used per point of discretization
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// zDensity2 = (optional) density distribution of layer 2, if absent, taken to be equal to 1st profile
// return value = array that contains the pseudopotentials

double* EvaluateFiniteWidthPseudoPotentials(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity, int points, double multiplier, double layerSeparation, AbstractZDensityProfile *zDensity2)
{
#ifdef HAVE_GSL
  gsl_spline *VmSpline=0;
  gsl_interp_accel *VmSplineAcc=0;
  double Z1min, Z1max, Width1, Z2min, Z2max, Width2;
  zDensity->GetSupport(Z1min, Z1max);
  if (zDensity2==0) zDensity2 = zDensity;
  zDensity2->GetSupport(Z2min, Z2max);
  Width1 = Z1max - Z1min;
  Width2 = Z2max - Z2min;

  double TrueMin = 0.0;
  if (layerSeparation>(Width1+Width2)/2.0)
    TrueMin = layerSeparation-(Width1+Width2)/2.0;
  double TrueMax = layerSeparation + (Width1+Width2)/2.0;
  
  
  // alternative choices
  const gsl_interp_type *t = gsl_interp_cspline;
  //const gsl_interp_type *t = gsl_interp_akima;  
  VmSpline = gsl_spline_alloc (t, points);
  VmSplineAcc = gsl_interp_accel_alloc ();  
  double * Distances = new double[points];
  double **PseudopotentialValues = new double*[points];

  for (int i = 0; i < points; ++i)
     {
       Distances[i]=TrueMin+(double)i*(TrueMax-TrueMin)/(points-1);
       PseudopotentialValues[i] = EvaluatePseudopotentials(nbrFlux, landauLevel, Distances[i], /* quiet */ true);
     }

  // number of PP coefficients:
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* TmpPseudopotentials = new double [points];
  double* FinalPseudopotentials = new double [points];

  // integration workspaces
  gsl_integration_workspace * IntW1
    = gsl_integration_workspace_alloc ((int)(multiplier*points));
  gsl_integration_workspace * IntW2
    = gsl_integration_workspace_alloc ((int)(multiplier*points));


  DiagPseudoPotentials::paramZ2 InnerParameters;
  
  InnerParameters.d = layerSeparation;
  InnerParameters.Rho2 = zDensity2;  
  InnerParameters.VmAcc = VmSplineAcc;
  InnerParameters.VmGen = NULL;

  DiagPseudoPotentials::paramZ1 OuterParameters;
  OuterParameters.Z2min = Z2min ;
  OuterParameters.Z2max = Z2max;
  OuterParameters.Rho1 = zDensity;  
  OuterParameters.epsAbs = 0;
  OuterParameters.epsRel = 1e-8;
  OuterParameters.innerVar = &InnerParameters;
  gsl_function TheInnerIntegrand;
  TheInnerIntegrand.function = &DiagPseudoPotentials::IntegrandZ2;
  TheInnerIntegrand.params = &InnerParameters;
  OuterParameters.innerIntegral = &TheInnerIntegrand;
  OuterParameters.sizeW = (int)(multiplier*points);
  OuterParameters.w=IntW2;

  gsl_function TheOuterIntegrand;
  TheOuterIntegrand.function = &DiagPseudoPotentials::IntegrandZ1;
  TheOuterIntegrand.params = &OuterParameters;

  double result, error;
  
  for (int l=0; l<=MaxMomentum; ++l)
    {
      // create interpolation of pseudopotential
      for (int p=0; p<points; ++p)
	TmpPseudopotentials[p]=PseudopotentialValues[p][l];      
      gsl_spline_init (VmSpline, Distances, TmpPseudopotentials, points);
      InnerParameters.Vm=VmSpline;
      
      // call the outer integration over z1
      gsl_integration_qag (&TheOuterIntegrand, Z1min, Z1max, OuterParameters.epsAbs, OuterParameters.epsRel,
			   OuterParameters.sizeW, /* KEY */ GSL_INTEG_GAUSS41, IntW1, &result, &error);        
      
      FinalPseudopotentials[l]=result;
      cout << "V_"<<l<<"="<<result<<" +/- "<< error<<endl;
    }

  // clean up
  
  gsl_integration_workspace_free (IntW1);
  gsl_integration_workspace_free (IntW2);
  
  gsl_spline_free (VmSpline);
  gsl_interp_accel_free (VmSplineAcc);

  for (int i = 0; i < points; ++i)
    delete [] PseudopotentialValues[i];
  delete [] PseudopotentialValues;
  delete [] TmpPseudopotentials;
  delete [] Distances;

  return FinalPseudopotentials;

#else
   cout << "EvaluateFiniteWidth requires linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}



// evalute pseudopotentials for coulomb interaction in a given Landau level with a given density profile,
// but without tabulating and interpolating the pseudopotentials
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// zDensity = density distribution of the layer 
// points = number of points where exact pseudopotentials are calculated
// multiplier = number of integration intervals used per point of discretization
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// zDensity2 = (optional) density distribution of layer 2, if absent, taken to be equal to 1st profile
// return value = array that contains the pseudopotentials

double* EvaluateFiniteWidthPseudoPotentialsNoInterpolation(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity, int points, double multiplier, double layerSeparation, AbstractZDensityProfile *zDensity2)
{
#ifdef HAVE_GSL
  cout.precision(14);
  double Z1min, Z1max, Width1, Z2min, Z2max, Width2;
  zDensity->GetSupport(Z1min, Z1max);
  if (zDensity2==0) zDensity2 = zDensity;
  zDensity2->GetSupport(Z2min, Z2max);
  Width1 = Z1max - Z1min;
  Width2 = Z2max - Z2min;

  PseudoPotentialGenerator VmGenerator(nbrFlux, landauLevel);
  
  // number of PP coefficients:
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* FinalPseudopotentials = new double [points];

  // integration workspaces
  gsl_integration_workspace * IntW1
    = gsl_integration_workspace_alloc ((int)(multiplier*points));
  gsl_integration_workspace * IntW2
    = gsl_integration_workspace_alloc ((int)(multiplier*points));


  DiagPseudoPotentials::paramZ2 InnerParameters;
  
  InnerParameters.d = layerSeparation;
  InnerParameters.Rho2 = zDensity2;  
  InnerParameters.VmAcc = NULL;
  InnerParameters.Vm = NULL;
  InnerParameters.VmGen = &VmGenerator;
  
  DiagPseudoPotentials::paramZ1 OuterParameters;
  OuterParameters.Z2min = Z2min ;
  OuterParameters.Z2max = Z2max;
  OuterParameters.Rho1 = zDensity;  
  OuterParameters.epsAbs = 0;
  OuterParameters.epsRel = 1e-8;
  OuterParameters.innerVar = &InnerParameters;
  gsl_function TheInnerIntegrand;
  TheInnerIntegrand.function = &DiagPseudoPotentials::IntegrandZ2NoInterpolation;
  TheInnerIntegrand.params = &InnerParameters;
  OuterParameters.innerIntegral = &TheInnerIntegrand;
  OuterParameters.sizeW = (int)(multiplier*points);
  OuterParameters.w=IntW2;

  gsl_function TheOuterIntegrand;
  TheOuterIntegrand.function = &DiagPseudoPotentials::IntegrandZ1;
  TheOuterIntegrand.params = &OuterParameters;

  double result, error;
  
  for (int l=0; l<=MaxMomentum; ++l)
    {
      InnerParameters.relativeM = MaxMomentum-l;
      // call the outer integration over z1
      gsl_integration_qag (&TheOuterIntegrand, Z1min, Z1max, OuterParameters.epsAbs, OuterParameters.epsRel,
			   OuterParameters.sizeW, /* KEY */ GSL_INTEG_GAUSS31, IntW1, &result, &error);        
      
      FinalPseudopotentials[l]=result;
      cout << "V_"<<l<<"="<<result<<" +/- "<< error<<endl;
    }

  // clean up
  
  gsl_integration_workspace_free (IntW1);
  gsl_integration_workspace_free (IntW2);
  
  return FinalPseudopotentials;

#else
   cout << "EvaluateFiniteWidth requires linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}
