#include "PseudoPotentials.h"

#include "Vector/RealVector.h"

#include "MathTools/Complex.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "PseudoPotentialGenerator.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Hamiltonian/ParticleOnSphereTwoLandauLevelDeltaHamiltonian.h"

#include <iostream>
#include <cstdlib>
#include <string>

using std::abs;
using std::cout;
using std::endl;
using std::string;
using std::sin;

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
      if ((nbrFlux & 1) ==0)
	--sizeRst;
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
			 parameters->epsAbs, parameters->epsRel/3.0, parameters->sizeW,
			 /* KEY */ GSL_INTEG_GAUSS15, parameters->w, &result, &error);        
//    gsl_integration_qags (parameters->innerIntegral, parameters->Z2min, parameters->Z2max,
//			 parameters->epsAbs, parameters->epsRel, parameters->sizeW,
//			 parameters->w, &result, &error);        


    return result*parameters->Rho1->GetValue(z1);
  }

  struct param1B
  {
    double TwoR;
    double NorthPoleCharge;
    double DSqrN;
    double SouthPoleCharge;
    double DSqrS;
    int OrbitalIndex;
    ParticleOnSphereGenericLLFunctionBasis *Basis;
    RealVector Value;
    Complex FctValue;
  };

  // integration of Coulomb interaction of a point charge
  double Integrand1B (double theta, void * params)
  {
    param1B *parameters=(param1B*)params;
    parameters->Value[0]=theta;
    parameters->Basis->GetFunctionValue(parameters->Value, parameters->FctValue, parameters->OrbitalIndex);
    double s=parameters->TwoR*sin((M_PI-theta)/2.);
    double s2=parameters->TwoR*sin(theta/2.);
    return 2.0*M_PI*sin(theta)*SqrNorm(parameters->FctValue)*(-parameters->NorthPoleCharge/sqrt(s*s+parameters->DSqrN)-parameters->SouthPoleCharge/sqrt(s2*s2+parameters->DSqrS));
  }
  
}

#endif


// evalute one body potentials for two charges +e located at given distance above the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPoleDistance = distance of charge located above the north pole
// southPoleDistance = distance of charge located above the south pole
// return value = array that conatins the pseudopotentials
double* EvaluateOneBodyCoulombPotentials(int nbrFlux, int landauLevel, double northPoleCharge, double northPoleDistance, double southPoleCharge, double southPoleDistance)
{
#ifdef HAVE_GSL
  ParticleOnSphereGenericLLFunctionBasis Basis (nbrFlux, landauLevel);

  DiagPseudoPotentials::param1B params;
  params.TwoR=sqrt((nbrFlux+2.0*landauLevel)*2.0);
  params.NorthPoleCharge = northPoleCharge;
  params.DSqrN = northPoleDistance*northPoleDistance;
  params.SouthPoleCharge = southPoleCharge;
  params.DSqrS = southPoleDistance*southPoleDistance;
  params.Basis=&Basis;
  params.Value.Resize(2);
  params.Value[1]=0.0;
  
  gsl_function TheIntegrand;
  TheIntegrand.function = &DiagPseudoPotentials::Integrand1B;
  TheIntegrand.params = &params;

  int NbrPoints=25*(nbrFlux+2*landauLevel);
  gsl_integration_workspace *Workspace = gsl_integration_workspace_alloc (NbrPoints);

  double result, error;

  double *Potentials = new double[nbrFlux+2*landauLevel+1];
    
  for (int i=0; i<=nbrFlux+2*landauLevel; ++i)
    {
      params.OrbitalIndex=i;
      // call the integration over theta
      gsl_integration_qag (&TheIntegrand, 0, M_PI, /*epsAbs*/ 0.0, /*epsRel*/ 1.0e-8,
      			   NbrPoints, /* KEY */ GSL_INTEG_GAUSS15, Workspace, &result, &error);
      Potentials[i]=result;
      //cout << "1B-Coulomb Potential["<<i<<"]="<<result<<endl;
    }
  gsl_integration_workspace_free (Workspace);
  return Potentials;
#else
  cout << "Need GSL libraries to evaluate 1-body pseudopotentials of charges."<<endl;
  return NULL;
#endif
}


// evalute pseudopotentials for coulomb interaction in a given Landau level with a given density profile
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// zDensity = density distribution of the layer 
// points = number of points where exact pseudopotentials are calculated
// multiplier = number of integration intervals used per point of discretization
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// zDensity2 = (optional) density distribution of layer 2, if absent, taken to be equal to 1st profile
// epsRel = tolerance given to integration routine
// return value = array that contains the pseudopotentials

double* EvaluateFiniteWidthPseudoPotentials(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity, int points, double multiplier, double layerSeparation, AbstractZDensityProfile *zDensity2, double epsRel)
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
  OuterParameters.epsRel = epsRel;
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
      			   OuterParameters.sizeW, /* KEY */ GSL_INTEG_GAUSS15, IntW1, &result, &error);        
//	gsl_integration_qags (&TheOuterIntegrand, Z1min, Z1max, OuterParameters.epsAbs, OuterParameters.epsRel,
  //                    OuterParameters.sizeW, IntW1, &result, &error);              

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
// epsRel = tolerance given to integration routine
// return value = array that contains the pseudopotentials

double* EvaluateFiniteWidthPseudoPotentialsNoInterpolation(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity, int points, double multiplier, double layerSeparation, AbstractZDensityProfile *zDensity2, double epsRel)
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
  OuterParameters.epsRel = epsRel;
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
			   OuterParameters.sizeW, /* KEY */ GSL_INTEG_GAUSS61, IntW1, &result, &error);        
      
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

// evaluate graphene bilayer pseudopotentials
//
// nbrflux = number of flux quanta for the lowest Landau level (i.e. twice the maximum momentum for a single particle
// l1, l2, l3, l4 = Landau level indices (0 for lowest LL, 1 for first excited LL)
//
double* EvaluateGrapheneBilayerPseudopotentials(int nbrFlux, int& nbrPseudopotentials, int llindex1, int llindex2, int llindex3, int llindex4, bool verbose)
{
  cout.precision(14); 
  double Q = 0.5 * nbrFlux;
  double l1 = Q + llindex1;
  double l2 = Q + llindex2;
  double l3 = Q + llindex3;
  double l4 = Q + llindex4;

  double Lmin = ((fabs(l1 - l2) > fabs(l3 - l4)) ? fabs(l1 - l2) : fabs(l3 - l4));
  double Lmax = ((l1 + l2 < l3 + l4) ? (l1 + l2) : (l3 + l4));

  nbrPseudopotentials = ((int)(Lmax - Lmin) + 1);

  cout<<"Lmin: "<<Lmin<<" LMax: "<<Lmax<<" Nbr= "<<nbrPseudopotentials<<endl;

  double* Pseudopotentials = new double [nbrPseudopotentials];

  //Prepare the tables of CG coeffs to speed up calculation

 
  ClebschGordanCoefficients Coefficients12((int)(2.0*l1), (int)(2.0*l2));
  ClebschGordanCoefficients Coefficients34((int)(2.0*l3), (int)(2.0*l4));

  double min12 = l1 < l2 ? l1 : l2;
  double min34 = l3 < l4 ? l3 : l4;
 
  double minj=l1;
  for (double m1 = -min12; m1<= min12; m1 += 1.0)
   for (double m2 = -min34; m2 <= min34; m2 += 1.0)
    if (fabs(m1-m2)<minj)
     minj = fabs(m1-m2);
  double maxj = ((l1 + l2) < (l3 + l4) ? (l1+l2) : (l3+l4));
  int length = (int) (maxj - minj) + 1;

  ClebschGordanCoefficients* Coefficients1j = new ClebschGordanCoefficients[length];
  ClebschGordanCoefficients* Coefficients2j = new ClebschGordanCoefficients[length];
  ClebschGordanCoefficients* Coefficients3j = new ClebschGordanCoefficients[length];
  ClebschGordanCoefficients* Coefficients4j = new ClebschGordanCoefficients[length];

  for (int j = (int)minj; j<= (int)maxj; j++)
   {
     Coefficients1j[j] = ClebschGordanCoefficients((int)(2.0*l1), 2*j); 
     Coefficients2j[j] = ClebschGordanCoefficients((int)(2.0*l2), 2*j); 
     Coefficients3j[j] = ClebschGordanCoefficients((int)(2.0*l3), 2*j); 
     Coefficients4j[j] = ClebschGordanCoefficients((int)(2.0*l4), 2*j); 
   }


  //Proceed to calculate pseudopotentials

  for (double l = 0.0; l <= (Lmax - Lmin); l += 1.0)
    {
     double TmpPseudopotentials = 0.0;
 
     for (double m1 = -min12; m1<= min12; m1 += 1.0)
      for (double m2 = -min34; m2 <= min34; m2 += 1.0)
       {
         double Sign = pow(-1.0, m2 - m1 + l1 + l2 + l3 + l4);
         double CGTerm1 = Coefficients12.GetCoefficient((int)(2.0*m1), (int)(-2.0*m1), (int)(2.0*(Lmax - l))) * Coefficients34.GetCoefficient((int)(2.0*m2), (int)(-2.0*m2), (int)(2.0*(Lmax - l)));
         
	 double abs_m12 = fabs(m1-m2);
         
         double min1234 = (l1 + l2) < (l3 + l4) ? l1 + l2 : l3 + l4;
         double CGTerm2 = 0.0;
         for (double j = abs_m12; j <= min1234; j+=1.0)
           if ( (l1 <= l3 + j) && (l1 >= fabs(l3-j)) && (l2 <= l4 + j) && (l2 >= fabs(l4-j)) && (l3 <= l1 + j) && (l3 >= fabs(l1-j)) && (l4 <= l2 + j) && (l4 >= fabs(l2-j))  )
              CGTerm2 += ( Coefficients1j[(int)j].GetCoefficient((int)(-2.0*m1), (int)(-2.0*(m2-m1)), (int)(2.0*l3)) * Coefficients2j[(int)j].GetCoefficient( (int)(2.0*m1), (int)(2.0*(m2-m1)), (int)(2.0*l4)) * Coefficients3j[(int)j].GetCoefficient( (int)(2.0*Q), 0, (int)(2.0*l1)) * Coefficients4j[(int)j].GetCoefficient( (int)(2.0*Q), 0, (int)(2.0*l2)) );

         TmpPseudopotentials += (Sign * CGTerm1 * CGTerm2);
       } 
      Pseudopotentials[(int)l] = TmpPseudopotentials / sqrt (Q); 
      if (verbose == true)
	cout << Pseudopotentials[(int)l] << "  ";
  }

  cout << endl;

  delete[] Coefficients1j;      
  delete[] Coefficients2j;      
  delete[] Coefficients3j;      
  delete[] Coefficients4j;      
  return Pseudopotentials;
}

// evalute pseudopotentials for dipolar interaction in the lowest Landau level Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double* EvaluateDipolarPseudopotentials(int nbrFlux, bool quiet)
{
  double* Pseudopotentials = new double [nbrFlux + 1];
  double Factor = 1.0 / sqrt(0.5 * ((double) nbrFlux));
  Factor *= Factor * Factor;
  Factor *= 0.5 * ((double) (nbrFlux + 1)) * ((double) (nbrFlux + 1));
  BinomialCoefficients Binomial (4 * nbrFlux);
  for (int j = 0; j < nbrFlux; ++j)
    {
      Pseudopotentials[nbrFlux - j] = (Factor / (((double) (nbrFlux + 1 + j)) * (((double) (nbrFlux - j))))) * 
	(Binomial.GetNumericalCoefficient(2 *(nbrFlux + j), nbrFlux + j) / Binomial(2 * nbrFlux, nbrFlux)) * 
	(Binomial.GetNumericalCoefficient(2 *(nbrFlux - j - 1), nbrFlux - j - 1) / Binomial.GetNumericalCoefficient(2 * nbrFlux, nbrFlux));
      if (quiet == false)
	cout << "V[" << (nbrFlux - j) << "] = " << Pseudopotentials[nbrFlux - j] << endl;
    }
  Pseudopotentials[0] = 0.0;
  return Pseudopotentials;
}

// ZP: Not very fast but it works... CAUTION: Pseudopotentials are not normalized in any particularly nice way, but they will produce a zero energy GS at nu=2/5 with cyclotron energy=0
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (2j+1) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for delta interaction with two Landau levels on Sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double** Evaluate2LLSphereDeltaPseudopotentials(int nbrFlux, bool quiet)
{
  int NbrFluxQuanta = nbrFlux;    
  int LzMaxUp = NbrFluxQuanta + 2;
  int LzMaxDown = NbrFluxQuanta;
  int TwoQ = NbrFluxQuanta;

  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[10] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
          "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
          "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown", "PseudopotentialsUpDownDownUp"};
  
  // these are the lenghts of the arrays corresponding to the labels above.           
  int PseudoLengths[10] = {TwoQ + 3, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ + 1}; 
  
  // these are the lenghts of the arrays corresponding to the labels above.           
  int PseudoMins[10] = {0, 0, 1, 0, 0, 1, 1, 1, 1, 1}; 
  
  int LLs[10][4] = {{1,1,1,1}, {1,1,0,0}, {1,1,1,0}, {0,0,1,1}, {0,0,0,0}, {0,0,1,0}, {1,0,1,1}, {1,0,0,0}, {1,0,1,0}, {1,0,0,1}}; 
    
  double **Pseudopotentials;
  Pseudopotentials = new double*[10];

  for(int pp = 0; pp < 10; pp++)
    {
  Pseudopotentials[pp] = new double[PseudoLengths[pp]]; 
        int n1 = LLs[pp][0], n2 = LLs[pp][1], n3 = LLs[pp][2], n4 = LLs[pp][3];
        if (!quiet) cout << "Pseudopotential: "<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<endl;
        int Twol1 = TwoQ + 2 * n1;
        int Twol2 = TwoQ + 2 * n2;
        int Twol3 = TwoQ + 2 * n3;
        int Twol4 = TwoQ + 2 * n4;
    
        ClebschGordanCoefficients Clebsch12 (Twol1, Twol2);
        ClebschGordanCoefficients Clebsch34 (Twol3, Twol4);
  
        int Twom1Min = -Twol1;
        if (Twom1Min < -Twol2)
          Twom1Min = -Twol2;

        int Twom1Max = Twol1;
        if (Twom1Max > Twol2)
          Twom1Max = Twol2;

        int Twom2Min = -Twol3;
        if (Twom2Min < -Twol4)
          Twom2Min = -Twol4;

        int Twom2Max = Twol3;
        if (Twom2Max > Twol4)
          Twom2Max = Twol4;

        int Twojmin = abs(Twol1-Twol3);
        if (Twojmin < abs(Twol2-Twol4))
          Twojmin = abs(Twol2-Twol4);

        int Twojmax = (Twol1+Twol3);
        if (Twojmax > (Twol2+Twol4))
          Twojmax = (Twol2+Twol4);

        //cout << "Lmax= "<<PseudoLengths[pp] - 1 + PseudoMins[pp]<<" Lmin= "<<PseudoMins[pp]<<endl;
        for ( int L = PseudoLengths[pp] - 1 + PseudoMins[pp] ; L >= PseudoMins[pp] ; L--)
     {
           int idx = PseudoLengths[pp] - 1 + PseudoMins[pp] - L;  
           double TmpVL = 0.0;
           for (int Twoj = Twojmin; Twoj <= Twojmax; Twoj+=2)
             {
               ClebschGordanCoefficients Clebsch1j (Twol1, Twoj);
               ClebschGordanCoefficients Clebsch2j (Twol2, Twoj);
               ClebschGordanCoefficients Clebsch3j (Twol3, Twoj);
               ClebschGordanCoefficients Clebsch4j (Twol4, Twoj);
               for(int Twom1 = Twom1Min; Twom1 <= Twom1Max; Twom1+=2)
                 {
                   for(int Twom2 = Twom2Min; Twom2 <= Twom2Max; Twom2+=2)
                     {

                       TmpVL += (Twoj + 1) * pow(-1.0, 0.5 * Twom2 - 0.5 * Twom1 + 0.5 * Twol1 + 0.5 * Twol2 + 0.5 * Twol3 + 0.5 * Twol4)
        * Clebsch12.CarefulGetCoefficient(Twom1, -Twom1, 2*L) * Clebsch34.CarefulGetCoefficient(Twom2, -Twom2, 2*L) 
                                * Clebsch1j.CarefulGetCoefficient(-Twom1, -(Twom2-Twom1), Twol3) * Clebsch2j.CarefulGetCoefficient(Twom1, (Twom2-Twom1), Twol4)
                                * Clebsch3j.CarefulGetCoefficient(TwoQ, 0, Twol1) * Clebsch4j.CarefulGetCoefficient(TwoQ, 0, Twol2);
                    }
                }  
              }
            if (!quiet) cout<<TmpVL<<" ";
            Pseudopotentials[pp][idx] = TmpVL;
          } 
         if (!quiet) cout<<endl; 
    }
  
  return Pseudopotentials;
}

// ZP: Not very fast but it works... CAUTION: Pseudopotentials are not normalized in any particularly nice way, but they will produce a zero energy GS at nu=2/5 with cyclotron energy=0
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (-j(j+1)/Q) * (2j+1) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for delta interaction with two Landau levels on Sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double** Evaluate2LLSphereLaplacianDeltaPseudopotentials(int nbrFlux, bool quiet)
{
  int NbrFluxQuanta = nbrFlux;    
  int LzMaxUp = NbrFluxQuanta + 2;
  int LzMaxDown = NbrFluxQuanta;
  int TwoQ = NbrFluxQuanta;

  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[10] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
			    "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
			    "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown", "PseudopotentialsUpDownDownUp"};
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoLengths[10] = {TwoQ + 3, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ + 1}; 
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoMins[10] = {0, 0, 1, 0, 0, 1, 1, 1, 1, 1}; 
  
  int LLs[10][4] = {{1,1,1,1}, {1,1,0,0}, {1,1,1,0}, {0,0,1,1}, {0,0,0,0}, {0,0,1,0}, {1,0,1,1}, {1,0,0,0}, {1,0,1,0}, {1,0,0,1}}; 
    
  double **Pseudopotentials;
  Pseudopotentials = new double*[10];

  for(int pp = 0; pp < 10; pp++)
    {
	Pseudopotentials[pp] = new double[PseudoLengths[pp]]; 
        int n1 = LLs[pp][0], n2 = LLs[pp][1], n3 = LLs[pp][2], n4 = LLs[pp][3];
        if (!quiet) cout << "Pseudopotential: "<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<endl;
        int Twol1 = TwoQ + 2 * n1;
        int Twol2 = TwoQ + 2 * n2;
        int Twol3 = TwoQ + 2 * n3;
        int Twol4 = TwoQ + 2 * n4;
    
        ClebschGordanCoefficients Clebsch12 (Twol1, Twol2);
        ClebschGordanCoefficients Clebsch34 (Twol3, Twol4);
  
        int Twom1Min = -Twol1;
        if (Twom1Min < -Twol2)
          Twom1Min = -Twol2;

        int Twom1Max = Twol1;
        if (Twom1Max > Twol2)
          Twom1Max = Twol2;

        int Twom2Min = -Twol3;
        if (Twom2Min < -Twol4)
          Twom2Min = -Twol4;

        int Twom2Max = Twol3;
        if (Twom2Max > Twol4)
          Twom2Max = Twol4;

        int Twojmin = abs(Twol1-Twol3);
        if (Twojmin < abs(Twol2-Twol4))
          Twojmin = abs(Twol2-Twol4);

        int Twojmax = (Twol1+Twol3);
        if (Twojmax > (Twol2+Twol4))
          Twojmax = (Twol2+Twol4);

        //cout << "Lmax= "<<PseudoLengths[pp] - 1 + PseudoMins[pp]<<" Lmin= "<<PseudoMins[pp]<<endl;
        for ( int L = PseudoLengths[pp] - 1 + PseudoMins[pp] ; L >= PseudoMins[pp] ; L--)
  	 {
           int idx = PseudoLengths[pp] - 1 + PseudoMins[pp] - L;	
           double TmpVL = 0.0;
           for (int Twoj = Twojmin; Twoj <= Twojmax; Twoj+=2)
             {
               ClebschGordanCoefficients Clebsch1j (Twol1, Twoj);
               ClebschGordanCoefficients Clebsch2j (Twol2, Twoj);
               ClebschGordanCoefficients Clebsch3j (Twol3, Twoj);
               ClebschGordanCoefficients Clebsch4j (Twol4, Twoj);
               for(int Twom1 = Twom1Min; Twom1 <= Twom1Max; Twom1+=2)
                 {
                   for(int Twom2 = Twom2Min; Twom2 <= Twom2Max; Twom2+=2)
                     {
                       TmpVL += (-0.5*Twoj*(0.5*Twoj+1)/(0.5*TwoQ)) * (Twoj + 1) * pow(-1.0, 0.5 * Twom2 - 0.5 * Twom1 + 0.5 * Twol1 + 0.5 * Twol2 + 0.5 * Twol3 + 0.5 * Twol4)
				* Clebsch12.CarefulGetCoefficient(Twom1, -Twom1, 2*L) * Clebsch34.CarefulGetCoefficient(Twom2, -Twom2, 2*L) 
                                * Clebsch1j.CarefulGetCoefficient(-Twom1, -(Twom2-Twom1), Twol3) * Clebsch2j.CarefulGetCoefficient(Twom1, (Twom2-Twom1), Twol4)
                                * Clebsch3j.CarefulGetCoefficient(TwoQ, 0, Twol1) * Clebsch4j.CarefulGetCoefficient(TwoQ, 0, Twol2);
                    }
                }  
              }
            if (!quiet) cout<<TmpVL<<" ";
            Pseudopotentials[pp][idx] = TmpVL;
          } 
         if (!quiet) cout<<endl; 
    }
	
  return Pseudopotentials;
}

// ZP: Not very fast but it works... 
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (1/sqrt(Q)) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for Coulomb interaction with two Landau levels on sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials

double** Evaluate2LLSphereCoulombPseudopotentials(int nbrFlux, bool quiet)
{
  int NbrFluxQuanta = nbrFlux;    
  int LzMaxUp = NbrFluxQuanta + 2;
  int LzMaxDown = NbrFluxQuanta;
  int TwoQ = NbrFluxQuanta;

  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[10] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
			    "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
			    "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown", "PseudopotentialsUpDownDownUp"};
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoLengths[10] = {TwoQ + 3, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ, TwoQ + 1, TwoQ + 1}; 
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoMins[10] = {0, 0, 1, 0, 0, 1, 1, 1, 1, 1}; 
  
  int LLs[10][4] = {{1,1,1,1}, {1,1,0,0}, {1,1,1,0}, {0,0,1,1}, {0,0,0,0}, {0,0,1,0}, {1,0,1,1}, {1,0,0,0}, {1,0,1,0}, {1,0,0,1}}; 
    
  double **Pseudopotentials;
  Pseudopotentials = new double*[10];

  for(int pp = 0; pp < 10; pp++)
    {
	Pseudopotentials[pp] = new double[PseudoLengths[pp]]; 
        int n1 = LLs[pp][0], n2 = LLs[pp][1], n3 = LLs[pp][2], n4 = LLs[pp][3];
        if (!quiet) cout << "Pseudopotential: "<<n1<<" "<<n2<<" "<<n3<<" "<<n4<<endl;
        int Twol1 = TwoQ + 2 * n1;
        int Twol2 = TwoQ + 2 * n2;
        int Twol3 = TwoQ + 2 * n3;
        int Twol4 = TwoQ + 2 * n4;
    
        ClebschGordanCoefficients Clebsch12 (Twol1, Twol2);
        ClebschGordanCoefficients Clebsch34 (Twol3, Twol4);
  
        int Twom1Min = -Twol1;
        if (Twom1Min < -Twol2)
          Twom1Min = -Twol2;

        int Twom1Max = Twol1;
        if (Twom1Max > Twol2)
          Twom1Max = Twol2;

        int Twom2Min = -Twol3;
        if (Twom2Min < -Twol4)
          Twom2Min = -Twol4;

        int Twom2Max = Twol3;
        if (Twom2Max > Twol4)
          Twom2Max = Twol4;

        int Twojmin = abs(Twol1-Twol3);
        if (Twojmin < abs(Twol2-Twol4))
          Twojmin = abs(Twol2-Twol4);

        int Twojmax = (Twol1+Twol3);
        if (Twojmax > (Twol2+Twol4))
          Twojmax = (Twol2+Twol4);

        //cout << "Lmax= "<<PseudoLengths[pp] - 1 + PseudoMins[pp]<<" Lmin= "<<PseudoMins[pp]<<endl;
        for ( int L = PseudoLengths[pp] - 1 + PseudoMins[pp] ; L >= PseudoMins[pp] ; L--)
  	 {
           int idx = PseudoLengths[pp] - 1 + PseudoMins[pp] - L;	
           double TmpVL = 0.0;
           for (int Twoj = Twojmin; Twoj <= Twojmax; Twoj+=2)
             {
               ClebschGordanCoefficients Clebsch1j (Twol1, Twoj);
               ClebschGordanCoefficients Clebsch2j (Twol2, Twoj);
               ClebschGordanCoefficients Clebsch3j (Twol3, Twoj);
               ClebschGordanCoefficients Clebsch4j (Twol4, Twoj);
               for(int Twom1 = Twom1Min; Twom1 <= Twom1Max; Twom1+=2)
                 {
                   for(int Twom2 = Twom2Min; Twom2 <= Twom2Max; Twom2+=2)
                     {

                       TmpVL += (1.0/sqrt(0.5 * TwoQ)) * pow(-1.0, 0.5 * Twom2 - 0.5 * Twom1 + 0.5 * Twol1 + 0.5 * Twol2 + 0.5 * Twol3 + 0.5 * Twol4)
				* Clebsch12.CarefulGetCoefficient(Twom1, -Twom1, 2*L) * Clebsch34.CarefulGetCoefficient(Twom2, -Twom2, 2*L) 
                                * Clebsch1j.CarefulGetCoefficient(-Twom1, -(Twom2-Twom1), Twol3) * Clebsch2j.CarefulGetCoefficient(Twom1, (Twom2-Twom1), Twol4)
                                * Clebsch3j.CarefulGetCoefficient(TwoQ, 0, Twol1) * Clebsch4j.CarefulGetCoefficient(TwoQ, 0, Twol2);
                    }
                }  
              }
            if (!quiet) cout<<TmpVL<<" ";
            Pseudopotentials[pp][idx] = TmpVL;
          } 
         if (!quiet) cout<<endl; 
    }
	
  return Pseudopotentials;
}

// evaluate pseudopotentials for triangular well on sphere using the midpoint method
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// width = thickness of the quantum well in units of the magnetic length
// bias = bias of the quantum well in magnetic energy unit
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// pseudospins = array that contains pseudospins (0 for up, 1 for down) ordered as written below
// return value = array that contains pseudopotentials

double* EvaluateTriangularWellPseudopotentials(int nbrFlux, int landauLevel, double width, double bias, int* pseudospins, int nbrPointsInteg, bool quiet)
{
  cout.precision(14);

  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* Pseudopotentials = new double [MaxMomentum+1];
  int m;
  for ( m=0; m <= MaxMomentum; ++m ) Pseudopotentials[m]=0; 

  // One has to compute \int_0^w dz1 \int_0^w dz2 Psi_p1(z1) Psi_p2(z2) Vm(z1-z2) Psi_p3(z1) Psi_p4(z2)
  // the index (1/2) being the pseudospin, i.e. the confinement subband number here
  // and Vm(z) is the mth pseudopotential of the bilayer with distance z
  double delta_z= width / nbrPointsInteg;
  double z,Vm;
  double integrand1,integrand2;
  
  double z1,z2;
  double* V;
  int i,j;
  if ( bias == 0 )  // Square quantum well: \psi_n(z)=\sqrt(2/w)\sin(n\pi x/w)
  {
    for ( i=0; i < nbrPointsInteg; ++i )
    {
      z1=(i+0.5)*delta_z;
      integrand1=sin(M_PI*z1/width*(1.+pseudospins[0]))*sin(M_PI*z1/width*(1.+pseudospins[2]));
      for ( j=0; j < nbrPointsInteg; ++j )
        {
          z2=(j+0.5)*delta_z;
	  integrand2=sin(M_PI*z2/width*(1.+pseudospins[1]))*sin(M_PI*z2/width*(1.+pseudospins[3]));
          z=abs(z1-z2);
          V = EvaluatePseudopotentials(nbrFlux, landauLevel, z, true);// Bilayer pseudopotentials
	  for ( m=0; m <= MaxMomentum; ++m) Pseudopotentials[m] += integrand1*integrand2*V[m] ;
        }
      //for ( m=0; m < nbrFlux; ++m) Pseudopotentials[m] *= integrand1 ;
    }
  
    for ( m=0; m <= MaxMomentum; ++m) Pseudopotentials[m] *= 4*delta_z*delta_z/width/width ;
  }
  else
  {
    TriangularWellEigenfunction Well(width, bias);

    for ( i=0; i < nbrPointsInteg; ++i )
      {
	z1=(i+0.5)*delta_z;
	integrand1=Well.Eigenfunction(z1,pseudospins[0])*Well.Eigenfunction(z1,pseudospins[2]);
	for ( j=0; j < nbrPointsInteg; ++j )
	  {
	    z2=(j+0.5)*delta_z;
	    integrand2=Well.Eigenfunction(z2,pseudospins[1])*Well.Eigenfunction(z2,pseudospins[3]);
	    z=abs(z1-z2);
	    V = EvaluatePseudopotentials(nbrFlux, landauLevel, z, true);// Bilayer pseudopotentials
	    for ( m=0; m <= MaxMomentum; ++m) Pseudopotentials[m] += integrand1*integrand2*V[m] ;
	  }
	//for ( m=0; m < nbrFlux; ++m) Pseudopotentials[m] *= integrand1 ;
      }
    
    for ( m=0; m <= MaxMomentum; ++m) Pseudopotentials[m] *= delta_z*delta_z ;
  }

  if (quiet == false) 
    {
       for ( m=0; m <= MaxMomentum; ++m) 
         cout << "V[" << m << "] = " << Pseudopotentials[m] << endl;	
    }

  return Pseudopotentials;
}
