#include "PseudoPotentialGenerator.h"

#include <cmath>
#include <iostream>

using std::sqrt;
using std::cout;
using std::endl;



// default constructor
PseudoPotentialGenerator::PseudoPotentialGenerator()
{
  this->NbrFlux=0;
}


// constructor
// nbrFlux = nbr Flux of sphere
// landauLevel = LL index
// layerSeparation = layer separation in magnetic lengths
PseudoPotentialGenerator::PseudoPotentialGenerator(int nbrFlux, int landauLevel, double layerSeparation)
{
  this->NbrFlux=nbrFlux;
  this->LandauLevel=landauLevel;
  this->MaxMomentum = nbrFlux + (landauLevel << 1);
  this->LayerSeparation=layerSeparation;
  this->MainCoefficients=new ClebschGordanCoefficients(MaxMomentum, MaxMomentum);
  this->Coefficients = new ClebschGordanCoefficients[MaxMomentum + 1];
  this->FormFactors = new double[MaxMomentum + 1];
  for (int l = 0; l <= MaxMomentum; ++l)
    this->Coefficients[l] = ClebschGordanCoefficients(MaxMomentum, l << 1);

}

// default destructor
PseudoPotentialGenerator::~PseudoPotentialGenerator()
{
  if (this->NbrFlux>0)
    {
      delete[] Coefficients;
      delete[] FormFactors;
      delete MainCoefficients;
    }
}

// auxiliary function to evaluate form-factors
void PseudoPotentialGenerator::EvaluateFormFactors()
{
  double dd = this->LayerSeparation*this->LayerSeparation;
  double base = ( sqrt(2.0*NbrFlux + dd) - LayerSeparation ) / ( sqrt(2.0*NbrFlux + dd) + LayerSeparation );
  FormFactors[0]=sqrt(base);
  for (int l = 1; l <= MaxMomentum; ++l)
    FormFactors[l] = FormFactors[l-1]*base;  
}



// evalute a pseudopotential for coulomb interaction in a given Landau level
//
// relativeM = relative angular momentum value
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// return value = V_m(d)
//
double PseudoPotentialGenerator::Evaluate(int relativeM, double layerSeparation)
{
  if (layerSeparation!=this->LayerSeparation)
    {
      this->LayerSeparation=layerSeparation;
      this->EvaluateFormFactors();
    }
  int l=relativeM;
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
	    TmpCoef2 = Coefficients[j].GetCoefficient(m1, m2 - m1, MaxMomentum) * Coefficients[j].GetCoefficient(NbrFlux, 0, MaxMomentum);
	    TmpCoef += FormFactors[j] * Sign * TmpCoef2 * TmpCoef2;
	    Sign *= -1.0;
	  }
	TmpPseudopotentials += (MainCoefficients->GetCoefficient(m1, -m1, l << 1) * 
				MainCoefficients->GetCoefficient(m2, -m2, l << 1) * TmpCoef);
      }
  return TmpPseudopotentials / sqrt (0.5 * NbrFlux);
}


// get a particular pseudopotential 
// 
double PseudoPotentialGenerator::Evaluate(int relativeM)
{
  return this->Evaluate(relativeM, this->LayerSeparation);
}
