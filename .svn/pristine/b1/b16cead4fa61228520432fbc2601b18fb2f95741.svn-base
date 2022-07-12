#include "config.h"
#include "MathTools/Complex.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// evaluate the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
//
// coefficients = bi-dimensional array where the numerical coefficients have to be stored
// nbrStates = number of states in the Lowest Landau Level
void EvaluateNumericalCoefficients (double*** coefficients, int nbrStates);

// evaluate the derivative of the Phi^4 contribution
//
// waveFunctionCoefficients = array containing the wave function coefficients
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = Phi^4 contribution
double EvaluatePhi4 (Complex* waveFunctionCoefficients, double*** numericalCoefficients, int nbrStates);

// evaluate the derivative of the Phi^4 contribution with respect to one of the wave function coefficients
//
// waveFunctionCoefficients = array containing the wave function coefficients
// index = index of the contribution with respect to one of the wave function coefficients for which derivative has to be evaluated
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = derivative of the Phi^4 contribution
Complex EvaluatePhi4Derivative (Complex* waveFunctionCoefficients, int index, double*** numericalCoefficients, int nbrStates);



int main(int argc, char** argv)
{
  cout.precision(14);

  int NbrStates = 11;
  int MaxNbrIteration = 1000;
  double Precision = 1e-8;
  int ResetLimit = 100;
  double Step = 0.01;

  double*** Coefficients = new double** [NbrStates];
  for (int i = 0; i < NbrStates; ++i)
    {
      Coefficients[i] = new double* [i + 1];
      for (int j = 0; j <= i ; ++j)
	Coefficients[i][j] = new double [NbrStates];
    }
  EvaluateNumericalCoefficients(Coefficients, NbrStates);
  Complex* CurrentX = new Complex[NbrStates];
  Complex* PreviousX = new Complex[NbrStates];
  CurrentX[0].Re = 2.0 * (drand48() - 0.5);
  double TmpNorm = CurrentX[0].Re;
  for (int i = 1; i < NbrStates; ++i)
    {
      CurrentX[i].Re = 2.0 * (drand48() - 0.5);
      CurrentX[i].Im = 2.0 * (drand48() - 0.5);
      TmpNorm += (CurrentX[i].Re * CurrentX[i].Re) + (CurrentX[i].Im * CurrentX[i].Re);
    }
  TmpNorm = 1.0 / sqrt(TmpNorm);
  CurrentX[0].Re *= TmpNorm;
  for (int i = 1; i < NbrStates; ++i)
    {
      CurrentX[i] *= TmpNorm;
    }

  double CurrentEnergy = EvaluatePhi4(CurrentX, Coefficients, NbrStates);
  cout << "energy: " << CurrentEnergy << endl;
  double PreviousEnergy = CurrentEnergy;
  Complex* CurrentDerivative = new Complex[NbrStates];
  Complex* PreviousDerivative = new Complex[NbrStates];
  Complex* CurrentDirection = new Complex[NbrStates];
  Complex* PreviousDirection = new Complex[NbrStates];
  for (int i = 0; i < NbrStates; ++i)
    {      
      CurrentDerivative[i] = EvaluatePhi4Derivative (CurrentX, i, Coefficients, NbrStates);
      PreviousDerivative[i] = CurrentDerivative[i];
      CurrentDirection[i] = -CurrentDerivative[i];
      PreviousDirection[i] = CurrentDirection[i];
    }
  double CurrentLagrangeMultiplier = 0.0;
  double PreviousLagrangeMultiplier = 0.0;
  double CurrentLagrangeMultiplierDerivative = 0.0;
  double PreviousLagrangeMultiplierDerivative = 0.0;
  double CurrentLagrangeMultiplierDirection = 0.0;


  double DirectionNorm = (CurrentDerivative[0].Re * CurrentDerivative[0].Re);
  for (int i = 1; i < NbrStates; ++i)
    DirectionNorm += (CurrentDerivative[i].Re * CurrentDerivative[i].Re) + (CurrentDerivative[i].Im * CurrentDerivative[i].Im);
  DirectionNorm += CurrentLagrangeMultiplierDerivative * CurrentLagrangeMultiplierDerivative;
  DirectionNorm = sqrt(DirectionNorm);

  int NbrIter = 1;
  double BetaFactor = 0.0;
  while ((NbrIter < MaxNbrIteration) && (DirectionNorm > Precision))
    {
      PreviousX[0].Re = CurrentX[0].Re;
      CurrentX[0].Re += (Step * CurrentDirection[0].Re);
      
      for (int i = 1; i < NbrStates; ++i)
	{
	  PreviousX[i] = CurrentX[i];
	  CurrentX[i] += (Step * CurrentDirection[i]);
	}
      PreviousLagrangeMultiplier = CurrentLagrangeMultiplier;
      CurrentLagrangeMultiplier += Step * CurrentLagrangeMultiplierDirection;

      PreviousLagrangeMultiplierDerivative = CurrentLagrangeMultiplierDerivative;
      CurrentLagrangeMultiplierDerivative = -1.0;
      for (int i = 0; i < NbrStates; ++i)
	{
	  PreviousDerivative[i] = CurrentDerivative[i];
	  CurrentDerivative[i] = EvaluatePhi4Derivative (CurrentX, i, Coefficients, NbrStates);
	  CurrentDerivative[i] += (2.0 * CurrentLagrangeMultiplier) * CurrentX[i];
	  cout << CurrentDerivative[i] << " ";
	  CurrentLagrangeMultiplierDerivative += (CurrentX[i].Re * CurrentX[i].Re) + (CurrentX[i].Im * CurrentX[i].Re);
	}
      cout << endl;

      PreviousEnergy = CurrentEnergy;
      CurrentEnergy = EvaluatePhi4(CurrentX, Coefficients, NbrStates);

      
      if ((NbrIter % ResetLimit) == 0)
	{
	  BetaFactor = 0.0;
	}
      else
	{
	  BetaFactor = CurrentDerivative[0].Re * CurrentDerivative[0].Re;
	  TmpNorm = CurrentDerivative[0].Re * CurrentDerivative[0].Re;
	  for (int i = 1; i < NbrStates; ++i)
	    {
	      BetaFactor += (CurrentDerivative[i].Re - PreviousDerivative[i].Re) * CurrentDerivative[i].Re;
	      BetaFactor += (CurrentDerivative[i].Im - PreviousDerivative[i].Im) * CurrentDerivative[i].Im;	      
	    }
	  cout << BetaFactor << " " << TmpNorm << endl;
	  BetaFactor += (CurrentLagrangeMultiplierDerivative - PreviousLagrangeMultiplierDerivative) * CurrentLagrangeMultiplierDerivative;
	  TmpNorm +=  PreviousLagrangeMultiplierDerivative * PreviousLagrangeMultiplierDerivative;
	  cout << BetaFactor << " " << TmpNorm << endl;
	  BetaFactor /= TmpNorm;
	}

      CurrentDirection[0].Re *= BetaFactor;
      CurrentDirection[0].Re -= CurrentDerivative[0].Re;
      for (int i = 1; i < NbrStates; ++i)
	{
	  CurrentDirection[i] *= BetaFactor;
	  CurrentDirection[i] -= CurrentDerivative[i];
	}    
      CurrentLagrangeMultiplierDirection *= BetaFactor;
      CurrentLagrangeMultiplierDirection -= CurrentLagrangeMultiplierDerivative;
      DirectionNorm = (CurrentDerivative[0].Re * CurrentDerivative[0].Re);
      for (int i = 1; i < NbrStates; ++i)
	DirectionNorm += (CurrentDerivative[i].Re * CurrentDerivative[i].Re) + (CurrentDerivative[i].Im * CurrentDerivative[i].Im);
      DirectionNorm += CurrentLagrangeMultiplierDerivative * CurrentLagrangeMultiplierDerivative;
      DirectionNorm = sqrt(DirectionNorm);
      
      cout << NbrIter << ": " << CurrentEnergy << " " << DirectionNorm << " " << BetaFactor << " " << CurrentLagrangeMultiplierDerivative << endl;
      ++NbrIter;
    }

  delete[] PreviousDirection;
  delete[] CurrentDirection;
  delete[] PreviousDerivative;  
  delete[] CurrentDerivative;  
  delete[] PreviousX;
  delete[] CurrentX;

}


// evaluate the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
//
// coefficients = bi-dimensional array where the numerical coefficients have to be stored
// nbrStates = number of states in the Lowest Landau Level

void EvaluateNumericalCoefficients (double*** coefficients, int nbrStates)
{
  int TwiceS = nbrStates - 1;
  double Factor = 2.0 / M_PI;
  FactorialCoefficient Coef;
  int i3;
  int i4;
  int Max;
  double* TmpArray;
  for (int i1 = 0; i1 < nbrStates; ++i1)
    {
      for (int i2 = 0; i2 <= i1; ++i2)
	{
	  TmpArray = coefficients[i1][i2];
	  Max = i1 + i2;
	  i3 = Max;
	  if ((i3 & 1) == 0)
	    i3 >>= 1;
	  else
	    {
	      i3 >>= 1;
	      i3 += 1;
	    }
	  if (Max >= nbrStates)
	    Max = nbrStates - 1;
	  for (; i3 <= Max; ++i3)
	    {
	      i4 =  i1 + i2 - i3;
	      Coef.SetToOne();
	      Coef.PartialFactorialMultiply(nbrStates - i1 , 2 * TwiceS - (i1 + i2));
	      Coef.PartialFactorialMultiply(nbrStates - i2 , 2 * TwiceS - (i1 + i2));
	      Coef.FactorialDivide(i1);
	      Coef.FactorialDivide(i2);
	      Coef.FactorialMultiply(i3 + i4);
	      Coef.FactorialMultiply(i3 + i4);
	      Coef.FactorialDivide(TwiceS - i3);
	      Coef.FactorialDivide(TwiceS - i4);
	      Coef.FactorialDivide(i3);
	      Coef.FactorialDivide(i4);
	      Coef.PartialFactorialDivide(nbrStates + 1, (2 * TwiceS) + 1);	
	      Coef.FactorialMultiply(nbrStates);
	      Coef.PartialFactorialDivide(nbrStates + 1, (2 * TwiceS) + 1);	
	      Coef.FactorialMultiply(nbrStates);
	      TmpArray[i3] = sqrt(Coef.GetNumericalValue() * Factor);
	      if (i1 == i2)
		TmpArray[i3] *= 0.5;
	      if (i3 == i4)
		TmpArray[i3] *= 0.5;	      
//	      cout << TmpArray[i3] << endl;
	    }
	}
   }
}

// evaluate the derivative of the Phi^4 contribution
//
// waveFunctionCoefficients = array containing the wave function coefficients
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = Phi^4 contribution

double EvaluatePhi4 (Complex* waveFunctionCoefficients, double*** numericalCoefficients, int nbrStates)
{
  double Phi4 = 0.0;
  int i3;
  int i4;
  Complex Tmp;
  Complex Tmp2;
  int Max;
  double* TmpArray;
  for (int i1 = 0; i1 < nbrStates; ++i1)
    for (int i2 = 0; i2 <= i1; ++i2)
      {
	TmpArray = numericalCoefficients[i1][i2];
	Tmp = (waveFunctionCoefficients[i1] * waveFunctionCoefficients[i2]);
	Max = i1 + i2;
	i3 = Max;
	if ((i3 & 1) == 0)
	  i3 >>= 1;
	else
	  {
	    i3 >>= 1;
	    i3 += 1;
	  }
	if (Max >= nbrStates)
	  Max = nbrStates - 1;
	for (; i3 <= Max; ++i3)
	  {
	    i4 =  i1 + i2 - i3;
	    Tmp2 = Tmp;
	    Tmp2.ConjugateProduct(waveFunctionCoefficients[i3]);
	    Phi4 += (Tmp2.Re * waveFunctionCoefficients[i4].Re + Tmp2.Im * waveFunctionCoefficients[i4].Im) * TmpArray[i3];
//	    cout << i1 << " " << i2 << " " << i3 << " " << i4 << " " << Phi4 << endl;
	  }
      }
  return Phi4;    
}

// evaluate the derivative of the Phi^4 contribution with respect to one of the wave function coefficients
//
// waveFunctionCoefficients = array containing the wave function coefficients
// index = index of the contribution with respect to one of the wave function coefficients for which derivative has to be evaluated
// numericalCoefficients = bi-dimensional array containing the numerical coefficients of arising for the integration of the product of four monopole harmonic in the LLL
// nbrStates = number of states in the Lowest Landau Level 
// return value = derivative of the Phi^4 contribution

Complex EvaluatePhi4Derivative (Complex* waveFunctionCoefficients, int index, double*** numericalCoefficients, int nbrStates)
{
  Complex DPhi4;
  int i3;
  int i4;
  Complex Tmp;
  Complex Tmp2;
  int Max;
  double* TmpArray;

  for (int i1 = index; i1 < nbrStates; ++i1)
    {
      Tmp2 = 4.0 * waveFunctionCoefficients[i1];
      TmpArray = numericalCoefficients[i1][index];
      Max = i1 + index;
      i3 = Max;
      if ((i3 & 1) == 0)
	i3 >>= 1;
      else
	{
	  i3 >>= 1;
	  i3 += 1;
	}
      if (Max >= nbrStates)
	Max = nbrStates - 1;
      for (; i3 <= Max; ++i3)
	{
	  i4 =  i1 + index - i3;
	  Tmp = (waveFunctionCoefficients[i4] * waveFunctionCoefficients[i3]);
	  DPhi4 += (Tmp.Re * Tmp2) * TmpArray[i3];
	}      
    }
  return DPhi4;    
}

