#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "MathTools/BinomialCoefficients.h"

#include "Polynomial/IntegerPolynomial.h"
#include "Polynomial/QDeformedBinomialCoefficients.h"

#include "Tools/FQHESpectrum/MomentumMultipletSet.h"


#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>

using std::cout;
using std::endl;


// get the total number of states corresponding to a given parafermionic quasihole system
// 
// nbrParticles = number of parafermions
// nbrQuasiholes = number of qusaiholes (must be a multiple of kType)
// kType = number of parafermion species
// binomialCoefficients = array that contains all usefull binomial coefficients
// return value = number of states
long GetNumberOfStates (int nbrParticles, int nbrQuasiholes, int kType, BinomialCoefficients& binomialCoefficients);

// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoeffients = reference on the binomial coefficients
// return value = number of way
long GetParafermionPartitionNumber (int nbrParafermions, int nbrBoxes, int kType, BinomialCoefficients& binomialCoeffients);

// get the number of ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// return value = number of way
long GetFixedLzFreeNumberBosonPartitionNumber (int momentum, int nbrStates);

// get all ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// permutations = reference on the array where permutations will be stored
// partitionNumber = reference on the integer where the corresponding partition will be stored
void GetFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, long& partitionNumber);

// recursive function associated to GetFixedLzFreeNumberBosonPermutations
// 
// momentum = total momemtum
// nbrStates = number of states
// permutations = reference on the array where permutations will be stored
// position = current position in the array
// return value = new current position in the array
int RecursiveFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, int position);

// get the multiplet decomposition corresponding to a given parafermionic quasihole system
// 
// nbrParticles = number of parafermions
// nbrQuasiholes = number of qusaiholes (must be a multiple of kType)
// kType = number of parafermion species
// binomialCoefficients = array that contains all usefull q-deformed binomial coefficients
// return value = corresponding set of multiplets
MomentumMultipletSet GetMomentumMultipletDecomposition (int nbrParticles, int nbrQuasiholes, int kType, QDeformedBinomialCoefficients& binomialCoefficients);


int main(int argc, char** argv)
{
  OptionManager Manager ("ParafermionQuasiholeDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "number of particles", 2);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('q', "nbr-quasiholes", "number of quasiholes", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "total-only", "only display the total number of states");
  (*SystemGroup) += new BooleanOption  ('\n', "lz-values", "display lz degeneracy instead of l degeneracy");
  (*SystemGroup) += new SingleStringOption  ('\n', "display-style", "set data output display style (raw, pretty, latex)", "raw");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ParafermionQuasiholeDimension -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int KValue = ((SingleIntegerOption*) Manager["k-value"])->GetInteger(); 
  int NbrParticles  = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrQuasiholes = ((SingleIntegerOption*) Manager["nbr-quasiholes"])->GetInteger(); 
  bool LzFlag = ((BooleanOption*) Manager["lz-values"])->GetBoolean();

  if ((NbrParticles % KValue) != 0)
    {
      cout << "the number of particles must be an multiple of k" << endl;
      return 1;
    }
  if ((NbrQuasiholes % KValue) != 0)
    {
      cout << "the number of quasiholes must be an multiple of k" << endl;
      return 1;
    }
  
  int MaxBinomial = 10;
  BinomialCoefficients TmpBinomialCoefficients(MaxBinomial);
  QDeformedBinomialCoefficients QDeformed(MaxBinomial);
  
  if (((BooleanOption*) Manager["total-only"])->GetBoolean() == false)
    {
      MomentumMultipletSet Multiplets(GetMomentumMultipletDecomposition(NbrParticles, NbrQuasiholes, KValue, QDeformed));
      if (strncmp("pretty", ((SingleStringOption*) Manager["display-style"])->GetString(), 6) == 0)
	if (LzFlag == false)
	  cout << Multiplets << endl;
	else
	  {
	    Multiplets.PrintLz(cout) << endl;
	  }
      else
	if (strncmp("latex", ((SingleStringOption*) Manager["display-style"])->GetString(), 5) == 0)
	  {
	    if (LzFlag == false)
	      {
		cout << "$" << Multiplets.GetNbrStates() << "$  & ";
		for (int i = (Multiplets.GetMaximumMomentum() & 1); i < Multiplets.GetMaximumMomentum(); i += 2)
		  cout << "$" << Multiplets[i] << "$  & ";
		cout << "$" << Multiplets[Multiplets.GetMaximumMomentum()] << "$" << endl;
	      }
	    else
	      {
		cout << "$" << Multiplets.GetNbrStates() << "$  & ";
		for (int i = (Multiplets.GetMaximumMomentum() & 1); i < Multiplets.GetMaximumMomentum(); i += 2)
		  {
		    int Value = 0;
		    for (int j = i; j <=  Multiplets.GetMaximumMomentum(); j += 2)
		      Value += Multiplets[j];
		    cout << "$" << Value << "$ & ";
		  }
		cout << "$" << Multiplets[Multiplets.GetMaximumMomentum()] << "$" << endl;
	      }
	  }       
	else
	  {
	    if (LzFlag == false)
	      {
		for (int i = (Multiplets.GetMaximumMomentum() & 1); i <= Multiplets.GetMaximumMomentum(); i += 2)
		  cout << Multiplets[i] << " ";
		cout << endl;
	      }
	    else
	      {
		for (int i = (Multiplets.GetMaximumMomentum() & 1); i <= Multiplets.GetMaximumMomentum(); i += 2)
		  {
		    int Value = 0;
		    for (int j = i; j <=  Multiplets.GetMaximumMomentum(); j += 2)
		      Value += Multiplets[j];
		    cout << Value << " ";
		  }
		cout << endl;
	      }
	  }
      long HilbertSpaceDimension = GetNumberOfStates(NbrParticles, NbrQuasiholes, KValue, TmpBinomialCoefficients);
      if (HilbertSpaceDimension == Multiplets.GetNbrStates())
	cout << HilbertSpaceDimension << endl;      
      else
	{
	  cout << "an error occured while evaluating momentum multiplet decomposition producing different number of states (" << Multiplets.GetNbrStates() 
	       << " and " << HilbertSpaceDimension << ")" << endl;
	  return -1;
	}
    }
  else
    {
      long HilbertSpaceDimension = GetNumberOfStates(NbrParticles, NbrQuasiholes, KValue, TmpBinomialCoefficients);
      cout << HilbertSpaceDimension << endl;      
    }
}



// get the total number of states corresponding to a given parafermionic quasihole system
// 
// nbrParticles = number of parafermions
// nbrQuasiholes = number of qusaiholes (must be a multiple of kType)
// kType = number of parafermion species
// binomialCoeffients = array that contains all usefull binomial coefficients
// return value = number of states

long GetNumberOfStates (int nbrParticles, int nbrQuasiholes, int kType, BinomialCoefficients& binomialCoefficients)
{
  long HilbertSpaceDimension = 0;
  int NbrUnclusteredParafermions = 0;
  while (NbrUnclusteredParafermions <= nbrParticles)
    {
      HilbertSpaceDimension += (GetParafermionPartitionNumber(NbrUnclusteredParafermions, nbrQuasiholes / kType, kType, binomialCoefficients) * 
				binomialCoefficients(nbrQuasiholes + ((nbrParticles - NbrUnclusteredParafermions) / kType), nbrQuasiholes));
      NbrUnclusteredParafermions += kType;
    }
  return HilbertSpaceDimension;
}


// get the number of ways to put F parafermions in n box
// 
// nbrParafermions = number of parafermions
// kType = number of parafermion species
// nbrBoxes = number of boxes
// binomialCoefficients = reference on the binomial coefficients
// return value = number of way

long GetParafermionPartitionNumber (int nbrParafermions, int nbrBoxes, int kType, BinomialCoefficients& binomialCoefficients)
{
  if ((nbrParafermions == 0) || (kType == 1))
    return 1l;
  if (kType == 2)
    {
      if (nbrParafermions <= nbrBoxes)
        return binomialCoefficients(nbrBoxes, nbrParafermions);
      else
	return 0l;
    }
  int** Permutations;
  long NbrPermutations;
  GetFixedLzFreeNumberBosonPermutations(nbrParafermions, kType - 1, Permutations, NbrPermutations);
  long NbrWays = 0l;
  long Tmp = 1l;
  int* TmpPermutation;
  int TmpCoefficient;
  for (int i = 0; i < NbrPermutations; ++i)
    {
      TmpPermutation = Permutations[i];
      Tmp = 1l;
      for (int j = 1; ((j < kType) && (Tmp != 0l)); ++j)
	{
	  int k = 1;
	  TmpCoefficient = (j * nbrBoxes * kType) + ((kType - (2 * (kType - j) * j)) * TmpPermutation[j - 1]);
	  for (; k < j; ++k)
	    TmpCoefficient -= (2 * (kType - j) * k) * TmpPermutation[k - 1];
	  ++k;
	  for (; k < kType; ++k)
	    TmpCoefficient -= (2 * (kType - k) * j) * TmpPermutation[k - 1];
	  if ((TmpCoefficient >= 0) && ((TmpCoefficient % kType) == 0) && ((TmpPermutation[j - 1] * kType) <= TmpCoefficient))
	    {
	      Tmp *= binomialCoefficients(TmpCoefficient / kType, TmpPermutation[j - 1]);
	    }
	  else
	    Tmp = 0l;
	}
      NbrWays += Tmp;
    }
  for (int i = 0; i < NbrPermutations; ++i)
    delete[] Permutations[i];
  delete[] Permutations;
  return NbrWays;
}

// get the number of ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// return value = number of way

long GetFixedLzFreeNumberBosonPartitionNumber (int momentum, int nbrStates)
{
  if ((nbrStates == 1) || (momentum == 0))
    return 1l;
  int Max = momentum / nbrStates;
  long Tmp = 0l;
  for (int i = 0; i <= Max; ++i)
    Tmp += GetFixedLzFreeNumberBosonPartitionNumber(momentum - (i * nbrStates), nbrStates - 1);
  return Tmp;
}

// get all ways to write sum_i=1^n i a_i = F where F and n (aka number of possibilities to have a state of total momentum F with a 
// free number of bosons and n one-body state carrying each a Lz momentum ranging from 1 to n)
// 
// momentum = total momemtum (aka F)
// nbrStates = number of states (aka n)
// permutations = reference on the array where permutations will be stored
// partitionNumber = reference on the integer where the corresponding partition will be stored

void GetFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, long& partitionNumber)
{
  partitionNumber = GetFixedLzFreeNumberBosonPartitionNumber(momentum, nbrStates);
  permutations = new int* [partitionNumber];
  for (int i = 0; i < partitionNumber; ++i)
    {
      permutations[i] = new int [nbrStates];
      for (int j = 0; j < nbrStates; ++j)
	permutations[i][j] = 0;
    }
  RecursiveFixedLzFreeNumberBosonPermutations(momentum, nbrStates, permutations, 0);
}

// recursive function associated to GetFixedLzFreeNumberBosonPermutations
// 
// momentum = total momemtum
// nbrStates = number of states
// permutations = reference on the array where permutations will be stored
// position = current position in the array
// return value = new current position in the array

int RecursiveFixedLzFreeNumberBosonPermutations (int momentum, int nbrStates, int**& permutations, int position)
{
  if (nbrStates == 1)
    {
      permutations[position][nbrStates - 1] = momentum;
      return position + 1;
    }
  if (momentum == 0)
    {
      return position + 1;
    }    
  int TmpPosition;
  int TmpMomentum = momentum;
  int Tmp = 0;
  while (TmpMomentum >= 0)
    {
      TmpPosition = RecursiveFixedLzFreeNumberBosonPermutations(TmpMomentum, nbrStates - 1, permutations, position);
      for (; position < TmpPosition; ++position)
	permutations[position][nbrStates - 1] = Tmp;
      TmpMomentum -= nbrStates;
      ++Tmp;
    }
  return position;
}

// get the multiplet decomposition corresponding to a given parafermionic quasihole system
// 
// nbrParticles = number of parafermions
// nbrQuasiholes = number of qusaiholes (must be a multiple of kType)
// kType = number of parafermion species
// binomialCoefficients = array that contains all usefull q-deformed binomial coefficients
// return value = corresponding set of multiplets

MomentumMultipletSet GetMomentumMultipletDecomposition (int nbrParticles, int nbrQuasiholes, int kType, QDeformedBinomialCoefficients& binomialCoefficients)
{
  int NbrUnclusteredParafermions = 0;
  int* TmpPermutation;
  int TmpCoefficient;
  int** Permutations;
  long NbrPermutations;
  long Constant = 1l;
  MomentumMultipletSet TotalMultiplet;
  while (NbrUnclusteredParafermions <= nbrParticles)
    {
      GetFixedLzFreeNumberBosonPermutations(NbrUnclusteredParafermions, kType - 1, Permutations, NbrPermutations);
      IntegerPolynomial TotalP;
      for (int i = 0; i < NbrPermutations; ++i)
	{
	  TmpPermutation = Permutations[i];
	  IntegerPolynomial Tmp (0, &Constant, false);
	  bool Flag = true;
	  int PowerShift = 0;
	  for (int j = 1; j < kType; ++ j)
	    {
	      int k = 1;
	      int Tmp = (kType - j) * j * TmpPermutation[j - 1];
	      for (; k < j; ++k)
		Tmp += ((kType - j) * k) * TmpPermutation[k - 1];
	      ++k;
	      for (; k < kType; ++k)
		Tmp += ((kType - k) * j) * TmpPermutation[k - 1];	      
	      PowerShift += Tmp * TmpPermutation[j - 1];
	    }
	  PowerShift /= kType;
	  for (int j = 1; ((j < kType) && (Flag == true)); ++j)
	    {
	      int k = 1;
	      TmpCoefficient = (j * nbrQuasiholes) + ((kType - (2 * (kType - j) * j)) * TmpPermutation[j - 1]);
	      for (; k < j; ++k)
		TmpCoefficient -= (2 * (kType - j) * k) * TmpPermutation[k - 1];
	      ++k;
	      for (; k < kType; ++k)
		TmpCoefficient -= (2 * (kType - k) * j) * TmpPermutation[k - 1];
	      if ((TmpCoefficient >= 0) && ((TmpCoefficient % kType) == 0) && ((TmpPermutation[j - 1] * kType) <= TmpCoefficient))
		{
		  Tmp *= binomialCoefficients(TmpCoefficient / kType, TmpPermutation[j - 1]);
		}
	      else
		Flag = false;
	    }
	  if (Flag == true)
	    {
	      Tmp.ShiftPowers(PowerShift);
	      TotalP += Tmp;
	    }
	}
      if (TotalP.Defined())
	{
	  MomentumMultipletSet Multiplet1(TotalP);
	  MomentumMultipletSet Multiplet2((nbrQuasiholes * (nbrParticles - NbrUnclusteredParafermions)) / kType);
	  Multiplet2.FindMultipletsForBosons(nbrQuasiholes, (nbrParticles - NbrUnclusteredParafermions) / kType);
	  TotalMultiplet.FuseAndAdd(Multiplet1, Multiplet2);
	}
      for (int i = 0; i < NbrPermutations; ++i)
	delete[] Permutations[i];
      delete[] Permutations;
      NbrUnclusteredParafermions += kType;
    }
  return TotalMultiplet;
}
