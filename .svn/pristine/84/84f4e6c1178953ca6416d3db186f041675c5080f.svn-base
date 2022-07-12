#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereLMinusOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("LMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lm", "number of time the L- operator has to be applied", 1);

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type LMinus -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz"])->GetInteger();
  int NbrLMinus = ((SingleIntegerOption*) Manager["nbr-lm"])->GetInteger();
  bool FermionFlag = false;
  if (((SingleStringOption*) Manager["statistics"])->GetString() == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromFileName(((SingleStringOption*) Manager["input-file"])->GetString(), NbrParticles, LzMax, FermionFlag) == false)
    {
      return -1;
    }
  if (((SingleStringOption*) Manager["statistics"])->GetString() != 0)
    if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
      {
	FermionFlag = true;
      }
    else
      if ((strcmp ("fermions", ((SingleStringOption*) Manager["statistics"])->GetString()) == 0))
	{
	  FermionFlag = false;
	}
      else
	{
	  cout << ((SingleStringOption*) Manager["statistics"])->GetString() << " is an undefined statistics" << endl;
	}  
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "error while reading " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphere* IntialSpace;
  ParticleOnSphere* TargetSpace;
  if (FermionFlag == true)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
    }
  else
    {
      IntialSpace = new BosonOnSphere(NbrParticles, Lz, LzMax);
    }
  for (int i = 1; i <= NbrLMinus; ++i)
    {
      if (FermionFlag == true)
	{
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz, LzMax - (2 * i), MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles, Lz, LzMax - (2 * i), MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), MemorySpace);
	    }
#endif
	}
      else
	{
	  TargetSpace = new BosonOnSphere(NbrParticles, Lz, LzMax - (2 * i));
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension());
      ParticleOnSphereLMinusOperator LMinus(IntialSpace, Lz, LzMax  - (2 * i) + 2);
      LMinus.Multiply(InitialVector, TargetVector);
      delete IntialSpace;
      IntialSpace = TargetSpace;
      InitialVector = TargetVector;
    }  
  if (InitialVector.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString()) == false)
    {
      cout << "error while writing " << ((SingleStringOption*) Manager["output-file"])->GetString() << endl;
      return -1;
    }
  delete IntialSpace;
  return 0;
}

