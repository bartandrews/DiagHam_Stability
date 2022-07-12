#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Operator/ParticleOnSphereWithSpinLMinusOperator.h"
#include "Operator/ParticleOnSphereWithSpinLPlusOperator.h"
#include "Operator/ParticleOnSphereWithSpinLPlusUpOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
//#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSpinLMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += DataGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "sz", "twice the total sz value of the system for the initial state", 0);
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "use Hilbert-space comprising all sz sectors", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "with all-sz: parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);

  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new SingleStringOption  ('S', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lm", "number of time the L- operator has to be applied", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "lplus", "apply L+ operator instead of L-");
  (*SystemGroup) += new BooleanOption  ('\n', "up-only", "apply operator affecting up-spins only");

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  if (Manager.GetString("input-file")==NULL)
    {
      cout << "An input file is required!"<<endl;
      exit(-1);
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  int TotalSz = Manager.GetInteger("sz");
  int PairParity = Manager.GetInteger("pair-parity");
  int NbrLMinus = Manager.GetInteger("nbr-lm");
  bool FermionFlag = false;
  if (Manager.GetBoolean("all-sz"))
    TotalSz = -1;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, Lz, TotalSz, FermionFlag) == false)
    {
      return -1;
    }
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }  
	}
    }
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }
  if ((Manager.GetBoolean("all-sz")==false)&&((NbrParticles&1) != (TotalSz&1)))
    {
      cout << "Sz and NbrParticles must have the same parity" << endl;
      return -1;
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  unsigned long MemorySpace = 9ul << 20;
  ParticleOnSphereWithSpin* InitialSpace;
  ParticleOnSphereWithSpin* TargetSpace;
  if (FermionFlag == true)
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
	    {
	      InitialSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	    }
	  else
	    {
	      // InitialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
	      exit(-1);
	    }
#else
	  if (LzMax <= 15)
	    {
	      InitialSpace = new FermionOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz, MemorySpace);
	    }
	  else
	    {
	      // InitialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
	      exit(-1);
	    }
#endif
	}
      else
	{
	  InitialSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, Lz, LzMax, MemorySpace);
	}
    }
  else
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
	  InitialSpace = new BosonOnSphereWithSpin(NbrParticles, Lz, LzMax, TotalSz);
	}
      else
	{
	  if (PairParity>=0)
	    InitialSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz, LzMax, PairParity, MemorySpace);
	  else
	    InitialSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz, LzMax, MemorySpace);
	}
    }
  if (Manager.GetBoolean("lplus")==false)
    {
      for (int i = 1; i <= NbrLMinus; ++i)
	{
	  if (FermionFlag == true)
	    {
	      if (Manager.GetBoolean("all-sz")==false)
		{
#ifdef __64_BITS__
		  if (LzMax <= 31)
		    {
		      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz - (2 * i), LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
		      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
		      exit(-1);
		    }
#else
		  if (LzMax <= 15)
		    {
		      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz - (2 * i), LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax - (2 * i), TotalSz, MemorySpace);
		      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
		      exit(-1);
		    }
#endif
		}
	      else
		{
		  TargetSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("all-sz")==false)
		TargetSpace = new BosonOnSphereWithSpin(NbrParticles, Lz - (2 * i), LzMax, TotalSz, MemorySpace);
	      else
		{
		  if (PairParity>=0)
		    TargetSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz - (2 * i), LzMax, PairParity, MemorySpace);
		  else
		    TargetSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz - (2 * i), LzMax, MemorySpace);
		}
	    }
	  InitialSpace->SetTargetSpace(TargetSpace);
	  TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension());
	  if (TargetSpace->GetHilbertSpaceDimension()!=InitialSpace->GetTargetHilbertSpaceDimension())
	    {
	      cout << "Problem with setting target space"<<endl;
	      exit(-1);
	    }
	  ParticleOnSphereWithSpinLMinusOperator LMinus(InitialSpace, Lz - (2 * i) + 2, LzMax);
	  LMinus.Multiply(InitialVector, TargetVector);
	  delete InitialSpace;
	  InitialSpace = TargetSpace;
	  InitialVector = TargetVector;
	  InitialVector/=InitialVector.Norm();
	}
    }
  else // LPlus - mode
    {
      for (int i = 1; i <= NbrLMinus; ++i)
	{
	  if (FermionFlag == true)
	    {
	      if (Manager.GetBoolean("all-sz")==false)
		{
		  
#ifdef __64_BITS__
		  if (LzMax <= 31)
		    {
		      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz + (2 * i), LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax + (2 * i), TotalSz, MemorySpace);
		      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
		      exit(-1);
		    }
#else
		  if (LzMax <= 15)
		    {
		      TargetSpace = new FermionOnSphereWithSpin(NbrParticles, Lz + (2 * i), LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      // TargetSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax + (2 * i), TotalSz, MemorySpace);
		      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
		      exit(-1);
		    }
#endif
		}
	      else
		{
		  TargetSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("all-sz")==false)
		TargetSpace = new BosonOnSphereWithSpin(NbrParticles, Lz + (2 * i), LzMax, TotalSz);
	      else
		{
		  if (PairParity>=0)
		    TargetSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz + (2 * i), LzMax, PairParity, MemorySpace);
		  else
		    TargetSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, Lz + (2 * i), LzMax, MemorySpace);
		}
	    }
	  InitialSpace->SetTargetSpace(TargetSpace);
	  TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension());
	  if (TargetSpace->GetHilbertSpaceDimension()!=InitialSpace->GetTargetHilbertSpaceDimension())
	    {
	      cout << "Problem with setting target space"<<endl;
	      exit(-1);
	    }
	  if (Manager.GetBoolean("up-only"))
	    {
	      ParticleOnSphereWithSpinLPlusUpOperator LPlusUp(InitialSpace, Lz + (2 * i) - 2, LzMax);
	      LPlusUp.Multiply(InitialVector, TargetVector);
	    }
	  else
	    {
	      ParticleOnSphereWithSpinLPlusOperator LPlus(InitialSpace, Lz + (2 * i) - 2, LzMax);
	      LPlus.Multiply(InitialVector, TargetVector);
	    }
	  delete InitialSpace;
	  InitialSpace = TargetSpace;
	  InitialVector = TargetVector;
	  double Norm = InitialVector.Norm();
	  cout << "Amplitude before normalization: "<<Norm<<endl;
	  if (Norm >0.0)
	    InitialVector/=Norm;
	  else
	    {
	      cout << "Cannot normalize target vector - zero vector found."<<endl;
	    }
	}
    }
  char *OutputName;
  if (Manager.GetString("output-file")!=NULL)
    {
      OutputName = new char[strlen(Manager.GetString("output-file"))+1];
      strcpy(OutputName,Manager.GetString("output-file"));
    }
  else
    {
      OutputName = new char[strlen(Manager.GetString("input-file"))+10];
      sprintf(OutputName,"%s_L-",Manager.GetString("input-file"));
    }
  if (InitialVector.WriteVector(OutputName) == false)
    {
      cout << "error while writing " << OutputName << endl;
      return -1;
    }
  delete InitialSpace;
  delete [] OutputName;
  return 0;
}

