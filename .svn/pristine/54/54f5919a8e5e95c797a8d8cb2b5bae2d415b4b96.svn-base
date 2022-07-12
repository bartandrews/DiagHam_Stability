#include "config.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/ParticleOnSphereWithSpinDensityDensityOperator.h"
#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSzLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"

#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"

#include <cstring>
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
  OptionManager Manager ("FQHESphereWithSpinCorrelation" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('b', "bilayer", "adjust normalization as a bilayer correlation function");
  (*SystemGroup) += new BooleanOption  ('\n', "shift", "calculate 'shift' as defined for bilayer states");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
//  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
//  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
  (*SystemGroup) += new BooleanOption  ('\n', "show-extracted", "show values extracted from file name");
//  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereCorrelation -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int PairParity = Manager.GetInteger("pair-parity");
  bool SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
  bool SzMinusParity = Manager.GetBoolean("minus-szparity");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool LzMinusParity = Manager.GetBoolean("minus-lzparity");
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  bool DensityFlag = Manager.GetBoolean("density");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  int TmpTotalSz=TotalSz;
  if (Manager.GetBoolean("all-sz"))
    TmpTotalSz=-1;
  if (NbrParticles==0)
    {
      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, TmpTotalSz, SzSymmetrizedBasis, SzMinusParity, 
							       LzSymmetrizedBasis, LzMinusParity, FermionFlag) == false)
	{
	  return -1;
	}
      else
	{
	  if (!Manager.GetBoolean("all-sz"))
	    TotalSz=TmpTotalSz;
	  if (Manager.GetBoolean("show-extracted") == true)
	    {
	      cout << "N=" << NbrParticles << "  LzMax=" << LzMax << "  TotalLz=" << TotalLz << "  TotalSz=" << TotalSz;
	      if (LzSymmetrizedBasis == true)
		{
		  cout << "  Lz symmetrized basis ";
		  if (LzMinusParity == true)
		    cout << "(minus parity) ";
		  else
		    cout << "(plus parity) ";
		}
	      if (SzSymmetrizedBasis == true)
		{
		  cout << "  Sz symmetrized basis ";
		  if (SzMinusParity == true)
		    cout << "(minus parity) ";
		  else
		    cout << "(plus parity) ";
		}
	      cout << endl;
	    }
	}
    }
  if (Manager.GetBoolean("lzsymmetrized-basis") == true)
    {
      LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
      LzMinusParity = Manager.GetBoolean("minus-lzparity");      
    }
  if (Manager.GetBoolean("szsymmetrized-basis") == true)
    {
      SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
      SzMinusParity = Manager.GetBoolean("minus-szparity");
    }
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("bosons", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }  
	}
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* StateFileName = Manager.GetString("state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector State;
  if (State.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }


  unsigned long MemorySpace = 9l << 20;
  ParticleOnSphereWithSpin* Space;
  if (FermionFlag == true)
    {
      if (Manager.GetBoolean("all-sz"))
	{
	  if (LzSymmetrizedBasis == false)
	    Space = new FermionOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, MemorySpace);
	  else
	    Space = new FermionOnSphereWithSpinAllSzLzSymmetry (NbrParticles, LzMax, LzMinusParity, MemorySpace);
	}
      else if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
	  {
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (LzMax <= 63)
#else
		  if (LzMax <= 31)
#endif
		    {
		      Space = new FermionOnSphereWithSpinLong(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return -1;
		    }	
	      }
	  }
      else
	{
#ifdef __128_BIT_LONGLONG__
	  if (LzMax >= 61)
#else
	    if (LzMax >= 29)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return -1;
	      }	
	  if (SzSymmetrizedBasis == true) 
	    if (LzSymmetrizedBasis == false)
	      {
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
		    }
		  else
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrParticles, TotalLz, LzMax, SzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);
		    }
		  }
	    else
#ifdef __64_BITS__
	      if (LzMax <= 28)
#else
		if (LzMax <= 13)
#endif
		  {
		    if (Manager.GetString("load-hilbert") == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, LzMax, SzMinusParity,
									LzMinusParity, MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);
		  }
		else
		  {
		    if (Manager.GetString("load-hilbert") == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrParticles, LzMax, SzMinusParity,
									    LzMinusParity, MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);
		    
		  }
	      else
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetry(Manager.GetString("load-hilbert"), MemorySpace);	      
		    }
		  else
		    {
		      if (Manager.GetString("load-hilbert") == 0)
			Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrParticles, LzMax, TotalSz, LzMinusParity, MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetryLong(Manager.GetString("load-hilbert"), MemorySpace);	      
		    }
	}
    }
  else
    {
      if (Manager.GetBoolean("all-sz"))
	{
	  if ( PairParity >=0 ) 
	    Space = new BosonOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, PairParity, MemorySpace);
	  else
	    Space = new BosonOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, MemorySpace);
	}
      else
	{
	  if (Manager.GetBoolean("use-alt") == false)
	    {
	      Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
	    }
	  else
	    {
	      Space = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, LzMax, TotalSz);
	    }
	}
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

    cout << "dim = " << Space->GetHilbertSpaceDimension() << endl;

  ParticleOnSphereFunctionBasis Basis(LzMax);

  Complex* Sum = new Complex [3];
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);
  Complex** PrecalculatedValues = new Complex* [3];
  for (int i = 0; i < 3; ++i)
    PrecalculatedValues[i] = new Complex [LzMax + 1];
   
  if (DensityFlag == true)
    {
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityOperator Operator (Space, i, 0, i, 0);
	  PrecalculatedValues[0][i] =   Operator.MatrixElement(State, State);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityOperator Operator (Space, i, 1, i, 1);
	  PrecalculatedValues[1][i] =   Operator.MatrixElement(State, State);
	}
    }
  else
    {
      Basis.GetFunctionValue(Value, TmpValue, LzMax);
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, i, 0, LzMax, 0, i, 0, LzMax, 0);
	  PrecalculatedValues[0][i] =   Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, i, 1, LzMax, 1, i, 1, LzMax, 1);
	  PrecalculatedValues[2][i] =  Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
      for (int i = 0; i <= LzMax; ++i)
	{
	  ParticleOnSphereWithSpinDensityDensityOperator Operator1 (Space, i, 0, LzMax, 1, i, 0, LzMax, 1);
	  PrecalculatedValues[1][i] =  Operator1.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	}
    }

  ofstream File;
  File.precision(14);
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("state") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  
  double Factor2 = sqrt (0.5 * LzMax );
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true) 
    Factor2 = 1.0;
  if (DensityFlag == true)
    {
      double Factor1 = 1.0;
      File << "# density coefficients for " << Manager.GetString("state") << endl;
      File << "#" << endl << "# (l+S) n_l^{u} n_l^{d}" << endl;
      double Sum2 = 0.0;
      for (int i = 0; i <= LzMax; ++i)
	{
	  File << "# " << i;
	  for (int j = 0; j < 2; ++j)
	    {
	      File << " " << PrecalculatedValues[j][i].Re;
	      Sum2 += PrecalculatedValues[j][i].Re;
	    }
	  File << endl;
	}
      File << "# sum = " << Sum2 << endl;
      File << "# dist (rad) rho_{u} rho_{d} rho_{u}+rho_{d}" << endl;
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  for (int j = 0; j < 3; ++j)
	    Sum[j] = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, i);
	      Sum[0] += PrecalculatedValues[0][i] * (Conj(TmpValue) * TmpValue);
	      Sum[1] += PrecalculatedValues[1][i] * (Conj(TmpValue) * TmpValue);
	    }
	  File << (X * Factor2);
	  File << " " << (Sum[0].Re  * Factor1) << " " << (Sum[1].Re  * Factor1) << " " << ((Sum[0].Re + Sum[1].Re) * Factor1) << endl;
	  X += XInc;
	}
      File.close();     
      return 0;
    }
  if (Manager.GetBoolean("coefficients-only"))
    {
      File << "# density-density correlation coefficients for " << Manager.GetString("state") << endl;
      File << "#" << endl << "# (l+S) n_l^{u,u} n_l^{u,d} n_l^{d,d}" << endl;
      for (int i = 0; i <= LzMax; ++i)
	{
	  File << i;
	  for (int j=0; j<3; ++j)
	    File << " " << PrecalculatedValues[j][i].Re;
	  File << endl;
	}
      File.close();
      for (int i = 0; i < 3; ++i)
	delete[] PrecalculatedValues[i];
      delete[] PrecalculatedValues;      
    }
  else
    {
      double Factor1;
      if (Manager.GetBoolean("bilayer"))
	Factor1 = (64.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      else
	Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      File << "# dist (rad) rho_{u,u} rho_{u,d} rho_{d,d}"<<endl;
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  for (int j = 0; j < 3; ++j)
	    Sum[j] = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, i);
	      for (int j = 0; j < 3; ++j)	    
		Sum[j] += PrecalculatedValues[j][i] * (Conj(TmpValue) * TmpValue);
	    }
	  File << (X * Factor2);
	  for (int j = 0; j < 3; ++j)
	    File << " " << (Norm(Sum[j])  * Factor1);
	  File << endl;
	  X += XInc;
	}
      File.close();
      for (int i = 0; i < 3; ++i)
	delete[] PrecalculatedValues[i];
      delete[] PrecalculatedValues;
      delete[] Sum;

      // CALCULATION OF "SHIFT OPERATOR":
      if (Manager.GetBoolean("shift"))
	{
	  char* OutputNameCorr = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "shift.dat");
	  Complex* PrecalculatedValues2 = new Complex [2*(LzMax + 1)];
	  Complex TmpValue2;      
	  int Pos=0;      
	  for (int m = 0; m <= LzMax; ++m)
	    {	    
	      ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, m, 0, LzMax, 1, m, 0, LzMax, 1);
	      PrecalculatedValues2[Pos] = Operator.MatrixElement(State, State);
	      ++Pos;
	      ParticleOnSphereWithSpinDensityDensityOperator Operator2 (Space, LzMax, 0, m, 1, m, 1, LzMax, 0);
	      PrecalculatedValues2[Pos] = Operator2.MatrixElement(State, State);
	      ++Pos;
	    }

	  ofstream File;
	  File.precision(14);
	  File.open(OutputNameCorr, ios::binary | ios::out);
	  X=0.0;
	  for (int x = 0; x < NbrPoints; ++x)
	    {
	      Value[0] = X;
	      Value[1] = 0.0;
	      Pos = 0;
	      Complex Sum3 = 0.0;
	      Sum2 = 0.0;	  
	      for (int m = 0; m <= LzMax; ++m)
		{
		  Basis.GetFunctionValue(Value, TmpValue, m);
		  double CommonFactor=SqrNorm(TmpValue);      
		  Sum3 += CommonFactor*PrecalculatedValues2[Pos];
		  ++Pos;
		  Sum2 += CommonFactor*PrecalculatedValues2[Pos];
		  ++Pos;
		}
	      if (Sum2!=0.0)
		File << (X * Factor2) << " " << Real(Sum3/Sum2) << " " << Imag(Sum3/Sum2) << endl;
	      else
		File << (X * Factor2) << " nan" << endl;
	      if (x==NbrPoints-1)
		cout << "Shift = " << Real(Sum3/Sum2) << endl; 
	      X += XInc;
	    }
	  File.close();
	  delete[] PrecalculatedValues;
	  delete [] OutputNameCorr;
	}
    }  
  return 0;
}

