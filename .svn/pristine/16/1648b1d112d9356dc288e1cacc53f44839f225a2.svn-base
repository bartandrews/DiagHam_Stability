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

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Operator/ParticleOnSphereWithSpinSzParityOperator.h"
#include "Operator/ParticleOnSphereSpinOperator.h"

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
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

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
  OptionManager Manager ("FQHESphereWithSpinLValue" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "show-all", "show all S^2_i average values");
  (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "no-spin", "do not compute the S^2 value of the state");
  (*SystemGroup) += new BooleanOption  ('\n', "no-szparity", "do not compute the parity under the Sz<->-Sz symmetry");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
//  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
//  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
//  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "show-extracted", "show values extracted from file name");
//  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSpinLValue -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereLValue -h" << endl;
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


  unsigned long MemorySpace = 9ul << 20;
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
	      Space = new BosonOnSphereWithSpin (NbrParticles, TotalLz, LzMax, TotalSz);
	    }
	  else
	    {
	      int LzSymmetry = 0;
	      int SzSymmetry = 0;
	      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, LzMax, TotalLz, TotalSz, SzSymmetry, LzSymmetry, FermionFlag) == false)
		{
		  cout << "error while retrieving system parameters from file name " << StateFileName << endl;
		  return -1;
		}
	      if (LzSymmetry == 0)
		{
		  if (SzSymmetry == 0)
		    {		      
		      Space = new BosonOnSphereWithSU2Spin (NbrParticles, TotalLz, LzMax, TotalSz);
		    }
		  else
		    {		      
		      Space = new BosonOnSphereWithSU2SpinSzSymmetry (NbrParticles, TotalLz, LzMax, TotalSz, (SzSymmetry == -1));
		    }		  
		}
	      else
		{
		  if (SzSymmetry == 0)
		    {		      
		      Space = new BosonOnSphereWithSU2SpinLzSymmetry (NbrParticles, LzMax, TotalSz, (LzSymmetry == -1));
		    }
		  else
		    {		      
		      Space = new BosonOnSphereWithSU2SpinLzSzSymmetry (NbrParticles, LzMax, TotalSz, (SzSymmetry == -1), (LzSymmetry == -1));
		    }
		}
	    }
	}
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  ParticleOnSphereWithSpinL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0);
  RealVector TmpState(Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &State, &TmpState);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  double L2Value = TmpState * State;
  double RawTmpAngularMomentum = 0.5 * (sqrt (((double)4.0 * L2Value) + (double)1.0) - 1.0);
  cout << "<L^2> = " << L2Value << endl
       << "<L> = " << RawTmpAngularMomentum << endl;
  if (Manager.GetBoolean("no-spin") == false)
    {
      bool fixedSz = !(Manager.GetBoolean("all-sz"));
      if(!(Manager.GetBoolean("show-all")))
      {
      ParticleOnSphereWithSpinS2Hamiltonian Hamiltonian2 (Space, NbrParticles, LzMax, TotalLz, TotalSz, Architecture.GetArchitecture(), 1.0, 0, fixedSz);
//	ParticleOnSphereWithSpinS2Hamiltonian Hamiltonian2 (Space, NbrParticles, LzMax, TotalLz, TotalSz, Architecture.GetArchitecture(), 1.0, -1, false, 0, fixedSz);
      VectorHamiltonianMultiplyOperation Operation2 (&Hamiltonian2, &State, &TmpState);
      Operation2.ApplyOperation(Architecture.GetArchitecture());
      L2Value = TmpState * State;
      RawTmpAngularMomentum = 0.5 * (sqrt (((double)4.0 * L2Value) + (double) 1.0) - 1.0);
      cout << "<S^2> = " << L2Value << endl
	   << "<S> = " << RawTmpAngularMomentum << endl;
      }
    }
//  if ((Manager.GetBoolean("all-sz"))||((Manager.GetBoolean("no-szparity") == false)&&(TotalSz==0)))
  if ( (Manager.GetBoolean("no-szparity") == false) && (TotalSz==0) && (!Manager.GetBoolean("all-sz")) )
    {
      Complex Tmp;
      if (SzSymmetrizedBasis == false)
	{
	  ParticleOnSphereWithSpinSzParityOperator Operator(Space);
	  Tmp = Operator.MatrixElement(State, State);
	}
      else
	{
	  if (SzMinusParity == false)
	    {
	      Tmp = 1.0;
	    }
	  else
	    {
	      Tmp = -1.0;
	    }
	}
      cout  << "<P_sz> = " << Tmp.Re << endl;      
    }
  if (Manager.GetBoolean("all-sz"))
    {
      Complex Tmp;
            
      cout << "====================================" << endl;

      ParticleOnSphereSpinOperator SxOperator(Space, 0, LzMax);
      Tmp = SxOperator.MatrixElement(State,State);
      cout << "<S_x> = " << Tmp.Re << endl;

      ParticleOnSphereSpinOperator SyOperator(Space, 1, LzMax);
      Tmp = SyOperator.MatrixElement(State,State);
      cout << "Im(<S_y>) = " << Tmp.Im << endl;
      
      ParticleOnSphereSpinOperator SzOperator(Space, 2, LzMax);
      Tmp = SzOperator.MatrixElement(State,State);
      cout << "<S_z> = " << Tmp.Re << endl;

      cout << "====================================" << endl;
      
      if(Manager.GetBoolean("show-all"))
      {
	Complex Sx2 = SxOperator.PartialMatrixElementSquare(State,State,0,Space->GetHilbertSpaceDimension());
	cout << "<S_x^2> = " << Sx2.Re << endl;
	
	Complex Sy2 = SyOperator.PartialMatrixElementSquare(State,State,0,Space->GetHilbertSpaceDimension());
	cout << "<S_y^2> = " << Sy2.Re << endl;
	
	Complex Sz2 = SzOperator.PartialMatrixElementSquare(State,State,0,Space->GetHilbertSpaceDimension());
	cout << "<S_z^2> = " << Sz2.Re << endl;
	
	cout << "====================================" << endl;

	double S2=Sx2.Re+Sy2.Re+Sz2.Re;
	cout << "<S^2> = " << S2 << endl;
	double RawS = 0.5 * (sqrt (((double)4.0 * S2) + (double) 1.0) - 1.0);
	cout << "<S> = " << RawS << endl;
      }
      
    }


  delete Space;
  return 0;
}

