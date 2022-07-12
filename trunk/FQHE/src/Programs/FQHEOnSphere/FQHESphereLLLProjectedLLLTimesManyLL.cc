#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"


#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsTimesSlaterProjectionOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereMultipleMonomialsTimesSlaterProjectionOperation.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

// convert a state such that its first nonzero component is equal to 1
//
// state = reference to the state to convert
// return value = converted state

RealVector PutFirstComponentToOne(RealVector& state);
LongRationalVector PutFirstComponentToOne(LongRationalVector& state);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereLLLProjectedLLLTimesTwoLL" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  
  (*SystemGroup) += new SingleStringOption  ('\n', "fermion", "name of the file corresponding to the fermion state in the Two Lowest Landau Level");
  (*SystemGroup) += new SingleStringOption  ('\n', "boson", "name of the file corresponding to the boson state in the LLL");
  (*SystemGroup) += new SingleStringOption  ('\n', "lll-state", "name of the file corresponding to the state in the LLL");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state (should be the one of the bosonic state)");
  (*SystemGroup) += new BooleanOption ('\n',"projection","the state will be projected into the LLL");
  (*SystemGroup) += new BooleanOption ('\n',"reverse-flux","the fluxes bind to each particle are in the opposite direction than the magnetic field");
  (*SystemGroup) += new BooleanOption  ('\n', "projected-haldane", "use an Haldane basis instead of the usual n-body basis for the projected state");
  (*SystemGroup) += new SingleStringOption  ('\n', "projected-referencefile", "use a file as the definition of the reference state for the projected state (should be the one of the bosonic state)");
  (*SystemGroup) += new BooleanOption ('\n', "resume", "the last calcul will be resumed from its last save step");
  (*SystemGroup) += new BooleanOption  ('\n', "3-ll", "consider particles within three Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "4-ll", "consider particles within four Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetric", "Lz->-Lz symmetry used");
  (*SystemGroup) += new SingleIntegerOption ('\n',"index", "create a temporary vector whose index component is put to 1, others to 0. No need of an input file in this case",-1);
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "energy-value", "value of the total effective cyclotron energy of the states that will be computed", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "1ll-nbrparticles", "constrain the number of particles in the lowest LL", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "2ll-nbrparticles", "constrain the number of particles in the second LL", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "3ll-nbrparticles", "constrain the number of particles in the third LL", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "4ll-nbrparticles", "constrain the number of particles in the fourth LL", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "resuming-index", "index from which the calculation will be resumed in case of interruption in the constain mode", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-component", "the min component", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-components", "the number of component computed", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n',"step", "number of time the Bosonic will be divided in",1);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "express the projected state in the normalized basis");
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int LLLNbrParticles = 0;
  int LLLLzMax = 0;
  int LLLTotalLz = 0;
  int NbrFermion = Manager.GetInteger("nbr-particles");
  int LzMaxFermion = Manager.GetInteger("lzmax");
  int TotalLzFermion = Manager.GetInteger("total-lz");
  bool FermionFlag = false;
  bool ReverseFluxFlag = Manager.GetBoolean("reverse-flux");
  int Step = Manager.GetInteger("step");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool Projection = Manager.GetBoolean("projection");
  bool Resume = Manager.GetBoolean("resume");
  int MinComponent = Manager.GetInteger("min-component");
  int NbrComponents = Manager.GetInteger("nbr-components");
  char* LLLFileName = Manager.GetString("lll-state");
  char* FermionFileName = Manager.GetString("fermion");
  int NbrLL = 2;
  if (Manager.GetBoolean("3-ll") == true)
    NbrLL = 3;
  if (Manager.GetBoolean("4-ll") == true)
    NbrLL = 4;
  cout <<"Number of Lambda Landau level = "<<NbrLL<<endl;
  
  bool Constraint = false;
  bool LLLFermionFlag = true;
  
  int * ParticleNumberConstrains = new int [4];
  ParticleNumberConstrains[0] = Manager.GetInteger("1ll-nbrparticles");
  ParticleNumberConstrains[1] = Manager.GetInteger("2ll-nbrparticles");
  ParticleNumberConstrains[2] = Manager.GetInteger("3ll-nbrparticles");
  ParticleNumberConstrains[3] = Manager.GetInteger("4ll-nbrparticles");
  int ResumingIndex = Manager.GetInteger("resuming-index");
  bool Symmetric = Manager.GetBoolean("symmetric");
  int Index = Manager.GetInteger("index");
  int EnergyValue = Manager.GetInteger("energy-value");
  
  
  if((EnergyValue != -1)||(ParticleNumberConstrains[0] != -1)||(ParticleNumberConstrains[1] != -1)||(ParticleNumberConstrains[2] != -1)||(ParticleNumberConstrains[3] != -1))
    Constraint = true;
  
  if ((Constraint == true) && ((NbrFermion == 0) || (LzMaxFermion == -1)))
    {
      cout <<"In the contraint mode, the particle number and the LzMax value must be given in option. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" <<endl;
      return -1; 
    }
  
  
  if (LLLFileName == 0)
    {
      cout << "error, a lowest Landau level state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }
  
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(LLLFileName, LLLNbrParticles, LLLLzMax, LLLTotalLz, LLLFermionFlag) == false)
    {
      return -1;
    }
  
  if((Index == -1) && (Constraint == false))
    {
      if (FermionFileName == 0)
	{
	  cout << "error, a fermionic state file should be provided. See man page for option syntax or type FQHESphereLLLProjectedLLLTimesManyLL -h" << endl;
	  return -1;
	}
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(FermionFileName, NbrFermion, LzMaxFermion, TotalLzFermion, FermionFlag) == false)
	{
	  return -1;
	}
      
    }
  if (LLLNbrParticles != NbrFermion)
    {
      cout << "The number of particles in the two states must be the same" <<endl;
      return -1;
    }
  
  int Parity = LLLTotalLz & 1;
  if (Parity != ((LLLNbrParticles * LLLLzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;
    }
  
  
  if (IsFile(LLLFileName) == false)
    {
      cout << "state " << LLLFileName << " does not exist or can't be opened" << endl;
      return -1;
    }
  
  
  ParticleOnSphere* LLLSpace = 0;
  
  if (LLLFermionFlag == false)
    {
      if(HaldaneBasisFlag == false)
        {
#ifdef  __64_BITS__
          if((LLLLzMax + LLLNbrParticles - 1) < 63)
#else
            if ((LLLLzMax + LLLNbrParticles - 1) < 31)
#endif
              {
                LLLSpace = new BosonOnSphereShort (LLLNbrParticles, LLLTotalLz, LLLLzMax);
              }
        }
      else
        {
          int* ReferenceState = 0;
          if (Manager.GetString("reference-file") == 0)
            {
              cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
              return -1;
            }
          ConfigurationParser ReferenceStateDefinition;
          if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
            {
              ReferenceStateDefinition.DumpErrors(cout) << endl;
              return -1;
            }
          if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == false) || (LLLNbrParticles <= 0))
            {
              cout << "NbrParticles is not defined or as a wrong value" << endl;
              return -1;
            }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LLLLzMax) == false) || (LLLLzMax < 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return -1;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
	      return -1;
	    }
	  
	  if (MaxNbrLz != (LLLLzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;
	    }
#ifdef  __64_BITS__
	  if (LLLLzMax  < 63)
#else
	    if (LLLLzMax  < 31)
#endif
	      LLLSpace = new BosonOnSphereHaldaneBasisShort (LLLNbrParticles,LLLTotalLz, LLLLzMax, ReferenceState);
        }
    }
  else
    {
      if (HaldaneBasisFlag == false)
        {
#ifdef __64_BITS__
          if (LLLLzMax <= 63)
#else
            if (LLLLzMax <= 31)
#endif
              {
                LLLSpace = new FermionOnSphere (LLLNbrParticles,LLLTotalLz,LLLLzMax);
              }
            else
              {
                cout << "This lowest Landau level space requires FermionOnSphereLong class for which this kind of calculation is not available"<<endl;
                return -1;
              }
        }
      else
        {
          int * ReferenceState = 0;
          ConfigurationParser ReferenceStateDefinition;
          if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
            {
              ReferenceStateDefinition.DumpErrors(cout) << endl;
              return -1;
            }
          if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == false) || (LLLNbrParticles <= 0))
            {
              cout << "NbrParticles is not defined or as a wrong value" << endl;
              return -1;
            }
          if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LLLLzMax) == false) || (LLLLzMax < 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return 0;
	    }
          int MaxNbrLz;
          if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
            {
              cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
              return -1;
            }
          if (MaxNbrLz != (LLLLzMax + 1))
            {
              cout << "wrong LzMax value in ReferenceState" << endl;
              return -1;
            }
#ifdef  __64_BITS__
          if (LLLLzMax  < 63)
#else
            if (LLLLzMax  < 31)
#endif
              LLLSpace = new FermionOnSphereHaldaneBasis (LLLNbrParticles,LLLTotalLz, LLLLzMax, ReferenceState);
        }
    }
  
  ParticleOnSphere * SpaceLL=0;
  int LzMaxUp = LzMaxFermion + 2;
  int LzMaxDown = LzMaxFermion;
  switch (NbrLL)
    {
    case 2:
      {
	SpaceLL = new FermionOnSphereTwoLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxUp, LzMaxDown);
      }
      break;
    case 3:
      SpaceLL = new FermionOnSphereThreeLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxFermion);
      break;
    case 4:
      SpaceLL = new FermionOnSphereFourLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxFermion);
      break;
    }
  
  
  ParticleOnSphere * FinalSpace = 0;
  if (LLLFermionFlag == false)
    {
      if (Projection == true)
	{
	  if (Manager.GetBoolean("projected-haldane") == false)
	    {
	      if (ReverseFluxFlag == false)
		{
		  FinalSpace = new FermionOnSphere (LLLNbrParticles, TotalLzFermion, LzMaxDown + LLLLzMax);
		}
	      else
		{
		  FinalSpace = new FermionOnSphere (LLLNbrParticles, TotalLzFermion,  LLLLzMax - LzMaxDown);
		}
	    }
	  else
	    {
	      int * ReferenceState = 0;
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(Manager.GetString("projected-referencefile")) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == false) || (LLLNbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int TmpLzMax = 0;
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", TmpLzMax) == false) || (TmpLzMax < 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return 0;
		}
	      if (TmpLzMax != (LzMaxDown + LLLLzMax))
		{
		  cout << "LzMax does not match the one of the projected state" << endl;
		  return 0;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
		  return -1;
		}
	      if (MaxNbrLz != (TmpLzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;
		}
#ifdef  __64_BITS__
	      if (LLLLzMax  < 63)
#else
		if (LLLLzMax  < 31)
#endif
		  FinalSpace = new FermionOnSphereHaldaneBasis (LLLNbrParticles, TotalLzFermion, TmpLzMax, ReferenceState);
	    }
	}
      else
	FinalSpace = new FermionOnSphereTwoLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxUp + LLLLzMax, LzMaxDown + LLLLzMax);
    }
  else
    {
      if (Projection == true)
	{
	  if (Manager.GetBoolean("projected-haldane") == false)
	    {
	      if (ReverseFluxFlag == false)
		{
		  FinalSpace = new BosonOnSphereShort (LLLNbrParticles, TotalLzFermion, LzMaxDown + LLLLzMax);
		}
	      else
		{
		  FinalSpace = new BosonOnSphereShort (LLLNbrParticles, TotalLzFermion, LLLLzMax - LzMaxDown);
		}
	    }
	  else
	    {
	      int * ReferenceState = 0;
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(Manager.GetString("projected-referencefile")) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return -1;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", LLLNbrParticles) == false) || (LLLNbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return -1;
		}
	      int TmpLzMax = 0;
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", TmpLzMax) == false) || (TmpLzMax < 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return 0;
		}
	      if (TmpLzMax != (LzMaxDown + LLLLzMax))
		{
		  cout << "LzMax does not match the one of the projected state" << endl;
		  return 0;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
		  return -1;
		}
	      if (MaxNbrLz != (TmpLzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return -1;
		}
#ifdef  __64_BITS__
	      if (LLLLzMax  < 63)
#else
		if (LLLLzMax  < 31)
#endif
		  FinalSpace = new BosonOnSphereHaldaneBasisShort (LLLNbrParticles, TotalLzFermion, TmpLzMax, ReferenceState);
	    }
	}
      else
	{
#ifdef  __64_BITS__
	  if ( 2 * (LzMaxDown + LLLLzMax) + 3 + LLLNbrParticles < 63)
#else
	    if ( 2 * (LzMaxDown + LLLLzMax) + 3 + LLLNbrParticles < 31)
#endif
	      {
		FinalSpace = new BosonOnSphereTwoLandauLevels (LLLNbrParticles, TotalLzFermion, LzMaxUp + LLLLzMax, LzMaxDown + LLLLzMax);
	      }
	}
    }
  
  
  char * OutputFileName = new char [100];
  
  if (Manager.GetBoolean("rational") == false)
    {
      if(Projection == false)
	{
	  if (LLLFermionFlag == true)
	    {
	      sprintf (OutputFileName, "bosons_%s_2ll_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
	    }
	  else
	    {
	      sprintf (OutputFileName, "fermions_%s_2ll_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
	    }
	}
      else
	{
	  if(ReverseFluxFlag == false)
	    {
	      if (LLLFermionFlag == true)
		{
		  sprintf (OutputFileName, "bosons_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
		}
	      else
		{
		  sprintf (OutputFileName, "fermions_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
		}
	    }
	  else
	    {
	      if (LLLFermionFlag == true)
		{
		  sprintf (OutputFileName, "bosons_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles, LLLLzMax - LzMaxDown, TotalLzFermion);
		}
	      else
		{
		  sprintf (OutputFileName, "fermions_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LLLLzMax - LzMaxDown, TotalLzFermion);
		}
	    }
	}
    }
  else
    {
      if(Projection == false)
	{
	  if (LLLFermionFlag == true)
	    {
	      sprintf (OutputFileName, "bosons_rational_%s_2ll_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
	    }
	  else
	    {
	      sprintf (OutputFileName, "fermions_rational_%s_2ll_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
	    }
	}
      else
	{
	  if(ReverseFluxFlag == false)
	    {
	      if (LLLFermionFlag == true)
		{
		  sprintf (OutputFileName, "bosons_rational_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
		}
	      else
		{
		  sprintf (OutputFileName, "fermions_rational_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LzMaxDown + LLLLzMax, TotalLzFermion);
		}
	    }
	  else
	    {
	      if (LLLFermionFlag == true)
		{
		  sprintf (OutputFileName, "bosons_rational_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles, LLLLzMax - LzMaxDown, TotalLzFermion);
		}
	      else
		{
		  sprintf (OutputFileName, "fermions_rational_%s_n_%d_2s_%d_lz_%d.0.vec",Manager.GetString("interaction-name"), LLLNbrParticles,LLLLzMax - LzMaxDown, TotalLzFermion);
		}
	    }
	}
    }
  
  int * MatchingConditionsIndex = 0;
  int NbrMatchingConditionsIndex = 0;
  
  if (Constraint == true)
    {
      MatchingConditionsIndex = new int [SpaceLL->GetHilbertSpaceDimension()];
      int * LLConfiguration = new int [NbrLL];
      cout << SpaceLL->GetHilbertSpaceDimension()<<endl;
      for(int k = 0 ; (k < NbrLL) ; k++)
	LLConfiguration[k] = 0;
      for(int i = 0 ; i < SpaceLL->GetHilbertSpaceDimension() ; i++)
	{
	  SpaceLL->LandauLevelOccupationNumber(i, LLConfiguration);
	  bool MatchingConstraint = true;
	  
	  if(EnergyValue == -1)
	    {
	      for(int k = 0 ; k < NbrLL ; k++)
		{
		  if ((ParticleNumberConstrains[k]!=-1)&&(LLConfiguration[k] != ParticleNumberConstrains[k]))
		    MatchingConstraint = false;
		  LLConfiguration[k] = 0;
		}
	    }
	  else
	    {
	      int Energy = 0;
	      for(int p = 0 ; p < NbrLL ; p++)
		{
		  Energy += p * LLConfiguration[p];
		  LLConfiguration[p] = 0;
		}
	      if(Energy != EnergyValue)
		{
		  MatchingConstraint = false;
		}
	    }
	  
	  if(MatchingConstraint == true)
	    {
	      MatchingConditionsIndex[NbrMatchingConditionsIndex] = i;
	      NbrMatchingConditionsIndex++;
	    }
	}
      if(NbrMatchingConditionsIndex == 0)
	{
	  cout <<"No matching conditions state";
	  return -1;
	}
      cout <<"Nbr States to be computed = "<<NbrMatchingConditionsIndex<<endl;
    }
  
  if (Manager.GetBoolean("rational") == false)
    {
      
      RealVector LLLState;
      if (LLLState.ReadVector (LLLFileName) == false)
	{
	  cout << "can't open vector file " << LLLFileName << endl;
	  return -1;
	}
      
      if (LLLSpace->GetHilbertSpaceDimension() != LLLState.GetVectorDimension())
	{
	  cout << "Error, the dimension of the LLL vector (" << LLLState.GetVectorDimension() 
	       << ") is not equal to the Hilbert space dimension (" << LLLSpace->GetHilbertSpaceDimension() << ")" << endl;;
	  return -1;
	}
      if (Step > LLLSpace->GetHilbertSpaceDimension())
	{
	  cout << " The LLL space cannot be divided in less than 1 dimension space"<<endl;
	  return -1;
	}
      
      RealVector * FermionState = 0;
      RealVector * OutputVector = 0;
      
      if (Constraint == false)
	{
	  if (Index == -1)
	    {
	      FermionState = new RealVector;
	      if (IsFile(FermionFileName) == false)
		{
		  cout << "state " << FermionFileName << " does not exist or can't be opened" << endl;
		  return -1;
		}
	      if (FermionState->ReadVector(FermionFileName) == false)
		{
		  cout << "can't open vector file " << FermionFileName << endl;
		  return -1;
		}
	      if (SpaceLL->GetHilbertSpaceDimension() != FermionState->GetVectorDimension())
		{
		  cout << "Error, the dimension of the fermionic vector " << FermionFileName << " (" 
		       << FermionState->GetVectorDimension() << ") is not equal to the Hilbert space dimension (" 
		       << SpaceLL->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	    }
	  else
	    {
	      FermionState = new RealVector(SpaceLL->GetHilbertSpaceDimension(),true);
	      (*FermionState)[Index]=1;
	    }
	  if (Resume == true)
	    {
	      char* ResumeVectorName = new char [256];
	      sprintf (ResumeVectorName, "temporary_projection_vector.vec");
	      char * LogFile = new char [256];
	      sprintf (LogFile, "projection.dat");
	      if (IsFile(LogFile) == false)
		{
		  cout << "The calcul cannot be resume as the LogFile doesn't exist." << endl;
		  return -1;
		}
	      ifstream File;
	      File.open(LogFile , ios::in);
	      
	      if (!File.is_open())
		{
		  cout << "Cannot open the file: " << LogFile<< endl;
		  return -1;
		}
	      File.seekg (0, ios::beg);
	      File >> MinComponent;
	      File.close();
	      if (IsFile(ResumeVectorName) == false)
		{
		  cout << "The calcul cannot be resume as there is no vector saved." << endl;
		  return -1;
		}
	      OutputVector = new RealVector;
	      if (OutputVector->ReadVector (ResumeVectorName) == false)
		{
		  cout << "can't open vector file " << ResumeVectorName << endl;
		  return -1;
		}
	      if (FinalSpace->GetHilbertSpaceDimension() != OutputVector->GetVectorDimension())
		{
		  cout << "Error, the dimension of the output vector " << ResumeVectorName << " (" 
		       << OutputVector->GetVectorDimension() << ") is not equal to the Hilbert space dimension (" 
		       << FinalSpace->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	    }
	  else
	    {
	      OutputVector = new RealVector(FinalSpace->GetHilbertSpaceDimension(),true);
	    }
	  
	  FQHESphereMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, LLLSpace, FinalSpace, FermionState, &LLLState, OutputVector, MinComponent, NbrComponents, Projection, Step, NbrLL, Symmetric, ReverseFluxFlag);
	  
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  if(NbrComponents+MinComponent != FinalSpace->GetHilbertSpaceDimension())
	    {
	      char* OutputFileName = new char [256];
	      sprintf (OutputFileName, "temporary_projection_vector.vec");
	      char * LogFile = new char [256];
	      sprintf (LogFile, "projection.dat");
	      (*OutputVector).WriteVector(OutputFileName);
	      ofstream File;
	      File.open(LogFile, ios::binary | ios::out);
	      File.precision(14);
	      File << NbrComponents+MinComponent;
	      File.close();
	    }
	  
	  if((Projection == true)||(NbrLL == 2))
	    {
	      if(Manager.GetBoolean("normalize"))
		{
		  FinalSpace->ConvertFromUnnormalizedMonomial((*OutputVector),0,true);
		}
	      else
		{
		  PutFirstComponentToOne((*OutputVector));
		}
	    }
	  
	  
	  
	  (*OutputVector).WriteVector(OutputFileName);
	  return 0;
	}
      else
	{
	  
	  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, LLLSpace, FinalSpace,MatchingConditionsIndex , &LLLState, NbrLL,NbrMatchingConditionsIndex ,Projection, Symmetric,Manager.GetBoolean("normalize"),OutputFileName,ResumingIndex,ReverseFluxFlag);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  return 0;
	}
    }
  else
    {
      
      LongRationalVector LLLState;
      if (LLLState.ReadVector (LLLFileName) == false)
	{
	  cout << "can't open vector file " << LLLFileName << endl;
	  return -1;
	}
      
      if (LLLSpace->GetHilbertSpaceDimension() != LLLState.GetVectorDimension())
	{
	  cout << "Error, the dimension of the LLL vector (" << LLLState.GetVectorDimension() 
	       << ") is not equal to the Hilbert space dimension (" << LLLSpace->GetHilbertSpaceDimension() << ")" << endl;;
	  return -1;
	}
      if (Step > LLLSpace->GetHilbertSpaceDimension())
	{
	  cout << " The LLL space cannot be divided in less than 1 dimension space"<<endl;
	  return -1;
	}
      
      LongRationalVector* FermionState = 0;
      LongRationalVector* OutputVector = 0;
      
      if (Constraint == false)
	{
	  if(Index == -1)
	    {
	      FermionState = new LongRationalVector;
	      if (IsFile(FermionFileName) == false)
		{
		  cout << "state " << FermionFileName << " does not exist or can't be opened" << endl;
		  return -1;
		}
	      if (FermionState->ReadVector(FermionFileName) == false)
		{
		  cout << "can't open vector file " << FermionFileName << endl;
		  return -1;
		}
	      if(SpaceLL->GetHilbertSpaceDimension() != FermionState->GetVectorDimension())
		{
		  cout << "Error, the dimension of the fermionic vector " << FermionFileName << " (" 
		       << FermionState->GetVectorDimension() << ") is not equal to the Hilbert space dimension (" 
		       << SpaceLL->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	    }
	  else
	    {
	      FermionState = new LongRationalVector(SpaceLL->GetHilbertSpaceDimension(),true);
	      (*FermionState)[Index] = 1l;
	    }
	  if(Resume == true)
	    {
	      char* ResumeVectorName = new char [256];
	      sprintf (ResumeVectorName, "temporary_projection_vector.vec");
	      char * LogFile = new char [256];
	      sprintf (LogFile, "projection.dat");
	      if (IsFile(LogFile) == false)
		{
		  cout << "The calcul cannot be resume as the LogFile doesn't exist." << endl;
		  return -1;
		}
	      ifstream File;
	      File.open(LogFile , ios::in);
	      
	      if (!File.is_open())
		{
		  cout << "Cannot open the file: " << LogFile<< endl;
		  return -1;
		}
	      File.seekg (0, ios::beg);
	      File >> MinComponent;
	      File.close();
	      if (IsFile(ResumeVectorName) == false)
		{
		  cout << "The calcul cannot be resume as there is no vector saved." << endl;
		  return -1;
		}
	      OutputVector = new LongRationalVector;
	      if (OutputVector->ReadVector (ResumeVectorName) == false)
		{
		  cout << "can't open vector file " << ResumeVectorName << endl;
		  return -1;
		}
	      if(FinalSpace->GetHilbertSpaceDimension() != OutputVector->GetVectorDimension())
		{
		  cout <<"Number of rows of the resume vector is not equal to the Hilbert space dimension!" <<endl;
		  return -1;
		}
	    }
	  else
	    {
	      OutputVector = new LongRationalVector (FinalSpace->GetHilbertSpaceDimension(),true);
	    }
	  FQHESphereMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, LLLSpace, FinalSpace, FermionState, &LLLState, OutputVector, MinComponent, NbrComponents, Projection, Step, NbrLL, Symmetric, ReverseFluxFlag);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  
	  if(NbrComponents+MinComponent!=FinalSpace->GetHilbertSpaceDimension())
	    {
	      char* OutputFileName = new char [256];
	      sprintf (OutputFileName, "temporary_projection_vector.vec");
	      char * LogFile = new char [256];
	      sprintf (LogFile, "projection.dat");
	      (*OutputVector).WriteVector(OutputFileName);
	      ofstream File;
	      File.open(LogFile, ios::binary | ios::out);
	      File.precision(14);
	      File << NbrComponents+MinComponent;
	      File.close();
	    }
	  
	  if(Projection == true)
	    {
	      PutFirstComponentToOne((*OutputVector));
	    }
	  
	  (*OutputVector).WriteVector(OutputFileName);
	  ofstream File;
	  return 0;
	}
      else
	{
	  FQHESphereMultipleMonomialsTimesSlaterProjectionOperation Operation(SpaceLL, LLLSpace, FinalSpace,MatchingConditionsIndex , &LLLState, NbrLL,NbrMatchingConditionsIndex ,Projection, Symmetric,Manager.GetBoolean("normalize"),OutputFileName,ResumingIndex,ReverseFluxFlag);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  return 0;
	}
    }
}

// convert a state such that its first nonzero component is equal to 1
//
// state = reference to the state to convert
// return value = converted state

RealVector PutFirstComponentToOne(RealVector& state)
{
  long Index = 0l; 
  while(fabs(state[Index])<1e-10)
    Index++;
  double Normalization = state[Index]; 
  for (long i = 0l ; i< state.GetVectorDimension() ; i++)
    {
      state[i] /= Normalization;
    }
  return state;
}

// convert a state such that its first nonzero component is equal to 1
//
// state = reference to the state to convert
// return value = converted state

LongRationalVector PutFirstComponentToOne(LongRationalVector& state)
{
  long Index = 0l; 
  while(state[Index].IsZero())
    Index++;
  LongRational Normalization = state[Index]; 
  for (long i = 0l ; i< state.GetVectorDimension() ; i++)
    {
      state[i] /= Normalization;
    }
  return state;
}
