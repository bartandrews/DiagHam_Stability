#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

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


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  ArchitectureManager Architecture;
  Architecture.AddOptionGroup(&Manager);
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-ll", "number of lambda levels", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "flux-quanta", "if positive or zero, for the number of flux quanta for the lambda levels", -1);
  (*SystemGroup) += new BooleanOption ('\n',"reverse-flux","the fluxes bind to each particle are in the opposite direction than the magnetic field");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "lambda-state", "provide a many-body state for the lambda levels, instead of using filled lambda levels");
  (*SystemGroup) += new BooleanOption ('\n',"disable-szsymmetry","disable the Sz<->-Sz symmetry for the Sz=0 sector");
  (*SystemGroup) += new BooleanOption ('\n',"disable-lzsymmetry","disable the Lz<->-Lz symmetry for the Lz=0 sector");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "cylinder-perimeter", "if non zero, fix the cylinder perimeter (in magnetic length unit) instead of the aspect ratio", 0);

  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the projected state assuming the sphere geometry");
  (*OutputGroup) += new BooleanOption ('\n', "cylinder-normalize", "normalize the projected state assuming the cylinder geometry");
  (*OutputGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*OutputGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*OutputGroup) += new SingleIntegerOption ('\n', "outputvector-index", "set the index of the output vector (i.e. the integer in the extention *.xxx.vec)", 0);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereWithSU2SpinLLLProjectedLLLTimesManyLL -h" << endl;
      return -1;
    }	
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles  = Manager.GetInteger("nbr-particles");
  int TotalSz = Manager.GetInteger("total-sz");;
  int NbrParticleUp = (NbrParticles + TotalSz) / 2;
  int NbrParticleDown = (NbrParticles - TotalSz) / 2;  
  int NbrLandauLevel = Manager.GetInteger("nbr-ll");
  int NbrFluxQuantumLambdaLevels = 0;
  if (Manager.GetInteger("flux-quanta") < 0)
    {
      if (NbrParticleUp >= NbrParticleDown)
	NbrFluxQuantumLambdaLevels = NbrParticleUp - (NbrLandauLevel * (NbrLandauLevel - 1));
      else
	NbrFluxQuantumLambdaLevels = NbrParticleDown - (NbrLandauLevel * (NbrLandauLevel - 1));
      if ((NbrFluxQuantumLambdaLevels % NbrLandauLevel) != 0)
	{
	  cout << "error, the number of particles is not compatible with " <<  NbrLandauLevel << " filled Landau levels" << endl;
	  return -1;
	}
      NbrFluxQuantumLambdaLevels /= NbrLandauLevel;
      --NbrFluxQuantumLambdaLevels;
    }
  else
    {
      if (Manager.GetString("lambda-state") == 0)
	{
	  cout << "error, --flux-quanta requires to provide an input vector using --lambda-state" << endl;
	  return -1;
	}
      NbrFluxQuantumLambdaLevels = Manager.GetInteger("flux-quanta");
    }
  int NbrFluxQuanta = (NbrParticles - 1);
  if (Manager.GetBoolean("reverse-flux") == false)
    NbrFluxQuanta += NbrFluxQuantumLambdaLevels;
  else
    NbrFluxQuanta -= NbrFluxQuantumLambdaLevels;
  int TotalLz = 0;
  double Ratio = Manager.GetDouble("aspect-ratio");
  double Perimeter = Manager.GetDouble("cylinder-perimeter");
  if (Perimeter != 0.0)
    {
      Ratio = 2.0 * M_PI * (NbrFluxQuanta + 1) / (Perimeter * Perimeter);
    }
  else
    {
      Perimeter = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) / Ratio);
    }    
  bool LzSymmmetryFlag = true;
  if ((TotalLz != 0) || (Manager.GetBoolean("disable-lzsymmetry") == true))
    {
      LzSymmmetryFlag = false;
    }
  bool SzSymmmetryFlag = true;
  if ((TotalSz != 0) || (Manager.GetBoolean("disable-szsymmetry") == true))
    {
      SzSymmmetryFlag = false;
    }
  char* DiscreteSymmetryName = new char[128];
  if (SzSymmmetryFlag == false)
    {
      if (LzSymmmetryFlag == false)
	{
	  sprintf (DiscreteSymmetryName, "");
	}
      else
	{
	  if (Manager.GetBoolean("minus-lzparity") == false)
	    {
	      sprintf (DiscreteSymmetryName, "_lzsym_1");
	    }
	  else
	    {
	      sprintf (DiscreteSymmetryName, "_lzsym_-1");
	    }
	}
    }
  else
    {
      if (LzSymmmetryFlag == false)
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_1");
	    }
	  else
	    {
	      sprintf (DiscreteSymmetryName, "_szsym_-1");
	    }
	}
      else
	{
	  if (Manager.GetBoolean("minus-szparity") == false)
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_1");
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("minus-lzparity") == false)
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_1_szsym_-1");
		}
	      else
		{
		  sprintf (DiscreteSymmetryName, "_lzsym_-1_szsym_-1");
		}
	    }
	}
    }


  BosonOnSphereWithSU2Spin* OutputSpace = 0;
  if (LzSymmmetryFlag == true)
    {
      if (SzSymmmetryFlag == true)
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"), 
								  Manager.GetBoolean("minus-lzparity"));
	}
      else
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinLzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-lzparity"));
	}
    }
  else
    {
      if (SzSymmmetryFlag == true)
	{
	  OutputSpace = new BosonOnSphereWithSU2SpinSzSymmetry(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"));
	}
      else
	{
	  OutputSpace = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	}
    }

  char* GeometryName = new char[128];
  if (Manager.GetBoolean("cylinder-normalize") == false)
    {
      sprintf (GeometryName, "sphere_su2", Ratio);
    }
  else
    {
      if (Manager.GetDouble("cylinder-perimeter") > 0.0)	
	{
	  sprintf (GeometryName, "cylinder_perimeter_%.6f_su2", Perimeter);
	}
      else
	{
	  sprintf (GeometryName, "cylinder_ratio_%.6f_su2", Ratio);
	}
    }
  char* OutputName = new char [512  + strlen(DiscreteSymmetryName)+ strlen(Manager.GetString("interaction-name")) + strlen(GeometryName)];
  sprintf (OutputName, "bosons_%s%s_%s_n_%d_2s_%d_sz_%d_lz_%d.%ld.vec", GeometryName, DiscreteSymmetryName, Manager.GetString("interaction-name"), 
	   NbrParticles, NbrFluxQuanta, TotalSz, TotalLz, Manager.GetInteger("outputvector-index"));

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
    {
      cout << "generating state " << OutputName << endl;
    }
  RealVector OutputVector (OutputSpace->GetHilbertSpaceDimension(), true);

  if (NbrLandauLevel == 1)
    {
      FermionOnSphereWithSpin* InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuantumLambdaLevels, TotalSz);  
      RealVector InputVector (InputSpace->GetHilbertSpaceDimension(), true);
      InputVector[0] = 1.0;
      OutputSpace->SlaterTimeSpinfulFermionicState(InputVector, OutputVector, InputSpace, 0, InputSpace->GetHilbertSpaceDimension(),
						   !(Manager.GetBoolean("normalize") | Manager.GetBoolean("cylinder-normalize")), 
						   Manager.GetBoolean("cylinder-normalize"), Perimeter, Architecture.GetArchitecture());
      delete InputSpace;
    }
  else
    {
      FermionOnSphereWithSpinTwoLandauLevels* InputSpace = new FermionOnSphereWithSpinTwoLandauLevels(NbrParticles, TotalLz, NbrFluxQuantumLambdaLevels, TotalSz);  
      RealVector InputVector (InputSpace->GetHilbertSpaceDimension(), true);
      if (Manager.GetInteger("flux-quanta") < 0)
	{
	  InputVector[0] = 1.0;
	}
      else
	{
	  if (InputVector.ReadVector(Manager.GetString("lambda-state")) == false)
	    {
	      cout << "error, can't read " << Manager.GetString("lambda-state") << endl;
	      return -1;
	    }
	}
      OutputSpace->SlaterTimeSpinfulFermionicState(InputVector, OutputVector, InputSpace, 0, InputSpace->GetHilbertSpaceDimension(),
						   !(Manager.GetBoolean("normalize") | Manager.GetBoolean("cylinder-normalize")), 
						   Manager.GetBoolean("cylinder-normalize"), Perimeter, Architecture.GetArchitecture());
      delete InputSpace;
    }

  if (Architecture.GetArchitecture()->CanWriteOnDisk())
    {
      if ((Manager.GetBoolean("normalize") == true) || (Manager.GetBoolean("cylinder-normalize") == true))
	OutputVector.Normalize();
      if (OutputVector.WriteVector(OutputName) == false)
	{
	  cout << "can't write " << OutputName << endl;
	}
    }

  delete[] OutputName;
  
  delete OutputSpace;

  return 0;
}

