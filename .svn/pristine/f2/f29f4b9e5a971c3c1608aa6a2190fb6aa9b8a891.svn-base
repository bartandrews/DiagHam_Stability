#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/FermionOnTorus.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
// #include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <limits>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMultipleU1ToU1" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
	
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('1', "state-1", "vector file that corresponds to the first component");
  (*SystemGroup) += new SingleStringOption  ('2', "state-2", "vector file that corresponds to the second component");
  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "text file that give the list of vectors to symmetrize");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('s', "single-state", "vector file that corresponds to the second component");
  (*SystemGroup) += new BooleanOption  ('a', "sym-y", "apply antiperiodic conditions with respect to Ly before symmetrizing");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals to group together when using the single-state option", 2);
  (*SystemGroup) += new BooleanOption  ('t', "twisted-torus", "apply single particle phases during symmetrization to compensate for the pi/4 tilting of the torus");
  (*SystemGroup) += new BooleanOption  ('\n', "subset-symmetrization", "symmetrize by picking equally space orbitals");
//   (*SystemGroup) += new SingleIntegerOption ('\n', "subset-periodicity", "distance between two consecutive orbitals to keep when using --subset-symmetrization", 2);
  (*SystemGroup) += new SingleIntegerOption ('\n', "subset-shift", "index of the first orbital to keep when using --subset-symmetrization", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "subset-twist", "add an additional phase (in pi units) for each kept and occupied orbital when using --subset-symmetrization", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ky-momentum", "compute the vector with given ky in mode sym-y", 0);
  
  (*SystemGroup) += new SingleDoubleOption ('\n', "precision", "if the norm of the symmetrized vector is below this threshold, it is considered to be zero", MACHINE_PRECISION);
//   (*SystemGroup) += new SingleBooleanOption  ('x', "magnetic-translation", "vector file that corresponds to the second component");
  
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file names");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("state-1") == 0) && (Manager.GetString("multiple-states") == 0))
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if ((Manager.GetString("state-1") != 0) && (IsFile(Manager.GetString("state-1")) == false))
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  

  bool HaveComplex = Manager.GetBoolean("complex");
  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalKy1 = 0;
  bool Statistics = true;
  ParticleOnTorus* Space1 = 0;
  bool UnnormalizedBasisFlag = false;
  bool SingleStateFlag = Manager.GetBoolean("single-state");
  int NbrStates = 1;
  double Precision = Manager.GetDouble("precision");

  if (Manager.GetString("state-1") != 0)
    {  
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						      NbrParticles1, NbrFluxQuanta1, TotalKy1, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  Space1 = new FermionOnTorus(NbrParticles1, NbrFluxQuanta1, TotalKy1);
	}
      else
	{ 
	  Space1 = new BosonOnTorusShort(NbrParticles1, NbrFluxQuanta1, TotalKy1);
	}
    }
  else
    {
      MultiColumnASCIIFile MultipleStateFile;
      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
	{
	  MultipleStateFile.DumpErrors(cout);
	  return -1;
	}
      NbrStates = MultipleStateFile.GetNbrLines();
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(MultipleStateFile(0, 0),
						      NbrParticles1, NbrFluxQuanta1, TotalKy1, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << MultipleStateFile(0, 0) << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  Space1 = new FermionOnTorus(NbrParticles1, NbrFluxQuanta1, TotalKy1);      
	}
      else
	{  
	  Space1 = new BosonOnTorusShort(NbrParticles1, NbrFluxQuanta1, TotalKy1);      
	}
    }

  if ((SingleStateFlag == true) && ((NbrFluxQuanta1 % Manager.GetInteger("nbr-orbitals")) != 0))
    {
      cout << "error: number of flux quanta should be a multiple of " << Manager.GetInteger("nbr-orbitals") << " to perform symmetrization operation on single state"<< endl;
      return -1;
    }
	  
  
  ParticleOnTorus** InputSpaces = 0;
  ParticleOnTorus* TargetSpace = 0; 
  int TotalNbrParticles = NbrParticles1;
  int TotalKy = TotalKy1;
  if (SingleStateFlag == false)    
    {      
      if ((Manager.GetString("state-2") == 0) && (Manager.GetString("multiple-states") == 0))
	{
	  cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHETorusMultipleU1ToU1 -h" << endl;
	  return -1;
	}
      if (Manager.GetString("state-2") != 0)
	{
	  NbrStates = 2;
	  if (IsFile(Manager.GetString("state-2")) == false)
	    {
	      cout << "can't open file " << Manager.GetString("state-2") << endl;
	    }
	  InputSpaces = new ParticleOnTorus*[NbrStates];
	  InputSpaces[0] = Space1;
	  int NbrParticles2 = 0; 
	  int NbrFluxQuanta2 = 0; 
	  int TotalKy2 = 0;
	  Statistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
							  NbrParticles2, NbrFluxQuanta2, TotalKy2, Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
	      return -1;
	    }
	  if (NbrFluxQuanta2 != NbrFluxQuanta1)
	    {
	      cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
	      return -1;
	    }
	  TotalNbrParticles += NbrParticles2;
	  TotalKy += TotalKy2;
	  if (Statistics == true)
	    {
	      InputSpaces[1] = new FermionOnTorus(NbrParticles2, NbrFluxQuanta1, TotalKy2);
	    }
	  else
	    {
	      InputSpaces[1] = new BosonOnTorusShort(NbrParticles2, NbrFluxQuanta1, TotalKy2);
	    }
	} 
      else
	{
	  MultiColumnASCIIFile MultipleStateFile;
	  if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
	    {
	      MultipleStateFile.DumpErrors(cout);
	      return -1;
	    }
	  InputSpaces = new ParticleOnTorus*[NbrStates];
	  InputSpaces[0] = Space1;
	  for (int i = 1; i < NbrStates; ++i)
	    {
	      int NbrParticles2 = 0; 
	      int NbrFluxQuanta2 = 0; 
	      int TotalKy2 = 0;
	      Statistics = true;
	      if (FQHEOnTorusFindSystemInfoFromVectorFileName(MultipleStateFile(0, i),
							      NbrParticles2, NbrFluxQuanta2, TotalKy2, Statistics) == false)
		{
		  cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
		  return -1;
		}
	      if (NbrFluxQuanta2 != NbrFluxQuanta1)
		{
		  cout << "error, " << Manager.GetString("state-1") << " and " << MultipleStateFile(0, i) << " don't have the same number of flux quanta" << endl;
		  return -1;
		}
	      TotalNbrParticles += NbrParticles2;
	      TotalKy += TotalKy2;
	      if (Statistics == true)
		{
		  InputSpaces[i] = new FermionOnTorus(NbrParticles2, NbrFluxQuanta1, TotalKy2);
		}
	      else
		{
		  InputSpaces[i] = new BosonOnTorusShort(NbrParticles2, NbrFluxQuanta1, TotalKy2);
		}
	    }
	}
      TotalKy %= NbrFluxQuanta1;
      cout << "target space N="<< TotalNbrParticles << " Nphi=" << NbrFluxQuanta1 << " Ky=" << TotalKy << endl;
      if (Statistics == true)
	{
	  TargetSpace = new FermionOnTorus(TotalNbrParticles, NbrFluxQuanta1, TotalKy);	       
	}
      else
	{
	  TargetSpace = new BosonOnTorusShort(TotalNbrParticles, NbrFluxQuanta1, TotalKy);	       
	}
    }

  if (HaveComplex == false)
    {
      if (SingleStateFlag == false)
	{	      
	  RealVector* States = new RealVector[NbrStates];
	  if ((Manager.GetString("state-1") != 0) && (Manager.GetString("state-2") != 0))
	    {
	      if (States[0].ReadVector (Manager.GetString("state-1")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
		  return -1;      
		}
	      if (States[1].ReadVector (Manager.GetString("state-2")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-2") << endl;
		  return -1;      
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile MultipleStateFile;
	      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
		{
		  MultipleStateFile.DumpErrors(cout);
		  return -1;
		}
	      for (int i = 0; i < NbrStates; ++i)
		{
		  if (States[i].ReadVector (MultipleStateFile(0, i)) == false)
		    {
		      cout << "can't open vector file " << MultipleStateFile(0, i) << endl;
		      return -1;      
		    }
		}
	    }
	  
	  char* OutputFileName = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
	    }
	  else
	    {
	      OutputFileName = new char [512];
	      if (Statistics == true)
		{
		  sprintf (OutputFileName, "fermions_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
		}
	      else
		{
		  sprintf (OutputFileName, "bosons_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
		}
	    }
	  
      	  RealVector OutputState;
	  if (NbrStates == 2)
	    OutputState = TargetSpace->SymmetrizeU1U1State (States[0], States[1], InputSpaces[0], InputSpaces[1], false, Architecture.GetArchitecture());	    
	  else
	    OutputState = TargetSpace->SymmetrizeU1U1State (States, InputSpaces, NbrStates, Precision, Architecture.GetArchitecture());
	  if (OutputState.Norm() > Precision)
	  {
	    OutputState /= OutputState.Norm();
	    if (OutputState.WriteVector(OutputFileName) == false)
	      {
		cout << "error while writing output state " << OutputFileName << endl;
		return -1;
	      }
	  }
	  
	  delete TargetSpace;
	}
      else
	{
	  RealVector State1;
	  if (State1.ReadVector (Manager.GetString("state-1")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	      return -1;      
	    }      
	  int NbrFluxQuanta = NbrFluxQuanta1 / Manager.GetInteger("nbr-orbitals");
	  char* OutputFileName = 0;
	  char* FullOutputFileName = 0;
	  RealVector* OutputStates = 0;
	  int* KySectors = 0;
	  int* NbrParticleSectors = 0;
	  int NbrKySectors = 0;
	  if (Manager.GetBoolean("sym-y") == false)
	    {
	      if (Manager.GetBoolean("subset-symmetrization") == true)
		{
		  if (Manager.GetDouble("subset-twist") == 0.0)
		    {
		      NbrKySectors = Space1->SymmetrizeSingleStatePeriodicSubsetOrbitals(State1, Manager.GetInteger("subset-shift"),
											 Manager.GetInteger("nbr-orbitals"),
											 OutputStates, NbrParticleSectors, KySectors, Architecture.GetArchitecture()); 
		    }
		  else
		    {
		      cout << "subset-twist requires complex vectors" << endl;
		      return -1;
		    }
		}
	      else
		{
		  NbrKySectors = Space1->SymmetrizeSingleStateGroupingDistantOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture(), Precision); 
		}
	    }
	  else
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingNeighbouringOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture(),  Precision); 	      
	    }
	  if ((NbrKySectors > 0) && (NbrParticleSectors == 0))
	    {
	      NbrParticleSectors = new int[NbrKySectors];
	      for (int i = 0; i < NbrKySectors; ++i)
		{
		  NbrParticleSectors[i] = NbrParticles1;
		}
	    }
	  int NbrGeneratedStates = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
	    }
	  else
	    {
	      OutputFileName = new char [512];
	      if (Manager.GetBoolean("sym-y") == false)
		{
		  if (Manager.GetBoolean("subset-symmetrization") == true)
		    {
		      if (Statistics == true)
			sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld", 
			       TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"));
		      else
			sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld", 
			       TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"));
		    }
		  else
		    {
		      if (Statistics == true)
			sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_xsymmetrized", TotalKy1);
		      else
			sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_xsymmetrized", TotalKy1);
		    }
		}
	      else
		{
		  if (Statistics == true)
		    sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_ysymmetrized", TotalKy1);
		  else
		    sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized", TotalKy1);
		}
	    }
	  for (int i = 0; i < NbrKySectors; ++i)
	    {
	      cout << "state generated in the N=" << NbrParticleSectors[i] << " Ky=" << KySectors[i] << " sector" << endl;
	      char* FullOutputFileName = new char [256 + strlen(OutputFileName)];
	      sprintf (FullOutputFileName , "%s_n_%d_2s_%d_ky_%d.0.vec", OutputFileName, NbrParticleSectors[i], NbrFluxQuanta, KySectors[i]);	      
	      if (OutputStates[i].WriteVector(FullOutputFileName) == false)
		{
		  cout << "error while writing output state " << FullOutputFileName << endl;
		  return -1;
		}
	      delete[] FullOutputFileName;
	    }
	  cout << "Symmetrization has generated " << NbrKySectors << " state(s)" << endl;
	}
    }
  else
    {
      if (SingleStateFlag == false)
	{	  
	  ComplexVector* States = new ComplexVector[NbrStates];
	  if ((Manager.GetString("state-1") != 0) && (Manager.GetString("state-2") != 0))
	    {
	      if (States[0].ReadVector (Manager.GetString("state-1")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
		  return -1;      
		}
	      if (States[1].ReadVector (Manager.GetString("state-2")) == false)
		{
		  cout << "can't open vector file " << Manager.GetString("state-2") << endl;
		  return -1;      
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile MultipleStateFile;
	      if (MultipleStateFile.Parse(Manager.GetString("multiple-states")) == false)
		{
		  MultipleStateFile.DumpErrors(cout);
		  return -1;
		}
	      for (int i = 0; i < NbrStates; ++i)
		{
		  if (States[i].ReadVector (MultipleStateFile(0, i)) == false)
		    {
		      cout << "can't open vector file " << MultipleStateFile(0, i) << endl;
		      return -1;      
		    }
		}
	    }
	  
	  char* OutputFileName = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
	    }
	  else
	    {
	      OutputFileName = new char [512];
	      if (Statistics == true)
		sprintf (OutputFileName, "fermions_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
	      else
		sprintf (OutputFileName, "bosons_torus_kysym_symmetrized_n_%d_2s_%d_ky_%d.0.vec", TotalNbrParticles, NbrFluxQuanta1, TotalKy);
	    }
	  
	  	  
	  ComplexVector OutputState = TargetSpace->SymmetrizeU1U1State (States, InputSpaces, NbrStates, Precision,  Architecture.GetArchitecture());

	  if (OutputState.Norm() > Precision)
	  {
	    OutputState /= OutputState.Norm();
	    if (OutputState.WriteVector(OutputFileName) == false)
	      {
		cout << "error while writing output state " << OutputFileName << endl;
		return -1;
	      }
	  }
	  
	  delete TargetSpace;
	}
      else
	{
	  ComplexVector State1;
	  if (State1.ReadVector (Manager.GetString("state-1")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	      return -1;      
	    }      
	  int NbrFluxQuanta = NbrFluxQuanta1 / Manager.GetInteger("nbr-orbitals");
	  char* OutputFileName = 0;
	  char* FullOutputFileName = 0;
	  ComplexVector* OutputStates = 0;
	  int* KySectors = 0;
	  int* NbrParticleSectors = 0;
	  int NbrKySectors = 0;
	  if (Manager.GetBoolean("sym-y") == false)
	    {
	      if (Manager.GetBoolean("subset-symmetrization") == true)
		{
		  if (Manager.GetDouble("subset-twist") == 0.0)
		    {
		      NbrKySectors = Space1->SymmetrizeSingleStatePeriodicSubsetOrbitals(State1, Manager.GetInteger("subset-shift"),
											 Manager.GetInteger("nbr-orbitals"),
											 OutputStates, NbrParticleSectors, KySectors, Architecture.GetArchitecture()); 
		    }
		  else
		    {
		      NbrKySectors = Space1->SymmetrizeSingleStatePeriodicSubsetOrbitals(State1, Manager.GetInteger("subset-shift"),
											 Manager.GetInteger("nbr-orbitals"), Manager.GetDouble("subset-twist"),
											 OutputStates, NbrParticleSectors, KySectors, Architecture.GetArchitecture()); 
		    }
		}
	      else
		{
		  NbrKySectors = Space1->SymmetrizeSingleStateGroupingDistantOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture(), Precision, Manager.GetBoolean("twisted-torus")); 
		}
	    }
	  else
	    {
	      NbrKySectors = Space1->SymmetrizeSingleStateGroupingNeighbouringOrbitals(State1, Manager.GetInteger("nbr-orbitals"), OutputStates, KySectors, Architecture.GetArchitecture(), Precision); 	      
	    }
	  if ((NbrKySectors > 0) && (NbrParticleSectors == 0))
	    {
	      NbrParticleSectors = new int[NbrKySectors];
	      for (int i = 0; i < NbrKySectors; ++i)
		{
		  NbrParticleSectors[i] = NbrParticles1;
		}
	    }
	  int NbrGeneratedStates = 0;
	  if (Manager.GetString("output-file") != 0)
	    {
	      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
	      strcpy (OutputFileName, Manager.GetString("output-file"));
		}
	  else
	    {
	      OutputFileName = new char [512];
	      if (Manager.GetBoolean("sym-y") == false)
		{
		  if (Manager.GetBoolean("subset-symmetrization") == true)
		    {
		      if (Manager.GetDouble("subset-twist") == 0.0)
			{
			  if (Statistics == true)
			    sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld", 
				     TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"));
			  else
			    sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld", 
				     TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"));
			}
		      else
			{
			  if (Statistics == true)
			    sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld_twist_%.6f", 
				     TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"), Manager.GetDouble("subset-twist"));
			  else
			    sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized_subset_shift_%ld_period_%ld_twist_%.6f", 
				     TotalKy1, Manager.GetInteger("subset-shift"), Manager.GetInteger("nbr-orbitals"), Manager.GetDouble("subset-twist"));
			}
		    }
		  else
		    {
		      if (Statistics == true)
			sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_xsymmetrized", TotalKy1);
		      else
			sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_xsymmetrized", TotalKy1);
		    }
		}
	      else
		{
		  if (Statistics == true)
		    sprintf (OutputFileName, "fermions_torus_kysym_sourceky_%d_ysymmetrized", TotalKy1);
		  else
		    sprintf (OutputFileName, "bosons_torus_kysym_sourceky_%d_ysymmetrized", TotalKy1);
		}
	    }
	  for (int i = 0; i < NbrKySectors; ++i)
	    {
	      cout << "state generated in the N=" << NbrParticleSectors[i] << " Ky=" << KySectors[i] << " sector" << endl;
	      char* FullOutputFileName = new char [256 + strlen(OutputFileName)];
	      sprintf (FullOutputFileName , "%s_n_%d_2s_%d_ky_%d.0.vec", OutputFileName, NbrParticleSectors[i], NbrFluxQuanta, KySectors[i]);	      
	      if (OutputStates[i].WriteVector(FullOutputFileName) == false)
		{
		  cout << "error while writing output state " << FullOutputFileName << endl;
		  return -1;
		}
	      delete[] FullOutputFileName;
	    }
	  cout << "Symmetrization has generated " << NbrKySectors << " state(s)" << endl;
	}
    }
  
  delete Space1;
}

