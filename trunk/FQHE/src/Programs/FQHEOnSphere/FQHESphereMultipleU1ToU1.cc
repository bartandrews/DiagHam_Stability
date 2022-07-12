#include "config.h"

#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereSymmetrizeU1U1StateOperation.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-1", "use squeezed basis instead of the usual n-body basis for the first component");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-2", "use squeezed basis instead of the usual n-body basis for the second component");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane-output", "use squeezed basis instead of the usual n-body basis for the symmetrized state");
  (*SystemGroup) += new BooleanOption  ('\n', "single-state", "symmetrize a unique state, groupings consecutive orbitals");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals to group together when using the single-state option -- distance between two consecutive orbitals to keep when using --subset-symmetrization", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "periodicity-symmetrization", "group orbitals that have the same momentum modulo a given periodicity, instead of neighboring orbitals");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals to group together when using the single-state option", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "subset-symmetrization", "symmetrize by picking equally space orbitals");
//   (*SystemGroup) += new SingleIntegerOption ('\n', "subset-periodicity", "distance between two consecutive orbitals to keep when using --subset-symmetrization", 2);
  (*SystemGroup) += new SingleIntegerOption ('\n', "subset-shift", "index of the first orbital to keep when using --subset-symmetrization", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "conserve-particles", "Only output the vectors with a number of particles conserved");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file1", "use a file as the definition of the reference state for the first component squeezed basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file2", "use a file as the definition of the reference state for the second component squeezed basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-fileoutput", "use a file as the definition of the reference state for the symmetrized state squeezed basis");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file names");
  (*OutputGroup) += new BooleanOption  ('u', "unnormalized-basis", "indicates that calculations and data are in the unnormalized basis");
  (*SystemGroup) += new BooleanOption  ('\n', "rational", "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new BooleanOption  ('\n', "fix-totallz", "fix the toal Lz value in in single-state mode");
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the target system (in single-state mode)", 0);
  
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("state-1") == 0)
    {
      cout << "error, an input file should be provided for the first component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state-1")) == false)
    {
      cout << "can't open file " << Manager.GetString("state-1") << endl;
    }
  if ((Manager.GetString("state-2") == 0) && (Manager.GetBoolean("single-state") == false))
    {
      cout << "error, an input file should be provided for the second component. See man page for option syntax or type FQHESphereMultipleU1ToU1 -h" << endl;
      return -1;
    }
  if ((Manager.GetBoolean("single-state") == false) && (IsFile(Manager.GetString("state-2")) == false))
    {
      cout << "can't open file " << Manager.GetString("state-2") << endl;
    }

  int NbrParticles1 = 0; 
  int NbrFluxQuanta1 = 0; 
  int TotalLz1 = 0;
  bool Statistics = true;
  
  int NbrParticles2 = 0; 
  int NbrFluxQuanta2 = 0; 
  int TotalLz2 = 0;
  
    
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-1"),
						   NbrParticles1, NbrFluxQuanta1, TotalLz1, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state-1") << endl;
      return -1;
    }
  
  
  RealVector State1;
  LongRationalVector RationalState1;
  RealVector State2;
  LongRationalVector RationalState2;
  if (Manager.GetBoolean("rational") == false)
    {
      if (State1.ReadVector (Manager.GetString("state-1")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	  return -1;      
	}
    }
  else
    {
      if (RationalState1.ReadVector (Manager.GetString("state-1")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("state-1") << endl;
	  return -1;      
	}
    }
  

  BosonOnSphereShort* Space1 = 0;
  BosonOnSphereShort* Space2 = 0;
  FermionOnSphere* FermionicSpace1 = 0;
  FermionOnSphere* FermionicSpace2 = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("haldane-1") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file1"), NbrParticles1, NbrFluxQuanta1, ReferenceState) == false)
	    return -1;
	  Space1 = new BosonOnSphereHaldaneBasisShort(NbrParticles1, TotalLz1, NbrFluxQuanta1, ReferenceState);	       
	}
      else
	{
	  Space1 = new BosonOnSphereShort(NbrParticles1, TotalLz1, NbrFluxQuanta1);	       
	}
    }
  else
    {
      if (Manager.GetBoolean("haldane-1") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file1"), NbrParticles1, NbrFluxQuanta1, ReferenceState) == false)
	    return -1;
	  FermionicSpace1 = new FermionOnSphereHaldaneBasis(NbrParticles1, TotalLz1, NbrFluxQuanta1, ReferenceState);	       
	}
      else
	{
	  FermionicSpace1 = new FermionOnSphere(NbrParticles1, TotalLz1, NbrFluxQuanta1);	       
	}
    }
  
  if (Manager.GetBoolean("single-state") == false)
    {
      Statistics = true;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state-2"),
						       NbrParticles2, NbrFluxQuanta2, TotalLz2, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state-2") << endl;
	  return -1;
	}
      /*  if (NbrParticles1 != NbrParticles2)
	  {
	  cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of particles" << endl;
	  return -1;
	  }*/
      if (NbrFluxQuanta2 != NbrFluxQuanta1)
	{
	  cout << "error, " << Manager.GetString("state-1") << " and " << Manager.GetString("state-2") << " don't have the same number of flux quanta" << endl;
	  return -1;
	}
      
      
      if (Manager.GetBoolean("rational") == false)
	{
	  if (State2.ReadVector (Manager.GetString("state-2")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
	      return -1;      
	    }
	}
      else
	{
	  if (RationalState2.ReadVector (Manager.GetString("state-2")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("state-2") << endl;
	      return -1;      
	    }
	}
      
      if (Manager.GetBoolean("haldane-2") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-file2"), NbrParticles2, NbrFluxQuanta2, ReferenceState) == false)
	    return -1;
	  if (Statistics == false)
	    {
	      Space2 = new BosonOnSphereHaldaneBasisShort(NbrParticles2, TotalLz2, NbrFluxQuanta2, ReferenceState);	       
	    }
	  else
	    {
	      FermionicSpace2 = new FermionOnSphereHaldaneBasis(NbrParticles2, TotalLz2, NbrFluxQuanta2, ReferenceState);	       
	    }
	}
      else
	{
	  if (Statistics == false)
	    {
	      Space2 = new BosonOnSphereShort(NbrParticles2, TotalLz2, NbrFluxQuanta2);	       
	    }
	  else
	    {
	      FermionicSpace2 = new FermionOnSphere(NbrParticles2, TotalLz2, NbrFluxQuanta2);
	    }
	}      
    }
  else
    {
//       if (((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals")) != 0) 
// 	{
// 	  cout << "Error: number of orbitals should be a multiple of " << Manager.GetInteger("nbr-orbitals") << " for symmetrization of a single state" << endl;
// 	  return -1;
// 	}
    }
  
  
  
  char* OutputFileName = 0;
  int NbrFluxQuanta = (((NbrFluxQuanta1 + 1) / ((int) Manager.GetInteger("nbr-orbitals"))) - 1);
  if ((Manager.GetBoolean("subset-symmetrization")) && (((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals")) != 0) && (Manager.GetInteger("subset-shift") < ((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals"))))
    NbrFluxQuanta += 1;
  
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [512];
      if (Manager.GetBoolean("single-state") == false)
	{
	  if (Statistics == false)
	    {
	      if (Manager.GetBoolean("unnormalized-basis") == false)
		{
		  if (Manager.GetBoolean("haldane-output") == true)
		    sprintf (OutputFileName, "bosons_haldane_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		  else
		    sprintf (OutputFileName, "bosons_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		}
	      else
		{
		  if (Manager.GetBoolean("haldane-output") == true)
		    sprintf (OutputFileName, "bosons_haldane_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		  else
		    sprintf (OutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("unnormalized-basis") == false)
		{
		  if (Manager.GetBoolean("haldane-output") == true)
		    sprintf (OutputFileName, "fermions_haldane_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		  else
		    sprintf (OutputFileName, "fermions_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		}
	      else
		{
		  if (Manager.GetBoolean("haldane-output") == true)
		    sprintf (OutputFileName, "fermions_haldane_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		  else
		    sprintf (OutputFileName, "fermions_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles1 + NbrParticles2, NbrFluxQuanta2, (TotalLz1 + TotalLz2));
		}
	    }
	}
      else
	{
	  if (Statistics == false)
	    {
	      if (Manager.GetBoolean("unnormalized-basis") == false)
		{
		  sprintf(OutputFileName, "bosons_symmetrized_n_%d_2s_%d", NbrParticles1, NbrFluxQuanta);
		}
	      else
		{
		  sprintf(OutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d", NbrParticles1, NbrFluxQuanta);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("unnormalized-basis") == false)
		{
		  sprintf(OutputFileName, "fermions_symmetrized_n_%d_2s_%d", NbrParticles1,NbrFluxQuanta);
		}
	      else
		{
		  sprintf(OutputFileName, "fermions_unnormalized_symmetrized_n_%d_2s_%d", NbrParticles1, NbrFluxQuanta);
		}
	    }
	}
    }
  
  
  BosonOnSphereShort* TargetSpace = 0;
  FermionOnSphere* FermionicTargetSpace = 0;
  if (Manager.GetBoolean("single-state") == false)
    {
      if (Manager.GetBoolean("haldane-output") == true)
	{
	  int* ReferenceState = 0;
	  int NbrParticles = 0;
	  int NbrFluxQuanta = 0;
	  if (FQHEGetRootPartition(Manager.GetString("reference-fileoutput"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
	    return -1;
	  if (NbrParticles != (NbrParticles1 + NbrParticles2))
	    {
	      cout << "error : number of particles for the output state ( " << NbrParticles << " ) is different from the one computed from the input states (" << (NbrParticles1 + NbrParticles2) << ")" << endl;
	      return -1;
	    }
	  if (NbrFluxQuanta != NbrFluxQuanta2)
	    {
	      cout << "error : number of flux quanta for the output state ( " << NbrFluxQuanta << " ) is different from the one computed from the input states (" << NbrFluxQuanta2 << ")" << endl;
	      return -1;
	    }
	  int TotalLz = 0;
	  if (Statistics == false)
	    {
	      TargetSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta2, ReferenceState);	       
	    }
	  else
	    {
	      FermionicTargetSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta2, ReferenceState);	       
	    }
	  if (TotalLz != (TotalLz1 + TotalLz2))
	    {
	      cout << "error : total Lz for the output state ( " << TotalLz << " ) is different from the one computed from the input states (" << (TotalLz1 + TotalLz2) << ")" << endl;
	      return -1;
	    }
	}
      else
	{
	  if (Statistics == false)
	    {
	      TargetSpace = new BosonOnSphereShort(NbrParticles1 + NbrParticles2, TotalLz1 + TotalLz2, NbrFluxQuanta2);	       
	    }
	  else
	    {
	      FermionicTargetSpace = new FermionOnSphere(NbrParticles1 + NbrParticles2, TotalLz1 + TotalLz2, NbrFluxQuanta2);	  
	    }
	}
      if (Manager.GetBoolean("rational") == false)
	{
	  RealVector OutputState = TargetSpace->SymmetrizeU1U1State (State1 , State2, Space1 , Space2 , Manager.GetBoolean("unnormalized-basis") , Architecture.GetArchitecture());
	  if (OutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }
	}
      else
	{
	  LongRationalVector RationalOutputState;
	  if (Statistics == false)
	    {
	      RationalOutputState = TargetSpace->SymmetrizeU1U1State (RationalState1, RationalState2, Space1, Space2, Architecture.GetArchitecture());
	    }
	  else
	    {
	      RationalOutputState = FermionicTargetSpace->AntiSymmetrizeU1U1State (RationalState1, RationalState2, FermionicSpace1, FermionicSpace2, Architecture.GetArchitecture());
	    }
	  if (RationalOutputState.WriteVector(OutputFileName) == false)
	    {
	      cout << "error while writing output state " << OutputFileName << endl;
	      return -1;
	    }
	}
    }
  else
    {
      if (Statistics == false)
	{
	  int NbrFluxQuanta = (NbrFluxQuanta1 - 1) / 2;
	  int TotalLz = Manager.GetInteger("total-lz");
	  char* FullOutputFileName = new char [strlen(OutputFileName) + 128];
	  sprintf (FullOutputFileName , "%s_lz_%d.0.vec", OutputFileName, (int) Manager.GetInteger("total-lz"));
	  TargetSpace = new BosonOnSphereShort(NbrParticles1, TotalLz, NbrFluxQuanta);
	  RealVector OutputState;
	  LongRationalVector RationalOutputState;
	  if (Manager.GetBoolean("rational") == false)
	    {
	      if (Manager.GetBoolean("fix-totallz") == false)
		{
		  RealVector* OutputStates;
		  int* LzSectors;
		  int NbrLzSectors;
		  int* NbrParticleSectors = 0;
		  int TargetNbrFluxQuanta = 0;
		  if (Manager.GetBoolean("periodicity-symmetrization") == false)
		    {
		      if (Manager.GetBoolean("subset-symmetrization") == true)
			{
			  NbrLzSectors = Space1->SymmetrizeSingleStatePeriodicSubsetOrbitals(State1, Manager.GetInteger("subset-shift"),
											     Manager.GetInteger("nbr-orbitals"),Manager.GetBoolean("unnormalized-basis"),
											     OutputStates, NbrParticleSectors, LzSectors);
			  TargetNbrFluxQuanta = ((NbrFluxQuanta1 + 1) / Manager.GetInteger("nbr-orbitals")) - 1;
			  if ((((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals")) != 0) && (Manager.GetInteger("subset-shift") < ((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals"))))
			    TargetNbrFluxQuanta += 1;
			}
		      else
			{
			  NbrLzSectors = Space1->SymmetrizeSingleStateOneIntoManyOrbital(State1, (int) Manager.GetInteger("nbr-orbitals"), OutputStates, Manager.GetBoolean("unnormalized-basis"), LzSectors);
			}
		    }
		  else
		    {
		      cout << "periodic symmetrization not implemented for real vectors" << endl;
		      return -1;
		    }
		  int NbrGeneratedStates = 0;
		  if ((NbrLzSectors > 0) && (NbrParticleSectors == 0))
		    {
		      NbrParticleSectors = new int[NbrLzSectors];
		      for (int i = 0; i < NbrLzSectors; ++i)
			{
			  NbrParticleSectors[i] = NbrParticles1;
			}
		    }
		  for (int i = 0; i < NbrLzSectors; ++i)
		    {
		      cout << "state generated in the N=" << NbrParticleSectors[i] << " 2*Lz=" << LzSectors[i] << " sector" << endl;
		      if ((Manager.GetBoolean("conserve-particles") == false) || (NbrParticleSectors[i] == NbrParticles1))
			{
			  RealVector& OutputState = OutputStates[i];
			  bool zeroFlag = true;
			  int RootPosition = 0;
			  if (Manager.GetBoolean("unnormalized-basis") == true)
			    {
			      for (int TmpPos = 0; (TmpPos < OutputState.GetVectorDimension()) && (zeroFlag == true); ++TmpPos)
				{
				  if (OutputState[TmpPos] > 1.0e-10)
				    zeroFlag = false;
				  else
				    ++RootPosition;
				}
			      if (zeroFlag)
				{	
				  cout << "this state is null." << endl;
				}
			      else
				{
				  OutputState /= OutputState[RootPosition]; 
				  sprintf (FullOutputFileName , "%s_lz_%d.0.vec", OutputFileName, LzSectors[i]);	      
				  if (OutputState.WriteVector(FullOutputFileName) == false)
				    {
				      cout << "error while writing output state " << FullOutputFileName << endl;
				      return -1;
				    }
				  ++NbrGeneratedStates;
				}	
			    }
			  else
			    {
			      if (OutputState.Norm() < 1.0e-10)
				{
				  cout << "this state is null." << endl;
				}
			      else
				{
				  OutputState /= OutputState.Norm(); 
				  if (Manager.GetBoolean("subset-symmetrization") == false)
				    {
				      if (Manager.GetString("output-file") == 0)
					{
					  if (Manager.GetBoolean("unnormalized-basis") == false)
					    {
					      sprintf(FullOutputFileName, "bosons_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
						      TargetNbrFluxQuanta, LzSectors[i]);
					    }
					  else
					    {
					      sprintf(FullOutputFileName, "bosons_unnormalized_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
						      TargetNbrFluxQuanta, LzSectors[i]);
					    }
					}
				      else
					{
					  sprintf(FullOutputFileName, "%s_lz_%d.0.vec", Manager.GetString("output-file"), LzSectors[i]);
					}
				    }
				  else
				    {
				      if (Manager.GetString("output-file") == 0)
					{
					  if (Manager.GetBoolean("unnormalized-basis") == false)
					    {
					      sprintf(FullOutputFileName, "bosons_subset_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
						      TargetNbrFluxQuanta, LzSectors[i]);
					    }
					  else
					    {
					      sprintf(FullOutputFileName, "bosons_unnormalized_subset_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
						      TargetNbrFluxQuanta, LzSectors[i]);
					    }
					}
				      else
					{
					  sprintf(FullOutputFileName, "%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("output-file"), NbrParticleSectors[i], 
						  TargetNbrFluxQuanta, LzSectors[i]);
					}
				    }      
				  if (OutputState.WriteVector(FullOutputFileName) == false)
				    {
				      cout << "error while writing output state " << FullOutputFileName << endl;
				      return -1;
				    }
				  ++NbrGeneratedStates;
				}
			    }
			}
		    }
		  cout << "Symmetrization has generated " << NbrGeneratedStates << " state(s)" << endl;
		}
	      else
		{
		  
		  if (((NbrParticles1 * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
		    {
		      cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
		      return -1;
		    } 
		  OutputState = TargetSpace->SymmetrizeU1U1SingleState(State1, Space1, Manager.GetBoolean("unnormalized-basis"));
		  if (OutputState.Norm() < 1.0e-10)
		    {
		      cout << "Symmetrized state is zero. No output." << endl;
		      return -1;
		    }
		  else
		    {		  
		      if (Manager.GetBoolean("unnormalized-basis") == true)
			{
			  int RootPosition = 0;
			  while (fabs(OutputState[RootPosition]) < 1.0e-12)
			    {
			      RootPosition += 1;
			    }
			  double RootCoef = OutputState[RootPosition];
			  OutputState /= RootCoef;
			}
		      else
			OutputState /= OutputState.Norm();
		      if (OutputState.WriteVector(FullOutputFileName) == false)
			{
			  cout << "error while writing output state " << FullOutputFileName << endl;
			  return -1;
			}
		    }
		}
	    }
	  else
	    {
	      LongRationalVector* RationalOutputStates;
	      int* LzSectors;
	      int* NbrParticleSectors = 0;
	      int NbrLzSectors;
	      int TargetNbrFluxQuanta = 0;
	      if (Manager.GetBoolean("periodicity-symmetrization") == false)
		{
		  if (Manager.GetBoolean("subset-symmetrization") == true)
		    {
		      NbrLzSectors = Space1->SymmetrizeSingleStatePeriodicSubsetOrbitals(RationalState1, Manager.GetInteger("subset-shift"),
											 Manager.GetInteger("nbr-orbitals"),
											 RationalOutputStates, NbrParticleSectors, LzSectors);
		      TargetNbrFluxQuanta = ((NbrFluxQuanta1 + 1) / Manager.GetInteger("nbr-orbitals")) - 1;
		      if ((((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals")) != 0) && (Manager.GetInteger("subset-shift") < ((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals"))))
			TargetNbrFluxQuanta += 1;
		    }
		  else
		    {
		      
		      NbrLzSectors = Space1->SymmetrizeSingleStateOneIntoManyOrbital(RationalState1, (int) Manager.GetInteger("nbr-orbitals"), RationalOutputStates, LzSectors);
		      TargetNbrFluxQuanta = (((NbrFluxQuanta1 + 1) / ((int) Manager.GetInteger("nbr-orbitals"))) - 1);
		    }
		}
	      else
		{
		  NbrLzSectors = Space1->SymmetrizeSingleStatePeriodicOrbitals(RationalState1, (NbrFluxQuanta1 + 1) / (int) Manager.GetInteger("nbr-orbitals"), RationalOutputStates,
									       LzSectors);
		  TargetNbrFluxQuanta = (((NbrFluxQuanta1 + 1) / ((int) Manager.GetInteger("nbr-orbitals"))) - 1);
		}
	      int NbrGeneratedStates = 0;
	      if ((NbrLzSectors > 0) && (NbrParticleSectors == 0))
		{
		  NbrParticleSectors = new int[NbrLzSectors];
		  for (int i = 0; i < NbrLzSectors; ++i)
		    {
		      NbrParticleSectors[i] = NbrParticles1;
		    }
		}
	      for (int i = 0; i < NbrLzSectors; ++i)
		{
		  cout << "state generated in the N=" << NbrParticleSectors[i] << " 2*Lz=" << LzSectors[i] << " sector" << endl;
		  if ((Manager.GetBoolean("conserve-particles") == false) || (NbrParticleSectors[i] == NbrParticles1))
		    {
		      LongRationalVector& RationalOutputState = RationalOutputStates[i];
		      bool zeroFlag = true;
		      int RootPosition = 0;
		      for (int TmpPos = 0; (TmpPos < RationalOutputState.GetVectorDimension()) && (zeroFlag == true); ++TmpPos)
			{
			  if (RationalOutputState[TmpPos].IsZero() == false)
			    zeroFlag = false;
			  else
			    ++RootPosition;
			}
		      if (zeroFlag)
			{
			  cout << "this state is null." << endl;
			}
		      else
			{
			  LongRational RootCoef = RationalOutputState[RootPosition];
			  RationalOutputState /= RootCoef; 
			  
			  if (Manager.GetBoolean("subset-symmetrization") == false)
			    {
			    if (Manager.GetString("output-file") == 0)
			      {
				if (Manager.GetBoolean("unnormalized-basis") == false)
				  {
				    sprintf(FullOutputFileName, "bosons_rational_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					    TargetNbrFluxQuanta, LzSectors[i]);
				  }
				else
				  {
				    sprintf(FullOutputFileName, "bosons_unnormalized_rational_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					  TargetNbrFluxQuanta, LzSectors[i]);
				  }
			      }
			    else
			      {
				sprintf(FullOutputFileName, "%s_lz_%d.0.vec", Manager.GetString("output-file"), LzSectors[i]);
			      }
			    }
			  else
			    {
			      if (Manager.GetString("output-file") == 0)
				{
				  if (Manager.GetBoolean("unnormalized-basis") == false)
				    {
				      sprintf(FullOutputFileName, "bosons_rational_subset_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					      TargetNbrFluxQuanta, LzSectors[i]);
				}
				  else
				    {
				  sprintf(FullOutputFileName, "bosons_unnormalized_rational_subset_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					  TargetNbrFluxQuanta, LzSectors[i]);
				    }
				}
			      else
				{
				  sprintf(FullOutputFileName, "%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("output-file"), NbrParticleSectors[i], 
					  TargetNbrFluxQuanta, LzSectors[i]);
				}
			    }
			}
		  
		      if (RationalOutputState.WriteVector(FullOutputFileName) == false)
			{
			  cout << "error while writing output state " << FullOutputFileName << endl;
			  return -1;
			}
		      ++NbrGeneratedStates;
		    }		  
		}
	      cout << "Symmetrization has generated " << NbrGeneratedStates << " state(s)" << endl;
	    }
	}
      else
	{
	  char* FullOutputFileName = new char [strlen(OutputFileName) + 128];
	  RealVector OutputState;
	  LongRationalVector RationalOutputState;
	  if (Manager.GetBoolean("rational") == false)
	    {
	      cout << "Single state symmetrization is not implemented for femionic non rational states" << endl;
	      return -1;
	    }
	  else
	    {
	      LongRationalVector* RationalOutputStates;
	      int* LzSectors;
	      int* NbrParticleSectors = 0;
	      int NbrLzSectors;
	      int TargetNbrFluxQuanta = 0;
	      if (Manager.GetBoolean("periodicity-symmetrization") == false)
		{
		  if (Manager.GetBoolean("subset-symmetrization") == true)
		    {
		      NbrLzSectors = FermionicSpace1->SymmetrizeSingleStatePeriodicSubsetOrbitals(RationalState1, Manager.GetInteger("subset-shift"),
												  Manager.GetInteger("nbr-orbitals"),
												  RationalOutputStates, NbrParticleSectors, LzSectors);
		      TargetNbrFluxQuanta = ((NbrFluxQuanta1 + 1) / Manager.GetInteger("nbr-orbitals")) - 1;
		      if ((((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals")) != 0) && (Manager.GetInteger("subset-shift") < ((NbrFluxQuanta1 + 1) % Manager.GetInteger("nbr-orbitals"))))
			TargetNbrFluxQuanta += 1;
		    }
		  else
		    {
		      
		      NbrLzSectors = FermionicSpace1->SymmetrizeSingleStateOneIntoManyOrbital(RationalState1, (int) Manager.GetInteger("nbr-orbitals"), RationalOutputStates, LzSectors);
		      TargetNbrFluxQuanta = (((NbrFluxQuanta1 + 1) / ((int) Manager.GetInteger("nbr-orbitals"))) - 1);
		    }
		}
	      else
		{
		  // 		  NbrLzSectors = FermionicSpace1->SymmetrizeSingleStatePeriodicOrbitals(RationalState1, (NbrFluxQuanta1 + 1) / (int) Manager.GetInteger("nbr-orbitals"), RationalOutputStates,
		  // 											LzSectors);
		  cout << "--periodicity-symmetrization is not implemented for fermions" << endl;
		  return 0;
		  TargetNbrFluxQuanta = (((NbrFluxQuanta1 + 1) / ((int) Manager.GetInteger("nbr-orbitals"))) - 1);
		}
	      int NbrGeneratedStates = 0;
	      if ((NbrLzSectors > 0) && (NbrParticleSectors == 0))
		{
		  NbrParticleSectors = new int[NbrLzSectors];
		  for (int i = 0; i < NbrLzSectors; ++i)
		    {
		      NbrParticleSectors[i] = NbrParticles1;
		    }
		}
	      for (int i = 0; i < NbrLzSectors; ++i)
		{
		  cout << "state generated in the N=" << NbrParticleSectors[i] << " 2*Lz=" << LzSectors[i] << " sector" << endl;
		  if (((Manager.GetBoolean("conserve-particles") == false) || (NbrParticleSectors[i] == NbrParticles1))
		      && ((Manager.GetBoolean("fix-totallz") == false) || (Manager.GetInteger("total-lz") == LzSectors[i])))
		    {
		      LongRationalVector& RationalOutputState = RationalOutputStates[i];
		      bool zeroFlag = true;
		      int RootPosition = 0;
		      for (int TmpPos = 0; (TmpPos < RationalOutputState.GetVectorDimension()) && (zeroFlag == true); ++TmpPos)
			{
			  if (RationalOutputState[TmpPos].IsZero() == false)
			    zeroFlag = false;
			  else
			    ++RootPosition;
			}
		      if (zeroFlag)
			{
			  cout << "this state is null." << endl;
			}
		      else
			{
			  LongRational RootCoef = RationalOutputState[RootPosition];
			  RationalOutputState /= RootCoef; 
			  if (Manager.GetBoolean("subset-symmetrization") == false)
			    {
			      if (Manager.GetString("output-file") == 0)
				{
				  if (Manager.GetBoolean("unnormalized-basis") == false)
				    {
				      sprintf(FullOutputFileName, "fermions_rational_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					      TargetNbrFluxQuanta, LzSectors[i]);
				    }
				  else
				    {
				      sprintf(FullOutputFileName, "fermions_unnormalized_rational_symmetrized_n_%d_2s_%d_lz_%d.0.vec", NbrParticleSectors[i], 
					      TargetNbrFluxQuanta, LzSectors[i]);
				    }
				}
			      else
				{
				  sprintf(FullOutputFileName, "%s_lz_%d.0.vec", Manager.GetString("output-file"), LzSectors[i]);
				}
			    }
			  else
			    {
			      if (Manager.GetString("output-file") == 0)
				{
				  sprintf(FullOutputFileName, "fermions_unnormalized_rational_subset_%d_%d_n_%d_2s_%d_lz_%d.0.vec", 
					  (int) Manager.GetInteger("subset-shift"), (int) Manager.GetInteger("nbr-orbitals"), NbrParticleSectors[i], TargetNbrFluxQuanta, LzSectors[i]);
				}
			      else
				{
				  sprintf(FullOutputFileName, "%s_n_%d_2s_%d_lz_%d.0.vec", Manager.GetString("output-file"), NbrParticleSectors[i], 
					  TargetNbrFluxQuanta, LzSectors[i]);
				}
			    }
			  if (RationalOutputState.WriteVector(FullOutputFileName) == false)
			    {
			      cout << "error while writing output state " << FullOutputFileName << endl;
			      return -1;
			    }
			  ++NbrGeneratedStates;
			}
		    }
		}
	      cout << "Symmetrization has generated " << NbrGeneratedStates << " state(s)" << endl;
	    }
	}
    }
  
  if (Space1 != 0)
    delete Space1;
  if (Space2 != 0)
    delete Space2;
  if (TargetSpace != 0)
    delete TargetSpace;
}

