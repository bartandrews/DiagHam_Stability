#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereWithSU4SpinAllEntanglement.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

//Permutations routine I need to test and its dependents
#include "GeneralTools/SmallIntegerArray.h"
#include "GeneralTools/List.h"
#include "GeneralTools/OrderedList.h"
#include "GeneralTools/Permutations.h"

#include "MathTools/FactorialCoefficient.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonsWithSpinLandauLevelLiftOperation.h"

#include <vector>
#include <set>

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
 
  OptionManager Manager ("FQHESphereBosonicStateLandauLevelLift" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
  ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += MiscGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
	
  (*SystemGroup) += new SingleStringOption ('\0', "state", "name of a an SU4 spin, 2N particle bosonic state");
  (*SystemGroup) += new SingleStringOption ('\n', "su2-state", "name of a an SU2 spin, N particle bosonic state (overrides su4)");
  (*SystemGroup) += new SingleStringOption ('p', "polarized-state", "name of a polarized N particle bosonic state");
  
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files, default is nbody)");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isospin", "twice the total isospin value of the system (0 if it has to be guessed from file name)", 0);
  
  (*SystemGroup) += new BooleanOption ('\n', "moore-read-fermions", "multiply by Laughlin factor to find state for fermionic Moore-Read 5/2 state", false);
  (*SystemGroup) += new SingleStringOption('\n', "jastrow-factor", "name of jastrow factor with which to multiply bosonic state to obtain fermionic 5/2 excitation");
  (*SystemGroup) += new SingleStringOption  ('\n', "resume-file", "use this file as the partial vector to resume from");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "resume-idx", "use this file as the partial vector to resume from", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "mpi-stages", "the number of stages divide into when using MPI  (default is 20)", 20);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "smp-stages", "the number of stages divide into when using SMP  (default is 20)", 20);
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "the output vector will be written in the normalized basis",false);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the product state into a binary file","product.vec");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "(optionally) output the product state into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "verbose", "print more output on screen");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);

  char* OutputFileName =  Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
 
  RealVector InitialState;

  bool SU2Mode=false;

  if ((Manager.GetString("state")==0)&&(Manager.GetString("su2-state")==0))
    {
      cout << "either state or su2-state need to be given as an input"<<endl;
      return -1;
    }
  if (Manager.GetString("su2-state")!=0)
    {
      SU2Mode = true;
      if (InitialState.ReadVector(Manager.GetString("su2-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("state") << endl;
	  return -1;
	}
    }
  else if (InitialState.ReadVector(Manager.GetString("state")) == false)
    {
      cout << "error while reading " << Manager.GetString("state") << endl;
      return -1;
    }

  RealVector PolarizedState;
  if (PolarizedState.ReadVector(Manager.GetString("polarized-state")) == false)
    {
      cout << "error while reading " << Manager.GetString("polarized-state") << endl;
      return -1;
    }
  //cout << PolarizedState << "\n";

  
  int NbrParticles = 0, NbrParticles2=0;
  int LzMaxIn = 0;
  int TotalLzIn = 0;
  int TotalSzIn = 0;
  int TotalIzIn = 0;
  int TotalEzIn = -1; // do not query entanglement
  bool FermionFlag = false;

  if (SU2Mode)
    {
      if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(Manager.GetString("su2-state"), NbrParticles, LzMaxIn, TotalLzIn, TotalSzIn, FermionFlag) == false)
	{
	  cout << "Error: Could not read su2 system information from input file name"<<endl;
	  return -1;
	}
      if (FermionFlag==true)
	{
	  cout << "Error: bosonic input file required!"<<endl;
	  return -1;
	}
    }
  else
    {
      if (FQHEOnSphereWithSU4SpinFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMaxIn, TotalLzIn, TotalSzIn, TotalIzIn, TotalEzIn, FermionFlag) == false)
	{
	  cout << "Error: Could not read su4 system information from input file name"<<endl;
	  return -1;
	}
      if (FermionFlag==true)
	{
	  cout << "Error: bosonic input file required!"<<endl;
	  return -1;
	}
      if (TotalIzIn!=0)
	{
	  cout << "Require states with Iz = 0, or equal populations in plus and minus for the moment."<<endl;
	  return -1;
	}

      cout << "Read system parameters: N="<<NbrParticles<<", LzMax="<<LzMaxIn<<", Lz="<<TotalLzIn<<", Sz="<<TotalSzIn<<", Iz="<<TotalIzIn<<endl;
    }


  int LzMaxPo = 0;
  int TotalLzPo = 0;

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("polarized-state"), NbrParticles2, LzMaxPo, TotalLzPo, FermionFlag) == false)
    {
      cout << "Error: Could not read system information from polarized file name" << endl;
      return -1;
    }
  if (FermionFlag==true)
    {
      cout << "Error: bosonic input file required!"<<endl;
      return -1;
    }

  if (SU2Mode)
    {
      if (NbrParticles2!=NbrParticles)
	{
	  cout << "Error: Initial state and polarized state need to have the same number of particles (su2-mode)"<<endl;
	  return -1;
	}
      cout << "LzMaxPo " << LzMaxPo << " TotalLzPo " << TotalLzPo << "\n";  

      BosonOnSphereWithSpin * InitialSpace;
      InitialSpace = new BosonOnSphereWithSpin (NbrParticles,TotalLzIn,LzMaxIn,TotalSzIn);
  
      if (InitialSpace->GetHilbertSpaceDimension() != InitialState.GetVectorDimension())
	{
	  cout << "dimension mismatch between the initial state (" << InitialState.GetVectorDimension() << ") and the Hilbert space (" << InitialSpace->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
  
      BosonOnSphereShort * PolarizedSpace;
      PolarizedSpace = new BosonOnSphereShort (NbrParticles,TotalLzPo,LzMaxPo);
      
      if (PolarizedSpace->GetHilbertSpaceDimension() != PolarizedState.GetVectorDimension())
	{
	  cout << "dimension mismatch between the polarized state (" << PolarizedState.GetVectorDimension() << ") and the Hilbert space (" << PolarizedSpace->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}

  
      int TotalLzFi = TotalLzIn + TotalLzPo;
      int LzMaxFi = LzMaxIn + LzMaxPo;
      BosonOnSphereWithSpin * FinalSpace;

      FinalSpace = new BosonOnSphereWithSpin (NbrParticles,TotalLzFi,LzMaxFi,TotalSzIn);

    
      RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
      int ResumeIdx = Manager.GetInteger("resume-idx");
		
      if (Manager.GetString("resume-file") != 0 )
	{
	  if ( Architecture.GetArchitecture()->ReadVector(OutputVector, Manager.GetString("resume-file")) == false )
	    {
	      cout << "error while reading " << Manager.GetString("resume-file") << endl;
	      return -1;
	    }		
	}
  

      // process operation tasks:
  


      cout << "LLL Polarized Space LzMax " << PolarizedSpace->GetLzValue() << "\n";

      cout << "Performing lift\n";
      //InitialSpace->BosonicStateWithSpinTimesBosonicState( InitialState, PolarizedState, OutputVector, PolarizedSpace, 0, InitialSpace->GetHilbertSpaceDimension(), FinalSpace);
      FQHESphereBosonsWithSpinLandauLevelLiftOperation MainOperation(InitialSpace, PolarizedSpace, FinalSpace, &InitialState, &PolarizedState, &OutputVector, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
      MainOperation.ApplyOperation(Architecture.GetArchitecture());

      // normalize output vector
      OutputVector.Normalize();

      if (Manager.GetBoolean("verbose"))
	{
	  cout << "Result="<<endl<<OutputVector << endl;

	  for(int i=0; i<FinalSpace->GetHilbertSpaceDimension(); i++)
	    {
	      FinalSpace->PrintState(cout, i);
	      cout << "\n";
	    }
	  cout << "\n";
	}
  
      if(OutputTxtFileName!=0)
	{
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputVector[i] << " ";
	      FinalSpace->PrintStateMonomial(File, i) << endl;
	      
	    }
	  File.close();
	}
      cout << "Result written to "<<Manager.GetString("bin-output")<<endl;
      Architecture.GetArchitecture()->WriteVector(OutputVector, Manager.GetString("bin-output"));
      
      delete InitialSpace;
      delete PolarizedSpace;
      delete FinalSpace;
      return 0; // end of SU2-mode here.
    }
  
  // full su4 lift
  if (2*NbrParticles2!=NbrParticles)
    {
      cout << "Error: Initial state must have twice the number of particles of the polarized state"<<endl;
      return -1;
    }
  cout << "LzMaxPo " << LzMaxPo << " TotalLzPo " << TotalLzPo << "\n";

  BosonOnSphereWithSU4SpinAllEntanglement * InitialSpace;

  cout << "Initialising the SU4 state with " << NbrParticles << " particles " << TotalLzIn << " TotalLzIn " << LzMaxIn << " LzMaxIn " << TotalSzIn << " TotalSzIn " << TotalIzIn << " TotalIzIn " << "\n";
  
  InitialSpace = new BosonOnSphereWithSU4SpinAllEntanglement (NbrParticles,TotalLzIn,LzMaxIn,TotalSzIn,TotalIzIn);
  if (InitialSpace->GetHilbertSpaceDimension() != InitialState.GetVectorDimension())
    {
      cout << "dimension mismatch between the initial state (" << InitialState.GetVectorDimension() << ") and the Hilbert space (" << InitialSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  BosonOnSphereShort * PolarizedSpace;
  PolarizedSpace = new BosonOnSphereShort (NbrParticles2,TotalLzPo,LzMaxPo);

  if (PolarizedSpace->GetHilbertSpaceDimension() != PolarizedState.GetVectorDimension())
    {
      cout << "dimension mismatch between the polarized state (" << PolarizedState.GetVectorDimension() << ") and the Hilbert space (" << PolarizedSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  
  int TotalLzFi = TotalLzIn + TotalLzPo;
  int LzMaxFi = LzMaxIn + LzMaxPo;
  BosonOnSphereWithSpin * FinalSpace;

  FinalSpace = new BosonOnSphereWithSpin (NbrParticles,TotalLzFi,LzMaxFi,TotalSzIn);

    
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
  int ResumeIdx = Manager.GetInteger("resume-idx");
		
  if (Manager.GetString("resume-file") != 0 )
    {
      if ( Architecture.GetArchitecture()->ReadVector(OutputVector, Manager.GetString("resume-file")) == false )
	{
	  cout << "error while reading " << Manager.GetString("resume-file") << endl;
	  return -1;
	}		
    }
  
  // process operation tasks:
  
  // FQHESphereBosonsWithSpinLandauLevelLiftOperation MainOperation(InitialSpace, PolarizedSpace, FinalSpace, &InitialState, &PolarizedState, &OutputVector, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
  //MainOperation.ApplyOperation(Architecture.GetArchitecture());    


  int NbrParticlesI = NbrParticles2;
  int NbrParticlesII = NbrParticles2;

  /*
  cout << "InitialSpace declared with states \n";
  for(int i=0; i<InitialSpace->GetHilbertSpaceDimension(); i++)
    {
      InitialSpace->PrintState(cout, i);
      cout << "\n";
    }

  */

  /*
  //bodge in romers state here to test
  InitialState[0] = 0.0;
  InitialState[1] = 2.0;
  InitialState[2] = -sqrt(2.0);
  InitialState[3] = -sqrt(2.0);
  InitialState[4] = 0.0;
  InitialState[5] = 2.0;
  */


  int TotalLzPlus, TotalLzMinus;
  int TotalSzPlus, TotalSzMinus;


  int NbrTotalLzGroupValues;
  int NbrTotalSzGroupValues;
  int MaxParticlesInGroup;
  if(NbrParticlesI > NbrParticlesII)
    {
      MaxParticlesInGroup = NbrParticlesI;
    }
  else
    {
      MaxParticlesInGroup = NbrParticlesII;
    }
 
  NbrTotalLzGroupValues = 2*MaxParticlesInGroup*LzMaxIn+1;
  NbrTotalSzGroupValues = 2*MaxParticlesInGroup+1;
  std::set< unsigned long > ** UniqueStatesToLiftAllLzSz = new std::set< unsigned long > * [NbrTotalLzGroupValues];
  std::vector< RealVector > ** LiftedStatesAllLzSz = new std::vector< RealVector > * [NbrTotalLzGroupValues];

  std::set< unsigned long > * UniqueStatesToLift;
  std::vector< RealVector > * LiftedStates;

  for(int lz=0; lz<NbrTotalLzGroupValues; lz++)
    {
      UniqueStatesToLiftAllLzSz[lz] = new std::set<unsigned long> [NbrTotalSzGroupValues];
      LiftedStatesAllLzSz[lz] = new std::vector< RealVector > [NbrTotalSzGroupValues];
    }

  std::pair< std::set< unsigned long >::iterator, bool > insertReturnValue;
  
  unsigned long * tempStatePlus = new unsigned long [2*(LzMaxIn+1)];
  unsigned long * tempStateMinus = new unsigned long [2*(LzMaxIn+1)];
  
  for(int i=0; i<InitialSpace->GetHilbertSpaceDimension(); i++)
    {
      if (fabs(InitialState[i]) < 1e-16)
	{
	  for(int lz=0; lz<=LzMaxIn; lz++)
	    {
	      tempStatePlus[lz] = 0x0ul;
	      tempStatePlus[lz+LzMaxIn+1] = 0x0ul;
	      tempStateMinus[lz] = 0x0ul;
	      tempStateMinus[lz+LzMaxIn+1] = 0x0ul;
	    }

	  unsigned long tempFermionPlus, tempFermionMinus;
	  unsigned long * tempStateUpPlus = &tempStatePlus[LzMaxIn+1];
	  unsigned long * tempStateDownPlus = &tempStatePlus[0];
	  unsigned long * tempStateUpMinus = &tempStateMinus[LzMaxIn+1];
	  unsigned long * tempStateDownMinus = &tempStateMinus[0];
	  
	  cout << "Calculating lift of group state \n";
	  InitialSpace->PrintState(cout, i);
	  cout << "\n";

	  InitialSpace->GetTotalLz(i, TotalLzPlus, TotalLzMinus);
	  InitialSpace->GetTotalSz(i, TotalSzPlus, TotalSzMinus);
	  
	  InitialSpace->GetBosonicStateDescription(i, tempStatePlus, tempStateMinus);
	  InitialSpace->GetFermionicStateDescription(i, tempFermionPlus, tempFermionMinus);	  

	  cout << "Declaring initial Hilbert spaces for each group\n";
	  cout << "Plus space TotalLz " << TotalLzPlus << " TotalSz " << TotalSzPlus << "\n";
	  cout << "Minus space TotalSz " << TotalLzMinus << " TotalSz " << TotalSzMinus << "\n";
	  BosonOnSphereWithSpin * InitialPlusSpace = new BosonOnSphereWithSpin(NbrParticlesI, TotalLzPlus, LzMaxIn, TotalSzPlus);
	  BosonOnSphereWithSpin * InitialMinusSpace = new BosonOnSphereWithSpin(NbrParticlesII, TotalLzMinus, LzMaxIn, TotalSzMinus);

	 	  
	  cout << "Declaring final Hilbert spaces for each group each with LzMax " << LzMaxIn + LzMaxPo << "\n";
	  cout << "Plus space TotalLz " << TotalLzPlus + TotalLzPo << " TotalSz " << TotalSzPlus << "\n";
	  cout << "Minus space TotalLz " << TotalLzMinus + TotalLzPo << " TotalSz " << TotalSzMinus << "\n";

	  BosonOnSphereWithSpin * FinalPlusSpace = new BosonOnSphereWithSpin(NbrParticlesI, TotalLzPlus + TotalLzPo, LzMaxIn + LzMaxPo, TotalSzPlus);
	  BosonOnSphereWithSpin * FinalMinusSpace = new BosonOnSphereWithSpin(NbrParticlesII, TotalLzMinus + TotalLzPo, LzMaxIn + LzMaxPo, TotalSzMinus);
	  
	  RealVector InitialPlusState(InitialPlusSpace->GetHilbertSpaceDimension(), true);
	  RealVector InitialMinusState(InitialMinusSpace->GetHilbertSpaceDimension(), true);
	  
	  cout << "InitialPlusState dimension " << InitialPlusSpace->GetHilbertSpaceDimension() << "\n";
	  cout << "InitialMinusState dimension " << InitialMinusSpace->GetHilbertSpaceDimension() << "\n";

	  cout << "tempStatePlus \n";
	  for(int lz=0; lz<=LzMaxIn; lz++)
	    {
	      cout << tempStateUpPlus[lz] << "u " << tempStateDownPlus[lz] << "d |";
	    }
	  cout << "\n";

	  cout << "tempStateMinus \n";
	  for(int lz=0; lz<=LzMaxIn; lz++)
	    {
	      cout << tempStateUpMinus[lz] << "u " << tempStateDownMinus[lz] << "d |";
	    }
	  cout << "\n";

	  int InitialPlusIndex = InitialPlusSpace->FindStateIndex( tempStateUpPlus, tempStateDownPlus );
	  int InitialMinusIndex = InitialMinusSpace->FindStateIndex( tempStateUpMinus, tempStateDownMinus );

	  /*
	  cout << "Index of plus state " << InitialPlusIndex << "\n";
	  cout << "Index of minus state " << InitialMinusIndex  << "\n";
	  */

	  InitialPlusState[ InitialPlusIndex ] = 1.0;
	  InitialMinusState[ InitialMinusIndex ] = 1.0;
	  
	  RealVector FinalPlusState( FinalPlusSpace->GetHilbertSpaceDimension(), true);
	  RealVector FinalMinusState( FinalMinusSpace->GetHilbertSpaceDimension(), true);
	  
	  int shiftedTotalLzPlus = TotalLzPlus + MaxParticlesInGroup*LzMaxIn;
	  int shiftedTotalSzPlus = TotalSzPlus + MaxParticlesInGroup;
	  int shiftedTotalLzMinus = TotalLzMinus + MaxParticlesInGroup*LzMaxIn;
	  int shiftedTotalSzMinus = TotalSzMinus + MaxParticlesInGroup;

	  //cout << "shiftedTotalLzPlus " << shiftedTotalLzPlus << " NbrTotalLzGroupValues " << NbrTotalLzGroupValues << " shiftedTotalSzPlus " << shiftedTotalSzPlus << " NbrTotalSzGroupValues " << NbrTotalSzGroupValues << " shiftedTotalLzMinus " << shiftedTotalLzMinus << "shiftedTotalSzMinus " << shiftedTotalSzMinus << "\n";

	  UniqueStatesToLift = &( UniqueStatesToLiftAllLzSz[shiftedTotalLzPlus][shiftedTotalSzPlus] );
	  LiftedStates = &( LiftedStatesAllLzSz[shiftedTotalLzPlus][shiftedTotalSzPlus] );

	  insertReturnValue = UniqueStatesToLift->insert(tempFermionPlus);
	  
	  if( insertReturnValue.second == true )
	    {
	      //new state to calculate lift
	      cout << "New state to calculate lift\n";
	      InitialPlusSpace->BosonicStateWithSpinTimesBosonicState(InitialPlusState, PolarizedState, FinalPlusState, PolarizedSpace, InitialPlusIndex, 1, FinalPlusSpace);

	      //FQHESphereBosonsWithSpinLandauLevelLiftOperation MainOperation(InitialPlusSpace, PolarizedSpace, FinalPlusSpace, &InitialPlusState, &PolarizedState, &FinalPlusState, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
	      //MainOperation.ApplyOperation(Architecture.GetArchitecture());

	      
	      /*
	      cout << "Lifted plus state\n";
	      for(int plusStateIndex=0; plusStateIndex<FinalPlusSpace->GetHilbertSpaceDimension(); plusStateIndex++)
		cout << FinalPlusState[plusStateIndex] << " ";
	      cout << "\n";
	      */
	      //save lifted state 
	      LiftedStates->insert( LiftedStates->begin() + std::distance(UniqueStatesToLift->begin(), insertReturnValue.first), FinalPlusState );
	      /*
	      cout << "Lifted states \n";
	      for(int tableIndex=0; tableIndex<LiftedStates->size(); tableIndex++)
		cout << tableIndex << "\n" << (*LiftedStates)[tableIndex] << "\n";
	      */
	    }
	  else
	    {
	      //look up result of lift
	      cout << "Looking up lift of state\n";
	      InitialPlusSpace->PrintState( cout, InitialPlusIndex );
	      cout << " fermion rep " << tempFermionPlus << " at position ";
	      cout << std::distance( UniqueStatesToLift->begin(), insertReturnValue.first ) << " in the table\n";
	      /*
	      for(std::set< unsigned long  >::iterator stateToLift = UniqueStatesToLift->begin(); stateToLift != UniqueStatesToLift->end(); stateToLift++)
		{
		  cout << (*stateToLift) << "\n";
		}
	      for(int tableIndex=0; tableIndex < LiftedStates->size(); tableIndex++)
		{
		  cout << tableIndex << "\n" << (*LiftedStates)[tableIndex] << "\n";
		}
	      */
	      FinalPlusState = (*LiftedStates)[ std::distance(UniqueStatesToLift->begin(), insertReturnValue.first) ];
	    }
	  
	  UniqueStatesToLift = &( UniqueStatesToLiftAllLzSz[shiftedTotalLzMinus][shiftedTotalSzMinus] );
	  LiftedStates = &( LiftedStatesAllLzSz[shiftedTotalLzMinus][shiftedTotalSzMinus] );

	  insertReturnValue = UniqueStatesToLift->insert(tempFermionMinus);
	  
	  if( insertReturnValue.second == true )
	    {
	      //new state to calculate lift
	      cout << "New state to lift\n";
	      InitialMinusSpace->BosonicStateWithSpinTimesBosonicState(InitialMinusState, PolarizedState, FinalMinusState, PolarizedSpace, InitialMinusIndex, 1, FinalMinusSpace);
	      cout << "Lifted minus state\n";
	      /*
	      for(int minusStateIndex=0; minusStateIndex<FinalMinusSpace->GetHilbertSpaceDimension(); minusStateIndex++)
		cout << FinalMinusState[minusStateIndex] << " ";
	      cout << "\n";
	      */
	      //save lifted state
	      LiftedStates->insert( LiftedStates->begin() + std::distance(UniqueStatesToLift->begin(), insertReturnValue.first), FinalMinusState );
	      /*
	      cout << "LiftedStates\n";
	      for(int tableIndex=0; tableIndex < LiftedStates->size(); tableIndex++)
		cout << tableIndex << "\n" << (*LiftedStates)[tableIndex] << "\n";
	      */
	    }
	  else 
	    {
	      //look up result of lift
	      cout << "Looking up lift of state\n";
	      InitialMinusSpace->PrintState( cout, InitialMinusIndex );
	      cout << " fermion rep " << tempFermionMinus << " at position ";
	      cout << std::distance( UniqueStatesToLift->begin(), insertReturnValue.first ) << " in the table \n";
	      /*
	      for(std::set< unsigned long  >::iterator stateToLift = UniqueStatesToLift->begin(); stateToLift != UniqueStatesToLift->end(); stateToLift++)
		{
		  cout << (*stateToLift) << "\n";
		}
	      for(int tableIndex=0; tableIndex < LiftedStates->size(); tableIndex++)
		{
		  cout << tableIndex << "\n" << (*LiftedStates)[tableIndex] << "\n";
		}
	      */
	      FinalMinusState = (*LiftedStates)[ std::distance(UniqueStatesToLift->begin(), insertReturnValue.first) ];
	    }
	  
	  unsigned long * FinalPlusStateUp = new unsigned long [LzMaxFi+1];
	  unsigned long * FinalPlusStateDown = new unsigned long [LzMaxFi+1];
	  unsigned long * FinalMinusStateUp = new unsigned long [LzMaxFi+1];
	  unsigned long * FinalMinusStateDown = new unsigned long [LzMaxFi+1];
	  
	  for(int plusIndex = 0; plusIndex < FinalPlusSpace->GetHilbertSpaceDimension(); plusIndex++)
	    {
	      FinalPlusSpace->GetBosonicDescription(plusIndex, FinalPlusStateUp, FinalPlusStateDown);
	      for(int minusIndex = 0; minusIndex < FinalMinusSpace->GetHilbertSpaceDimension(); minusIndex++)
		  {
		    FinalMinusSpace->GetBosonicDescription(minusIndex, FinalMinusStateUp, FinalMinusStateDown);
		    /*
		    cout << "Final state pre-symmetrisation\n";
		    for(int lz=0; lz<=LzMaxFi; lz++)
		      {
			cout << FinalPlusStateUp[lz] << "u+ " << FinalPlusStateDown[lz] << "d+ " << FinalMinusStateUp[lz] << "u- " << FinalMinusStateDown[lz] << "d-| ";
		      }
		    */
		    cout << "\n";
		    cout << "Symmetrising over groups\n";
		    cout << "FinalPlusState coefficient " << FinalPlusState[plusIndex] << "\n";
		    cout << "FinalMinusState coefficient " << FinalMinusState[minusIndex] << "\n";
		    cout << "InitialState component " << InitialState[i] << "\n";
		    FinalSpace->SymmetriseOverGroupsAndAddToVector(FinalPlusStateUp, FinalMinusStateUp, FinalPlusStateDown, FinalMinusStateDown, InitialState[i]*FinalPlusState[plusIndex]*FinalMinusState[minusIndex], OutputVector);
		    cout << "Symmetrised over term\n";
		  }
	    }

	  delete [] FinalPlusStateUp;
	  delete [] FinalPlusStateDown;
	  delete [] FinalMinusStateUp;
	  delete [] FinalMinusStateDown;
	  
	  delete InitialPlusSpace;
	  delete InitialMinusSpace;
	  delete FinalPlusSpace;
	  delete FinalMinusSpace;
	}
    }

  delete [] tempStatePlus;
  delete [] tempStateMinus;


  if(Manager.GetBoolean("moore-read-fermions"))
    {
      cout << "Multiplying by Jastrow factor to map onto excitation over the Moore-Read 5/2 state\n";
      RealVector JastrowState;
      if (JastrowState.ReadVector(Manager.GetString("jastrow-factor")) == false)
	{
	  cout << "error while reading " << Manager.GetString("jastrow-factor") << endl;
	  return -1;
	}
      
      int NbrJastrowParticles=0;
      int LzMaxJas = 0;
      int TotalLzJas = 0;
      bool JastrowFermionFlag;

      if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("jastrow-factor"), NbrJastrowParticles, LzMaxJas, TotalLzJas, JastrowFermionFlag) == false)
	{
	  cout << "Error: Could not read system information from polarized file name" << endl;
	  return -1;
	}
      if (JastrowFermionFlag==true)
	{
	  cout << "Error: bosonic input file required!"<<endl;
	  return -1;
	}    
      
      if(NbrJastrowParticles != NbrParticles)
	{
	  cout << "Error: Jastrow factor has " << NbrJastrowParticles << " != " << NbrParticles << " particles of texture\n";
	  return -1;
	}

      BosonOnSphereShort * JastrowSpace = new BosonOnSphereShort(NbrParticles, TotalLzJas, LzMaxJas );

      int LzMaxFiveHalves = LzMaxFi + LzMaxJas;
      int TotalLzFiveHalves = TotalLzFi + TotalLzJas;

      BosonOnSphereWithSpin * FiveHalvesSpace = new BosonOnSphereWithSpin(NbrParticles, TotalLzFiveHalves, LzMaxFiveHalves, TotalSzIn);
      RealVector FiveHalvesOutputVector( FiveHalvesSpace->GetHilbertSpaceDimension(), true);
      //FinalSpace->BosonicStateWithSpinTimesBosonicState( OutputVector, JastrowState, FiveHalvesOutputVector, JastrowSpace, 0, FinalSpace->GetHilbertSpaceDimension(), FiveHalvesSpace );

      FQHESphereBosonsWithSpinLandauLevelLiftOperation MainOperation2(FinalSpace, JastrowSpace, FiveHalvesSpace, &OutputVector, &JastrowState, &FiveHalvesOutputVector, Manager.GetInteger("mpi-stages"), Manager.GetInteger("smp-stages"), ResumeIdx);
      MainOperation2.ApplyOperation(Architecture.GetArchitecture());
      
      
      if(Manager.GetBoolean("normalize"))
	FiveHalvesOutputVector.Normalize();
 
      Architecture.GetArchitecture()->WriteVector(FiveHalvesOutputVector, OutputFileName);

      ofstream File;
      if(OutputTxtFileName!=0)
	{
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < FiveHalvesSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << FiveHalvesOutputVector[i] << " ";
	      FiveHalvesSpace->PrintStateMonomial(File, i) << endl;
	    }
	}
      File.close();

      delete JastrowSpace;
      delete FiveHalvesSpace;

    }
  else
    {
      
      if(Manager.GetBoolean("normalize"))
	OutputVector.Normalize();
      
      cout << "Result of the Bosonic Landau lift:\n";
      for(int i=0; i<FinalSpace->GetHilbertSpaceDimension(); i++)
	{
	  cout << OutputVector[i] << " ";
	  FinalSpace->PrintState(cout, i);
	  cout << "\n";
	}
      
      Architecture.GetArchitecture()->WriteVector(OutputVector, OutputFileName);

      ofstream File;
      if(OutputTxtFileName!=0)
	{
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < FinalSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputVector[i] << " ";
	      FinalSpace->PrintStateMonomial(File, i) << endl;
	    }
	}
      File.close();
    

    }

  delete InitialSpace;

  for(int lz=0; lz<NbrTotalLzGroupValues; lz++)
    {
      delete [] UniqueStatesToLiftAllLzSz[lz];
      delete [] LiftedStatesAllLzSz[lz];
    }
  delete [] UniqueStatesToLiftAllLzSz;
  delete [] LiftedStatesAllLzSz;
  
  delete PolarizedSpace;
  delete FinalSpace;
  
      
}

  
  
