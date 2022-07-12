#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeKy.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHELatticeConvertSymmetrizedState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ProcessGroup = new OptionGroup ("process options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += ProcessGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");  
  (*SystemGroup) += new SingleIntegerOption  ('k', "ky", "constraint of momentum in y-direction (-1=all)", 0);

  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");

  (*ProcessGroup) += new SingleIntegerOption  ('s', "split", "split vector in multiple segments", 0);
  (*ProcessGroup) += new SingleIntegerOption  ('i', "segment-index", "number of segment to be calculated by this process", 0);
  
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the default name un-symmetrized.vec");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  if (Manager.GetString("input-file") == 0)
    {
      cout << "error, one input file should be provided. See man page for option syntax or type FQHELatticeConvertSymmetrizedState -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-file")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-file") << endl;
    }
  if ((Manager.GetInteger("split")!=0)&&(Manager.GetInteger("segment-index")>Manager.GetInteger("split")-1))
    {
      cout << "Segment index has to be chosen in interval [0,nbr-split]"<<endl;
      exit(1);
    }
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int Ky = Manager.GetInteger("ky");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = 9 << 20;

  char* InputFileName = Manager.GetString("input-file");
  bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");

  double Interaction=0.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  if (FQHEOnLatticeFindSystemInfoWithKyFromVectorFileName(InputFileName, NbrParticles, Lx, Ly, Ky, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }  
  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;

  if (Manager.GetBoolean("boson") == true)
    {
      if (Manager.GetBoolean("boson") == true)
	Statistics = false;
    }

  ComplexVector State;
  if (State.ReadVector (InputFileName) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-file") << endl;
      return -1;      
    }

  cout << "Using the following system info: "<<endl;
  cout << "N="<<NbrParticles<<", Lx="<<Lx<<", Ly="<<Ly<<", N_phi="<<NbrFluxQuanta<<", Ky="<<Ky<<endl;

  char *OutputName = Manager.GetString("output-file");
  if (OutputName==NULL)
    {
      OutputName = new char[strlen(Manager.GetString("input-file"))+20];
      if (Manager.GetInteger("split")!=0)	
	sprintf(OutputName,"%s.seg_%ld-%ld",Manager.GetString("input-file"),Manager.GetInteger("segment-index"),
		Manager.GetInteger("split"));
      else sprintf(OutputName,"%s.full",Manager.GetString("input-file"));
    }

  if (Statistics == true)
    {
      cout << "No fermionic systems implemented, yet..."<<endl;
      return -1;

    }
  else
    {
      ComplexVector OutputState;
      BosonOnLatticeKy InitialSpace(NbrParticles, Lx, Ly, Ky, NbrFluxQuanta, MemorySpace); 
      BosonOnLattice TargetSpace(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
      if (SymmetrizeFlag)
	{
	  if (TargetSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  cout << "Symmetrization not implemented, yet!"<<endl;
	  // OutputState = InitialSpace.ConvertToNBodyKyMomentumBasis(State, TargetSpace);
	}
      else
	{
	  if (InitialSpace.GetHilbertSpaceDimension() != State.GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  if (Manager.GetInteger("split")>0)
	    {
	      int SegSize=InitialSpace.GetHilbertSpaceDimension()/Manager.GetInteger("split");
	      int StartD=Manager.GetInteger("segment-index")*SegSize;
	      int StopD=(Manager.GetInteger("segment-index")+1)*SegSize;
	      if (Manager.GetInteger("segment-index")==Manager.GetInteger("split")-1)
		StopD=InitialSpace.GetHilbertSpaceDimension();
	      cout<<"Converting segment ["<<StartD<<", "<<StopD<<"]"<<endl;
	      OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace, StartD, StopD-StartD);
	    }
	  else
	    OutputState = InitialSpace.ConvertToNbodyBasis(State, TargetSpace);
	}
      if (OutputState.WriteVector(OutputName) == false)
	{
	  cout << "error while writing output state " << OutputName << endl;
	  return -1;
	}
    }
  delete [] OutputName;
}

