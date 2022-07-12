#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeKy.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
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
  OptionManager Manager ("FQHELatticeConvertSymmetrizedState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ProcessGroup = new OptionGroup ("process options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  ArchitectureManager Architecture;
  
  Manager += SystemGroup;
  Manager += ProcessGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  Architecture.AddOptionGroup(&Manager);
  (*SystemGroup) += new MultipleStringOption  ('\0', "input-file", "input state file name");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");  
  (*SystemGroup) += new SingleIntegerOption  ('k', "ky", "constraint of momentum in y-direction (-1=all)", 0);

  (*SystemGroup) += new BooleanOption  ('b', "boson", "use bosonic statistics (override autodetection from input file name)");
  (*SystemGroup) += new BooleanOption  ('r', "symmetrize", "symmetrize state (instead of unsymmetrizing it)");
  (*SystemGroup) += new BooleanOption  ('i', "inverse", "symmetrize state (instead of unsymmetrizing it)");
  (*ProcessGroup) += new SingleIntegerOption  ('\n', "nbr-component", "compute only a certain number of component", 0) ;
 
  (*ProcessGroup) += new SingleIntegerOption  ('s', "split", "split vector in multiple segments", 0) ;
  (*ProcessGroup) += new SingleIntegerOption  ('i', "segment-index", "number of segment to be calculated by this process", 0);
  
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the default name un-symmetrized.vec");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

if ((Manager.GetInteger("split")!=0)&&(Manager.GetInteger("segment-index")>Manager.GetInteger("split")-1))
    {
      cout << "Segment index has to be chosen in interval [0,nbr-split]"<<endl;
      exit(1);
    }
      int NbrVectors;
  char** VectorFiles = Manager.GetStrings("input-file", NbrVectors);

  int *NbrParticles = new int [NbrVectors];
  int *Lx = new int [NbrVectors];
  int *Ly = new int [NbrVectors];
  int Ky = Manager.GetInteger("ky");
  int *NbrFluxQuanta = new int [NbrVectors];
  unsigned long MemorySpace = 9 << 20;


  
  bool SymmetrizeFlag = Manager.GetBoolean("symmetrize");

  bool Statistics=false;
  bool HardCore=false;
  bool InverseFlag = Manager.GetBoolean("inverse");
  int NbrComponent = Manager.GetInteger("nbr-component");

  for(int i = 0; i < NbrVectors; i++)
    {
      NbrParticles[i] = 0;
      Lx[i] = 0;
      Ly[i]= 0;
      NbrFluxQuanta[i] = 0;
  if (FQHEOnLatticeFindSystemInfoFromFileName1(VectorFiles[i], NbrParticles[i], Lx[i], Ly[i], NbrFluxQuanta[i], Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }
    }
    
    for(int i = 1; i < NbrVectors; i++)
    {
      if((NbrParticles[i] != NbrParticles[i - 1]) || (Lx[i] != Lx[i - 1]) || (Ly[i] != Ly[i - 1]) || 
	 (NbrFluxQuanta[i] != NbrFluxQuanta[i - 1]))
	{ 
	  cout << "Vectors are not all in the same HilbertSpace." << endl;
	  return -1;
	}
    }
    
  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;

  if (Manager.GetBoolean("boson") == true)
    {
      if (Manager.GetBoolean("boson") == true)
	Statistics = false;
    }

  ComplexVector * InputState = new ComplexVector[NbrVectors];
  
  for(int i = 0; i < NbrVectors; i++)
    {
      if (IsFile(VectorFiles[i]) == false)
	{
	  cout << "state " << VectorFiles[i] << " does not exist or can't be opened" << endl;
	  return -1;
	}
      if(InputState[i].ReadVector(VectorFiles[i]) == false)
	{
	  cout << "error while reading " << VectorFiles[i] << endl;
	  return -1;
	}
    }

  cout << "Using the following system info: "<<endl;
  cout << "N="<<NbrParticles[0]<<", Lx="<<Lx[0]<<", Ly="<<Ly[0]<<", N_phi="<<NbrFluxQuanta[0]<<", Ky="<<Ky<<endl;

  char** VectorFilesOut = new char*[NbrVectors];
  VectorFilesOut[0] = Manager.GetString("output-file");
  if (VectorFilesOut[0]==NULL)
    {
      VectorFilesOut[0] = new char[strlen(VectorFiles[0])+20];
      if (Manager.GetInteger("split")!=0)	
	sprintf(VectorFilesOut[0],"%s.seg_%ld-%ld",VectorFiles[0],Manager.GetInteger("segment-index"),
		Manager.GetInteger("split"));
      else sprintf(VectorFilesOut[0],"%s.full",VectorFiles[0]);
      if(InverseFlag == true)
      {
	
      for(int i = 0; i < NbrVectors; i++)
	{
	  VectorFilesOut[i] = new char [strlen (VectorFiles[i])+ 32];
	  char* TmpPos = VectorFilesOut[i];
	  char* TmpPos2 = strstr(VectorFiles[i], "_q_");
	  strncpy (TmpPos,VectorFilesOut[i],TmpPos2-VectorFilesOut[i]);
	  sprintf (TmpPos, "_k_%d%s", Ky, TmpPos2);
	}
      }
    }

  if (Statistics == true)
    {
      cout << "No fermionic systems implemented, yet..."<<endl;
      return -1;

    }
  else
    {
      ComplexVector * OutputState = new ComplexVector[NbrVectors];
      BosonOnLatticeKy InitialSpace(NbrParticles[0], Lx[0], Ly[0], Ky, NbrFluxQuanta[0], MemorySpace); 
      BosonOnLattice TargetSpace (NbrParticles[0], Lx[0], Ly[0], NbrFluxQuanta[0], MemorySpace);
      if(NbrComponent == 0)
	NbrComponent = InitialSpace.GetHilbertSpaceDimension();
      if (SymmetrizeFlag)
	{
	  if (TargetSpace.GetHilbertSpaceDimension() != InputState[0].GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	  cout << "Symmetrization not implemented, yet!"<<endl;
	  // OutputState = InitialSpace.ConvertToNBodyKyMomentumBasis(State, TargetSpace);
	}
      else
	{
	  if(InverseFlag == true)
	  {
	    if (TargetSpace.GetHilbertSpaceDimension() != InputState[0].GetVectorDimension())
	    {
	      cout << "dimension mismatch between Hilbert space and input state" << endl;
	      return -1;
	    }
	    
	    OutputState = InitialSpace.ConvertFromNbodyBasis(InputState, TargetSpace,NbrVectors,NbrComponent,Architecture.GetArchitecture());
	  }
	  else
	  {
	  if (InitialSpace.GetHilbertSpaceDimension() != InputState[0].GetVectorDimension())
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
	      OutputState[0] = InitialSpace.ConvertToNbodyBasis(InputState[0], TargetSpace, StartD, StopD-StartD);
	    }
	  else
	    OutputState[0] = InitialSpace.ConvertToNbodyBasis(InputState[0], TargetSpace);
	  }
	}
	for(int i=0; i < NbrVectors; i++)
	{
      if (OutputState[i].WriteVector(VectorFilesOut[i]) == false)
	{
	  cout << "error while writing output state " << VectorFilesOut[i] << endl;
	  return -1;
	}
	}
    }
  delete [] VectorFilesOut;
}

