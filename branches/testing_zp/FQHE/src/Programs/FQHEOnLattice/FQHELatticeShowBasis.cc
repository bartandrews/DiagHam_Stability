#include "HilbertSpace/BosonOnLatticeKy.h"
#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "GeneralTools/FilenameTools.h"
#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "Options/Options.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;






int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeShowBasis" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "input-files", "input state file name(s)");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles of target many-body state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons (~Gutzwiller projection)");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");
  (*SystemGroup) += new BooleanOption('r',"first-real","Multiply each vector with a phase such that the first non-zero coefficient is real");
  (*SystemGroup) += new SingleIntegerOption('e',"real-element","Index of element to be made real with option 'first-real'",0);

  (*SystemGroup) += new SingleIntegerOption  ('k', "ky", "constraint of momentum in y-direction", 0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int Ky = Manager.GetInteger("ky");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = 9 << 20;


  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("input-files",NbrVectors);
  int VectorDimension=0;
  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  Complex *Phases = new Complex[NbrVectors];
  for (int i=0; i<NbrVectors; ++i) Phases[i]=1.0;

  double Interaction=-1.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  bool KySymmetry=false;
  if (NbrVectors > 0 )
    {
      if (FQHEOnLatticeFindSystemInfoWithKyFromVectorFileName(VectorFiles[0], NbrParticles, Lx, Ly, Ky, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	{
	  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[0], NbrParticles, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	    {
	      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
	      exit(1);
	    }
	  KySymmetry=false;
	  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
	  if (Manager.GetBoolean("no-hard-core"))
	    HardCore=false;
	}
      else
	{
	  KySymmetry=true;
	}
      
      bool tmpB, haveVector=false;
      for (int i=0; i<NbrVectors; ++i)
	{
	  if ((tmpB=Vectors[i].ReadVector(VectorFiles[i]))==false)
	    exit(1);
	  if (!haveVector)
	    VectorDimension=Vectors[i].GetVectorDimension();
	  if (haveVector && (Vectors[i].GetVectorDimension()!=VectorDimension))
	    {
	      cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
	      exit(1);
	    }
	  haveVector=haveVector | tmpB;
	}

      if (!haveVector)
	{
	  cout << "No valid vector files found: printing basis only!"<<endl;
	  NbrVectors=0;
	}
    }
  
  ParticleOnLattice *Space;
  if (KySymmetry)
    {
      cout << "N="<<NbrParticles<<", Lx="<<Lx<<", Ly="<<Ly<<", Q="<<NbrFluxQuanta<<", Ky="<<Ky<<endl;
      Space = new BosonOnLatticeKy(NbrParticles, Lx, Ly, Ky, NbrFluxQuanta, MemorySpace);      
    }
  else
    {
      cout << "N="<<NbrParticles<<", Lx="<<Lx<<", Ly="<<Ly<<", Q="<<NbrFluxQuanta<<endl;
      if (HardCore)
	Space =new HardCoreBosonOnLattice(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
      else Space = new BosonOnLattice(NbrParticles, Lx, Ly, NbrFluxQuanta, MemorySpace);
    }

  if (Space->GetHilbertSpaceDimension() != VectorDimension)
    {
      cout << "Dimension of vectors does not match the size of the Hilbert-Space!"<<endl;
      exit(1);	
    }

  if (Manager.GetBoolean("first-real"))
    {      
      for (int k=0; k<NbrVectors; ++k)
	{
	  int e=(Manager.GetInteger("real-element")<Space->GetHilbertSpaceDimension()?Manager.GetInteger("real-element"):0);
	  while ((Norm(Vectors[k][e])<1e-10)&&(e<Space->GetHilbertSpaceDimension()))
	    ++e;
	  Phases[k]=Polar(1.0,-Arg(Vectors[k][e]));
	}
	  
    }

  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
    {
      Space->PrintState(cout, i);
      if (NbrVectors>0)
	cout << " : " << "("<<Real(Vectors[0][i]*Phases[0])<<"+I*"<<Imag(Vectors[0][i]*Phases[0])<<")";
      for (int k=1; k<NbrVectors; ++k)
	cout << "  " << "("<<Real(Vectors[k][i]*Phases[k])<<"+I*"<<Imag(Vectors[k][i]*Phases[k])<<")";
      cout << endl;
    }

  
  delete Space;  
  return 0;
}
