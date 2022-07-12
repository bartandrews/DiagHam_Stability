#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


void SymmetrizeVectors(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
		       ComplexVector *Vectors, int k, int nbrK, ArchitectureManager &Architecture)
{
  Complex Tmp;
  for (int i=0; i<NbrVectors; ++i)
    {
      ComplexVector TmpState(Vectors[i].GetVectorDimension(), true);
      ComplexVector TranslatedState(Vectors[i], true);
      ComplexVector ResultingState(Vectors[i], true);
      for (int n=1; n<nbrK; ++n)
	{
	  VectorOperatorMultiplyOperation Operation (Operator, &TranslatedState, &TmpState);      
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  ResultingState.AddLinearCombination (Polar(1.0,2.0*M_PI*(double)k/(double)nbrK*n), TmpState);
	  TranslatedState.Copy(TmpState);
	}
      if (ResultingState.SqrNorm()<1e-30)
	{
	  cout<<"Attention, vector "<<i<<" has no component in this momentum"<<endl;	  
	}
      else
	{
	  ResultingState/=ResultingState.Norm();
	  Vectors[i].Copy(ResultingState);
	}
    }
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeProjectMomentum" , "0.01");  
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");

  (*SystemGroup) += new SingleIntegerOption  ('k', "target-kx", "target momentum in x-direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "target-ky", "target momentum in y-direction", 0);
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*OutputGroup) += new SingleStringOption ('o', "output", "write resulting vector to this file",NULL);
  (*OutputGroup) += new BooleanOption  ('V', "verbose", "give additional output");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);
  
  if (NbrVectors==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }
  double Interaction=-1.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[0], NbrBosons, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }  
  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;
  
  int VectorDimension=0;
  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  bool tmpB, haveVector=false;
  for (int i=0; i<NbrVectors; ++i)
    {
      tmpB = Vectors[i].ReadVector(VectorFiles[i]);
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
      cout << "No valid vector files found!"<<endl;
      exit(1);
    }  

  ParticleOnLattice* Space;
  if (HardCore)
    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);

  if (VectorDimension != Space->GetHilbertSpaceDimension())
    {
      cout<<"Dimension of vectors does not match size of Hilbert-space!"<<endl;
	  exit(1);
    }
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  ParticleOnLatticeTranslationOperator *TranslationOperator= new ParticleOnLatticeTranslationOperator(Space);

  cout<< "====== Constructing momentum eigenvector ====="<<endl;

  int Degeneracy=1;
  int n1=1, n2=1;
  int FluxModulo = FindGCD(NbrFluxQuanta, Lx*Ly);
  int r=NbrFluxQuanta/FluxModulo;
  int t=Lx*Ly/FluxModulo;

  while ((((Ly*n1)%t)!=0) && (n1<Lx)) ++n1;
  while ((((Lx*n2)%t)!=0) && (n2<Ly)) ++n2;

  if ((Lx%n1)!=0)
    cout << "Extending range of n1 to Lx"<<endl;
  if ((Ly%n2)!=0)
    cout << "Extending range of n2 to Ly"<<endl;

  if (((n1*n2*NbrFluxQuanta)%t) != 0)
    {
      cout << "Cannot resolve translations: Brillouin zone trivial?"<<endl;
      n1=Lx;
      n2=Ly;
    }

  while ((r*NbrBosons*n1*n2*Degeneracy)%t != 0) ++Degeneracy;
  
  cout << "N_phi = "<<r<<"/"<<t<<endl;
  cout << "n1="<<n1<<", n2="<<n2<<", global degeneracy: "<<Degeneracy<<endl;

  int RemainingDegeneracy=Degeneracy;

  if ((Ly/n2)<Degeneracy)
    {
      int GCD = FindGCD(Lx/n1, Degeneracy);
      RemainingDegeneracy/=GCD;
      if (GCD!=1) cout << "Multiplying factor "<<GCD<<" of degeneracy onto n1"<<endl;
      n1*=GCD;      
    }

  if ((Ly/n2)%RemainingDegeneracy!=0)
    {
      cout<<"Did not treat degeneracy properly -> need to put onto n1?"<<endl;
      exit(1);      
    }
  else
    {
      if (RemainingDegeneracy!=1)
	cout << "Multiplying factor "<<RemainingDegeneracy<<" of degeneracy onto n2"<<endl;
      n2*=RemainingDegeneracy;
      RemainingDegeneracy=1;
    }
  
  int kx=Manager.GetInteger("target-kx");
  int ky=Manager.GetInteger("target-ky");
  
  TranslationOperator->SetTranslationComponents(n1,0);
  SymmetrizeVectors(TranslationOperator, NbrVectors, Vectors, kx, Lx/n1, Architecture);
  
  TranslationOperator->SetTranslationComponents(0,n2);
  SymmetrizeVectors(TranslationOperator, NbrVectors, Vectors, ky, Ly/n2, Architecture);
  

  if (NbrVectors>1)
    {
      char *NewExtension=new char[20];
      sprintf(NewExtension,"kx_%d_ky_%d.vec",kx,ky);
      char *OldExtension=GetExtensionFromFileName(VectorFiles[0]);
      for (int i=0; i<NbrVectors; ++i)
	{
	  char *OutputName=0;
	  if ((OldExtension==0) || (strcmp(OldExtension,".vec")!=0))
	    OutputName = AddExtensionToFileName(VectorFiles[i], NewExtension);
	  else
	    {
	      OutputName = ReplaceExtensionToFileName(VectorFiles[i], "vec", NewExtension);
	      if (OutputName==NULL)
		OutputName = AddExtensionToFileName(VectorFiles[i], NewExtension);
	    }
	  Vectors[i].WriteVector(OutputName);
	  delete [] OutputName;
	}
      delete [] NewExtension;
      delete [] OldExtension;
    }
  else
    {
      char *OutputName=Manager.GetString("output");
      if (OutputName==NULL)
	{
	  char *NewExtension=new char[20];
	  sprintf(NewExtension,"kx_%d_ky_%d.vec",kx,ky);
	  char *OldExtension=GetExtensionFromFileName(VectorFiles[0]);
	  if ((OldExtension==0) || (strcmp(OldExtension,".vec")!=0))
	    OutputName = AddExtensionToFileName(VectorFiles[0], NewExtension);
	  else
	    {
	      OutputName = ReplaceExtensionToFileName(VectorFiles[0], "vec", NewExtension);
	      if (OutputName==NULL)
		OutputName = AddExtensionToFileName(VectorFiles[0], NewExtension);
	    }
	  delete [] NewExtension;
	  delete [] OldExtension;
	}
      Vectors[0].WriteVector(OutputName);
      delete [] OutputName;
    }
  
  delete Space;  
  delete [] Vectors;
  delete TranslationOperator;
}
