#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Operator/ParticleOnLatticeDensityDensityOperator.h"
#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

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


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");

  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);  
  (*MiscGroup) += new BooleanOption  ('g', "gnuplot", "format output for gnuplot, and make plot");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  bool DensityFlag = ((BooleanOption*) Manager["density"])->GetBoolean();

  if (NbrVectors==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }
  ParticleOnLattice* Space = 0;
  bool Plot = Manager.GetBoolean("gnuplot");
  for (int i=0; i<NbrVectors; ++i)
    {      
      double Interaction=0.0;
      int TmpI=-1;
      bool Statistics=false;
      bool HardCore=false;
      if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[i], NbrBosons, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	{
	  cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
	  exit(-1);
	}
      HardCore=(HardCore||Manager.GetBoolean("hard-core"));
  
      int VectorDimension=0;
      ComplexVector State;
      if (State.ReadVector (VectorFiles[i]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[i] << endl;
	  return -1;      
	}
      VectorDimension=State.GetVectorDimension();
      char* OutputNameCorr = new char [10 + strlen (VectorFiles[i])];
      char* ending;
      if ((ending=strstr(VectorFiles[i],".vec"))!=0)
	strncpy(OutputNameCorr, VectorFiles[i], ending-VectorFiles[i]);
      else
	strcpy(OutputNameCorr, VectorFiles[i]);
      char* PlotNameCmd=0;
      char* PlotNamePS=0;
      if (Plot)
	{
	  PlotNameCmd = new char  [10 + strlen (VectorFiles[i])];
	  PlotNamePS = new char  [10 + strlen (VectorFiles[i])];
	  sprintf (PlotNameCmd, "%s.gp", OutputNameCorr);
	  sprintf (PlotNamePS, "%s.ps", OutputNameCorr);
	}
      if (DensityFlag == false)
	sprintf (OutputNameCorr, "%s.rho_rho.dat", OutputNameCorr);
      else
	sprintf (OutputNameCorr, "%s.rho.dat", OutputNameCorr);
      
	
      cout << "<in  "<<VectorFiles[i]<<endl<<">out "<<OutputNameCorr<<endl;
      if ((Space == 0) || (Space->GetHilbertSpaceDimension()!=VectorDimension))
	{
	  if (Space == 0) delete Space;
	  if (HardCore)
	    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
	  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
	  
	  if (VectorDimension != Space->GetHilbertSpaceDimension())
	    {
	      cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of Hilbert-space!"<<endl;
	      exit(-11);
	    }
	}
      ofstream File;
      File.precision(14);
      File.open(OutputNameCorr, ios::binary | ios::out);      
      Complex Tmp;
      if (DensityFlag == false)
	{
	  File << "# density-density correlation for " << VectorFiles[i]<< endl;
	  File << "# x\ty\tg"<< endl;
	  ParticleOnLatticeOneBodyOperator Operator0 (Space, 0, 0);
	  Complex Zero=Operator0.MatrixElement(State, State);
	  for (int x = 0; x < Lx; ++x)
	    {
	      for (int y = 0; y < Ly; ++y)
		{
		  int q=Space->EncodeQuantumNumber(x, y, 0, Tmp);
		  ParticleOnLatticeDensityDensityOperator Operator (Space, q, 0, q, 0);
		  double ME=Real(Operator.MatrixElement(State, State));
		  if (q==0) ME = Real(Zero)*Real(Zero);
		  File << x << "\t" << y << "\t" << ME << endl;		  
		}
	      if (Plot) File << endl;
	    }
	}
      else
	{
	  File << "# density-profile for " << VectorFiles[i]<< endl;
	  File << "# x\t<\tg"<< endl;
	  for (int x = 0; x < Lx; ++x)
	    {
	      for (int y = 0; y < Ly; ++y)
		{
		  int q=Space->EncodeQuantumNumber(x, y, 0, Tmp);
		  ParticleOnLatticeOneBodyOperator Operator (Space, q, q);
		  File << x << "\t" << y << "\t" << Real(Operator.MatrixElement(State, State)) << endl;		  
		}
	      if (Plot) File << endl;
	    }
	}
      File.close();
      if (Plot)
	{
	  File.open(PlotNameCmd, ios::binary | ios::out);      
	  File<< "set pm3d at b" << endl <<"set xlabel \"x\""<<endl<<"set ylabel \"y\""<<endl;
	  File<<"set terminal postscript eps enhanced \"Helvetica\" 14 color"<<endl;
	  File<<"set palette rgbformulae 4,9,15"<<endl;
	  File<<"set output \""<<PlotNamePS<<"\""<<endl;
	  File<<"splot \""<<OutputNameCorr<<"\" u 1:2:3 title \"correlations\" w lines 1"<<endl;
	  File.close();
	  char tmpC[255];
	  sprintf(tmpC,"gnuplot %s",PlotNameCmd);
	  system(tmpC);
	  cout<< "created graph "<<PlotNamePS<<endl;
	  delete [] PlotNameCmd;
	  delete [] PlotNamePS;
	}
      delete[] OutputNameCorr;	  

    }

  return 0;
}


