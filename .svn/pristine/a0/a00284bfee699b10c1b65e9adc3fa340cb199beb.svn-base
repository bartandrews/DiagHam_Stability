#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"
#include "HilbertSpace/BosonOnLatticeGeneric.h"

#include "Tools/FQHESpectrum/LatticePhases.h"

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
#include <cstring>
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
  OptionManager Manager ("FQHELatticeBosonsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('r', "reference", "reference site for two-particle correlations", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");

  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "fluctuations", "plot local density fluctuations", false);
  (*SystemGroup) += new BooleanOption  ('\n', "all-ops", "plot all available operators", false);
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
  bool CorrelationFlag = Manager.GetBoolean("density");
  bool DensityFlag = Manager.GetBoolean("density");
  bool FluctuationFlag = Manager.GetBoolean("fluctuations");
  bool AllOps = Manager.GetBoolean("all-ops");

  if (DensityFlag||FluctuationFlag)
    CorrelationFlag=false;
  
  if (NbrVectors==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }

  double Interaction=-1.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  bool GenericLattice=false;
  int NbrSites=0;
  int NbrSubLattices=1;
  
  LatticePhases *Lattice = NULL;

  ParticleOnLattice* Space = 0;
  bool Plot = Manager.GetBoolean("gnuplot");
  for (int i=0; i<NbrVectors; ++i)
    {
      if (FQHEOnLatticeHaveGeneralLattice(VectorFiles[i]))
	{
	  GenericLattice=true;
	  if (Manager.GetString("lattice-definition")==NULL)
	    {
	      cout << "Please indicate the file with the lattice-definition for this vector"<<endl;
	      exit(1);
	    }
	  // get the lattice geometry
	  Lattice = new LatticePhases();
	  NbrSites = Lattice->GetNbrSites();
	  Lx = Lattice->GetLatticeLength(0);
	  Ly = Lattice->GetLatticeLength(1);
	  NbrSubLattices = Lattice->GetNbrSubLattices();
	  char* LatticeName = Lattice->GeometryString();
	  if (strstr(VectorFiles[0], LatticeName)==0)
	    {
	      cout << "The given lattice parameters do not coincide with the filename, verify lattice definition, and repetition of unit cells"<<endl;
	    }
	  bool HaveContFlux;
	  double ContFlux;
	  delete [] LatticeName;
	  if (FQHEOnLatticeFindSystemInfoFromGeneralVectorFileName(VectorFiles[0], NbrBosons, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore, HaveContFlux, ContFlux) == false)
	    {
	      cout<<"Please use standard file-names, or indicate all necessary system parameters!"<<endl;
	      exit(1);
	    }
	  
	}
      else
	{
	  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[0], NbrBosons, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	    {
	      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
	      exit(1);
	    }
	  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
	  if (Manager.GetBoolean("no-hard-core"))
	    HardCore=false;
	  NbrSites = Lx*Ly;
	}

      int ReferenceSite = Manager.GetInteger("reference")%NbrSites;
      
      
      int VectorDimension=0;
      ComplexVector State;
      if (State.ReadVector (VectorFiles[i]) == false)
	{
	  cout << "can't open vector file " << VectorFiles[i] << endl;
	  return -1;      
	}
      VectorDimension=State.GetVectorDimension();
      char* OutputNameCorr = new char [64 + strlen (VectorFiles[i])];
      char* OutputNameBase = new char [64 + strlen (VectorFiles[i])];
      char* ending;
      if ((ending=strstr(VectorFiles[i],".vec"))!=0)
	{
	  char *last;
	  while ((last=strstr(ending+3,".vec"))!=0)
	    {
	      ending = last;
	    }
	  strncpy(OutputNameBase, VectorFiles[i], ending-VectorFiles[i]);
	  OutputNameBase[ending-VectorFiles[i]]='\0';
	}
      else
	strcpy(OutputNameBase, VectorFiles[i]);
      if (ReferenceSite!=0)
	sprintf (OutputNameBase, "%s_ref_%d", OutputNameBase, ReferenceSite);
      char* PlotNameCmd=0;
      char* PlotNamePS=0;
      if (Plot)
	{
	  PlotNameCmd = new char  [10 + strlen (VectorFiles[i])];
	  PlotNamePS = new char  [10 + strlen (VectorFiles[i])];
	  sprintf (PlotNameCmd, "%s.gp", OutputNameBase);
	  sprintf (PlotNamePS, "%s.ps", OutputNameBase);
	}
      if (AllOps == false)
	{
	  if (FluctuationFlag == false)
	    {
	      if (DensityFlag == false)
		sprintf (OutputNameCorr, "%s.rho_rho.dat", OutputNameBase);
	      else
		sprintf (OutputNameCorr, "%s.rho.dat", OutputNameBase);
	    }
	  else sprintf (OutputNameCorr, "%s.sigma_rho.dat", OutputNameBase);
	}
      else
	sprintf (OutputNameCorr, "%s_CHANNEL.dat", OutputNameBase);
      
	
      cout << "<in  "<<VectorFiles[i]<<endl<<">out "<<OutputNameCorr<<endl;
      if ((Space == 0) || (Space->GetHilbertSpaceDimension()!=VectorDimension))
	{
	  if (Space != 0) delete Space;
	  if (GenericLattice)
	    {
	      if (HardCore)
		Space = new HardCoreBosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace);
	      else Space = new BosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace);
	    }
	  else
	    {
	      if (HardCore)
		Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
	      else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
	    }
	  
	  if (VectorDimension != Space->GetHilbertSpaceDimension())
	    {
	      cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of Hilbert-space!"<<endl;
	      exit(-11);
	    }
	}
      ofstream File;
      File.precision(14);
      Complex Tmp;
      int CellPosition[2];
      RealVector SitePosition;
      if ((CorrelationFlag == true)||(AllOps == true))
	{
	  sprintf (OutputNameCorr, "%s.rho_rho.dat", OutputNameBase);
	  File.open(OutputNameCorr, ios::binary | ios::out);
	  File << "# density-density correlation for " << VectorFiles[i]<< endl;
	  File << "# x\ty\tg"<< endl;
	  ParticleOnLatticeOneBodyOperator Operator0 (Space, ReferenceSite, ReferenceSite);
	  Complex Zero=Operator0.MatrixElement(State, State);
	  for (int x = 0; x < Lx; ++x)
	    {
	      for (int y = 0; y < Ly; ++y)
		{
		  for (int Sub=0; Sub<NbrSubLattices; ++Sub)
		    {
		      int q=Space->EncodeQuantumNumber(x, y, Sub, Tmp);
		      ParticleOnLatticeDensityDensityOperator Operator (Space, q, ReferenceSite);
		      double ME;
		      if (q==ReferenceSite)
			ME = Real(Zero)*(Real(Zero)-1);
		      else
			{
			  ParticleOnLatticeOneBodyOperator Operator1 (Space, q, q);
			  Complex One=Operator1.MatrixElement(State, State);
			  ME=Real(Operator.MatrixElement(State, State))-Real(Zero)*Real(One);
			}
		      double X,Y;
		      if (GenericLattice)
			{
			  CellPosition[0]=x;
			  CellPosition[1]=y;
			  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
			  X=SitePosition[0];
			  Y=SitePosition[1];
			}
		      else
			{
			  X=(double)x;
			  Y=(double)y;
			}
		      
		      File << X << "\t" << Y << "\t" << ME << endl;
		    }
		}
	      if (Plot) File << endl;
	    }
	  File.close();
	}
      if ((DensityFlag == true)||(AllOps==true))
	{
	  sprintf (OutputNameCorr, "%s.rho.dat", OutputNameBase);
	  File.open(OutputNameCorr, ios::binary | ios::out);
	  File << "# density-profile for " << VectorFiles[i]<< endl;
	  File << "# x\ty\tg"<< endl;
	  for (int x = 0; x < Lx; ++x)
	    {
	      for (int y = 0; y < Ly; ++y)
		for (int Sub=0; Sub<NbrSubLattices; ++Sub)
		  {
		    int q=Space->EncodeQuantumNumber(x, y, Sub, Tmp);
		    ParticleOnLatticeOneBodyOperator Operator (Space, q, q);
		    double X,Y;
		    if (GenericLattice)
		      {
			CellPosition[0]=x;
			CellPosition[1]=y;
			SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
			X=SitePosition[0];
			Y=SitePosition[1];
		      }
		    else
		      {
			X=(double)x;
			Y=(double)y;
		      }
		    File << X << "\t" << Y << "\t" << Real(Operator.MatrixElement(State, State)) << endl;		  
		  }
	      if (Plot) File << endl;
	    }
	  File.close();
	}
      if ((FluctuationFlag == true)||(AllOps==true))
	{
	  sprintf (OutputNameCorr, "%s.sigma_rho.dat", OutputNameBase);
	  File.open(OutputNameCorr, ios::binary | ios::out);
	  File << "# local density fluctuations for " << VectorFiles[i]<< endl;
	  File << "# x\ty\tg"<< endl;
	  for (int x = 0; x < Lx; ++x)
	    {
	      for (int y = 0; y < Ly; ++y)
		{
		  for (int Sub=0; Sub<NbrSubLattices; ++Sub)
		    {
		      int q=Space->EncodeQuantumNumber(x, y, Sub, Tmp);
		      ParticleOnLatticeOneBodyOperator OperatorN (Space, q, q);
		      ParticleOnLatticeDensityDensityOperator OperatorNSqr (Space, q, q);
		      double N=Real(OperatorN.MatrixElement(State, State));
		      double NSqr=Real(OperatorNSqr.MatrixElement(State, State));
		      double X,Y;
		      if (GenericLattice)
			{
			  CellPosition[0]=x;
			  CellPosition[1]=y;
			  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
			  X=SitePosition[0];
			  Y=SitePosition[1];
			}
		      else
			{
			  X=(double)x;
			  Y=(double)y;
			}
		      
		      File << X << "\t" << Y << "\t" << NSqr+N*(1-N) << endl;
		    }
		}
	      if (Plot) File << endl;
	    }
	  File.close();
	}
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


