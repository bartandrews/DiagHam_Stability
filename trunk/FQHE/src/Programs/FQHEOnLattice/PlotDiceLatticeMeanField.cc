#include "Vector/ComplexVector.h"

#include "Tools/FQHESpectrum/LatticePhases.h"
#include "HilbertSpace/SingleParticleOnLatticeGeneric.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

int GetSiteIndex(int SiteIndex)
{
  if ( SiteIndex % 2 == 0)
    return (6*(SiteIndex/2)+1);
  else
    return (6*((SiteIndex-1)/2)+5);
}



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("PlotDiceLatticeMeanField" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  Manager += PrecalculationGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "wavefunction", "parameters of a mean-field GP wavefunction");

  (*SystemGroup) += new MultipleIntegerOption  ('\n',"hub-offset", "offset of index for 6x connected site",',',0,"1,5");

  (*MiscGroup) += new BooleanOption  ('g', "gnuplot", "format output for gnuplot, and make plot");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
    
  char* WaveFunctionFile = Manager.GetString("wavefunction");


  if (WaveFunctionFile==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }

  int NbrSites=0;
  int NbrSubLattices=6;
  
  LatticePhases *Lattice = NULL;

  bool Plot = Manager.GetBoolean("gnuplot");

  if (Manager.GetString("lattice-definition")==NULL)
    {
      cout << "Please indicate the file with the lattice-definition for the correponding dice lattice"<<endl;
      exit(1);
    }
  // get the lattice geometry
  Lattice = new LatticePhases();
  NbrSites = Lattice->GetNbrSites();
  int Lx = Lattice->GetLatticeLength(0);
  int Ly = Lattice->GetLatticeLength(1);
  NbrSubLattices = Lattice->GetNbrSubLattices();
  char* LatticeGeometryStr = new char[10];
  sprintf(LatticeGeometryStr,"_%dx%d_",Lx,Ly);
  if (strstr(WaveFunctionFile, LatticeGeometryStr)==0)
    {
      cout << "The given lattice parameters do not coincide with the filename, verify lattice definition, and repetition of unit cells"<<endl;
    }
  delete [] LatticeGeometryStr;
  
      
  int WaveFunctionDimension=0;
  ComplexVector State;
  if (State.ReadVector (WaveFunctionFile) == false)
    {
      cout << "can't open vector file " << WaveFunctionFile << endl;
      return -1;      
    }
  WaveFunctionDimension=State.GetVectorDimension();
  char* OutputNamePlot = new char [64 + strlen (WaveFunctionFile)];
  char* ending;
  if ((ending=strstr(WaveFunctionFile,".par"))!=0)
    {
      char *last;
      while ((last=strstr(ending+3,".par"))!=0)
	{
	  ending = last;
	}
      strncpy(OutputNamePlot, WaveFunctionFile, ending-WaveFunctionFile);
      OutputNamePlot[ending-WaveFunctionFile]='\0';
    }
  else
    strcpy(OutputNamePlot, WaveFunctionFile);

  char* PlotNameCmd=0;
  char* PlotNamePS=0;
  if (Plot)
    {
      PlotNameCmd = new char  [10 + strlen (WaveFunctionFile)];
      PlotNamePS = new char  [10 + strlen (WaveFunctionFile)];
      sprintf (PlotNameCmd, "%s.gp", OutputNamePlot);
      sprintf (PlotNamePS, "%s.ps", OutputNamePlot);
    }
  sprintf (OutputNamePlot, "%s.psi", OutputNamePlot);
      
	
  cout << "<in  "<<WaveFunctionFile<<endl<<">out "<<OutputNamePlot<<endl;
	  
  if (WaveFunctionDimension != 2*Lx*Ly)
    {
      cout<<"Dimension of wavefunction "<<WaveFunctionFile<<" does not match expected size!"<<endl;
      exit(-11);
    }

  int NbrSite, NbrSiteEffective;
  int NbrNeighbors;
  int* Neighbors;
  double *Phases, *Amplitudes;
  int** PeriodicTranslations;
  int Length;
  int *HubOffset;
  HubOffset = Manager.GetIntegers("hub-offset",Length);
  if (Length!=2)
    {
      cout << "Please indicate the sublattice indices of the two hub sites in each unit cell"<<endl;
      exit(1);
    }
  

  SingleParticleOnLatticeGeneric *Space = new SingleParticleOnLatticeGeneric(Lattice, 0, 0.0, 0.0);

  ComplexVector ActualWaveFunction;
  ActualWaveFunction.ResizeAndClean(6*Lx*Ly);
  ComplexVector LocalWaveFunction;
  LocalWaveFunction.ResizeAndClean(6*Lx*Ly);  
  Complex Tmp;
  
  for (int x = 0; x < Lx; ++x)
    {
      for (int y = 0; y < Ly; ++y)
	{
	  for (int Sub=0; Sub<2; ++Sub)
	    {
	      NbrSiteEffective = 2*(x+y*Lx)+Sub;
	      NbrSite=Space->EncodeQuantumNumber(x, y, HubOffset[Sub], Tmp);
	      //cout << "Same? "<<NbrSiteEffective<<", "<< (NbrSite-HubOffset[Sub])/3+Sub <<endl;
	      cout << "Same? "<<NbrSite<<", "<< GetSiteIndex(NbrSiteEffective) << endl;
	      Lattice->GetNeighbors(NbrSite, NbrNeighbors, Neighbors, Phases, PeriodicTranslations, Amplitudes);
	      if (NbrNeighbors!=6)
		{
		  cout << "Did not find 6x connected hub on site "<< NbrSite <<" actual neighbors: "<< NbrNeighbors<<endl;
		  exit(1);
		}
	      LocalWaveFunction.ClearVector();
	      LocalWaveFunction[NbrSites-NbrSite-1].Re=1.0/sqrt(2.0);  
	      for (int n=0; n<NbrNeighbors; ++n)
		{
		  LocalWaveFunction[NbrSites-Neighbors[n]-1].Re=std::cos(M_PI*Phases[n])/sqrt(12.0);
		  LocalWaveFunction[NbrSites-Neighbors[n]-1].Im=std::sin(M_PI*Phases[n])/sqrt(12.0);
		}
	      // cout << "Local state on effective site "<<NbrSiteEffective<<endl<<LocalWaveFunction;
	      ActualWaveFunction.AddLinearCombination(State[NbrSiteEffective], LocalWaveFunction);
	    }
	}
    }  
  
  // now, go on to plot things
  
  ofstream File;
  File.precision(14);
  File.open(OutputNamePlot, ios::binary | ios::out);

  int CellPosition[2];
  int CellPosition2[2];
  RealVector SitePosition, SitePosition2;
  double X,Y;

  for (int x = 0; x < Lx; ++x)
    for (int y = 0; y < Ly; ++y)
      for (int Sub=0; Sub<6; ++Sub)
	{
	  NbrSite=Space->EncodeQuantumNumber(x, y, Sub, Tmp);
	  CellPosition[0]=x;
	  CellPosition[1]=y;
	  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	  X=SitePosition[0];
	  Y=SitePosition[1];
	  
	  File << X << "\t" << Y << "\t" << ActualWaveFunction[NbrSites-NbrSite-1].Re << "\t" << ActualWaveFunction[NbrSites-NbrSite-1].Im << "\t" << Norm(ActualWaveFunction[NbrSites-NbrSite-1]) << endl;
	}

  File.close();
  
  char * OutputNameTheta = ReplaceExtensionToFileName(OutputNamePlot, "psi", "theta");

  File.open(OutputNameTheta, ios::binary | ios::out);
  int Sub2;
  double CentralPhase, OtherPhase, Theta;
  for (int x = 0; x < Lx; ++x)
    {
      for (int y = 0; y < Ly; ++y)
	{
	  for (int Sub=0; Sub<2; ++Sub)
	    {
	      NbrSiteEffective = 2*(x+y*Lx)+Sub;
	      NbrSite=Space->EncodeQuantumNumber(x, y, HubOffset[Sub], Tmp);
	      Lattice->GetNeighbors(NbrSite, NbrNeighbors, Neighbors, Phases, PeriodicTranslations, Amplitudes);
	      if (NbrNeighbors!=6)
		{
		  cout << "Did not find 6x connected hub"<<endl;
		  exit(1);
		}
	      CellPosition[0]=x;
	      CellPosition[1]=y;
	      SitePosition = Lattice->GetSitePosition(CellPosition,HubOffset[Sub]);
	      CentralPhase = Arg(ActualWaveFunction[NbrSites-NbrSite-1]);

	      cout << "central position:"<<endl<< SitePosition;
	      for (int n=0; n<NbrNeighbors; ++n)
		{
		  Lattice->GetSiteCoordinates(Neighbors[n],CellPosition2, Sub2);
		  SitePosition2 = Lattice->GetSitePosition(CellPosition2,Sub2);
		  SitePosition2.AddLinearCombination(-(double)PeriodicTranslations[n][0]*Lx,Lattice->GetLatticeVector(0));
		  SitePosition2.AddLinearCombination(-(double)PeriodicTranslations[n][1]*Ly,Lattice->GetLatticeVector(1));
		  cout << "connection "<<NbrSite<<"-"<<Neighbors[n]<<endl<<"with other position "<<endl<<SitePosition2;
		  cout << "from translations "<<PeriodicTranslations[n][0]<<","<<PeriodicTranslations[n][1]<<endl;
		  SitePosition2 += SitePosition;
		  SitePosition2*=0.5;
		  cout << "average"<<endl<<SitePosition2;
		  OtherPhase = Arg(ActualWaveFunction[NbrSites-Neighbors[n]-1]);
		  Theta = CentralPhase - OtherPhase + M_PI*Phases[n]; // check conventions against FQHELatticeBosonsMeanField
		  while( Theta < -M_PI)
		    Theta+=2.0*M_PI;
		  while( Theta > M_PI)
		    Theta-=2.0*M_PI;
		  File << SitePosition2[0] << "\t" << SitePosition2[1] << "\t" << (Theta>0.0?Theta:0.0) << "\t";
		  File << (Theta<0.0?-Theta:0.0) << "\t" << -2.0*Imag(Conj(ActualWaveFunction[NbrSites-NbrSite-1])*ActualWaveFunction[NbrSites-Neighbors[n]-1]*Polar(1.0, -M_PI*Phases[n])) << endl;
		}
	    }
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
      File<<"splot \""<<OutputNamePlot<<"\" u 1:2:3 title \"correlations\" w lines 1"<<endl;
      File.close();
      char tmpC[255];
      sprintf(tmpC,"gnuplot %s",PlotNameCmd);
      system(tmpC);
      cout<< "created graph "<<PlotNamePS<<endl;
      delete [] PlotNameCmd;
      delete [] PlotNamePS;
    }
  delete[] OutputNamePlot;
  delete[] OutputNameTheta;
  
  delete Space;
  return 0;
}


