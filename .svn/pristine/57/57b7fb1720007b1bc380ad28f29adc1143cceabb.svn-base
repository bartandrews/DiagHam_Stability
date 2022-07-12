#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEMonteCarlo/AbstractMCSamplingFunction.h"
#include "Tools/FQHEMonteCarlo/QHESamplingFunctionManager.h"
#include "Tools/FQHEMonteCarlo/SimpleMonteCarloAlgorithm.h"
//#include "Tools/FQHEMonteCarlo/DiskCoulombEnergy.h"
//#include "Tools/FQHEMonteCarlo/DiskGeneralEnergy.h"
#include "Tools/FQHEMonteCarlo/SimpleTwoBodyCorrelatorOnDisk.h"
#include "Tools/FQHEMonteCarlo/SimpleDensityOnDisk.h"

#include "MCObservables/RealObservable.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskMonteCarlo" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
    
  QHEWaveFunctionManager WaveFunctionManager(QHEWaveFunctionManager::DiskGeometry, QHEWaveFunctionManager::DiskWithBackground);
  QHESamplingFunctionManager SamplingFunctionManager(QHESamplingFunctionManager::DiskGeometry);
  
  Manager += SystemGroup;
  SimpleMonteCarloAlgorithm::AddOptionGroup(&Manager);
  WaveFunctionManager.AddOptionGroup(&Manager);
  SamplingFunctionManager.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum angular momentum per particle", 15);
  (*SystemGroup) += new SingleStringOption('\n',"interaction-params","File containing parameters of interaction");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");

  (*SystemGroup) += new BooleanOption ('b', "background-only", "print background energy, and exit");
  
  (*SystemGroup) += new SingleStringOption('o',"output","file to write a log of results");  
  (*SystemGroup) += new SingleStringOption('\n',"plot","plot effective interaction to file and exit",NULL);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-points", "number points in plot", 250);

  (*SystemGroup) += new BooleanOption  ('c', "correlations", "sample correlation function");

  (*SystemGroup) += new SingleIntegerOption ('\n', "corr-resolution", "number of bins for correlation function", 512);
  (*SystemGroup) += new SingleIntegerOption ('\n', "corr-highres-range", "number of bins to magnify at small r", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "corr-highres-factor", "ratio by which to inrease density of bins in high res area", 3);

  (*SystemGroup) += new BooleanOption  ('d', "density", "sample correlation function");

  (*SystemGroup) += new SingleIntegerOption ('\n', "rho-resolution", "number of bins for density observable", 512);
  (*SystemGroup) += new SingleIntegerOption ('\n', "rho-highres-range", "number of bins to magnify at small r", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "rho-highres-factor", "ratio by which to inrease density of bins in high res area", 3);


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout); // new feature to process options and display help in 1 line!

  if (Manager.GetBoolean("list-wavefunctions") == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }
  
  SamplingFunctionManager.TestForShow();
  

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrFlux = Manager.GetInteger("lzmax");
  
  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  AbstractMCSamplingFunction* SamplingFunction = SamplingFunctionManager.GetSamplingFunction();

  if (TestWaveFunction==0 && SamplingFunction==0)
    {
      cout << "Invalid wavefunction requested - list available states with --list-wavefunctions"<<endl;
      exit(-1);
    }
  if (TestWaveFunction!=0)
    cout << "Function: " << WaveFunctionManager.GetDescription()<<endl;
  if (SamplingFunction!=0)
    cout << "Sampler:  " << SamplingFunctionManager.GetDescription()<<endl;

  SimpleMonteCarloAlgorithm MonteCarloRoutine(AbstractParticleCollection::OnDiskCollection, NbrParticles, TestWaveFunction, SamplingFunction,
						      &Manager);
  
/*  
  if (Manager.GetString("interaction-params")==0)
    {
      // add Coulomb energy in LLL
      DiskCoulombEnergy *Energy=new DiskCoulombEnergy(NbrFlux);
      MonteCarloRoutine.AddObservable(Energy);
    }
  else
    {
      // add a general energy function 
      DiskGeneralEnergy *Energy=new DiskGeneralEnergy(NbrFlux, Manager.GetString("interaction-params"));
      if (Manager.GetString("plot")!=NULL)
	{
	  ofstream MyFile(Manager.GetString("plot"),std::ios::out);
	  Energy->PlotPotential(MyFile,Manager.GetInteger("nbr-points"));
	  MyFile.close();
	  cout << "Interaction plotted in file "<<Manager.GetString("plot")<<endl;
	  exit(0);
	}
      
      MonteCarloRoutine.AddObservable(Energy);
      if (Manager.GetBoolean("background-only"))
	{
	  cout << "# Background energy"<<endl;
	  cout << "bg " << Energy->GetTotalBackgroundEnergy()<<endl;
	  if (Manager.GetString("output")!=NULL)
	    {
	      ofstream LogFile(Manager.GetString("output"),std::ios::out);
	      LogFile.precision(12);
	      LogFile << "# Background energy"<<endl;
	      LogFile << Energy->GetTotalBackgroundEnergy()<<endl;
	      LogFile.close();
	    }
	  exit(0);
	}
    }
  */
  
  
  SimpleTwoBodyCorrelatorOnDisk *Correlations=NULL;
  if (Manager.GetBoolean("correlations"))
    {
      int highResRange=Manager.GetInteger("corr-highres-range");
      int highResFactor=Manager.GetInteger("corr-highres-factor");
      std::cout << "NbrFlux="<<NbrFlux<<endl;
      Correlations = new SimpleTwoBodyCorrelatorOnDisk(2.5*std::sqrt((double)NbrFlux), Manager.GetInteger("corr-resolution"), highResFactor * highResRange, highResRange);
      MonteCarloRoutine.AddObservable(Correlations, 20);
      Correlations->IncludeInPrint(false); // do not write output to logfile.
    }

  SimpleDensityOnDisk *Density=NULL;
  if (Manager.GetBoolean("density"))
    {
      int highResRange=Manager.GetInteger("rho-highres-range");
      int highResFactor=Manager.GetInteger("rho-highres-factor");
      std::cout << "NbrFlux="<<NbrFlux<<endl;
      Density = new SimpleDensityOnDisk((1.25*std::sqrt((double)NbrFlux)+2)/(1.0-Manager.GetDouble("defect-angle")), Manager.GetInteger("rho-resolution"), highResFactor * highResRange, highResRange);
      MonteCarloRoutine.AddObservable(Density, 20);
      Density->IncludeInPrint(false); // do not write output to logfile.
    }

  // run simulation
  MonteCarloRoutine.Simulate();

  if (Manager.GetString("output")!=NULL)
    {
      ofstream LogFile(Manager.GetString("output"),std::ios::out);
      MonteCarloRoutine.WriteObservations(LogFile);
      LogFile.close();
    }

  // print final results:
  cout << "Final results:" << endl;
  MonteCarloRoutine.WriteObservations(cout);  
  
  
  if(Correlations!=NULL)
    {
      ofstream LogFile("correlations.dat",std::ios::out);
      Correlations->WriteDataFile(LogFile);
    }
  if(Correlations!=NULL)
    {
      ofstream LogFile("density_radial.dat",std::ios::out);
      Density->WriteDataFile(LogFile);
    }

}

