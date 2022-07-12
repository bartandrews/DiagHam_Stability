#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEMonteCarlo/AbstractMCSamplingFunction.h"
#include "Tools/FQHEMonteCarlo/QHESamplingFunctionManager.h"
#include "Tools/FQHEMonteCarlo/SimpleMonteCarloOnSphereAlgorithm.h"
#include "Tools/FQHEMonteCarlo/SphereCoulombEnergy.h"
#include "Tools/FQHEMonteCarlo/SphereGeneralEnergy.h"

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
  OptionManager Manager ("FQHESphereMonteCarlo" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
    
  QHEWaveFunctionManager WaveFunctionManager(QHEWaveFunctionManager::SphereGeometry);
  QHESamplingFunctionManager SamplingFunctionManager(QHESamplingFunctionManager::SphereGeometry);
  
  Manager += SystemGroup;
  SimpleMonteCarloOnSphereAlgorithm::AddOptionGroup(&Manager);
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

  if (TestWaveFunction==0)
    {
      cout << "Invalid wavefunction requested - list available states with --list-wavefunctions"<<endl;
      exit(-1);
    }
  cout << "Function: " << WaveFunctionManager.GetDescription()<<endl;
  if (SamplingFunction!=0)
    cout << "Sampler:  " << SamplingFunctionManager.GetDescription()<<endl;

  SimpleMonteCarloOnSphereAlgorithm MonteCarloRoutine(NbrParticles, TestWaveFunction, SamplingFunction,
						      &Manager);
  
  if (Manager.GetString("interaction-params")==0)
    {
      // add Coulomb energy in LLL
      SphereCoulombEnergy *Energy=new SphereCoulombEnergy(/*Flux */ NbrFlux);
      MonteCarloRoutine.AddObservable(Energy);
    }
  else
    {
      // add a general energy function 
      SphereGeneralEnergy *Energy=new SphereGeneralEnergy(NbrFlux, Manager.GetString("interaction-params"));
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
  
}

