#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEMonteCarlo/AbstractMCSamplingFunction.h"
#include "Tools/FQHEMonteCarlo/QHESamplingFunctionManager.h"
#include "Tools/FQHEMonteCarlo/SimpleMonteCarloOnSphereAlgorithm.h"
#include "Tools/FQHEMonteCarlo/SphereBilayerCoulombEnergy.h"
#include "Tools/FQHEMonteCarlo/SphereWithSpinGeneralEnergy.h"

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
  OptionManager Manager ("FQHESphereWithSpinMonteCarlo" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");      
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
    
  QHEWaveFunctionManager WaveFunctionManager(QHEWaveFunctionManager::SphereWithSpinGeometry);
  QHESamplingFunctionManager SamplingFunctionManager(QHESamplingFunctionManager::SphereWithSpinGeometry);
  
  Manager += SystemGroup;
  SimpleMonteCarloOnSphereAlgorithm::AddOptionGroup(&Manager);
  WaveFunctionManager.AddOptionGroup(&Manager);
  SamplingFunctionManager.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum angular momentum", 9);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the total spin", 0);
  
  (*SystemGroup) += new SingleDoubleOption  ('D', "lowest-d", "smallest layer separation d where the energy is evaluated", 0.0, true, /* minimum value */ 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "spacing-d", "spacing of layer separation d where the energy is evaluated", 0.25, true, /* minimum value */ 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-d", "number of layer separation d where the energy is evaluated", 13);
  (*SystemGroup) += new SingleStringOption('\n',"all-params","File containing all parameters of the interaction");
  (*SystemGroup) += new SingleStringOption('\n',"params-inter","File containing parameters of inter-spin interaction");
  (*SystemGroup) += new SingleStringOption('\n',"params-intra","File containing parameters of intra-spin interaction");
  (*SystemGroup) += new SingleStringOption('\n',"params-other","File containing parameters of second intra-spin interaction for down spins (if different from up channel)");
  
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "list-samplingfunctions", "list all available sampling-fuctions");

  (*SystemGroup) += new BooleanOption ('b', "background-only", "print background energy, and exit");  
  (*SystemGroup) += new SingleIntegerInternalOption  ('n', "SzTotal", "number of layer separation d where the energy is evaluated", 0);

  (*SystemGroup) += new SingleStringOption('o',"output","file to write a log of results");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout); // new feature to process options and display help in 1 line!

  if (Manager.GetBoolean("list-wavefunctions") == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }
  
  SamplingFunctionManager.TestForShow();
  

  int NbrParticles = Manager.GetInteger("nbr-particles");
  double LowestD = Manager.GetDouble("lowest-d");
  double Spacing = Manager.GetDouble("spacing-d");
  int NbrSeparations = Manager.GetInteger("nbr-d");
  int SzTotal = Manager.GetInteger("total-sz");
  int NbrFlux = Manager.GetInteger("lzmax");
  int NbrUp = (NbrParticles + SzTotal)/2;
  
  
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

  // add observables
  if ((Manager.GetString("params-inter")==0)&&(Manager.GetString("all-params")==0))
    {
      // add Coulomb energy for bilayer      
      SphereBilayerCoulombEnergy *Energy = new SphereBilayerCoulombEnergy(/*Flux */ NbrParticles-1, NbrSeparations, LowestD, Spacing);
      MonteCarloRoutine.AddObservable(Energy);
    }
  else
    {
      // add a general energy function
      SphereWithSpinGeneralEnergy *Energy;
      if (Manager.GetString("all-params")!=0)
	{
	  Energy=new SphereWithSpinGeneralEnergy(NbrUp, NbrFlux, Manager.GetString("all-params"));
	}
      else
	{
	  if (Manager.GetString("params-intra")==0)
	    {
	      cout << "Attention: both inter- and intra-spin interactions are required!"<<endl;
	      exit(-1);
	    }
	  Energy=new SphereWithSpinGeneralEnergy(NbrUp, NbrFlux,
						 Manager.GetString("params-inter"),
						 Manager.GetString("params-intra"),
						 Manager.GetString("params-other"));
	}
      MonteCarloRoutine.AddObservable(Energy);
      if (Manager.GetBoolean("background-only"))
	{
	  cout << "# Background energy"<<endl;
	  cout << "bg "<<Energy->GetTotalBackgroundEnergy()<<endl;
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

  // print final results:
  cout << "Final results:" << endl;
  MonteCarloRoutine.WriteObservations(cout);

  if (Manager.GetString("output")!=NULL)
    {
      ofstream LogFile(Manager.GetString("output"),std::ios::out);
      MonteCarloRoutine.WriteObservations(LogFile);
      LogFile.close();
    }

  if (TestWaveFunction!=NULL) delete TestWaveFunction;
  if (SamplingFunction!=NULL) delete SamplingFunction;
}

