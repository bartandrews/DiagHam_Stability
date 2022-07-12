#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEMonteCarlo/AbstractMCSamplingFunction.h"
#include "Tools/FQHEMonteCarlo/QHESamplingFunctionManager.h"
#include "Tools/FQHEMonteCarlo/SimpleMonteCarloOnSphereAlgorithm.h"
#include "Tools/FQHEMonteCarlo/SphereBilayerCoulombEnergy.h"

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
  (*SystemGroup) += new SingleDoubleOption  ('D', "lowest-d", "smallest layer separation d where the energy is evaluated", 0.0, true, /* minimum value */ 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('s', "spacing-d", "spacing of layer separation d where the energy is evaluated", 0.25, true, /* minimum value */ 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-d", "number of layer separation d where the energy is evaluated", 13);
  
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "list-samplingfunctions", "list all available sampling-fuctions");
  (*SystemGroup) += new SingleIntegerInternalOption  ('n', "SzTotal", "number of layer separation d where the energy is evaluated", 0);
  
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
  
  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  AbstractMCSamplingFunction* SamplingFunction = SamplingFunctionManager.GetSamplingFunction();

  cout << "Function: " << WaveFunctionManager.GetDescription()<<endl;
  cout << "Sampler:  " << SamplingFunctionManager.GetDescription()<<endl;

  SimpleMonteCarloOnSphereAlgorithm MonteCarloRoutine(NbrParticles, TestWaveFunction, SamplingFunction,
						      &Manager);

  // add observables
  SphereBilayerCoulombEnergy Energy(/*Flux */ NbrParticles-1, NbrSeparations, LowestD, Spacing);
  
  MonteCarloRoutine.AddObservable(&Energy);  

  // run simulation
  MonteCarloRoutine.Simulate();

  // print final results:
  cout << "Final results:" << endl;
  Energy.WriteDataFile(cout);  
  
}

