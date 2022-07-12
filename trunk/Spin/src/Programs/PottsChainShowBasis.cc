#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Potts3Chain.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PottsChainShowBasis" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('n', "zn", "Zn value of the Potts model", 3);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 2);
  (*SystemGroup) += new  SingleIntegerOption ('q', "q-value", "Zn charge", 0);
  (*SystemGroup) += new SingleStringOption  ('e', "state", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PottsChainShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int ZnValue = Manager.GetInteger("zn");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  int QValue = Manager.GetInteger("q-value");
  double Error = Manager.GetDouble("hide-component");
  if (Manager.GetString("state") != 0)
    {
       if (PottsFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrSpins, QValue) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	  return -1;
	}    
    }

  AbstractSpinChain* Space = 0;
  switch (ZnValue)
    {
    case 2 :
      Space = new Spin1_2ChainFixedParity (NbrSpins, QValue);
      break;
    case 3 :
      Space = new Potts3Chain (NbrSpins, QValue, 1000000);
      break;
    default :
      {
	cout << "Zn with n=" << ZnValue  << " is not available" << endl;
	return -1;
      }
    }
  
  if (Manager.GetString("state") == 0)
    {
      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
	Space->PrintState(cout, i) << endl;
    }
  else
    {
      RealVector State;
      if (State.ReadVectorTest(Manager.GetString("state")) == true)
	{
	  if (State.ReadVector(Manager.GetString("state")) == false)
	    {
	      cout << "error while reading " << Manager.GetString("state") << endl;
	      return -1;
	    }
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    if (fabs(State[i]) > Error)
	      Space->PrintState(cout, i) << " : "  << State[i] << endl;
	}
      else
	{
	  ComplexVector ComplexState;
	  if (ComplexState.ReadVector(Manager.GetString("state")) == false)
	    {
	      cout << "error while reading " << Manager.GetString("state") << endl;
	      return -1;
	    }
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    if (Norm(ComplexState[i]) > Error)
	      Space->PrintState(cout, i) << " : "  << ComplexState[i] << endl;
	}
    }
  return 0;
}
