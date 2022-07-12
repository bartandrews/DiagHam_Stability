#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFixedParity.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetryLong.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetriesLong.h"

#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"


#include "HilbertSpace/DoubledSpin1_2_Chain.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslations.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslationsStaggered.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsStaggered.h"
//#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslations.h"
// #include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers.h"
#include "HilbertSpace/DoubledSpin1_2_ChainWithTranslations.h"
#include "HilbertSpace/DoubledSpin1_2_ChainWithTranslations_alternative.h"

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
  OptionManager Manager ("GenericSpinChainShowBasis" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('z', "sz-value", "twice the value of sz", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "fixed-parity", "the spin chain has a fixed total parity");
  (*SystemGroup) += new  BooleanOption ('\n', "doubled-spinchain", "use the double spin chain HilbertSpace");
  (*SystemGroup) += new  BooleanOption ('\n', "staggered", "use the double spin chain HilbertSpace");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "parity", "parity of the spin chain", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "periodic-chain", "consider periodic instead of open chain", false);
  (*SystemGroup) += new  BooleanOption ('\n', "no-sz", "don't consider sz quantum number (only available for certain chains up to now)", false);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "max-momentum", "max momentum in the periodic case (should be a divisor of the number of sites)", -1);  
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sz-symmetry", "set the Sz<->-Sz symmetry sector (0 if it should not be used)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "inversion-symmetry", "set the inversion symmetry sector (0 if it should not be used)", 0);

  (*SystemGroup) += new  BooleanOption ('\n', "zero-half", "consider the 0 +1/2 spin chain", false);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetry", "use Hilbert space with ZZ symmetry in the case of the 0+1/2 doubled chain");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "zket", "value of the Z operator in the ket layer", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "zbra", "value of the Z operator in the bra layer", 0);

  (*SystemGroup) += new BooleanOption  ('\n', "sublatticeQN", "use Hilbert space with sublattice quantum numbers (implies the use of ZZ symmetry) in the case of the 0+1/2 doubled chain");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sket", "value of the Sublattice QN in the ket layer", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sbra", "value of the Sublattice QN the bra layer", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sproduct", "value of the product of the Sublattice QNs", 0);
  
  (*SystemGroup) += new  SingleIntegerOption ('k', "momentum", "momentum sector (for periodic chain)", 0);
  (*SystemGroup) += new SingleStringOption  ('e', "state", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new BooleanOption ('\n', "complex-vector", "the eigenstate to be  displayed is a complex vector");

  (*SystemGroup) += new BooleanOption  ('\n', "save-normalized-vector", "Save the phase-fixed vector for eigenvectors of PEPS transfer matrix");

  (*SystemGroup) += new BooleanOption  ('\n', "unnormalize", "remove the symmetry factors from the basis (i.e. the sqrt(orbit size) factors)");
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericSpinChainShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");
  int SzValue = Manager.GetInteger("sz-value");
  int Momentum = Manager.GetInteger("momentum");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  if (MaxMomentum == -1) 
    MaxMomentum = NbrSpins;
  double Error = Manager.GetDouble("hide-component");
  bool SymmetryFlag = Manager.GetBoolean("symmetry"); 
  bool SubLatticeQuantumNumberFlag =  Manager.GetBoolean("sublatticeQN"); 
  bool NoSzQuantumNumberFlag =  Manager.GetBoolean("no-sz"); 
  if  (SubLatticeQuantumNumberFlag == true )
    SymmetryFlag = true;
  
  if (Manager.GetBoolean("doubled-spinchain"))
    {
      if( (SpinValue==1 ) && ( Manager.GetBoolean("zero-half") == false ) )
	{
	  if (Manager.GetBoolean("periodic-chain") == false)
	    {
	      DoubledSpin1_2_Chain Space (NbrSpins, SzValue,  1000000, 1000000);
	      
	      if (Manager.GetString("state") == 0)
		{
		  for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
		    Space.PrintState(cout, i) << endl;
		}
	      else
		{
		  RealVector State;
		  if (State.ReadVector(Manager.GetString("state")) == false)
		    {
		      cout << "error while reading " << Manager.GetString("state") << endl;
		      return -1;
		    }
		  
		  if (Space.GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
		    {
		      cout << "dimension mismatch between Hilbert space and ground state" << endl;
		      return 0;
		    }
		  
		  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
		    if (fabs(State[i]) > Error)
		      Space.PrintState(cout, i) << " : "  << State[i] << endl;
		}
	    }
	  else
	    {
	      DoubledSpin1_2_ChainWithTranslations  Space (NbrSpins,  Momentum,  SzValue,  1000000, 1000000);
	      DoubledSpin1_2_ChainWithTranslations_alternative  Space1 (NbrSpins,  Momentum,  SzValue,  1000000, 1000000);
	      cout << "dim = " <<  Space1.GetHilbertSpaceDimension() << endl ;
	      if (Manager.GetString("state") == 0)
		{
		  for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
		  {
		    Space.PrintState(cout, i) << endl;
// 		    cout << i << " " ;
		    Space1.PrintState(cout, i) << endl;
		    cout << endl;
		  }
		}
	      else
		{
		  RealVector State;
		  if (State.ReadVector(Manager.GetString("state")) == false)
		    {
		      cout << "error while reading " << Manager.GetString("state") << endl;
		      return -1;
		    }
		  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
		    if (fabs(State[i]) > Error)
		      Space.PrintState(cout, i) << " : "  << State[i] << endl;
		  
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("zero-half")==true)
	    {
	      if (Manager.GetBoolean("periodic-chain") == false)
		{
		  AbstractSpinChain* Space = 0;
		  if (Manager.GetBoolean( "staggered") == true)
		    {
		      if (SymmetryFlag == false )
			{
			  Space = new  DoubledSpin0_1_2_ChainWithTranslationsStaggered ( NbrSpins, SzValue,  1000000, 1000000);
			}
		      else
			{
			  Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry(NbrSpins, SzValue ,Manager.GetInteger("zbra"),Manager.GetInteger("zket"),10000,10000);
			}
		    }
		  else
		    {
		      if (SymmetryFlag == false )
			{
			  Space = new DoubledSpin0_1_2_ChainWithTranslations ( NbrSpins, SzValue,  1000000, 1000000);
			}
		      else
			{
			  if (SubLatticeQuantumNumberFlag == false ) 
			    Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry(NbrSpins,  SzValue, Manager.GetInteger("zbra"),Manager.GetInteger("zket"),10000,10000);
			  else
			    Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (NbrSpins,  SzValue, Manager.GetInteger("zbra"),Manager.GetInteger("zket"),  Manager.GetInteger("sbra") *  Manager.GetInteger("sbra"), Manager.GetInteger("sket") *  Manager.GetInteger("sket"),  Manager.GetInteger("sket")* Manager.GetInteger("sbra")* Manager.GetInteger("sproduct"), 100000,100000);
			}
		    }
		  
		  if (Manager.GetString("state") == 0)
		    {
		      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
			Space->PrintState(cout, i) << endl;
		    }
		  else
		    {
		      if (Manager.GetBoolean("complex-vector") == false)
			{
			  RealVector State;
			  if (State.ReadVector(Manager.GetString("state")) == false)
			    {
			      cout << "error while reading " << Manager.GetString("state") << endl;
			      return -1;
			    }
			  if (Space->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
			    {
			      cout << "dimension mismatch between Hilbert space and ground state" << endl;
			      return 0;
			    }
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (fabs(State[i]) > Error)
			      Space->PrintState(cout, i) << " : "  << State[i] << endl;
			}
		      else
			{
			  ComplexVector State;
			  if (State.ReadVector(Manager.GetString("state")) == false)
			    {
			      cout << "error while reading " << Manager.GetString("state") << endl;
			      return -1;
			    }
			  if (Space->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
			    {
			      cout << "dimension mismatch between Hilbert space and ground state" << endl;
			      return 0;
			    }
			  ((DoubledSpin0_1_2_ChainWithTranslations *) Space)->NormalizeDensityMatrix(State);
			  if ( Manager.GetBoolean("save-normalized-vector") == true ) 
			    {
			      State.WriteVector(Manager.GetString("state"));
			    }
			  
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (Norm(State[i]) > Error)
			      Space->PrintState(cout, i) << " : "  << State[i] << endl;
			}
		    }
		}
	      else
		{
		  AbstractSpinChain* Space = 0;
		  if (Manager.GetBoolean( "staggered") == true)
		    {
		      if (SymmetryFlag == false )
			{
			  Space = new  DoubledSpin0_1_2_ChainWithTranslationsStaggered (NbrSpins,Momentum, SzValue,  1000000, 1000000);
			}
		      else
			{
			  Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry(NbrSpins, Momentum,  SzValue, Manager.GetInteger("zbra"),Manager.GetInteger("zket"),10000,10000);
			}
		    }
		  else
		    {
		      if (SymmetryFlag == false )
			{
			  Space = new DoubledSpin0_1_2_ChainWithTranslations (NbrSpins, Momentum, NbrSpins /MaxMomentum  ,SzValue,  1000000, 1000000);
			}
		      else
			{
			  if (SubLatticeQuantumNumberFlag == false ) 
			    Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry(NbrSpins, Momentum, NbrSpins /MaxMomentum ,  SzValue, Manager.GetInteger("zbra"),Manager.GetInteger("zket"),10000,10000);
			  else
			    {
			      Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (NbrSpins, Momentum, NbrSpins /MaxMomentum  , SzValue,  Manager.GetInteger("sket") *  Manager.GetInteger("sket"), Manager.GetInteger("sbra") *  Manager.GetInteger("sbra"), Manager.GetInteger("sket")* Manager.GetInteger("sbra")* Manager.GetInteger("sproduct"), 100000,100000);
			    }
			}
		    }
		  
		  
		  if (Manager.GetString("state") == 0)
		    {
		      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
			Space->PrintState(cout, i) << endl;
		    }
		  else
		    {
		      if (Manager.GetBoolean("complex-vector") == false)
			{
			  RealVector State;
			  if (State.ReadVector(Manager.GetString("state")) == false)
			    {
			      cout << "error while reading " << Manager.GetString("state") << endl;
			      return -1;
			    }
			  if (Space->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
			    {
			      cout << "dimension mismatch between Hilbert space and ground state" << endl;
			      return 0;
			    }
 
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (fabs(State[i]) > Error)
			      Space->PrintState(cout, i) << " : "  << State[i] << endl;
			}
		      else
			{
			  ComplexVector State;
			  if (State.ReadVector(Manager.GetString("state")) == false)
			    {
			      cout << "error while reading " << Manager.GetString("state") << endl;
			      return -1;
			    }
			  if (Space->GetLargeHilbertSpaceDimension() != State.GetLargeVectorDimension())
			    {
			      cout << "dimension mismatch between Hilbert space and ground state" << endl;
			      return 0;
			    }
			  ((DoubledSpin0_1_2_ChainWithTranslations *) Space)->NormalizeDensityMatrix(State);
			  
			  if ( Manager.GetBoolean("save-normalized-vector") == true ) 
			    {
			      State.WriteVector(Manager.GetString("state"));
			    }
			  
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (Norm(State[i]) > Error)
			      Space->PrintState(cout, i) << " : "  << State[i] << endl;
			}
		    }
		}
	      return 0;
	    }
	}
    }
  
  if (Manager.GetBoolean("zero-half") == true)
    {
      if (Manager.GetBoolean("periodic-chain") == false)
	{
	  AbstractSpinChain* Space = 0;

	  if (Manager.GetBoolean( "staggered") == true)
	    {
	      Space = new Spin0_1_2_ChainWithTranslationsStaggered ( NbrSpins, SzValue,  1000000, 1000000);
	    }
	  else
	    {
	      if (NoSzQuantumNumberFlag == true )
		Space = new Spin0_1_2_ChainWithTranslations ( NbrSpins,  1000000, 1000000);
	      else
		Space = new Spin0_1_2_ChainWithTranslations ( NbrSpins, SzValue,  1000000, 1000000);
	    }
	  
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("complex-vector") == false)
		{
		  RealVector State;
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
		  ComplexVector State;
		  if (State.ReadVector(Manager.GetString("state")) == false)
		    {
		      cout << "error while reading " << Manager.GetString("state") << endl;
		      return -1;
		    }
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (Norm(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;
		}
	    }
	  return 0;
	}
      else
	{
	  AbstractSpinChain* Space = 0;
	  if (	  SubLatticeQuantumNumberFlag == false)
	    {
	      if (Manager.GetBoolean( "staggered") == true)
		{
		  Space = new Spin0_1_2_ChainWithTranslationsStaggered ( NbrSpins,  Momentum, SzValue,  1000000, 1000000);
		}
	      else
		{
		  Space = new Spin0_1_2_ChainWithTranslations ( NbrSpins,  Momentum,  NbrSpins / MaxMomentum, SzValue,  1000000, 1000000);
		}
	    }
	  else
	    {
	      Space = new Spin0_1_2_ChainWithTranslationsAndSublatticeQuantumNumbers  ( NbrSpins,  Momentum,  NbrSpins / MaxMomentum, SzValue,  Manager.GetInteger("sbra"),  1000000, 1000000);
	    }
	      
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("complex-vector") == false)
		{
		  RealVector State;
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
		  ComplexVector State;
		  if (State.ReadVector(Manager.GetString("state")) == false)
		    {
		      cout << "error while reading " << Manager.GetString("state") << endl;
		      return -1;
		    }
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (Norm(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;
		}
	    }
	}
      return 0;
    }
  
  
  if (Manager.GetBoolean("periodic-chain") == false)
    {
      AbstractSpinChain* Space = 0;
      switch (SpinValue)
	{
	case 1 :
	  {
	    if (NoSzQuantumNumberFlag == false)
	      {
		if (Manager.GetBoolean("fixed-parity") == false)
		  Space = new Spin1_2Chain (NbrSpins, SzValue, 1000000);
		else
		  Space = new Spin1_2ChainFixedParity (NbrSpins, Manager.GetInteger("parity"));
	      }
	    else
	      {
		Space = new Spin1_2ChainFull (NbrSpins);
	      }
	  }
	  break;
	case 2 :
	  Space = new Spin1Chain (NbrSpins, SzValue, 1000000);
	  break;
	default :
	  {
	    if ((SpinValue & 1) == 0)
	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
	    else 
	      cout << "spin " << SpinValue << "/2 are not available" << endl;
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
	 if (Manager.GetBoolean("complex-vector") == false)
	   {
	     RealVector State;
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
	     ComplexVector State;
	     if (State.ReadVector(Manager.GetString("state")) == false)
	       {
	      cout << "error while reading " << Manager.GetString("state") << endl;
	      return -1;
	       }
	     for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	       if (Norm(State[i]) > Error)
		 Space->PrintState(cout, i) << " : "  << State[i] << endl;
	   }
        }
     }
   else //periodic bc
    {
      AbstractSpinChainWithTranslations* Space = 0;
      switch (SpinValue)
	{
	case 1 :
	  {
	    if ((Manager.GetInteger("inversion-symmetry") != 0) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins / 2)))))
	      {
		if ((Manager.GetInteger("sz-symmetry") != 0) && (SzValue == 0))
		  {
		    Space = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momentum, 1, Manager.GetInteger("inversion-symmetry"), 
										      Manager.GetInteger("sz-symmetry"), SzValue, 1000000, 1000000);
		  }
		else
		  {
		    if (NbrSpins < 32)
		      {
			Space = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momentum, 1, Manager.GetInteger("inversion-symmetry"), SzValue, 1000000, 1000000);
		      }
		    else
		      {
			Space = new Spin1_2ChainWithTranslationsAndInversionSymmetryLong (NbrSpins, Momentum, 1, Manager.GetInteger("inversion-symmetry"), SzValue, 1000000, 1000000);
		      }
		  }
	      }
	    else
	      {
		if ((Manager.GetInteger("sz-symmetry") != 0) && (SzValue == 0))
		  {
		    Space = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSpins, Momentum, 1, Manager.GetInteger("sz-symmetry"), SzValue, 1000000, 1000000);
		  }
		else
		  {
		    if (NbrSpins < 32)
		      {
			Space = new Spin1_2ChainWithTranslations (NbrSpins, Momentum, 1, SzValue, 1000000, 1000000);
		      }
		    else
		      {
			Space = new Spin1_2ChainWithTranslationsLong (NbrSpins, Momentum, 1, SzValue, 1000000, 1000000);
		      }
		  }
	      }
	  }
	  break;
	case 2 :
	  {
	    if ((Manager.GetInteger("sz-symmetry") == 0) || (SzValue != 0))
		{
		  if ((Manager.GetInteger("inversion-symmetry") != 0) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins / 2)))))
		    {
		      Space = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, Momentum, Manager.GetInteger("inversion-symmetry"), SzValue);
		    }
		  else
		    {
		      Space = new Spin1ChainWithTranslations (NbrSpins, Momentum, SzValue);
		    }
		}
	      else
		{
		  if ((Manager.GetInteger("inversion-symmetry") != 0) && ((Momentum == 0) || (((NbrSpins & 1) == 0) && (Momentum == (NbrSpins / 2)))))
		    {
		      Space = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, Momentum, Manager.GetInteger("inversion-symmetry"), 
										      Manager.GetInteger("sz-symmetry"), SzValue);
		    }
		  else
		    {
		      Space = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, Momentum, Manager.GetInteger("sz-symmetry"), SzValue);
		    }
		}
	    }
	    break;
	  default :
	    {
	        if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		else 
		  cout << "spin " << SpinValue << "/2 are not available" << endl;
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
 	 if (Manager.GetBoolean("complex-vector") == false)
	   {
	     RealVector State;
	     if (State.ReadVector(Manager.GetString("state")) == false)
	       {
		 cout << "error while reading " << Manager.GetString("state") << endl;
		 return -1;
	       }
	     if (Manager.GetBoolean("unnormalize") == false)
	       {
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   {
		     if (fabs(State[i]) > Error)
		       Space->PrintState(cout, i) << " : "  << State[i] << endl;
		   }
	       }
	     else
	       {
		 double* NormalizationFactors = Space->GetBasisNormalization();
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   {
		     if (fabs(State[i]) > Error)
		       Space->PrintState(cout, i) << " : "  << State[i] << " orb size=" << 1.0 / (NormalizationFactors[i] * NormalizationFactors[i]) << endl;
		   }
	       }
	   }
	 else
	   {
	     ComplexVector State;
	     if (State.ReadVector(Manager.GetString("state")) == false)
	       {
		 cout << "error while reading " << Manager.GetString("state") << endl;
		 return -1;
	       }
	     if (Manager.GetBoolean("unnormalize") == false)
	       {
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   {
		     if (Norm(State[i]) > Error)
		       Space->PrintState(cout, i) << " : "  << State[i] << endl;
		   }
	       }
	     else
	       {
		 double* NormalizationFactors = Space->GetBasisNormalization();
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   {
		     if (Norm(State[i]) > Error)
		       Space->PrintState(cout, i) << " : "  << (State[i] * NormalizationFactors[i]) << endl;
		   }
	       }
	   }
       }


    }
      
  return 0;
}
