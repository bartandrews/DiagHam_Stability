#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeNonPeriodicMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnHyperCubicLatticeWithSpinMomentumSpace.h"

#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWannierSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace.h"

#include "Vector/ComplexVector.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("FQHETopInsulatorShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('t', "nbr-sitet", "number of sites along the t direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kx", "total momentum along the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "total momentum along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kz", "total momentum along the z direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kt", "total momentum along the t direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "4d", "consider a 4d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption  ('\n', "spin-conserved", "assume that the spin is conserved in the two band model");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz", "twice the spin Sz value (only useful when --spin-conserved is on)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "non-periodic", "look at the non-periodic hilbert space with a cut in momentum space described by nbr-allowed-site and min-k");
  (*SystemGroup) += new SingleIntegerOption  ('X', "nbr-allowed-sitex", "number of x momenta allowed for a single particle", 3);
  (*SystemGroup) += new SingleIntegerOption  ('Y', "nbr-allowed-sitey", "number of y momenta allowed for a single particle", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-kx", "minimal x momentum allowed for a single particle", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-ky", "minimal y momentum allowed for a single particle", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-kz", "minimal z momentum allowed for a single particle", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-kt", "minimal t momentum allowed for a single particle", 4);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "wannier", "use wannier wavefunction basis");
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-autodetect", "do not autdetect system parameter from state file name");
  (*SystemGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (not available for the non-periodic momentum space or the case with spin)",0);
  //  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  //  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_topinsulator_nbrsubbands_n_nbrparticles_x_nbrsitex_y_nbrsitey.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  int NbrSiteZ = Manager.GetInteger("nbr-sitez"); 
  int NbrSiteT = Manager.GetInteger("nbr-sitet"); 
  int NbrAllowedKx = Manager.GetInteger("nbr-allowed-sitex"); 
  int NbrAllowedKy = Manager.GetInteger("nbr-allowed-sitey"); 
  int MinKx = Manager.GetInteger("min-kx"); 
  int MinKy = Manager.GetInteger("min-ky");
  int MinKz = Manager.GetInteger("min-kz");
  int MinKt = Manager.GetInteger("min-kt");
  int TotalKx = Manager.GetInteger("kx"); 
  int TotalKy = Manager.GetInteger("ky");
  int TotalKz = Manager.GetInteger("kz");
  int TotalKt = Manager.GetInteger("kt");
  int Sz = Manager.GetInteger("sz");

  if ((Manager.GetString("state") != 0) && (Manager.GetBoolean("no-autodetect") == false))
    {
      double Mass = 0.0;
      bool Statistics = false;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							      NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy, Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	  return -1;
	}	  
    }
 
  AbstractQHEParticle* Space;
  if (Manager.GetBoolean("boson") == false)
    {
      if (Manager.GetBoolean("non-periodic") == true)
	{
	  if (Manager.GetInteger("nbr-subbands") == 1)
	    Space = new FermionOnSquareLatticeNonPeriodicMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrAllowedKx, NbrAllowedKy, MinKx, MinKy, TotalKx, TotalKy);
	  else
	    return 0;
	}
      else
	{
	  if (Manager.GetInteger("nbr-subbands") == 1)
	    {
              if (Manager.GetBoolean("3d") == false)
                {
                  Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
                  if (Manager.GetString("save-hilbert") != 0)
                    {
                      Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
                      return 0;
                    }
                }
              else
                {
                  Space = new FermionOnCubicLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx, TotalKy, TotalKz);
                }
	    }
	  else
	    {
	      if (Manager.GetInteger("nbr-subbands") == 2)
		{
		  if (Manager.GetBoolean("3d") == false)
		    {
		      if (Manager.GetBoolean("4d") == false)
			{
			  if (Manager.GetBoolean("spin-conserved") == false)
			    {
			      Space = new FermionOnSquareLatticeWithSpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
			    }
			  else
			    {
			      Space = new FermionOnSquareLatticeWithSpinMomentumSpace(NbrParticles, (Sz + NbrParticles) / 2, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
			    }
			}
		      else
			{
			  Space = new FermionOnHyperCubicLatticeWithSpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, NbrSiteT, TotalKx, TotalKy, TotalKz, TotalKt);
			}
		    }
		  else
		    {
		      Space = new FermionOnCubicLatticeWithSpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx, TotalKy, TotalKz);
		    }		
		}
	      else
		{
		  if (Manager.GetInteger("nbr-subbands") == 4)
		    {
		      if (Manager.GetBoolean("3d") == true)
			{
			  if (NbrSiteX*NbrSiteY*NbrSiteZ <= 15)
			  {
			    Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx, TotalKy, TotalKz);
			  }
			  else
			  {
			     Space = new FermionOnCubicLatticeWithSU4SpinMomentumSpaceLong(NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx, TotalKy, TotalKz);
			  }
			}
		      else 
			{
			  Space = new FermionOnSquareLatticeWithSU4SpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
			}
		    }
		  else
		    {
		      return 0;
		    }
		}
	    }
	}
    }
  else
    {
      if (Manager.GetInteger("nbr-subbands") == 1)
	{
	  if (Manager.GetBoolean("3d") == false)
	    {
              if (Manager.GetBoolean("wannier"))
                {
                  Space = new BosonOnSquareLatticeWannierSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKy);
                }
              else
                {
                  Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
                }
              if (Manager.GetString("save-hilbert") != 0)
                {
                  Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
                  return 0;
                }
            }
          else
            {
	      Space = new BosonOnCubicLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, TotalKx, TotalKy, TotalKz);
            }
	}
      if (Manager.GetInteger("nbr-subbands") == 2)
	{
	  if (Manager.GetBoolean("3d") == false)
	    {
	      if (Manager.GetBoolean("4d") == false)
		{
		  if (Manager.GetBoolean("spin-conserved") == false)
		    {
		      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
		    }
		  else
		    {
		      Space = new BosonOnSquareLatticeWithSU2SpinMomentumSpace (NbrParticles, (Sz + NbrParticles) / 2, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
		    }
		}
	      else
		{
		  Space = new BosonOnHyperCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, NbrSiteT,
										TotalKx, TotalKy, TotalKz, TotalKt);
		}
	    }
	  else
	    {
	      Space = new BosonOnCubicLatticeWithSU2SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ,
								       TotalKx, TotalKy, TotalKz);
	    }
	}
      if (Manager.GetInteger("nbr-subbands") == 3)
	{
	  Space = new BosonOnSquareLatticeWithSU3SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
	}
      if (Manager.GetInteger("nbr-subbands") == 4)
	{
	  if (Manager.GetBoolean("3d") == false)
	    {
	      Space = new BosonOnSquareLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
	    }
	  else
	    {
	      if (NbrParticles + NbrSiteX*NbrSiteY*NbrSiteZ <= 64)
	      {
		Space = new BosonOnCubicLatticeWithSU4SpinMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ,
								       TotalKx, TotalKy, TotalKz);
	      }
	      else
	      {
		Space = new BosonOnCubicLatticeWithSU4SpinMomentumSpaceLong (NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ,
								       TotalKx, TotalKy, TotalKz);
	      }
	    }
	}
      if (Manager.GetString("save-hilbert") != 0)
	{
	  Space->WriteHilbertSpace(Manager.GetString("save-hilbert"));
	  return 0;
	}
    }

  if (Manager.GetString("state") == 0)
    {
          cout<< Space->GetHilbertSpaceDimension()<<endl;
      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
	Space->PrintState(cout, i) << endl;;
      cout << endl;
    }
  else
    {
      int NbrHiddenComponents = 0;
      double WeightHiddenComponents = 0.0;
      double Normalization = 0.0;
      ComplexVector State;
      if (State.ReadVector(Manager.GetString("state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("state") << endl;
	  return -1;
	}
      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
	{
	  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
      if (Manager.GetDouble("hide-component") > 0.0)
	{
	  double Error = Manager.GetDouble("hide-component");
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      if (Norm(State[i]) > Error)
		Space->PrintState(cout, i) << " : "  << State[i] << " Norm : "<< Norm(State[i]) << " Phase : " <<  Arg(State[i]) << endl;
	      else
		{
		  WeightHiddenComponents += SqrNorm(State[i]);
		  NbrHiddenComponents++;
		}
	      Normalization += SqrNorm(State[i]);
	    }
	  cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
	}
      else
      {
	for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	  Space->PrintState(cout, i) << " : "  << State[i] << endl;;
      }
    }
  delete Space;
  return 0;
}

