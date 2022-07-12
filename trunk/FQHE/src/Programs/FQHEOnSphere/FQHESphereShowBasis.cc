#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinPartialPolarization.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinPartialPolarization.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinAllLz.h"
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereWithSpinTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOn4DSphere.h"
#include "HilbertSpace/BosonOn4DSphereLong.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"
#include "HilbertSpace/FermionOnCP2.h"
#include "HilbertSpace/FermionOnCP2Long.h"
#include "HilbertSpace/FermionOnS2xS2.h"
#include "HilbertSpace/BosonOnS2xS2.h"
#include "HilbertSpace/BosonOnS2xS2Long.h"
#include "HilbertSpace/BosonOnS2xS2HardcoreNoNearestNeighbors.h"
#include "HilbertSpace/BosonOnS2xS2HardcoreNoNearestNeighborsLong.h"
#include "HilbertSpace/FermionOnS2xS2WithExclusionPrinciple.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/FilenameTools.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrflux-2", "number of flux quanta for the second sphere on the S2XS2 geometry", 20);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "jz-value", "twice the total value of jz (only useful in the bosonic 4D mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kz-value", "twice the total value of kz (only useful in the bosonic 4D mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "tz-value", "the total value of Tz (only useful in the bosonic CP2 mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-value", "the total value of Y (only useful in the bosonic CP2 mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-lz", "consider particles with all Lz components");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for the SU(2) spinful bosonic states");
  (*SystemGroup) += new BooleanOption  ('\n', "4-D", "consider particles on the 4D sphere (only available in the bosonic mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "cp2", "consider particles on the CP2 ");
  (*SystemGroup) += new BooleanOption  ('\n', "s2xs2", "consider particles on the S2xS2 geometry");
  (*SystemGroup) += new BooleanOption  ('\n', "s2s2-hardcorenonn", "consider particles on the S2xS2 geometry, with hardcore contraint and no nearest neighbor");
  (*SystemGroup) += new BooleanOption  ('\n', "truncated-cp2", "consider particles on a truncated CP2 geometry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-y", "minimum value of y for a truncated CP2 geometry", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "tzZ3symmetrized-basis", "use Tz <-> -Tz and Z3 permutations symmetrized version of the CP2 basis (only valid if total-tz=0 and total-y = 0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2)/su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "all-sz", "consider particles with SU(2) spin all Sz components");
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "add index of the Hilbert space vectors");
  (*SystemGroup) += new BooleanOption  ('\n', "add-szvalue", "add Sz value to each Hilbert space vector (valid only for all-sz)");
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbrspin-polarized", "number of orbitals which ar fully spin up polarized (from the one with the lowest momentum)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-pairing", "only fix Sz and not the number of particles");
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(3) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "quasiholes", "consider fermions with spin, where only the (k, r) partitions are admitted for each spin -- only available in use-pairing fermionic mode");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles within two Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "restrict-polarization", "restrict number of particles in each Landau level (provided by nbrparticles-up and nbrparticles-down)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-up", "number of particles in N=1 LL", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrparticles-down", "number of particles in N=0 LL", 0);  
  (*SystemGroup) += new BooleanOption  ('\n', "3-ll", "consider particles within three Landau levels");
  (*SystemGroup) += new BooleanOption  ('\n', "4-ll", "consider particles within four Landau levels");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('c', "complex-vector" , "vector state is complex instead of real");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*SystemGroup) += new SingleStringOption ('\n', "get-index", "find the index of a given n-body state");
  (*SystemGroup) += new MultipleIntegerOption ('\n', "pauli", "print only states obeying a general Pauli exclusion principle",',');
  (*OutputGroup) += new BooleanOption  ('\n', "variance", "show state variance");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_suN_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  int NbrFluxQuanta2 = Manager.GetInteger("nbrflux-2"); 
  int TotalLz = Manager.GetInteger("lz-value");
  int TotalTz = Manager.GetInteger("total-tz");
  int TotalY = Manager.GetInteger("total-y");
  int TotalJz = Manager.GetInteger("jz-value");
  int TotalKz = Manager.GetInteger("kz-value");
  int TzValue = Manager.GetInteger("tz-value");
  int YValue = Manager.GetInteger("y-value");
  bool AllLzFlag = Manager.GetBoolean("all-lz");
  bool SU2SpinFlag = Manager.GetBoolean("su2-spin");
  bool LzSymmetrizedBasis = Manager.GetBoolean("lzsymmetrized-basis");
  bool SzSymmetrizedBasis = Manager.GetBoolean("szsymmetrized-basis");
  bool AllSzFlag = Manager.GetBoolean("all-sz");
  bool AddIndex = Manager.GetBoolean("add-index");
  bool AddSzValue = Manager.GetBoolean("add-szvalue");
  bool SU3SpinFlag = Manager.GetBoolean("su3-spin");
  bool SU4SpinFlag = Manager.GetBoolean("su4-spin");
  bool TwoLLFlag = Manager.GetBoolean("2-ll");
  bool ThreeLLFlag = Manager.GetBoolean("3-ll");
  bool FourLLFlag = Manager.GetBoolean("4-ll");
  bool ComplexFlag = Manager.GetBoolean("complex-vector");
  bool FourDFlag = Manager.GetBoolean("4-D");
  bool CP2Flag = Manager.GetBoolean("cp2");
  bool S2xS2Flag = Manager.GetBoolean("s2xs2");
  bool TruncatedCP2Flag = Manager.GetBoolean("truncated-cp2");
  int MinY = Manager.GetInteger("min-y");
  bool SymFlagTzZ3 = Manager.GetBoolean("tzZ3symmetrized-basis");
  bool TzMinusParity = Manager.GetBoolean("minus-tzparity");
  int TotalSz = Manager.GetInteger("total-sz");
  int TotalIz = Manager.GetInteger("total-isosz");
  int TotalPz = Manager.GetInteger("total-entanglement");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  int PauliK=0, PauliR=0;
  int NbrOrbitals;
  if (FourDFlag == true)
    NbrOrbitals = (NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)*(NbrFluxQuanta + 3) / 6;
  if (CP2Flag == true)
    {
      NbrOrbitals = (NbrFluxQuanta + 1)*(NbrFluxQuanta + 2) / 2;
      if (TruncatedCP2Flag == false)      
	MinY = -2*NbrFluxQuanta;
      else
	NbrOrbitals -= (2*NbrFluxQuanta + MinY) * (2*NbrFluxQuanta + MinY + 3) / 18;
    }
  
  int TmpL=0;
  int *TmpIs=Manager.GetIntegers("pauli",TmpL);
  if (TmpL>0)
    {
      if (TmpL!=2)
	{
	  cout << "--pauli takes two arguments k,r describing the exclusion statistics"<<endl;
	  exit(1);
	}
      PauliK=TmpIs[0];
      PauliR=TmpIs[1];
      cout << "applying ("<<PauliK<<", "<<PauliR<<") exclusion statistics"<<endl;
    }
  
  if (FourDFlag == true)
    {
      if (Manager.GetBoolean("boson") == false)
	cout << "Warning : --4-D option only implemented in bosonic mode" << endl;
      
      if (- (NbrFluxQuanta*(NbrFluxQuanta+1)*(2*NbrFluxQuanta+1))/6 + NbrFluxQuanta*(NbrFluxQuanta*(NbrFluxQuanta+1))/2 + (NbrFluxQuanta + 1)*(NbrFluxQuanta + 1) + NbrParticles > 129)
	{
	  cout << "number of fermionic orbitals is too big to allow the storage of one state in one long integer" << endl;
	  return -1; 
	} 
      if (((NbrParticles * NbrFluxQuanta) & 1) != ((TotalJz + TotalKz) & 1)) 
	{
	  cout << "incompatible values for the number of particles, the number of flux quanta and twice the total jz and kz value (nbr-particles * nbr-flux and jz + kz should have the same parity)" << endl;
	  return -1;
	}
    }
  
  if (CP2Flag == true)
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  if ((NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)/2 + NbrParticles > 64)
	    {
	      cout  << "number of fermionic orbitals is too big to allow the storage of one state in one long integer" << endl;
	      return -1;
	    }
	}
      else
	if ((NbrFluxQuanta + 1)*(NbrFluxQuanta + 2)/2 > 128)
	  {
	    cout  << "number of fermionic orbitals is too big to allow the storage of one state in one long integer" << endl;
	    return -1;
	  }
      if ((YValue + 3*TzValue + 2*NbrParticles*NbrFluxQuanta) % 6 != 0)
	{
	  cout << "Y + 3Tz + 2N*Nphi should  be a multiple of 6" << endl;
	  return -1;
	}
      if ((YValue - 3*TzValue + 2*NbrParticles*NbrFluxQuanta) % 6 != 0)
	{
	  cout << "Y - 3Tz + 2N*Nphi should  be a multiple of 6" << endl;
	  return -1;
	}
    }
  
  if ((AllLzFlag == false) && (FourDFlag == false) && (CP2Flag == false))
    if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
      {
        cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
        return -1;
      }
  
      
  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      if (TwoLLFlag == false)
	{
	  if ((SU2SpinFlag == false) && (SU3SpinFlag == false))
	    {
	      if (HaldaneBasisFlag == false)
		{
		  if (FourDFlag == false)
		    {
		      if (CP2Flag == false)
			{
			  if (S2xS2Flag == false)
			    {
			      Space = new BosonOnSphereShort(NbrParticles, TotalLz, NbrFluxQuanta);
			    }
			  else
			    {
			      if (Manager.GetBoolean("s2s2-hardcorenonn") == true)
				{
				  if (((NbrFluxQuanta + 1 ) * (NbrFluxQuanta2 + 1)) < 64)
				    {
				      Space = new BosonOnS2xS2HardcoreNoNearestNeighbors(NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);
				    }
				  else
				    {
				      Space = new BosonOnS2xS2HardcoreNoNearestNeighborsLong(NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);

				    }
				}
			      else
				{
				  Space = new BosonOnS2xS2(NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);
				}
			    }
			}
		      else
			if (SymFlagTzZ3 == false)
			  Space = new BosonOnCP2(NbrParticles, NbrFluxQuanta, TzValue, YValue);
			else
			  Space = new BosonOnCP2TzZ3Symmetry(NbrParticles, NbrFluxQuanta, TzValue, YValue, TzMinusParity);
		    }
		  else
		  {
		    if (NbrOrbitals + NbrParticles < 65)
		      Space = new BosonOn4DSphere(NbrParticles, NbrFluxQuanta, TotalJz, TotalKz);
		    else 
		      Space = new BosonOn4DSphereLong(NbrParticles, NbrFluxQuanta, TotalJz, TotalKz);
		   }
		}
	      else
		{
		  int* ReferenceState = 0;
		  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
		    return -1;
		  Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
		}
	    }
	  else
	    {
	      if (SU2SpinFlag == true)
		{
		  if (Manager.GetInteger("nbrspin-polarized") > 0)
		    {
		      Space = new BosonOnSphereWithSU2SpinPartialPolarization(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, 
									      (int) Manager.GetInteger("nbrspin-polarized"));
		    }
		  else
		    {
		      if (Manager.GetBoolean("use-alt") == false)
			{		       
			  Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
			}
		      else
			{		       
			  if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
			    {
			      if (AllLzFlag == true)
				{
				  Space = new BosonOnSphereWithSU2SpinAllLz(NbrParticles, NbrFluxQuanta, TotalSz);
				}
			      else
				{
				  if (AllSzFlag == true)
				    {
				      Space = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, NbrFluxQuanta);
				    }
				  else
				    {
				      Space = new BosonOnSphereWithSU2Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
				    }
				}
			    }
			  else
			    {
			      if ((SzSymmetrizedBasis == true)  && (TotalSz == 0) && (LzSymmetrizedBasis == true) && (TotalLz == 0))
				{
				  if (AllSzFlag == true)
				    {
				      Space = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrParticles, NbrFluxQuanta, Manager.GetBoolean("minus-szparity"),
										       Manager.GetBoolean("minus-lzparity"));
				    }
				  else
				    {
				      Space = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"),
										       Manager.GetBoolean("minus-lzparity"));					  
				    }
				}
			      else 
				{
				  if ((SzSymmetrizedBasis == true)  && (TotalSz == 0))
				    {
				      if (AllSzFlag == true)
					{
					  Space = new BosonOnSphereWithSU2SpinSzSymmetry(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetBoolean("minus-szparity"));
					}
				      else
					{
					  Space = new BosonOnSphereWithSU2SpinSzSymmetry(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-szparity"));
					}
				    }
				  else
				    {
				      if (AllSzFlag == true)
					{
					  Space = new BosonOnSphereWithSU2SpinLzSymmetry(NbrParticles, NbrFluxQuanta, Manager.GetBoolean("minus-lzparity"));
					}
				      else
					{
					  Space = new BosonOnSphereWithSU2SpinLzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-lzparity"));
					}
				    }
				}
			    }
			}
		    }
		}
	      else
		{
		  if (SU3SpinFlag == true)
		    {
		      Space = new BosonOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);
		    }
		  else
		    {
		      if (SU4SpinFlag == true)
			Space = new BosonOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, TotalIz, TotalPz);
		    }
		}
	    }
	}
      else
        {
	  Space = new BosonOnSphereTwoLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta + 2, NbrFluxQuanta);    
        }

    }
  else
    {
      if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false) && (AllSzFlag == false) && (TwoLLFlag == false) && (ThreeLLFlag == false) && (FourLLFlag == false) && (CP2Flag == false) && (S2xS2Flag == false))
	{
	  if (HaldaneBasisFlag == false)
	    {
              if (AllLzFlag == false) 
                {
#ifdef __64_BITS__
	            if (NbrFluxQuanta <= 63)
		      Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	            else
		      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#else
  	            if (NbrFluxQuanta <= 31)
		      Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	            else
		      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#endif
                }
              else //Consider all Lz values 
                Space = new FermionOnSphereFull(NbrParticles, NbrFluxQuanta);
	    }
 	  else
 	    {
 	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
		return -1;
 	      Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
 	    }
	}
      else
 	if (SU2SpinFlag == true)
	  {
	    if (HaldaneBasisFlag == false)
	      {                 
		if (Manager.GetBoolean("use-pairing") == false)
		  {
		    if (TwoLLFlag == false)
		      {
			if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
			  {
			    if (Manager.GetInteger("nbrspin-polarized") > 0)
			      {
				Space = new FermionOnSphereWithSpinPartialPolarization(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, 
											(int) Manager.GetInteger("nbrspin-polarized"));
			      }
			    else
			      {
				Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
			      }
			  }
			else //either Lz or Sz symmetrized basis
			  {
			    if ((SzSymmetrizedBasis == true)  && (TotalSz == 0) && (LzSymmetrizedBasis == true) && (TotalLz == 0))
			      {
				Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrParticles, NbrFluxQuanta, Manager.GetBoolean("minus-szparity"),
										Manager.GetBoolean("minus-lzparity"));
			      }
			    else 
			      if ((SzSymmetrizedBasis == true)  && (TotalSz == 0))
				{
				  Space = new FermionOnSphereWithSpinSzSymmetry(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetBoolean("minus-szparity"));
				}
			      else
				Space = new FermionOnSphereWithSpinLzSymmetry(NbrParticles, NbrFluxQuanta, TotalSz, Manager.GetBoolean("minus-lzparity"));
			  }
		      }
		    else
		      {
			Space = new FermionOnSphereWithSpinTwoLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
		      }
		  }
		else
		  {
		    if (Manager.GetBoolean("quasiholes") == false)
		      {
			if ((SzSymmetrizedBasis == false) && (LzSymmetrizedBasis == false))
			  Space = new FermionOnSphereWithSpinAndPairing(TotalLz, NbrFluxQuanta, TotalSz);
		      }
		    else
		      {
			if (PauliK == 0 or PauliR == 0)
			  {
			    PauliK = 1;
			    PauliR = 2;
			  }
			Space = new QuasiholeOnSphereWithSpinAndPairing(PauliK, PauliR, TotalLz, NbrFluxQuanta, TotalSz, 0, "fermions");
		      }
		  }
	      }
	    else
	      {
		 int LzMax = NbrFluxQuanta;
		 int** ReferenceStates = 0;
		 int NbrReferenceStates;
		 bool TexturelessFlag; 
		 if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates, TexturelessFlag) == false)
		  {
		    cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
		    return 0;
		  }
		 if ( TexturelessFlag == false ) 
		  {
		      Space = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates);
		  }
		else
		  {
		    int **TexturelessReferenceState = new int*[NbrReferenceStates];
		    for ( int j = 0 ; j < NbrReferenceStates ; j++ ) 
		      {
			TexturelessReferenceState[j] = new int[LzMax+1];
			for ( int i = 0 ; i < (LzMax + 1) ; i++ )
			  {
			    if ( ReferenceStates[j][i] == 3 ) 
			      {
				  TexturelessReferenceState[j][i] = 2;
			      }
			    else if ( (ReferenceStates[j][i] == 1) || (ReferenceStates[j][i] == 2) ) 
			      {
				  TexturelessReferenceState[j][i] = 1;
			      }
			    else
			      {
				  TexturelessReferenceState[j][i] = 0;
			      }
			  }
		      }	
		    Space = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, TexturelessReferenceState, NbrReferenceStates, true);
		  }
	      }	    
	  }
	else 
	  if (AllSzFlag == true)
	    Space = new FermionOnSphereWithSpinAllSz(NbrParticles, TotalLz, NbrFluxQuanta);
	  else
	    if (SU3SpinFlag == true)
	      Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);
	    else
	      if (SU4SpinFlag == true)
		Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, TotalIz, TotalPz);	    
	      else
		if (TwoLLFlag == true)
                  {
                    if (Manager.GetBoolean("restrict-polarization")) 
                      Space = new FermionOnSphereTwoLandauLevels( Manager.GetInteger("nbrparticles-up"),  Manager.GetInteger("nbrparticles-down"), TotalLz, NbrFluxQuanta + 2, NbrFluxQuanta);	    
                    else
	  	      Space = new FermionOnSphereTwoLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta + 2, NbrFluxQuanta);	    
                  }
		else
		  if (ThreeLLFlag == true)
		    Space = new FermionOnSphereThreeLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta);	 
		  else
		    if (FourLLFlag == true)
		      Space = new FermionOnSphereFourLandauLevels(NbrParticles, TotalLz, NbrFluxQuanta);
		    else
		      {
			if (CP2Flag == true)
			  {
			    if (TruncatedCP2Flag == false)
			      {
				if (NbrOrbitals < 64)
				  Space = new FermionOnCP2(NbrParticles, NbrFluxQuanta, TzValue, YValue);
				else
				  Space = new FermionOnCP2Long(NbrParticles, NbrFluxQuanta, TzValue, YValue);
			      }
			    else
			      {
				Space = new FermionOnCP2(NbrParticles, NbrFluxQuanta, MinY, TzValue, YValue);
			      }
			  }
			if (S2xS2Flag == true)
			  {
			      if (Manager.GetBoolean("s2s2-hardcorenonn") == true)
				{
				  Space = new FermionOnS2xS2WithExclusionPrinciple(NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);
				}
			      else
				{
				  Space = new FermionOnS2xS2(NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);
				}
			  }
		      }
    }
  
  if (Manager.GetString("get-index") != 0)
    {
      long TmpIndex = Space->FindStateIndex(Manager.GetString("get-index"));
      if (TmpIndex == Space->GetHilbertSpaceDimension())
	{
	  cout << "state " << Manager.GetString("get-index") << " not found" << endl;
	}
      else
	{
	  cout << TmpIndex << " : ";
	  Space->PrintState(cout, TmpIndex) << endl;	   
	}
      return 0;
    }


  ofstream File;
  char* OutputFileName = 0;
  if (Manager.GetBoolean("save-disk") == true)
    {
      if (Manager.GetString("output-file") == 0)
	{
	  OutputFileName = new char[512];
	  if (Manager.GetBoolean("boson") == true)
	    if ((SU2SpinFlag == false) && (SU4SpinFlag == false) && (FourDFlag == false) && (CP2Flag == false) && (S2xS2Flag == false))
	      sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
	    else
	      if (SU2SpinFlag == true)
		sprintf (OutputFileName, "bosons_sphere_su2_n_%d_2s_%d_lz_%d_sz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz);
	      else
		if (SU3SpinFlag == true)
		  sprintf (OutputFileName, "bosons_sphere_su3_n_%d_2s_%d_lz_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalTz, TotalY);
		else
		  if (SU4SpinFlag == true)
		    sprintf (OutputFileName, "bosons_sphere_su4_n_%d_2s_%d_lz_%d_sz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz);
		  else
		    if (S2xS2Flag == true)
		      sprintf (OutputFileName, "bosons_s2xs2_n_%d_2s1_%d_2s2_%d_lz_%d_kz_%d.basis", NbrParticles, NbrFluxQuanta, NbrFluxQuanta2, TotalLz, TotalKz);
		    else
		      if (FourDFlag == true)
			sprintf (OutputFileName, "bosons_sphere4d_n_%d_2s_%d_jz_%d_kz_%d.basis", NbrParticles, NbrFluxQuanta, TotalJz, TotalKz);
		      else
			{
			  if (SymFlagTzZ3 == false)
			    sprintf (OutputFileName, "bosons_sphereCP2_n_%d_2s_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TzValue, YValue);
			  else
			    sprintf (OutputFileName, "bosons_sphereCP2_tzZ3sym_n_%d_2s_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TzValue, YValue);
			}
	  else
	    if (SU2SpinFlag == true)
	      sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d_lz_%d_sz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz);
	    else
	      if (AllSzFlag == true)
		sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d_lz_%d_allsz.basis", NbrParticles, NbrFluxQuanta, TotalLz);
	      else
		if (SU3SpinFlag == true)
		  sprintf (OutputFileName, "fermions_sphere_su3_n_%d_2s_%d_lz_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalTz, TotalY);
		else
		  if (SU4SpinFlag == true)
		    sprintf (OutputFileName, "fermions_sphere_su4_n_%d_2s_%d_lz_%d_sz_%d_iz_%d_pz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz, TotalIz, TotalPz);
		  else
		    if (TwoLLFlag == true)
		      sprintf (OutputFileName, "fermions_sphere_2ll_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
		    else
		      if (ThreeLLFlag == true)
			sprintf (OutputFileName, "fermions_sphere_3ll_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
		      else
			if (FourLLFlag == true)
			  sprintf (OutputFileName, "fermions_sphere_4ll_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
		      else
			if (CP2Flag == true)
			  sprintf(OutputFileName, "fermions_spherecp2_n_%d_2s_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TzValue, YValue);
			else
			  sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
	  if (PauliK != 0)
	    {
	      char* TmpNewExtension = new char[32];
	      sprintf (TmpNewExtension, "pauli_%d_%d.basis", PauliK, PauliR);
	      char* OutputFileName2 = ReplaceExtensionToFileName(OutputFileName, "basis", TmpNewExtension);
	      delete[] OutputFileName;
	      delete[] TmpNewExtension;
	      OutputFileName = OutputFileName2;
	    }
	}
      else
	{
	  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFileName, Manager.GetString("output-file"));
	}
      File.open(OutputFileName, ios::binary | ios::out);
      if (!File.is_open())
	{
	  cout << "Cannot create file " << OutputFileName << endl;
	  return -1;
	}
      File.precision(14);
    }
  else
    cout.precision(14);
  if (Manager.GetString("state") == 0)
    {
      if (Manager.GetBoolean("save-disk") == true)
	{
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      if ((PauliK == 0) || (Space->HasPauliExclusions(i, PauliK, PauliR)))
		{
		  if (AddIndex == true) 
		    File << i << " ";
		  Space->PrintState(File, i);
		  if (AddSzValue == true) File<<" Sz= "<<Space->GetSzValue(i)<< endl;
		else File<<endl;
		}
	    }
	}
      else
	{
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      if ((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR)) || (Manager.GetBoolean("quasiholes")))
		{
		  if (AddIndex == true) 
		    cout << i <<" ";
		  Space->PrintState(cout, i);
		  if (AddSzValue == true) cout<<" Sz= "<<Space->GetSzValue(i)<< endl;
		  else cout<<endl;
		}
	    }
	}
    }
  else
   if (ComplexFlag == false)
    {
      if (Manager.GetBoolean("rational") == false)
	{
	  int NbrHiddenComponents = 0;
	  double WeightHiddenComponents = 0.0;
	  double Normalization = 0.0;
	  RealVector State;

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
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if ((fabs(State[i]) > Error)&&((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
			{
			  if (SymFlagTzZ3 == false)
			  {
			    if (AddIndex == true) 
			      File << i << " ";	
			    Space->PrintState(File, i) << " : " << State[i];
			    if (AddSzValue == true) File<<" Sz= "<<Space->GetSzValue(i)<< endl;
			    else File<<endl;
			  }
			  else
			  {
			    Space->PrintState(File, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    File << endl;
			  }
			}
		      else		     
			{
			  WeightHiddenComponents += State[i] * State[i];
			  ++NbrHiddenComponents; 
			}
		      Normalization += State[i] * State[i];
		    }
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  {
		    if ((fabs(State[i]) > Error)&&((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
		      {
			if (SymFlagTzZ3 == false)
			{
			  if (AddIndex == true) 
			    cout << i <<" ";
			  Space->PrintState(cout, i) << " : " << State[i];
			  if (AddSzValue == true) cout<<" Sz= "<<Space->GetSzValue(i)<< endl;
			  else cout<<endl;
			}
			else
			  {
			    Space->PrintState(cout, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    cout << endl;
			  }
		      }
		    else		     
		      {
			WeightHiddenComponents += State[i] * State[i];
			++NbrHiddenComponents; 
		      }
		    Normalization += State[i] * State[i];
		  }
	      cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
	    }
	  else
	    {             
	      if (Manager.GetBoolean("variance"))
		{
		  if (Manager.GetBoolean("save-disk") == true)
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(File, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(File, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;
		    }
		  else
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(cout, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(cout, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("save-disk") == true)
		    {
		      if (PauliK==0)
		      {
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			{
			if (SymFlagTzZ3 == false)
			    Space->PrintState(File, i) << " : " << State[i] << endl;
			else
			  {
			    Space->PrintState(File, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    File << endl;
			  }
			}
		      }
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(File, i) << " : " << State[i] << endl;
		    }
		  else
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			{
			if (SymFlagTzZ3 == false)
			    Space->PrintState(cout, i) << " : " << State[i] << endl;
			else
			  {
			    Space->PrintState(cout, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    cout << endl;
			  }
			}
		      else
			{
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if ((Space->HasPauliExclusions(i,PauliK,PauliR)) || Manager.GetBoolean("quasiholes"))
			      Space->PrintState(cout, i) << " : " << State[i] << endl;	   
			}
		    }
		}
	    }
	}
      else
	{
	  int NbrHiddenComponents = 0;
	  double WeightHiddenComponents = 0.0;
	  double Normalization = 0.0;
	  LongRationalVector State;
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
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if ((State[i] != 0l) && ((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
			{
			  if (SymFlagTzZ3 == false)
			  {
			    if (AddIndex == true) 
			      File << i << " ";	
			    Space->PrintState(File, i) << " : " << State[i];
			    if (AddSzValue == true) 
			      File <<" Sz= "<<Space->GetSzValue(i)<< endl;
			    else 
			      File << endl;
			  }
			  else
			  {
			    Space->PrintState(File, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    File << endl;
			  }
			}
		      else		     
			{
			  ++NbrHiddenComponents; 
			}
		    }
		}
	      else
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if ((State[i] != 0l) && ((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
			{
			  if (SymFlagTzZ3 == false)
			  {
			    if (AddIndex == true) 
			      cout << i <<" ";
			    Space->PrintState(cout, i) << " : " << State[i];
			    if (AddSzValue == true)
			      cout<< " Sz= "<< Space->GetSzValue(i)<< endl;
			    else 
			      cout << endl;
			  }
			  else
			  {
			    Space->PrintState(cout, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    cout << endl;
			  }
			}
		      else		     
			{
			  ++NbrHiddenComponents; 
			}
		    }
		}
	      cout << NbrHiddenComponents << " hidden components" << endl;
	    }
	  else
	    {
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  if (PauliK == 0)
		    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		      if (SymFlagTzZ3 == false)
			Space->PrintState(File, i) << " : " << State[i] << endl;
		      else
			  {
			    Space->PrintState(File, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    File << endl;
			  }
		  else
		    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		      if (Space->HasPauliExclusions(i,PauliK,PauliR))
			Space->PrintState(File, i) << " : " << State[i] << endl;
		}
	      else
		{
		  if (PauliK == 0)
		    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		      if (SymFlagTzZ3 == false)
			Space->PrintState(cout, i) << " : " << State[i] << endl;
		      else
			  {
			    Space->PrintState(cout, i) << " : " << State[i] <<  " : " << Space->GetSymmetryDimension(i);
			    cout << endl;
			  }
		  else
		    {
		      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			if (Space->HasPauliExclusions(i,PauliK,PauliR))
			  Space->PrintState(cout, i) << " : " << State[i] << endl;	   
		    }
		}
	    }
	}
    }
  else //complex vector
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
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if ((Norm(State[i]) > Error)&&((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
			{
			  if (AddIndex == true) 
			    File << i << " ";	
			  Space->PrintState(File, i) << " : " << State.Re(i)<<" "<<State.Im(i);
			  if (AddSzValue == true) File<<" Sz= "<<Space->GetSzValue(i)<< endl;
			  else File<<endl;
			}
		      else		     
			{
			  WeightHiddenComponents += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
			  ++NbrHiddenComponents; 
			}
		      Normalization += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
		    }
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  {
		    if ((Norm(State[i]) > Error)&&((PauliK==0)||(Space->HasPauliExclusions(i,PauliK,PauliR))))
		      {
			if (AddIndex == true) 
			  cout << i <<" ";
			Space->PrintState(cout, i) << " : " << State.Re(i)<<" "<<State.Im(i);
			if (AddSzValue == true) cout<<" Sz= "<<Space->GetSzValue(i)<< endl;
			else cout<<endl;
		      }
		    else		     
		      {
			WeightHiddenComponents += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
			++NbrHiddenComponents; 
		      }
		    Normalization += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
		  }
	      cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
	    }
	  else
	    {             
	      if (Manager.GetBoolean("variance"))
		{
		  if (Manager.GetBoolean("save-disk") == true)
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(File, i) << " : " << Space->StateVariance(i) << " " << State.Re(i) << " " << State.Im(i) << " " << i <<endl;
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(File, i) << " : " << Space->StateVariance(i) << " " << State.Re(i) << " " << State.Im(i) << " " << i <<endl;
		    }
		  else
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(cout, i) << " : " << Space->StateVariance(i) << " " << State.Re(i) << " " << State.Im(i) << " " << i <<endl;
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(cout, i) << " : " << Space->StateVariance(i) << " " << State.Re(i) << " " << State.Im(i) << " " << i <<endl;
		    }
		}
	      else
		{
		  if (Manager.GetBoolean("save-disk") == true)
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(File, i) << " : " << State.Re(i) << " " << State.Im(i) << endl;
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  if (Space->HasPauliExclusions(i,PauliK,PauliR))
			    Space->PrintState(File, i) << " : " << State.Re(i) << " " << State.Im(i) << endl;
		    }
		  else
		    {
		      if (PauliK==0)
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(cout, i) << " : " << State.Re(i) << " " << State.Im(i) << endl;
		      else
			{
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (Space->HasPauliExclusions(i,PauliK,PauliR))
			      Space->PrintState(cout, i) << " : " << State.Re(i) << " " << State.Im(i) << endl;	   
			}
		    }
		}
	    }
    }

 
  if (OutputFileName != 0)
    {
      File.close();
      delete[] OutputFileName;
    }
}

