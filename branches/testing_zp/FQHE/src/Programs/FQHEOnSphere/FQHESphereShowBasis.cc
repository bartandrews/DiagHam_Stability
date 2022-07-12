#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
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
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2)/su(4) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-sz", "consider particles with SU(2) spin all Sz components");
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "add index of the Hilbert space vectors");
  (*SystemGroup) += new BooleanOption  ('\n', "add-szvalue", "add Sz value to each Hilbert space vector (valid only for all-sz)");
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*OutputGroup) += new BooleanOption  ('\n', "variance", "show state variance");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_suN_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereShowBasis -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  int TotalTz = ((SingleIntegerOption*) Manager["total-tz"])->GetInteger();
  int TotalY = ((SingleIntegerOption*) Manager["total-y"])->GetInteger();
  bool SU2SpinFlag = ((BooleanOption*) Manager["su2-spin"])->GetBoolean();
  bool AllSzFlag = ((BooleanOption*) Manager["all-sz"])->GetBoolean();
  bool AddIndex = ((BooleanOption*) Manager["add-index"])->GetBoolean();
  bool AddSzValue = ((BooleanOption*) Manager["add-szvalue"])->GetBoolean();
  bool SU3SpinFlag = ((BooleanOption*) Manager["su3-spin"])->GetBoolean();
  bool SU4SpinFlag = ((BooleanOption*) Manager["su4-spin"])->GetBoolean();
  int TotalSz = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  int TotalIz = ((SingleIntegerOption*) Manager["total-isosz"])->GetInteger();
  int TotalPz = ((SingleIntegerOption*) Manager["total-entanglement"])->GetInteger();
  bool BernevigFlag = false;
    
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
    {
      cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
      return -1;
    }

  ParticleOnSphere* Space = 0;
  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      if (SU2SpinFlag == false)
	{
	  if (BernevigFlag == false)
	    {
	      Space = new BosonOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	    }
	  else
	    {
	      int* ReferenceState = 0;
	      ReferenceState = new int[NbrFluxQuanta + 1];
	      for (int i = 0; i <= NbrFluxQuanta; ++i)
		ReferenceState[i] = 0;
	      for (int i = 0; i <= NbrFluxQuanta; i += 2)
		ReferenceState[i] = 1;
	      //	  Space = new BosonOnSphereBernevigBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
	      return 0;
	    }
	}
      else
	{
	  Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	}
    }
  else
    {
      if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false) && (AllSzFlag == false))
	{
	  if (BernevigFlag == false)
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
 	  else
 	    {
 	      int* ReferenceState = 0;
 	      ReferenceState = new int[NbrFluxQuanta + 1];
 	      for (int i = 0; i <= NbrFluxQuanta; ++i)
 		ReferenceState[i] = 0;
 	      for (int i = 0; i <= NbrFluxQuanta; i += 3)
 		ReferenceState[i] = 1;
 	      Space = 0;
	      //new FermionOnSphereBernevigBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, 3);
 	    }
	}
      else
 	if (SU2SpinFlag == true)
	 Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	else 
	  if (AllSzFlag == true)
	    Space = new FermionOnSphereWithSpinAllSz(NbrParticles, TotalLz, NbrFluxQuanta);
	  else
	   if (SU3SpinFlag == true)
	     Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);
	   else
	    if (SU4SpinFlag == true)
	      Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, TotalIz, TotalPz);	    
    }
  

  ofstream File;
  char* OutputFileName = 0;
  if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
    {
      if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
	{
	  OutputFileName = new char[512];
	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	    if ((SU2SpinFlag == false) && (SU4SpinFlag == false))
	      sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
	    else
	      if (SU2SpinFlag == true)
		sprintf (OutputFileName, "bosons_sphere_su2_n_%d_2s_%d_lz_%d_sz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz);
	      else
		if (SU3SpinFlag == true)
		  sprintf (OutputFileName, "bosons_sphere_su3_n_%d_2s_%d_lz_%d_tz_%d_y_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalTz, TotalY);
		else
		  sprintf (OutputFileName, "bosons_sphere_su4_n_%d_2s_%d_lz_%d_sz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz, TotalSz);
	  else
	    if ((SU2SpinFlag == false) && (SU4SpinFlag == false))
	      sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d_lz_%d.basis", NbrParticles, NbrFluxQuanta, TotalLz);
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
	}
      else
	{
	  OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
	  strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
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
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	{
	  if (AddIndex == true) File<<i<<" ";
	  Space->PrintState(File, i);
	  if (AddSzValue == true) File<<" Sz= "<<Space->GetSzValue(i)<< endl;
	   else File<<endl;
        }
    else
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	{
	  if (AddIndex == true) cout<<i<<" ";
	  Space->PrintState(cout, i);
	  if (AddSzValue == true) cout<<" Sz= "<<Space->GetSzValue(i)<< endl;
           else cout<<endl;
        }
      
   else
     {
       int NbrHiddenComponents = 0;
       double WeightHiddenComponents = 0.0;
       double Normalization = 0.0;
       RealVector State;
       if (State.ReadVector(((SingleStringOption*) Manager["state"])->GetString()) == false)
	 {
	   cout << "error while reading " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
	   return -1;
	 }
       if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
	 {
	   cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	   return -1;
	 }
       if (((SingleDoubleOption*) Manager["hide-component"])->GetDouble() > 0.0)
	 {
	   double Error = ((SingleDoubleOption*) Manager["hide-component"])->GetDouble();
	   if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
	     {
	       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		 {
		   if (fabs(State[i]) > Error)
                     {
		       if (AddIndex == true) File<<i<<" ";	
		       Space->PrintState(File, i) << " : " << State[i];
		       if (AddSzValue == true) File<<" Sz= "<<Space->GetSzValue(i)<< endl;
                        else File<<endl;
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
		 if (fabs(State[i]) > Error)
		   {
                     if (AddIndex == true) cout<<i<<" ";
		     Space->PrintState(cout, i) << " : " << State[i];
	             if (AddSzValue == true) cout<<" Sz= "<<Space->GetSzValue(i)<< endl;
                      else cout<<endl;
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
	       if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   Space->PrintState(File, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;	   
	       else
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   Space->PrintState(cout, i) << " : " << Space->StateVariance(i) << " " << State[i] << " " << i <<endl;	   
	     }
	   else
	     {
	       if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   Space->PrintState(File, i) << " : " << State[i] << endl;	   
	       else
		 for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		   Space->PrintState(cout, i) << " : " << State[i] << endl;	   
	     }
	 }
     }

  if (OutputFileName != 0)
    {
      File.close();
      delete[] OutputFileName;
    }
}

