////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of manager for particles on sphere                 //
//                                                                            //
//                        last modification : 05/10/2008                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "HilbertSpace/ParticleOnSphereManager.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereUnnormalizedBasis.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSqueezedBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"

#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinZ3Symmetry.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzZ3Symmetry.h"

#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include "GeneralTools/ConfigurationParser.h"


#include <stdlib.h>


// default constructor
//
// fermionFlag = use fermionic statistics if true
// bosonFlag = use bosonic statistics if true, if both are fermionFlag and bosonFlag are equal, add an option to handle statistic selection
// sUKIndex = indicate which SU(K) internal degree of freedom as to be used
// filenameFlag = true if some option values can/will be retrieved from a file name

ParticleOnSphereManager::ParticleOnSphereManager(bool fermionFlag, bool bosonFlag, int sUKIndex, bool filenameFlag)
{
  this->FermionFlag = fermionFlag;
  this->BosonFlag = bosonFlag;
  if ((this->FermionFlag == false) && (this->BosonFlag == false))
    {
      this->FermionFlag = true;
      this->BosonFlag = true;
    }
  this->SUKIndex = sUKIndex;
  this->FilenameFlag = filenameFlag;
  this->Options = 0;
}

// destructor
//

ParticleOnSphereManager::~ParticleOnSphereManager()
{
}
 
// add an option group containing all options related to the Hilbert space construction 
//
// manager = pointer to the option manager

void ParticleOnSphereManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* SystemGroup  = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  (*(this->Options)) += SystemGroup;
  (*(this->Options)) += PrecalculationGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  if ((this->FermionFlag == true) && (this->BosonFlag == true))
    (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  switch (this->SUKIndex)
    {
    case 1:
      {
	(*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)"); 
	(*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
	if (this->BosonFlag == false)
	  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
	(*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
	if (this->BosonFlag == false)
	  (*SystemGroup) += new BooleanOption  ('\n', "unnormalized-basis", "do not normalized Fock states"); 
	if (this->FermionFlag == true)
	  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);	
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
      }
      break;
    case 2:
      {
	(*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
	(*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
	(*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
      }
      break;
    case 3:
      {
	(*SystemGroup) += new SingleIntegerOption  ('t', "total-tz", "twice the quantum number of the system associated to the Tz generator", 0);
	(*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "three time the quantum number of the system associated to the Y generator", 0);
	(*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "tzsymmetrized-basis", "use Tz <-> -Tz symmetrized version of the basis (only valid if total-tz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "z3symmetrized-basis", "use Z3 symmetrized version of the basis (only valid if total-y=0 and total-tz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-tzparity", "select the  Tz <-> -Tz symmetric sector with negative parity");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
	(*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
      }
      break;
    case 4:
      {
	(*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
	(*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin (i.e valley SU(2) degeneracy) of the system", 0);
	(*SystemGroup) += new BooleanOption  ('\n', "use-entanglement", "use a define value for the spin-isopsin entanglement of the system");
	(*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system", 0);
	(*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
      }
      break;
    }
}

// get the Hilbert space defined by the running options
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpace(int totalLz)
{
  if ((this->FermionFlag == true) && (this->BosonFlag == true) && (((SingleStringOption*) (*(this->Options))["statistics"])->GetString() != 0))
    {
      if ((strcmp ("fermions", ((SingleStringOption*) (*(this->Options))["statistics"])->GetString()) == 0))
	{
	  this->FermionFlag = true;
	  this->BosonFlag = false;
	}
      else
	if ((strcmp ("bosons", ((SingleStringOption*) (*(this->Options))["statistics"])->GetString()) == 0))
	  {
	    this->FermionFlag = false;
	    this->BosonFlag = true;
	  }
	else
	  {
	    cout << ((SingleStringOption*) (*(this->Options))["statistics"])->GetString() << " is an undefined statistics" << endl;
	    return 0;
	}  
    }
  switch (this->SUKIndex)
    {
    case 1:
      return this->GetHilbertSpaceU1(totalLz);
    case 2:
      return this->GetHilbertSpaceSU2(totalLz);
    case 3:
      return this->GetHilbertSpaceSU3(totalLz);
    case 4:
      return this->GetHilbertSpaceSU4(totalLz);
    }
  return 0;
}

// get the Hilbert space defined by the running options and a given Total Lz value and for the U(1) case
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpaceU1(int totalLz)
{
  if (this->BosonFlag == false)
    {
      int NbrParticles = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      bool HaldaneBasisFlag = ((BooleanOption*) (*(this->Options))["haldane"])->GetBoolean();
      bool SymmetrizedBasis = ((BooleanOption*) (*(this->Options))["symmetrized-basis"])->GetBoolean();
      bool UnnormalizedBasis = ((BooleanOption*) (*(this->Options))["unnormalized-basis"])->GetBoolean();
      unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) (*(this->Options))["fast-search"])->GetInteger()) << 20;
      ParticleOnSphere* Space = 0;
      if (UnnormalizedBasis == true)
	{
	  Space = new  FermionOnSphereUnnormalizedBasis(NbrParticles, totalLz, LzMax, MemorySpace);
	  return Space;
	}
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	  if (LzMax <= 30)
#endif
	    if ((SymmetrizedBasis == false) || (totalLz != 0))
	      Space = new FermionOnSphere(NbrParticles, totalLz, LzMax, MemorySpace);
	    else
	      {
		if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
		  Space = new FermionOnSphereSymmetricBasis(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		else
		  Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
		if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
		  {
		    ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
		    return 0;
		  }
	      }
	  else
#ifdef __128_BIT_LONGLONG__
	    if (LzMax <= 126)
#else
	      if (LzMax <= 62)
#endif
		{
		  if ((SymmetrizedBasis == false) || (totalLz != 0))
		    Space = new FermionOnSphereLong(NbrParticles, totalLz, LzMax, MemorySpace);
		  else
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
			Space = new FermionOnSphereSymmetricBasisLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		      else
			Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		      if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
			{
			  ((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			  return 0;
			}
		    }
		}
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, totalLz, LzMax, MemorySpace);
	}
      else
	{
	  int* ReferenceState = 0;
	  if (((SingleStringOption*) (*(this->Options))["reference-file"])->GetString() == 0)
	    {
	      ReferenceState = new int[LzMax + 1];
	      for (int i = 0; i <= LzMax; ++i)
		ReferenceState[i] = 0;
	      if (strcasecmp(((SingleStringOption*) (*(this->Options))["reference-state"])->GetString(), "laughlin") == 0)
		for (int i = 0; i <= LzMax; i += 3)
		  ReferenceState[i] = 1;
	      else
		if (strcasecmp(((SingleStringOption*) (*(this->Options))["reference-state"])->GetString(), "pfaffian") == 0)
		  for (int i = 0; i <= LzMax; i += 4)
		    {
		      ReferenceState[i] = 1;
		      ReferenceState[i + 1] = 1;
		    }
		else
		  if (strcasecmp(((SingleStringOption*) (*(this->Options))["reference-state"])->GetString(), "readrezayi3") == 0)
		    for (int i = 0; i <= LzMax; i += 5)
		      {
			ReferenceState[i] = 1;
			ReferenceState[i + 1] = 1;
			ReferenceState[i + 2] = 1;
		      }
		  else
		    {
		      cout << "unknown reference state " << ((SingleStringOption*) (*(this->Options))["reference-state"])->GetString() << endl;
		      return 0;
		    }
	    }
	  else
	    {
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(((SingleStringOption*) (*(this->Options))["reference-file"])->GetString()) == false)
		{
		  ReferenceStateDefinition.DumpErrors(cout) << endl;
		  return 0;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		{
		  cout << "NbrParticles is not defined or as a wrong value" << endl;
		  return 0;
		}
	      if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
		{
		  cout << "LzMax is not defined or as a wrong value" << endl;
		  return 0;
		}
	      int MaxNbrLz;
	      if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		{
		  cout << "error while parsing ReferenceState in " << ((SingleStringOption*) (*(this->Options))["reference-file"])->GetString() << endl;
		  return 0;     
		}
	      if (MaxNbrLz != (LzMax + 1))
		{
		  cout << "wrong LzMax value in ReferenceState" << endl;
		  return 0;     
		}
	    }
	  if (SymmetrizedBasis == false)
	     {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
		   {
		     if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
		       Space = new FermionOnSphereHaldaneBasis(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneBasis(NbrParticles, totalLz, LzMax, ReferenceState, MemorySpace);
		     if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
		       {
			 ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			 return 0;
		       }
		   }
	       else
#ifdef __128_BIT_LONGLONG__
		 if (LzMax <= 126)
#else
		   if (LzMax <= 62)
#endif
		     {
		       if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
			 Space = new FermionOnSphereHaldaneBasisLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		       else
			 Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, totalLz, LzMax, ReferenceState, MemorySpace);
		       if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
			 {
			   ((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			   return 0;
			 }
		     }	       
	     }
	  else
	    {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
		   {
		     if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
		       Space = new FermionOnSphereHaldaneSymmetricBasis(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		     if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
		       {
			 ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			 return 0;
		       }
		   }
		 else
#ifdef __128_BIT_LONGLONG__
		   if (LzMax <= 126)
#else
		     if (LzMax <= 62)
#endif
		       {
			 if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() != 0)
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
			 else
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
			 if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
			   {
			     ((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			     return 0;
			   }
		       }
	    }
	}
      return Space;
    }
  else
    {
      int NbrBosons = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      bool HaldaneBasisFlag = ((BooleanOption*) (*(this->Options))["haldane"])->GetBoolean();
      ParticleOnSphere* Space = 0;
      if (HaldaneBasisFlag == false)
	{
	  bool SymmetrizedBasis = ((BooleanOption*) (*(this->Options))["symmetrized-basis"])->GetBoolean();
#ifdef  __64_BITS__
	  if ((LzMax + NbrBosons - 1) < 63)
#else
	    if ((LzMax + NbrBosons - 1) < 31)	
#endif
	      {
		if ((SymmetrizedBasis == false) || (totalLz != 0))
		  Space = new BosonOnSphereShort(NbrBosons, totalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasisShort(NbrBosons, LzMax);
	      }
	    else
	      {
		if ((SymmetrizedBasis == false) || (totalLz != 0))
		  Space = new BosonOnSphere (NbrBosons, totalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasis(NbrBosons, LzMax);
	      }
	}
      else 
	{
	  int* ReferenceState = 0;
	  if (((SingleStringOption*) (*(this->Options))["reference-file"])->GetString() == 0)
	    {
	      cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
	      return 0;
	    }
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(((SingleStringOption*) (*(this->Options))["reference-file"])->GetString()) == false)
	    {
	      ReferenceStateDefinition.DumpErrors(cout) << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrBosons) == false) || (NbrBosons <= 0))
	    {
	      cout << "NbrParticles is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
	    {
	      cout << "LzMax is not defined or as a wrong value" << endl;
	      return 0;
	    }
	  int MaxNbrLz;
	  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
	    {
	      cout << "error while parsing ReferenceState in " << ((SingleStringOption*) (*(this->Options))["reference-file"])->GetString() << endl;
	      return 0;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return 0;     
	    }
#ifdef  __64_BITS__
	  if ((LzMax + NbrBosons - 1) < 63)
#else
	    if ((LzMax + NbrBosons - 1) < 31)	
#endif
	      Space = new BosonOnSphereHaldaneBasisShort(NbrBosons, totalLz, LzMax, ReferenceState);	  
	}
      return Space;
    }
}

// get the Hilbert space defined by the running options and a given Total Lz value and for the SU(2) case
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpaceSU2(int totalLz)
{
  if (this->BosonFlag == false)
    {
      ParticleOnSphereWithSpin* Space = 0;
      int NbrFermions = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      int SzTotal = ((SingleIntegerOption*) (*(this->Options))["total-sz"])->GetInteger();
      bool LzSymmetrizedBasis = ((BooleanOption*) (*(this->Options))["lzsymmetrized-basis"])->GetBoolean();
      bool SzSymmetrizedBasis = ((BooleanOption*) (*(this->Options))["szsymmetrized-basis"])->GetBoolean();
      unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) (*(this->Options))["fast-search"])->GetInteger()) << 20;
      if (((SzSymmetrizedBasis == false) || (SzTotal != 0)) && ((LzSymmetrizedBasis == false) || (totalLz != 0)))
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	    if (LzMax <= 15)
#endif
	      {
		Space = new FermionOnSphereWithSpin(NbrFermions, totalLz, LzMax, SzTotal, MemorySpace);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (LzMax <= 63)
#else
		  if (LzMax <= 31)
#endif
		    {
		      Space = new FermionOnSphereWithSpinLong(NbrFermions, totalLz, LzMax, SzTotal, MemorySpace);
		    }
		  else
		    {
		      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		      return 0;
		    }	
	      }
	}
      else
	{
#ifdef __128_BIT_LONGLONG__
	  if (LzMax >= 61)
#else
	    if (LzMax >= 29)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return 0;
	      }	
	  if ((SzSymmetrizedBasis == true)  && (SzTotal == 0))
	    if ((LzSymmetrizedBasis == false) || (totalLz != 0))
	      {
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinSzSymmetry(NbrFermions, totalLz, LzMax, ((BooleanOption*) (*(this->Options))["minus-szparity"])->GetBoolean(), MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		    }
		  else
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrFermions, totalLz, LzMax, ((BooleanOption*) (*(this->Options))["minus-szparity"])->GetBoolean(), MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinSzSymmetryLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		    }
	      }
	    else
#ifdef __64_BITS__
	      if (LzMax <= 28)
#else
		if (LzMax <= 13)
#endif
		  {
		    if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrFermions, LzMax, ((BooleanOption*) (*(this->Options))["minus-szparity"])->GetBoolean(),
									((BooleanOption*) (*(this->Options))["minus-lzparity"])->GetBoolean(), MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		  }
		else
		  {
		    if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
		      {
			Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrFermions, LzMax, ((BooleanOption*) (*(this->Options))["minus-szparity"])->GetBoolean(),
									    ((BooleanOption*) (*(this->Options))["minus-lzparity"])->GetBoolean(), MemorySpace);
		      }
		    else
		      Space = new FermionOnSphereWithSpinLzSzSymmetryLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		    
		  }
	      else
#ifdef __64_BITS__
		if (LzMax <= 28)
#else
		  if (LzMax <= 13)
#endif
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinLzSymmetry(NbrFermions, LzMax, SzTotal, ((BooleanOption*) (*(this->Options))["minus-lzparity"])->GetBoolean(), MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);	      
		    }
		  else
		    {
		      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
			Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrFermions, LzMax, SzTotal, ((BooleanOption*) (*(this->Options))["minus-lzparity"])->GetBoolean(), MemorySpace);
		      else
			Space = new FermionOnSphereWithSpinLzSymmetryLong(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);	      
		    }
	  if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
	    {
	      ((FermionOnSphereWithSpinLzSzSymmetry*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
	      return 0;
	    }
	}
      return Space;
    }
  else
    {
      ParticleOnSphereWithSpin* Space = 0;
      int NbrBosons = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      int SzTotal = ((SingleIntegerOption*) (*(this->Options))["total-sz"])->GetInteger();
      Space = new BosonOnSphereWithSpin(NbrBosons, totalLz, LzMax, SzTotal);
      return Space;
    }
}

// get the Hilbert space defined by the running options and a given Total Lz value and for the SU(3) case
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpaceSU3(int totalLz)
{
  if (this->BosonFlag == false)
    {
      int NbrFermions = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      int TotalTz = ((SingleIntegerOption*) (*(this->Options))["total-tz"])->GetInteger();
      int TotalY = ((SingleIntegerOption*) (*(this->Options))["total-y"])->GetInteger();
//      bool LzSymmetrizedBasis = ((BooleanOption*) (*(this->Options))["lzsymmetrized-basis"])->GetBoolean();
      bool TzSymmetrizedBasis = ((BooleanOption*) (*(this->Options))["tzsymmetrized-basis"])->GetBoolean();
      bool Z3SymmetrizedBasis = ((BooleanOption*) (*(this->Options))["z3symmetrized-basis"])->GetBoolean();
      unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) (*(this->Options))["fast-search"])->GetInteger()) << 20;      
      ParticleOnSphereWithSU3Spin* Space = 0;
      if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == false))
	{
#ifdef __64_BITS__
	  if (LzMax <= 20)
#else
	    if (LzMax <= 9)
#endif
	      {
		Space = new FermionOnSphereWithSU3Spin(NbrFermions, totalLz, LzMax, TotalTz, TotalY, MemorySpace);
	      }
	    else
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return 0;
	      }	
	}
      else
	{
#ifdef __64_BITS__
	  if (LzMax > 20)
#else
	    if (LzMax > 9)
#endif
	      {
		cout << "States of this Hilbert space cannot be represented in a single word." << endl;
		return 0;
	      }	
	  if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == false))
	    {
	      if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(NbrFermions, totalLz, LzMax, TotalY, ((BooleanOption*) (*(this->Options))["minus-tzparity"])->GetBoolean(), MemorySpace);
		  if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
		    {
		      ((FermionOnSphereWithSU3SpinTzSymmetry*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
		      return 0;
		    }
		}
	      else
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		}
	    }
	  else
	    if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == true))
	      {
		if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(NbrFermions, totalLz, LzMax, TotalTz, MemorySpace);
		    if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
		      {
			((FermionOnSphereWithSU3SpinZ3Symmetry*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			return 0;
		      }
		  }
		else
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		  }
	      }
	    else
	      if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == true))
		{
		  if (((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString() == 0)
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(NbrFermions, totalLz, LzMax, ((BooleanOption*) (*(this->Options))["minus-tzparity"])->GetBoolean(), MemorySpace);
		      if (((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString() != 0)
			{
			  ((FermionOnSphereWithSU3SpinTzZ3Symmetry*) Space)->WriteHilbertSpace(((SingleStringOption*) (*(this->Options))["save-hilbert"])->GetString());
			  return 0;
			}
		    }
		  else
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(((SingleStringOption*) (*(this->Options))["load-hilbert"])->GetString(), MemorySpace);
		    }
		}
	}
      return Space;
    }
  else
    {
      return 0;
    }
}

// get the Hilbert space defined by the running options and a given Total Lz value and for the SU(4) case
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpaceSU4(int totalLz)
{
  if (this->BosonFlag == false)
    {
      int NbrFermions = ((SingleIntegerOption*) (*(this->Options))["nbr-particles"])->GetInteger();
      int LzMax = ((SingleIntegerOption*) (*(this->Options))["lzmax"])->GetInteger();
      int SzTotal = ((SingleIntegerOption*) (*(this->Options))["total-sz"])->GetInteger();
      int IsoSzTotal = ((SingleIntegerOption*) (*(this->Options))["total-isosz"])->GetInteger();
      int TotalEntanglement = ((SingleIntegerOption*) (*(this->Options))["total-entanglement"])->GetInteger();
      unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) (*(this->Options))["fast-search"])->GetInteger()) << 20;
      ParticleOnSphereWithSU4Spin* Space;
#ifdef __64_BITS__
      if (LzMax <= 15)
#else
      if (LzMax <= 7)
#endif
        {
	  if (((BooleanOption*) (*(this->Options))["use-entanglement"])->GetBoolean())
	    Space = new FermionOnSphereWithSU4Spin(NbrFermions, totalLz, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
	  else
	    Space = new FermionOnSphereWithSU4Spin(NbrFermions, totalLz, LzMax, SzTotal, IsoSzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return 0;
	}	
      if ((((BooleanOption*) (*(this->Options))["use-entanglement"])->GetBoolean()) && (Space->GetHilbertSpaceDimension() == 0))
	{
	  cout << "zero dimension Hilbert space" << endl;
	  return 0;	  
	}
      return Space;
    }
  else
    {
      return 0;
    }
}

