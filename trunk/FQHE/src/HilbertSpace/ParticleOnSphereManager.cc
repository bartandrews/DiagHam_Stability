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

#include "Options/Options.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereEfficient.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinAndPairing.h"
#include "HilbertSpace/FermionOnSphereWithSpinPartialPolarization.h"

#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLzSzSymmetry.h"

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
#include "HilbertSpace/BosonOnSphereHaldaneSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinLzSzSymmetry.h"
#include "HilbertSpace/BosonOnSphereWithSpinOld.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/BosonOnSphereWithSU2SpinPartialPolarization.h"

#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU4SpinAllEntanglement.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/ConfigurationParser.h"


#include <stdlib.h>
#include <cstring>
#include <iostream>


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

void ParticleOnSphereManager::AddOptionGroup(OptionManager* manager, const char* comment)
{
  this->Options = manager;
  char tmpC[255];
  if (comment==NULL)
    sprintf(tmpC,"system options");
  else
    sprintf(tmpC,"system options (%s)",comment);  
  OptionGroup* SystemGroup  = new OptionGroup (tmpC);
  if (comment==NULL)
    sprintf(tmpC,"precalculation options");
  else
    sprintf(tmpC,"precalculation options (%s)",comment);    
  OptionGroup* PrecalculationGroup = new OptionGroup (tmpC);
  (*(this->Options)) += SystemGroup;
  (*(this->Options)) += PrecalculationGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  if ((this->FermionFlag == true) && (this->BosonFlag == true))
    (*SystemGroup) += new SingleStringOption  ('\n', "statistics", "particle statistics (bosons or fermions)", "fermions");
  switch (this->SUKIndex)
    {
    case 1:
      {
	(*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)"); 
	(*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
	(*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
	if (this->BosonFlag == false)
	  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
	(*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
	if (this->FermionFlag == true)
	  {
	    (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
	    (*PrecalculationGroup) += new BooleanOption  ('\n', "use-efficient", "use exhaustive search table");
	  }
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
	(*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
	(*PrecalculationGroup) += new SingleIntegerOption  ('\n', "huge-memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
      }
      break;
    case 2:
      {
	(*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
	(*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0)");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
	(*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");
	(*SystemGroup) += new SingleIntegerOption ('\n', "nbrspin-polarized", "number of orbitals which ar fully spin up polarized (from the one with the lowest momentum)", 0);
	if (this->FermionFlag == true) // symmetrized bases not defined for bosons at the moment
	  {
	    (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
	    (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
	  }
	if (this->BosonFlag == true)
	  {
	    // boson options
	    (*SystemGroup) += new BooleanOption  ('\n', "all-sz", "use Hilbert-space with all values of sz");
	    (*SystemGroup) += new SingleIntegerOption  ('\n', "pair-parity", "parity for N_up as compared to int(N/2) (0=same, 1=different, -1=none)", -1);
	    (*PrecalculationGroup) += new BooleanOption  ('\n', "use-old", "use full integer representation of bosonic states (slow)");
	    (*PrecalculationGroup) += new BooleanOption  ('\n', "use-alt", "use alternative Hilbert space for  bosonic states");
	  }
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
	(*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the total Sz projection", 0);
	(*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the total isospin operator", 0);
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
  if ((this->FermionFlag == true) && (this->BosonFlag == true) && (this->Options->GetString("statistics") != 0))
    {
      if ((strcmp ("fermions", this->Options->GetString("statistics")) == 0))
	{
	  this->FermionFlag = true;
	  this->BosonFlag = false;
	}
      else
	if ((strcmp ("bosons", this->Options->GetString("statistics")) == 0))
	  {
	    this->FermionFlag = false;
	    this->BosonFlag = true;
	  }
	else
	  {
	    cout << this->Options->GetString("statistics") << " is an undefined statistics" << endl;
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
  bool SymmetrizedBasis = this->Options->GetBoolean("symmetrized-basis");
  int NbrParticles = this->Options->GetInteger("nbr-particles");
  int LzMax = this->Options->GetInteger("lzmax");
  bool HaldaneBasisFlag = this->Options->GetBoolean("haldane");
  if (this->BosonFlag == false)
    {
      unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;
      ParticleOnSphere* Space = 0;
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	  if (LzMax <= 30)
#endif
	    if ((SymmetrizedBasis == false) || (totalLz != 0))
	      {
		if (this->Options->GetBoolean("use-efficient"))
		  Space = new FermionOnSphereEfficient(NbrParticles, totalLz, LzMax, MemorySpace);
		else
		  Space = new FermionOnSphere(NbrParticles, totalLz, LzMax, MemorySpace);
	      }
	    else
	      {
		if (this->Options->GetString("load-hilbert") != 0)
		  Space = new FermionOnSphereSymmetricBasis(this->Options->GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
		if (this->Options->GetString("save-hilbert") != 0)
		  {
		    ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
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
		      if (this->Options->GetString("load-hilbert") != 0)
			Space = new FermionOnSphereSymmetricBasisLong(this->Options->GetString("load-hilbert"), MemorySpace);
		      else
			Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		      if (this->Options->GetString("save-hilbert") != 0)
			{
			  ((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
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
	  if (this->Options->GetString("reference-file") == 0)
	    {
	      ReferenceState = new int[LzMax + 1];
	      for (int i = 0; i <= LzMax; ++i)
		ReferenceState[i] = 0;
	      if (strcasecmp(this->Options->GetString("reference-state"), "laughlin") == 0)
		for (int i = 0; i <= LzMax; i += 3)
		  ReferenceState[i] = 1;
	      else
		if (strcasecmp(this->Options->GetString("reference-state"), "pfaffian") == 0)
		  for (int i = 0; i <= LzMax; i += 4)
		    {
		      ReferenceState[i] = 1;
		      ReferenceState[i + 1] = 1;
		    }
		else
		  if (strcasecmp(this->Options->GetString("reference-state"), "readrezayi3") == 0)
		    for (int i = 0; i <= LzMax; i += 5)
		      {
			ReferenceState[i] = 1;
			ReferenceState[i + 1] = 1;
			ReferenceState[i + 2] = 1;
		      }
		  else
		    {
		      cout << "unknown reference state " << this->Options->GetString("reference-state") << endl;
		      return 0;
		    }
	    }
	  else
	    {
	      ConfigurationParser ReferenceStateDefinition;
	      if (ReferenceStateDefinition.Parse(this->Options->GetString("reference-file")) == false)
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
		  cout << "error while parsing ReferenceState in " << this->Options->GetString("reference-file") << endl;
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
	       if (this->Options->GetBoolean("huge-basis") == false)
		 {
#ifdef __64_BITS__
		   if (LzMax <= 62)
#else
		     if (LzMax <= 30)
#endif
		       {
			 if (this->Options->GetString("load-hilbert") != 0)
			   Space = new FermionOnSphereHaldaneBasis(this->Options->GetString("load-hilbert"), MemorySpace);
			 else
			   Space = new FermionOnSphereHaldaneBasis(NbrParticles, totalLz, LzMax, ReferenceState, MemorySpace);
			 if (this->Options->GetString("save-hilbert") != 0)
			   {
			     ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			     cout << "Wrote hilbert space "<<this->Options->GetString("save-hilbert")<<" of dimension d="<<Space->GetHilbertSpaceDimension()<<endl;
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
			     if (this->Options->GetString("load-hilbert") != 0)
			       Space = new FermionOnSphereHaldaneBasisLong(this->Options->GetString("load-hilbert"), MemorySpace);
			     else
			       Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, totalLz, LzMax, ReferenceState, MemorySpace);
			     if (this->Options->GetString("save-hilbert") != 0)
			       {
				 ((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
				 return 0;
			       }
			   }	       
		 }
	       else
		 {
		   if (this->Options->GetString("load-hilbert") == 0)
		     {
		       cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		       return 0;
		     }
		   Space = new FermionOnSphereHaldaneHugeBasis (this->Options->GetString("load-hilbert"), this->Options->GetInteger("huge-memory"));
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
		     if (this->Options->GetString("load-hilbert") != 0)
		       Space = new FermionOnSphereHaldaneSymmetricBasis(this->Options->GetString("load-hilbert"), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		     if (this->Options->GetString("save-hilbert") != 0)
		       {
			 ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
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
			 if (this->Options->GetString("load-hilbert") != 0)
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(this->Options->GetString("load-hilbert"), MemorySpace);
			 else
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
			 if (this->Options->GetString("save-hilbert") != 0)
			   {
			     ((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			     return 0;
			   }
		       }
	    }
	}
      return Space;
    }
  else
    {
      ParticleOnSphere* Space = 0;
      if (HaldaneBasisFlag == false)
	{
	  bool SymmetrizedBasis = this->Options->GetBoolean("symmetrized-basis");
#ifdef  __64_BITS__
	  if ((LzMax + NbrParticles - 1) < 63)
#else
	    if ((LzMax + NbrParticles - 1) < 31)	
#endif
	      {
		if ((SymmetrizedBasis == false) || (totalLz != 0))
		  Space = new BosonOnSphereShort(NbrParticles, totalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasisShort(NbrParticles, LzMax);
	      }
	    else
	      {
		if ((SymmetrizedBasis == false) || (totalLz != 0))
		  Space = new BosonOnSphere (NbrParticles, totalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasis(NbrParticles, LzMax);
	      }
	}
      else 
	{
	  int* ReferenceState = 0;
	  if (this->Options->GetString("reference-file") == 0)
	    {
	      cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
	      return 0;
	    }
	  ConfigurationParser ReferenceStateDefinition;
	  if (ReferenceStateDefinition.Parse(this->Options->GetString("reference-file")) == false)
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
	      cout << "error while parsing ReferenceState in " << this->Options->GetString("reference-file") << endl;
	      return 0;     
	    }
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return 0;     
	    }
	  if (SymmetrizedBasis == false)
	     {
	       if (this->Options->GetBoolean("huge-basis") == false)
		 {
#ifdef  __64_BITS__
		   if ((LzMax + NbrParticles - 1) < 63)
#else
		     if ((LzMax + NbrParticles - 1) < 31)	
#endif
		       {
			 if (this->Options->GetString("load-hilbert") != 0)
			   Space = new BosonOnSphereHaldaneBasisShort(this->Options->GetString("load-hilbert"));
			 else
			   Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, totalLz, LzMax, ReferenceState);
			 if (this->Options->GetString("save-hilbert") != 0)
			   {
			     ((BosonOnSphereHaldaneBasisShort*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			     cout << "Wrote hilbert space "<<this->Options->GetString("save-hilbert")<<" of dimension d="<<Space->GetHilbertSpaceDimension()<<endl;
			     return 0;
			   }
		       }
		 }
	       else
		 {
		   if (this->Options->GetString("load-hilbert") == 0)
		     {
		       cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		       return 0;
		     }
		   Space = new BosonOnSphereHaldaneHugeBasisShort (this->Options->GetString("load-hilbert"), this->Options->GetInteger("huge-memory"));
		 }
	     }
	  else
	    {
#ifdef  __64_BITS__
	      if ((LzMax + NbrParticles - 1) < 62)
#else
		if ((LzMax + NbrParticles - 1) < 30)	
#endif
		  {
		    if (this->Options->GetString("load-hilbert") != 0)
		      Space = new BosonOnSphereHaldaneSymmetricBasisShort(this->Options->GetString("load-hilbert"));
		    else
		      Space = new BosonOnSphereHaldaneSymmetricBasisShort(NbrParticles, LzMax, ReferenceState);
		    if (this->Options->GetString("save-hilbert") != 0)
		      {
			((BosonOnSphereHaldaneSymmetricBasisShort*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			cout << "Wrote hilbert space "<<this->Options->GetString("save-hilbert")<<" of dimension d="<<Space->GetHilbertSpaceDimension()<<endl;
			return 0;
		      }
		  }
	    }
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
      int NbrFermions = this->Options->GetInteger("nbr-particles");
      int LzMax = this->Options->GetInteger("lzmax");
      int SzTotal = this->Options->GetInteger("total-sz");
      bool LzSymmetrizedBasis = this->Options->GetBoolean("lzsymmetrized-basis");
      bool SzSymmetrizedBasis = this->Options->GetBoolean("szsymmetrized-basis");
      unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;
      if (this->Options->GetBoolean("haldane") == false)
	{
	  if (((SzSymmetrizedBasis == false) || (SzTotal != 0)) && ((LzSymmetrizedBasis == false) || (totalLz != 0)))
	    {
#ifdef __64_BITS__
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		      if (this->Options->GetInteger("nbrspin-polarized") > 0)
			{
			  Space = new FermionOnSphereWithSpinPartialPolarization(NbrFermions, totalLz, LzMax, SzTotal, 
										 (int) this->Options->GetInteger("nbrspin-polarized"), MemorySpace);
			}
		      else
			{
			  Space = new FermionOnSphereWithSpin(NbrFermions, totalLz, LzMax, SzTotal, MemorySpace);
			}
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
			  if (this->Options->GetString("load-hilbert") == 0)
			    Space = new FermionOnSphereWithSpinSzSymmetry(NbrFermions, totalLz, LzMax, this->Options->GetBoolean("minus-szparity"), MemorySpace);
			  else
			    Space = new FermionOnSphereWithSpinSzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);
			}
		      else
			{
			  if (this->Options->GetString("load-hilbert") == 0)
			    Space = new FermionOnSphereWithSpinSzSymmetryLong(NbrFermions, totalLz, LzMax, this->Options->GetBoolean("minus-szparity"), MemorySpace);
			  else
			    Space = new FermionOnSphereWithSpinSzSymmetryLong(this->Options->GetString("load-hilbert"), MemorySpace);
			}
		  }
		else
#ifdef __64_BITS__
		  if (LzMax <= 28)
#else
		    if (LzMax <= 13)
#endif
		      {
			if (this->Options->GetString("load-hilbert") == 0)
			  {
			    Space = new FermionOnSphereWithSpinLzSzSymmetry(NbrFermions, LzMax, this->Options->GetBoolean("minus-szparity"),
									    this->Options->GetBoolean("minus-lzparity"), MemorySpace);
			  }
			else
			  Space = new FermionOnSphereWithSpinLzSzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);
		      }
		    else
		      {
			if (this->Options->GetString("load-hilbert") == 0)
			  {
			    Space = new FermionOnSphereWithSpinLzSzSymmetryLong(NbrFermions, LzMax, this->Options->GetBoolean("minus-szparity"),
										this->Options->GetBoolean("minus-lzparity"), MemorySpace);
			  }
			else
			  Space = new FermionOnSphereWithSpinLzSzSymmetryLong(this->Options->GetString("load-hilbert"), MemorySpace);
			
		      }
		  else
#ifdef __64_BITS__
		    if (LzMax <= 28)
#else
		      if (LzMax <= 13)
#endif
			{
			  if (this->Options->GetString("load-hilbert") == 0)
			    Space = new FermionOnSphereWithSpinLzSymmetry(NbrFermions, LzMax, SzTotal, this->Options->GetBoolean("minus-lzparity"), MemorySpace);
			  else
			    Space = new FermionOnSphereWithSpinLzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);	      
			}
		      else
			{
			  if (this->Options->GetString("load-hilbert") == 0)
			    Space = new FermionOnSphereWithSpinLzSymmetryLong(NbrFermions, LzMax, SzTotal, this->Options->GetBoolean("minus-lzparity"), MemorySpace);
			  else
			    Space = new FermionOnSphereWithSpinLzSymmetryLong(this->Options->GetString("load-hilbert"), MemorySpace);	      
			}
	      if (this->Options->GetString("save-hilbert") != 0)
		{
		  ((FermionOnSphereWithSpinLzSzSymmetry*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
		  return 0;
		}
	    }
	  return Space;
	}
      else
	{
	  int** ReferenceStates = 0;
	  int NbrReferenceStates;
	  bool TexturelessFlag; 
	  if (this->Options->GetString("reference-file") == 0)
	    {
	      cout << "error, a reference file is needed for fermions in Haldane basis" << endl;
	      return 0;
	    }
	  if (FQHEGetRootPartitionSU2(this->Options->GetString("reference-file"), NbrFermions, LzMax, ReferenceStates, NbrReferenceStates, TexturelessFlag) == false)
	    {
	      cout << "error while parsing " << this->Options->GetString("reference-file") << endl;	      
	      return 0;
	    }
	  if ((SzSymmetrizedBasis == true) && (SzTotal == 0) && (LzSymmetrizedBasis == true) && (totalLz == 0))
	    {	      	      
#ifdef __64_BITS__
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		    if ( TexturelessFlag == false )
		      {		    
			Space = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrFermions, LzMax, this->Options->GetBoolean("minus-szparity"), this->Options->GetBoolean("minus-lzparity"), ReferenceStates, NbrReferenceStates);
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
			Space = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrFermions, LzMax, this->Options->GetBoolean("minus-szparity"), this->Options->GetBoolean("minus-lzparity"), TexturelessReferenceState, NbrReferenceStates, true);
		      }
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
	      if (LzMax <= 31)
#else
		if (LzMax <= 15)
#endif
		  {
		    if ( TexturelessFlag == false ) 
		      {
			  Space = new FermionOnSphereWithSpinHaldaneBasis(NbrFermions, totalLz, LzMax, SzTotal, ReferenceStates, NbrReferenceStates);
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
			Space = new FermionOnSphereWithSpinHaldaneBasis(NbrFermions, totalLz, LzMax, SzTotal, TexturelessReferenceState, NbrReferenceStates, true);
		      }
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (LzMax <= 63)
#else
		      if (LzMax <= 31)
#endif
			{
			  Space = new FermionOnSphereWithSpinHaldaneBasisLong (NbrFermions, totalLz, LzMax, SzTotal, ReferenceStates, NbrReferenceStates);
			}
		      else
			{
			  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
			  return 0;
			}	
		  }
	    }
	  return Space;
	}	
    }
  else
    {
      ParticleOnSphereWithSpin* Space = 0;
      int NbrBosons = this->Options->GetInteger("nbr-particles");
      int LzMax = this->Options->GetInteger("lzmax");
      bool SzSymmetrizedBasis = this->Options->GetBoolean("szsymmetrized-basis");
      bool LzSymmetrizedBasis = this->Options->GetBoolean("lzsymmetrized-basis");
      int SzTotal = this->Options->GetInteger("total-sz");
      unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;
      if  ((LzSymmetrizedBasis == false) || (totalLz != 0))
	{
	  if (SzSymmetrizedBasis == true)
	    {
	      if (this->Options->GetBoolean("all-sz"))
		{
		  if (this->Options->GetString("load-hilbert") == 0)
		    {
		      Space = new BosonOnSphereWithSU2SpinSzSymmetry(NbrBosons, totalLz, LzMax, this->Options->GetBoolean("minus-szparity"));
		    }
		  else
		    {
		      Space = new BosonOnSphereWithSU2SpinSzSymmetry(this->Options->GetString("load-hilbert"));
		    }
		}
	      else
		{
		  Space = new BosonOnSphereWithSU2SpinSzSymmetry(NbrBosons, totalLz, LzMax, SzTotal, this->Options->GetBoolean("minus-szparity"), MemorySpace);
		}
	    }
	  else
	    {
	      if (this->Options->GetBoolean("all-sz"))
		{
		  if (this->Options->GetBoolean("use-alt"))
		    {
		      Space = new BosonOnSphereWithSU2Spin(NbrBosons, totalLz, LzMax);
		    }
		  else
		    {
		      if (this->Options->GetInteger("pair-parity") >= 0)
			Space = new BosonOnSphereWithSpinAllSz(NbrBosons, totalLz, LzMax, this->Options->GetInteger("pair-parity"), MemorySpace);
		      else
			Space = new BosonOnSphereWithSpinAllSz(NbrBosons, totalLz, LzMax, MemorySpace);
		    }
		}
	      else
		{
		  unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;
		  if (this->Options->GetBoolean("use-old"))
		    {
		      Space = new BosonOnSphereWithSpinOld(NbrBosons, totalLz, LzMax, SzTotal);
		    }
		  else
		    {
		      if (this->Options->GetInteger("nbrspin-polarized") > 0)
			{
			  if (this->Options->GetString("load-hilbert") == 0)
			    {
			      Space = new BosonOnSphereWithSU2SpinPartialPolarization(NbrBosons, totalLz, LzMax, SzTotal, 
										      (int) this->Options->GetInteger("nbrspin-polarized"), MemorySpace);
			    }
			  else
			    {
			      Space = new BosonOnSphereWithSU2SpinPartialPolarization(this->Options->GetString("load-hilbert"), MemorySpace);
			    }
			}
		      else
			{
			  if (this->Options->GetBoolean("use-alt"))
			    {
			      if (this->Options->GetString("load-hilbert") == 0)
				{
				  Space = new BosonOnSphereWithSU2Spin(NbrBosons, totalLz, LzMax, SzTotal, MemorySpace);
				}
			      else
				{
				  Space = new BosonOnSphereWithSU2Spin(this->Options->GetString("load-hilbert"), MemorySpace);
				}
			    }
			  else
			    {
			      Space = new BosonOnSphereWithSpin(NbrBosons, totalLz, LzMax, SzTotal, MemorySpace);
			    }
			}
		    }
		}
	    }
	}
      else
	{
	  if (SzSymmetrizedBasis == true)
	    {
	      if (this->Options->GetBoolean("all-sz"))
		{
		  Space = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrBosons, LzMax, 
								   this->Options->GetBoolean("minus-szparity"), this->Options->GetBoolean("minus-lzparity"));
		}
	      else
		{
		  if (this->Options->GetString("load-hilbert") == 0)
		    {
		      Space = new BosonOnSphereWithSU2SpinLzSzSymmetry(NbrBosons, LzMax, SzTotal, this->Options->GetBoolean("minus-szparity"), this->Options->GetBoolean("minus-lzparity"), MemorySpace);
		    }
		  else
		    {
		      Space = new BosonOnSphereWithSU2SpinLzSzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);
		    }
		}
	    }
	  else
	    {
	      if (this->Options->GetBoolean("all-sz"))
		{
		  Space = new BosonOnSphereWithSU2SpinLzSymmetry(NbrBosons, LzMax, this->Options->GetBoolean("minus-lzparity"));
		}
	      else
		{
		  if (this->Options->GetString("load-hilbert") == 0)
		    {
		      Space = new BosonOnSphereWithSU2SpinLzSymmetry(NbrBosons, LzMax, SzTotal, this->Options->GetBoolean("minus-lzparity"));
		    }
		  else
		    {
		      Space = new BosonOnSphereWithSU2SpinLzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);
		    }
		}
	    }
	}
      if (this->Options->GetString("save-hilbert") != 0)
	{
	  Space->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
	}
      return Space;
    }
  return 0;
}

// get the Hilbert space defined by the running options and a given Total Lz value and for the SU(3) case
//
// totalLz = twice the system total Lz value
// return value = pointer to the Lanczos algorithm

ParticleOnSphere* ParticleOnSphereManager::GetHilbertSpaceSU3(int totalLz)
{
  if (this->BosonFlag == false)
    {
      int NbrFermions = this->Options->GetInteger("nbr-particles");
      int LzMax = this->Options->GetInteger("lzmax");
      int TotalTz = this->Options->GetInteger("total-tz");
      int TotalY = this->Options->GetInteger("total-y");
//      bool LzSymmetrizedBasis = this->Options->GetBoolean("lzsymmetrized-basis");
      bool TzSymmetrizedBasis = this->Options->GetBoolean("tzsymmetrized-basis");
      bool Z3SymmetrizedBasis = this->Options->GetBoolean("z3symmetrized-basis");
      unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;      
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
	      if (this->Options->GetString("load-hilbert") == 0)
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(NbrFermions, totalLz, LzMax, TotalY, this->Options->GetBoolean("minus-tzparity"), MemorySpace);
		  if (this->Options->GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereWithSU3SpinTzSymmetry*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
		      return 0;
		    }
		}
	      else
		{
		  Space = new FermionOnSphereWithSU3SpinTzSymmetry(this->Options->GetString("load-hilbert"), MemorySpace);
		}
	    }
	  else
	    if ((TzSymmetrizedBasis == false) && (Z3SymmetrizedBasis == true))
	      {
		if (this->Options->GetString("load-hilbert") == 0)
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(NbrFermions, totalLz, LzMax, TotalTz, MemorySpace);
		    if (this->Options->GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereWithSU3SpinZ3Symmetry*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			return 0;
		      }
		  }
		else
		  {
		    Space = new FermionOnSphereWithSU3SpinZ3Symmetry(this->Options->GetString("load-hilbert"), MemorySpace);
		  }
	      }
	    else
	      if ((TzSymmetrizedBasis == true) && (Z3SymmetrizedBasis == true))
		{
		  if (this->Options->GetString("load-hilbert") == 0)
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(NbrFermions, totalLz, LzMax, this->Options->GetBoolean("minus-tzparity"), MemorySpace);
		      if (this->Options->GetString("save-hilbert") != 0)
			{
			  ((FermionOnSphereWithSU3SpinTzZ3Symmetry*) Space)->WriteHilbertSpace(this->Options->GetString("save-hilbert"));
			  return 0;
			}
		    }
		  else
		    {
		      Space = new FermionOnSphereWithSU3SpinTzZ3Symmetry(this->Options->GetString("load-hilbert"), MemorySpace);
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
  int NbrParticles = this->Options->GetInteger("nbr-particles");
  int LzMax = this->Options->GetInteger("lzmax");
  int SzTotal = this->Options->GetInteger("total-sz");
  int IsoSzTotal = this->Options->GetInteger("total-isosz");
  int TotalEntanglement = this->Options->GetInteger("total-entanglement");
  unsigned long MemorySpace = ((unsigned long) this->Options->GetInteger("fast-search")) << 20;
  if (this->BosonFlag == false)
    {
      ParticleOnSphereWithSU4Spin* Space;
#ifdef __64_BITS__
      if (LzMax <= 15)
#else
      if (LzMax <= 7)
#endif
        {
	  if (this->Options->GetBoolean("use-entanglement"))
	    Space = new FermionOnSphereWithSU4Spin(NbrParticles, totalLz, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
	  else
	    Space = new FermionOnSphereWithSU4Spin(NbrParticles, totalLz, LzMax, SzTotal, IsoSzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return 0;
	}	
      if ((this->Options->GetBoolean("use-entanglement")) && (Space->GetHilbertSpaceDimension() == 0))
	{
	  cout << "zero dimension Hilbert space" << endl;
	  return 0;	  
	}
      return Space;
    }
  else
    {
      ParticleOnSphereWithSU4Spin* Space;
      if (this->Options->GetBoolean("use-entanglement"))
	Space = new BosonOnSphereWithSU4Spin(NbrParticles, totalLz, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
      else
	Space = new BosonOnSphereWithSU4SpinAllEntanglement(NbrParticles, totalLz, LzMax, SzTotal, IsoSzTotal, MemorySpace);
      if (Space->GetHilbertSpaceDimension()==0)
	{
	  cout << "zero dimension Hilbert space" << endl;
	  return 0;
	}
      return Space;
    }
}

