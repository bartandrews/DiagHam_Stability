
#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Operator/ParticleOnSphereCreationOperator.h"
#include "Operator/ParticleOnSphereAnnihilationOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereLMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += DataGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-lz", "use full Hilbert space (all Lz sectors)");
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "orbital-number", "number of the orbital where a particle will be created (0, 1, 2, ..., lzmax)", 0);

  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file");
  (*SystemGroup) += new BooleanOption  ('\n', "annihilate-particle", "annihilate particle instead of creating it");

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");
  (*DataGroup) += new SingleStringOption  ('\n', "interaction-name", "interaction name for the output files when computing more than one state", "lminus");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLMinus -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  int OrbitalNumber = Manager.GetInteger("orbital-number");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, Lz, FermionFlag) == false)
    {
      return -1;
    }
  if ((Manager.GetString("statistics")) != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	  {
	    FermionFlag = false;
	  }
	else
	  {
	    cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	  }  
    }
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphere* IntialSpace;
  ParticleOnSphere* TargetSpace;
  cout << "Creating initial space"<<endl;
  if (FermionFlag == true)
    {
     if (Manager.GetBoolean("all-lz")) 
          IntialSpace = new FermionOnSphereFull(NbrParticles, LzMax);
     else
      {
         cout << "initial space Lz="<<Lz<<"/2"<<endl;
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
     }
    }
  else
    {
      IntialSpace = new BosonOnSphereShort(NbrParticles, Lz, LzMax);
    }


  if (Manager.GetBoolean("annihilate-particle") == false)
   {
         cout << "Creating target space"<<endl;

         if (FermionFlag == true)
	 {
           if (Manager.GetBoolean("all-lz"))              
              TargetSpace = new FermionOnSphereFull(NbrParticles + 1, LzMax);
          else
           {
             cout << "target Lz="<<Lz + (2 * OrbitalNumber - LzMax)<<"/2"<<endl;
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#endif
         }
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}

      cout << "Initial dim= " << IntialSpace->GetHilbertSpaceDimension() << " Target dim= " << TargetSpace->GetHilbertSpaceDimension() << endl;


      Architecture.GetArchitecture()->SetDimension(IntialSpace->GetHilbertSpaceDimension());

      cout << "computing c^+_"<< OrbitalNumber << " |Psi>" << endl;
      ParticleOnSphereCreationOperator TmpOperator(IntialSpace, OrbitalNumber);
      VectorOperatorMultiplyOperation Operation(&TmpOperator, &InitialVector, &TargetVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());

    } 	 
    else //annihilate particle
    {

      cout << "Creating target space "<<endl;

      if (FermionFlag == true)
	{

           if (Manager.GetBoolean("all-lz"))              
              TargetSpace = new FermionOnSphereFull(NbrParticles - 1, LzMax);
           else
           {
              cout << "target space Lz="<<Lz - (2 * OrbitalNumber - LzMax)<<"/2"<<endl;
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#endif
         }  
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}

      cout << "Initial dim= " << IntialSpace->GetHilbertSpaceDimension() << " Target dim= " << TargetSpace->GetHilbertSpaceDimension() << endl;

      Architecture.GetArchitecture()->SetDimension(IntialSpace->GetHilbertSpaceDimension());

      cout << "computing c_"<< OrbitalNumber << " |Psi>" << endl;
      ParticleOnSphereAnnihilationOperator TmpOperator(IntialSpace, OrbitalNumber);
      VectorOperatorMultiplyOperation Operation(&TmpOperator, &InitialVector, &TargetVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());


    }	

    char *OutputName;
    if (Manager.GetString("output-file")!=NULL)
	{
	  OutputName = new char[strlen(Manager.GetString("output-file"))+1];
	  strcpy(OutputName,Manager.GetString("output-file"));
	}
    else
	{
     OutputName = new char [256];
     int FinalLz, FinalNbrParticles;
     if (Manager.GetBoolean("annihilate-particle") == false)
        { 
          FinalNbrParticles = NbrParticles + 1; 
          FinalLz = Lz + (2 * OrbitalNumber - LzMax);
        }  
     else
        {
          FinalNbrParticles = NbrParticles - 1;  
          FinalLz = Lz - (2 * OrbitalNumber - LzMax);
        } 
     sprintf (OutputName, "fermions_addremoveparticle_n_%d_2s_%d_lz_%d.0.vec", FinalNbrParticles, LzMax, FinalLz);
	}  
    if (TargetVector.WriteVector(OutputName) == false)
	{
	  cout << "error while writing " << OutputName << endl;
	  return -1;
	}
    cout<<"Norm= " << TargetVector.Norm() << endl;
    delete [] OutputName;

  delete IntialSpace;
  return 0;
}


/* //Old code, working but redundant
#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Operator/ParticleOnSphereCreationOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereLMinus" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* DataGroup = new OptionGroup ("data options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += DataGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz", "twice the total lz value of the system for the initial state", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-lz", "use full Hilbert space (all Lz sectors)");
  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (boson or fermion, try to guess it from file name if not defined)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "orbital-number", "number of the orbital where a particle will be created (0, 1, 2, ..., lzmax)", 0);

  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a haldane basis with the given reference file for the input file");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a haldane basis with the given reference file for the output file");
  (*SystemGroup) += new BooleanOption  ('\n', "annihilate-particle", "annihilate particle instead of creating it");

  (*DataGroup) += new SingleStringOption  ('i', "input-file", "input vector file name");
  (*DataGroup) += new SingleStringOption  ('o', "output-file", "output vector file name");
  (*DataGroup) += new SingleStringOption  ('\n', "interaction-name", "interaction name for the output files when computing more than one state", "lminus");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLMinus -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz");
  int OrbitalNumber = Manager.GetInteger("orbital-number");
  bool FermionFlag = false;
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-file"), NbrParticles, LzMax, Lz, FermionFlag) == false)
    {
      return -1;
    }
  if ((Manager.GetString("statistics")) != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	  {
	    FermionFlag = false;
	  }
	else
	  {
	    cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	  }  
    }
  int Parity = Lz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }

  RealVector InitialVector; 
  RealVector TargetVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
  long MemorySpace = 9l << 20;
  ParticleOnSphere* IntialSpace;
  ParticleOnSphere* TargetSpace;
  cout << "Creating initial space"<<endl;
  if (FermionFlag == true)
    {
     if (Manager.GetBoolean("all-lz")) 
          IntialSpace = new FermionOnSphereFull(NbrParticles, LzMax);
     else
      {
         cout << "initial space Lz="<<Lz<<"/2"<<endl;
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#else
      if (LzMax <= 31)
	{
	  IntialSpace = new FermionOnSphere(NbrParticles, Lz, LzMax, MemorySpace);
	}
      else
	{
	  IntialSpace = new FermionOnSphereUnlimited(NbrParticles, Lz, LzMax, MemorySpace);
	}
#endif
     }
    }
  else
    {
      IntialSpace = new BosonOnSphereShort(NbrParticles, Lz, LzMax);
    }


  if (Manager.GetBoolean("annihilate-particle") == false)
   {
         cout << "Creating target space"<<endl;

         if (FermionFlag == true)
	 {
           if (Manager.GetBoolean("all-lz"))              
              TargetSpace = new FermionOnSphereFull(NbrParticles + 1, LzMax);
          else
           {
             cout << "target Lz="<<Lz + (2 * OrbitalNumber - LzMax)<<"/2"<<endl;
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#endif
         }
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles + 1, Lz + (2 * OrbitalNumber - LzMax), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}

      cout << "Initial dim= " << IntialSpace->GetHilbertSpaceDimension() << " Target dim= " << TargetSpace->GetHilbertSpaceDimension() << endl;

      double TmpCoefficient = 0.0;
      int TmpIndex;
      int* TmpOrbitals = new int[LzMax];
      for (int i = 0; i < IntialSpace->GetHilbertSpaceDimension(); ++i)
       {
         //Get the orbital occupancies of the initial vector and construct a representative word
	 IntialSpace->GetOccupied(i, TmpOrbitals);
	 //for (int k = 0; k < NbrParticles; k++) cout << TmpOrbitals[k] << " ";

	 unsigned long TmpState = 0x0ul;
         for (int k = 0; k < NbrParticles; ++k)
            TmpState |= 0x1ul << (TmpOrbitals[k]);

         //int TmpLzMax = LzMax;
         //while ((TmpState >> TmpLzMax) == 0x0ul)
         //  --TmpLzMax;

 	 //Act with c_M^+
         TmpCoefficient = 0.0;
	 TmpState = IntialSpace->Ad(TmpState, OrbitalNumber, TmpCoefficient);
	 if (TmpCoefficient != 0.0)
	  {
             //Get occupied orbitals of new state TmpState
             int k = 0;
             for (int l = 0; l <= LzMax; ++l)
               if ((TmpState >> l) & 0x1ul)
                  TmpOrbitals[k++] = l;

            //Find the new vector in a target Hilbert space 
            TmpIndex = TargetSpace->FindStateIndex(TmpOrbitals);
            //cout<<"  i= "<<i<<" ;  ->    Index = " << TmpIndex<<"  ; ";
	    //TargetSpace->GetOccupied(TmpIndex, TmpOrbitals); 
  	    //for (int k = 0; k < (NbrParticles+1); k++) cout << TmpOrbitals[k] << " ";
            //cout<<" coeff= "<<InitialVector[i] * TmpCoefficient;

	    if (TmpIndex < TargetVector.GetVectorDimension())
	       {
	          TargetVector[TmpIndex] += InitialVector[i] * TmpCoefficient;
	       }
          }
         //cout<<endl;
       }
      delete[] TmpOrbitals; 
	} 	 
    else //annihilate particle
    {

      cout << "Creating target space "<<endl;

      if (FermionFlag == true)
	{

           if (Manager.GetBoolean("all-lz"))              
              TargetSpace = new FermionOnSphereFull(NbrParticles - 1, LzMax);
           else
           {
              cout << "target space Lz="<<Lz - (2 * OrbitalNumber - LzMax)<<"/2"<<endl;
#ifdef __64_BITS__
	  if (LzMax <= 63)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#else
	  if (LzMax <= 31)
	    {
	      TargetSpace = new FermionOnSphere(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
	  else
	    {
	      TargetSpace = new FermionOnSphereUnlimited(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax, MemorySpace);
	    }
#endif
         }  
	}
      else
	{
	  TargetSpace = new BosonOnSphereShort(NbrParticles - 1, Lz - (2 * OrbitalNumber - LzMax), LzMax);
	}
      IntialSpace->SetTargetSpace(TargetSpace);
      TargetVector = RealVector(TargetSpace->GetHilbertSpaceDimension(), true);
      if (TargetSpace->GetHilbertSpaceDimension()!=IntialSpace->GetTargetHilbertSpaceDimension())
	{
	  cout << "Problem with setting target space"<<endl;
	  exit(-1);
	}

      cout << "Initial dim= " << IntialSpace->GetHilbertSpaceDimension() << " Target dim= " << TargetSpace->GetHilbertSpaceDimension() << endl;

      double TmpCoefficient = 0.0;
      int TmpIndex;
      int* TmpOrbitals = new int[LzMax];
      for (int i = 0; i < IntialSpace->GetHilbertSpaceDimension(); ++i)
       {
             //Get the orbital occupancies of the initial vector and construct a representative word
	     //IntialSpace->GetOccupied(i, TmpOrbitals);
	     //for (int k = 0; k < NbrParticles; k++) cout << TmpOrbitals[k] << " ";

	     //unsigned long TmpState = 0x0ul;
             //for (int k = 0; k < NbrParticles; ++k)
             //  TmpState |= 0x1ul << (TmpOrbitals[k]);

            //int TmpLzMax = LzMax;
            //while ((TmpState >> TmpLzMax) == 0x0ul)
            //  --TmpLzMax;

 	    //Act with c_M
             TmpCoefficient = 0.0;
	     TmpIndex = IntialSpace->A(i, OrbitalNumber, TmpCoefficient);
             //cout<<"Tmp index="<<TmpIndex<<endl;
	     if (TmpCoefficient != 0.0)
	      {
                //Get occupied orbitals of new state TmpState
                //int k = 0;
                //for (int l = 0; l <= LzMax; ++l)
                //   if ((TmpState >> l) & 0x1ul)
                //      TmpOrbitals[k++] = l;

               //Find the new vector in a target Hilbert space 
               //TmpIndex = TargetSpace->FindStateIndex(TmpOrbitals);
               //cout<<"  i= "<<i<<" ;  ->    Index = " << TmpIndex<<"  ; ";
	       //TargetSpace->GetOccupied(TmpIndex, TmpOrbitals); 
  	       //for (int k = 0; k < (NbrParticles-1); k++) cout << TmpOrbitals[k] << " ";
               //   cout<<" coeff= "<<InitialVector[i] * TmpCoefficient;

	       if (TmpIndex < TargetVector.GetVectorDimension())
	        {
	          TargetVector[TmpIndex] += InitialVector[i] * TmpCoefficient;
	        }
            }
            //cout<<endl;
       }
      delete[] TmpOrbitals; 

    }	

    char *OutputName;
    if (Manager.GetString("output-file")!=NULL)
	{
	  OutputName = new char[strlen(Manager.GetString("output-file"))+1];
	  strcpy(OutputName,Manager.GetString("output-file"));
	}
    else
	{
     OutputName = new char [256];
     int FinalLz, FinalNbrParticles;
     if (Manager.GetBoolean("annihilate-particle") == false)
        { 
          FinalNbrParticles = NbrParticles + 1; 
          FinalLz = Lz + (2 * OrbitalNumber - LzMax);
        }  
     else
        {
          FinalNbrParticles = NbrParticles - 1;  
          FinalLz = Lz - (2 * OrbitalNumber - LzMax);
        } 
     sprintf (OutputName, "fermions_addremoveparticle_n_%d_2s_%d_lz_%d.0.vec", FinalNbrParticles, LzMax, FinalLz);
	}  
    if (TargetVector.WriteVector(OutputName) == false)
	{
	  cout << "error while writing " << OutputName << endl;
	  return -1;
	}
    delete [] OutputName;

  delete IntialSpace;
  return 0;
}

*/
