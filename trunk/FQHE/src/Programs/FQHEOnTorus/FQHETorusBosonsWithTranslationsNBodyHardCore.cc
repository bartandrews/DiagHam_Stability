#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"


#include "Hamiltonian/ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/FQHEOnTorusMainTask.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "QuantumNumber/AbstractQuantumNumber.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHETorusBosonsWithTranslationsNBodyHardCore" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption ('l', "max-momentum", "maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-nbody", "number of particle that can interact simultaneously through the n-body hard-core interaction", 3);
  (*SystemGroup) += new SingleDoubleOption ('r', "ratio", "ratio between the two torus lengths", 1.0);
  (*SystemGroup) += new BooleanOption  ('\n', "all-points", "calculate all points", false);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new BooleanOption  ('\n', "full-reducedbz", "calculate all points within the full reduced Brillouin zone", false);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");

  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsWithTranslationsNBodyHardCore -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles");
  int MaxMomentum = Manager.GetInteger("max-momentum");
  int XMomentum = Manager.GetInteger("x-momentum");
  int YMomentum = Manager.GetInteger("y-momentum");
  double XRatio = Manager.GetDouble("ratio");
  int NbrNBody = Manager.GetInteger("nbr-nbody");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  bool FirstRun = true;
  
  char* OutputName = new char [256];
  sprintf (OutputName, "bosons_torus_%dbody_hardcore_n_%d_2s_%d_ratio_%f.dat", NbrNBody, NbrParticles, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);

  int MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  bool GenerateMomenta = false;
  if ((XMomentum < 0)||(YMomentum < 0))
    GenerateMomenta = true;
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
  int YMaxMomentum = (MaxMomentum - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum;

  int NbrMomenta;
  int* XMomenta;
  int* YMomenta;
  int* Multiplicities = NULL;
  int CenterX = 0;
  int CenterY = 0;
    
  if (GenerateMomenta == false)
    {
      NbrMomenta=1;
      XMomenta = new int[1];
      YMomenta = new int[1];
      XMomenta[0] = XMomentum;
      YMomenta[0] = YMomentum;
    }
  else
    {
      if (Manager.GetString("selected-points") == 0)
	{
	  if ((YMomentum != 0) || (XMomentum != 0))
	    {
	      if (YMomentum != 0)
		{
		  NbrMomenta = MomentumModulo;
		  YMomenta = new int[NbrMomenta];
		  XMomenta = new int[NbrMomenta];	  
		  for (int Kx = 0; Kx < NbrMomenta; ++Kx)
		    {
		      YMomenta[Kx] = YMomentum;
		      XMomenta[Kx] = Kx;
		    }
		}
	      if (XMomentum != 0)
		{
		  NbrMomenta = MomentumModulo;
		  YMomenta = new int[NbrMomenta];
		  XMomenta = new int[NbrMomenta];
		  for (int Ky = 0; Ky < NbrMomenta; ++Ky)
		    {
		      XMomenta [Ky] = XMomentum;
		      YMomenta[Ky] = Ky;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("all-points"))
		{
		  int Pos=0;
		  NbrMomenta = (XMaxMomentum - XMomentum+1) * (YMaxMomentum - YMomentum+1);
		  XMomenta = new int[NbrMomenta];
		  YMomenta = new int[NbrMomenta];
		  for (; XMomentum <= XMaxMomentum; ++XMomentum)
		    for (int YMomentum2 = YMomentum; YMomentum2<= YMaxMomentum; ++YMomentum2)
		      {
			XMomenta[Pos] = XMomentum;
			YMomenta[Pos] = YMomentum2;
			++Pos;
		      }
		}
	      else // determine inequivalent states in BZ
		{
		  if (Manager.GetBoolean("full-reducedbz"))
		    {
		      int Pos=0;
		      XMaxMomentum = MomentumModulo;
		      YMaxMomentum = MomentumModulo;
		      NbrMomenta = MomentumModulo * MomentumModulo;
		      XMomenta = new int[NbrMomenta];
		      YMomenta = new int[NbrMomenta];
		      for (; XMomentum < XMaxMomentum; ++XMomentum)
			for (int YMomentum2 = YMomentum; YMomentum2 < YMaxMomentum; ++YMomentum2)
			  {
			    XMomenta[Pos] = XMomentum;
			    YMomenta[Pos] = YMomentum2;
			    ++Pos;
			  }
		    }
		  else
		    {
		      CenterX=0;
		      CenterY=0;
		      if (XRatio == 1.0)
			{
			  NbrMomenta=0;
			  for (int Kx = CenterX; Kx <= CenterX+MomentumModulo/2; ++Kx)
			    for (int Ky= (Kx-CenterX) + CenterY; Ky <= CenterY+MomentumModulo/2; ++Ky)
			      {
				++NbrMomenta;
			      }
			  int Pos=0;
			  XMomenta = new int[NbrMomenta];
			  YMomenta = new int[NbrMomenta];
			  Multiplicities = new int[NbrMomenta];
			  for (int Kx = 0; Kx <= MomentumModulo/2; ++Kx)
			    for (int Ky = Kx; Ky <= MomentumModulo/2; ++Ky, ++Pos)
			      {
				XMomenta[Pos] = CenterX + Kx;
				YMomenta[Pos] = CenterY + Ky;
				if (Kx==0)
				  {
				    if (Ky==0)
				      Multiplicities[Pos] = 1; // BZ center
				    else if (Ky == MomentumModulo/2)
				      Multiplicities[Pos] = 2;
				    else Multiplicities[Pos] = 4;
				  }
				else if (Kx == MomentumModulo/2)
				  {
				    Multiplicities[Pos] = 1; // BZ corner
				  }
				else
				  {
				    if (Ky == Kx) // diagonal ?
				      {
					Multiplicities[Pos] = 4; 
				      }
				    else
				      {
					if (Ky == MomentumModulo/2)
					  Multiplicities[Pos] = 4;
					else
					  Multiplicities[Pos] = 8;
				      }
				  }
			      }
			}
		      else // rectangular torus
			{
			  NbrMomenta=(MomentumModulo/2+1) * (MomentumModulo/2+1);
			  int Pos = 0;
			  XMomenta = new int[NbrMomenta];
			  YMomenta = new int[NbrMomenta];
			  Multiplicities = new int[NbrMomenta];
			  for (int Kx = 0; Kx<=MomentumModulo/2; ++Kx)
			    for (int Ky= 0; Ky<=MomentumModulo/2; ++Ky, ++Pos)
			      {
				XMomenta[Pos] = CenterX + Kx;
				YMomenta[Pos] = CenterY + Ky;
				if (Kx == 0)
				  {
				    if (Ky == 0)
				      Multiplicities[Pos] = 1; // BZ center
				    else // on Gamma->X]
				      Multiplicities[Pos] = 2;
				  }
				else
				  {
				    if (Ky == 0)
				      Multiplicities[Pos] = 2;
				    else
				      {
					if (Kx == MomentumModulo/2)
					  {
					    if (Ky==MomentumModulo/2) // BZ corner?
					      Multiplicities[Pos] = 1;
					    else
					      Multiplicities[Pos] = 2;
					  }
					else
					  {
					    if (Ky == MomentumModulo/2) // on edge?
					      Multiplicities[Pos] = 2;
					    else
					      Multiplicities[Pos] = 4;
					  }
				      }
				  }
			      }
			}
		    }
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile MomentumFile;
	  if (MomentumFile.Parse(Manager.GetString("selected-points")) == false)
	    {
	      MomentumFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrMomenta = MomentumFile.GetNbrLines();
	  XMomenta = MomentumFile.GetAsIntegerArray(0);
	  YMomenta = MomentumFile.GetAsIntegerArray(1);
	}
    }
  
  for (int Pos = 0;Pos < NbrMomenta; ++Pos)
    {
      XMomentum = XMomenta[Pos];
      YMomentum = YMomenta[Pos];
      cout << "----------------------------------------------------------------" << endl;
      cout << " Ratio = " << XRatio << endl;

      BosonOnTorusWithMagneticTranslationsShort* Space = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      
      AbstractQHEHamiltonian* Hamiltonian = 0;
      Hamiltonian = new ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian(Space, NbrParticles, MaxMomentum, XMomentum, XRatio, 
											NbrNBody, Architecture.GetArchitecture(), Memory);

      double Shift = -1.0;
      Hamiltonian->ShiftHamiltonian(Shift);
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [512];
	  char* TmpName = RemoveExtensionFromFileName(OutputName, ".dat");
	  sprintf (EigenvectorName, "%s_kx_%d_ky_%d", TmpName, XMomentum, YMomentum);
	  delete [] TmpName;
	}
      FQHEOnTorusMainTask Task (&Manager, Space, &Lanczos, Hamiltonian, YMomentum, Shift, OutputName, FirstRun, EigenvectorName, XMomentum);
      MainTaskOperation TaskOperation (&Task);
      Task.SetKxValue(XMomentum);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      if (FirstRun == true)
	FirstRun = false;
      delete Hamiltonian;
      delete Space;
    }
  File.close();
  delete[] XMomenta;
  delete[] YMomenta;
  delete[] OutputName;
  
  return 0;
}
