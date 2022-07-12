#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneLargeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (n=0 being the LLL). WARNING: if using this option, l+1 should be the total number of orbitals that defines your input state (i.e., the flux Q will be determined from l=2Q+2n)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "p-truncated", "use a p-truncated basis instead of the full squeezed basis");
  (*SystemGroup) += new SingleIntegerOption ('\n', "p-truncation", "p-truncation for the p-truncated basis (if --p-truncated is used)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new BooleanOption  ('\n', "large-basis", "use large Hilbert space support (i.e. handle non-squeezed Hilbert space larger than 2^31 without hard-drive storage)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "pair-amplitude", "evaluate pair amplitudes", false);
  (*SystemGroup) += new BooleanOption  ('\n', "energy-expectation", "evaluate energy expectation value from pair amplitudes", false);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new BooleanOption  ('\n', "structure-factor", "evaluate the projected structure factor instead of (density-)density (use with radians option)", false);
  (*SystemGroup) += new BooleanOption  ('\n', "guidingcenter-structurefactor", "evaluate the guiding center structure factor instead of (density-)density", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "huge-memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "large-memory", "maximum memory (in kBytes) that can allocated for precalculations when using huge mode", 1);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int LandauLevel = Manager.GetInteger("landau-level");
  int NbrPoints = Manager.GetInteger("nbr-points");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool DensityFlag = Manager.GetBoolean("density");
  bool StructureFactorFlag = Manager.GetBoolean("structure-factor");
  bool GuidingCenterStructureFactorFlag = Manager.GetBoolean("guidingcenter-structurefactor");
  bool PairAmplitudeFlag = Manager.GetBoolean("pair-amplitude");
  bool EnergyExpectationFlag = Manager.GetBoolean("energy-expectation");
  bool ChordFlag = Manager.GetBoolean("chord");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  bool Statistics = true;
  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHESphereFermionsCorrelation requires a state" << endl;
      return -1;
    }

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("eigenstate") << endl;
      return -1;
    }

  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
  if (HaldaneBasisFlag == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 62)
#else
	if (LzMax <= 30)
#endif
	  if ((SymmetrizedBasis == false) || (TotalLz != 0))
	    {
	      Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	    }
	  else
	    {
	      if (Manager.GetString("load-hilbert") != 0)
		Space = new FermionOnSphereSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
	      else
		Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
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
		if ((SymmetrizedBasis == false) || (TotalLz != 0))
		  Space = new FermionOnSphereLong(NbrParticles, TotalLz, LzMax, MemorySpace);
		else
		  {
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new FermionOnSphereSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
	      }
	    else
	      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax, MemorySpace);
    }
  else //Haldane basis
    {
      int* ReferenceState = 0;
      if (Manager.GetString("reference-file") == 0)
	{
	  return -1;
	}
      else
	{
	  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles,LzMax, ReferenceState) == false)
	    return -1;
	}
      if (Manager.GetBoolean("p-truncated") == true)
	{
	    Space = new FermionOnSpherePTruncated(NbrParticles, TotalLz, LzMax, Manager.GetInteger("p-truncation"), ReferenceState);
        }
      else
       { 
         if (SymmetrizedBasis == false)
	   {
	     if (Manager.GetBoolean("large-basis") == true)
	       {
	         if (Manager.GetString("load-hilbert") != 0)
		   Space = new FermionOnSphereHaldaneLargeBasis(Manager.GetString("load-hilbert"), Manager.GetInteger("large-memory") << 10);
	         else
		   {
		      Space = new FermionOnSphereHaldaneLargeBasis(NbrParticles, TotalLz, LzMax, ReferenceState, Manager.GetInteger("large-memory") << 10);	  
		      if (Manager.GetString("save-hilbert") != 0)
		       {
		         ((FermionOnSphereHaldaneLargeBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		         return 0;
		       }
		   }
	       }
	     else
	       {
	         if (Manager.GetBoolean("huge-basis") == true)
		   {
		     if (Manager.GetString("save-hilbert") != 0)
		       {
		         Space = new FermionOnSphereHaldaneHugeBasis (NbrParticles, TotalLz, LzMax, Manager.GetInteger("file-size"), ReferenceState, ((unsigned long) Manager.GetInteger("huge-memory")) << 20, false);
		         ((FermionOnSphereHaldaneHugeBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		         return 0;
		       }
		     if (Manager.GetString("load-hilbert") == 0)
		       {
		         cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		         return -1;
		       }
		     Space = new FermionOnSphereHaldaneHugeBasis (Manager.GetString("load-hilbert"), Manager.GetInteger("huge-memory"));
		   }
	         else
		   {
#ifdef __64_BITS__
		     if (LzMax <= 62)
#else
		        if (LzMax <= 30)
#endif
		          {
			    if (Manager.GetString("load-hilbert") != 0)
			      Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
			    else
			      Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
			    if (Manager.GetString("save-hilbert") != 0)
			     {
			       ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
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
			        if (Manager.GetString("load-hilbert") != 0)
			          Space = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
			        else
			          Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
			        if (Manager.GetString("save-hilbert") != 0)
			         {
				   ((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
				   return 0;
			         }
			      }	       
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
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new FermionOnSphereHaldaneSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
		        ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
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
		         if (Manager.GetString("load-hilbert") != 0)
		           Space = new FermionOnSphereHaldaneSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		         else
		           Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
		         if (Manager.GetString("save-hilbert") != 0)
		          {
			    ((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			    return 0;
		          }
		      }
	    }
      }
   }
  cout << Space->GetHilbertSpaceDimension() << endl;

  if (GuidingCenterStructureFactorFlag == true)
   {
     cout<<"Evaluating the guiding-center structure factor: " << endl;
     cout<<"S_L = <rho_{L,0}rho_{L,0}>, where rho_{L,M}=\\sum_m \\sqrt{2L+1}C_{m0m}^{SLS} c_m^+ c_m" << endl;
     cout<<"Normalization 1/N_{orb}, for L->infty S->nu-1/nu^2" << endl;
     ofstream File;
     File.precision(14);
     if (Manager.GetString("output-file") != 0)
       File.open(Manager.GetString("output-file"), ios::binary | ios::out);
     else
      {
        cout << "Enter output file! " << endl;
        exit(1);
      }
  
     RealVector State;
     if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
      {
        if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	    return -1;      
	  }
	if (Space->GetLargeHilbertSpaceDimension()!=State.GetLargeVectorDimension())
	  {
	    cout << "Dimension mismatch between state and Hilbert space!"<<endl;
	    return -1;
	  }

        Complex* DensityMatEl;
        DensityMatEl = new Complex[LzMax + 1];

	//cout<<"Density: ";
        for (int i = 0; i <= LzMax; ++i)
         {
           ParticleOnSphereDensityOperator Operator (Space, i);
           DensityMatEl[i] = Operator.MatrixElement(State, State);
	   //cout<<"i= "<<i<<" "<<DensityMatEl[i]<<" ";
         }
        cout<<endl;

       OperatorMatrixElementOperation* Operation;
       ParticleOnSphereDensityDensityOperator* DensityDensityOperator;

	//double S = 0.5 * (double)LzMax;
       ClebschGordanCoefficients CoeffLLS(LzMax, LzMax);  
       for (int L = 0; L <= LzMax; ++L)
        {
          Complex SumIJ(0.0,0.0);
          for (int i = 0; i <= LzMax; ++i)
           {
            for (int j = 0; j <= LzMax; ++j)
	     {        
               double Factor = (LzMax + 1.0) * pow(-1.0, 2*LzMax - i - j) * CoeffLLS.GetCoefficient(((j << 1) - LzMax), -((j << 1) - LzMax), L << 1) * CoeffLLS.GetCoefficient(((i << 1) - LzMax), -((i << 1) - LzMax), L << 1);

               DensityDensityOperator = new ParticleOnSphereDensityDensityOperator(Space, i, j, i, j);
               Operation = new OperatorMatrixElementOperation (DensityDensityOperator, State, State);          
               Operation->ApplyOperation(Architecture.GetArchitecture());               

               SumIJ +=  Factor * Operation->GetScalar();
               if (i != j)
                 SumIJ += Factor * DensityMatEl[j]; 

               delete DensityDensityOperator;
               delete Operation;
	     }
           }     
         cout << L <<" "<< -(1.0/(double)(LzMax + 1)) * SumIJ.Re << " " << -(1.0/(double)(LzMax + 1)) * SumIJ.Im <<endl;
         File << L <<" "<< -(1.0/(double)(LzMax + 1)) * SumIJ.Re << " " << -(1.0/(double)(LzMax + 1)) * SumIJ.Im <<endl;
        }
       delete[] DensityMatEl;
       return 0;
     } 
   }
  else if (StructureFactorFlag == true)
   {
     cout<<"Projected structure factor is evaluated for LLL; " << endl;
     cout<<"use normalization convention from He, Simon & Halperin." << endl;
     ofstream File;
     File.precision(14);
     if (Manager.GetString("output-file") != 0)
       File.open(Manager.GetString("output-file"), ios::binary | ios::out);
     else
      {
        cout << "Enter output file! " << endl;
        exit(1);
      }
  
     RealVector State;
     if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
      {
        if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	    return -1;      
	  }
	if (Space->GetLargeHilbertSpaceDimension()!=State.GetLargeVectorDimension())
	  {
	    cout << "Dimension mismatch between state and Hilbert space!"<<endl;
	    return -1;
	  }

        Complex* DensityMatEl;
        DensityMatEl = new Complex[LzMax + 1];

        for (int i = 0; i <= LzMax; ++i)
         {
           ParticleOnSphereDensityOperator Operator (Space, i);
           DensityMatEl[i] = Operator.MatrixElement(State, State);
         }

       double S = 0.5 * (double)LzMax;
       ClebschGordanCoefficients CoeffLLS(LzMax, LzMax);  
       for (int L = 0; L <= LzMax; ++L)
        {
          double Factor = pow(-1.0, 3.0 * LzMax + (L << 1)) * pow(LzMax + 1.0, 2.0) /(4.0 * M_PI * (2.0 * L + 1.0));
          Factor *=  pow(CoeffLLS.GetCoefficient(-LzMax, LzMax, L << 1), 2.0);

          Complex SumIJ(0.0,0.0);
          for (int i = 0; i <= LzMax; ++i)
           {
            double FactorI = Factor * CoeffLLS.GetCoefficient(((i << 1) - LzMax), -((i << 1) - LzMax), L << 1);
            for (int j = 0; j <= LzMax; ++j)
	     {        
               double FactorJ = FactorI * CoeffLLS.GetCoefficient(((j << 1) - LzMax), -((j << 1) - LzMax), L << 1) * pow(-1.0, -(i - S)-(j-S));
               ParticleOnSphereDensityDensityOperator Operator (Space, i, j, i, j);
               SumIJ -=  FactorJ * Operator.MatrixElement(State, State);
               if (i == j)
                 SumIJ += FactorJ * DensityMatEl[i]; 
                
	     }
           }     
         cout << L <<" "<< (4.0 * M_PI/(double)NbrParticles) * SumIJ.Re << " " << (4.0 * M_PI/(double)NbrParticles) * SumIJ.Im <<endl;
         File << L <<" "<< (4.0 * M_PI/(double)NbrParticles) * SumIJ.Re << " " << (4.0 * M_PI/(double)NbrParticles) * SumIJ.Im <<endl;
        }
       delete[] DensityMatEl;
       return 0;
     } 
   }

   if (PairAmplitudeFlag == true)
   {
     cout<<"Evaluating the pair amplitudes <xi_{J,M}^+ xi_{J,M}>" << endl;
     cout<<"assuming M can be set to zero"<<endl;
     cout<<"xi_{J,M}=\\sum_{m1,m2} C_{m2m1M}^{SSJ} c_m1 c_m2 "<<endl;
     cout<<"Total energy = \\sum_{J,M} V_J <xi_{J,M}^+ xi_{J,M}>"<<endl;
     ofstream File;
     File.precision(14);
     if (Manager.GetString("output-file") != 0)
       File.open(Manager.GetString("output-file"), ios::binary | ios::out);
     else
      {
        cout << "Enter output file! " << endl;
        exit(1);
      }
  
     RealVector State;
     if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
      {
        if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	    return -1;      
	  }

       ClebschGordanCoefficients CoeffLLS(LzMax, LzMax);
       Complex TotalEnergy (0.0, 0.0);
       double* PseudoPotentials = 0;
       if (EnergyExpectationFlag)
        {
          ConfigurationParser InteractionDefinition;
          if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
 	   {
	     InteractionDefinition.DumpErrors(cout) << endl;
 	     return -1;
 	   }
          int TmpNbrPseudoPotentials;
          if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == false)
 	   {
	     cout << "Weights is not defined or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	     return -1;
	   }
          cout << "LzMax= " << TmpNbrPseudoPotentials << " " << LzMax << endl;
          if (TmpNbrPseudoPotentials != (LzMax +1))
 	   {
	     cout << "Invalid number of pseudo-potentials" << endl;
	     return -1;	  
	   }
         }

       OperatorMatrixElementOperation* Operation;
       ParticleOnSphereDensityDensityOperator* DensityDensityOperator;

       for (int L = 0; L <= LzMax; ++L)
        {
          Complex PairAmpL(0.0,0.0);
          for (int m1 = 0; m1 <= LzMax; ++m1)
           {
            for (int m3 = 0; m3 <= LzMax; ++m3)
	     {        
               double Factor = CoeffLLS.GetCoefficient(((m1 << 1) - LzMax), -((m1 << 1) - LzMax), 2*LzMax - (L << 1)) * CoeffLLS.GetCoefficient(-((m3 << 1) - LzMax), ((m3 << 1) - LzMax), 2*LzMax - (L << 1));

               DensityDensityOperator = new ParticleOnSphereDensityDensityOperator(Space, m1, LzMax - m1, m3, LzMax - m3);
               Operation = new OperatorMatrixElementOperation (DensityDensityOperator, State, State);          
               Operation->ApplyOperation(Architecture.GetArchitecture());               

               PairAmpL += Factor * Operation->GetScalar(); 

               delete DensityDensityOperator;
               delete Operation;
	     }
           }     
         cout << L <<" "<< PairAmpL.Re << " " << PairAmpL.Im <<endl;
         File << L <<" "<< PairAmpL.Re << " " << PairAmpL.Im <<endl;
         if (EnergyExpectationFlag)
   	   TotalEnergy += (2*(LzMax-L)+1) * PairAmpL * PseudoPotentials[L];
        }

      if (EnergyExpectationFlag)
 	cout << "Twice total energy: " <<TotalEnergy<<endl;

       delete [] PseudoPotentials;
       return 0;
     } 
   }

  int TwoQ = LzMax - 2 * LandauLevel;

  AbstractFunctionBasis* Basis;
  if (LandauLevel == 0)
    Basis = new ParticleOnSphereFunctionBasis(TwoQ);
  else
    Basis = new ParticleOnSphereGenericLLFunctionBasis(TwoQ, LandauLevel);

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);


  Complex* PrecalculatedMatrixElements = new Complex [LzMax + 1];
  Complex* PrecalculatedValues = new Complex [LzMax + 1];
  RealVector State;
  if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
    {
      if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (Space->GetLargeHilbertSpaceDimension()!=State.GetLargeVectorDimension())
	  {
	    cout << "Dimension mismatch between state and Hilbert space!"<<endl;
	    return -1;
	  }
      
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax - LandauLevel);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax - LandauLevel, i, LzMax - LandauLevel);
	    PrecalculatedMatrixElements[i] = Operator.MatrixElement(State, State);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedMatrixElements[i] = Operator.MatrixElement(State, State);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State);
	  }
    }
  else
    {
      ComplexVector ComplexState;
      if (ComplexState.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (Space->GetLargeHilbertSpaceDimension()!= ComplexState.GetLargeVectorDimension())
	  {
	    cout << "Dimension mismatch between state and Hilbert space!"<<endl;
	    return -1;
	  }
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax - LandauLevel);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax - LandauLevel, i, LzMax - LandauLevel);
	    PrecalculatedMatrixElements[i] = Operator.MatrixElement(ComplexState, ComplexState);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedMatrixElements[i] = Operator.MatrixElement(ComplexState, ComplexState);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState);
	  }
    }

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("eigenstate") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  if (DensityFlag == true)      
    File << "# density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  else
    File << "# density-density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  File << "#" << endl << "# (l+S)    n_l" << endl;
  if (CoefficientOnlyFlag == false)
    {
      for (int i = 0; i <= LzMax; ++i)
	File << "# " << i << " " << PrecalculatedMatrixElements[i].Re<< endl;
    }
  else
    {
      for (int i = 0; i <= LzMax; ++i)
	File << i << " " << PrecalculatedMatrixElements[i].Re<< endl;
    }
  if (CoefficientOnlyFlag == false)
    {
      double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      if (DensityFlag == true)
	Factor1 = 1.0;//4.0 * M_PI;
      double Factor2;
      if (Manager.GetBoolean("radians") == true)
	Factor2 = 1.0;
      else
	Factor2 = sqrt (0.5 * TwoQ);
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  Sum = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis->GetFunctionValue(Value, TmpValue, i);
	      Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	    }
	  if (ChordFlag == false)
	    File << (X * Factor2) << " " << (Norm(Sum)  * Factor1) << endl;
	  else
	    File << (2.0 * Factor2 * sin (X * 0.5)) << " " << Norm(Sum)  * Factor1 << endl;
	  X += XInc;
	}
    }
  File.close();
 
  delete[] PrecalculatedValues;
  delete[] PrecalculatedMatrixElements;

  return 0;
}


