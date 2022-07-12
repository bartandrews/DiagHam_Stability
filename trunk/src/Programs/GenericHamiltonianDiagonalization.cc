#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/FileBasedHamiltonian.h"
#include "Hamiltonian/FileBasedHermitianHamiltonian.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("GenericHamiltonianDiagonalization" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('\n', "hamiltonian", "text file where the hamiltonian matrix elements are stored");
  (*SystemGroup) += new  SingleStringOption ('\n', "bin-hamiltonian", "use a binary encoded hamiltonian");
  (*SystemGroup) += new  SingleStringOption ('\n', "multiple-binhamiltonians", "use a series of binary encoded hamiltonians with weights (if defined, the third column is a complex phase)");
  (*SystemGroup) += new  BooleanOption ('\n', "fortran", "assume indices are 1-based instead of 0-based (i.e. fortran index convention)");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "skip-lines", "skip the first n-tf lines of the input file", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "data-column", "index of the column that contains the matrix elements (or their real part)", 0);
  (*SystemGroup) += new BooleanOption  ('c', "complex", "indicate that the Hamiltonian is complex");
  (*SystemGroup) += new SingleDoubleOption ('\n', "shift-spectrum", "shift the spectrum by a constant value during the diagonalization", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "rescaling-factor", "apply a rescaling factor to the hamiltonian", 1.0);
  (*SystemGroup) += new SingleStringOption ('\n', "add-diagonal", "add a diagonal contribution to the hamiltonian");
  (*SystemGroup) += new SingleIntegerOption ('\n', "column-diagonal", "indicates which column has to be used in --add-diagonal", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");  
  (*SystemGroup) += new BooleanOption  ('s', "nosort-rowindices", "do not sort the row indices, assuming they are already sorted from the smallest to the largest"); 
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
 (*OutputGroup) += new SingleStringOption ('o', "output-file", "prefix to use for output file names", "dummy");
  (*OutputGroup) += new SingleStringOption ('\n', "eigenstate-file", "prefix to use for the eigenstate output file names", "dummy");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "friendlyshowhamiltonian-error", "when using --friendlyshow-hamiltonian, value below which a matrix element is considered equal to zero", 0.0);
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericHamiltonianDiagonalization -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("hamiltonian") == 0) && (Manager.GetString("bin-hamiltonian") == 0) && (Manager.GetString("multiple-binhamiltonians") == 0))
    {
      cout << "no hamiltonian provided" << endl; 
      return -1;
    }
  if ((Manager.GetString("hamiltonian") != 0) && (IsFile(Manager.GetString("hamiltonian")) == false))
    {
      cout << "can't open hamiltonian " << Manager.GetString("hamiltonian") << endl; 
      return -1;
    }

  if ((Manager.GetString("bin-hamiltonian") != 0) && (IsFile(Manager.GetString("bin-hamiltonian")) == false))
    {
      cout << "can't open hamiltonian " << Manager.GetString("bin-hamiltonian") << endl; 
      return -1;
    }
  if ((Manager.GetString("multiple-binhamiltonians") != 0) && (IsFile(Manager.GetString("multiple-binhamiltonians")) == false))
    {
      cout << "can't open hamiltonian lst " << Manager.GetString("multiple-binhamiltonians") << endl; 
      return -1;
    }

  AbstractHamiltonian* Hamiltonian = 0;
  bool ComplexFlag = Manager.GetBoolean("complex");
  if (Manager.GetString("hamiltonian") != 0)
    {
      if (ComplexFlag == false)
	{
	  Hamiltonian  = new FileBasedHamiltonian(Manager.GetString("hamiltonian"), Manager.GetInteger("data-column"), false, Manager.GetBoolean("fortran"), Manager.GetInteger("skip-lines"), Manager.GetBoolean("nosort-rowindices"));
	}
      else
	{
	  Hamiltonian  = new FileBasedHermitianHamiltonian(Manager.GetString("hamiltonian"), Manager.GetInteger("data-column"), false, Manager.GetBoolean("fortran"), Manager.GetInteger("skip-lines"), Manager.GetBoolean("nosort-rowindices"));
	}
      
      Architecture.GetArchitecture()->SetDimension(Hamiltonian->GetHilbertSpaceDimension());	
      
      char* CommentLine = new char [strlen(Manager.GetString("hamiltonian")) + 256];
      sprintf (CommentLine, "eigenvalues of %s\n#", Manager.GetString("hamiltonian"));
      
      char* EigenvectorFileName = Manager.GetString("eigenstate-file");
      
      Hamiltonian->ShiftHamiltonian(Manager.GetDouble("shift-spectrum"));
      
      if (ComplexFlag == false)
	{
	  GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, "", CommentLine, Manager.GetDouble("shift-spectrum"),  Manager.GetString("output-file"), true, EigenvectorFileName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	}
      else
	{
	  Lanczos.SetComplexAlgorithms();
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, "", CommentLine, Manager.GetDouble("shift-spectrum"),  Manager.GetString("output-file"), true, EigenvectorFileName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	}
    }
  else
    {
      UndescribedHilbertSpace* DummyHilbertSpace = 0;
      RealSymmetricMatrix HRepReal;
      HermitianMatrix HRep;
      if (ComplexFlag == false)
	{
	  if (Manager.GetString("bin-hamiltonian") != 0)
	    {
	      if (HRepReal.ReadMatrix(Manager.GetString("bin-hamiltonian")) == false)
		{
		  cout << "can't read " << Manager.GetString("bin-hamiltonian") << endl;
		  return -1;
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile HamiltonianFile;
	      if (HamiltonianFile.Parse(Manager.GetString("multiple-binhamiltonians")) == false)
		{
		  return -1;
		}
	      if (HRepReal.ReadMatrix(HamiltonianFile(0, 0)) == false)
		{
		  cout << "can't read " << HamiltonianFile(0, 0) << endl;
		  return -1;
		}
	      double* TmpWeigths;
	      if (HamiltonianFile.GetNbrColumns() > 1)
		{
		  TmpWeigths = HamiltonianFile.GetAsDoubleArray(1);
		}
	      else
		{
		  TmpWeigths = new double[HamiltonianFile.GetNbrLines()];
		  for (int i = 0; i < HamiltonianFile.GetNbrLines(); ++i)
		    {
		      TmpWeigths[i] = 1.0;
		    }
		}	      
	      HRepReal *= TmpWeigths[0];
	      double* TmpPhases = 0;
	      if (HamiltonianFile.GetNbrColumns() > 2)
		{
		  TmpPhases = HamiltonianFile.GetAsDoubleArray(2);
		  HRep = HermitianMatrix (HRepReal, TmpPhases[0]);		  
		  ComplexFlag = true;
		}
	      for (int i = 1; i < HamiltonianFile.GetNbrLines(); ++i)
		{
		  RealSymmetricMatrix TmpH;
		  if (TmpH.ReadMatrix(HamiltonianFile(0, i)) == false)
		    {
		      cout << "can't read " << HamiltonianFile(0, i) << endl;
		      return -1;
		    }
		  if (TmpH.GetNbrRow() != HRepReal.GetNbrRow())
		    {
		      cout << "error, " << HamiltonianFile(0, i) << " and " << HamiltonianFile(0, 0) 
			   << " do not have the same size" << endl;
		    }
		  TmpH *= TmpWeigths[i];
		  if (TmpPhases == 0)
		    HRepReal += TmpH;
		  else
		    {
		      HermitianMatrix TmpHermitianH (TmpH, TmpPhases[i]);
		      HRep += TmpHermitianH;
		    }
		}
	    }
	  if (Manager.GetDouble("rescaling-factor") != 1.0)
	    {
	      if (ComplexFlag == false)
		HRepReal *= Manager.GetDouble("rescaling-factor");
	      else
		HRep *= Manager.GetDouble("rescaling-factor");
	    }
	  if (Manager.GetString("add-diagonal") != 0)
	    {
	      MultiColumnASCIIFile DiagonalFile;
	      if (DiagonalFile.Parse(Manager.GetString("add-diagonal")) == false)
		{
		  DiagonalFile.DumpErrors(cout);
		  return -1;
		}
	      if (DiagonalFile.GetNbrLines() != HRepReal.GetNbrRow())
		{
		  cout << "error, " << Manager.GetString("add-diagonal") << " does not have the coorect number of lines (is " << DiagonalFile.GetNbrLines() 
		       << ", should be " << HRepReal.GetNbrRow() << ")" << endl;
		}
	      double* TmpDiagonalElements = DiagonalFile.GetAsDoubleArray(Manager.GetInteger("column-diagonal"));
	      for (int i = 0; i < HRepReal.GetNbrRow(); ++i)
		{
		  if (ComplexFlag == false)
		    HRepReal.AddToMatrixElement(i, i, TmpDiagonalElements[i]);
		  else
		    HRep.AddToMatrixElement(i, i, TmpDiagonalElements[i]);
		}
	    }
	  DummyHilbertSpace = new UndescribedHilbertSpace(HRepReal.GetNbrRow());
	  if (ComplexFlag == false)
	    Hamiltonian  = new ExplicitHamiltonian(DummyHilbertSpace, &HRepReal);
	  else
	    Hamiltonian = new ExplicitHamiltonian(DummyHilbertSpace, &HRep);
	}
      else
	{
	  if (Manager.GetString("bin-hamiltonian") != 0)
	    {
	      if (HRep.ReadMatrix(Manager.GetString("bin-hamiltonian")) == false)
		{
		  cout << "can't read " << Manager.GetString("bin-hamiltonian") << endl;
		  return -1;
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile HamiltonianFile;
	      if (HamiltonianFile.Parse(Manager.GetString("multiple-binhamiltonians")) == false)
		{
		  return -1;
		}
	      if (HRep.ReadMatrix(HamiltonianFile(0, 0)) == false)
		{
		  cout << "can't read " << HamiltonianFile(0, 0) << endl;
		  return -1;
		}
	      double* TmpWeigths;
	      if (HamiltonianFile.GetNbrColumns() > 1)
		{
		  TmpWeigths = HamiltonianFile.GetAsDoubleArray(1);
		}
	      else
		{
		  TmpWeigths = new double[HamiltonianFile.GetNbrLines()];
		  for (int i = 0; i < HamiltonianFile.GetNbrLines(); ++i)
		    {
		      TmpWeigths[i] = 1.0;
		    }
		}	      
	      HRep *= TmpWeigths[0];
	      for (int i = 1; i < HamiltonianFile.GetNbrLines(); ++i)
		{
		  HermitianMatrix TmpH;
		  if (TmpH.ReadMatrix(HamiltonianFile(0, i)) == false)
		    {
		      cout << "can't read " << HamiltonianFile(0, i) << endl;
		      return -1;
		    }
		  if (TmpH.GetNbrRow() != HRep.GetNbrRow())
		    {
		      cout << "error, " << HamiltonianFile(0, i) << " and " << HamiltonianFile(0, 0) 
			   << " do not have the same size" << endl;
		    }
		  TmpH *= TmpWeigths[i];
		  HRep += TmpH;
		}
	    }
	  if (Manager.GetDouble("rescaling-factor") != 1.0)
	    {
	      HRep *= Manager.GetDouble("rescaling-factor");
	    }
	  if (Manager.GetString("add-diagonal") != 0)
	    {
	      MultiColumnASCIIFile DiagonalFile;
	      if (DiagonalFile.Parse(Manager.GetString("add-diagonal")) == false)
		{
		  DiagonalFile.DumpErrors(cout);
		  return -1;
		}
	      if (DiagonalFile.GetNbrLines() != HRep.GetNbrRow())
		{
		  cout << "error, " << Manager.GetString("add-diagonal") << " does not have the coorect number of lines (is " << DiagonalFile.GetNbrLines() 
		       << ", should be " << HRep.GetNbrRow() << ")" << endl;
		}
	      double* TmpDiagonalElements = DiagonalFile.GetAsDoubleArray(Manager.GetInteger("column-diagonal"));
	      for (int i = 0; i < HRep.GetNbrRow(); ++i)
		{
		  HRep.AddToMatrixElement(i, i, TmpDiagonalElements[i]);
		}
	    }
	  DummyHilbertSpace = new UndescribedHilbertSpace(HRep.GetNbrRow());
	  Hamiltonian = new ExplicitHamiltonian(DummyHilbertSpace, &HRep);
	}
      Architecture.GetArchitecture()->SetDimension(Hamiltonian->GetHilbertSpaceDimension());	
      
      char* CommentLine;
      if (Manager.GetString("bin-hamiltonian") != 0)
	{
	  CommentLine = new char [strlen(Manager.GetString("bin-hamiltonian")) + 256];
	  sprintf (CommentLine, "eigenvalues of %s\n#", Manager.GetString("bin-hamiltonian"));
	}
      else
	{
	  CommentLine = new char [strlen(Manager.GetString("multiple-binhamiltonians")) + 256];
	  sprintf (CommentLine, "eigenvalues of %s\n#", Manager.GetString("multiple-binhamiltonians"));
	}
      
      char* EigenvectorFileName = Manager.GetString("eigenstate-file");

      if (ComplexFlag == false)
	{
	  GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, "", CommentLine, Manager.GetDouble("shift-spectrum"),  Manager.GetString("output-file"), true, EigenvectorFileName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	}
      else
	{
	  Lanczos.SetComplexAlgorithms();
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, "", CommentLine, Manager.GetDouble("shift-spectrum"),  Manager.GetString("output-file"), true, EigenvectorFileName);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	}
    }
  return 0;
}
