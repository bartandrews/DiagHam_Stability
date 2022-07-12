#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <cstring>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("MatrixElement" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of the vector files obtained using exact diagonalization (4 files required)");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('g', "gauge", "Take phase such that largest coefficient is real and positive");
  (*SystemGroup) += new SingleStringOption  ('i', "interaction", "File defining the interaction matrix elements");
  (*SystemGroup) += new BooleanOption  ('\n', "quiet", "discard any output except the matrix element (suitable for scripting)");
  (*SystemGroup) += new BooleanOption  ('\n', "verbose", "give a lot of output");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool QuietFlag = Manager.GetBoolean("quiet");
  bool VerboseFlag = Manager.GetBoolean("verbose");

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  if (NbrVectors!=4)
    {
      cout << "Precisely four vector are required for a two-body interaction!"<<endl;
      exit(1);
    }

  if (QuietFlag == false)
    for (int i=0; i<NbrVectors; ++i)    
      cout << "File "<<i<<"  "<<VectorFiles[i]<<endl;

  // read configuration file
  if (Manager.GetString("interaction")==NULL)
    {
      cout << "A definition of the interaction is required!"<<endl;
      exit(1);
    }
  if (VerboseFlag)
    cout << "Parsing interaction definition"<<endl;
  ConfigurationParser LatticeDefinition;
  if (LatticeDefinition.Parse(Manager.GetString("interaction")) == false)
    {
      LatticeDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  int HilbertSpaceDimension;
  if ((LatticeDefinition.GetAsSingleInteger("HilbertSpaceDimension", HilbertSpaceDimension) == false) || (HilbertSpaceDimension <= 0))
    {
      cout << "HilbertSpaceDimension is not defined or has invalid value" << endl;
      exit(-1);
    }
  int NbrMatrixElements;
  if ((LatticeDefinition.GetAsSingleInteger("NbrMatrixElements", NbrMatrixElements) == false) || (NbrMatrixElements <= 0))
    {
      cout << "NbrMatrixElements is not defined or has invalid value" << endl;
      exit(-1);
    }
  int TmpLength;
  double *MatrixElements;
  int *DiagonalEntries;
	
  if (LatticeDefinition.GetAsDoubleArray("MatrixElements",' ',MatrixElements, TmpLength) == false)
    {
      cout << "MatrixElements is not defined or has invalid value" << endl;
      exit(-1);
    }
  if ((LatticeDefinition["Sparse"]!=NULL)&&
      ( (strcmp(LatticeDefinition["Sparse"],"yes")==0) || (strcmp(LatticeDefinition["Sparse"],"YES")==0)
	|| (strcmp(LatticeDefinition["Sparse"],"true")==0) || (strcmp(LatticeDefinition["Sparse"],"TRUE")==0) ))
    {
      cout << "Sparse matrix not yet implemented"<<endl;
      exit(-1);
    }
  else
    {
      if (LatticeDefinition.GetAsIntegerArray("DiagonalEntries",' ',DiagonalEntries, TmpLength) == false)
	{
	  cout << "DiagonalEntries is not defined or has invalid value" << endl;
	  exit(-1);
	}
      if (TmpLength!=HilbertSpaceDimension)
	{
	  cout << "Number of diagonal matrix entries has to equal Hilbert-space dimension"<<endl;
	  exit(-1);
	}
      for (int i=0; i<HilbertSpaceDimension; ++i)
	if ((DiagonalEntries[i]<0) || (DiagonalEntries[i]>=NbrMatrixElements))
	  {
	    cout << "Entry no. "<<i<<" in DiagonalEntries is out of bounds" << endl;
	    exit(-1);
	  }
      if ((LatticeDefinition["HaveOffDiagonal"]!=NULL)&&
	  ( (strcmp(LatticeDefinition["HaveOffDiagonal"],"yes")==0) || (strcmp(LatticeDefinition["HaveOffDiagonal"],"YES")==0)
	    || (strcmp(LatticeDefinition["HaveOffDiagonal"],"true")==0) || (strcmp(LatticeDefinition["HaveOffDiagonal"],"TRUE")==0) ))
	{
	  cout << "OffDiagonal matrix elements not yet implementated"<<endl;
	  exit(-1);
	}
    }
  
  if (Manager.GetBoolean("complex"))
    {      
      ComplexVector* States= new ComplexVector[4];
      Complex *Results=new Complex[NbrMatrixElements];
      for (int i=0; i<NbrMatrixElements; ++i)
	Results[i]=0.0;
      for (int i=0; i<4; ++i)
	{
	  if (States[i].ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  if (States[i].GetVectorDimension()!=HilbertSpaceDimension)
	    {
	      cout << "Vector size and dimension defined in interaction file have to coindice"<<endl;
	      exit(-1);
	    }
	  if (Manager.GetBoolean("gauge"))
	    {
	      int maxI=0;
	      double maxN=0.0;
	      for (int n=0; n<HilbertSpaceDimension; ++n)
		if (Norm(States[i][n])>maxN)
		  {
		    maxI=n;
		    maxN=Norm(States[i][n]);
		  }
	      Complex Phase=Polar(1.0,-Arg(States[i][maxI]));
	      if (VerboseFlag)
		cout << "State "<<i<<": gauging with phase: "<<Phase<<" for element "<<maxI<<endl;
	      States[i]*=Phase;
	    }
	}
      for (int i=0; i<HilbertSpaceDimension; ++i)
	{
	  if (VerboseFlag)
	    cout << "i="<<i<<" DiagonalEntry="<<DiagonalEntries[HilbertSpaceDimension-1-i] <<" term "<<Conj(States[0][i])*Conj(States[1][i])*States[2][i]*States[3][i]<<" t_0="<<States[0][i]<<endl;
	  Results[DiagonalEntries[HilbertSpaceDimension-1-i]] += Conj(States[0][i])*Conj(States[1][i])*States[2][i]*States[3][i];
	}

      if (QuietFlag)
	{
	  cout <<Results[0].Re<<" "<<Results[0].Im<<endl;
	  for (int i=1; i<NbrMatrixElements; ++i)
	    cout <<Results[i].Re<<" "<<Results[i].Im<<endl;
	}
      else
	{
	  Complex Result=0.0;
	  for (int i=0; i<NbrMatrixElements; ++i)
	    Result+=MatrixElements[i]*Results[i];
	  if (fabs(Result.Im)>1e-12)
	    cout << "M="<<Norm(Result)<<"*Exp["<<Arg(Result)/M_PI<<"*I*Pi]"<<endl;
	  else
	    cout << "M="<<Result<<endl;
	  cout << "General M="<<Results[0]<<"*m(0)";
	  for (int i=1; i<NbrMatrixElements; ++i)
	    cout << " + "<<Results[i]<<"*m("<<i<<")";
	  cout << endl;
	}
      delete [] Results;
    }
  else // real vectors
    {
      RealVector* States= new RealVector[4];
      double *Results=new double[NbrMatrixElements];
      for (int i=0; i<NbrMatrixElements; ++i)
	Results[i]=0.0;
      for (int i=0; i<4; ++i)
	{
	  if (States[i].ReadVector (VectorFiles[i]) == false)
	    {
	      cout << "can't open vector file " << VectorFiles[i] << endl;
	      return -1;      
	    }
	  if (States[i].GetVectorDimension()!=HilbertSpaceDimension)
	    {
	      cout << "Vector size and dimension defined in interaction file have to coindice"<<endl;
	      exit(-1);
	    }
	  for (int i=0; i<HilbertSpaceDimension; ++i)
	    Results[DiagonalEntries[i]] += States[0][i]*States[1][i]*States[2][i]*States[3][i];
	}
      double Result=0.0;
      for (int i=0; i<NbrMatrixElements; ++i)
	Result+=MatrixElements[i]*Results[i];
      if (QuietFlag)
	{
	  cout <<Results[0]<<endl;
	  for (int i=1; i<NbrMatrixElements; ++i)
	    cout <<Results[i]<<endl;
	}
      else
	{
	  cout << "M="<<Result<<endl;
	  cout << "General M="<<Results[0]<<"*m(0)";
	  for (int i=1; i<NbrMatrixElements; ++i)
	    cout << " + "<<Results[i]<<"*m("<<i<<")";
	  cout << endl;
	}
      delete [] Results;
    }
}
