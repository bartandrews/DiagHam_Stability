#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("BuildSuperPosition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of all vector files that should be superposed");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new SingleDoubleOption  ('r', "random-component", "amplitude of a random component to be added",0.5);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "names of output filename","superposition.vec");
  (*SystemGroup) += new SingleStringOption  ('f', "description-file", "build the superposition for a linear combination of a vector set described in a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  if (Manager.GetString("description-file") == 0)
    {
      int NbrVectors, ValidVectors=0, VectorDimension=0;
      char** VectorFiles = Manager.GetStrings("states",NbrVectors);
      
      
      Vector **Vectors = new Vector*[NbrVectors], *Result, *Random;
      
      bool tmpB, haveVector=false;
      
      if (Manager.GetBoolean("complex"))
	{
	  Result = new ComplexVector();
	  for (int i=0; i<NbrVectors; ++i)
	    {
	      Vectors[i]= new ComplexVector();
	      tmpB = ((ComplexVector*)Vectors[i])->ReadVector(VectorFiles[i]);
	      if (!haveVector)
		VectorDimension=((ComplexVector*)Vectors[i])->GetVectorDimension();
	      if (haveVector && (((ComplexVector*)Vectors[i])->GetVectorDimension()!=VectorDimension))
		{
		  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
		  exit(1);
		}
	      haveVector=haveVector | tmpB;
	      if (tmpB) ++ValidVectors;
	    }
	  if (ValidVectors==0)
	    {
	      cout <<"At least one valid vector is required!"<<endl;
	      exit(1);
	    }
	  Random = new ComplexVector(VectorDimension);
	  for (int i = 0; i < VectorDimension; ++i)
	    {
	      ((ComplexVector*)Random)->Re(i) = (rand() - 32767) * 0.5;
	      ((ComplexVector*)Random)->Im(i) = (rand() - 32767) * 0.5;
	    }
	  *((ComplexVector*)(Random)) /= Random->Norm();
	}
      else
	{
	  Result = new RealVector();
	  for (int i=0; i<NbrVectors; ++i)
	    {
	      Vectors[i] = new RealVector();
	      tmpB = ((RealVector*)Vectors[i])->ReadVector(VectorFiles[i]);
	      if (!haveVector)
		VectorDimension=((RealVector*)Vectors[i])->GetVectorDimension();
	      if (haveVector && (((RealVector*)Vectors[i])->GetVectorDimension()!=VectorDimension))
		{
		  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
		  exit(1);
		}
	      haveVector=haveVector | tmpB;
	      if (tmpB) ++ValidVectors;
	    }
	  if (ValidVectors==0)
	    {
	      cout <<"At least one valid vector is required!"<<endl;
	      exit(1);
	    }
	  Random = new RealVector(VectorDimension);
	  for (int i = 0; i < VectorDimension; ++i)
	    (*(RealVector*)(Random))[i] = (rand() - 32767) * 0.5;
	  *((RealVector*)(Random)) /= Random->Norm();
	}	
      
      double RandomComponent = Manager.GetDouble("random-component");  
      double GivenComponent = (1.0-RandomComponent)/NbrVectors;
      
      Result->Resize(VectorDimension);
      Result->ClearVector();
      
      for (int i=0; i<NbrVectors; ++i)
	{
	  cout<<"Adding Vector "<<i<<endl;
	  Result->AddLinearCombination(GivenComponent,*(Vectors[i]));
	}
      
      Result->AddLinearCombination(RandomComponent,*Random);
      
      char *OutputName = Manager.GetString("output");
      
      if (Manager.GetBoolean("complex"))
	{
	  *((ComplexVector*)(Result)) /= Result->Norm();
	  ((ComplexVector*)Result)->WriteVector(OutputName);
	}
      else
	{
	  *((RealVector*)(Result)) /= Result->Norm();
	  ((RealVector*)Result)->WriteVector(OutputName);
	}
    }
  else
    {
      MultiColumnASCIIFile Description;
      if (Description.Parse(Manager.GetString("description-file")) == false)
	{
	  Description.DumpErrors(cout);
	  return -1;
	}
      if (Description.GetNbrColumns() < 2)
	{
	  cout << "wrong number of columns in " << Manager.GetString("description-file") << endl;
	  return -1;
	}
      if (Description.GetNbrLines() < 2)
	{
	  cout << "not enough vectors in " << Manager.GetString("description-file") << endl;
	  return -1;
	}
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector Result;
	  if (Result.ReadVector(Description(0, 0)) == false)
	    {
	      cout << "can't open file " << Description(0, 0) << endl;
	      return -1;	      
	    }
	  double* Coefficients = Description.GetAsDoubleArray(1);
	  Result *= Coefficients[0];
 	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      ComplexVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}	      
	      Result.AddLinearCombination(Coefficients[i], TmpVector);
	    }
	  Result /= Result.Norm();
	  Result.WriteVector(Manager.GetString("output"));
	}
      else
	{
	  RealVector Result;
	  if (Result.ReadVector(Description(0, 0)) == false)
	    {
	      cout << "can't open file " << Description(0, 0) << endl;
	      return -1;	      
	    }
	  double* Coefficients = Description.GetAsDoubleArray(1);
	  Result *= Coefficients[0];
 	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      RealVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}	      
	      Result.AddLinearCombination(Coefficients[i], TmpVector);
	    }
	  Result /= Result.Norm();
	  Result.WriteVector(Manager.GetString("output"));
	}
    }
}
  
