#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RationalVector.h"
#include "Vector/LongRationalVector.h"

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
  (*SystemGroup) +=  new BooleanOption  ('\n', "rational", "input vectors are rational vectors");
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif
  (*SystemGroup) += new SingleDoubleOption  ('r', "random-component", "amplitude of a random component to be added",0.0);
  (*SystemGroup) += new SingleIntegerOption  ('R', "random-only", "generate a pure random vector, argument is dimension",0);
  (*SystemGroup) += new SingleStringOption  ('\n', "random-orthogonal", "build a pure random vector orthgonal to the one provided as an argument");
  (*SystemGroup) += new BooleanOption  ('n', "no-normalize", "do NOT normalize the final vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "names of output filename","superposition.vec");
  (*SystemGroup) += new SingleStringOption  ('f', "description-file", "build the superposition for a linear combination of a vector set described in a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  if (Manager.GetInteger("random-only") > 0)
    {
      int VectorDimension=Manager.GetInteger("random-only");
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector Vector(VectorDimension);
	  for (int i = 0; i < VectorDimension; ++i)
	    {
	      Vector.Re(i) = (rand() - 32767) * 0.5;
	      Vector.Im(i) = (rand() - 32767) * 0.5;
	    }
	  Vector /= Vector.Norm();
	  Vector.WriteVector(Manager.GetString("output"));
	}
      else
	{
	  RealVector Vector(VectorDimension);
	  for (int i = 0; i < VectorDimension; ++i)
	    Vector[i] = (rand() - 32767) * 0.5;
	  Vector /= Vector.Norm();
	  Vector.WriteVector(Manager.GetString("output"));
	}
      cout << "Generated random vector of dimension "<<VectorDimension<<" (file: "<<
	Manager.GetString("output")<<" )"<<endl;
      exit(0);
    }
  
  if (Manager.GetString("random-orthogonal") != 0)
    {
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector InputVector;
	  if (InputVector.ReadVector(Manager.GetString("random-orthogonal")) == false)
	    {
	      cout << "can't open " << Manager.GetString("random-orthogonal") << endl;
	      return -1;
	    }
	  ComplexVector Vector(InputVector.GetLargeVectorDimension());
	  for (long i = 0l; i < InputVector.GetLargeVectorDimension(); ++i)
	    {
	      Vector.Re(i) = (rand() - 32767) * 0.5;
	      Vector.Im(i) = (rand() - 32767) * 0.5;
	    }
	  Vector /= Vector.Norm();	  
	  Vector.AddLinearCombination((InputVector* Vector) / (- InputVector.Norm()), InputVector);
	  Vector.WriteVector(Manager.GetString("output"));
	}
      else
	{
	  RealVector InputVector;
	  if (InputVector.ReadVector(Manager.GetString("random-orthogonal")) == false)
	    {
	      cout << "can't open " << Manager.GetString("random-orthogonal") << endl;
	      return -1;
	    }
	  RealVector Vector(InputVector.GetLargeVectorDimension());
	  for (long i = 0l; i < InputVector.GetLargeVectorDimension(); ++i)
	    {
	      Vector[i] = (rand() - 32767) * 0.5;
	    }
	  Vector.AddLinearCombination((InputVector * Vector) / (-Vector.Norm() * InputVector.Norm()), InputVector);
	  Vector /= Vector.Norm();	  
	  Vector.WriteVector(Manager.GetString("output"));
	}
      cout << "Generated random vector orthogonal to " << Manager.GetString("random-orthogonal") << endl;
      exit(0);
    }

  if (Manager.GetString("description-file") == 0)
    {
      int NbrVectors, ValidVectors=0, VectorDimension=0;
      char** VectorFiles = Manager.GetStrings("states",NbrVectors);
      
      
      Vector **Vectors = new Vector*[NbrVectors], *Result, *Random=NULL;
      
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
	  if (Manager.GetDouble("random-component")!=0.0)
	    {
	      Random = new ComplexVector(VectorDimension);
	      for (int i = 0; i < VectorDimension; ++i)
		{
		  ((ComplexVector*)Random)->Re(i) = (rand() - 32767) * 0.5;
		  ((ComplexVector*)Random)->Im(i) = (rand() - 32767) * 0.5;
		}
	      *((ComplexVector*)(Random)) /= Random->Norm();
	    }
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
	  if (Manager.GetDouble("random-component")!=0.0)
	    {
	      Random = new RealVector(VectorDimension);
	      for (int i = 0; i < VectorDimension; ++i)
		(*(RealVector*)(Random))[i] = (rand() - 32767) * 0.5;
	      *((RealVector*)(Random)) /= Random->Norm();
	    }
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
      if (RandomComponent!=0.0)
	Result->AddLinearCombination(RandomComponent,*Random);
      
      char *OutputName = Manager.GetString("output");
      
      if (Manager.GetBoolean("complex"))
	{
	  if (!Manager.GetBoolean("no-normalize"))
	    *((ComplexVector*)(Result)) /= Result->Norm();
	  ((ComplexVector*)Result)->WriteVector(OutputName);
	}
      else
	{
	  if (!Manager.GetBoolean("no-normalize"))
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
	
      if ((Manager.GetBoolean("complex") == false) && (Description.GetNbrColumns() == 3))
	{
	  RealVector Source;
	  if (Source.ReadVector(Description(0, 0)) == false)
	    {
	      cout << "can't open file " << Description(0, 0) << endl;
	      return -1;	      
	    }
	  double* CoefficientsRe = Description.GetAsDoubleArray(1);
	  double* CoefficientsIm = Description.GetAsDoubleArray(2);
	  Complex* Coefficients = new Complex[Description.GetNbrLines()];
	  if (CoefficientsIm!=NULL)
	    {
	      for (int i = 0; i < Description.GetNbrLines(); ++i)
		Coefficients[i]=Complex(CoefficientsRe[i],CoefficientsIm[i]);
	    }
	  else
	    {
	      for (int i = 0; i < Description.GetNbrLines(); ++i)
		Coefficients[i]=Complex(CoefficientsRe[i]);
	    }
	    
	  ComplexVector Result(Source.GetVectorDimension());
	  Result.AddLinearCombination(Coefficients[0], Source);
// 	  cout << "d*Vector1 = "<<Result<<endl;
 	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      RealVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}
	      Result.AddLinearCombination(Coefficients[i], TmpVector);
// 	      cout << "+d*Vector"<<i<<" = "<<Result<<endl;
	    }
	  if (!Manager.GetBoolean("no-normalize"))
	    Result /= Result.Norm();
	  Result.WriteVector(Manager.GetString("output"));
	  delete [] Coefficients;
	  return 0;
	}
	
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector Result;
	  if (Result.ReadVector(Description(0, 0)) == false)
	    {
	      cout << "can't open file " << Description(0, 0) << endl;
	      return -1;	      
	    }
	  double* CoefficientsRe = Description.GetAsDoubleArray(1);
	  double* CoefficientsIm = Description.GetAsDoubleArray(2);
	  Complex* Coefficients = new Complex[Description.GetNbrLines()];
	  if (CoefficientsIm!=NULL)
	    {
	      for (int i = 0; i < Description.GetNbrLines(); ++i)
		Coefficients[i]=Complex(CoefficientsRe[i],CoefficientsIm[i]);
	    }
	  else
	    {
	      for (int i = 0; i < Description.GetNbrLines(); ++i)
		Coefficients[i]=Complex(CoefficientsRe[i]);
	    }
	  Result *= Coefficients[0];
// 	  cout << "d*Vector1 = "<<Result<<endl;
 	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      ComplexVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}
	      Result.AddLinearCombination(Coefficients[i], TmpVector);
// 	      cout << "+d*Vector"<<i<<" = "<<Result<<endl;
	    }
	  if (!Manager.GetBoolean("no-normalize"))
	    Result /= Result.Norm();
	  Result.WriteVector(Manager.GetString("output"));
	  delete [] Coefficients;
	  return 0;
	}
      if (Manager.GetBoolean("rational"))
	{
#ifdef __GMP__
	  if (Manager.GetBoolean("use-gmp") == false)
#else
	    if (Manager.GetBoolean("use-longlong") == false)	
#endif
	      {
		RationalVector Result;
		if (Result.ReadVector(Description(0, 0)) == false)
		  {
		    cout << "can't open file " << Description(0, 0) << endl;
		    return -1;	      
		  }
		Rational* Coefficients = Description.GetAsRationalArray(1);
		Result *= Coefficients[0];
		for (int i = 1; i < Description.GetNbrLines(); ++i)
		  {
		    RationalVector TmpVector;
		    if (TmpVector.ReadVector(Description(0, i)) == false)
		      {
			cout << "can't open file " << Description(0, i) << endl;
			return -1;	      	    
		      }	      
		    Result.AddLinearCombination(Coefficients[i], TmpVector);
		  }
		Result.WriteVector(Manager.GetString("output"));
		delete [] Coefficients;
	      }
	    else
	      {
		LongRationalVector Result;
		if (Result.ReadVector(Description(0, 0)) == false)
		  {
		    cout << "can't open file " << Description(0, 0) << endl;
		    return -1;	      
		  }
		LongRational* Coefficients = Description.GetAsLongRationalArray(1);
		Result *= Coefficients[0];
		for (int i = 1; i < Description.GetNbrLines(); ++i)
		  {
		    LongRationalVector TmpVector;
		    if (TmpVector.ReadVector(Description(0, i)) == false)
		      {
			cout << "can't open file " << Description(0, i) << endl;
			return -1;	      	    
		      }	      
		    Result.AddLinearCombination(Coefficients[i], TmpVector);
		  }
		Result.WriteVector(Manager.GetString("output"));
		delete [] Coefficients;
	      }
	  return 0;
	}
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
      if (!Manager.GetBoolean("no-normalize"))
	Result /= Result.Norm();
      Result.WriteVector(Manager.GetString("output"));
      delete [] Coefficients;
    }
}
  
