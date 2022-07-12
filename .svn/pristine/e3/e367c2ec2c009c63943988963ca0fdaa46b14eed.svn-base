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
  OptionManager Manager ("MergeVectors" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  
  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of all vector files that should be merged");
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) +=  new BooleanOption  ('\n', "rational", "input vectors are rational vectors");
#ifdef __GMP__
  (*SystemGroup) += new BooleanOption  ('\n', "use-gmp", "use arbitrary precision integers instead of fixed precision integers in rational mode");
#else
  (*SystemGroup) += new BooleanOption  ('\n', "use-longlong", "use 128bit(64bits) integers instead of 64bits(32bits) integers in rational mode");
#endif
  (*SystemGroup) += new SingleStringOption  ('o', "output", "names of output filename","merge.vec");
  (*SystemGroup) += new SingleStringOption  ('f', "description-file", "build the merge from a vector set described in a text file");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  
  if (Manager.GetString("description-file") == 0)
    {
      int NbrVectors, VectorDimension = 0;
      char** VectorFiles = Manager.GetStrings("states",NbrVectors);
      
      if (Manager.GetBoolean("complex"))
	{
	  ComplexVector Result;
	  ComplexVector TmpVector;
	  if (Result.ReadVector(VectorFiles[0]) == false)
	    {
	      cout << "can't open file " << VectorFiles[0] << endl;
	      return -1;	      	    
	    }	
	  VectorDimension = Result.GetVectorDimension();
	  for (int i = 1; i < NbrVectors; ++i)
	    {
	      if (TmpVector.ReadVector(VectorFiles[i]) == false)
		{
		  cout << "can't open file " <<VectorFiles[i] << endl;
		  return -1;	      	    
		}      
	      if (TmpVector.GetVectorDimension() != VectorDimension)
		{
		  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
		  exit(1);
		}
	      
	      for (int k = 0 ; k < VectorDimension; k++)
		{
		  if (TmpVector[k] != 0)
		    {
		      if(Result[k] == 0)
			Result[k] = TmpVector[k];
		      else
			{
			  if(Result[k] != TmpVector[k])
			    {
			      cout <<"Two non zero components are different"<<endl;
			      exit(1);
			    }
			}
		    }
		}
	    }
	  Result.WriteVector(Manager.GetString("output"));
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
		RationalVector TmpVector;
		if (Result.ReadVector(VectorFiles[0]) == false)
		  {
		    cout << "can't open file " << VectorFiles[0] << endl;
		    return -1;	      	    
		  }	
		VectorDimension = Result.GetVectorDimension();
		for (int i = 1; i < NbrVectors; ++i)
		  {
		    if (TmpVector.ReadVector(VectorFiles[i]) == false)
		      {
			cout << "can't open file " <<VectorFiles[i] << endl;
			return -1;	      	    
		      }	      
		    if (TmpVector.GetVectorDimension() != VectorDimension)
		      {
			cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
			exit(1);
		      }
		    
		    for (int k = 0 ; k < VectorDimension; k++)
		      {
			if (TmpVector[k] != 0l)
			  {
			    if(Result[k] == 0l)
			      Result[k] = TmpVector[k];
			    else
			      {
				if(Result[k] != TmpVector[k])
				  {
				    cout <<"Two non zero components are different"<<endl;
				    exit(1);
				  }
			      }
			  }
		      }
		  }
		Result.WriteVector(Manager.GetString("output"));
		return 0;
	      }
	  
	  LongRationalVector Result;
	  LongRationalVector TmpVector;
	  if (Result.ReadVector(VectorFiles[0]) == false)
	    {
	      cout << "can't open file " << VectorFiles[0] << endl;
	      return -1;	      	    
	    }	
	  cout <<"Vecteur : "<<VectorFiles[0]<<endl;
	  VectorDimension = Result.GetVectorDimension();
	  for (int i = 1; i<NbrVectors; ++i)
	    {
	      cout <<"Vecteur : "<<VectorFiles[i]<<endl;
	      if (TmpVector.ReadVector(VectorFiles[i]) == false)
		{
		  cout << "can't open file " <<VectorFiles[i] << endl;
		  return -1;	      	    
		}	      
	      if (TmpVector.GetVectorDimension() != VectorDimension)
		{
		  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
		  exit(1);
		}
	      for (int k = 0 ; k < TmpVector.GetVectorDimension(); k++)
		{
		  if (TmpVector[k].IsZero() == false)
		    {
		      if(Result[k].IsZero() == true)
			Result[k] = TmpVector[k];
		      else
			{
			  if(Result[k] != TmpVector[k])
			    {
			      cout <<"Two non zero components are different"<<endl;
			      cout << k <<" " << Result[k] << " "<<TmpVector[k]<<endl;
			      exit(1);
			    }
			}
		    }
		}
	    }
	  Result.WriteVector(Manager.GetString("output"));
	  return 0;
	}
      RealVector Result;
      RealVector TmpVector;
      if (Result.ReadVector(VectorFiles[0]) == false)
	{
	  cout << "can't open file " << VectorFiles[0] << endl;
	  return -1;	      	    
	}	
      VectorDimension = Result.GetVectorDimension();
      for (int i = 1; i<NbrVectors; ++i)
	{
	  if (TmpVector.ReadVector(VectorFiles[i]) == false)
	    {
	      cout << "can't open file " <<VectorFiles[i] << endl;
	      return -1;	      	    
	    }	      
	  if (TmpVector.GetVectorDimension() != VectorDimension)
	    {
	      cout<<"Dimension of vector "<< VectorFiles[i] << " does not match size of previous vectors!" << endl;
	      exit(1);
	    }
	  
	  for (int k = 0 ; k < VectorDimension; k++)
	    {
	      if (TmpVector[k] != 0)
		{
		  if(Result[k] == 0)
		    Result[k] = TmpVector[k];
		  else
		    {
		      if(Result[k] != TmpVector[k])
			{
			  cout <<"Two non zero components are different"<<endl;
			  exit(1);
			}
		    }
		}
	    }
	}
      Result.WriteVector(Manager.GetString("output"));
      return 0;	
    }	
  else
    {
      MultiColumnASCIIFile Description;
      if (Description.Parse(Manager.GetString("description-file")) == false)
	{
	  Description.DumpErrors(cout);
	  return -1;
	}
      if (Description.GetNbrColumns() < 1)
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
	  
 	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      ComplexVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}
	      
	      for (int k = 0 ; k < TmpVector.GetVectorDimension(); k++)
		{
		  if (TmpVector[k] != 0)
		    {
		      if(Result[k] == 0)
			Result[k] = TmpVector[k];
		      else
			{
			  if(Result[k] != TmpVector[k])
			    {
			      cout <<"Two non zero components are different"<<endl;
			      exit(1);
			    }
			}
		    }
		}
	    }
	  Result.WriteVector(Manager.GetString("output"));
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
		
		for (int i = 1; i < Description.GetNbrLines(); ++i)
		  {
		    RationalVector TmpVector;
		    if (TmpVector.ReadVector(Description(0, i)) == false)
		      {
			cout << "can't open file " << Description(0, i) << endl;
			return -1;	      	    
		      }
		    
		    for (int k = 0 ; k < TmpVector.GetVectorDimension(); k++)
		      {
			if (TmpVector[k] != 0)
			  {
			    if(Result[k] == 0)
			      Result[k] = TmpVector[k];
			    else
			      {
				if(Result[k] != TmpVector[k])
				  {
				    cout <<"Two non zero components are different"<<endl;
				    exit(1);
				  }
			      }
			  }
		      }
		  }
		Result.WriteVector(Manager.GetString("output"));
	      }
	    else
	      {
		LongRationalVector Result;
		if (Result.ReadVector(Description(0, 0)) == false)
		  {
		    cout << "can't open file " << Description(0, 0) << endl;
		    return -1;	      
		  }
		for (int i = 1; i < Description.GetNbrLines(); ++i)
		  {
		    LongRationalVector TmpVector;
		    if (TmpVector.ReadVector(Description(0, i)) == false)
		      {
			cout << "can't open file " << Description(0, i) << endl;
			return -1;	      	    
		      }	      
		    for (int k = 0 ; k < TmpVector.GetVectorDimension(); k++)
		      {
			if (TmpVector[k].IsZero() == false)
			  {
			    if(Result[k].IsZero() == true)
			      Result[k] = TmpVector[k];
			    else
			      {
				if(Result[k] != TmpVector[k])
				  {
				    cout <<"Two non zero components are different"<<endl;
				    exit(1);
				  }
			      }
			  }
			
		      }
		  }
		Result.WriteVector(Manager.GetString("output"));
		return 0;
	      }
	}
      else
	{
	  RealVector Result;
	  if (Result.ReadVector(Description(0, 0)) == false)
		{
		  cout << "can't open file " << Description(0, 0) << endl;
		  return -1;	      
		}
	  for (int i = 1; i < Description.GetNbrLines(); ++i)
	    {
	      RealVector TmpVector;
	      if (TmpVector.ReadVector(Description(0, i)) == false)
		{
		  cout << "can't open file " << Description(0, i) << endl;
		  return -1;	      	    
		}	      
	      for (int k = 0 ; k < TmpVector.GetVectorDimension(); k++)
		{
		  if (TmpVector[k] != 0)
		    {
		      if(TmpVector[k] == 0)
			TmpVector[k] = TmpVector[k];
		      else
			{
			  if(fabs(Result[k]-TmpVector[k]) > 1e-10)
			    {
			      cout <<"Two non zero components are different"<<endl;
			      exit(1);
			    }
			}
		    }
		}
	    }
	  Result.WriteVector(Manager.GetString("output"));
	}
    }
}
