#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

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

  // some running options and help
  OptionManager Manager ("GenericOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "names of the vector files obtained using exact diagonalization");
  (*SystemGroup) += new SingleStringOption  ('\n', "state-list", "provide the names of the vector using a single column formatted text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "states-matrix", "use a series of vectors stored as c");
  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "Assume vectors consist of complex numbers");
  (*SystemGroup) += new BooleanOption  ('s', "scalar-product", "Get the scalar product, not the overlap");
  (*SystemGroup) += new BooleanOption  ('p', "polar", "print complex numbers in polar coordinates");
  (*SystemGroup) += new BooleanOption  ('\n', "conjugate", "Conjugate the second (complex) number");
  (*SystemGroup) += new BooleanOption  ('\n', "discard-sign", "compute sum_i |v1_i * v2_i| instead of sum_i v1_i * v2_i");
  (*SystemGroup) += new BooleanOption  ('x', "no-cross", "calculate only overlap of 1st vector with all others");
  (*SystemGroup) += new SingleDoubleOption  ('t', "threshold", "apply threshold and display only those overlaps exceeding this value",0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "sum", "sum all computed overlaps");
  (*SystemGroup) += new BooleanOption  ('\n', "no-square", "calculate only the scalar products");
  (*SystemGroup) += new BooleanOption  ('n', "normalize", "normalize vectors before calculating any overlaps");
  (*SystemGroup) += new BooleanOption  ('d', "dimension", "show vector dimension");
  (*SystemGroup) += new BooleanOption  ('\n', "quiet", "discard any output except the overlaps");
  (*SystemGroup) += new SingleStringOption ('\n', "matrix", "export the overlap as a binary matrix");
  
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
  bool Scalar = Manager.GetBoolean("scalar-product") | Manager.GetBoolean("no-square");

  double Threshold = Manager.GetDouble("threshold");
  
  int NbrVectors;
  char** VectorFiles = 0;

  if (Manager.GetString("state-list") == 0)
    {
      VectorFiles = Manager.GetStrings("states", NbrVectors);
    }
  else
    {
      MultiColumnASCIIFile Description;
      if (Description.Parse(Manager.GetString("state-list")) == false)
	{
	  Description.DumpErrors(cout);
	  return -1;
	}      
      if (Description.GetNbrColumns() < 1)
	{
	  cout << "wrong number of columns in " << Manager.GetString("state-list") << endl;
	  return -1;
	}
      if (Manager.GetStrings("states", NbrVectors) != 0)
	{
	  int TmpNbrVectors1 = 0;
	  char** TmpVectorFiles1 =  Manager.GetStrings("states", TmpNbrVectors1);
	  char** TmpVectorFiles2 = Description.GetAsStringArray(0);
	  int TmpNbrVectors2 = Description.GetNbrLines();
	  NbrVectors = TmpNbrVectors1 + TmpNbrVectors2;
	  VectorFiles = new char* [NbrVectors];
	  for (int i = 0; i < TmpNbrVectors1; ++i)
	    VectorFiles[i] = TmpVectorFiles1[i];
	  for (int i = 0; i < TmpNbrVectors2; ++i)
	    VectorFiles[i + TmpNbrVectors1] = TmpVectorFiles2[i];
	  delete[] TmpVectorFiles2;
	}
      else
	{
	  VectorFiles = Description.GetAsStringArray(0);
	  NbrVectors = Description.GetNbrLines();
	}
    }

  if ((NbrVectors < 2) && (Manager.GetString("states-matrix")) == 0)
    {
      cout << "At least two vector files are required!"<<endl;
      exit(1);
    }

  if (QuietFlag == false)
    for (int i = 0; i < NbrVectors; ++i)    
      cout << "File "<<i<<"  "<<VectorFiles[i]<<endl;


  int MaxVectors= (Manager.GetBoolean("no-cross") ? 1:NbrVectors);

  bool HaveComplex = Manager.GetBoolean("complex");


  if (Manager.GetString("states-matrix") == 0)
    {
      if (HaveComplex)
	{      
	  Complex sp = 0.0;
	  Complex TotalOverlap = 0.0;
	  ComplexVector State1, State2;
	  ComplexMatrix OverlapMatrix;
	  if (MaxVectors == 1) 
	    OverlapMatrix = ComplexMatrix(MaxVectors, NbrVectors - 1, true);
	  else
	    OverlapMatrix = ComplexMatrix(MaxVectors, NbrVectors, true);
	  for (int i=0; i < MaxVectors; ++i)
	    {
	      if (State1.ReadVector (VectorFiles[i]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[i] << endl;
		  return -1;      
		}
	      if ((i==0)&&(Manager.GetBoolean("dimension")))
		{
		  cout << "Vector dimension = "<<State1.GetVectorDimension() <<endl;
		}
	      if (Manager.GetBoolean("normalize"))
		State1 /= State1.Norm();
	      for (int j= i + 1; j < NbrVectors; ++j)
		{	      
		  if (State2.ReadVector (VectorFiles[j]) == false)
		    {
		      cout << "can't open vector file " << VectorFiles[j] << endl;
		      return -1;      
		    }
		  if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		    {
		      cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		      return -2;
		    }
		  if (Manager.GetBoolean("normalize"))
		    State2 /= State2.Norm();
		  sp = 0.0;
		  if (Manager.GetBoolean("discard-sign"))
		    for (int k = 0; k < State1.GetVectorDimension(); ++k)
		      sp += Norm(State1[k] * State2[k]);
		  else
		    if (Manager.GetBoolean("conjugate"))
		      sp = EuclidianScalarProduct(State1, State2);
		    else
		      sp = State1 * State2;
		  if (Scalar == false)
		    {
		      TotalOverlap += SqrNorm(sp);
		    }
		  else
		    {
		      TotalOverlap += sp;
		    }
		  if (MaxVectors == 1)
		    {
		      OverlapMatrix.SetMatrixElement(0, j - 1, sp);
		    }
		  else
		    {
		      OverlapMatrix.SetMatrixElement(j, i, sp);
		      OverlapMatrix.SetMatrixElement(i, j, sp);
		      OverlapMatrix.SetMatrixElement(i, i, 1.0);
		    }
		  if (Scalar == false)
		    {
		      if (QuietFlag == false)
			{
			  if (SqrNorm(sp)>=Threshold)
			    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
			}
		      else
			cout << SqrNorm(sp) << endl;
		    }
		  else
		    {
		      if (QuietFlag == false)
			{
			  if (Norm(sp)>=Threshold)
			    {
			      if (Manager.GetBoolean("polar"))
				cout << "<"<<i<<"|"<<j<<"> = " << Norm(sp)<<"*Exp("<< Arg(sp)/M_PI<< " I Pi)" << endl;
			      else
				cout << "<"<<i<<"|"<<j<<"> = " << sp << endl;
			    }
			}
		      else
			cout << sp.Re << " " << sp.Im << endl;
		    }
		}
	    }
	  if (Manager.GetBoolean("sum") == true)
	    {
	      if (Scalar == false)
		{
		  if (QuietFlag == false)
		    cout << "Total overlap = " << TotalOverlap.Re << endl;
		  else
		    cout << TotalOverlap.Re << endl;
		}
	      else
		{
		  if (QuietFlag == false)
		    cout << "Total overlap = " << TotalOverlap << endl;
		  else
		    cout << TotalOverlap << endl;
		}
	    }
	  if (Manager.GetString("matrix") != 0)
	    {
	      OverlapMatrix.WriteMatrix(Manager.GetString("matrix"));
	    }
	}
      else // real vectors
	{
	  double TotalOverlap = 0.0;
	  double sp = 0.0;
	  RealVector State1, State2;
	  RealMatrix OverlapMatrix;
	  if (MaxVectors == 1) 
	    OverlapMatrix = RealMatrix(MaxVectors, NbrVectors - 1, true);
	  else
	    OverlapMatrix = RealMatrix(MaxVectors, NbrVectors, true);
	  for (int i = 0; i < MaxVectors; ++i)
	    {
	      if (State1.ReadVector (VectorFiles[i]) == false)
		{
		  cout << "can't open vector file " << VectorFiles[i] << endl;
		  return -1;      
		}
	      if ((i==0)&&(Manager.GetBoolean("dimension")))
		{
		  cout << "Vector dimension = "<<State1.GetVectorDimension() <<endl;
		}
	      if (Manager.GetBoolean("normalize"))
		State1 /= State1.Norm();
	      for (int j = i + 1; j < NbrVectors; ++j)
		{	      
		  if (State2.ReadVector (VectorFiles[j]) == false)
		    {
		      cout << "can't open vector file " << VectorFiles[j] << endl;
		      return -1;      
		    }
		  if (State1.GetVectorDimension() != State2.GetVectorDimension() )
		    {
		      cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
		      return -2;
		    }
		  if (Manager.GetBoolean("normalize"))
		    State2/=State2.Norm();
		  sp = 0.0;
		  if (Manager.GetBoolean("discard-sign"))
		    {
		      for (int k = 0; k < State1.GetVectorDimension(); ++k)
			sp += fabs(State1[k] * State2[k]);
		    }
		  else
		    {
		      sp = State1 * State2 ;
		    }
		  if (MaxVectors == 1)
		    {
		      OverlapMatrix.SetMatrixElement(0, j - 1, sp);
		    }
		  else
		    {
		      OverlapMatrix.SetMatrixElement(j, i, sp);
		      OverlapMatrix.SetMatrixElement(i, j, sp);
		      OverlapMatrix.SetMatrixElement(i, i, 1.0);
		    }
		  if (Scalar == false)
		    {
		      TotalOverlap += SqrNorm(sp);
		    }
		  else
		    {
		      TotalOverlap += sp;
		    }
		  if (Scalar == false)
		    {
		      if (QuietFlag == false)
			{
			  if (SqrNorm(sp) >= Threshold)
			    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
			}
		      else
			{
			  cout << SqrNorm(sp) << endl;
			}
		    }
		  else
		    {
		      if (QuietFlag == false)
			{
			  if (fabs(sp)>=Threshold)
			    cout << "<"<<i<<"|"<<j<<"> = " << sp << endl;
			}
		      else
			{
			  cout << sp << endl;
			}
		    }
		  
		}
	    }
	  if (Manager.GetBoolean("sum") == true)
	    {
	      if (QuietFlag == false)
		cout << "Total overlap = " << TotalOverlap << endl;
	      else
		cout << TotalOverlap << endl;
	    }
	  if (Manager.GetString("matrix") != 0)
	    {
	      OverlapMatrix.WriteMatrix(Manager.GetString("matrix"));
	    }
	}
    }
  else
    {
      if (HaveComplex)
	{
	  ComplexMatrix VectorsAsAMatrix;
	  if (VectorsAsAMatrix.ReadMatrix(Manager.GetString("states-matrix")) == false)
	    {
	      cout << "error while reading matrix " << Manager.GetString("states-matrix") << endl;
	      return -1;
	    }
	  ComplexVector* StateVectors = new ComplexVector [NbrVectors + VectorsAsAMatrix.GetNbrColumn()];
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      if (StateVectors[i].ReadVector (VectorFiles[i]) == false)
		{
		  cout << "error while reading vector " << Manager.GetString(VectorFiles[i]) << endl;
		  return -1;		  
		}
	      if (StateVectors[i].GetVectorDimension() != VectorsAsAMatrix.GetNbrRow())
		{
		  cout << "dimension mismacth bewteen vector " << Manager.GetString(VectorFiles[i]) << " (" << StateVectors[i].GetVectorDimension() << ") and "
		       << Manager.GetString("states-matrix") << " (" << VectorsAsAMatrix.GetNbrRow() << ")" << endl;		  
		}
	    }
	  for (int i = 0; i < VectorsAsAMatrix.GetNbrColumn(); ++i)
	    {
	      StateVectors[i + NbrVectors] = VectorsAsAMatrix[i];
	    }
	  NbrVectors += VectorsAsAMatrix.GetNbrColumn();
	  Complex sp = 0.0;
	  Complex TotalOverlap = 0.0;
	  ComplexMatrix OverlapMatrix;
	  if (MaxVectors == 1) 
	    OverlapMatrix = ComplexMatrix(MaxVectors, NbrVectors - 1, true);
	  else
	    OverlapMatrix = ComplexMatrix(MaxVectors, NbrVectors, true);
	  for (int i=0; i < MaxVectors; ++i)
	    {
	      if (Manager.GetBoolean("normalize"))
		StateVectors[i] /= StateVectors[i].Norm();
	      for (int j= i + 1; j < NbrVectors; ++j)
		{	      
		  if (Manager.GetBoolean("normalize"))
		    StateVectors[j] /= StateVectors[j].Norm();
		  sp = 0.0;
		  if (Manager.GetBoolean("discard-sign"))
		    for (int k = 0; k < StateVectors[i].GetVectorDimension(); ++k)
		      sp += Norm(StateVectors[i][k] * StateVectors[j][k]);
		  else
		    if (Manager.GetBoolean("conjugate"))
		      sp = EuclidianScalarProduct(StateVectors[i], StateVectors[j]);
		    else
		      sp = StateVectors[i] * StateVectors[j];
		  if (Scalar == false)
		    {
		      TotalOverlap += SqrNorm(sp);
		    }
		  else
		    {
		      TotalOverlap += sp;
		    }
		  if (MaxVectors == 1)
		    {
		      OverlapMatrix.SetMatrixElement(0, j - 1, sp);
		    }
		  else
		    {
		      OverlapMatrix.SetMatrixElement(j, i, sp);
		      OverlapMatrix.SetMatrixElement(i, j, sp);
		      OverlapMatrix.SetMatrixElement(i, i, 1.0);
		    }
		  if (Scalar == false)
		    {
		      if (QuietFlag == false)
			{
			  if (SqrNorm(sp)>=Threshold)
			    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
			}
		      else
			cout << SqrNorm(sp) << endl;
		    }
		  else
		    {
		      if (QuietFlag == false)
			{
			  if (Norm(sp)>=Threshold)
			    {
			      if (Manager.GetBoolean("polar"))
				cout << "<"<<i<<"|"<<j<<"> = " << Norm(sp)<<"*Exp("<< Arg(sp)/M_PI<< " I Pi)" << endl;
			      else
				cout << "<"<<i<<"|"<<j<<"> = " << sp << endl;
			    }
			}
		      else
			cout << sp.Re << " " << sp.Im << endl;
		    }
		}
	    }
	  if (Manager.GetBoolean("sum") == true)
	    {
	      if (Scalar == false)
		{
		  if (QuietFlag == false)
		    cout << "Total overlap = " << TotalOverlap.Re << endl;
		  else
		    cout << TotalOverlap.Re << endl;
		}
	      else
		{
		  if (QuietFlag == false)
		    cout << "Total overlap = " << TotalOverlap << endl;
		  else
		    cout << TotalOverlap << endl;
		}
	    }
	  if (Manager.GetString("matrix") != 0)
	    {
	      OverlapMatrix.WriteMatrix(Manager.GetString("matrix"));
	    }
	  delete[] StateVectors;
	}
      else
	{
	  RealMatrix VectorsAsAMatrix;
	  if (VectorsAsAMatrix.ReadMatrix(Manager.GetString("states-matrix")) == false)
	    {
	      cout << "error while reading matrix " << Manager.GetString("states-matrix") << endl;
	      return -1;
	    }
	  RealVector* StateVectors = new RealVector [NbrVectors + VectorsAsAMatrix.GetNbrColumn()];
	  for (int i = 0; i < NbrVectors; ++i)
	    {
	      if (StateVectors[i].ReadVector (VectorFiles[i]) == false)
		{
		  cout << "error while reading vector " << Manager.GetString(VectorFiles[i]) << endl;
		  return -1;		  
		}
	      if (StateVectors[i].GetVectorDimension() != VectorsAsAMatrix.GetNbrRow())
		{
		  cout << "dimension mismacth bewteen vector " << Manager.GetString(VectorFiles[i]) << " (" << StateVectors[i].GetVectorDimension() << ") and "
		       << Manager.GetString("states-matrix") << " (" << VectorsAsAMatrix.GetNbrRow() << ")" << endl;		  
		}
	    }
	  for (int i = 0; i < VectorsAsAMatrix.GetNbrColumn(); ++i)
	    {
	      StateVectors[i + NbrVectors] = VectorsAsAMatrix[i];
	    }
	  NbrVectors += VectorsAsAMatrix.GetNbrColumn();
	  double TotalOverlap = 0.0;
	  double sp = 0.0;
	  RealMatrix OverlapMatrix;
	  if (MaxVectors == 1) 
	    OverlapMatrix = RealMatrix(MaxVectors, NbrVectors - 1, true);
	  else
	    OverlapMatrix = RealMatrix(MaxVectors, NbrVectors, true);
	  for (int i = 0; i < MaxVectors; ++i)
	    {
	      if (Manager.GetBoolean("normalize"))
		StateVectors[i] /= StateVectors[i].Norm();
	      for (int j = i + 1; j < NbrVectors; ++j)
		{	      
		  if (Manager.GetBoolean("normalize"))
		    StateVectors[j] /= StateVectors[j].Norm();
		  sp = 0.0;
		  if (Manager.GetBoolean("discard-sign"))
		    {
		      for (int k = 0; k < StateVectors[i].GetVectorDimension(); ++k)
			sp += fabs(StateVectors[i][k] * StateVectors[j][k]);
		    }
		  else
		    {
		      sp = StateVectors[i] * StateVectors[j] ;
		    }
		  if (MaxVectors == 1)
		    {
		      OverlapMatrix.SetMatrixElement(0, j - 1, sp);
		    }
		  else
		    {
		      OverlapMatrix.SetMatrixElement(j, i, sp);
		      OverlapMatrix.SetMatrixElement(i, j, sp);
		      OverlapMatrix.SetMatrixElement(i, i, 1.0);
		    }
		  if (Scalar == false)
		    {
		      TotalOverlap += SqrNorm(sp);
		    }
		  else
		    {
		      TotalOverlap += sp;
		    }
		  if (Scalar == false)
		    {
		      if (QuietFlag == false)
			{
			  if (SqrNorm(sp) >= Threshold)
			    cout << "Overlap |<"<<i<<"|"<<j<<">|^2 = " << SqrNorm(sp) << endl;
			}
		      else
			{
			  cout << SqrNorm(sp) << endl;
			}
		    }
		  else
		    {
		      if (QuietFlag == false)
			{
			  if (fabs(sp)>=Threshold)
			    cout << "<"<<i<<"|"<<j<<"> = " << sp << endl;
			}
		      else
			{
			  cout << sp << endl;
			}
		    }
		  
		}
	    }
	  if (Manager.GetBoolean("sum") == true)
	    {
	      if (QuietFlag == false)
		cout << "Total overlap = " << TotalOverlap << endl;
	      else
		cout << TotalOverlap << endl;
	    }
	  if (Manager.GetString("matrix") != 0)
	    {
	      OverlapMatrix.WriteMatrix(Manager.GetString("matrix"));
	    }
	  delete[] StateVectors;
	}
    }
  return 0;
}
