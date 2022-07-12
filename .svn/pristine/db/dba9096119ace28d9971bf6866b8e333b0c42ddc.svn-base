#include "MaximallyCondensedStateOnLattice.h"

#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"

#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include "Tools/NewUnconstrainedOptimizsation.h"

#include <cstdlib>
#include <sys/time.h>
#include <iostream>
using std::cout;
using std::endl;


// flag for testing mode
//#define TESTING

// constructor for contact interactions on a square lattice
//
// nbrStates = number of quantum states
// states = state vectors to superpose
// space = Hilbert-space of states
// lx = Lx dimension
// ly = Ly dimension
// sublattices = number of sublattices
// randomGenerator = external random number generator
MaximallyCondensedStateOnLattice::MaximallyCondensedStateOnLattice(AbstractArchitecture *architecture, int nbrStates, ComplexVector *states, ParticleOnLattice* space, int lx, int ly, int sublattices, AbstractRandomNumberGenerator *randomGenerator)
{
  this->Architecture=architecture;
  this->NbrVectors = nbrStates;
  this->Vectors = states;
  this->Space = space;
  // lattice dimensions
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSubLattices = sublattices;
  this->LastMaximumEV = 0.0;
  this->SphereParametrization= NSphereParameters(nbrStates,true);
  this->VariationalParameters.Resize(this->SphereParametrization.GetNbrParameters());
  this->VariationalParameters.ClearVector();

  this->DensityMatrixDimension = Lx * Ly * NbrSubLattices;

  // storage for density matrices
  this->DiagonalDensityMatrices = new HermitianMatrix[NbrVectors];
  this->NbrOffDiagonal = (this->NbrVectors * (this->NbrVectors - 1));
  this->OffDiagonalDensityMatrices = new ComplexMatrix[NbrOffDiagonal];
  for (int i=0; i<NbrVectors; ++i)
    this->DiagonalDensityMatrices[i].Resize(DensityMatrixDimension,DensityMatrixDimension);
  for (int i=0; i<NbrOffDiagonal; ++i)
    this->OffDiagonalDensityMatrices[i].Resize(DensityMatrixDimension,DensityMatrixDimension);
  this->CurrentDensityMatrix.Resize(DensityMatrixDimension, DensityMatrixDimension);
  this->CurrentHermitianMatrix.Resize(DensityMatrixDimension, DensityMatrixDimension);
  
  // matrices as temporary space for calculations
  this->M.Resize(DensityMatrixDimension,DensityMatrixDimension);
  this->Q.Resize(DensityMatrixDimension,DensityMatrixDimension); 

  this->EvaluateDensityMatrices();
  
  if (randomGenerator!=NULL)
    {
      this->RandomNumbers = randomGenerator;
      ExternalGenerator=true;
    }
  else
    {
      timeval RandomTime;
      gettimeofday (&(RandomTime), 0);
      this->RandomNumbers = new NumRecRandomGenerator(RandomTime.tv_sec);
      ExternalGenerator=false;
    }
}

  
 

// destructor
//
MaximallyCondensedStateOnLattice::~MaximallyCondensedStateOnLattice()
{
  if (NbrVectors>0)
    {
      delete [] DiagonalDensityMatrices;
      delete [] OffDiagonalDensityMatrices;
    }
  if (ExternalGenerator==false)
    delete this->RandomNumbers;
}

// get the parameters of the Many-Body state that was last calculated
// return = state
ComplexVector & MaximallyCondensedStateOnLattice::GetVariationalParameters()
{
  return this->SphereParametrization.GetComplexCoordinates();
}

// set trial parameters
void MaximallyCondensedStateOnLattice::SetVariationalParameters(RealVector &variationalParameters)
{
  this->SphereParametrization.SetParameters(&(variationalParameters[0]));
}

// randomize trial parameters
void MaximallyCondensedStateOnLattice::RandomizeVariationalParameters()
{
  int i=0;
  for (; i<VariationalParameters.GetVectorDimension()/2; ++i)
    VariationalParameters[i]=M_PI*RandomNumbers->GetRealRandomNumber();
  for (; i<VariationalParameters.GetVectorDimension(); ++i)
    VariationalParameters[i]=2.0*M_PI*RandomNumbers->GetRealRandomNumber();
  this->SphereParametrization.SetParameters(&(VariationalParameters[0]));
}


// get the wavefunction corresponding to the current parameters
// return = complex vector of local amplitudes and phases
ComplexVector MaximallyCondensedStateOnLattice::GetWaveFunction()
{
  ComplexVector TmpParameters = this->SphereParametrization.GetComplexCoordinates();
  ComplexVector Result;
  Result.Copy(Vectors[0],TmpParameters[0]);
  for (int i=1; i<NbrVectors; ++i)
    Result.AddLinearCombination(TmpParameters[i],Vectors[i]);
  return ComplexVector(Result);
}
  
  
// optimize wavefunction starting from present settings of VariationalParameters
// nbrEigenvals = number of eigenvalues that should be summed up for the optimization
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double MaximallyCondensedStateOnLattice::Optimize(int nbrEigenvals, double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = SphereParametrization.GetNbrParameters();
  this->NbrEvaluations=0;
  this->NbrEigenvalues=nbrEigenvals;
  if (this->NbrEigenvalues>DensityMatrixDimension)
    this->NbrEigenvalues=DensityMatrixDimension-1;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1;
  int rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = &(this->VariationalParameters[0]);
  double (MaximallyCondensedStateOnLattice::*TargetFunction)(int, double*)=&MaximallyCondensedStateOnLattice::EvaluateCondensateFraction;
  MaximallyCondensedStateOnLattice *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << "total: "<<NbrEvaluations<< " evaluations"<<endl;
  delete [] Work;
  return Result;
}

double MaximallyCondensedStateOnLattice::SimplexOptimize(double targetSize, int maxIter, double initialStep)
{
  cout << "Attention: MaximallyCondensedStateOnLattice::SimplexOptimize not implemented"<<endl;
  return 0.0;
}



// target function for optimizer routine:
double MaximallyCondensedStateOnLattice::EvaluateCondensateFraction(int nbrParameters, double *x)
{
  for (int i=0; i<this->SphereParametrization.GetNbrParameters(); ++i)
    if (this->VariationalParameters[i]!=x[i])
      this->VariationalParameters[i]=x[i];
  this->SetVariationalParameters(this->VariationalParameters);
  Complex TmpC;
  Complex TmpC2;
  this->CurrentDensityMatrix.ClearMatrix();
  this->ResultingParameters = SphereParametrization.GetComplexCoordinates();
  for (int n=0; n<NbrVectors; ++n)
    {
      for (int m=0; m<n; ++m)
	CurrentDensityMatrix.AddLinearCombination(this->ResultingParameters[n]*Conj(this->ResultingParameters[m]),
					   OffDiagonalDensityMatrices[(NbrVectors-1)*n+m]);
      CurrentDensityMatrix.AddLinearCombination(this->ResultingParameters[n]*Conj(this->ResultingParameters[n]),DiagonalDensityMatrices[n]);
      for (int m=n+1; m<NbrVectors; ++m)
	CurrentDensityMatrix.AddLinearCombination(this->ResultingParameters[n]*Conj(this->ResultingParameters[m]),
					   OffDiagonalDensityMatrices[(NbrVectors-1)*n+m-1]);
    }
  for (int i=0; i<DensityMatrixDimension; ++i)
    {
      CurrentDensityMatrix.GetMatrixElement(i,i,TmpC);
      if (fabs(TmpC.Im)>1e-13)
	{
	  cout << "Error: Imaginary part on diagonal in element ("<<i<<","<<i<<"): "<< TmpC <<" !"<<endl;
	  exit(1);
	}
      CurrentHermitianMatrix.SetMatrixElement(i,i,TmpC.Re);
      for (int j=i+1; j<DensityMatrixDimension; ++j)
	{
	  CurrentDensityMatrix.GetMatrixElement(i,j,TmpC);
	  CurrentDensityMatrix.GetMatrixElement(j,i,TmpC2);
	  if (Norm(TmpC-Conj(TmpC2))>1e-13)
	    {
	      cout << "Error: Matrix not hermitian in elements ("<<i<<","<<j<<"): "<< TmpC << " vs "<< TmpC2 <<" !"<<endl;
	      exit(1);
	    }
	  CurrentHermitianMatrix.SetMatrixElement(i,j,TmpC);
	}
    }
  ++this->NbrEvaluations;
  CurrentHermitianMatrix.Diagonalize(M, Q, 1e-12, 1000);  
  this->LastMaximumEV = M[DensityMatrixDimension-1];
  for (int i=1; i<this->NbrEigenvalues; ++i)
    this->LastMaximumEV += M[DensityMatrixDimension-1-i];
  return -this->LastMaximumEV;
}

// evaluate all interaction factors
//   
void MaximallyCondensedStateOnLattice::EvaluateDensityMatrices()
{
  ComplexVector TargetVector;
  Complex Tmp;
  Complex *ScalarProducts = new Complex[NbrVectors];
  int CreationIndex, AnnihilationIndex;
  TargetVector.Resize(Vectors[0].GetVectorDimension());
  ParticleOnLatticeOneBodyOperator *DensityOperator= new ParticleOnLatticeOneBodyOperator(this->Space);
#ifdef TESTING
  // inefficient implementation, for testing
  for (int n=0; n<NbrVectors; ++n)
    for (int m=0; m<NbrVectors; ++m)
      {
	for (int CreationX=0; CreationX<this->Lx; ++CreationX)
	  for (int CreationY=0; CreationY<this->Ly; ++CreationY)
	    for (int CreationSub=0; CreationSub<this->NbrSubLattices; ++CreationSub)
	      {
		CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, CreationSub, Tmp);	
		for (int AnnihilationX=0; AnnihilationX<this->Lx; ++AnnihilationX)
		  for (int AnnihilationY=0; AnnihilationY<this->Ly; ++AnnihilationY)
		    for (int AnnihilationSub=0; AnnihilationSub<this->NbrSubLattices; ++AnnihilationSub)
		      {
			AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, AnnihilationSub, Tmp);
			DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
			Tmp=DensityOperator->MatrixElement(Vectors[m], Vectors[n]);
			CurrentDensityMatrix.SetMatrixElement(CreationIndex, AnnihilationIndex, Tmp);
		      }
	      }
	cout << "<"<<m<<"|rho|"<<n<<">="<<endl<<CurrentDensityMatrix;
      }

  
  ComplexVector TmpV1 (Vectors[0].GetVectorDimension(), true);
  ComplexVector TmpV2 (Vectors[0].GetVectorDimension(), true);

  ComplexMatrix TmpMatrix(Vectors[0].GetVectorDimension(),Vectors[0].GetVectorDimension());
  for (int CreationX=0; CreationX<this->Lx; ++CreationX)
    for (int CreationY=0; CreationY<this->Ly; ++CreationY)
      for (int CreationSub=0; CreationSub<this->NbrSubLattices; ++CreationSub)
	{
	  CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, CreationSub, Tmp);	
	  for (int AnnihilationX=0; AnnihilationX<this->Lx; ++AnnihilationX)
	    for (int AnnihilationY=0; AnnihilationY<this->Ly; ++AnnihilationY)
	      for (int AnnihilationSub=0; AnnihilationSub<this->NbrSubLattices; ++AnnihilationSub)
		{
		  AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, AnnihilationSub, Tmp);
		  DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
		  for (int i = 0; i < Vectors[0].GetVectorDimension(); i++)
		    {
		      TmpV1[i] = Complex(1.0, 0.0);
		      DensityOperator->Multiply(TmpV1, TmpV2);
		      for (int j = 0; j < Vectors[0].GetVectorDimension(); j++)
			TmpMatrix.SetMatrixElement(i, j, TmpV2[j]);
		      TmpV1[i] = Complex(0.0, 0.0);
		    }
		  cout << "Density operator (A^+_"<<CreationIndex<<" A_"<<AnnihilationIndex<<")="<<endl<<TmpMatrix<<endl;
		}
	}
  for (int n=0; n<NbrVectors; ++n)
    cout << "Vector["<<n<<"]="<<endl<<Vectors[n]<<endl;
#endif
  
  for (int n=0; n<NbrVectors; ++n)
    {
      for (int CreationX=0; CreationX<this->Lx; ++CreationX)
	for (int CreationY=0; CreationY<this->Ly; ++CreationY)
	  for (int CreationSub=0; CreationSub<this->NbrSubLattices; ++CreationSub)
	    {
	      CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, CreationSub, Tmp);	
	      for (int AnnihilationX=0; AnnihilationX<this->Lx; ++AnnihilationX)
		for (int AnnihilationY=0; AnnihilationY<this->Ly; ++AnnihilationY)
		  for (int AnnihilationSub=0; AnnihilationSub<this->NbrSubLattices; ++AnnihilationSub)
		    {
		      AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, AnnihilationSub, Tmp);
		      DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
		      // calculate possible matrix elements in subspace of vectors
		      VectorOperatorMultiplyOperation Operation(DensityOperator,&(Vectors[n]),&TargetVector);
		      Operation.ApplyOperation(this->Architecture);
		      MultipleComplexScalarProductOperation Operation2(&TargetVector, Vectors, NbrVectors, ScalarProducts);
		      Operation2.ApplyOperation(this->Architecture);
#ifdef TESTING
		      cout << "Checking operator multiplication: "<<endl;
		      cout << TargetVector<<endl;
		      for (int m=0; m<NbrVectors; ++m)
			cout << "A^+_"<<CreationIndex<<" A_"<<AnnihilationIndex<<": Scalar Product ("<<n<<","<<m<<")="
			     <<ScalarProducts[m]<<endl;
#endif
		      for (int m=0; m<n; ++m)
			OffDiagonalDensityMatrices[(NbrVectors-1)*n+m].
			  SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[m]));
		      DiagonalDensityMatrices[n].SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[n]));
		      for (int m=n+1; m<NbrVectors; ++m)
			OffDiagonalDensityMatrices[(NbrVectors-1)*n+m-1].
			  SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[m]));
		    }
	    }
    }
#ifdef TESTING
  cout << "Diagonal matrices: "<<endl;
  for (int n=0; n<NbrVectors; ++n)
    cout << "<"<<n<<"|rho|"<<n<<">="<<endl<<DiagonalDensityMatrices[n];
  cout << "Off-diagonal matrices: "<<endl;
  for (int n=1; n<NbrVectors; ++n)
    for (int m=0; m<n; ++m)
      {
	cout << "<"<<m<<"|rho|"<<n<<">="<<endl<<OffDiagonalDensityMatrices[(NbrVectors-1)*n+m];
	cout << "<"<<n<<"|rho|"<<m<<">="<<endl<<OffDiagonalDensityMatrices[(NbrVectors-1)*m+n-1];
      }
#endif
  delete [] ScalarProducts;
  delete DensityOperator;
}
