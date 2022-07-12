#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

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


void GetTranslationMatrix(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
			  ComplexVector *Vectors, ComplexMatrix &MatrixRepresentation,
			  ComplexVector &TmpState, ArchitectureManager &Architecture)
{
  Complex Tmp;
  for (int i=0; i<NbrVectors; ++i)
    {
      VectorOperatorMultiplyOperation Operation (Operator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());           
      for (int j=0; j<NbrVectors; ++j)
	{
	  Tmp = Vectors[j] * TmpState;
	  MatrixRepresentation.SetMatrixElement(i,j,Tmp);
	}
    }
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeDensityMatrix" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "redirect output to this file",NULL);
  (*OutputGroup) += new BooleanOption ('\n', "plot-density", "plot the density matrix eigenstates");
  (*OutputGroup) += new BooleanOption ('\n', "save-vectors", "write vectors, which yield maximum density matrix EV");  
  (*OutputGroup) += new BooleanOption  ('v', "momentum-vectors", "writes the basis of momentum eigenstates");  
  (*OutputGroup) += new SingleDoubleOption  ('r',"dynamic-range","range of density operator eigenvalues to be displayed",1e-5);
  (*OutputGroup) += new BooleanOption  ('\n', "show-translation", "display the matrix defining the translation operator");
  (*OutputGroup) += new BooleanOption  ('\n', "show-basis", "show elements of vector in basis and exit");
  (*MiscGroup) += new SingleIntegerOption ('\n', "nbr-density", "number of density matrix eigenstates to be written out",1);
  (*MiscGroup) += new SingleIntegerOption ('s',"superpositions","in case of two input vectors, number of values for phase in superpositions",12);
  (*MiscGroup) += new BooleanOption  ('V', "verbose", "give additional output");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  int NbrVectors;
  char** VectorFiles = Manager.GetStrings("states",NbrVectors);

  if (NbrVectors==0)
    {
      cout << "At least one vector file is required!"<<endl;
      exit(1);
    }
  double Interaction=-1.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[0], NbrBosons, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }  
  HardCore=(HardCore||Manager.GetBoolean("hard-core"));
  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;
  
  int NbrSites = Lx*Ly;
  int VectorDimension=0;
  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  bool tmpB, haveVector=false;
  for (int i=0; i<NbrVectors; ++i)
    {
      tmpB = Vectors[i].ReadVector(VectorFiles[i]);
      if (!haveVector)
	VectorDimension=Vectors[i].GetVectorDimension();
      if (haveVector && (Vectors[i].GetVectorDimension()!=VectorDimension))
	{
	  cout<<"Dimension of vector "<<VectorFiles[i]<<" does not match size of previous vectors!"<<endl;
	  exit(1);
	}
      haveVector=haveVector | tmpB;
	    
//       cout << "Vector "<<i<<":"<<endl;
//       for (int j=0; j<VectorDimension; ++j)
// 	{
// 	  cout<<Vectors[i][j]<<" ( ";
// 	  Space->PrintState(cout,j);
// 	  cout << " )"<<endl;
// 	}
//       int dX, dY;
//       for (int fi=0; fi<VectorDimension; ++fi)
// 	if (Space->IsTranslation(0, fi, dX, dY))
// 	    cout << "Potential translation " << 0 << "->"<<fi<<endl;
    }

  if (!haveVector)
    {
      cout << "No valid vector files found!"<<endl;
      exit(1);
    }  

  ParticleOnLattice* Space;
  if (HardCore)
    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);

  if (VectorDimension != Space->GetHilbertSpaceDimension())
    {
      cout<<"Dimension of vectors does not match size of Hilbert-space!"<<endl;
	  exit(1);
    }

  if (Manager.GetBoolean("show-basis"))
    {
      for (int v=0; v<NbrVectors; ++v)
	{
	  cout << "Components of vector "<<v<<":"<<endl;
	  for (int i=0; i<VectorDimension; ++i)
	    {
	      Space->PrintState(cout, i);
	      cout <<" :  "<<Vectors[v][i]<<endl;
	    }
	}
      exit(0);
    }
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  ParticleOnLatticeOneBodyOperator *DensityOperator= new ParticleOnLatticeOneBodyOperator(Space);
  ParticleOnLatticeTranslationOperator *TranslationOperator= new ParticleOnLatticeTranslationOperator(Space);

  cout<< "========= Analysis of density matrix ========"<<endl;
  
  int DensityMatrixDimension = NbrSites*NbrVectors;
  HermitianMatrix Rho(DensityMatrixDimension);  

  Complex Tmp;
  int CreationIndex, AnnihilationIndex, TotalIndexI, TotalIndexJ;  
  for (int CreationX=0; CreationX<Lx; ++CreationX)
    for (int CreationY=0; CreationY<Ly; ++CreationY)
      {
	CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, 0, Tmp);	
	for (int AnnihilationX=0; AnnihilationX<Lx; ++AnnihilationX)
	  for (int AnnihilationY=0; AnnihilationY<Ly; ++AnnihilationY)
	    {
	      AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, 0, Tmp);
	      DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
	      // calculate possible matrix elements in subspace of vectors
	      for (int numVector=0; numVector<NbrVectors; ++numVector)
		for (int numVector2=0; numVector2<NbrVectors; ++numVector2)
		  {
		    TotalIndexI = CreationIndex+numVector*NbrSites;
		    TotalIndexJ = AnnihilationIndex+numVector2*NbrSites;
		    if (TotalIndexI<=TotalIndexJ)
		      {
			Tmp=DensityOperator->MatrixElement(Vectors[numVector], Vectors[numVector2]);
			Rho.SetMatrixElement(TotalIndexI,TotalIndexJ,Tmp);
		      }
		  }	      
	    }
      }  
  
  //cout << "Matrix="<<endl<<Rho<<endl;
  // calculate eigenvalues & vectors of Rho
  double dynamics = Manager.GetDouble("dynamic-range");
  RealDiagonalMatrix M;
  ComplexMatrix Q;  
  Rho.Diagonalize(M, Q, 1e-12, 1000);
  for (int i=0; i<DensityMatrixDimension; ++i)
    if (fabs(M[DensityMatrixDimension-1-i])>dynamics*M[DensityMatrixDimension-1])
      cout << "EV["<<i<<"] = " << M[DensityMatrixDimension-1-i] << endl;
  //cout << "Transition Matrix: "<<endl<<Q<<endl;
  cout << "First Eigenvector: "<<endl;
  Complex TmpC;
  for (int i=0; i<DensityMatrixDimension; ++i)
    {      
      Q.GetMatrixElement(i,DensityMatrixDimension-1,TmpC);
      cout << TmpC.Re << "+I*" << TmpC.Im << endl;
    }

  if (NbrVectors==1)
    {
      // write eigenstate of density matrix in basis of Hilbert-space
      char *RhoVecOut = new char[strlen(VectorFiles[0])+10];
      ComplexVector EigenState(DensityMatrixDimension);
      for (int s=0; s<Manager.GetInteger("nbr-density"); ++s)
	{
	  sprintf(RhoVecOut,"%s.dm%d",VectorFiles[0],s);
      
	  for (int i=0; i<DensityMatrixDimension; ++i)
	    {
	      Q.GetMatrixElement(i,DensityMatrixDimension-1-s,EigenState[DensityMatrixDimension-1-i]);
	    }
	  EigenState.WriteVector(RhoVecOut);
	  if (Manager.GetBoolean("plot-density"))
	    {
	      sprintf(RhoVecOut,"%s.dm%d.dat",VectorFiles[0],s);
	      ofstream DataFile(RhoVecOut);
	      DataFile << "# X\tY\tv_x\tv_y"<<endl;
	      for (int x=0; x<Lx; ++x)
		for (int y=0; y<Ly; ++y)
		  {
		    int q = Space->EncodeQuantumNumber(x, y, 0, Tmp);	
		    DataFile << x<<"\t"<<y<<"\t"<<EigenState[DensityMatrixDimension-1-q].Re<<"\t"<<EigenState[DensityMatrixDimension-1-q].Im<<endl;
		  }
	      DataFile.close();
	    }
		
	}
      delete [] RhoVecOut;
    }  

  if (NbrVectors==2)
    {
      int DensityMatrixDimension2 = NbrSites;
      RealDiagonalMatrix M2;	  
      HermitianMatrix Rho2(DensityMatrixDimension2);
      cout << "====== Analysing superpositions of form |1> + e^(i phi) |2> ======" << endl;
      ComplexVector Superposition = ComplexVector(Vectors[0].GetVectorDimension());
      ComplexVector *EigenStates = new ComplexVector[Manager.GetInteger("superpositions")];      
      int Max1=-1, Max2=-1;
      double MaxVal1=0.0, MaxVal2=0.0;      
      
      for (int k=0; k<Manager.GetInteger("superpositions");++k)
	{	  
	  Complex Phase = Polar(sqrt(0.5),(2.0*M_PI*k)/Manager.GetInteger("superpositions"));
	  Superposition.Copy(Vectors[0],sqrt(0.5));
	  Superposition.AddLinearCombination (Phase, Vectors[1]);
	  for (int CreationX=0; CreationX<Lx; ++CreationX)
	    for (int CreationY=0; CreationY<Ly; ++CreationY)
	      {
		CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, 0, Tmp);	
		for (int AnnihilationX=0; AnnihilationX<Lx; ++AnnihilationX)
		  for (int AnnihilationY=0; AnnihilationY<Ly; ++AnnihilationY)
		    {
		      AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, 0, Tmp);
		      DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
		      // calculate possible matrix elements in subspace of vectors
		      if (CreationIndex <= AnnihilationIndex)
			{
			  Tmp=DensityOperator->MatrixElement(Superposition, Superposition);
			  Rho2.SetMatrixElement(CreationIndex, AnnihilationIndex, Tmp);
			}
		    }
	      }
	  Rho2.Diagonalize(M2, Q, 1e-10, 250);
	  cout << "EV's["<<2*k<<"/"<<Manager.GetInteger("superpositions")<<"pi] = " << M2[DensityMatrixDimension2-1] << ", "
	       <<M2[DensityMatrixDimension2-2] <<", "<<M2[DensityMatrixDimension2-3]<<endl;
	  EigenStates[k].Resize(DensityMatrixDimension2);
	  for (int i=0; i<DensityMatrixDimension2; ++i)
	    Q.GetMatrixElement(i,DensityMatrixDimension2-1,EigenStates[k][DensityMatrixDimension2-1-i]);
	  if (M2[DensityMatrixDimension2-1]>=MaxVal1-1e-8)
	    {
	      MaxVal2=MaxVal1;
	      Max2=Max1;
	      MaxVal1 = M2[DensityMatrixDimension2-1];
	      Max1=k;
	    }
	}
      
      char *RhoVecOut = new char[strlen(VectorFiles[0])+10];
      
      if (Max1>=0)
	{
	  sprintf(RhoVecOut,"%s-2a.dm",VectorFiles[0]);
	  cout << "Writing density matrix eigenstate for EV["<<2*Max1<<"/"<<Manager.GetInteger("superpositions")<<"pi] to "<<RhoVecOut<<endl;
	  EigenStates[Max1].WriteVector(RhoVecOut);
	  if (Manager.GetBoolean("save-vectors"))
	    {	      
	      sprintf(RhoVecOut,"%s-2a.vec",VectorFiles[0]);
	      cout << "Writing superposition Vec["<<2*Max1<<"/"<<Manager.GetInteger("superpositions")<<"pi] to "<<RhoVecOut<<endl;
	      Complex Phase = Polar(sqrt(0.5),(2.0*M_PI*Max1)/Manager.GetInteger("superpositions"));
	      Superposition.Copy(Vectors[0],sqrt(0.5));
	      Superposition.AddLinearCombination (Phase, Vectors[1]);
	      Superposition/=Superposition.Norm();
	      Superposition.WriteVector(RhoVecOut);
	    }
	}
      if (Max2>=0)
	{
	  sprintf(RhoVecOut,"%s-2b.dm",VectorFiles[0]);
	  cout << "Writing density matrix eigenstate for EV["<<2*Max2<<"/"<<Manager.GetInteger("superpositions")<<"pi] to "<<RhoVecOut<<endl;
	  EigenStates[Max2].WriteVector(RhoVecOut);
	  if (Manager.GetBoolean("save-vectors"))
	    {	      
	      sprintf(RhoVecOut,"%s-2b.vec",VectorFiles[0]);
	      cout << "Writing superposition Vec["<<2*Max2<<"/"<<Manager.GetInteger("superpositions")<<"pi] to "<<RhoVecOut<<endl;
	      Complex Phase = Polar(sqrt(0.5),(2.0*M_PI*Max2)/Manager.GetInteger("superpositions"));
	      Superposition.Copy(Vectors[0],sqrt(0.5));
	      Superposition.AddLinearCombination (Phase, Vectors[1]);
	      Superposition/=Superposition.Norm();
	      Superposition.WriteVector(RhoVecOut);
	    }
	}
      delete [] EigenStates;
      delete [] RhoVecOut;
    }

  if (NbrVectors>1)
    {
      int DensityMatrixDimension2 = NbrSites;
      RealDiagonalMatrix M2;	  
      HermitianMatrix Rho2(DensityMatrixDimension2);
      cout << "====== Analysing sum of density matrices ======" << endl;
      for (int CreationX=0; CreationX<Lx; ++CreationX)
	for (int CreationY=0; CreationY<Ly; ++CreationY)
	  {
	    CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, 0, Tmp);	
	    for (int AnnihilationX=0; AnnihilationX<Lx; ++AnnihilationX)
	      for (int AnnihilationY=0; AnnihilationY<Ly; ++AnnihilationY)
		{
		  AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, 0, Tmp);
		  DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
		  // calculate possible matrix elements in subspace of vectors
		  // if (CreationIndex <= AnnihilationIndex)
		  Tmp=0.0;
		  for (int i=0; i<NbrVectors; ++i)
		    Tmp+=DensityOperator->MatrixElement(Vectors[i], Vectors[i]);
		  Rho2.SetMatrixElement(CreationIndex, AnnihilationIndex, Tmp);
		}
	  }
      Rho2.Diagonalize(M2, 1e-10, 250);      
      for (int i=0; i<DensityMatrixDimension2; ++i)
	if (fabs(M2[DensityMatrixDimension2-1-i])
	    >dynamics*M2[DensityMatrixDimension2-1])
	  cout << "Sum-EV["<<i<<"] = " << M2[DensityMatrixDimension2-1-i] << endl;
    }


  cout<< "====== Analysis of momentum eigenvalues ====="<<endl;

  ComplexVector TmpState(VectorDimension);
  ComplexVector TmpState2(VectorDimension);


  if (Manager.GetBoolean("show-translation"))
  {
    // testing unitarity of translation operator matrix and display it:
    ComplexMatrix TrRep(VectorDimension, VectorDimension);  
    
    TranslationOperator->SetTranslationComponents(1,0);
    for (int i=0; i<VectorDimension; ++i)
      {
	TmpState2.ClearVector();
	TmpState2.Re(i)=1.0;
	VectorOperatorMultiplyOperation Operation (TranslationOperator, &TmpState2, &TmpState);      
	Operation.ApplyOperation(Architecture.GetArchitecture());      
	for (int j=0; j<VectorDimension; ++j)
	  TrRep.SetMatrixElement(j,i,TmpState[j]);
      }
    
    cout << "Representation of T_x"<<endl<<TrRep<<endl;
  }

  
  ComplexMatrix XTranslationMatrix(NbrVectors, NbrVectors);
  ComplexMatrix YTranslationMatrix(NbrVectors, NbrVectors);
  ComplexVector TmpState3(VectorDimension);

  int Degeneracy=1;
  int n1=1, n2=1;
  int FluxModulo = FindGCD(NbrFluxQuanta, Lx*Ly);
  int r=NbrFluxQuanta/FluxModulo;
  int t=Lx*Ly/FluxModulo;

  while ((((Ly*n1)%t)!=0) && (n1<Lx)) ++n1;
  while ((((Lx*n2)%t)!=0) && (n2<Ly)) ++n2;

  if ((Lx%n1)!=0)
    cout << "Extending range of n1 to Lx"<<endl;
  if ((Ly%n2)!=0)
    cout << "Extending range of n2 to Ly"<<endl;

  if (((n1*n2*NbrFluxQuanta)%t) != 0)
    {
      cout << "Cannot resolve translations: Brillouin zone trivial?"<<endl;
      n1=Lx;
      n2=Ly;
    }

  while ((r*NbrBosons*n1*n2*Degeneracy)%t != 0) ++Degeneracy;
  
  cout << "N_phi = "<<r<<"/"<<t<<endl;
  cout << "n1="<<n1<<", n2="<<n2<<", global degeneracy: "<<Degeneracy<<endl;

  int RemainingDegeneracy=Degeneracy;

  ComplexMatrix EVecX(NbrVectors, NbrVectors);
  ComplexMatrix EVecY(NbrVectors, NbrVectors);  
  ComplexDiagonalMatrix EValX(NbrVectors, NbrVectors);
  ComplexDiagonalMatrix EValY(NbrVectors, NbrVectors);
  
  TranslationOperator->SetTranslationComponents(n1,0);
  GetTranslationMatrix(TranslationOperator, NbrVectors, Vectors, XTranslationMatrix, TmpState, Architecture);

  XTranslationMatrix.Diagonalize(EValX);  
  if ((fabs(Norm(EValX[0])-1.0)>1e-10)||((Ly/n2)<Degeneracy))
    {
      int GCD = FindGCD(Lx/n1, Degeneracy);
      RemainingDegeneracy/=GCD;
      if (GCD!=1) cout << "Multiplying factor "<<GCD<<" of degeneracy onto n1"<<endl;
      n1*=GCD;
      TranslationOperator->SetTranslationComponents(n1,0);
      GetTranslationMatrix(TranslationOperator, NbrVectors, Vectors, XTranslationMatrix, TmpState, Architecture);
    }
  
  if ((Ly/n2)%RemainingDegeneracy!=0)
    {
      cout<<"Did not treat degeneracy properly -> need to put onto n1?"<<endl;
      exit(1);      
    }
  else
    {
      if (RemainingDegeneracy!=1)
	cout << "Multiplying factor "<<RemainingDegeneracy<<" of degeneracy onto n2"<<endl;
      n2*=RemainingDegeneracy;
      RemainingDegeneracy=1;
    }
      
  TranslationOperator->SetTranslationComponents(0,n2);
  GetTranslationMatrix(TranslationOperator, NbrVectors, Vectors, YTranslationMatrix, TmpState, Architecture);

  if (Manager.GetBoolean("show-translation"))
    {
      cout << "XTranslationMatrix="<<endl<<XTranslationMatrix<<endl;
      cout << "YTranslationMatrix="<<endl<<YTranslationMatrix<<endl;
    }

  ComplexMatrix EVecXY(NbrVectors, NbrVectors);


  // form linear superposition of Tx and Ty to diagonalize:
  ComplexMatrix Z((Matrix&)XTranslationMatrix);
  Z*=log(91.0); // scale with some random number > 1
  Z+=YTranslationMatrix;
  Z.Diagonalize(EValX,EVecXY);
  ComplexMatrix QH=EVecXY.GetAdjoint();

  bool IsDiagonal;
      
  ComplexDiagonalMatrix XEV(EVecXY.GetAdjoint()*(XTranslationMatrix*EVecXY),IsDiagonal, 1e-6);
  
  if (IsDiagonal)
    {
      if (Manager.GetBoolean("verbose"))
	cout << "EigenValues(Tx)="<<endl<<XEV<<endl;
    }
  else
    cout << "EigenValues(Tx)=  !!! Attention, was not fully diagonal !!!"
	 <<endl<<EVecXY.GetAdjoint()*(XTranslationMatrix*EVecXY)<<endl;

  ComplexDiagonalMatrix YEV(EVecXY.GetAdjoint()*(YTranslationMatrix*EVecXY),IsDiagonal, 1e-6);
  if (IsDiagonal)
    {
      if (Manager.GetBoolean("verbose"))
	cout << "EigenValues(Ty)="<<endl<<YEV<<endl;
    }
  else
    cout << "EigenValues(Ty)=  !!! Attention, was not fully diagonal !!!"
	 <<endl<<EVecXY.GetAdjoint()*(YTranslationMatrix*EVecXY)<<endl;
  
  if (Manager.GetBoolean("verbose"))
    cout << "Eigenvectors="<<endl<<EVecXY<<endl;

  cout << "#i\tKx\tKy"<<endl;
  for (int i=0; i<NbrVectors; ++i)
    {
      cout <<i<<"\t"<<Arg(XEV[i])/M_PI<<"\t"<<Arg(YEV[i])/M_PI;
      if (fabs(Norm(XEV[i])-1.0)>1e-10) cout << "\t!!abs(Tx)="<<Norm(XEV[i]);
      if (fabs(Norm(YEV[i])-1.0)>1e-10) cout << "\t!!abs(Ty)="<<Norm(YEV[i]);
      cout << endl;
    }

  if (Manager.GetBoolean("momentum-vectors"))
    {
      char *vectorName=new char [strlen(VectorFiles[0])+20];
      strcpy(vectorName,VectorFiles[0]);
      int endBase=strlen(vectorName)-1;
      int countDot=0;
      for (;(endBase>=0)&&(countDot<2);--endBase)
	if (vectorName[endBase]=='.') ++countDot;
      if (countDot<2)
	{
	  countDot=0;
	  for (;(endBase>=0)&&(countDot<1);--endBase)
	    if (vectorName[endBase]=='.') ++countDot;
	}
      endBase++;
      int nbrVec;
      int minNbrVec=1000;
      int maxNbrVec=-1;
      for (int i=0;i<NbrVectors;++i)
	{
	  sscanf(VectorFiles[i]+endBase+1,"%d.vec",&nbrVec);
	  if (nbrVec>maxNbrVec) maxNbrVec = nbrVec;
	  if (nbrVec<minNbrVec) minNbrVec = nbrVec;
	  //cout << "Number of vector="<<nbrVec<<" char " <<(char)('A'+nbrVec)<<endl;
	}
      if ((minNbrVec>=0) && (minNbrVec<26)&&(maxNbrVec<26)&&(maxNbrVec>=0))
	VectorFiles[0][endBase]='\0';
      else
	{
	  minNbrVec=0;
	  VectorFiles[0][strlen(VectorFiles[0])-4]='\0';
	}
      for (int i=0;i<NbrVectors;++i)
	{
	  sprintf(vectorName,"%s.%c.vec",VectorFiles[0],'a'+minNbrVec+i);
	  TmpState.ClearVector();
	  for (int j=0; j<NbrVectors;++j)
	    TmpState.AddLinearCombination(Conj(EVecXY[i][j]),Vectors[j]);
	  TmpState/=TmpState.Norm();
	  cout << "Vector-"<<i<<"="<<vectorName<<endl;
	  TmpState.WriteVector(vectorName);
	}
    }
      
  delete Space;
  delete [] Vectors;
  delete DensityOperator;
  delete TranslationOperator;
}
