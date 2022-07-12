#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"

#include "HilbertSpace/GenericSUNSpinCollection.h"

#include "Hamiltonian/SUNSpinOnLatticeQuadraticHamiltonian.h"

#include "Tools/LatticeConnections.h"


#include "Operator/ChiralExchangeOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MathTools/h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
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

  OptionManager Manager ("SUNSpinsOnLAtticeCorrelations" , "0.01");  
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);  
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "state", "filename of state vector to be processed");
  (*SystemGroup) += new SingleIntegerOption  ('N', "level-n", "level of SU(N) symmetry group (0=parse from file-name)", 0);
  (*SystemGroup) += new MultipleIntegerOption  ('t', "cartan", "eigenvalues of the generators of the cartan algebra (omit = parse)",',');
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "redirect output to this file",NULL);

  (*MiscGroup) += new BooleanOption  ('\n', "show-basis", "show elements of input vector in given many-body basis");  
  (*MiscGroup) += new BooleanOption  ('V', "verbose", "give additional output");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);

  int LevelN=Manager.GetInteger("level-n");  

  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

    if (Manager.GetString("lattice-definition")==NULL)
    {
      cout << "Please indicate the file with the lattice-definition for this vector"<<endl;
      exit(1);
    }
  // get the lattice geometry
  LatticeConnections *Lattice = new LatticeConnections();

  char *InputFile = Manager.GetString("state");
  
  NbrSpins = Lattice->GetNbrSites();
  int Lx = Lattice->GetLatticeLength(0);
  int Ly = Lattice->GetLatticeLength(1);
  NbrSubLattices = Lattice->GetNbrSubLattices();
  char* LatticeName = Lattice->GeometryString();
  if (strstr(InputFile, LatticeName)==0)
	{
	  cout << "The given lattice parameters do not coincide with the filename, verify lattice definition, and repetition of unit cells"<<endl;
	  exit(1);
	}
  delete [] LatticeName;
  char *Descriptor=0;
  int FileNbrSpins=0;
  int NbrCartan;
  int *CartanQuantumNumbers;
  CartanQuantumNumbers = Manager.GetIntegers("cartan",NbrCartan);

  if (SUNSpinFindSystemInfoFromFileName(InputFile, LevelN, Descriptor, FileNbrSpins, CartanQuantumNumbers) == false)
    {
      cout<<"Please use standard file-names, or indicate all necessary system parameters!"<<endl;
      exit(1);
    }
  
  int VectorDimension=0;
  RealVector InputVector;
  
  InputVector.ReadVector(InputFile);
  VectorDimension=InputVector.GetVectorDimension();

  GenericSUNSpinCollection *Space = new GenericSUNSpinCollection(LevelN, NbrSpins, CartanQuantumNumbers, MemorySpace);

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

  ChiralEchangeOperator *ChiralExchange= new ChiralEchangeOperator(Space);
  

  cout<< "========= Analysis of chiral spin exchange ========"<<endl;
  
  int DensityMatrixDimension = NbrSites*NbrVectors;
  HermitianMatrix Rho(DensityMatrixDimension);  

  Complex Tmp;
  int CreationIndex, AnnihilationIndex, TotalIndexI, TotalIndexJ;  
  for (int CreationX=0; CreationX<Lx; ++CreationX)
    for (int CreationY=0; CreationY<Ly; ++CreationY)
      for (int CreationSub=0; CreationSub<NbrSubLattices; ++CreationSub)
	{
	  CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, CreationSub, Tmp);
	  for (int AnnihilationX=0; AnnihilationX<Lx; ++AnnihilationX)
	    for (int AnnihilationY=0; AnnihilationY<Ly; ++AnnihilationY)
	      for (int AnnihilationSub=0; AnnihilationSub<NbrSubLattices; ++AnnihilationSub)
		{
		  AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, AnnihilationSub, Tmp);
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
      RealVector EigenState(DensityMatrixDimension);
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
      RealVector Superposition = RealVector(Vectors[0].GetVectorDimension());
      RealVector *EigenStates = new RealVector[Manager.GetInteger("superpositions")];      
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

  
  // stop here if we have a generic lattice -> translations as of yet not implemented
  if (GenericLattice) exit(0);    

  cout<< "====== Analysis of momentum eigenvalues ====="<<endl;

  RealVector TmpState(VectorDimension);
  RealVector TmpState2(VectorDimension);


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
  RealVector TmpState3(VectorDimension);

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
