#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeGeneric.h"
#include "HilbertSpace/HardCoreBosonOnLatticeGeneric.h"

#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHESpectrum/LatticePhases.h"
#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"
#include "Tools/FQHEWaveFunction/MaximallyCondensedStateOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"
#include "Operator/ParticleOnLatticeMomentumOperator.h"
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
  LatticePhases::AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);  
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new MultipleStringOption  ('\0', "states", "filenames of state vectors to be processed");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice", 0);
  (*SystemGroup) += new SingleDoubleOption  ('Q', "cont-flux", "multiples of flux quanta piercing the lattice", 0.0);
  (*SystemGroup) += new MultipleDoubleOption  ('s', "solenoid-flux", "twist in periodic boundary conditions phi_x[,phi_y])",',');
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
  (*OutputGroup) += new BooleanOption  ('\n', "currents", "calculate expectation values of current operators and exit");
  (*SystemGroup) += new BooleanOption  ('\n', "normalize-phase", "change phase to make largest absolute value real and positive");
  (*OutputGroup) += new SingleDoubleOption  ('\n', "display-threshold", "only display values larger than threshold",0.0);
  (*MiscGroup) += new SingleIntegerOption ('\n', "nbr-density", "number of density matrix eigenstates to be written out",1);
  (*MiscGroup) += new BooleanOption  ('\n', "joint", "evaluate joint density matrix for multiple vectors");
  (*MiscGroup) += new BooleanOption  ('\n', "equal", "only use equal weight superpositions of two vectors");
  (*MiscGroup) += new SingleIntegerOption ('u',"superpositions","in case of two input vectors, number of values for phase in superpositions",12);
  (*MiscGroup) += new SingleIntegerOption ('\n',"max-iter","maximum number of iterations for optimizing condensate fraction",250);
  (*MiscGroup) += new SingleDoubleOption ('\n',"opt-tolerance","tolerance for optimizing condensate fraction",1e-4);
  (*MiscGroup) += new BooleanOption  ('\n', "opt-random", "start optimization from a randomized initial condition");
  (*MiscGroup) += new SingleIntegerOption  ('\n', "opt-target", "target for optimization on sum over first N eigenvalues",1);
  (*MiscGroup) += new BooleanOption  ('\n', "expansion", "obtain expansion image of state(s)");
  (*MiscGroup) += new SingleIntegerOption ('\n',"apply-gauge","Apply a gauge-transform before calculating expansion image - 1:doubled->Dali");
  (*MiscGroup) += new MultipleDoubleOption  ('\n', "exp-offset", "apply offset of momentum values for Fourier-Transform O_x,O_y)",',');
  (*MiscGroup) += new BooleanOption  ('\n', "no-momenta", "do not calculate momentum of states");
  (*MiscGroup) += new BooleanOption  ('V', "verbose", "give additional output");  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  double SolenoidX=0.0, SolenoidY=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-flux", tmpI);
    if (tmpI>0) SolenoidX=Fluxes[0];
    if (tmpI>1) SolenoidY=Fluxes[1];
    if (tmpI>0) delete [] Fluxes;
  }

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
  bool GenericLattice=false;
  int NbrSites=0;
  int NbrSubLattices=1;
  bool HaveContFlux = false;
  double ContFlux =0.0;
  LatticePhases *Lattice = NULL;
  if ((Manager.GetString("lattice-definition")!=NULL)||(FQHEOnLatticeHaveGeneralLattice(VectorFiles[0])))
    {
      GenericLattice=true;
      if (Manager.GetString("lattice-definition")==NULL)
	{
	  cout << "Please indicate the file with the lattice-definition for this vector"<<endl;
	  exit(1);
	}
      // get the lattice geometry
      Lattice = new LatticePhases();
      NbrSites = Lattice->GetNbrSites();
      Lx = Lattice->GetLatticeLength(0);
      Ly = Lattice->GetLatticeLength(1);
      NbrSubLattices = Lattice->GetNbrSubLattices();
      char* LatticeName = Lattice->GeometryString();
      if (strstr(VectorFiles[0], LatticeName)==0)
	{
	  cout << "The given lattice parameters do not coincide with the filename, verify lattice definition, and repetition of unit cells"<<endl;
	}
      delete [] LatticeName;
      if (FQHEOnLatticeFindSystemInfoFromGeneralVectorFileName(VectorFiles[0], NbrBosons, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore, HaveContFlux, ContFlux) == false)
	{
	  cout<<"Please use standard file-names, or indicate all necessary system parameters!"<<endl;
	  exit(1);
	}
      
    }
  else
    {
      if (FQHEOnLatticeFindSystemInfoFromVectorFileName(VectorFiles[0], NbrBosons, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
	{
	  cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
	  exit(1);
	}
      HardCore=(HardCore||Manager.GetBoolean("hard-core"));
      if (Manager.GetBoolean("no-hard-core"))
	HardCore=false;
      NbrSites = Lx*Ly;
    }
  
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
    }

  if (!haveVector)
    {
      cout << "No valid vector files found!"<<endl;
      exit(1);
    }  

  ParticleOnLattice* Space;
  if (GenericLattice)
    {
      if (HardCore)
	Space = new HardCoreBosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
      else Space = new BosonOnLatticeGeneric(NbrBosons, Lattice, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
    }
  else
    {
      if (HardCore)
	Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
      else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace, SolenoidX, SolenoidY);
    }

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
	  if (Manager.GetBoolean("normalize-phase"))
	    {
	      double MaxNorm=0.0;
	      int MaxIndex=0;
	      for (int i = 0; i < Vectors[v].GetVectorDimension();++i)
		{
		  if (Norm(Vectors[v][i])>MaxNorm)
		    {
		      MaxNorm=Norm(Vectors[v][i]);
		      MaxIndex=i;
		    }
		}
	      Vectors[v] *= Polar(1.0,-Arg(Vectors[v][MaxIndex]));
	    }
	  if (Manager.GetDouble("display-threshold")>0.0)
	    {
	      double Threshold=Manager.GetDouble("display-threshold");
	      for (int i=0; i<VectorDimension; ++i)
		{
		  if (Norm(Vectors[v][i])>Threshold)
		    {
		      Space->PrintState(cout, i);
		      cout.precision(5);
		      cout <<" :  "<<Vectors[v][i];
		      cout <<" -- "<<Norm(Vectors[v][i])<<" "<<Arg(Vectors[v][i])/M_PI<<endl;
		      cout.precision(14);
		    }
		}
	      }
	  else
		
	    for (int i=0; i<VectorDimension; ++i)
	      {
		Space->PrintState(cout, i);
		cout.precision(5);
		cout <<" :  "<<Vectors[v][i].Re;
		if (Vectors[v][i].Im>0)
		  cout << "+"<<Vectors[v][i].Im<<"I";
		else
		  cout << "-"<<-Vectors[v][i].Im<<"I";
		cout <<" -- "<<Norm(Vectors[v][i])<<" "<<Arg(Vectors[v][i])/M_PI<<endl;
		cout.precision(14);
	      }
	}
      exit(0);
    }
  
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  if (Manager.GetBoolean("currents"))
    {
      if (GenericLattice==false)
	{
	  cout << "Calculation of currents not implemented for default square lattice - use a configuration file"<<endl;
	  exit(1);
	}

      char* CurrentFieldName = ReplaceExtensionToFileName(VectorFiles[0],"vec","j.dat");
	    
      ofstream VectorField;
      VectorField.open(CurrentFieldName,ios::out);
      delete [] CurrentFieldName;
      double FluxDensity;
      if (NbrFluxQuanta>Lx*Ly)
	{
	  FluxDensity=0.5;
	  cout << "Assuming flux density of n_phi="<<FluxDensity<<endl;
	}
      else
	{
	  FluxDensity=(double)NbrFluxQuanta/(double)NbrSites;
	}
      // calculate current operator
      RealVector SitePosition, SitePosition2, RelativePosition;
      int CellPosition[2];
      int Sub;

      int NbrNeighbors;
      int *Neighbors;
      double *Phases;
      double *Amplitudes;
      int **PeriodicTranslations;
      ParticleOnLatticeOneBodyOperator *DensityOperatorOut= new ParticleOnLatticeOneBodyOperator(Space);
      ParticleOnLatticeOneBodyOperator *DensityOperatorIn= new ParticleOnLatticeOneBodyOperator(Space);
      for (int s=0; s<NbrSites; ++s)
	{
	  Lattice->GetSiteCoordinates(s, CellPosition, Sub);
	  SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
	  Lattice->GetNeighbors(s, NbrNeighbors, Neighbors, Phases, PeriodicTranslations, Amplitudes);
	  ComplexVector TargetState(Space->GetHilbertSpaceDimension()), TargetState2(Space->GetHilbertSpaceDimension());
	  for (int n=0; n<NbrNeighbors; ++n)
	    if (Neighbors[n]>s)
	      {
		DensityOperatorOut->SetCreationAnnihilationIndex(Neighbors[n],s);
		DensityOperatorIn->SetCreationAnnihilationIndex(s,Neighbors[n]);

		double *Current = new double[NbrNeighbors];
		for (int v=0; v<NbrVectors; ++v)
		  {
		    VectorOperatorMultiplyOperation Operation (DensityOperatorIn, &(Vectors[v]), &TargetState);
		    Operation.ApplyOperation(Architecture.GetArchitecture());
		    VectorOperatorMultiplyOperation Operation2 (DensityOperatorOut, &(Vectors[v]), &TargetState2);
		    Operation2.ApplyOperation(Architecture.GetArchitecture());

		    Complex PhaseIn = Vectors[v] * TargetState * Polar(1.0, -2.0*M_PI*FluxDensity*Phases[n]);

		    Complex PhaseOut = Vectors[v] * TargetState2 * Polar(1.0, 2.0*M_PI*FluxDensity*Phases[n]);

		    // correct expression for trial states:
		    // double Current = -2.0*Imag(Conj(ResultingState[NbrSites-s-1])*ResultingState[NbrSites-Neighbors[n]-1]*Polar(1.0, -2.0*M_PI*FluxDensity*Phases[n]));
		    Current[n] = - Imag(PhaseIn + PhaseOut);
		    if (Amplitudes!=NULL)
		      Current[n] *= Amplitudes[n];
		  }
		//cout << "J_"<<s<<","<<Neighbors[n]<<" = "<<Current<<endl;
		Lattice->GetSiteCoordinates(Neighbors[n], CellPosition, Sub);
		SitePosition2 = Lattice->GetSitePosition(CellPosition,Sub);
		
		for (int d=0; d<2; ++d)
		  SitePosition2.AddLinearCombination((double)(-PeriodicTranslations[n][d]*Lattice->GetLatticeLength(d)),Lattice->GetLatticeVector(d));
		//cout << "Connection "<<endl<<SitePosition<<"to"<<endl<<SitePosition2<< "(shift "<<PeriodicTranslations[n][0]<<","<< PeriodicTranslations[n][1]<<")"<<endl;
		RelativePosition.Copy(SitePosition2);
		RelativePosition.AddLinearCombination(-1.0,SitePosition);
		RelativePosition/=RelativePosition.Norm();
		SitePosition2.AddLinearCombination(1.0,SitePosition);
		SitePosition2*=0.5;
		
		VectorField << SitePosition2[0] << "\t" << SitePosition2[1]
			    << "\t" << RelativePosition[0]*Current[0] << "\t" << RelativePosition[1]*Current[0];

		for (int v=1; v<NbrVectors; ++v)
		  VectorField << "\t" << RelativePosition[0]*Current[v] << "\t" << RelativePosition[1]*Current[v];
		VectorField << endl;
		delete [] Current;
	      }
	}
      delete DensityOperatorOut;
      delete DensityOperatorIn;
    }

  ParticleOnLatticeOneBodyOperator *DensityOperator= new ParticleOnLatticeOneBodyOperator(Space);

  cout<< "========= Analysis of ";
  if (NbrVectors>1) cout << "extended ";
  cout <<"density matrix ========"<<endl;
 
  int DensityMatrixDimension = NbrSites*NbrVectors;
  HermitianMatrix Rho(DensityMatrixDimension);  
  Complex Tmp;
  int CreationIndex, AnnihilationIndex, TotalIndexI, TotalIndexJ;  
  double dynamics = Manager.GetDouble("dynamic-range");
  RealDiagonalMatrix M;
  ComplexMatrix Q;  

  if ((NbrVectors==1)||(Manager.GetBoolean("joint")))
    {
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

      if (NbrVectors > 1)
	{
	  ComplexMatrix Block1(NbrSites, NbrSites);
	  ComplexMatrix Block2(NbrSites, NbrSites);
	  ComplexMatrix Rst1(NbrSites, NbrSites);
	  ComplexMatrix Rst2(NbrSites, NbrSites);
	  for (int i=0; i<NbrVectors; ++i)
	    {
	      for (int k=0; k<NbrSites; ++k)
		for (int l=0; l<NbrSites; ++l)
		  Block1.SetMatrixElement(k,l,Rho(i*NbrSites+k,i*NbrSites+l));
	      for (int j=0; j<NbrVectors; ++j)
		{
		  if (j==i) ++j;
		  if (j>=NbrVectors) break;
		  for (int k=0; k<NbrSites; ++k)
		    for (int l=0; l<NbrSites; ++l)
		      Block2.SetMatrixElement(k,l,Rho(i*NbrSites+k,j*NbrSites+l));
		  // 		  cout << "Block1="<<endl<<Block1;
		  // 		  cout << "Block2="<<endl<<Block2;
		  Rst1 = Block1*Block2;
		  //		  cout << "1*2="<<endl<<Rst1;
		  Rst2 = Block2*Block1;
		  //		  cout << "2*1="<<endl<<Rst2;
		  Rst1 -= Rst2;
		  //		  cout << "1*2-2*1="<<endl<<Rst1;
		  bool NonZero = false;
		  for (int k=0; k<NbrSites; ++k)
		    for (int l=0; l<NbrSites; ++l)
		      if (Norm(Rst1[k][l])>1e-12)
			NonZero = true;
		  if (NonZero == true)
		    cout << "Blocks ("<<i<<"," <<i<<"), ("<<i<<", "<<j<<") do NOT commute"<<endl;
		  else
		    cout << "Blocks ("<<i<<"," <<i<<"), ("<<i<<", "<<j<<") commute"<<endl;
		}
	    }
	}
      
      //cout << "Matrix="<<endl<<Rho<<endl;
      // calculate eigenvalues & vectors of Rho
      Rho.Diagonalize(M, Q, 1e-12, 1000);
      for (int i=0; i<DensityMatrixDimension; ++i)
	if (fabs(M[DensityMatrixDimension-1-i])>dynamics*M[DensityMatrixDimension-1])
	  cout << "EV["<<i<<"] = " << M[DensityMatrixDimension-1-i] << endl;
      //cout << "Transition Matrix: "<<endl<<Q<<endl;
      if (Manager.GetInteger("opt-target")>1)
	{
	  double Sum=M[DensityMatrixDimension-1];
	  for (int i=1; i<Manager.GetInteger("opt-target"); ++i)
	    Sum-=M[DensityMatrixDimension-1-i];
	  cout << "EV[0-"<<Manager.GetInteger("opt-target")-1<<"] = " <<Sum<<endl;
	}
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
	      if (Manager.GetBoolean("normalize-phase"))
		{
		  double MaxNorm=0.0;
		  int MaxIndex=0;
		  for (int i = 0; i < EigenState.GetVectorDimension();++i)
		    {
		      if (Norm(EigenState[i])>MaxNorm)
			{
			  MaxNorm=Norm(EigenState[i]);
			  MaxIndex=i;
			}
		    }
		  EigenState *= Polar(1.0,-Arg(EigenState[MaxIndex]));
		}
	      EigenState.WriteVector(RhoVecOut);
	      if (Manager.GetBoolean("plot-density"))
		{
		  sprintf(RhoVecOut,"%s.dm%d.dat",VectorFiles[0],s);
		  ofstream DataFile(RhoVecOut);
		  DataFile << "# X\tY\tv_x\tv_y"<<endl;
		  int CellPosition[2];
		  RealVector SitePosition;
		  for (int x=0; x<Lx; ++x)
		    for (int y=0; y<Ly; ++y)
		      for (int Sub=0; Sub<NbrSubLattices; ++Sub)
			{
			  int q = Space->EncodeQuantumNumber(x, y, Sub, Tmp);
			  double X,Y;
			  if (GenericLattice)
			    {
			      CellPosition[0]=x;
			      CellPosition[1]=y;
			      SitePosition = Lattice->GetSitePosition(CellPosition,Sub);
			      X=SitePosition[0];
			      Y=SitePosition[1];
			    }
			  else
			    {
			      X=(double)x;
			      Y=(double)y;
			    }
			  DataFile << X<<"\t"<<Y<<"\t"<<EigenState[DensityMatrixDimension-1-q].Re<<"\t"<<EigenState[DensityMatrixDimension-1-q].Im<<endl;
			}
		  DataFile.close();
		}
	      
	    }
	  delete [] RhoVecOut;
	}  
    }
      
  if (NbrVectors>=2)
    {
      if ((NbrVectors==2)&&(Manager.GetBoolean("equal")))
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
	      Superposition /= Superposition.Norm();
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
      
	  char *RhoVecOut = new char[strlen(VectorFiles[0])+50];
      
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
      else  // study general superpositions of n vectors
	{
	  MaximallyCondensedStateOnLattice * BestCondensate
	    = new MaximallyCondensedStateOnLattice(Architecture.GetArchitecture(), NbrVectors, Vectors, Space, Lx, Ly, NbrSubLattices);
	  if (Manager.GetBoolean("opt-random"))
	    BestCondensate->RandomizeVariationalParameters();
	  cout << "Optimization in subspace"<<endl;
	  for (int n=0; n<NbrVectors; ++n)
	    cout << "|"<<n<<" > = "<<VectorFiles[n]<<endl;
	  BestCondensate->Optimize(Manager.GetInteger("opt-target"), Manager.GetDouble("opt-tolerance"), Manager.GetInteger("max-iter"));
	  ComplexVector Parameters = BestCondensate->GetVariationalParameters();
	  cout << "Found condensate: "<< Norm(Parameters[0]) << " |0>";
	  for (int i=1; i<NbrVectors; ++i)
	    cout << " + " << Norm(Parameters[i]) << " exp(I*"<<Arg(Parameters[i]) <<") |"<<i<<" >";
	  cout << endl;
	  if (Manager.GetInteger("opt-target")==1)
	    cout << "with EV[0]_max = "<<BestCondensate->GetDensityMatrixEigenvalue()<<endl;
	  else
	    cout << "with EV[0-"<<Manager.GetInteger("opt-target")-1<<"]_max = "<<BestCondensate->GetDensityMatrixEigenvalue()<<endl;
	  if (Manager.GetBoolean("save-vectors"))
	    {
	      char *OldBase = RemoveExtensionFromFileName(VectorFiles[0],".vec");
	      char *OutputName = new char[strlen(OldBase)+20];
	      sprintf(OutputName,"%s.opt_v%d.vec",OldBase,NbrVectors);
	      char *ActualName = GetUniqueFileName(OutputName);
	      ComplexVector CondensedState = BestCondensate->GetWaveFunction();
	      cout << "Writing optimized condensate state: "<<ActualName<<endl;
	      CondensedState.WriteVector(ActualName);
	      delete [] OldBase;
	      delete [] OutputName;
	      delete [] ActualName;
	      delete BestCondensate;
	    }
	}
    }

  if (NbrVectors>1)
    {
      int DensityMatrixDimension2 = NbrSites;
      RealDiagonalMatrix M2;	  
      HermitianMatrix Rho2(DensityMatrixDimension2);
      cout << "====== Analysing sum of density matrices ======" << endl;
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

  if (Manager.GetBoolean("expansion"))
    {
      ComplexVector TmpState(VectorDimension);
      ParticleOnLatticeMomentumOperator *MomentumOperator = new ParticleOnLatticeMomentumOperator(Space, Lx, Ly, NbrSubLattices);
      double *GaugeTransform = NULL;
      int Q;
      // apply gauge transformations
      if (Manager.GetInteger("apply-gauge")!=0)
	{
	  GaugeTransform = new double[NbrSites];
	  //	  TranslationsX = new double[Ly*NbrSubLattices];
	  //	  TranslationsY = new double[Lx*NbrSubLattices];
	  
	  switch (Manager.GetInteger("apply-gauge"))
	    {    
	    case 1:
	      if (!HaveContFlux)
		{
		  cout << "Gauge #1 is written for a continuous staggered flux square lattice!"<<endl;
		  exit(1);
		}
	      for (int x=0; x<Lx; ++x)
		for (int y=0; y<Ly; ++y)
		  {
		    Q = Space->EncodeQuantumNumber(x, y, 0, Tmp);
		    GaugeTransform[Q] = M_PI*ContFlux*x;
		    Q = Space->EncodeQuantumNumber(x, y, 1, Tmp);
		    GaugeTransform[Q] = -M_PI*ContFlux*x;
		  }
// 	      for (int i=0; i<Ly; ++i)
// 		{
// 		  TranslationsX[2*i]=Lx*M_PI*ContFlux;
// 		  TranslationsX[2*i+1]=-Lx*M_PI*ContFlux;
// 		}
// 	      for (int i=0; i<Lx*NbrSubLattices; ++i) TranslationsY[i]=0.0;
	      break;
	    default:
	      cout << "No gauge transform is associated with this code"<<endl;
	      for (int i=0; i<NbrSites; ++i) GaugeTransform[i]=0.0;
// 	      for (int i=0; i<Ly*NbrSubLattices; ++i) TranslationsX[i]=0.0;
// 	      for (int i=0; i<Lx*NbrSubLattices; ++i) TranslationsY[i]=0.0;
	      break;
	    }
	  for (int n=0; n<NbrVectors; ++n)
	    Space->GaugeTransformVector(GaugeTransform, Vectors[n]);
	}

      double OffsetX=0.0, OffsetY=0.0;
      if (Manager.GetDoubles("exp-offset")!=NULL)
	{
	  int Length;
	  double *TmpD = Manager.GetDoubles("exp-offset", Length);
	  if (Length!=2)
	    {
	      cout << "--exp-offset needs to momentum offsets OffsetX,OffsetY"<<endl;
	      exit(1);
	    }
	  OffsetX=TmpD[0];
	  OffsetY=TmpD[1];
	}
      char *DataFileName = new char[24+strlen(VectorFiles[0])];
      for (int n=0; n<NbrVectors; ++n)
	{
	  sprintf(DataFileName,"%s.exp",VectorFiles[n]);
	  ofstream DataFile(DataFileName);
	  DataFile << "# kx\tky\trho"<<endl;
	  for (int kx=0; kx<Lx; ++kx)
	    for (int ky=0; ky<Ly; ++ky)
	      {
		MomentumOperator->SetMomentum(kx,ky,OffsetX,OffsetY);
		VectorOperatorMultiplyOperation Operation (MomentumOperator, &(Vectors[n]), &TmpState);
		Operation.ApplyOperation(Architecture.GetArchitecture());
		if (kx<=Lx/2)
		  DataFile << (2.0*kx*M_PI)/Lx;
		else
		  DataFile << (2.0*(kx-Lx)*M_PI)/Lx;
		if (ky<=Ly/2)
		  DataFile << "\t" << (2.0*ky*M_PI)/Ly;
		else
		  DataFile << "\t" << (2.0*(ky-Ly)*M_PI)/Ly;
		DataFile << "\t" << (Vectors[n]*TmpState).Re << endl;
		// DataFile << "\t" << (Vectors[n]*TmpState).Im << endl;
	      }
	  DataFile.close();
	}
      delete [] DataFileName;
    }
  
  // stop here if we have a generic lattice -> translations as of yet not implemented
  if (Manager.GetBoolean("no-momenta")) exit(0);

  cout<< "====== Analysis of momentum eigenvalues ====="<<endl;

  ParticleOnLatticeTranslationOperator *TranslationOperator= new ParticleOnLatticeTranslationOperator(Space);
    
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
      
      if (GenericLattice)
	cout << "Representation of T_(1)"<<endl<<TrRep<<endl;
      else
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

  if (!GenericLattice)
    {
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
    }

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
	  if (Manager.GetBoolean("normalize-phase"))
	    {
//	      double MaxNorm=0.0;
	      int MaxIndex=0;
// 	      for (int i = 0; i < TmpState.GetVectorDimension();++i)
// 		{
// 		  if (Norm(TmpState[i])>MaxNorm)
// 		    {
// 		      MaxNorm=Norm(TmpState[i]);
// 		      MaxIndex=i;
// 		    }
// 		}
	      TmpState *= Polar(1.0,-Arg(TmpState[MaxIndex]));
	    }

	  cout << "Vector-"<<i<<"="<<vectorName<<endl;
	  TmpState.WriteVector(vectorName);
	}
    }
      
  delete Space;
  delete [] Vectors;
  delete DensityOperator;
  delete TranslationOperator;
  if (Lattice!=NULL) delete Lattice;
}
