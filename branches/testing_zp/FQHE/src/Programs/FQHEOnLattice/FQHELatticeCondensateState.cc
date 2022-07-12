#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/QHEOnLatticeFileTools.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

#include <bitset>
using std::bitset;


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



class GutzwillerWaveFunction
{
public:
  // constructor
  // nbrParticles = particles in the condensate (should match space)
  // condensateWF = condensate wavefunction, in notation as single-particle hilbert-space
  // unoccupiedTerms = norms of empty states, allowing form \prod_i (Norm_i + \psi_i a^\dagger_i) |state>
  // space = target-space of many-body state
  GutzwillerWaveFunction(int nbrParticles, ComplexVector& condensateWF, ComplexVector& unoccupiedTerms, ParticleOnLattice *space);

  // destructor
  ~GutzwillerWaveFunction();

  // get Many-Body state
  // return = resultingState
  ComplexVector & GetGutzwillerWaveFunction();    

protected:

  // main recursion to calculate State (\sum \psi_i a^\dagger_i)^N |state>
  // exponent = power N remaining to be applied
  // state = state to be acted upon
  // prefactor = previous coefficients applied to state
  // in last stage of recursion, writes to this->TargetVector
  void Product (int exponent, unsigned long state, Complex prefactor);

  // main recursion to calculate State \prod_i(\sum \chi_i + \psi_i a^\dagger_i) |state>
  // exponent = power N remaining to be applied
  // state = state to be acted upon
  // prefactor = previous coefficients applied to state
  // in last stage of recursion, writes to this->TargetVector
  void Product2 (int nextQ, int nbrBosons, unsigned long state, Complex prefactor);

  // target state, internal use
  ComplexVector TargetVector;

  // target condensate wavefunction
  ComplexVector CondensateState;

  // additional terms for wavefunction
  ComplexVector UnoccupiedTerms;

  // number of particles in many-body state
  int NbrParticles;

  // target Hilbert space
  ParticleOnLattice *Space;

  // number of states on lattice
  int NbrStates;

  // dimension of Space
  int Dim;
};


// constructor
// nbrParticles = particles in the condensate (should match space)
// condensateWF = condensate wavefunction, in notation as single-particle hilbert-space
// space = target-space of many-body state
GutzwillerWaveFunction::GutzwillerWaveFunction(int nbrParticles, ComplexVector &condensateWF, ComplexVector& unoccupiedTerms, ParticleOnLattice *space)
{
  this->NbrParticles = nbrParticles;
  this->Space = space;
  this->CondensateState = ComplexVector(condensateWF,true);
  this->UnoccupiedTerms = ComplexVector(unoccupiedTerms,true);
  // believe this is the way it should be done to get CondensateState[q] associated with quantum no. q.
//   cout << "initial condensate: "<<CondensateState<<endl;
  this->CondensateState.ReverseVector();
  this->UnoccupiedTerms.ReverseVector();
//   cout << "reversed condensate: "<<CondensateState<<endl;
  this->Dim = Space->GetHilbertSpaceDimension();
  this->NbrStates = CondensateState.GetVectorDimension();
}

// destructor
GutzwillerWaveFunction::~GutzwillerWaveFunction()
{
}

// get Many-Body state
// return = resultingState
ComplexVector & GutzwillerWaveFunction::GetGutzwillerWaveFunction()
{
  this->TargetVector.Resize(Dim);
  this->TargetVector.ClearVector();
  // call main recursion
  if (this->UnoccupiedTerms.Norm()<1e-6)
    {
      cout << "Using b^+^N|0>"<<endl;
      this->Product(NbrParticles, 0x0ul, 1.0);
    }
  else
    {
      cout << "Using prod_i (a_o + b^+_i) |0>"<<endl;
      this->Product2(NbrStates-1, 0, 0x0ul, 1.0);
    }
  this->TargetVector/=this->TargetVector.Norm();
//   cout <<"Test norm: "<<TargetVector.Norm()<<endl;
  return this->TargetVector;
}

// main recursion to calculate State (\sum \psi_i a^\dagger_i)^N |state>
// exponent = power N remaining to be applied
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerWaveFunction::Product (int exponent, unsigned long state, Complex prefactor)
{
  int Index;
  unsigned long ResultingState;
  double AdFactor;
  if (exponent>1)
    {
      for (int q=0; q<this->NbrStates; ++q)
	{
	  ResultingState = Space->Ad(state, q, AdFactor);
	  if (ResultingState!=0x0ul)
	    Product(exponent-1, ResultingState, prefactor*AdFactor*CondensateState[q]);
	}
    }
  else
    {      
      for (int q=0; q<this->NbrStates; ++q)
	{
	  ResultingState = Space->Ad(state, q, AdFactor);
	  if (ResultingState!=0x0ul)
	    {	      
	      if ((Index=Space->CarefulFindStateIndex(ResultingState,-1))<Dim)
		{
		  TargetVector[Index]+= prefactor*AdFactor*CondensateState[q];
		}
	    }

	}
    }
}


// main recursion to calculate State \prod_i (\sum \chi_i + \psi_i a^\dagger_i) |state>
// nextQ = value quantum number in next operator to be applied
// nbrBosons = number of bosons already in state
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerWaveFunction::Product2 (int nextQ, int nbrBosons, unsigned long state, Complex prefactor)
{
  int Index;
  unsigned long ResultingState;
  double AdFactor;
  if (nextQ>0)
    {
      if (nbrBosons<this->NbrParticles)
	{
	  ResultingState = Space->Ad(state,nextQ,AdFactor);
	  Product2(nextQ-1, nbrBosons+1, ResultingState, prefactor*CondensateState[nextQ]);
	}
      Product2(nextQ-1, nbrBosons, state, prefactor*AdFactor*UnoccupiedTerms[nextQ]);
    }
  else
    {
      if (nbrBosons==this->NbrParticles-1)
	{
	  ResultingState = Space->Ad(state,nextQ, AdFactor);
	  if ((Index=Space->CarefulFindStateIndex(ResultingState,-1))<Dim)
	    {
	      TargetVector[Index]+= prefactor*AdFactor*CondensateState[nextQ];
	    }
	}
      else if (nbrBosons==this->NbrParticles)
	{
	  if ((Index=Space->CarefulFindStateIndex(state,-1))<Dim)
	    {
	      TargetVector[Index]+= prefactor*UnoccupiedTerms[nextQ];
	    }
	}
    }
}


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeCondensateState" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "condensate", "filename of vector describing condensate WF");

  (*SystemGroup) += new SingleStringOption  ('a', "alt-input", "alternative input for filename of vector describing condensate WF in format provided by Nigel");
  

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles of target many-body state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (grabbed from condensate)", 0);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons (~Gutzwiller projection)");
  (*SystemGroup) += new BooleanOption('n',"no-hard-core","Do not use Hilbert-space of hard-core bosons (overriding detection from filename)");
  (*SystemGroup) += new BooleanOption('\n',"alt-ignore-amplitude","Ignore amplitude of condensate functions theta_i");
  (*SystemGroup) += new SingleIntegerOption('\n',"alt-format","0: Nigel's format (theta, phi), 1: my format (psi_x, psi_y)",0);

  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
    
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  unsigned long MemorySpace = 9ul << 20;

  char* CondensateFile;  

  if (Manager.GetString("alt-input")!=NULL)
    {
      CondensateFile = Manager.GetString("alt-input");
    }
  else
    {
      CondensateFile = Manager.GetString("condensate");
      if (Manager.GetString("condensate")==0)
	{
	  cout << "A vector for the condensate is required" << endl;
	  exit(1);
	}
    }  

  double Interaction=10.0;
  int TmpI=-1;
  bool Statistics=false;
  bool HardCore=false;
  int NbrBosonsCondensate=0;
  
  
  if (FQHEOnLatticeFindSystemInfoFromVectorFileName(CondensateFile, NbrBosonsCondensate, Lx, Ly, Interaction, NbrFluxQuanta, TmpI, Statistics, HardCore) == false)
    {
      cout<<"Please use standard file-names, or indicate all system parameters!"<<endl;
      exit(1);
    }

  HardCore=(HardCore||Manager.GetBoolean("hard-core"));

  if (Manager.GetBoolean("no-hard-core"))
    HardCore=false;

  if (NbrBosons==0) NbrBosons=NbrBosonsCondensate;

  cout << "Calculating Gutzwiller state for "<<NbrBosons<<" bosons on a lattice with geometry"<<endl;
  cout << "Lx="<<Lx<<", Ly="<<Ly<<", N_phi="<<NbrFluxQuanta;
  if (HardCore)
    cout <<" and hardcore interactions."<<endl;
  else
    cout << endl;

  ComplexVector CondensateState;
  ComplexVector UnoccupiedNorms;

  if (Manager.GetString("alt-input")!=NULL)
    {
      MultiColumnASCIIFile Parser(' ');

      if (Parser.Parse(CondensateFile)==false)
	{
	  cout<<"Could not read condensate file "<<CondensateFile<<" !"<<endl;
	  exit(1);
	}

      if (Parser.GetNbrColumns()!=4)
	{
	  cout<<"Error: Four columns, [ x y theta phi ] expected in condensate file !"<<endl;
	  exit(1);
	}
      if (Parser.GetNbrLines()!=Lx*Ly)
	{
	  cout<<"Error: Wrong number of lines in condensate file !"<<endl;
	  exit(1);
	}

      CondensateState.Resize(Lx*Ly);
      UnoccupiedNorms.Resize(Lx*Ly);

      int *x = Parser.GetAsIntegerArray(0);
      int *y = Parser.GetAsIntegerArray(1);

      double *theta = Parser.GetAsDoubleArray(2);
      double *phi = Parser.GetAsDoubleArray(3);

      ParticleOnLattice *Space = new BosonOnLattice(1, Lx, Ly, NbrFluxQuanta, MemorySpace);

      int Dim=Lx*Ly;
      Complex Tmp;
      if (Manager.GetInteger("alt-format")==0)
	{
	  if (Manager.GetBoolean("alt-ignore-amplitude"))
	    {
	      for (int i=0; i<Dim; ++i)
		{      
		  int q = Space->EncodeQuantumNumber(x[i]-1, y[i]-1, 0, Tmp);
		  CondensateState[Dim-1-q]=Polar(1.0, phi[i]); // CondensateState[Dim-1-q]
		  UnoccupiedNorms[Dim-1-q].Re=0.0;
		  UnoccupiedNorms[Dim-1-q].Im=0.0;
		}
	    }
	  else
	    {
	      cout << "# x \t y \t vx \t vy"<<endl;
	      for (int i=0; i<Dim; ++i)
		{      
		  int q = Space->EncodeQuantumNumber(x[i]-1, y[i]-1, 0, Tmp);
		  CondensateState[Dim-1-q]=Complex(cos(theta[i]/2.0)*cos(phi[i]), cos(theta[i]/2.0)*sin(phi[i])); // CondensateState[Dim-1-q]
		  UnoccupiedNorms[Dim-1-q].Re=sin(theta[i]/2.0);
		  UnoccupiedNorms[Dim-1-q].Im=0.0;
		  cout <<x[i]-1<<"\t"<<y[i]-1<<"\t"<<CondensateState[Dim-1-q].Re<<"\t"<<CondensateState[Dim-1-q].Im<<endl;
		}
	      cout <<"CondensateState:"<<endl<<CondensateState<<endl;
	    }
	}
      else if (Manager.GetInteger("alt-format")==1)
	{
	  for (int i=0; i<Dim; ++i)
	    {      
	      int q = Space->EncodeQuantumNumber(x[i], y[i], 0, Tmp);
	      CondensateState[Dim-1-q].Re= theta[i]; 
	      CondensateState[Dim-1-q].Im= phi[i]; 
	      UnoccupiedNorms[Dim-1-q].Re=0.0;
	      UnoccupiedNorms[Dim-1-q].Im=0.0;
	    }
	}
      delete Space;      
      
    }
  else
    {
      if (CondensateState.ReadVector(Manager.GetString("condensate"))==false)
	{
	  cout<<"Could not read condensate state!"<<endl;
	  exit(1);
	}
      if (CondensateState.GetVectorDimension()!=Lx*Ly)
	{
	  cout<<"Number of sites in condensate does not match the lattice Lx="<<Lx<<", Ly="<<Ly<<"!"<<endl;
	  exit(1);
	}
      UnoccupiedNorms.Resize(Lx*Ly);
      for (int i=0; i<Lx*Ly; ++i) UnoccupiedNorms[i]=0.0;
    }
  
  

  ParticleOnLattice* Space;
  if (HardCore)
    Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  else Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);

  GutzwillerWaveFunction GutzwillerState(NbrBosons, CondensateState, UnoccupiedNorms, Space);

  ComplexVector GutzwillerStateVector = GutzwillerState.GetGutzwillerWaveFunction();
  char extension[20];
  sprintf(extension,"N_%d.gw",NbrBosons);
  char *TmpFileName = AddExtensionToFileName(CondensateFile, extension);

  GutzwillerStateVector.WriteVector(TmpFileName);

  delete [] TmpFileName;
}
