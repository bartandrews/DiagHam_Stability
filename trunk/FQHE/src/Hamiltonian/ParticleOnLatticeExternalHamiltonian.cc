////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//   to particles with contact interactions on a lattice in magnetic field    //
//                                                                            //
//                      last modification : 13/02/2008                        //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeExternalHamiltonian.h"
#include "Output/MathematicaOutput.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
//#define DEBUG_OUTPUT



// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrStates = number of quantum states
// oneParticleTerms = file describing single particle terms
// twoParticleTerms = file describing two-particle terms
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
// hermitianFlag = flag indicating whether to use hermitian symmetry
ParticleOnLatticeExternalHamiltonian::ParticleOnLatticeExternalHamiltonian(ParticleOnLattice* particles, int nbrParticles, int nbrStates, const char* oneParticleTerms, const char* twoParticleTerms, AbstractArchitecture* architecture, unsigned long memory, char* precalculationFileName, bool hermitianFlag, LatticePhases *latticeGeometry)
{
  this->Particles=particles;
  this->NbrParticles=nbrParticles;
  if (oneParticleTerms!=NULL)
    {
      this->OneParticleTerms = new char[strlen(oneParticleTerms)+2];
      strcpy(this->OneParticleTerms,oneParticleTerms);
    }
  else this->OneParticleTerms=NULL;
  if (twoParticleTerms!=NULL)
    {
      this->TwoParticleTerms = new char[strlen(twoParticleTerms)+2];
      strcpy(this->TwoParticleTerms,twoParticleTerms);
    }
  else this->TwoParticleTerms=NULL;
  this->LatticeGeometry=latticeGeometry;
  this->HaveKySymmetry=false;
  this->KyMax=0;  
  this->NbrSites = nbrStates;
  this->NbrFluxQuanta=0;
  this->HamiltonianShift=0.0;
  this->FluxDensity=0.0;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  if (hermitianFlag)
    this->HermitianSymmetrizeInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  PrintMemorySize(cout, TmpMemory)<< endl;
	  if (memory > 0)
	    this->EnableFastMultiplication();
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//
ParticleOnLatticeExternalHamiltonian::~ParticleOnLatticeExternalHamiltonian()
{
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->Q3PerQ12;
      delete [] this->Q4PerQ12;
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }
  if (OneParticleTerms!=NULL)
    delete [] this->OneParticleTerms;
  if (this->TwoParticleTerms!=NULL)
    delete [] this->TwoParticleTerms;

}


// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream
ostream& operator << (ostream& Str, ParticleOnLatticeExternalHamiltonian& H)
{
  Str << "Need to implement ostream& operator << for ParticleOnLatticeExternalHamiltonian!" << endl;
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream
MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnLatticeExternalHamiltonian& H)
{
  Str << "Need to implement MathematicaOutput& operator << for ParticleOnLatticeExternalHamiltonian!\n";
  return Str;
}


// evaluate all interaction factors
//   
void ParticleOnLatticeExternalHamiltonian::EvaluateInteractionFactors()
{  
  if (this->OneParticleTerms!=0)
    {
      if (!IsFile(this->OneParticleTerms))
	{
	  cout << "Could not read file with single particle interactions "<<this->OneParticleTerms<<endl;
	  exit(1);
	}
      MultiColumnASCIIFile Parser;
      if ((Parser.Parse(this->OneParticleTerms))&&((Parser.GetNbrColumns()==4)||(Parser.GetNbrColumns()==3)))
	{
	  this->NbrHoppingTerms = Parser.GetNbrLines();
	  this->HoppingTerms = new Complex[NbrHoppingTerms];
	  this->KineticQf = Parser.GetAsIntegerArray (0);
	  this->KineticQi = Parser.GetAsIntegerArray (1);
	  double *TmpRe = Parser.GetAsDoubleArray (2);
	  if (Parser.GetNbrColumns()==4)
	    {
	      double *TmpIm = Parser.GetAsDoubleArray (3);
	      for (int i=0; i<NbrHoppingTerms; ++i)
		this->HoppingTerms[i]=Complex(TmpRe[i],TmpIm[i]);
	      delete [] TmpIm; 
	    }
	  else
	    {
	      for (int i=0; i<NbrHoppingTerms; ++i)
		this->HoppingTerms[i]=Complex(TmpRe[i],0.0);
	    }
	  delete [] TmpRe;
	}
      else
	{
	  cout << "Error parsing single particle interactions "<<this->OneParticleTerms<<
	    " (3 or 4 columns required: Qf, Qi, Re(M), [Im(M)])"<<endl;
	  exit(1);
	}
    }
  else
    {
      this->NbrHoppingTerms = 0;
    }


  if ((this->LatticeGeometry!=NULL)&&(this->LatticeGeometry->HaveOneParticlePotentials()))
    {
      int OldNumberTerms=this->NbrHoppingTerms;
      this->NbrHoppingTerms+=LatticeGeometry->GetNbrLocalPotentials();
      int *NewQi=new int[this->NbrHoppingTerms];
      int *NewQf=new int[this->NbrHoppingTerms];
      Complex *NewHoppingTerms=new Complex[this->NbrHoppingTerms];
      if (OldNumberTerms>0)
	{
	  for (int i=0; i<OldNumberTerms; ++i)
	    {
	      NewQi[i] = KineticQi[i];
	      NewQf[i] = KineticQi[i];
	      NewHoppingTerms[i] = HoppingTerms[i];
	    }
	  delete[]KineticQi;
	  delete[]KineticQf;
	  delete[]HoppingTerms;
	}
      this->KineticQi=NewQi;
      this->KineticQf=NewQf;
      this->HoppingTerms=NewHoppingTerms;
      cout << "Adding one particle potentials in Hamiltonian"<<endl;
      int NbrPotentials;
      int *PotentialPositions;
      double *Potentials = LatticeGeometry->GetOneParticlePotentials(NbrPotentials, PotentialPositions);
      for (int n=0; n<NbrPotentials; ++n)
	{
	  KineticQi[OldNumberTerms] = PotentialPositions[n];
	  KineticQf[OldNumberTerms] = KineticQi[OldNumberTerms];
	  HoppingTerms[OldNumberTerms] = Potentials[n];
#ifdef DEBUG_OUTPUT
	  cout << "H["<<KineticQi[OldNumberTerms]<<"->"<<KineticQf[OldNumberTerms]<<"]="<<HoppingTerms[OldNumberTerms]<<endl;
#endif
	  ++OldNumberTerms;
	}
    }

  this->NbrInteractionFactors=0;
  
  if (this->TwoParticleTerms!=NULL)
    {
      if (!IsFile(this->TwoParticleTerms))
	{
	  cout << "Could not read file with two-particle interactions "<<this->TwoParticleTerms<<endl;
	  exit(1);
	}
      MultiColumnASCIIFile Parser;
      if ((Parser.Parse(this->TwoParticleTerms))&&(Parser.GetNbrColumns()>=6))
	{
	  int TmpNbrLines = Parser.GetNbrLines();
	  this->NbrQ12Indices = 0;
	  this->NbrDiagonalInteractionFactors = 0;
	  int TmpNbrInteractionFactors=0;
	  int *TmpQ1 = Parser.GetAsIntegerArray (0);
	  int *TmpQ2 = Parser.GetAsIntegerArray (1);
	  int *TmpQ3 = Parser.GetAsIntegerArray (2);
	  int *TmpQ4 = Parser.GetAsIntegerArray (3);
	  double *TmpRe = Parser.GetAsDoubleArray (4);
	  double *TmpIm = Parser.GetAsDoubleArray (5);
	  // count pairs of Q1,Q2 first
	  int oldQ1=TmpQ1[0], oldQ2=TmpQ2[0];
	  int *TmpNbrQ34Values = new int[TmpNbrLines];
	  bool HaveLargerQ1=false, HaveSmallerQ1=false;
	  bool HaveLargerQ3=false, HaveSmallerQ3=false;
	  TmpNbrQ34Values[0]=0;
	  int Pos=0;
	  bool HaveZeroIndex=false;
	  while (Pos<TmpNbrLines)
	    {
	      if ((!HaveZeroIndex)&&((TmpQ1[Pos]==0)||(TmpQ2[Pos]==0)||(TmpQ3[Pos]==0)||(TmpQ4[Pos]==0)))
		HaveZeroIndex=true;
	      while ((Pos<TmpNbrLines)&&(TmpQ1[Pos]==oldQ1)&&(TmpQ2[Pos]==oldQ2))
		{
		  if (TmpQ1[Pos]>TmpQ2[Pos]) HaveLargerQ1=true;
		  else if (TmpQ1[Pos]<TmpQ2[Pos]) HaveSmallerQ1=true;
		  if (TmpQ3[Pos]>TmpQ4[Pos]) HaveLargerQ3=true;
		  else if (TmpQ3[Pos]<TmpQ4[Pos]) HaveSmallerQ3=true;
		  // have diagonal element?
		  if ((TmpQ1[Pos]==TmpQ2[Pos])&&(TmpQ1[Pos]==TmpQ3[Pos])&&(TmpQ1[Pos]==TmpQ4[Pos]))
		    {
		      ++this->NbrDiagonalInteractionFactors;
		      //cout << Pos<<" : "<< TmpQ1[Pos]<<" " << TmpQ2[Pos]<<" " << TmpQ3[Pos]<<" " << TmpQ4[Pos]<<" : diagonal"<<endl;
		    }
		  else
		    {
		      ++TmpNbrQ34Values[NbrQ12Indices];
		      ++TmpNbrInteractionFactors;
		      //cout << Pos<<" : "<< TmpQ1[Pos]<<" " << TmpQ2[Pos]<<" " << TmpQ3[Pos]<<" " << TmpQ4[Pos]<<" : off-diagonal"<<endl;
		    }
		  ++Pos;
		}
	      if (Pos<TmpNbrLines)
		{
		  if (TmpQ1[Pos]>TmpQ2[Pos]) HaveLargerQ1=true;
		  else if (TmpQ1[Pos]<TmpQ2[Pos]) HaveSmallerQ1=true;
		  if (TmpQ3[Pos]>TmpQ4[Pos]) HaveLargerQ3=true;
		  else if (TmpQ3[Pos]<TmpQ4[Pos]) HaveSmallerQ3=true;
		  ++NbrQ12Indices;	     
		  oldQ1=TmpQ1[Pos];
		  oldQ2=TmpQ2[Pos];
		  // have diagonal element?
		  if ((TmpQ1[Pos]==TmpQ2[Pos])&&(TmpQ1[Pos]==TmpQ3[Pos])&&(TmpQ1[Pos]==TmpQ4[Pos]))
		    {
		      ++this->NbrDiagonalInteractionFactors;
		      TmpNbrQ34Values[NbrQ12Indices]=0;
		      //cout << Pos<<" : "<< TmpQ1[Pos]<<" " << TmpQ2[Pos]<<" " << TmpQ3[Pos]<<" " << TmpQ4[Pos]<<" : diagonal"<<endl;
		    }
		  else
		    {
		      ++TmpNbrInteractionFactors;
		      TmpNbrQ34Values[NbrQ12Indices]=1;
		      //cout << Pos<<" : "<< TmpQ1[Pos]<<" " << TmpQ2[Pos]<<" " << TmpQ3[Pos]<<" " << TmpQ4[Pos]<<" : off-diagonal"<<endl;
		    }
		  ++Pos;
		}
	    }
	  if (HaveZeroIndex==false)
	    {
	      cout << "Warning: quantum number q=0 not found - need to shift indices?"<<endl;
	    }
	  ++NbrQ12Indices;
#ifdef DEBUG_OUTPUT
	  int sum=NbrDiagonalInteractionFactors;
	  cout << "Count of matrix elements:"<<endl
	       << "Diagonal: "<<NbrDiagonalInteractionFactors<<endl
	       << "NbrQ12Indices= "<<NbrQ12Indices<<endl;
	  for (int i=0; i<NbrQ12Indices; ++i)
	    {
	      cout << "  NbrQ34Values["<<i<<"]="<<TmpNbrQ34Values[i]<<endl;
	      sum+=TmpNbrQ34Values[i];
	    }
	  cout << "total elements: "<<sum<<" (lines: "<<TmpNbrLines<<")"<<endl;
#endif
	  // assign memory
	  this->NbrQ34Values = new int[NbrQ12Indices];
	  for (int i=0; i<NbrQ12Indices; ++i)
	    this->NbrQ34Values[i]=TmpNbrQ34Values[i];
	  delete [] TmpNbrQ34Values;
	  this->InteractionFactors = new Complex[TmpNbrInteractionFactors];
	  this->Q1Value = new int[NbrQ12Indices];
	  this->Q2Value = new int[NbrQ12Indices];
	  this->Q3PerQ12 = new int*[NbrQ12Indices];
	  this->Q4PerQ12 = new int*[NbrQ12Indices];
	  for (int i=0; i<NbrQ12Indices; ++i)
	    {
	      this->Q3PerQ12[i] = new int[NbrQ34Values[i]];
	      this->Q4PerQ12[i] = new int[NbrQ34Values[i]];
	    }
	  this->DiagonalQValues = new int[NbrDiagonalInteractionFactors];
	  this->DiagonalInteractionFactors = new double[NbrDiagonalInteractionFactors];
	  // test whether we need to apply symmetry factors
	  bool Q12Symmetry=false, Q34Symmetry=false;
	  if (HaveSmallerQ1^HaveLargerQ1)
	    Q12Symmetry=true;
	  if (HaveSmallerQ3^HaveLargerQ3)
	    Q34Symmetry=true;
	  if ((Q12Symmetry)||(Q34Symmetry))
	    cout << "Assuming symmetry in";
	  if (Q12Symmetry) cout << " Q12";
	  if (Q34Symmetry) cout << " Q34";
	  cout<<endl;
	  // read matrix elements into storage structure
	  int Q12Index=0;
	  int Q34Index=0;
	  int GeneralIndex=0;
	  int DiagonalIndex=0;
	  Q1Value[0]=TmpQ1[0];
	  Q2Value[0]=TmpQ2[0];
	  oldQ1=TmpQ1[0];
	  oldQ2=TmpQ2[0];
	  Pos=0;
	  while (Pos<TmpNbrLines)
	    {
	      while ((Pos<TmpNbrLines)&&(TmpQ1[Pos]==oldQ1)&&(TmpQ2[Pos]==oldQ2))
		{
		  // have diagonal element?
		  if ((TmpQ1[Pos]==TmpQ2[Pos])&&(TmpQ1[Pos]==TmpQ3[Pos])&&(TmpQ1[Pos]==TmpQ4[Pos]))
		    {
		      // cout << "Diagonal Interaction "<<TmpRe[Pos]<<" v_"<<TmpQ1[Pos]<<endl;
		      this->DiagonalQValues[DiagonalIndex]=TmpQ1[Pos];
		      this->DiagonalInteractionFactors[DiagonalIndex++]=TmpRe[Pos];
		    }
		  else
		    {
		      Q3PerQ12[Q12Index][Q34Index]=TmpQ3[Pos];
		      Q4PerQ12[Q12Index][Q34Index++]=TmpQ4[Pos];
		      double SymFactor=1.0;
		      if ((Q12Symmetry)&&(TmpQ1[Pos]!=TmpQ2[Pos])) SymFactor*=2.0;
		      if ((Q34Symmetry)&&(TmpQ3[Pos]!=TmpQ4[Pos])) SymFactor*=2.0;
		      // cout << "Interaction "<<TmpRe[Pos]<<"+I*"<<TmpIm[Pos]<<" with Sym "<<SymFactor<<"  v_"<<TmpQ1[Pos]<<" "
		      //      << TmpQ2[Pos]<<" "<<TmpQ3[Pos]<<" "<<TmpQ4[Pos]<<" "<<endl;
		      this->InteractionFactors[GeneralIndex++]=SymFactor*Complex(TmpRe[Pos],TmpIm[Pos]);
		    }
		  ++Pos;
		}
	      if (Pos<TmpNbrLines)
		{
		  ++Q12Index;
		  Q34Index=0;
		  Q1Value[Q12Index]=TmpQ1[Pos];
		  Q2Value[Q12Index]=TmpQ2[Pos];
		  oldQ1=TmpQ1[Pos];
		  oldQ2=TmpQ2[Pos];
		  // have diagonal element?
		  if ((TmpQ1[Pos]==TmpQ2[Pos])&&(TmpQ1[Pos]==TmpQ3[Pos])&&(TmpQ1[Pos]==TmpQ4[Pos]))
		    {
		      // cout << "Diagonal Interaction "<<TmpRe[Pos]<<" v_"<<TmpQ1[Pos]<<endl;
		      this->DiagonalQValues[DiagonalIndex]=TmpQ1[Pos];
		      this->DiagonalInteractionFactors[DiagonalIndex++]=TmpRe[Pos];
		    }
		  else
		    {
		      Q3PerQ12[Q12Index][Q34Index]=TmpQ3[Pos];
		      Q4PerQ12[Q12Index][Q34Index++]=TmpQ4[Pos];
		      double SymFactor=1.0;
		      if ((Q12Symmetry)&&(TmpQ1[Pos]!=TmpQ2[Pos])) SymFactor*=2.0;
		      if ((Q34Symmetry)&&(TmpQ3[Pos]!=TmpQ4[Pos])) SymFactor*=2.0;
		      // cout << "Interaction "<<TmpRe[Pos]<<"+I*"<<TmpIm[Pos]<<" with Sym "<<SymFactor<<"  v_"<<TmpQ1[Pos]<<" "
		      //      <<TmpQ2[Pos]<<" "<<TmpQ3[Pos]<<" "<<TmpQ4[Pos]<<" "<<endl;		      
		      this->InteractionFactors[GeneralIndex++]=SymFactor*Complex(TmpRe[Pos],TmpIm[Pos]);
		    }
		  ++Pos;
		}
	    }
	  if (GeneralIndex!=TmpNbrInteractionFactors)
	    {
	      cout << "Inconsistency in count of matrix elements for ParticleOnLatticeExternalHamiltonian"<<endl;
	      exit(1);
	    }
	  if (DiagonalIndex!=NbrDiagonalInteractionFactors)
	    {
	      cout << "Inconsistency in count of diagonal matrix elements for ParticleOnLatticeExternalHamiltonian"<<endl;
	      exit(1);
	    }
	  delete [] TmpQ1;
	  delete [] TmpQ2;
	  delete [] TmpQ3;
	  delete [] TmpQ4;
	  delete [] TmpRe;
	  delete [] TmpIm;
	}
      else
	{
	  cout << "Error parsing two-particle interactions "<<this->TwoParticleTerms<<
	    " [6 columns required: Q1, Q2, Q3, Q4, Re(M), Im(M)]"<<endl;
	  exit(-1);
	}
    }
  else
    {
      // we have no general four-particle interactions:     
      this->NbrQ12Indices=0;
      cout << "No two-particle interactions"<<endl;
    }
  
}
