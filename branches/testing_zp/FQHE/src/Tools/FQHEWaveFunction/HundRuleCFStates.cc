////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2008 Gunnar Moeller                   //
//                                                                            //
//                                                                            //
//           class implementing composite fermion state with partially        //
//             filled highest CF shell for a wave function on sphere          //
//                                                                            //
//                        last modification : 16/01/2008                      //
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


#include "HundRuleCFStates.h"
#include "Tools/FQHEMonteCarlo/ParticleOnSphereCollection.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include <iostream>

// switch debugging output
#define DEBUG

using std::cout;
using std::endl;


// default constructor
HundRuleCFStates::HundRuleCFStates()
{
  this->NbrParticles = 0;
}

// standard constructor
// nbrParticles = number of particles in system
// nbrEffectiveFlux = effective flux seen by composite fermions
// jastrowP = power of jastrow factors div 2
// overrideK = calculate for value of angular momentum < maximumL (only active for 2 or more particles)
//
HundRuleCFStates::HundRuleCFStates(int nbrParticles, int nbrEffectiveFlux, int jastrowP, int overrideL)
{
  this->NbrParticles = nbrParticles;
  if (nbrParticles<0)
    {
      cout << "The number of particles has to be a positive integer" << endl;
      exit(1);
    }
  this->TwoSEff = abs(nbrEffectiveFlux);
  if (nbrEffectiveFlux<0) ReverseFluxFlag = true;
  else ReverseFluxFlag = false;
  this->TwoS = 2*jastrowP*(nbrParticles-1) + nbrEffectiveFlux;
  this->JastrowP = jastrowP;
  this->NumShells=0;
  this->ElementNorm=1.0;
  this->Flag.Initialize();
  int particlesLeft=NbrParticles;
  while (particlesLeft>0)
    {
      particlesLeft-=TwoSEff+1+2*NumShells;
      ++NumShells;	
    }
  this->NbrOrbitalsInHighestShell=TwoSEff+1+2*(NumShells-1);
  this->NbrParticlesInHighestShell = particlesLeft+NbrOrbitalsInHighestShell;
  this->HighestShellLzMax = NbrOrbitalsInHighestShell - 1;
  this->TotalL=0;
  for (int i=0; i<NbrParticlesInHighestShell;++i)
    {
      this->TotalL += 2*(HighestShellLzMax-i) - HighestShellLzMax;
    }
  this->SelectMPosition=this->TotalL;
  this->Orbitals = new JainCFOnSphereOrbitals(nbrParticles, NumShells, nbrEffectiveFlux, 2*JastrowP);
  this->Ji = new Complex[this->NbrParticles];
#ifdef __USE_LAPACK_HERE__
  this->SlaterDeterminant.Resize(NbrParticles);
#else
  this->SlaterDeterminant.Resize(NbrParticles, NbrParticles);
#endif
  cout << "Configuration of CF's calculated:"<<endl;
  cout << NbrParticles << " electrons filling " << NumShells << " CF shells"<<endl;
  cout << NbrParticlesInHighestShell << " electrons in highest shell with total L = ";  
  if (TotalL&1) cout << TotalL<<"/2";
  else cout << TotalL/2.0;
  cout << endl;

  if (TotalL==0) // very trivial case...
    {
      cout << "This is a filled shell state -> can use more efficient class to generate" << endl;
      // just reserve some space to avoid different case in destructor
      this->NbrTermsPerLz = new int[1];
      this->TermsPerLz = new SlaterComponent*[1];
      this->NbrTermsPerLz[0]=1;
      this->TermsPerLz[0]=new SlaterComponent[1];
      this->TermsPerLz[0][0]=SlaterComponent(NbrParticlesInHighestShell,HighestShellLzMax);
    }  
  else if (NbrParticlesInHighestShell==1) // trivial case...
    {
      this->NbrTermsPerLz = new int[TotalL+1];
      this->TermsPerLz = new SlaterComponent*[TotalL+1];
      for (int i=0; i<=TotalL; ++i)
	{
	  this->NbrTermsPerLz[i]=1;
	  this->TermsPerLz[i]=new SlaterComponent[1];
	  int *tmpI=new int[1];
	  tmpI[0]=i;
	  this->TermsPerLz[i][0]=SlaterComponent(1.0,1,tmpI,HighestShellLzMax);
	}
      this->LEigenstates = 0;
    }
  else if (NbrParticlesInHighestShell==HighestShellLzMax) // still trivial: one hole...
    {
      this->NbrTermsPerLz = new int[TotalL+1];
      this->TermsPerLz = new SlaterComponent*[TotalL+1];
      for (int i=0; i<=TotalL; ++i)
	{
	  this->NbrTermsPerLz[i]=1;
	  this->TermsPerLz[i]=new SlaterComponent[1];
	  int *tmpI=new int[NbrParticlesInHighestShell];
	  int k=0;
	  for(int l=0; l<NbrParticlesInHighestShell; ++l)
	    {
	      if (k==TotalL-i) ++k;
	      tmpI[l]=k++;
	    }
	  this->TermsPerLz[i][0]=SlaterComponent(1.0,NbrParticlesInHighestShell,tmpI,HighestShellLzMax);
	}
      this->LEigenstates = 0;
    }
  else if (NbrParticlesInHighestShell==2) // case entirely given by Clebsch-Gordan coefficients
    {
      if ((overrideL>=0) && (this->TotalL>overrideL))
	this->TotalL=overrideL;
      ClebschGordanCoefficients VectorCoupling(HighestShellLzMax,HighestShellLzMax);
      this->NbrTermsPerLz = new int[TotalL+1];
      this->TermsPerLz = new SlaterComponent*[TotalL+1];
      int *TmpM1=new int[TotalL+1];
      int *TmpM2=new int[TotalL+1];
      for (int i=0; i<=TotalL; ++i)
	{
	  this->NbrTermsPerLz[i]=0;
	  int M = 2*i - TotalL;
	  for (int m1=-HighestShellLzMax+2; m1<=HighestShellLzMax; m1+=2)
	    {
	      int m2 = M - m1;
	      if ((m2 < m1) && (m2>=-HighestShellLzMax) && (m2<=HighestShellLzMax))
		{
		  if ((VectorCoupling.GetCoefficient (m1, m2, this->TotalL) != 0.0) &&
		      (fabs(VectorCoupling.GetCoefficient (m1, m2, this->TotalL) +
			    VectorCoupling.GetCoefficient (m2, m1, this->TotalL)) < 1e-13))
		    {
		      TmpM1[this->NbrTermsPerLz[i]] = m1;
		      TmpM2[this->NbrTermsPerLz[i]] = m2;		      
		      ++this->NbrTermsPerLz[i];
		    }
		}
	    }
	  this->TermsPerLz[i]=new SlaterComponent[this->NbrTermsPerLz[i]];
	  for (int k=0; k<this->NbrTermsPerLz[i]; ++k)
	    {
	      int *tmpI=new int[2];
	      tmpI[0]=(TmpM2[k]+HighestShellLzMax)/2;
	      tmpI[1]=(TmpM1[k]+HighestShellLzMax)/2;
	      this->TermsPerLz[i][k]=SlaterComponent(VectorCoupling.GetCoefficient (TmpM1[k], TmpM2[k],
						       this->TotalL), 2, tmpI, HighestShellLzMax);
	    }
	}
      delete [] TmpM1;
      delete [] TmpM2;
      this->LEigenstates = 0;
    }
  else if (NbrParticlesInHighestShell==HighestShellLzMax-1) // case still needs to be checked...
    // 2 holes in upper shell: case still entirely given by Clebsch-Gordan coefficients
    {
      cout << "2 holes in upper shell"<<endl;
      ClebschGordanCoefficients VectorCoupling(HighestShellLzMax,HighestShellLzMax);
      this->NbrTermsPerLz = new int[TotalL+1];
      this->TermsPerLz = new SlaterComponent*[TotalL+1];
      int *TmpM1=new int[TotalL+1];
      int *TmpM2=new int[TotalL+1];
      for (int i=0; i<=TotalL; ++i)
	{
	  this->NbrTermsPerLz[i]=0;
	  int M = 2*i - TotalL;
	  for (int m1=-HighestShellLzMax+2; m1<=HighestShellLzMax; m1+=2)
	    {
	      int m2 = M - m1;
	      if ((m2 < m1) && (m2>=-HighestShellLzMax) && (m2<=HighestShellLzMax))
		{
		  if ((VectorCoupling.GetCoefficient (m1, m2, this->TotalL) != 0.0) &&
		      (fabs(VectorCoupling.GetCoefficient (m1, m2, this->TotalL) +
			    VectorCoupling.GetCoefficient (m2, m1, this->TotalL)) < 1e-13))
		    {
		      TmpM1[this->NbrTermsPerLz[i]] = m1;
		      TmpM2[this->NbrTermsPerLz[i]] = m2;		      
		      ++this->NbrTermsPerLz[i];
		    }
		}
	    }
	  this->TermsPerLz[i]=new SlaterComponent[this->NbrTermsPerLz[i]];
	  for (int k=0; k<this->NbrTermsPerLz[i]; ++k)
	    {
	      int *tmpI=new int[NbrParticlesInHighestShell];
	      int avoid1 = (-TmpM1[k]+HighestShellLzMax)/2;
	      int avoid2 = (-TmpM2[k]+HighestShellLzMax)/2;
	      int q=0;
	      for(int l=0; l<NbrParticlesInHighestShell; ++l)
		{
		  if ((q==avoid1)||(q==avoid2))
		    {
		      ++q;
		      if ((q==avoid1)||(q==avoid2))
			++q;
		    }
		  tmpI[l]=q++;
		}
	      this->TermsPerLz[i][k]=SlaterComponent(VectorCoupling.GetCoefficient (TmpM1[k], TmpM2[k],
						       this->TotalL), NbrParticlesInHighestShell,
						     tmpI, HighestShellLzMax);
	    }
	}
      delete [] TmpM1;
      delete [] TmpM2;
      this->LEigenstates = 0;
    }  
  else // NbrParticlesInHighestShell > 2:
    {
      
      this->LEigenstates = new SlaterSuperposition[TotalL+1];
      this->LEigenstates[TotalL] = SlaterSuperposition(NbrParticlesInHighestShell, HighestShellLzMax);
      for (int i=0; i<TotalL; ++i)
	this->LEigenstates[TotalL-i-1]=this->LEigenstates[TotalL-i].ApplyLMinus();
      
      for (int i=0; i<=TotalL; ++i)
	cout << this->LEigenstates[i] << endl;
    }
#ifdef DEBUG
  for (int i=0; i<=TotalL; ++i)
    {
      cout << "State[M="<<(2*i-TotalL)/2.0<<"]= ";
       cout << TermsPerLz[i][0];
      for (int k=1; k<this->NbrTermsPerLz[i]; ++k)
	cout << " + "<<TermsPerLz[i][k];
      cout << endl;
    }
#endif  
}


// copy constructor
HundRuleCFStates::HundRuleCFStates(HundRuleCFStates &toCopy)
{
  this->NbrParticles=toCopy.NbrParticles;
  this->NumShells=toCopy.NumShells;
  this->NbrOrbitalsInHighestShell=toCopy.NbrOrbitalsInHighestShell;
  this->NbrParticlesInHighestShell=toCopy.NbrParticlesInHighestShell;
  this->HighestShellLzMax=toCopy.HighestShellLzMax;
  this->TotalL=toCopy.TotalL;
  this->ReverseFluxFlag=toCopy.ReverseFluxFlag;
  this->TwoS=toCopy.TwoS;
  this->TwoSEff=toCopy.TwoSEff;
  this->JastrowP=toCopy.JastrowP;
  this->Orbitals=toCopy.Orbitals;
  this->Ji = toCopy.Ji;
  this->LEigenstates=toCopy.LEigenstates;
  this->NbrTermsPerLz=toCopy.NbrTermsPerLz;
  this->TermsPerLz=toCopy.TermsPerLz;
  this->ElementNorm=toCopy.ElementNorm;
  this->SelectMPosition=toCopy.SelectMPosition;
#ifdef __USE_LAPACK_HERE__
  this->SlaterDeterminant.Resize(NbrParticles);
#else
  this->SlaterDeterminant.Resize(NbrParticles, NbrParticles);
#endif
  this->Flag=toCopy.Flag;
}

HundRuleCFStates& HundRuleCFStates::operator = (HundRuleCFStates &toCopy)
{
  if ((this->NbrParticles != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] this->NbrTermsPerLz;
      for (int i=0; i<=TotalL; ++i)
	delete [] this->TermsPerLz[i];
      delete [] this->TermsPerLz;
      delete [] this->Ji;
      delete this->Orbitals;
      if (this->LEigenstates!=0)
	delete [] this->LEigenstates;
    }
  
  this->NbrParticles=toCopy.NbrParticles;
  this->NumShells=toCopy.NumShells;
  this->NbrOrbitalsInHighestShell=toCopy.NbrOrbitalsInHighestShell;
  this->NbrParticlesInHighestShell=toCopy.NbrParticlesInHighestShell;
  this->HighestShellLzMax=toCopy.HighestShellLzMax;
  this->TotalL=toCopy.TotalL;
  this->ReverseFluxFlag=toCopy.ReverseFluxFlag;
  this->TwoS=toCopy.TwoS;
  this->TwoSEff=toCopy.TwoSEff;
  this->JastrowP=toCopy.JastrowP;
  this->Orbitals=toCopy.Orbitals;
  this->Ji=toCopy.Ji;
  this->LEigenstates=toCopy.LEigenstates;
  this->ElementNorm=toCopy.ElementNorm;
  this->NbrTermsPerLz=toCopy.NbrTermsPerLz;
  this->TermsPerLz=toCopy.TermsPerLz;
  this->SelectMPosition=toCopy.SelectMPosition;
#ifdef __USE_LAPACK_HERE__
  this->SlaterDeterminant.Resize(NbrParticles);
#else
  this->SlaterDeterminant.Resize(NbrParticles, NbrParticles);
#endif  
  this->Flag=toCopy.Flag;
  
  return *this;
}

  // destructor
HundRuleCFStates::~HundRuleCFStates()
{
  if ((this->NbrParticles != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] this->NbrTermsPerLz;
      for (int i=0; i<=TotalL; ++i)
	delete [] this->TermsPerLz[i];
      delete [] this->TermsPerLz;
      delete this->Orbitals;
      delete [] this->Ji;
      if (this->LEigenstates!=0)
	delete [] this->LEigenstates;
    }
}
  

// clone function 
//
// return value = clone of the function
// attention, copy constructor not properly defined, yet.
Abstract1DComplexFunction* HundRuleCFStates::Clone ()
{
  return (Abstract1DComplexFunction*)( new HundRuleCFStates(*this));
}


// select Lz of state returned with operator()
// return value=true upon success (M in bounds)
bool HundRuleCFStates::SelectMValue(int newMValue)
{
  //cout << "TotalL="<<TotalL<<", newM="<<newMValue<<endl;
  if ((newMValue>=-this->TotalL)&&(newMValue<=this->TotalL)&&(((newMValue+this->TotalL)&1)==0))
    {
      cout << "Selected Hund's rule state with M=";
      if (newMValue&1)
	cout<<newMValue<<"/2"<<endl;
      else
	cout<<newMValue/2<<endl;
      this->SelectMPosition=(newMValue+this->TotalL)/2;
      return true;
    }
  else
    {
      return false;
    }
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  
Complex HundRuleCFStates::operator ()(RealVector& x)
{
  this->OrbitalValues = (*(this->Orbitals))(x);
  this->EvaluateTables();
  return EvaluateAState(SelectMPosition);
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex HundRuleCFStates::CalculateFromSpinorVariables(ComplexVector& uv)
{
  this->OrbitalValues = this->Orbitals->CalculateFromSpinorVariables(uv);
  this->EvaluateTables();
  return EvaluateAState(SelectMPosition);
}


void HundRuleCFStates::GetValues(RealVector& x, Complex *result)
{
  this->OrbitalValues = (*(this->Orbitals))(x);
  this->EvaluateTables();  
  for (int m=0; m<=this->TotalL; ++m)
    result[m]=EvaluateAState(m);
}

void HundRuleCFStates::GetValuesFromSpinorVariables(ComplexVector& uv, Complex *result)
{
  this->OrbitalValues = this->Orbitals->CalculateFromSpinorVariables(uv);
  this->EvaluateTables();  
  for (int m=0; m<=this->TotalL; ++m)
    result[m]=EvaluateAState(m);
}

// add up the Slater determinants for a single state:
Complex HundRuleCFStates::EvaluateAState(int LzPosition)
{
  Complex Result = 0.0;
  SlaterComponent *TmpTerms;
  TmpTerms = TermsPerLz[LzPosition];
  //cout << "Calculating M="<<(2*LzPosition-HighestShellLzMax)/2.0<<" (LzPosition="<<LzPosition<<")"<<endl;
  for (int t=0; t<this->NbrTermsPerLz[LzPosition]; ++t)
    {
      for (int alpha=0; alpha<NbrParticles-NbrParticlesInHighestShell; ++alpha)
	for (int n=0; n<NbrParticles; ++n)
	  {
	    //cout << "Setting element of lower shells (" << alpha<<", "<<n<<")"<<endl;
#ifdef __USE_LAPACK_HERE__
	    this->SlaterDeterminant.SetMatrixElement(alpha,n,Ji[n]*this->ElementNorm*OrbitalValues[alpha][n]);
#else
	    this->SlaterDeterminant[alpha].Re(n) = Real(Ji[n]*this->ElementNorm*OrbitalValues[alpha][n]);
	    this->SlaterDeterminant[alpha].Im(n) = Imag(Ji[n]*this->ElementNorm*OrbitalValues[alpha][n]);
#endif
	  }
      int *TmpLzPositions=TmpTerms[t].GetLzPositionPtr();
      int alpha0 = NbrParticles-NbrParticlesInHighestShell; 
      int alpha1, alpha2 = NbrParticles-NbrParticlesInHighestShell;
      for (int i=0; i<NbrParticlesInHighestShell; ++i)
	{
	  alpha1=alpha0 + TmpLzPositions[i];
	  for (int n=0; n<NbrParticles; ++n)
	    {
	      //cout << "Setting element of highest shells (" << alpha2<<", "<<n<<")"<<endl;
#ifdef __USE_LAPACK_HERE__
	      this->SlaterDeterminant.SetMatrixElement(alpha2,n,Ji[n]*this->ElementNorm*OrbitalValues[alpha1][n]);
#else
	      this->SlaterDeterminant[alpha2].Re(n) = Real(Ji[n]*this->ElementNorm*OrbitalValues[alpha1][n]);
	      this->SlaterDeterminant[alpha2].Im(n) = Imag(Ji[n]*this->ElementNorm*OrbitalValues[alpha1][n]);
#endif
	    }
	  ++alpha2;
	}
      Complex SlaterDet = this->SlaterDeterminant.Determinant();
      SlaterDet*=TmpTerms[t].GetPrefactor()*this->Interpolation;
      Result+=SlaterDet;
    }
  return Result;
}


// set wavefunction to one for a given set of particle coordinates
void HundRuleCFStates::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)1.0/this->NbrParticles);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)1.0/this->NbrParticles);
      else 
	this->ElementNorm*= pow(det,(double)-1.0/this->NbrParticles);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}



// utility function to set the right dynamic interval for Monte-Carlo
void HundRuleCFStates::AdaptAverageMCNorm(int thermalize, int average)
{
  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(2*this->NbrParticles);
  this->AdaptNorm(Particles->GetPositions());
  Complex TmpMetropolis, TrialValue = (*this)(Particles->GetPositions());  
  double PreviousSamplingAmplitude = SqrNorm(TrialValue);
  double CurrentSamplingAmplitude = PreviousSamplingAmplitude;
  int NextCoordinates=0;
  // do some MC moves: accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
  for (int i = 0; i < thermalize; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticles) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticles) --NextCoordinates;      
    }
  this->AdaptNorm(Particles->GetPositions());
  double SumTrialValues=0.0;
  for (int i = 0; i < average; ++i)
    {
      Particles->Move(NextCoordinates);
      TmpMetropolis = (*this)(Particles->GetPositions());
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((Particles->GetRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      NextCoordinates = (int) (((double) NbrParticles) * Particles->GetRandomNumber());
      if (NextCoordinates == NbrParticles) --NextCoordinates;            
      SumTrialValues+=Norm(TmpMetropolis);
    }  
  this->ElementNorm*= pow(SumTrialValues/average,(double)-2.0/this->NbrParticles);
  delete Particles;
}


// evaluate all precalculations
void HundRuleCFStates::EvaluateTables()
// assumes OrbitalValues initialized
//
{
  Complex tmp;

  // evaluate single particle Jastrow factors
  this->Interpolation=1.0;
  if (Orbitals->TestCriticality(Interpolation) == 0)
    {
      for (int i=0;i<this->NbrParticles;i++)
	{
	  Ji[i]=1.0;
	  for(int j=0;j<i;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
	  for(int j=i+1;j<NbrParticles;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
	  tmp = Ji[i];
	  for (int i=1;i<JastrowP;++i)
	    Ji[i]*=tmp;
	}
    }
  else // if some interpolation occurred, the true values of the Ji's have to be recalculated:
    {
      for (int i=0;i<this->NbrParticles;i++)
	{
	  Ji[i]=1.0;
	  for(int j=0;j<i;j++) Ji[i] *= ((Orbitals->SpinorU(i) * Orbitals->SpinorV(j)) - (Orbitals->SpinorU(j) * Orbitals->SpinorV(i)));
	  for(int j=i+1;j<NbrParticles;j++) Ji[i] *= ((Orbitals->SpinorU(i) * Orbitals->SpinorV(j)) - (Orbitals->SpinorU(j) * Orbitals->SpinorV(i)));
	  for (int i=1;i<JastrowP;++i)
	    Ji[i]*=tmp;
	}
    }
}
