#include "ParticleOnSphereCollectionSouthPole.h"

#include <iostream>
using std::cout;
using std::endl;

// switch for debugging output:
//#define DEBUG_OUTPUT

// numerical constants: square root of 1/2
#define INV_SQRT2 0.707106781186548

ParticleOnSphereCollectionSouthPole::ParticleOnSphereCollectionSouthPole()
{
  this->NbrParticles = 0;
}

ParticleOnSphereCollectionSouthPole::ParticleOnSphereCollectionSouthPole(int N, long seed)
{
  this->NbrParticles = N;
  this->SpinorUCoordinates = new Complex[N];
  this->SpinorVCoordinates = new Complex[N];
  this->N1 = new double[N];
  this->N2 = new double[N];
  this->N3 = new double[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = new NumRecRandomGenerator(seed);
  ExternalGenerator = false;
  this->Theta0 = (0.2*M_PI/sqrt((double)N/10.0));
  this->ThetaPhi.Resize(2*N);
  this->Randomize();

  this->Distances = new RealSymmetricMatrix(NbrParticles, true);
  this->DistancesUpToDate = false;
}

ParticleOnSphereCollectionSouthPole::ParticleOnSphereCollectionSouthPole(int N, AbstractRandomNumberGenerator *generator)
{
  this->NbrParticles = N;
  this->SpinorUCoordinates = new Complex[N];
  this->SpinorVCoordinates = new Complex[N];
  this->N1 = new double[N];
  this->N2 = new double[N];
  this->N3 = new double[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = generator;
  this->ExternalGenerator = true;
  this->Theta0 = (0.2*M_PI/sqrt((double)N/10.0));
  this->ThetaPhi.Resize(2*N);  
  this->Randomize();

  this->Distances = new RealSymmetricMatrix(NbrParticles, true);
  this->DistancesUpToDate = false;

}


ParticleOnSphereCollectionSouthPole::ParticleOnSphereCollectionSouthPole(const ParticleOnSphereCollectionSouthPole &tocopy)
{
  this->NbrParticles = tocopy.NbrParticles;
  this->LastMoved = tocopy.LastMoved;
  this->LastU = tocopy.LastU;
  this->LastV = tocopy.LastV;
  this->LastTheta = tocopy.LastTheta;
  this->LastPhi = tocopy.LastPhi;
  this->SpinorUCoordinates = tocopy.SpinorUCoordinates;
  this->SpinorVCoordinates = tocopy.SpinorVCoordinates;
  this->N1 = tocopy.N1;
  this->N2 = tocopy.N2;
  this->N3 = tocopy.N3;
  this->Flag  = tocopy.Flag;
  this->Generator = tocopy.Generator;
  this->ExternalGenerator = tocopy.ExternalGenerator;
  this->ThetaPhi = tocopy.ThetaPhi;
  this->Theta0 = tocopy.Theta0;
  this->Distances = tocopy.Distances;
  this->DistancesUpToDate = tocopy.DistancesUpToDate;
}

ParticleOnSphereCollectionSouthPole::~ParticleOnSphereCollectionSouthPole()
{
  if ((this->Flag.Used() == true)&&(this->Flag.Shared() == false))
    {
      delete [] this->SpinorUCoordinates;
      delete [] this->SpinorVCoordinates;
      delete [] this->N1;
      delete [] this->N2;
      delete [] this->N3;
      if (!this->ExternalGenerator)
	delete Generator;
      delete [] Distances;
    }
}

// randomly moves particle number nbrParticle
void ParticleOnSphereCollectionSouthPole::Move(int nbrParticle)
{
  // store old positions
  this->LastMoved = nbrParticle;
  this->LastU = SpinorUCoordinates[nbrParticle].Re;
  this->LastV = SpinorVCoordinates[nbrParticle];
  this->LastN1 = N1[nbrParticle];
  this->LastN2 = N2[nbrParticle];
  this->LastN3 = N3[nbrParticle];
  this->LastTheta = this->ThetaPhi[nbrParticle<<1];
  this->LastPhi = this->ThetaPhi[(nbrParticle<<1)+1];
  // make a random move with gaussian distribution around the initial position
  double theta = this->Theta0*fabs(Generator->GetGaussianRandomNumber());
  double phi = 2.0*M_PI*Generator->GetRealRandomNumber();
  
  double st = sin(theta/2.0);
  double ct = cos(theta/2.0);
  double sp = sin(phi);
  double cp = cos(phi);
  double z = st*cp;
  Complex U(LastN2, LastN3);
  Complex W(ct,-st*cp);

  N1[nbrParticle] = LastN1*(ct*ct + st*st*(2.0*sp*sp-1.0));
  Complex Tmp=U*W*W - (2.0*LastN1*z)*Complex(W.Im,-W.Re) + (z*z)*Complex(LastN2,-LastN3);
  N2[nbrParticle] = Tmp.Re;
  N3[nbrParticle] = Tmp.Im;
  z=Tmp.Im*Tmp.Im;
  if(z < 0.99999)
    {
      double root=sqrt((1.-z)/(1.-Tmp.Im));
      SpinorUCoordinates[nbrParticle].Re = INV_SQRT2 * root;
      SpinorVCoordinates[nbrParticle] = (INV_SQRT2 / root) * Complex(N1[nbrParticle],Tmp.Re);
    }
  else
    {
      if (Tmp.Im>0)
	{
	  double eps = 1.0-Tmp.Im;
	  SpinorUCoordinates[nbrParticle].Re = sqrt(1.0-0.5*eps);
	  SpinorVCoordinates[nbrParticle] = Complex(N1[nbrParticle],Tmp.Re) * (0.5/SpinorUCoordinates[nbrParticle].Re);
	}
      else
	{
	  double eps = -1.0-Tmp.Im;
	  SpinorUCoordinates[nbrParticle].Re = sqrt(0.5*(1.-z)/(1.-Tmp.Im));
	  Tmp=Complex(N1[nbrParticle],Tmp.Re);
	  Tmp/=Norm(Tmp);
	  SpinorVCoordinates[nbrParticle] = sqrt(1.0-(eps*(1.0-0.5*eps)/(2.0-eps))) * Tmp;
	}
    }
  this->ThetaPhi[nbrParticle<<1] = 2.0*acos(fabs(SpinorUCoordinates[nbrParticle].Re));
  this->ThetaPhi[(nbrParticle<<1)+1] = Arg(SpinorVCoordinates[nbrParticle]);
  cout << "Checking norm: "<<N1[nbrParticle]*N1[nbrParticle]+N2[nbrParticle]*N2[nbrParticle]+N3[nbrParticle]*N3[nbrParticle]<<" "<<SqrNorm(SpinorUCoordinates[nbrParticle])+SqrNorm(SpinorVCoordinates[nbrParticle])<<endl;
  // cout << "new coordinates: ("<<this->ThetaPhi[nbrParticle<<1]<<", " << this->ThetaPhi[(nbrParticle<<1)+1] << ")" <<endl;
}

// randomly select a particle and move it
int ParticleOnSphereCollectionSouthPole::Move()
{
  this->LastMoved = (int) (((double) NbrParticles) * Generator->GetRealRandomNumber());
  if (LastMoved == NbrParticles) --LastMoved;
  this->Move(LastMoved);
  return LastMoved;
}


// randomly moves all particles 
void ParticleOnSphereCollectionSouthPole::RotateAll()
{
  double theta = this->Theta0*fabs(Generator->GetGaussianRandomNumber());
  double phi = 2.0*M_PI*Generator->GetRealRandomNumber();
  this->RotateAll(theta, phi);
}

// moves all particles by angles theta and phi
void ParticleOnSphereCollectionSouthPole::RotateAll(double theta, double phi)
{
  cout << "Attention: ParticleOnSphereCollectionSouthPole::RotateAll needs to be implemented!"<<endl;
  /*
  Complex move_v = sin(theta/2.0)*Polar( phi/2.0);
  Complex move_u = cos(theta/2.0)*Polar(-phi/2.0);
  Complex lastU, lastV;
  for (int nbrParticle=0; nbrParticle<NbrParticles; ++nbrParticle)
    {
      lastU=SpinorUCoordinates[nbrParticle];
      lastV=SpinorVCoordinates[nbrParticle];
      SpinorUCoordinates[nbrParticle] =lastU*Conj(move_u)+lastV*Conj(move_v);
      SpinorVCoordinates[nbrParticle] =-lastU*move_v+lastV*move_u;
      this->ThetaPhi[nbrParticle<<1] = (2.0*acos(Norm(SpinorUCoordinates[nbrParticle])));
      this->ThetaPhi[(nbrParticle<<1)+1] = Arg(SpinorVCoordinates[nbrParticle]);
      // cout << "new coordinates: ("<<this->ThetaPhi[nbrParticle<<1]<<", " << this->ThetaPhi[(nbrParticle<<1)+1] << ")" <<endl;
    }

  // rotate old positions, as well
  lastU=this->LastU;
  lastV=this->LastV;
  this->LastU=lastU*Conj(move_u)+lastV*Conj(move_v);
  this->LastV=-lastU*move_v+lastV*move_u;
  this->LastTheta = (2.0*acos(Norm(LastU)));
  this->LastPhi = (Arg(LastV)-Arg(LastU));
  */
}


// restore last move
void ParticleOnSphereCollectionSouthPole::RestoreMove()
{
  if (LastMoved>-1)
    {      
      this->SpinorUCoordinates[LastMoved].Re = this->LastU;
      this->SpinorVCoordinates[LastMoved] = this->LastV;
      this->N1[LastMoved] = this->LastN1;
      this->N2[LastMoved] = this->LastN2;
      this->N3[LastMoved] = this->LastN3;
      this->ThetaPhi[LastMoved<<1] = this->LastTheta;
      this->ThetaPhi[(LastMoved<<1)+1] = this->LastPhi;
      this->LastMoved = -1;
    }
}


void ParticleOnSphereCollectionSouthPole::MultiplyStepLength(double multiplier)
{
  this->Theta0*=multiplier;
}

// get single position theta
double ParticleOnSphereCollectionSouthPole::Theta(int nbrParticle)
{
  return this->ThetaPhi[nbrParticle<<1];
}


// get single position phi
double ParticleOnSphereCollectionSouthPole::Phi(int nbrParticle)
{
  return this->ThetaPhi[(nbrParticle<<1)+1];
}

// set single position
void ParticleOnSphereCollectionSouthPole::SetPosition(int nbrParticle, double theta, double phi)
{
  double s2=sin(theta/2.);
  double s=sin(theta);
  double c=cos(theta);
  double sp=sin(phi);
  double cp=cos(phi);
  this->SpinorUCoordinates[nbrParticle].Re = cos(theta/2.0);
  this->SpinorUCoordinates[nbrParticle].Im = 0.0;
  this->SpinorVCoordinates[nbrParticle].Re = s2*cp;
  this->SpinorVCoordinates[nbrParticle].Im = s2*sp;
  this->N1[nbrParticle]=s*cp;
  this->N2[nbrParticle]=s*sp;
  this->N3[nbrParticle]=c;
  this->ThetaPhi[nbrParticle<<1] = theta;   // (2.0*acos(Norm(SpinorUCoordinates[nbrParticle])));
  this->ThetaPhi[(nbrParticle<<1)+1] = phi; // (Arg(SpinorVCoordinates[nbrParticle])-Arg(SpinorUCoordinates[nbrParticle]));
}

double ParticleOnSphereCollectionSouthPole::GetRandomNumber()
{
  return Generator->GetRealRandomNumber();
}

// randomize particle positions
void ParticleOnSphereCollectionSouthPole::Randomize()
{
  double phi, theta,  c, s, s2, sp, cp;
  for (int i=0; i<NbrParticles; ++i)
    {
      phi=Generator->GetRealRandomNumber()*2.0*M_PI;
      theta=acos(1.0-2.0*Generator->GetRealRandomNumber());
      s2=sin(theta/2.);
      s=sin(theta);
      c=cos(theta);
      sp=sin(phi);
      cp=cos(phi);
      this->SpinorUCoordinates[i].Re = cos(theta/2.);
      this->SpinorUCoordinates[i].Im = 0.0;    // gauge choice!
      this->SpinorVCoordinates[i].Re = s2*cp;
      this->SpinorVCoordinates[i].Im = s2*sp;
      this->N1[i]=s*cp;
      this->N2[i]=s*sp;
      this->N3[i]=c;
      this->ThetaPhi[i<<1] = theta;   // (2.0*acos(Norm(SpinorUCoordinates[i])));
      this->ThetaPhi[(i<<1)+1] = phi; // (Arg(SpinorVCoordinates[i])-Arg(SpinorUCoordinates[i]));
    }
}

// get absolute values of all relative distances
// distances = matrix in which to return the distances
void ParticleOnSphereCollectionSouthPole::GetDistances(RealSymmetricMatrix &distances)
{
  if ((distances.GetNbrRow()!=NbrParticles)||(distances.GetNbrColumn()!=NbrParticles))
    distances.Resize(NbrParticles,NbrParticles);
  for (int i=0; i<NbrParticles; ++i)
    {
      distances(i,i)=0.0;
      for (int j=i+1; j<NbrParticles; ++j)
	{
	  distances(i,j)=(*this->Distances)(i,j)=Norm(this->SpinorUCoordinates[i].Re*this->SpinorVCoordinates[j]-this->SpinorUCoordinates[j].Re*this->SpinorVCoordinates[i]);
	}
    }
  this->DistancesUpToDate = true;
  return;
}

// get reference to internal matrix with particle distances
const RealSymmetricMatrix& ParticleOnSphereCollectionSouthPole::GetDistances()
{
  if (!DistancesUpToDate)
    {
      for (int i=0; i<NbrParticles; ++i)
	{
	  // this->Distances(i,i)=0.0;
	  for (int j=i+1; j<NbrParticles; ++j)
	    {
	      (*this->Distances)(i,j)=Norm(this->SpinorUCoordinates[i]*this->SpinorVCoordinates[j]-this->SpinorUCoordinates[j]*this->SpinorVCoordinates[i]);
	    }
	}
    }
  this->DistancesUpToDate = false;
  return *this->Distances;
}



// toggle positions of first N/2 particles with the remaining N/2 positions
//
void ParticleOnSphereCollectionSouthPole::ToggleHalfHalf()
{
  if (NbrParticles&1) return;
  double TmpD;
  Complex TmpC;
  int NUp = this->NbrParticles/2;
  // exchange spin up and spin down
  for (int j = 0; j < NUp; ++j)
    {
      TmpD = N1[j];
      N1[j] = N1[j+NUp];
      N1[j+NUp] = TmpD;
      
      TmpD = N2[j];
      N2[j] = N2[j+NUp];
      N2[j+NUp] = TmpD;
      
      TmpD = N3[j];
      N3[j] = N3[j+NUp];
      N3[j+NUp] = TmpD;
      
      TmpD = ThetaPhi[j << 1];      
      ThetaPhi[j << 1] = ThetaPhi[(j+NUp) << 1];
      ThetaPhi[(j+NUp) << 1] = TmpD;
      
      TmpD = ThetaPhi[1 + (j << 1)];
      ThetaPhi[1+(j <<1)] = ThetaPhi[1+ ((j+NUp) << 1)];
      ThetaPhi[1+ ((j+NUp) << 1)] = TmpD;
      
      TmpD = SpinorUCoordinates[j].Re;
      SpinorUCoordinates[j].Re = SpinorUCoordinates[j+NUp].Re;
      SpinorUCoordinates[j+NUp].Re = TmpD;
      
      TmpC = SpinorVCoordinates[j];
      SpinorVCoordinates[j] = SpinorVCoordinates[j+NUp];
      SpinorVCoordinates[j+NUp] = TmpC;
    }
}
