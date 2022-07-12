#include "ParticleOnDiskCollection.h"

// switch for debugging output:
//#define DEBUG_OUTPUT


ParticleOnDiskCollection::ParticleOnDiskCollection()
{
  this->NbrParticles = 0;
}

/// create a collection of particles for a circular droplet of N particles at a density nu
/// seed = optional random seed
ParticleOnDiskCollection::ParticleOnDiskCollection(int N, double nu, long seed)
{
  if(nu==0.0)
    {
      std::cout << "Warning, assuming nu=1 for initialisation!"<<std::endl;
      nu=1.0;
    }
  this->NbrParticles = N;
  this->Positions.Resize(2*N); // real-vector holding positions as (x1,y1, x2, y2, ...)
  this->CoordinatesZ = (Complex*)(&(this->Positions[0])); // new Complex[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = new NumRecRandomGenerator(seed);
  ExternalGenerator = false;
  if (nu<=0.0)
    this->Nu = 1.0;
  else
    this->Nu = nu;
  this->R0 = (1.0/sqrt(nu));
  this->Randomize();

}

ParticleOnDiskCollection::ParticleOnDiskCollection(int N, double nu, AbstractRandomNumberGenerator *generator)
{
  if(nu==0.0)
    {
      std::cout << "Warning, assuming nu=1 for initialisation!"<<std::endl;
      nu=1.0;
    }
  this->NbrParticles = N;
  this->Positions.Resize(2*N); // real-vector holding positions as (x1,y1, x2, y2, ...)
  this->CoordinatesZ = (Complex*)(&(this->Positions[0])); // new Complex[N];
  this->LastMoved=-1;
  this->Flag.Initialize();
  this->Generator = generator;
  ExternalGenerator = true;
  if (nu<=0.0)
    this->Nu = 1.0;
  else
    this->Nu = nu;
  this->R0 = (0.5/sqrt(nu));
  this->Randomize();
}



ParticleOnDiskCollection::ParticleOnDiskCollection(const ParticleOnDiskCollection &tocopy)
{
  this->NbrParticles = tocopy.NbrParticles;
  this->LastMoved = tocopy.LastMoved;
  this->LastZ = tocopy.LastZ;
  // this->LastR = tocopy.LastR;
  // this->LastPhi = tocopy.LastPhi;
  this->CoordinatesZ = tocopy.CoordinatesZ;
  this->Flag  = tocopy.Flag;
  this->Generator = tocopy.Generator;
  this->ExternalGenerator = tocopy.ExternalGenerator;
  this->Positions = tocopy.Positions;
  this->R0 = tocopy.R0;
}

ParticleOnDiskCollection::~ParticleOnDiskCollection()
{
  if ((this->Flag.Used() == true)&&(this->Flag.Shared() == false))
    {
      // delete [] CoordinatesZ;
      if (!this->ExternalGenerator)
	delete Generator;
    }
}

// randomly moves particle number nbrParticle
void ParticleOnDiskCollection::Move(int nbrParticle)
{
  // store old positions
  this->LastMoved = nbrParticle;
  this->LastZ = CoordinatesZ[nbrParticle];
  // make a random move with gaussian distribution around the initial position
  Complex move_dz = Polar( this->R0*Generator->GetGaussianRandomNumber(), M_PI*Generator->GetRealRandomNumber());
  CoordinatesZ[nbrParticle] += move_dz;
  // std::cout << "new coordinates after move: "<<this->CoordinatesZ[nbrParticle]<<std::endl;
}

// randomly select a particle and move it
int ParticleOnDiskCollection::Move()
{
  this->LastMoved = (int) (this->NbrParticles * Generator->GetRealRandomNumber());
  if (this->LastMoved == this->NbrParticles) --LastMoved;
  this->Move(LastMoved);
  return LastMoved;
}



// restore last move
void ParticleOnDiskCollection::RestoreMove()
{
  if (LastMoved>-1)
    {
      Complex tmp = CoordinatesZ[LastMoved];
      CoordinatesZ[LastMoved] = this->LastZ;
      this->LastZ = tmp;
      // this->RPhi[LastMoved<<1] = this->LastR;
      // this->RPhi[(LastMoved<<1)+1] = this->LastPhi;
    }
}


void ParticleOnDiskCollection::MultiplyStepLength(double multiplier)
{
  this->R0*=multiplier;
}

// get single position theta
double ParticleOnDiskCollection::GetR(int nbrParticle)
{
  return Norm(this->CoordinatesZ[nbrParticle]);
}


// get single position phi
double ParticleOnDiskCollection::Phi(int nbrParticle)
{
  return Arg(this->CoordinatesZ[nbrParticle]);
  //  return this->RPhi[(nbrParticle<<1)+1];
}

void ParticleOnDiskCollection::SetPosition(int nbrParticle, double x, double y)
{
  if(nbrParticle < 0 || nbrParticle >= this->NbrParticles) return;
  this->Positions[2*nbrParticle] = x;
  this->Positions[2*nbrParticle+1] = y;
}

double ParticleOnDiskCollection::GetRandomNumber()
{
  return Generator->GetRealRandomNumber();
}

// randomize particle positions
void ParticleOnDiskCollection::Randomize()
{
  double Sqr, Rmax = sqrt((double)this->NbrParticles/(M_PI*this->Nu)), SumSqr=0.0;
  for (int i=0; i<this->NbrParticles; ++i)
    {      
      do
	{
	  this->Positions[2*i] = (1.0-2.0*Generator->GetRealRandomNumber());
	  this->Positions[2*i+1] = (1.0-2.0*Generator->GetRealRandomNumber());
	  Sqr=SqrNorm(this->CoordinatesZ[i]);
	} while (Sqr>1.0);
      this->CoordinatesZ[i]*=Rmax;
      SumSqr+=SqrNorm(this->CoordinatesZ[i]);
    }

  // std::cout << "Randomized coordinates: \n";
  // std::cout << "Rmax = "<<Rmax<< ", nu=" << this->Nu << ", sum_i |z|**2 = "<<SumSqr<<"\n";
  // this->PrintPositions(std::cout);
}

// get absolute values of all relative distances
// distances = matrix in which to return the distances
void ParticleOnDiskCollection::GetDistances(RealSymmetricMatrix &distances)
{
  if ((distances.GetNbrRow()!=NbrParticles)||(distances.GetNbrColumn()!=NbrParticles))
    distances.Resize(NbrParticles,NbrParticles);
  for (int i=0; i<NbrParticles; ++i)
    {
      distances(i,i)=0.0;
      for (int j=i+1; j<NbrParticles; ++j)
	{
	  distances(i,j)=Norm(this->CoordinatesZ[i]-this->CoordinatesZ[j]);
	}
    }
  return;
}


// toggle positions of first N/2 particles with the remaining N/2 positions
//
void ParticleOnDiskCollection::ToggleHalfHalf()
{
  if (NbrParticles&1) return;
  double TmpD;
  Complex TmpC;
  int NUp = this->NbrParticles/2;
  // exchange spin up and spin down
  for (int j = 0; j < NUp; ++j)
    {
      // TmpD = RPhi[j << 1];      
      // RPhi[j << 1] = RPhi[(j+NUp) << 1];
      // RPhi[(j+NUp) << 1] = TmpD;
      // TmpD = RPhi[1 + (j << 1)];
      // RPhi[1+(j <<1)] = RPhi[1+ ((j+NUp) << 1)];
      // RPhi[1+ ((j+NUp) << 1)] = TmpD;
      TmpC = CoordinatesZ[j];
      CoordinatesZ[j] = CoordinatesZ[j+NUp];
      CoordinatesZ[j+NUp] = TmpC;
    }
}


// print positions (mainly for testing)
ostream &ParticleOnDiskCollection::PrintPositions(ostream &Str)
{
  for (int i=0; i<NbrParticles; ++i)
    {
      Str << "r_" << i << " = ("<<Positions[2*i]<<","<<Positions[2*i+1] << ") = " << CoordinatesZ[i] << "\n";
     }
  return Str;
}
