#include "Options/Options.h"

#include "MathTools/RandomNumber/AbstractRandomNumberGenerator.h"
#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEDiskMonteCarloFluctuations" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-samples", "number of Monte Carlo samples", 1e8);
  (*SystemGroup) += new SingleIntegerOption ('\n', "init-samples", "number of Monte Carlo samples before starting to record", 1e5);
  (*SystemGroup) += new SingleIntegerOption ('\n', "measure-frequency", "number of steps between measures", 100);
  (*SystemGroup) += new SingleIntegerOption ('q', "laughlin-q", "invert of the Laughlin filling fraction", 1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "discretize-disk", "measure fluctuation of number of particles in all sections of angle theta = i / (p 2pi), i < p", 100);
  (*SystemGroup) += new SingleDoubleOption ('\n', "ellipse-t", "number of particles", 0.0);
   (*SystemGroup) += new SingleStringOption ('\n', "file-name", "string of character to characterize the run", "0");
   (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
   
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinTimeReversalSymmetryAndPairing -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LaughlinInvertFillingFraction = Manager.GetInteger("laughlin-q");
  int Shift = 0;
  int NbrOccupiedOrbitals = NbrParticles * LaughlinInvertFillingFraction - Shift;
  int NbrStepInitialization = Manager.GetInteger("init-samples");
  int MeasurementPeriod = Manager.GetInteger("measure-frequency");
  int NbrSampleSizes = Manager.GetInteger("discretize-disk");
  double SampleSize = 2.0 * M_PI / ((double) NbrSampleSizes);
  int Tmp;
  long NumberSamples = Manager.GetInteger("nbr-samples");
  double EllipseParameter = Manager.GetDouble("ellipse-t");
  double RadiusSquareInsideDroplet = (NbrParticles * LaughlinInvertFillingFraction * (1.0 - EllipseParameter))/ (2.0 * (1.0 + EllipseParameter));
  
  
  double*** MomentOrderNbrParticles = new double** [2];
  for (int l = 0; l < 2; ++l)
  {
    MomentOrderNbrParticles[l] = new double*[4];
    for (int i = 0; i < 4; ++i)
    {
      MomentOrderNbrParticles[l][i] = new double [NbrSampleSizes];
      for (int j = 0; j < NbrSampleSizes; ++j)
	MomentOrderNbrParticles[l][i][j] = 0;
    }
  }
  
  int** CurrentRunNumberParticles = new int*[2];
  for (int l = 0; l < 2; ++l)
    CurrentRunNumberParticles[l] = new int[NbrSampleSizes];
    
  
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);
  StdlibRandomNumberGenerator* RandomNumberGenerator = 0;
  RandomNumberGenerator = new StdlibRandomNumberGenerator(TotalStartingTime.tv_sec);
  
  double* ParticleCoordinateX = new double[NbrParticles];
  double* ParticleCoordinateY = new double[NbrParticles];
  for (int i = 0; i < NbrParticles; ++i)
  {
    double TmpRadius = sqrt(-2.0 * log (1.0 - RandomNumberGenerator->GetRealRandomNumber())) + 0.3;
    double TmpPhi = 2.0 * M_PI * RandomNumberGenerator->GetRealRandomNumber();
    double XCoordinate = TmpRadius * cos(TmpPhi) * sqrt(NbrParticles);
    double YCoordinate = TmpRadius * sin(TmpPhi) * sqrt(NbrParticles);
    ParticleCoordinateX[i] = XCoordinate;
    ParticleCoordinateY[i] = YCoordinate;    
  }
  
  int TmpParticleIndex;
  double TransferR;
  double TransferPhi;
  double RadiusSquare;
  double Phi;
  double TmpNextX;
  double TmpNextY;
  double TmpRatio;
  double TmpNextSquareNorm;
  double TmpCurrentSquareNorm;
  
  for (long i = 0 ; i < (NbrStepInitialization + Manager.GetInteger("nbr-samples")); ++i)
  {
    TmpParticleIndex = (RandomNumberGenerator->GetIntegerRandomNumber()) % NbrParticles;
    
    TransferR = - 0.5 * log(1.0 - RandomNumberGenerator->GetRealRandomNumber());
    TransferPhi = 2.0 * M_PI * RandomNumberGenerator->GetRealRandomNumber();
    
    TmpNextX = ParticleCoordinateX[TmpParticleIndex] + TransferR * cos(TransferPhi);
    TmpNextY = ParticleCoordinateY[TmpParticleIndex] + TransferR * sin(TransferPhi);
    
    TmpRatio = 1.0;
    for (int j = 0; j < NbrParticles; ++j)
    {
      if (j != TmpParticleIndex)
      {
	TmpNextSquareNorm = (TmpNextX - ParticleCoordinateX[j])*(TmpNextX - ParticleCoordinateX[j]) + (TmpNextY - ParticleCoordinateY[j])*(TmpNextY - ParticleCoordinateY[j]);
	TmpCurrentSquareNorm = (ParticleCoordinateX[TmpParticleIndex] - ParticleCoordinateX[j])*(ParticleCoordinateX[TmpParticleIndex] - ParticleCoordinateX[j]) + (ParticleCoordinateY[TmpParticleIndex] - ParticleCoordinateY[j])*(ParticleCoordinateY[TmpParticleIndex] - ParticleCoordinateY[j]);
	
	TmpRatio *= (TmpNextSquareNorm / TmpCurrentSquareNorm);
      }
    }
   
    TmpRatio = pow(TmpRatio, LaughlinInvertFillingFraction) * exp(-((1.0 + EllipseParameter) * (TmpNextX*TmpNextX - ParticleCoordinateX[TmpParticleIndex]*ParticleCoordinateX[TmpParticleIndex]) + (1.0 - EllipseParameter) * (TmpNextY*TmpNextY  - ParticleCoordinateY[TmpParticleIndex]*ParticleCoordinateY[TmpParticleIndex])));
    
    //Metropolis algorithm: decide or not to accept move
    if (RandomNumberGenerator->GetRealRandomNumber() < TmpRatio)
    {
      ParticleCoordinateX[TmpParticleIndex] = TmpNextX;
      ParticleCoordinateY[TmpParticleIndex] = TmpNextY;
    }
    
    if ((i >= NbrStepInitialization) && ((i % MeasurementPeriod) == 0))
    {
      for (int l = 0; l < 2; ++l)
	for (int p = 0; p < NbrSampleSizes; ++p)
	  CurrentRunNumberParticles[l][p] = 0;
      for (int j = 0; j < NbrParticles; ++j)
      {
	RadiusSquare = ParticleCoordinateX[j]*ParticleCoordinateX[j] + ParticleCoordinateY[j]*ParticleCoordinateY[j];
	Phi = atan2(ParticleCoordinateY[j], ParticleCoordinateX[j]) + M_PI;
	Tmp = 1;
	while ((Tmp < NbrSampleSizes) && (Phi  > ((double) Tmp) * SampleSize))
	  ++Tmp;
	for (int p = Tmp - 1; p < NbrSampleSizes; ++p)
	{
	  ++CurrentRunNumberParticles[0][p];
	  if (RadiusSquare < RadiusSquareInsideDroplet)
	    ++CurrentRunNumberParticles[1][p];
	}
      }
      for (int p = 0; p < NbrSampleSizes; ++p)
      {
	for (int l = 0; l < 2; ++l)
	{
	  double currentNbrParticles  = ((double) CurrentRunNumberParticles[l][p]);
	  for (int order = 0; order < 4; ++order)
	  {
	    MomentOrderNbrParticles[l][order][p] += currentNbrParticles;
	    currentNbrParticles *= ((double) CurrentRunNumberParticles[l][p]);
	  }
	}
      }
    }
    
  }
  for (int l = 0; l < 2; ++l)
    for (int p = 0; p < NbrSampleSizes; ++p)
      for (int order = 0; order < 4; ++order)
	MomentOrderNbrParticles[l][order][p] /= ((double) (Manager.GetInteger("nbr-samples") / MeasurementPeriod));
    
  char* FileName = new char[512];
  if (EllipseParameter == 0.0)
    sprintf(FileName, "fluctuation_nbr_particles_laughlin_%d_disk_n_%d_MCsamples_%.0e_period_%d_%s.dat", LaughlinInvertFillingFraction, NbrParticles, ((double) Manager.GetInteger("nbr-samples")), ((int) Manager.GetInteger("discretize-disk")), Manager.GetString("file-name"));
  else
    sprintf(FileName, "fluctuation_nbr_particles_laughlin_%d_ellipse_%.6f_n_%d_MCsamples_%.0e_period_%d_%s.dat", LaughlinInvertFillingFraction, EllipseParameter, NbrParticles, ((double) Manager.GetInteger("nbr-samples")), ((int) Manager.GetInteger("discretize-disk")), Manager.GetString("file-name"));
  
  double** Variance = new double* [2];
  for (int l = 0; l < 2; ++l)
    Variance[l] = new double [NbrSampleSizes];
  
  ofstream File;
  File.open(FileName, ios::binary | ios::out);
  File.precision(14);
  File << "# l <N_A>_c  <N_A^2>_c <N_A^3>_c <N_A^4>_c" << endl;
  for (int p = 0; p < NbrSampleSizes; ++p)
  {
    File << p << " ";
    for (int l = 0; l < 2; ++l)
    {
      double SquareAverage = MomentOrderNbrParticles[l][0][p] * MomentOrderNbrParticles[l][0][p];
      double CubeAverage = MomentOrderNbrParticles[l][0][p] * SquareAverage;
      double FourthAverage = MomentOrderNbrParticles[l][0][p] * CubeAverage;
      Variance[l][p] = MomentOrderNbrParticles[l][1][p] - SquareAverage;
      double ThirdCumulant = MomentOrderNbrParticles[l][2][p] - 3 * MomentOrderNbrParticles[l][1][p] * MomentOrderNbrParticles[l][0][p] + 2 * CubeAverage;
      double FourthCumulant = MomentOrderNbrParticles[l][3][p] - 4 * MomentOrderNbrParticles[l][2][p] * MomentOrderNbrParticles[l][0][p] - 3 * MomentOrderNbrParticles[l][1][p] * MomentOrderNbrParticles[l][1][p] + 12 * MomentOrderNbrParticles[l][1][p] * SquareAverage - 6 * FourthAverage;
       File << MomentOrderNbrParticles[l][0][p] << " " << Variance[l][p] << " " <<  ThirdCumulant << " " << FourthCumulant << " ";
    }
    if ((p >= NbrSampleSizes / 2) && (p < NbrSampleSizes - 1))
    {
      double RescaledFluctuations = (Variance[0][p] + Variance[0][NbrSampleSizes - 2 - p]) - (Variance[1][p] + Variance[1][NbrSampleSizes - 2 - p]);
      File << (RescaledFluctuations / 2.0);
    }
    File << endl;
  }
  File.close();
  
  
  delete[] FileName;
  
  delete[] ParticleCoordinateX;
  delete[] ParticleCoordinateY;
  for (int l = 0; l < 2; ++l)
  {
    for (int i = 0; i < 4; ++i)
      delete[] MomentOrderNbrParticles[l][i];
    delete[] CurrentRunNumberParticles[l];
    delete[] MomentOrderNbrParticles[l];
    delete[] Variance[l];
  }
  delete[] MomentOrderNbrParticles;
  delete[] CurrentRunNumberParticles;
  delete[] Variance;
}