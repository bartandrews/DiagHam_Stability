#include "TriangularWellEigenfunction.h"

#ifndef PREC_GSL
#define PREC_GSL GSL_PREC_DOUBLE
#endif

// Uncomment this two lines to ignore GSL errors
//#include <gsl/gsl_errno.h>
//gsl_set_error_handler_off ();

using namespace std;

// Constructor
//
// width = thickness of the quantum well
// bias = bias of the triangular quantum well in energy per length
//
TriangularWellEigenfunction::TriangularWellEigenfunction(double width, double bias)
{
  // Initial parameters
  this->w=width;
  this->V=bias;
  
  // Useful constants
  this->c1=(this->V)/(this->w);  // (V/w)
  this->c2=1./pow((this->c1),2/3.);  // (V/w)^(-2/3)

  // Compute eigenenergies
  this->NewtonMaxStepsNbr=10;
  this->NewtonEpsRelative=1e-9;
  this->E[0]=Eigenenergy(0);
  this->E[1]=Eigenenergy(1);

  // Compute norms
  this->integration_absolute_error=0;
  this->integration_relative_error=1e-9;
  this->norm[0]=this->Norm(0);
  this->norm[1]=this->Norm(1);

  // Useful constants
  this->c3[0]=this->c2*(this->c1*this->w-this->E[0]);
  this->c3[1]=this->c2*(this->c1*this->w-this->E[1]);
}

//Destructor
//
TriangularWellEigenfunction::~TriangularWellEigenfunction()
{
}

// Returns the energy of the nth eigenstate of the triangular quantum well
//
// n = number of the eigenstate, n=0 (ground state) or n=1 (first excited state)
//
// return value = nth eigenenergy 
//
double TriangularWellEigenfunction::GetEnergy(int n)
{
  return this->E[n];
}

// General solution of the triangular-well potential Schroedinger equation with energy E
//
// x = position of the confinement potential where the function is evaluated
// E = energy 
//
// return value = value of the function  at x
//
double TriangularWellEigenfunction::Psi(double x,double E)
{
#ifdef __GSL__
  double TmpConst=this->c2*(this->c1*this->w-E);
  double y=this->c2*(this->c1*x-E); 
  return (gsl_sf_airy_Bi(TmpConst, PREC_GSL)*gsl_sf_airy_Ai(y, PREC_GSL)-gsl_sf_airy_Ai(TmpConst, PREC_GSL)*gsl_sf_airy_Bi(y, PREC_GSL));
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

// Derivative with respect to E of Psi (defined above) at x=0
//
// x = position of the confinement potential where the derivative is evaluated
// E = energy 
//
// return value = value of the derivative at E
//
double TriangularWellEigenfunction::PsiPrime(double E)
{
#ifdef __GSL__
  double u=(this->V-E)*this->c2, v=-E*this->c2;
  
  return this->c2*(gsl_sf_airy_Ai_deriv(u, PREC_GSL)*gsl_sf_airy_Bi(v, PREC_GSL)-gsl_sf_airy_Ai_deriv(v, PREC_GSL)*gsl_sf_airy_Bi(u, PREC_GSL)+gsl_sf_airy_Ai(u, PREC_GSL)*gsl_sf_airy_Bi_deriv(v, PREC_GSL)-gsl_sf_airy_Ai(v, PREC_GSL)*gsl_sf_airy_Bi_deriv(u, PREC_GSL));
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

// nth eigenfunction of the triangular well confinement potential
//
// x = position of the confinement potential where the function is evaluated
// n = number of the eigenstate (can be either 0 for the ground state or 1)
//
// return value = value of the eigenfunction  at x
//
double TriangularWellEigenfunction::Eigenfunction(double x, int n)
{ 
#ifdef __GSL__
  double TmpConst=this->c3[n];
  double y=this->c2*( this->c1 * x - this->E[n]);
  double N=(this->norm[n]);

  return N*(gsl_sf_airy_Bi(TmpConst, PREC_GSL)*gsl_sf_airy_Ai(y, PREC_GSL)-gsl_sf_airy_Ai(TmpConst, PREC_GSL)*gsl_sf_airy_Bi(y, PREC_GSL));
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

// Define parameters given to a GSL_function, the first one is a pointer which has to be "cast" in order to get a pointer to an object
struct GSL_function_parameters { void * AutoPointer ; double Energy ; } ;

// Define a static GSL_function to integrate the square of the eigenfunction Psi, and then get its norm
//
// x = position
// p = pointer, which adress points toward a structure
//
// return value = Psi^2 in a GSL_function form
//
double TriangularWellEigenfunction::GSL_f(double x, void * p) 
{
#ifdef __GSL__
  // parameters, i.e. pointer to the class and energy are given in the pointer p
  // The "cast" transforms p back to a pointer to a structure
  struct GSL_function_parameters * params = (struct GSL_function_parameters *) p ;
  // Another cast is used to 
  TriangularWellEigenfunction* Well = (TriangularWellEigenfunction*) (params->AutoPointer);
  // Energy is given as the second parameter
  double E = params->Energy;
  double y=Well->Psi(x,E);
  return y*y;
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

// Norm of the nth eigenfunction, norm=1/sqrt(\int_0^w Psi^2(x,E) dx)
//
// n = number of the eigenstate (can be either 0 for the ground state or 1)
//
// return value = norm of the nth eigenstate
//
double TriangularWellEigenfunction::Norm(int n)
{
#ifdef __GSL__
  double E=this->E[n];
  double result,error;
  size_t nbr_eval;
  
  struct GSL_function_parameters params = { this, E };
  gsl_function GSL_PsiSquare;
  GSL_PsiSquare.function = &GSL_f;
  GSL_PsiSquare.params = &params;
  
  gsl_integration_qng (&GSL_PsiSquare, 0, w, this->integration_absolute_error, this->integration_relative_error, &result, &error, &nbr_eval); 

  if ( n == 1) return -1./sqrt(result);
  else return 1./sqrt(result);
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

// Energy of the nth level, given by a WKB-type approximation, used for root finding of psi(x) as the initial value, which gives accurates energies
//
// n = number of the eigenstate (can be either 0 for the ground state or 1)
//
// return value = wkb approximation of the nth eigenenergy
//
double TriangularWellEigenfunction::wkb(int n)
{
  double dummy1=3.*M_PI/2.*(n+11./12.)/(this->w);
  double cut=2.*dummy1*dummy1;
  double dummy2=(n+1.)*M_PI/(this->w);
  if(V < cut) return (this->V)/2.+dummy2*dummy2;
  else return pow(3.*M_PI*(this->V)/(this->w)/2.*(n+0.75),2./3.);
}

double TriangularWellEigenfunction::GSL_Newton_f(double E, void* p)
{
#ifdef __GSL__
   struct GSL_function_parameters * params = (struct GSL_function_parameters *) p ;
 
   TriangularWellEigenfunction* Well = (TriangularWellEigenfunction*) (params->AutoPointer);
 
   double TmpConst=Well->c2*(Well->c1*Well->w-E);
   double y=Well->c2*(-E); 
   
   return (gsl_sf_airy_Bi(TmpConst, PREC_GSL)*gsl_sf_airy_Ai(y, PREC_GSL)-gsl_sf_airy_Ai(TmpConst, PREC_GSL)*gsl_sf_airy_Bi(y, PREC_GSL));
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

double TriangularWellEigenfunction::GSL_Newton_df(double E, void* p)
{
#ifdef __GSL__
   struct GSL_function_parameters * params = (struct GSL_function_parameters *) p ;
 
 TriangularWellEigenfunction* Well = (TriangularWellEigenfunction*) (params->AutoPointer);

   double u=(Well->V-E)*Well->c2, v=-E*Well->c2;
  
   return Well->c2*(gsl_sf_airy_Ai_deriv(u, PREC_GSL)*gsl_sf_airy_Bi(v, PREC_GSL)-gsl_sf_airy_Ai_deriv(v, PREC_GSL)*gsl_sf_airy_Bi(u, PREC_GSL)+gsl_sf_airy_Ai(u, PREC_GSL)*gsl_sf_airy_Bi_deriv(v, PREC_GSL)-gsl_sf_airy_Ai(v, PREC_GSL)*gsl_sf_airy_Bi_deriv(u, PREC_GSL));
#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}
 
void TriangularWellEigenfunction::GSL_Newton_fdf(double E, void* p, double *g, double *dg)
{
#ifdef __GSL__
   struct GSL_function_parameters * params = (struct GSL_function_parameters *) p ;
 
   TriangularWellEigenfunction* Well = (TriangularWellEigenfunction*) (params->AutoPointer);

   double TmpConst=Well->c2*(Well->c1*Well->w-E);
   double y=-E*(Well->c2); 
   double u=(Well->V-E)*Well->c2, v=-E*Well->c2;
  
   *g = (gsl_sf_airy_Bi(TmpConst, PREC_GSL)*gsl_sf_airy_Ai(y, PREC_GSL)-gsl_sf_airy_Ai(TmpConst, PREC_GSL)*gsl_sf_airy_Bi(y, PREC_GSL));
   *dg = Well->c2*(gsl_sf_airy_Ai_deriv(u, PREC_GSL)*gsl_sf_airy_Bi(v, PREC_GSL)-gsl_sf_airy_Ai_deriv(v, PREC_GSL)*gsl_sf_airy_Bi(u, PREC_GSL)+gsl_sf_airy_Ai(u, PREC_GSL)*gsl_sf_airy_Bi_deriv(v, PREC_GSL)-gsl_sf_airy_Ai(v, PREC_GSL)*gsl_sf_airy_Bi_deriv(u, PREC_GSL));

#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return;
#endif
}



double TriangularWellEigenfunction::Eigenenergy(int n)
{
#ifdef __GSL__
  double root=wkb(n),oldroot,delta;
  int i = 0, status;

  struct GSL_function_parameters params = {this};
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  gsl_function_fdf FDF;
 
   FDF.f = &GSL_Newton_f;
   FDF.df = &GSL_Newton_df;
   FDF.fdf = &GSL_Newton_fdf;
   FDF.params = &params;
 
   T = gsl_root_fdfsolver_newton;
   s = gsl_root_fdfsolver_alloc (T);
   gsl_root_fdfsolver_set (s, &FDF, root);
 
   do
     {
       i++;
       status = gsl_root_fdfsolver_iterate (s);
       oldroot = root;
       root = gsl_root_fdfsolver_root (s);
       delta=fabs((root-oldroot)/root);
       //status = gsl_root_test_delta (root, oldroot, 0, NewtonEpsRelative);
       //cout << status << endl;
       //cout << "Step " << i << ", root = " << root << ", delta = " << delta << endl; // testing purpose
     }
   while ( (delta > this->NewtonEpsRelative) && (i < this->NewtonMaxStepsNbr));
 
   gsl_root_fdfsolver_free (s);
  
  // Newton root-finding method written manually
  //for ( int i=0 ; i < this->NewtonStepsNbr ; i++ ) 
  //{
      //cout << " Step " << i << ", E(" << n << ") = " << root << endl;
  //    root -= Psi(0,root)/PsiPrime(root);
  //}

  return root;

#else
   cout << "TriangularWell subroutines require linking to the Gnu Scientific Library!"<<endl;
   return 0;
#endif
}

