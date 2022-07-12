////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class for calculation of Pseudopotential Coefficients  //
//                                                                            //
//                        last modification : 19/11/2007                      //
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


#ifndef PSEUDOPOTENTIALS_H
#define PSEUDOPOTENTIALS_H


#include "config.h"
#include "AbstractZDensityProfile.h"
#include "TriangularWellEigenfunction.h"

// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials
double* EvaluatePseudopotentials(int nbrFlux, int landauLevel, double layerSeparation=0.0, bool quiet=false);

// evalute pseudopotentials for dipolar interaction in the lowest Landau level Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials
double* EvaluateDipolarPseudopotentials(int nbrFlux, bool quiet = false);

// evalute pseudopotentials for coulomb interaction in a given Landau level with a given density profile
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// zDensity = density distribution of the layer 
// points = number of points where exact pseudopotentials are calculated
// multiplier = number of integration intervals used per point of discretization
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// zDensity2 = (optional) density distribution of layer 2, if absent, taken to be equal to 1st profile
// epsRel = tolerance given to integration routine
// return value = array that conatins the pseudopotentials
double* EvaluateFiniteWidthPseudoPotentials(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity,
					   int points=200, double multiplier=5.0, double layerSeparation=0.0,
					   AbstractZDensityProfile *zDensity2=0, double epsRel=1e-8);


// evalute pseudopotentials for coulomb interaction in a given Landau level with a given density profile,
// but without tabulating and interpolating the pseudopotentials
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// zDensity = density distribution of the layer 
// points = number of points where exact pseudopotentials are calculated
// multiplier = number of integration intervals used per point of discretization
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// zDensity2 = (optional) density distribution of layer 2, if absent, taken to be equal to 1st profile
// epsRel = tolerance given to integration routine
// return value = array that contains the pseudopotentials

double* EvaluateFiniteWidthPseudoPotentialsNoInterpolation(int nbrFlux, int landauLevel, AbstractZDensityProfile *zDensity, int points, double multiplier, double layerSeparation, AbstractZDensityProfile *zDensity2, double epsRel=1e-8);


// evalute one body potentials for two impurities located at the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPolePotential = potential of the impurity located at the north pole
// southPolePotential = potential of the impurity located at the south pole
// return value = array that conatins the pseudopotentials
double* EvaluateOneBodyPotentials(int nbrFlux, int landauLevel, double northPolePotential, double southPolePotential);


// evalute one body potentials for two charges +e located at given distance above the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPoleDistance = distance of charge located above the north pole
// southPoleDistance = distance of charge located above the south pole
// return value = array that conatins the pseudopotentials
double* EvaluateOneBodyCoulombPotentials(int nbrFlux, int landauLevel, double northPoleCharge, double northPoleDistance, double southPoleCharge, double southPoleDistance);


// evaluate pseudopotentials coefficients of the monomials r^n in units of 1/R
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle of the corresponding LLL problem)
// exponentN = exponent of the monomial
// onlyOdd = boolean indidicating whether it's sufficient to reproduce only the odd pseudopotentials V_(2m+1)
// bool verbose = print pseudopotentials on screen
// return value = array that conatins the coefficients V_m(r^n)
// where m runs over 0,...,nbrFlux, or if option onlyOdd given, from 0 to nbrFlux/2 with entries V_(2m+1)(r^n)
double* GetMonomialPseudopotentials(int nbrFlux, int exponentN, bool onlyOdd, bool verbose = false);
  
// evaluate graphene bilayer pseudopotentials
//
// nbrflux = number of flux quanta for the lowest Landau level (i.e. twice the maximum momentum for a single particle
// l1, l2, l3, l4 = Landau level indices (0 for lowest LL, 1 for first excited LL)
double* EvaluateGrapheneBilayerPseudopotentials(int nbrFlux, int& nbrPseudopotentials, int l1, int l2, int l3, int l4, bool verbose = true);

// ZP: Not very fast but it works... CAUTION: Pseudopotentials are not normalized in any particularly nice way, but they will produce a zero energy GS at nu=2/5 with cyclotron energy=0
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (2j+1) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for delta interaction with two Landau levels on Sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials
double** Evaluate2LLSphereDeltaPseudopotentials(int nbrFlux, bool quiet);

// ZP: Not very fast but it works... CAUTION: Pseudopotentials are not normalized in any particularly nice way, but they will produce a zero energy GS at nu=2/5 with cyclotron energy=0
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (-j(j+1)/Q) * (2j+1) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for delta interaction with two Landau levels on Sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that conatins the pseudopotentials
double** Evaluate2LLSphereLaplacianDeltaPseudopotentials(int nbrFlux, bool quiet);

// ZP: Not very fast but it works... 
//
// Uses the formula: V_L = sum_{j=max(|l1-l3|,|l2-l4|)}^{min(l1+l3,l2+l4)} sum_{m1=max(-l1,-l2)}^{min(l1,l2)} sum_{m2=max(-l3,-l4)}^{min(l3,l4)} 
//                           (1/sqrt(Q)) * (-1)^(m2-m1+l1+l2+l3+l4)
//                         <l1,m1; l2,-m1| L,0> <l3,m2; l4,-m2| L,0> <l3,Q; j,0| l1,Q> <l4,Q; j,0| l2,Q> <l1,-m1; j,-(m2-m1)| l3,-m2> <l2,m1; j,(m2-m1)| l4,m2>
//
// evalute pseudopotentials for Coulomb interaction with two Landau levels on sphere
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that contains the pseudopotentials
double** Evaluate2LLSphereCoulombPseudopotentials(int nbrFlux, bool quiet);

// evalute pseudopotentials for Coulomb interaction, with triangular quantum well confinement, on Sphere 
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum on LLL)
// width = width of the well, in units of the magnetic length
// slope = slope of the triangular well, in energy unit per length unit
// quiet = indicate whether Coulomb Pseudopotentials should be printed on screen
// return value = array that contains the pseudopotentials
double* EvaluateTriangularWellPseudopotentials(int nbrFlux, int landauLevel, double bias, double slope, int* pseudospins, int nbrPointsInteg, bool quiet);


#endif //PSEUDOPOTENTIALS_H
