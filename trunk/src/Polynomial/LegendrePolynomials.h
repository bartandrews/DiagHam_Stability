////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                       adopted from the BDMC library for                    //
//                    Bold Diagrammatic Monte Carlo Project                   //
//                                                                            //
//                    Copyright (C) 2013-2021 Gunnar Möller                   //
//                                                                            //
//                       Class author: Gunnar Möller                          //
//                                                                            //
//            Class implementing the Legendre Polynomials by recursion        //
//                                                                            //
//                         File ported : 29/04/2022                           //
//                                                                            //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef LEGENDREPOLYNOMIALS_H
#define LEGENDREPOLYNOMIALS_H

#include "config.h"

#include <cmath>

#include "Vector/RealVector.h"

class LegendrePolynomials
{

public:
    // precision below which abscissae are treated as equal on consecutive calls
    static double Precision;

    // constructor
    LegendrePolynomials(int MaxOrder);

    //
    double GetPolynomial(int n, double x);

    // Update values for shifted Legendre polynomials and return reference to internal storage
    // @param x argument in [-1,1]
    // @param[out] output array of output elements
    // @param requiredOrder order up to which it should be evaluated
    RealVector const & GetPolynomials(double x);

    // Calculate values for shifted Legendre polynomials up to the required order and return to a supplied array
    // @param x argument in [-1,1]
    // @param[out] output array of output elements
    // @param requiredOrder order up to which it should be evaluated
    void GetPolynomials(double x, double *output, int requiredOrder=-1);

protected:
    // highest degree of Legendre polynomials to evaluate
    int NbrLegendrePolynomials;
    // storage of current values
    RealVector Values;
    // last argument used to evaluate the polynomials
    double LastArgument;
    // last order up to which evaluation is complete
    int LastOrder;


private:
    // Update values for Legendre polynomials up to the maximum order
    // @param x argument in [-1,1]
    // @param requiredOrder order up to which it should be evaluated
    void UpdateAllPolynomials(double x) { return UpdatePolynomials(x, NbrLegendrePolynomials-1); }

    // Update values for Legendre polynomials up to the given order
    // @param x argument in [-1,1]
    // @param requiredOrder order up to which it should be evaluated
    void UpdatePolynomials(double x, int requiredOrder);

    /// Update values for Legendre polynomials up to the given order, if recalculation required
    /// @param x argument in [-1,1]
    /// @param requiredOrder order up to which it should be evaluated
    void ConditionalUpdatePolynomials(double x, int requiredOrder);

    /// check
    bool IsStored(int n, double x)
    {
      return (this->LastOrder>=n && (std::fabs(x-this->LastArgument) < Precision));
    }

};

#endif // LEGENDREPOLYNOMIALS_H
