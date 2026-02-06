/**
 * @file HybridizationFunction.hpp
 * @brief Header for bath hybridization (Delta) calculations.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#ifndef HybridizationFunction_hpp
#define HybridizationFunction_hpp

#include <complex>
#include <vector>
#include <iostream>

/* ----------------------------------------------------------------------------------------------------
                                      Component Functions
 ----------------------------------------------------------------------------------------------------*/

/**
 * Calculations for specific spinor components (alpha, beta).
 * W: Bandwidth/Frequency, Gam: Coupling strength, mu: Chemical potential.
 */

std::complex<double> Delta_11(double W, double Gam, double potential_mu, double t, double t_prime);
std::complex<double> Delta_12(double W, double Gam, double potential_mu, double t, double t_prime);
std::complex<double> Delta_21(double W, double Gam, double potential_mu, double t, double t_prime);
std::complex<double> Delta_22(double W, double Gam, double potential_mu, double t, double t_prime);

/* ----------------------------------------------------------------------------------------------------
                                      Dispatch Interface
 ----------------------------------------------------------------------------------------------------*/

/**
 * @brief Main interface to compute the hybridization function Delta_{alpha, alpha'}(t - t').
 * @param alpha Spinor index of the first operator.
 * @param alpha_prime Spinor index of the second operator.
 * @param ParametreBath Vector containing bath parameters: {W, Gamma, mu}.
 */
std::complex<double> Hybridization_Funct(int const& alpha, int const& alpha_prime,
                                         double t, double t_prime,
                                         std::vector<double> ParametreBath);

#endif /* HybridizationFunction_hpp */
