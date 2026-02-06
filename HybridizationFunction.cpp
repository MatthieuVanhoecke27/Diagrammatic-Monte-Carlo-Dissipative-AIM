/**
 * @file HybridizationFunction.cpp
 * @brief Implementation of bath hybridization functions (Delta) for DMC.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "HybridizationFunction.hpp"

using namespace std;

// Constant for complex i to improve readability
const complex<double> I(0.0, 1.0);

/* ----------------------------------------------------------------------------------------------------
                                      Hybridization Components (Delta)
 ----------------------------------------------------------------------------------------------------*/

/**
 * @brief Delta with alpha = 0, beta = 0
 */
complex<double> Delta_11(double W, double Gam, double potential_mu, double t, double t_prime) {
    double dt = t - t_prime;
    
    if (std::abs(dt) < 1e-15) {
        return -I * Gam * (potential_mu + W);
    }
    
    if (t > t_prime) {
        return -Gam * (1.0 - std::exp(-I * W * dt)) / dt;
    } else {
        return Gam * (std::exp(I * W * dt) - 1.0) / dt;
    }
}

/**
 * @brief Delta with alpha = 1, beta = 1
 */
complex<double> Delta_22(double W, double Gam, double potential_mu, double t, double t_prime) {
    double dt = t - t_prime;
    
    if (std::abs(dt) < 1e-15) {
        return -Gam * I * W;
    }
    
    if (t > t_prime) {
        if (potential_mu > W) return 0.0;
        return -Gam * (std::exp(I * W * dt) - 1.0) / dt;
    } else {
        return Gam * (1.0 - std::exp(-I * W * dt)) / dt;
    }
}

/**
 * @brief Delta with alpha = 0, beta = 1
 */
complex<double> Delta_12(double W, double Gam, double potential_mu, double t, double t_prime) {
    double dt = t - t_prime;
    
    if (std::abs(dt) < 1e-15) {
        return -I * Gam * W;
    }
    
    // Expression identical for t > t_prime and t < t_prime
    return -I * Gam * (std::exp(I * W * dt) - 1.0) / dt;
}

/**
 * @brief Delta with alpha = 1, beta = 0
 */
complex<double> Delta_21(double W, double Gam, double potential_mu, double t, double t_prime) {
    double dt = t - t_prime;
    
    if (std::abs(dt) < 1e-15) {
        return I * Gam * W;
    }
    
    // Expression identical for t > t_prime and t < t_prime
    return I * Gam * (1.0 - std::exp(-I * W * dt)) / dt;
}

/* ----------------------------------------------------------------------------------------------------
                                      Main Interface Function
 ----------------------------------------------------------------------------------------------------*/

/**
 * @brief Dispatches the hybridization calculation based on alpha indices.
 * @param ParametreBath Vector containing {W, Gamma, mu}
 */
complex<double> Hybridization_Funct(int alpha, int alpha_prime, double t, double t_prime, vector<double> ParametreBath) {
    // Safety check for parameters
    if (ParametreBath.size() < 3) {
        cerr << "Error: ParametreBath must contain at least 3 elements {W, Gam, mu}" << endl;
        return 0.0;
    }

    double W = ParametreBath[0];
    double Gam = ParametreBath[1];
    double mu = ParametreBath[2];

    if (alpha == 0) {
        return (alpha_prime == 0) ? Delta_11(W, Gam, mu, t, t_prime)
                                 : Delta_12(W, Gam, mu, t, t_prime);
    }
    else if (alpha == 1) {
        return (alpha_prime == 0) ? Delta_21(W, Gam, mu, t, t_prime)
                                 : Delta_22(W, Gam, mu, t, t_prime);
    }
    else {
        cerr << "Error Hybridization Function: Invalid Alpha index " << alpha << endl;
        return 0.0;
    }
}
