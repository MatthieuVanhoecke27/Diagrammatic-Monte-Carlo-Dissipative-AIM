/**
 * @file MoveRemoving.hpp
 * @brief Declarations for operator removal moves in Diagrammatic Monte Carlo.
 * @author Matthieu Vanhoecke
 * @date 2023-03-27
 */

#ifndef MoveRemoving_hpp
#define MoveRemoving_hpp

#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <complex>
#include <random>
#include <vector>
#include "FonctionPrincipale.hpp"

/* ------------------------------------------------------------------------------------------------------------------------------------
                                      Move Selection and Logic
 ------------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Declaration of removal-related methods for the clist class.
 * Note: These methods are implemented in MoveRemoving.cpp.
 */



void propose_removing(node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prev_alpha);

void Move_Removing(node *& Ope_removeEarly, node *& Ope_removeLatter,
                  node *& prevAlphaEarly, node *& prevAlphaLatter,
                  node *& prevTimeEarly, node *& prevTimeLatter);

void Trace_Removing(double T_max, node *& Ope_removeEarly, node *& Ope_removeLatter,
                   node *& prev_alphaEarly, node *& prev_alphaLatter,
                   std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig);


#endif /* MoveRemoving_hpp */
