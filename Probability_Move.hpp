/**
 * @file Probability_Move.hpp
 * @brief Header file for proposal probabilities and Metropolis-Hastings ratios.
 * @author Matthieu Vanhoecke
 * @date 2023-03-27
 */

#ifndef Probability_Move_hpp
#define Probability_Move_hpp

#include <iostream>
#include <complex>
#include <random>
#include <vector>
#include <cmath>
#include <armadillo>
#include "FonctionPrincipale.hpp"

/* ------------------------------------------------------------------------------------------------------------------------------------
                                      Helper Function Declarations
 ------------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Checks if both alpha channels (0 and 1) are empty for a given configuration.
 * @return 0 if empty, 1 if at least one operator exists.
 */
int EmptyChannel_NbreOpe(int & NbreOpe_Alpha0, int & NbreOpe_Alpha1);

/**
 * @brief Checks if a specific alpha is empty across both physical channels (0 and 1).
 * @return 0 if empty, 1 if not.
 */
int EmptyAlpha_NbreOpe(int & NbreOpe_Channel0, int & NbreOpe_Channel1);

/* ------------------------------------------------------------------------------------------------------------------------------------
                                      Class clist: Probability Methods
 ------------------------------------------------------------------------------------------------------------------------------------*/


double Proba_adding(int channel, double T_max, node *& upper, node *& lower,
                   kink & OperatorEarly, kink & OperatorLatter,
                   std::vector<std::vector<int>> & NombreOpe);

double Proba_Removing(kink & Ope_removeEarly, kink & Ope_removeLatter,
                     std::vector<std::vector<int>> & NombreOpe);

double Proba_Shifting(node *& NodeShift, node *& prev_alpha, double & New_time,
                     std::vector<std::vector<int>> & NombreOpe);

double Sign_TimeOrdering(node * startNode, node *& EndNode);

double Probability_move(int & move, kink OperatorEarly, kink OperatorLatter,
                       node * upper, node * lower, int channel, double T_max,
                       std::vector<std::vector<int>> & NombreOpe);


#endif /* Probability_Move_hpp */
