/**
 * @file MoveShifting.hpp
 * @brief Declarations for operator shifting moves in the DMC algorithm.
 * @author Matthieu Vanhoecke
 * @date 2023-03-27
 */

#ifndef MoveShifting_hpp
#define MoveShifting_hpp

#include <stdio.h>
#include <iostream>
#include <cstdio>
#include <complex>
#include <random>
#include <vector>

#include "FonctionPrincipale.hpp"

/* ------------------------------------------------------------------------------------------------------------------------------------
                                      Shift Move Declarations
 ------------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Method declarations for the clist class related to time-shifting operators.
 * These functions handle the selection, physical relocation in the list, and weight evaluation.
 */



void propose_shift(node *& NodeShift, node *& prev_alpha, double & New_time, double T_max);

void Move_shifting(node *& NodeShift, node *& prev_alpha, double & New_time);

void Trace_Shifting(double T_max, node *& NodeShift, node *& prev_alpha, double & New_time,
                   std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig);

#endif /* MoveShifting_hpp */
