/**
 * @file Probability_Move.cpp
 * @brief Implementation of proposal probability calculations and time-ordering signs.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "Probability_Move.hpp"

/*----------------------------------------------------------------------------------------------------
                                      Proposal Probability: ADDING Move
 -----------------------------------------------------------------------------------------------------*/

/**
 * @brief Calculates the probability of proposing an 'adding' move.
 * @details Computes the product of selection probabilities for channel, vertex type, and time sampling.
 */
double clist::Proba_adding(int channel, double T_max, node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter, std::vector<std::vector<int>> & NombreOpe) {
    
    double Probabilite_move = 1.0;
    
    // Probability of selecting the 'Add' move among available moves
    if ((NombreOpe[0][0] + NombreOpe[1][0] + NombreOpe[0][1] + NombreOpe[1][1]) != 0) {
        Probabilite_move = Probabilite_move * (1.0 / 4.0);
    }

    // Probability of channel selection
    Probabilite_move = Probabilite_move * 0.5;
    
    // Probability of choosing between Diagonal and Off-Diagonal vertex
    Probabilite_move = Probabilite_move * 0.5;
    
    if (OperatorEarly.alpha == OperatorLatter.alpha) {
        // --- Diagonal Vertex Case ---
        // Probability of selecting alpha (0 or 1)
        Probabilite_move = Probabilite_move * (1.0 / 2.0);
        
        // Probability of selecting the insertion index
        Probabilite_move = Probabilite_move * (1.0 / (NombreOpe[OperatorEarly.alpha][OperatorEarly.channel] + 1.0));
        
        // --- Time Sampling Probability ---
        if (lower == nullptr) {
            if (upper == nullptr) {
                // List is empty for this specific alpha/channel
                if (OperatorEarly.alpha == OperatorEarly.flag) {
                    // Initial full density matrix: Early is annihilation operator
                    Probabilite_move = Probabilite_move * (1.0 / T_max) * (1.0 / OperatorLatter.time);
                }
                else {
                    // Initial empty density matrix: Early is creation operator
                    Probabilite_move = Probabilite_move * (1.0 / T_max) * (1.0 / (T_max - OperatorEarly.time));
                }
            }
            else {
                // Case: Insertion at the beginning of the list
                Probabilite_move = Probabilite_move * (1.0 / (upper->data.time)) * (1.0 / (OperatorLatter.time));
            }
        }
        else {
            if (upper == nullptr) {
                // Case: Insertion at the end of the list
                Probabilite_move = Probabilite_move * (1.0 / (T_max - lower->data.time)) * (1.0 / (OperatorLatter.time - lower->data.time));
            }
            else {
                // Case: Insertion between two existing nodes
                Probabilite_move = Probabilite_move * (1.0 / (upper->data.time - lower->data.time)) * (1.0 / (OperatorLatter.time - lower->data.time));
            }
        }
    }
    else {
        // --- Off-Diagonal Vertex Case ---
        // Operators already selected, determine time sampling probability
        if (upper == nullptr) {
            if (lower == nullptr) {
                Probabilite_move = Probabilite_move * (1.0 / T_max) * (1.0 / T_max);
            }
            else {
                Probabilite_move = Probabilite_move * (1.0 / (T_max - lower->data.time)) * (1.0 / T_max);
            }
        }
        else {
            if (lower == nullptr) {
                Probabilite_move = Probabilite_move * (1.0 / T_max) * (1.0 / (T_max - upper->data.time));
            }
            else {
                Probabilite_move = Probabilite_move * (1.0 / (T_max - lower->data.time)) * (1.0 / (T_max - upper->data.time));
            }
        }
    }
    return Probabilite_move;
}

/*----------------------------------------------------------------------------------------------------
                                     Proposal Probability: REMOVING Move
 -----------------------------------------------------------------------------------------------------*/

/**
 * @brief Helper function to check if a channel is empty.
 * @return 0 if empty, 1 if not empty.
 */
int EmptyChannel_NbreOpe(int & NbreOpe_Alpha0, int & NbreOpe_Alpha1) {
    return (NbreOpe_Alpha0 == 0 && NbreOpe_Alpha1 == 0) ? 0 : 1;
}

/**
 * @brief Calculates the probability of proposing a 'removing' move.
 */
double clist::Proba_Removing(kink & Ope_removeEarly, kink & Ope_removeLatter, std::vector<std::vector<int>> & NombreOpe) {
    
    double Probabilite_move = 1.0;
    Probabilite_move = Probabilite_move * (1.0 / 4.0);

    // Channel selection probability
    if (EmptyChannel_NbreOpe(NombreOpe[0][0], NombreOpe[1][0]) == 1 && EmptyChannel_NbreOpe(NombreOpe[0][1], NombreOpe[1][1]) == 1) {
        Probabilite_move = Probabilite_move * (1.0 / 2.0);
    }

    // Diagonal vs Off-Diagonal move selection probability
    if (NombreOpe[0][Ope_removeEarly.channel] == 0 || NombreOpe[1][Ope_removeEarly.channel] == 0) {
        Probabilite_move = Probabilite_move * 1.0;
    }
    else if (NombreOpe[0][Ope_removeEarly.channel] == 1 && NombreOpe[1][Ope_removeEarly.channel] == 1) {
        Probabilite_move = Probabilite_move * 1.0;
    }
    else {
        Probabilite_move = Probabilite_move * (1.0 / 2.0);
    }

    // Vertex selection probability (for Diagonal moves)
    if (Ope_removeEarly.alpha == Ope_removeLatter.alpha) {
        if (NombreOpe[1][Ope_removeEarly.channel] > 1 && NombreOpe[0][Ope_removeEarly.channel] > 1) {
            Probabilite_move = Probabilite_move * (1.0 / 2.0);
        }
        // Index selection probability
        Probabilite_move = Probabilite_move * (1.0 / (static_cast<double>(NombreOpe[Ope_removeEarly.alpha][Ope_removeEarly.channel]) - 1.0));
    }
    return Probabilite_move;
}

/*----------------------------------------------------------------------------------------------------
                                     Proposal Probability: SHIFTING Move
 -----------------------------------------------------------------------------------------------------*/

int EmptyAlpha_NbreOpe(int & NbreOpe_Channel0, int & NbreOpe_Channel1) {
    return (NbreOpe_Channel0 == 0 && NbreOpe_Channel1 == 0) ? 0 : 1;
}

/**
 * @brief Calculates the probability of proposing a 'shifting' move.
 */
double clist::Proba_Shifting(node *& NodeShift, node *& prev_alpha, double & New_time, std::vector<std::vector<int>> & NombreOpe) {
    
    double Probabilite_move = 1.0;
    Probabilite_move = Probabilite_move * (1.0 / 3.0);

    // Alpha selection probability
    if (EmptyAlpha_NbreOpe(NombreOpe[0][0], NombreOpe[0][1]) == 1 && EmptyAlpha_NbreOpe(NombreOpe[1][0], NombreOpe[1][1]) == 1) {
        Probabilite_move = Probabilite_move * (1.0 / 2.0);
    }

    // Channel selection probability
    if (NombreOpe[NodeShift->data.alpha][0] != 0 && NombreOpe[NodeShift->data.alpha][1] != 0) {
        Probabilite_move = Probabilite_move * (1.0 / 2.0);
    }
    
    // Operator index selection probability
    Probabilite_move = Probabilite_move * (1.0 / static_cast<double>(NombreOpe[NodeShift->data.alpha][NodeShift->data.channel]));
    
    // Time sampling probability within the interval [prev_alpha, fwd_alpha]
    Probabilite_move = Probabilite_move * (1.0 / (NodeShift->fwd_alpha->data.time - prev_alpha->data.time));
    
    return Probabilite_move;
}

/*----------------------------------------------------------------------------------------------------
                                          Time Ordering Sign
 -----------------------------------------------------------------------------------------------------*/

/**
 * @brief Determines the sign from the time-ordering operator.
 * @details Accounts for the fermionic statistics when swapping spin degrees of freedom.
 */
double clist::Sign_TimeOrdering(node * startNode, node *& EndNode) {
    double sign = 1.0;
    node * ReadNode = startNode;
    int targetChannel = std::abs(startNode->data.channel - 1);

    if (EndNode == nullptr) {
        while (ReadNode != nullptr) {
            if (ReadNode->data.channel == targetChannel) {
                sign *= -1.0;
            }
            ReadNode = ReadNode->p_next;
        }
    }
    else {
        while (ReadNode->data.time != EndNode->data.time) {
            if (ReadNode->data.channel == targetChannel) {
                sign *= -1.0;
            }
            ReadNode = ReadNode->p_next;
        }
    }
    return sign;
}

/*----------------------------------------------------------------------------------------------------
                                      Metropolis Acceptance Ratio
 -----------------------------------------------------------------------------------------------------*/

/**
 * @brief Computes the ratio of proposal probabilities for the Metropolis-Hastings acceptance step.
 * @return The ratio P(backwards) / P(forwards).
 */
double clist::Probability_move(int & move, kink OperatorEarly, kink OperatorLatter, node * upper, node * lower, int channel, double T_max, std::vector<std::vector<int>> & NombreOpe) {
    if (move == 0) {
        // --- Move: Adding a vertex ---
        double numerator = Proba_adding(channel, T_max, upper, lower, OperatorEarly, OperatorLatter, NombreOpe);
        
        // Temporarily update counts to calculate the reverse probability
        ++NombreOpe[OperatorEarly.alpha][OperatorEarly.channel];
        ++NombreOpe[OperatorLatter.alpha][OperatorLatter.channel];
        
        double denominator = Proba_Removing(OperatorEarly, OperatorLatter, NombreOpe);
        
        // Restore counts
        --NombreOpe[OperatorEarly.alpha][OperatorEarly.channel];
        --NombreOpe[OperatorLatter.alpha][OperatorLatter.channel];
        
        return (denominator / numerator);
    }
    else if (move == 1) {
        // --- Move: Removing a vertex ---
        double numerator = Proba_Removing(OperatorEarly, OperatorLatter, NombreOpe);
        
        --NombreOpe[OperatorEarly.alpha][OperatorEarly.channel];
        --NombreOpe[OperatorLatter.alpha][OperatorLatter.channel];
        
        double denominator = Proba_adding(OperatorEarly.channel, T_max, upper, lower, OperatorEarly, OperatorLatter, NombreOpe);
        
        // Restore counts
        ++NombreOpe[OperatorEarly.alpha][OperatorEarly.channel];
        ++NombreOpe[OperatorLatter.alpha][OperatorLatter.channel];
        
        return (denominator / numerator);
    }
    else {
        return 1.0;
    }
}
