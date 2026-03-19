/**
 * @file MoveRemoving.cpp
 * @brief Implementation of vertex removal moves for Diagrammatic Monte Carlo.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "MoveRemoving.hpp"
#include <random>

/* ------------------------------------------------------------------------------------------------------------------------------------
                                    Propose Removing & Execution
 ------------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Selects a vertex (pair of operators) to be proposed for removal.
 * @param Ope_removeEarly Pointer to the early operator of the vertex to remove.
 * @param Ope_removeLatter Pointer to the latter operator of the selected vertex.
 * @param prev_alpha Pointer to the node preceding Ope_removeEarly in the alpha-linked list.
 */
void clist::propose_removing(node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prev_alpha) {
    
    std::random_device rd;
    std::default_random_engine eng(rd());
    
    if (empty_list() == 0) {
        std::cout << "Error: Propose Removing called but no vertex available." << std::endl;
    }
    else {
        Ope_removeEarly = nullptr;
        Ope_removeLatter = nullptr;
        
        int channel;
        
        // --- Determine Channel Selection ---
        if (nb_operator[0][0] == 0 && nb_operator[1][0] == 0) {
            channel = 1;
        }
        else if (nb_operator[0][1] == 0 && nb_operator[1][1] == 0) {
            channel = 0;
        }
        else {
            channel = rand() % 2;
        }
        
        // --- Determine Diagonal vs Off-Diagonal move type ---
        int Rand_Diag_OffDiag;
        if (nb_operator[0][channel] == 0 || nb_operator[1][channel] == 0) {
            Rand_Diag_OffDiag = 0;
        }
        else if (nb_operator[0][channel] < 2 && nb_operator[1][channel] < 2) {
            Rand_Diag_OffDiag = 1;
        }
        else {
            Rand_Diag_OffDiag = rand() % 2;
        }
        
        if (Rand_Diag_OffDiag == 0) {
            int alpha;
            
            // --- Determine Alpha label within selected Channel ---
            if (nb_operator[0][channel] == 0) {
                alpha = 1;
            }
            else if (nb_operator[1][channel] == 0) {
                alpha = 0;
            }
            else {
                if (nb_operator[0][channel] == 1) {
                    alpha = 1;
                }
                else if (nb_operator[1][channel] == 1) {
                    alpha = 0;
                }
                else {
                    alpha = rand() % 2;
                }
            }
            
            // --- Select Ope_removeEarly and Ope_removeLatter for Diagonal vertex ---
            int index = rand() % (nb_operator[alpha][channel] - 1);
            
            Ope_removeEarly = alph_head[alpha][channel];
            prev_alpha = nullptr;
            
            for (int i = 0; i < index; i = i + 1) {
                if (i == index - 1) {
                    prev_alpha = Ope_removeEarly;
                }
                Ope_removeEarly = Ope_removeEarly->fwd_alpha;
            }
            Ope_removeLatter = Ope_removeEarly->fwd_alpha;
        }
        else {
            // --- Select operators for Off-Diagonal vertex removal ---
            if (tail[0][channel]->data.time > tail[1][channel]->data.time) {
                Ope_removeEarly = tail[1][channel];
                Ope_removeLatter = tail[0][channel];
            }
            else {
                Ope_removeEarly = tail[0][channel];
                Ope_removeLatter = tail[1][channel];
            }
        }
    }
}

/**
 * @brief Unlinks selected operators from the time and alpha-ordered lists.
 * @param Ope_removeEarly The early operator to be unlinked.
 * @param Ope_removeLatter The latter operator to be unlinked.
 * @param prevAlphaEarly Previous node for early operator (alpha-linked).
 * @param prevAlphaLatter Previous node for latter operator (alpha-linked).
 * @param prevTimeEarly Previous node for early operator (time-linked).
 * @param prevTimeLatter Previous node for latter operator (time-linked).
 */
void clist::Move_Removing(node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prevAlphaEarly, node *& prevAlphaLatter, node *& prevTimeEarly, node *& prevTimeLatter) {

    if (Ope_removeEarly->data.alpha == Ope_removeLatter->data.alpha) {
        // --- Process Diagonal Vertex Removal ---
        // First remove the latter node, then the early node
        remove_node(Ope_removeEarly, prevTimeLatter, Ope_removeLatter);
        remove_node(prevAlphaEarly, prevTimeEarly, Ope_removeEarly);
    }
    else {
        // --- Process Off-Diagonal Vertex Removal ---
        // Remove both nodes using their respective alpha and time predecessors
        remove_node(prevAlphaEarly, prevTimeEarly, Ope_removeEarly);
        remove_node(prevAlphaLatter, prevTimeLatter, Ope_removeLatter);
    }
}

/*--------------------------------------------------------------------------------------------------------------------------------
                                                    Trace Evaluation
 ---------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Computes the trace ratio Tr[C_new] / Tr[C_old] for the removal move.
 * @param T_max Maximum simulation time.
 * @param Ope_removeEarly Node of the early operator.
 * @param Ope_removeLatter Node of the latter operator.
 * @param prev_alphaEarly Predecessor of early operator in alpha order.
 * @param prev_alphaLatter Predecessor of latter operator in alpha order.
 * @param l_sig Vector for trace calculation results.
 * @param W_sig Vector for weight calculation results.
 * @param O_sig Vector for operator calculation results.
 */
void clist::Trace_Removing(double T_max, node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prev_alphaEarly, node *& prev_alphaLatter, std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig) {
    
    int AddRemove = 1;
    
    if (Ope_removeEarly->data.alpha == Ope_removeLatter->data.alpha) {
        // --- Diagonal Vertex Trace Calculation ---
        Trace_Adding(T_max, Ope_removeLatter->fwd_alpha, prev_alphaEarly, Ope_removeEarly->data, Ope_removeLatter->data, l_sig, W_sig, O_sig, AddRemove);
        
        // Invert signs for removal logic
        l_sig[0] = -1.0 * l_sig[0];
        l_sig[1] = -1.0 * l_sig[1];
        W_sig[0] = -1.0 * W_sig[0];
        W_sig[1] = -1.0 * W_sig[1];
        O_sig[0] = -1.0 * O_sig[0];
        O_sig[1] = -1.0 * O_sig[1];
    }
    else {
        // --- Off-Diagonal Vertex Trace Calculation ---
        prev_alphaEarly = nullptr;
        prev_alphaLatter = nullptr;
        
        // Search for the predecessor of Ope_removeEarly
        node * ReadNode = Ope_removeEarly;
        if (ReadNode->p_prev != nullptr) {
            while (prev_alphaEarly == nullptr) {
                ReadNode = ReadNode->p_prev;
                if ((ReadNode->data.alpha == Ope_removeEarly->data.alpha) && (ReadNode->data.channel == Ope_removeEarly->data.channel)) {
                    prev_alphaEarly = ReadNode;
                    break;
                }
                if (ReadNode->p_prev == nullptr) {
                    prev_alphaEarly = nullptr;
                    break;
                }
            }
        }
        
        // Search for the predecessor of Ope_removeLatter
        ReadNode = Ope_removeLatter;
        if (ReadNode->p_prev != nullptr) {
            while (prev_alphaLatter == nullptr) {
                ReadNode = ReadNode->p_prev;
                if ((ReadNode->data.alpha == Ope_removeLatter->data.alpha) && (ReadNode->data.channel == Ope_removeLatter->data.channel)) {
                    prev_alphaLatter = ReadNode;
                    break;
                }
                if (ReadNode->p_prev == nullptr) {
                    prev_alphaLatter = nullptr;
                    break;
                }
            }
        }
        
        // Final trace calculation using the located predecessors
        Trace_Adding(T_max, prev_alphaLatter, prev_alphaEarly, Ope_removeEarly->data, Ope_removeLatter->data, l_sig, W_sig, O_sig, AddRemove);
        
        // Invert signs for removal logic
        l_sig[0] = -1.0 * l_sig[0];
        l_sig[1] = -1.0 * l_sig[1];
        W_sig[0] = -1.0 * W_sig[0];
        W_sig[1] = -1.0 * W_sig[1];
        O_sig[0] = -1.0 * O_sig[0];
        O_sig[1] = -1.0 * O_sig[1];
    }
}
