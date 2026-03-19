/**
 * @file MoveShifting.cpp
 * @brief Implementation of the shift moves in the Diagrammatic Monte Carlo algorithm.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "MoveShifting.hpp"

/*------------------------------------------------------------------------------------------------------------------------------------
                                         Propose Shift, Move Execution & Trace
 ------------------------------------------------------------------------------------------------------------------------------------*/

/**
 * @brief Proposes a new time for a randomly selected operator within the configuration.
 * @param NodeShift Output: pointer to the node selected to be shifted.
 * @param prev_alpha Output: pointer to the previous node in the alpha-linked list.
 * @param New_time Output: the newly sampled time for the operator.
 * @param T_max Maximum simulation time (boundary for sampling).
 */
void clist::propose_shift(node *& NodeShift, node *& prev_alpha, double & New_time, double T_max) {
    short int channel;
    short int alpha;
    
    if (empty_list() == 0) {
        std::cout << "Error Shifting -> List is empty" << std::endl;
    }
    else {
        // --- Determine Alpha selection ---
        if (empty_alpha(0) == 0) {
            alpha = 1;
        }
        else if (empty_alpha(1) == 0) {
            alpha = 0;
        }
        else {
            alpha = rand() % 2;
        }
        
        // --- Determine Channel selection ---
        if (nb_operator[alpha][0] == 0) {
            channel = 1;
        }
        else if (nb_operator[alpha][1] == 0) {
            channel = 0;
        }
        else {
            channel = rand() % 2;
        }
        
        std::random_device rd;
        std::default_random_engine eng_shift(rd());
        
        // --- Sampling Shifting Parameters ---
        int index = rand() % nb_operator[alpha][channel];
        double t_lower = 0.0;
        double t_upper = T_max;
        
        prev_alpha = nullptr;
        NodeShift = alph_head[alpha][channel];

        for (int i = 0; i < index; i = i + 1) {
            if (i == index - 1) {
                prev_alpha = NodeShift;
            }
            NodeShift = NodeShift->fwd_alpha;
        }
        
        // Boundaries for time sampling based on neighboring nodes in alpha-list
        if (NodeShift->fwd_alpha != nullptr) {
            t_upper = (NodeShift->fwd_alpha)->data.time;
        }
        if (prev_alpha != nullptr) {
            t_lower = prev_alpha->data.time;
        }
        
        std::uniform_real_distribution<double> distr_dag(t_lower, t_upper);
        New_time = distr_dag(eng_shift);
    }
}

/**
 * @brief Executes the shifting move by updating the node's position in the linked lists.
 * @param NodeShift Pointer to the node being shifted.
 * @param prev_alpha Pointer to the previous node in the alpha-list.
 * @param New_time The target time for the shift.
 */
void clist::Move_shifting(node *& NodeShift, node *& prev_alpha, double & New_time) {
    if (NodeShift == nullptr) {
        std::cout << "Error Move_shifting -> No node provided for shift" << std::endl;
    }
    else {
        // --- Locate previous_time for removal logic ---
        node * readNode = NodeShift;
        node * previous_time = nullptr;
        kink NewNode(NodeShift->data.alpha, NodeShift->data.channel, NodeShift->data.flag, New_time);
        short int test = 0;

        while (test == 0) {
            if (readNode->p_prev == nullptr) break;
            readNode = readNode->p_prev;
            if ((readNode->data.flag == NodeShift->data.flag) && (readNode->data.channel == NodeShift->data.channel)) {
                previous_time = readNode;
                test = 1;
                break;
            }
        }
        
        // Remove node from current position
        remove_node(prev_alpha, previous_time, NodeShift);
        
        // --- Search for new insertion point (previous_time) ---
        if (previous_time != nullptr && (previous_time->data.time < NewNode.time)) {
            if (previous_time->fwd_time != nullptr) {
                while ((previous_time->fwd_time->data.time) < New_time) {
                    previous_time = previous_time->fwd_time;
                    if (previous_time->fwd_time == nullptr) break;
                }
            }
        }
        else {
            previous_time = time_head[NewNode.flag][NewNode.channel];
            if (previous_time != nullptr) {
                if (previous_time->data.time > NewNode.time) {
                    previous_time = nullptr;
                }
                else {
                    if (previous_time->fwd_time != nullptr) {
                        while ((previous_time->fwd_time->data.time) < New_time) {
                            previous_time = previous_time->fwd_time;
                            if (previous_time->fwd_time == nullptr) break;
                        }
                    }
                }
            }
        }
        
        // --- Search for new insertion point (previous_node) ---
        node * previous_node = nullptr;
        if (previous_time != nullptr && (previous_time->data.time < NewNode.time)) {
            readNode = previous_time;
        }
        else {
            readNode = head;
        }

        if (readNode != nullptr && readNode->data.time > New_time) {
            readNode = nullptr;
        }

        if (readNode != nullptr) {
            if (readNode->p_next != nullptr) {
                while ((readNode->p_next->data.time) < New_time) {
                    readNode = readNode->p_next;
                    if (readNode->p_next == nullptr) break;
                }
            }
        }
        previous_node = readNode;

        // Re-insert node with updated time
        add_node(prev_alpha, previous_time, previous_node, NewNode);
    }
}

/**
 * @brief Calculates the trace ratio components for a shifting move.
 * @param T_max Maximum simulation time.
 * @param NodeShift The node being shifted.
 * @param prev_alpha Previous node in alpha-linked list.
 * @param New_time Target time for the shift.
 * @param l_sig Weight adjustment for the same Hilbert space.
 * @param W_sig Overlap change between H and H_tilde.
 * @param O_sig Overlap change between Sigma and Sigma_bar.
 */
void clist::Trace_Shifting(double T_max, node *& NodeShift, node *& prev_alpha, double & New_time, std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig) {
    
    l_sig[NodeShift->data.alpha] = 0.0;

    // --- Trace ratio for same Hilbert space (l_sigma) ---
    if (NodeShift->data.flag == NodeShift->data.alpha) {
        l_sig[NodeShift->data.alpha] = New_time - NodeShift->data.time;
    }
    else {
        l_sig[NodeShift->data.alpha] = NodeShift->data.time - New_time;
    }
    
    // --- Overlap H and H_tilde (W_sigma) ---
    W_sig[NodeShift->data.channel] = 0.0;
    if (NodeShift->data.flag == NodeShift->data.alpha) {
        if (New_time < NodeShift->data.time) {
            double overlap = 0.0;
            Overlap_HHtilde(overlap, T_max, NodeShift->data.alpha, NodeShift->data.channel, New_time, NodeShift->data.time);
            W_sig[NodeShift->data.channel] = -1.0 * overlap;
        }
        else {
            Overlap_HHtilde(W_sig[NodeShift->data.channel], T_max, NodeShift->data.alpha, NodeShift->data.channel, NodeShift->data.time, New_time);
        }
    }
    else {
        if (New_time < NodeShift->data.time) {
            Overlap_HHtilde(W_sig[NodeShift->data.channel], T_max, NodeShift->data.alpha, NodeShift->data.channel, New_time, NodeShift->data.time);
        }
        else {
            double overlap = 0.0;
            Overlap_HHtilde(overlap, T_max, NodeShift->data.alpha, NodeShift->data.channel, NodeShift->data.time, New_time);
            W_sig[NodeShift->data.channel] = -1.0 * overlap;
        }
    }
    
    // --- Overlap Sigma and Sigma_bar (O_sigma) ---
    O_sig[NodeShift->data.alpha] = 0.0;
    
    if (NodeShift->data.flag == NodeShift->data.alpha) {
        if (New_time < NodeShift->data.time) {
            double overlap = 0.0;
            Overlap_SigSigBar(overlap, T_max, NodeShift->data.alpha, NodeShift->data.channel, New_time, NodeShift->data.time);
            O_sig[NodeShift->data.alpha] = -1.0 * overlap;
        }
        else {
            Overlap_SigSigBar(O_sig[NodeShift->data.alpha], T_max, NodeShift->data.alpha, NodeShift->data.channel, NodeShift->data.time, New_time);
        }
    }
    else {
        if (New_time < NodeShift->data.time) {
            Overlap_SigSigBar(O_sig[NodeShift->data.alpha], T_max, NodeShift->data.alpha, NodeShift->data.channel, New_time, NodeShift->data.time);
        }
        else {
            double overlap = 0.0;
            Overlap_SigSigBar(overlap, T_max, NodeShift->data.alpha, NodeShift->data.channel, NodeShift->data.time, New_time);
            O_sig[NodeShift->data.alpha] = -1.0 * overlap;
        }
    }
}
