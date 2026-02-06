/**
 * @file FonctionPrincipale.hpp
 * @brief Core structures and clist class definition for DMC Dephasing.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#ifndef FonctionPrincipale_hpp
#define FonctionPrincipale_hpp

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <armadillo>

/* ----------------------------------------------------------------------------------------------------
                                      Global Configuration
 ----------------------------------------------------------------------------------------------------*/

// Number of physical channels (e.g., spin up/down)
constexpr int nbr_channel = 2;

/**
 * @brief Initial state density matrix configuration.
 * rho_initial[alpha][channel] = 0 (empty), 1 (full)
 */
const int rho_initial[2][nbr_channel] = {{0, 0}, {0, 0}};

/* ----------------------------------------------------------------------------------------------------
                                      Class: kink (Data Unit)
 ----------------------------------------------------------------------------------------------------*/

class kink {
public:
    int alpha;   ///< Spinor component (0 or 1)
    int channel; ///< Physical channel (e.g., 0 for Spin Up, 1 for Spin Down)
    int flag;    ///< Time type: 0 for t, 1 for t-bar
    double time; ///< Physical time in [0, T_max]

    // Modern Constructor
    kink(int alph = 0, int chan = 0, int flg = 0, double t = 0.0)
        : alpha(alph), channel(chan), flag(flg), time(t) {}

    // Overloaded assignment operator
    kink& operator=(const kink& other) {
        if (this != &other) {
            alpha = other.alpha;
            channel = other.channel;
            flag = other.flag;
            time = other.time;
        }
        return *this;
    }

    void print() const {
        std::cout << "[Kink] Alpha: " << alpha << " | Chan: " << channel
                  << " | Flag: " << flag << " | Time: " << time << std::endl;
    }
};

/* ----------------------------------------------------------------------------------------------------
                                      Class: node (Linked List Element)
 ----------------------------------------------------------------------------------------------------*/



class node {
    friend class clist;
public:
    // Constructors
    node(kink d)
        : data(d), p_next(nullptr), p_prev(nullptr), fwd_alpha(nullptr), fwd_time(nullptr) {}

    node(kink d, node* next, node* prev, node* n_alph, node* n_time)
        : data(d), p_next(next), p_prev(prev), fwd_alpha(n_alph), fwd_time(n_time) {}

    kink data;
    node *p_next;    ///< Linear forward pointer (global time order)
    node *p_prev;    ///< Linear backward pointer (global time order)
    node *fwd_alpha; ///< Pointer to next node with same alpha/channel
    node *fwd_time;  ///< Pointer to next node with same flag/channel
};

/* ----------------------------------------------------------------------------------------------------
                                      Class: clist (Main Controller)
 ----------------------------------------------------------------------------------------------------*/

class clist {
public:
    clist();
    ~clist();

    // Node Management
    void add_beginning(kink & data);
    void add_end(node *& previousNode, node *& prev_alpha, node *& prev_time, kink & data);
    void add_node(node *& prev_alph, node *& prev_time, node *& previous_node, kink data);
    
    void remove_beginning();
    void remove_node(node *& previousAlpha, node *& previous_time, node *& NodeDel);
    void empty();

    // State Queries
    int empty_alpha(int const alpha);
    int empty_channel(int const channel);
    int empty_list();
    int empty_channel_and_alpha(int const alph, int const chan);

    // Diagnostics & Printing
    void print_config();
    void print_config_alpha(int & alpha);
    void print_config_channel(int & channel);
    void print_config_flag(int & flag, int & channel);
    void saveConfig();
    void test_ordonneList();
    void TestHead();
    void Test_previous();
    void test_boucle(int alpha, int channel);

    // Monte Carlo Move Proposals
    void propose_adding(int channel, double T_max, node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter);
    void propose_removing(node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prev_alpha);
    void propose_shift(node *& NodeShift, node *& prev_alpha, double & New_time, double T_max);

    // Acceptance Actions
    void Move_Adding(node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter);
    void Move_Removing(node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prevAlphaEarly, node *& prevAlphaLatter, node *& prevTimeEarly, node *& prevTimeLatter);
    void Move_shifting(node *& NodeShift, node *& prev_alpha, double & New_time);

    // Search Helpers
    void Rch_Previous_Alpha_time(node *& startNode, node *& prev_time, node *& prev_alpha, int alpha, int flag, int channel);
    void RecherchePreviousNode(kink & KinkAdd, node *& previousNode, node *& previousAlpha, node *& previous_time);
    void RecherchePreviousTime_and_index(kink & KinkAdd, int & Index, node *& previous_time);

    // Trace & Physics Calculations
    void Trace_Adding(double T_max, node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter, std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig, int AddRemove);
    void Overlap_HHtilde(double & W_sig, double T_max, int alpha, int channel, double time_lower, double time_upper);
    void Overlap_SigSigBar(double & O_sig, double T_max, int alpha, int channel, double time_lower, double time_upper);
    void Trace_Shifting(double T_max, node *& NodeShift, node *& prev_alpha, double & New_time, std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig);
    void Trace_Removing(double T_max, node *& Ope_removeEarly, node *& Ope_removeLatter, node *& prev_alphaEarly, node *& prev_alphaLatter, std::vector<double> & l_sig, std::vector<double> & W_sig, std::vector<double> & O_sig);

    // Probability & Signs
    double Proba_adding(int channel, double T_max, node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter, std::vector<std::vector<int>> & NombreOpe);
    double Proba_Removing(kink & Ope_removeEarly, kink & Ope_removeLatter, std::vector<std::vector<int>> & NombreOpe);
    double Proba_Shifting(node *& NodeShift, node *& prev_alpha, double & New_time, std::vector<std::vector<int>> & NombreOpe);
    double Probability_move(int & move, kink Ope_removeEarly, kink Ope_removeLatter, node * upper, node * lower, int channel, double T_max, std::vector<std::vector<int>> & NombreOpe);
    double Sign_TimeOrdering(node * startNode, node *& EndNode);

private:
    node * head;                         ///< Global head (minimum time)
    node * tail[2][nbr_channel];         ///< Tail per alpha/channel
    node * alph_head[2][nbr_channel];    ///< Head per alpha/channel
    node * time_head[2][nbr_channel];    ///< Head per flag/channel
    int nb_operator[2][nbr_channel];     ///< Operator counts
    int nb_node;                         ///< Total node count
};

#endif /* FonctionPrincipale_hpp */
