/**
 * @file FonctionPrincipale.cpp
 * @brief Implementation of the clist class for Diagrammatic Monte Carlo data management.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "FonctionPrincipale.hpp"
#include <fstream>
#include <iomanip>

/* ----------------------------------------------------------------------------------------------------
                                      Constructor & Destructor
 ----------------------------------------------------------------------------------------------------*/

clist::clist() {
    head = nullptr;
    nb_node = 0;

    // Initialize multi-dimensional tracking arrays
    for(int a = 0; a < 2; ++a) {
        for(int c = 0; c < 2; ++c) {
            tail[a][c] = nullptr;
            alph_head[a][c] = nullptr;
            time_head[a][c] = nullptr;
            nb_operator[a][c] = 0;
        }
    }
}

clist::~clist() {
    empty();
}

/* ----------------------------------------------------------------------------------------------------
                                      State Checking Methods
 ----------------------------------------------------------------------------------------------------*/

int clist::empty_alpha(int const alpha) {
    return (alph_head[alpha][0] == nullptr && alph_head[alpha][1] == nullptr) ? 0 : 1;
}

int clist::empty_channel(int const chan) {
    return (alph_head[0][chan] == nullptr && alph_head[1][chan] == nullptr) ? 0 : 1;
}

int clist::empty_list() {
    return (nb_node == 0) ? 0 : 1;
}

int clist::empty_channel_and_alpha(int const alph, int const chan) {
    return (alph_head[alph][chan] == nullptr) ? 0 : 1;
}

/* ----------------------------------------------------------------------------------------------------
                                      Configuration Printing
 ----------------------------------------------------------------------------------------------------*/

void clist::print_config() {
    if (nb_node == 0) {
        std::cout << "List is empty." << std::endl;
        return;
    }
    
    node* read_node = head;
    int i = 0;
    while (read_node != nullptr) {
        std::cout << "Node #" << i++ << " | Alpha: " << read_node->data.alpha
                  << " | Chan: " << read_node->data.channel
                  << " | Flag: " << read_node->data.flag
                  << " | Time: " << read_node->data.time << std::endl;
        read_node = read_node->p_next;
    }
}

/* ----------------------------------------------------------------------------------------------------
                                      Node Insertion Logic
 ----------------------------------------------------------------------------------------------------*/



void clist::add_beginning(kink & data) {
    node* new_node = new node(data);
    new_node->p_next = head;

    if (head != nullptr) {
        head->p_prev = new_node;
    }

    // Update global head
    head = new_node;

    // Handle tail update for specific branch
    if (empty_channel_and_alpha(data.alpha, data.channel) == 0) {
        tail[data.alpha][data.channel] = new_node;
    }

    // Update branch-specific head (alpha-channel)
    new_node->fwd_alpha = alph_head[data.alpha][data.channel];
    alph_head[data.alpha][data.channel] = new_node;

    // Update flag-specific head (flag-channel)
    new_node->fwd_time = time_head[data.flag][data.channel];
    time_head[data.flag][data.channel] = new_node;

    ++nb_node;
    ++nb_operator[data.alpha][data.channel];
}

void clist::add_end(node *& previousNode, node *& prev_alpha, node *& prev_time, kink & data) {
    node* new_node = new node(data);
    new_node->p_prev = previousNode;
    tail[data.alpha][data.channel] = new_node;

    if (nb_node == 0) {
        head = new_node;
    } else {
        previousNode->p_next = new_node;

        // Link specific alpha sequence
        if (prev_alpha == nullptr) alph_head[data.alpha][data.channel] = new_node;
        else prev_alpha->fwd_alpha = new_node;

        // Link specific time sequence
        if (prev_time == nullptr) time_head[data.flag][data.channel] = new_node;
        else prev_time->fwd_time = new_node;
    }

    ++nb_node;
    ++nb_operator[data.alpha][data.channel];
}

void clist::add_node(node *& prev_alph, node *& prev_time, node *& previous_node, kink data) {
    if (previous_node == nullptr) {
        add_beginning(data);
    } else if (previous_node->p_next == nullptr) {
        add_end(previous_node, prev_alph, prev_time, data);
    } else {
        node* new_node = new node(data);
        node* next_node = previous_node->p_next;

        // Linear links
        previous_node->p_next = new_node;
        new_node->p_prev = previous_node;
        new_node->p_next = next_node;
        next_node->p_prev = new_node;

        // Alpha-specific links
        if (prev_alph == nullptr) {
            new_node->fwd_alpha = alph_head[data.alpha][data.channel];
            alph_head[data.alpha][data.channel] = new_node;
            if (new_node->fwd_alpha == nullptr) tail[data.alpha][data.channel] = new_node;
        } else {
            new_node->fwd_alpha = prev_alph->fwd_alpha;
            if (prev_alph->fwd_alpha == nullptr) tail[data.alpha][data.channel] = new_node;
            prev_alph->fwd_alpha = new_node;
        }

        // Flag-specific links
        if (prev_time == nullptr) {
            new_node->fwd_time = time_head[data.flag][data.channel];
            time_head[data.flag][data.channel] = new_node;
        } else {
            new_node->fwd_time = prev_time->fwd_time;
            prev_time->fwd_time = new_node;
        }

        ++nb_node;
        ++nb_operator[data.alpha][data.channel];
    }
}

/* ----------------------------------------------------------------------------------------------------
                                      Node Removal Logic
 ----------------------------------------------------------------------------------------------------*/

void clist::remove_beginning() {
    if (nb_node == 0) return;

    node* del = head;
    int a = del->data.alpha;
    int c = del->data.channel;
    int f = del->data.flag;

    head = head->p_next;
    if (head != nullptr) head->p_prev = nullptr;

    // Update secondary heads
    alph_head[a][c] = del->fwd_alpha;
    time_head[f][c] = del->fwd_time;

    delete del;
    --nb_node;
    --nb_operator[a][c];

    if (alph_head[a][c] == nullptr) tail[a][c] = nullptr;
}

void clist::remove_node(node *& previousAlpha, node *& previous_time, node *& NodeDel) {
    int a = NodeDel->data.alpha;
    int c = NodeDel->data.channel;
    int f = NodeDel->data.flag;

    // 1. Update Alpha Chain
    if (previousAlpha == nullptr) alph_head[a][c] = NodeDel->fwd_alpha;
    else previousAlpha->fwd_alpha = NodeDel->fwd_alpha;
    
    if (NodeDel->fwd_alpha == nullptr) tail[a][c] = previousAlpha;

    // 2. Update Time/Flag Chain
    if (previous_time == nullptr) time_head[f][c] = NodeDel->fwd_time;
    else previous_time->fwd_time = NodeDel->fwd_time;

    // 3. Update Linear Chain
    if (NodeDel->p_prev == nullptr && NodeDel->p_next == nullptr) {
        head = nullptr;
    } else {
        if (NodeDel->p_prev == nullptr) {
            head = NodeDel->p_next;
            head->p_prev = nullptr;
        } else {
            NodeDel->p_prev->p_next = NodeDel->p_next;
        }

        if (NodeDel->p_next != nullptr) {
            NodeDel->p_next->p_prev = NodeDel->p_prev;
        }
    }

    --nb_operator[a][c];
    delete NodeDel;
    --nb_node;
}

void clist::empty() {
    while (nb_node > 0) remove_beginning();
}

/* ----------------------------------------------------------------------------------------------------
                                      Diagnostics & I/O
 ----------------------------------------------------------------------------------------------------*/

void clist::test_ordonneList() {
    if (nb_node < 2) return;
    node* res = head;
    while (res->p_next != nullptr) {
        if (res->data.time > res->p_next->data.time) {
            std::cout << "[ERROR] List not ordered at t=" << res->data.time << std::endl;
            return;
        }
        res = res->p_next;
    }
    std::cout << "[OK] Time order verified." << std::endl;
}

void clist::saveConfig() {
    std::string path = "/Users/adm.matthieu.vanhoecke/Desktop/Dia_MonteCarlo/Config_N=100.txt";
    std::ofstream ffile(path);

    if (!ffile.is_open()) {
        std::cerr << "Error: Could not open " << path << " for writing." << std::endl;
        return;
    }

    node* read = head;
    while (read != nullptr) {
        ffile << std::fixed << std::setprecision(8)
              << read->data.time << " "
              << (double)read->data.alpha << " "
              << (double)read->data.channel << " "
              << (double)read->data.flag << "\n";
        read = read->p_next;
    }
    ffile.close();
}
