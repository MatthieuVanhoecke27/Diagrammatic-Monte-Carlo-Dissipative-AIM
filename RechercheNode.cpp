/**
 * @file RechercheNode.cpp
 * @brief Fonctions de navigation et de recherche dans la clist.
 * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#include "RechercheNode.hpp"

/*-------------------------------------------------------------------------------------------------------------------------------------
    Recherche des nœuds précédents pour Alpha et Time
    Utilisé pour reconnecter la liste lors d'une suppression (Remove).
 ------------------------------------------------------------------------------------------------------------------------------------*/

void clist::Rch_Previous_Alpha_time(node *& startNode, node *& prev_time, node *& prev_alpha, int alpha, int flag, int channel) {
    if (startNode == nullptr) {
        prev_time = nullptr;
        prev_alpha = nullptr;
        return;
    }

    node* readNode = startNode->p_prev; // On commence par l'élément juste avant startNode
    prev_alpha = nullptr;
    prev_time = nullptr;

    bool foundAlpha = false;
    bool foundTime = false;

    // On remonte la liste globale pour trouver les voisins directs
    while (readNode != nullptr && (!foundAlpha || !foundTime)) {
        if (readNode->data.channel == channel) {
            // Recherche du précédent avec le même alpha
            if (!foundAlpha && readNode->data.alpha == alpha) {
                prev_alpha = readNode;
                foundAlpha = true;
            }
            // Recherche du précédent avec le même flag (t ou t-bar)
            if (!foundTime && readNode->data.flag == flag) {
                prev_time = readNode;
                foundTime = true;
            }
        }
        readNode = readNode->p_prev;
    }
}

/*-------------------------------------------------------------------------------------------------------------------------------------
    Recherche du point d'insertion pour un nouveau Kink
    Détermine la position temporelle correcte (previousNode) et les liens secondaires.
 ------------------------------------------------------------------------------------------------------------------------------------*/



void clist::RecherchePreviousNode(kink & KinkAdd, node *& previousNode, node *& previousAlpha, node *& previous_time) {
    previous_time = nullptr;
    previousNode = nullptr;
    
    // On commence la recherche à partir de previousAlpha s'il est fourni, sinon du head
    node* readNode = (previousAlpha != nullptr) ? previousAlpha : head;

    // Si la liste est totalement vide
    if (readNode == nullptr) return;

    // Phase 1 : Si le readNode est déjà après le temps cible, on doit reculer
    if (readNode->data.time > KinkAdd.time) {
        while (readNode != nullptr && readNode->data.time > KinkAdd.time) {
            readNode = readNode->p_prev;
        }
    }

    // Phase 2 : On avance jusqu'à trouver la position juste avant le temps cible
    while (readNode != nullptr) {
        // Est-ce le nœud immédiatement avant dans le temps global ?
        if (readNode->data.time < KinkAdd.time) {
            previousNode = readNode;
        } else {
            // On a dépassé le temps cible
            break;
        }

        // Est-ce le nœud précédent pour le même flag/channel ?
        if (readNode->data.flag == KinkAdd.flag && readNode->data.channel == KinkAdd.channel) {
            previous_time = readNode;
        }
        
        // Sécurité pour ne pas boucler à l'infini
        if (readNode->p_next == nullptr) break;
        readNode = readNode->p_next;
    }
    
    // Note : previousAlpha est généralement passé en argument car il est
    // souvent déjà connu par l'algorithme de proposition de mouvement.
}

/*-------------------------------------------------------------------------------------------------------------------------------------
    Recherche de l'index et du temps précédent
    Essentiel pour calculer le ratio des déterminants (Metropolis-Hastings).
 ------------------------------------------------------------------------------------------------------------------------------------*/

void clist::RecherchePreviousTime_and_index(kink & Kinkremove, int & Index, node *& previous_time) {
    Index = 0;
    previous_time = nullptr;
    
    // On utilise la tête de liste spécifique au flag/channel pour aller plus vite
    node* readNode = time_head[Kinkremove.flag][Kinkremove.channel];
    
    if (readNode == nullptr || readNode->data.time >= Kinkremove.time) {
        return; // L'élément à supprimer est le premier de sa chaîne
    }

    while (readNode != nullptr) {
        // Si le prochain nœud dans la chaîne est celui qu'on cherche ou est après
        if (readNode->fwd_time == nullptr || readNode->fwd_time->data.time >= Kinkremove.time) {
            previous_time = readNode;
            Index++;
            break;
        }
        readNode = readNode->fwd_time;
        Index++;
    }
}
