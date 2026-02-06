//#include "MoveAdding.hpp"

/**
 * @brief Proposes the addition of a new vertex (pair of kinks) in the diagram.
 * * This function handles the stochastic generation of both diagonal and off-diagonal
 * vertices. It identifies the temporal intervals and selects appropriate time
 * coordinates for the operator pair based on the current configuration and
 * boundary conditions (initial density matrix).
 */
void clist::propose_adding(int channel, double T_max, node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter) {
    
    if (T_max == 0) {
        cout << "ERREUR T_max" << endl;
    }
    
    std::random_device rd;
    std::default_random_engine eng(rd());
    
    double OPtime = 0.0;
    double OPtime_dag = 0.0;
    
    OperatorEarly.channel = channel;
    OperatorLatter.channel = channel;

    int Rand_Diag_OffDiag = rand() % 2;

    if (Rand_Diag_OffDiag == 0) {
        /** --- Diagonal Vertex Addition --- **/
        int alpha = rand() % 2;
        
        OperatorEarly.alpha = alpha;
        OperatorLatter.alpha = alpha;

        if (empty_channel_and_alpha(alpha, channel) == 0) {
            /** Handle case where the worldline is currently empty **/
            upper = NULL;
            lower = NULL;
            
            if (rho_initial[alpha][channel] == 0) {
                /** Initial state is vacuum: must insert d_dagger then d to create a segment **/
                std::uniform_real_distribution<double> distr_dag(0.0, T_max);
                OPtime_dag = distr_dag(eng);
                std::uniform_real_distribution<double> distr(OPtime_dag, T_max);
                OPtime = distr(eng);
            }
            else {
                /** Initial state is occupied: must insert d then d_dagger to create a hole **/
                std::uniform_real_distribution<double> distr_dag(0.0, T_max);
                OPtime_dag = distr_dag(eng);
                std::uniform_real_distribution<double> distr(0.0, OPtime_dag);
                OPtime = distr(eng);
            }
        }
        else {
            /** Handle insertion into an existing worldline **/
            double Tlower = 0.0;
            double Tupper = T_max;
            
            int index = rand() % (nb_operator[alpha][channel] + 1);
            
            upper = NULL;
            lower = NULL;

            if (index == 0) {
                /** Insertion at the start of the list **/
                upper = alph_head[alpha][channel];
                Tlower = 0.0;
                Tupper = upper->data.time;
            }
            else if (index == nb_operator[alpha][channel]) {
                /** Insertion at the end of the list **/
                upper = NULL;
                lower = tail[alpha][channel];
                Tlower = lower->data.time;
                Tupper = T_max;
            }
            else {
                /** Insertion between two existing nodes **/
                upper = alph_head[alpha][channel];
                for (int i = 0; i < index; i = i + 1) {
                    if (i == index - 1) {
                        lower = upper;
                    }
                    upper = upper->fwd_alpha;
                }
                Tlower = lower->data.time;
                Tupper = upper->data.time;
            }

            /** Maintain worldline continuity by checking the neighbor operator types **/
            if (upper != NULL) {
                if (upper->data.flag == abs(upper->data.alpha - 1)) {
                    /** Upper node is a creation operator: sampling creation-annihilation sequence **/
                    std::uniform_real_distribution<double> distr(Tlower, Tupper);
                    OPtime = distr(eng);
                    std::uniform_real_distribution<double> distr_dag(Tlower, OPtime);
                    OPtime_dag = distr_dag(eng);
                }
                else {
                    /** Upper node is an annihilation operator **/
                    std::uniform_real_distribution<double> distr_dag(Tlower, Tupper);
                    OPtime_dag = distr_dag(eng);
                    std::uniform_real_distribution<double> distr(Tlower, OPtime_dag);
                    OPtime = distr(eng);
                }
            }
            else {
                if (abs(lower->data.flag) == abs(lower->data.alpha - 1)) {
                    std::uniform_real_distribution<double> distr_dag(Tlower, Tupper);
                    OPtime_dag = distr_dag(eng);
                    std::uniform_real_distribution<double> distr(Tlower, OPtime_dag);
                    OPtime = distr(eng);
                }
                else {
                    std::uniform_real_distribution<double> distr(Tlower, Tupper);
                    OPtime = distr(eng);
                    std::uniform_real_distribution<double> distr_dag(Tlower, OPtime);
                    OPtime_dag = distr_dag(eng);
                }
            }
        }

        /** Finalize Kink properties based on temporal order **/
        if (OPtime > OPtime_dag) {
            OperatorEarly.time = OPtime_dag;
            OperatorEarly.flag = abs(alpha - 1);
            OperatorLatter.time = OPtime;
            OperatorLatter.flag = alpha;
        }
        else {
            OperatorLatter.time = OPtime_dag;
            OperatorLatter.flag = abs(alpha - 1);
            OperatorEarly.time = OPtime;
            OperatorEarly.flag = alpha;
        }
    }
    else {
        /** --- Off-Diagonal Vertex Addition --- **/
        /** Involves simultaneous operators on H and H-tilde branches **/
        double time_OpeH = 0.0;
        double time_OpeHtilde = 0.0;
        
        node * OpeH = tail[0][OperatorEarly.channel];
        node * OpeHtilde = tail[1][OperatorLatter.channel];

        if (OpeH == NULL && OpeHtilde == NULL) {
            std::uniform_real_distribution<double> distr_dag(0., T_max);
            time_OpeH = distr_dag(eng);
            std::uniform_real_distribution<double> distr(0., T_max);
            time_OpeHtilde = distr(eng);
            OpeH = NULL;
            OpeHtilde = NULL;
        }
        else {
            if (OpeH == NULL) {
                std::uniform_real_distribution<double> distr_dag(0., T_max);
                time_OpeH = distr_dag(eng);
                std::uniform_real_distribution<double> distr(OpeHtilde->data.time, T_max);
                time_OpeHtilde = distr(eng);
                OpeH = NULL;
            }
            else if (OpeHtilde == NULL) {
                std::uniform_real_distribution<double> distr_dag(OpeH->data.time, T_max);
                time_OpeH = distr_dag(eng);
                std::uniform_real_distribution<double> distr(0., T_max);
                time_OpeHtilde = distr(eng);
                OpeHtilde = NULL;
            }
            else {
                std::uniform_real_distribution<double> distr_dag(OpeH->data.time, T_max);
                time_OpeH = distr_dag(eng);
                std::uniform_real_distribution<double> distr(OpeHtilde->data.time, T_max);
                time_OpeHtilde = distr(eng);
            }
        }

        /** Assign operator roles and flags based on the sampled times for off-diagonal updates **/
        if (OpeH != NULL && OpeHtilde != NULL) {
            if (time_OpeHtilde > time_OpeH) {
                OperatorEarly.time = time_OpeH;
                OperatorEarly.alpha = 0;
                OperatorEarly.flag = abs(OpeH->data.flag - 1);
                OperatorLatter.time = time_OpeHtilde;
                OperatorLatter.alpha = 1;
                OperatorLatter.flag = abs(OpeHtilde->data.flag - 1);
                lower = OpeH;
                upper = OpeHtilde;
            }
            else {
                OperatorLatter.time = time_OpeH;
                OperatorLatter.alpha = 0;
                OperatorLatter.flag = abs(OpeH->data.flag - 1);
                OperatorEarly.time = time_OpeHtilde;
                OperatorEarly.alpha = 1;
                OperatorEarly.flag = abs(OpeHtilde->data.flag - 1);
                lower = OpeHtilde;
                upper = OpeH;
            }
        }
        else if ((OpeH == NULL && OpeHtilde != NULL) || (OpeH != NULL && OpeHtilde == NULL)) {
            if (OpeH != NULL) {
                if (time_OpeHtilde > time_OpeH) {
                    OperatorEarly.time = time_OpeH;
                    OperatorEarly.alpha = 0;
                    OperatorEarly.flag = abs(OpeH->data.flag - 1);
                    OperatorLatter.time = time_OpeHtilde;
                    OperatorLatter.alpha = 1;
                    OperatorLatter.flag = OpeH->data.flag;
                    lower = OpeH;
                    upper = NULL;
                }
                else {
                    OperatorLatter.time = time_OpeH;
                    OperatorLatter.alpha = 0;
                    OperatorLatter.flag = abs(OpeH->data.flag - 1);
                    OperatorEarly.time = time_OpeHtilde;
                    OperatorEarly.alpha = 1;
                    OperatorEarly.flag = OpeH->data.flag;
                    lower = NULL;
                    upper = OpeH;
                }
            }
            if (OpeHtilde != NULL) {
                if (time_OpeHtilde > time_OpeH) {
                    OperatorEarly.time = time_OpeH;
                    OperatorEarly.alpha = 0;
                    OperatorEarly.flag = OpeHtilde->data.flag;
                    OperatorLatter.time = time_OpeHtilde;
                    OperatorLatter.alpha = 1;
                    OperatorLatter.flag = abs(OpeHtilde->data.flag - 1);
                    lower = NULL;
                    upper = OpeHtilde;
                }
                else {
                    OperatorLatter.time = time_OpeH;
                    OperatorLatter.alpha = 0;
                    OperatorLatter.flag = OpeHtilde->data.flag;
                    OperatorEarly.time = time_OpeHtilde;
                    OperatorEarly.alpha = 1;
                    OperatorEarly.flag = abs(OpeHtilde->data.flag - 1);
                    lower = OpeHtilde;
                    upper = NULL;
                }
            }
        }
        else {
            /** Handle case for empty list based on initial density matrix density **/
            if (rho_initial[0][channel] == 0 && rho_initial[1][channel] == 0) {
                if (time_OpeHtilde > time_OpeH) {
                    OperatorEarly.time = time_OpeH;
                    OperatorEarly.alpha = 0;
                    OperatorEarly.flag = 1;
                    OperatorLatter.time = time_OpeHtilde;
                    OperatorLatter.alpha = 1;
                    OperatorLatter.flag = 0;
                    lower = OpeH;
                    upper = OpeHtilde;
                }
                else {
                    OperatorLatter.time = time_OpeH;
                    OperatorLatter.alpha = 0;
                    OperatorLatter.flag = 1;
                    OperatorEarly.time = time_OpeHtilde;
                    OperatorEarly.alpha = 1;
                    OperatorEarly.flag = 0;
                    lower = OpeHtilde;
                    upper = OpeH;
                }
            }
            else {
                if (time_OpeHtilde > time_OpeH) {
                    OperatorEarly.time = time_OpeH;
                    OperatorEarly.alpha = 0;
                    OperatorEarly.flag = 0;
                    OperatorLatter.time = time_OpeHtilde;
                    OperatorLatter.alpha = 1;
                    OperatorLatter.flag = 1;
                    lower = OpeH;
                    upper = OpeHtilde;
                }
                else {
                    OperatorLatter.time = time_OpeH;
                    OperatorLatter.alpha = 0;
                    OperatorLatter.flag = 0;
                    OperatorEarly.time = time_OpeHtilde;
                    OperatorEarly.alpha = 1;
                    OperatorEarly.flag = 1;
                    lower = OpeHtilde;
                    upper = OpeH;
                }
            }
        }
    }
}

/**
 * @brief Commits the proposed vertex addition by inserting new nodes into the linked lists.
 */
void clist::Move_Adding(node *& upper, node *& lower, kink & OperatorEarly, kink & OperatorLatter) {
    if (OperatorEarly.alpha == OperatorLatter.alpha) {
        node * previous_time = NULL;
        node * previous_node = NULL;
        
        /** Insert the first kink of the diagonal vertex **/
        RecherchePreviousNode(OperatorEarly, previous_node, lower, previous_time);
        add_node(lower, previous_time, previous_node, OperatorEarly);
        
        /** Update traversal pointer and insert the second kink **/
        if (lower == NULL) {
            lower = alph_head[OperatorEarly.alpha][OperatorEarly.channel];
        }
        else {
            lower = lower->fwd_alpha;
        }
        RecherchePreviousNode(OperatorLatter, previous_node, lower, previous_time);
        add_node(lower, previous_time, previous_node, OperatorLatter);
    }
    else {
        node * previous_time = NULL;
        node * previous_node = NULL;
        
        /** Insert off-diagonal kinks into their respective alpha-specific lists **/
        RecherchePreviousNode(OperatorEarly, previous_node, lower, previous_time);
        add_node(lower, previous_time, previous_node, OperatorEarly);
        
        RecherchePreviousNode(OperatorLatter, previous_node, upper, previous_time);
        add_node(upper, previous_time, previous_node, OperatorLatter);
    }
}

/**
 * @brief Calculates the segment overlap between different spin channels.
 */
void clist::Overlap_SigSigBar(double & O_sig, double T_max, int alpha, int channel, double time_lower, double time_upper) {
    Overlap_HHtilde(O_sig, T_max, abs(alpha - 1), abs(channel - 1), time_lower, time_upper);
}

/**
 * @brief Integrates the temporal overlap between H and H-tilde configurations for a given channel.
 * * Scans the worldline of the opposite alpha branch to find periods of simultaneous occupation.
 */
void clist::Overlap_HHtilde(double & W_sig, double T_max, int alpha, int channel, double time_lower, double time_upper) {
    node * readNode = alph_head[abs(alpha - 1)][channel];
    
    if (readNode == NULL) {
        if (rho_initial[abs(alpha - 1)][channel] == 1) {
            W_sig = W_sig + time_upper - time_lower;
        }
        else {
            W_sig = W_sig + 0.0;
        }
    }
    else {
        /** Synchronize readNode position with the start of the temporal segment **/
        if (readNode->data.alpha == readNode->data.flag) {
            if (readNode->data.time < time_lower) {
                readNode = readNode->fwd_alpha;
            }
            else {
                if (time_upper < readNode->data.time) {
                    W_sig = W_sig + time_upper - time_lower;
                }
                else {
                    W_sig = W_sig - time_lower + readNode->data.time;
                }
                readNode = readNode->fwd_alpha;
            }
        }
        
        /** Iterate through creation/annihilation pairs to accumulate overlap durations **/
        if (readNode != NULL) {
            while (readNode->data.time < time_upper) {
                if (readNode->data.time < time_lower) {
                    if (readNode->fwd_alpha == NULL) {
                        W_sig = W_sig + time_upper - time_lower;
                        break;
                    }
                    else {
                        if (((readNode->fwd_alpha)->data.time) <= time_upper && ((readNode->fwd_alpha)->data.time) > time_lower) {
                            W_sig = W_sig + (readNode->fwd_alpha)->data.time - time_lower;
                        }
                        if (((readNode->fwd_alpha)->data.time) > time_upper) {
                            W_sig = W_sig + time_upper - time_lower;
                        }
                    }
                }
                else {
                    if (readNode->fwd_alpha == NULL) {
                        W_sig = W_sig + time_upper - readNode->data.time;
                        break;
                    }
                    else {
                        if (((readNode->fwd_alpha)->data.time) < time_upper || ((readNode->fwd_alpha)->data.time) == time_upper) {
                            W_sig = W_sig + (readNode->fwd_alpha)->data.time - readNode->data.time;
                        }
                        if (((readNode->fwd_alpha)->data.time) > time_upper) {
                            W_sig = W_sig + time_upper - readNode->data.time;
                        }
                    }
                }
                
                if ((readNode->fwd_alpha)->fwd_alpha != NULL) {
                    readNode = (readNode->fwd_alpha)->fwd_alpha;
                }
                else {
                    break;
                }
            }
        }
    }
}






void clist::Trace_Adding(double T_max,node *& upper ,node *& lower, kink & OperatorEarly, kink & OperatorLatter  , std::vector<double> & l_sig , std::vector<double> & W_sig, std::vector<double> & O_sig , int AddRemove ){
    
    if (OperatorEarly.alpha == OperatorLatter.alpha){
        
        if (lower==NULL){
            
            if (OperatorEarly.flag == OperatorEarly.alpha){
                
                if (upper !=NULL){
                    
                    l_sig[OperatorEarly.alpha] =  OperatorEarly.time - OperatorLatter.time;
                }
                else{
                    
                    l_sig[OperatorEarly.alpha] = OperatorEarly.time  - OperatorLatter.time;
                }
            }
            else{
                
                l_sig[OperatorEarly.alpha] = OperatorLatter.time - OperatorEarly.time;
            }
        
        }
        else{
            if (lower ->data.alpha == lower ->data.flag){
                
                l_sig[OperatorEarly.alpha] = OperatorLatter.time - OperatorEarly.time;
            }
            else{
                
                if (upper==NULL){
                    
                    l_sig[OperatorEarly.alpha] = OperatorEarly.time - lower->data.time + T_max - OperatorLatter.time - (T_max- lower->data.time );
                }
                else{
                    
                    l_sig[OperatorEarly.alpha] = OperatorEarly.time - lower->data.time + upper-> data.time- OperatorLatter.time - ( upper-> data.time - lower->data.time );
                }
            }
        }
        
        if (lower==NULL){
            if (upper==NULL){
                
                if (OperatorEarly.flag == OperatorEarly.alpha){
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time ,OperatorLatter.time );
                    
                    W_sig[OperatorEarly.channel] = -1.* W_sig[OperatorEarly.channel] ;
                    
                }
                else{
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time  ,OperatorLatter.time );
                }
            }
            else{
                
                if (OperatorEarly.flag == OperatorEarly.alpha){
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time  ,OperatorLatter.time );
                    W_sig[OperatorEarly.channel] = -1.*W_sig[OperatorEarly.channel];
                }
                else{
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time  ,OperatorLatter.time );
 
                }
            }
        }
        else{
            if (lower ->data.alpha == lower ->data.flag){
                
                Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, lower->data.alpha , lower->data.channel,OperatorEarly.time ,OperatorLatter.time );
            }
            else{
                
                Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max, lower->data.alpha , lower->data.channel,OperatorEarly.time ,OperatorLatter.time );
                W_sig[OperatorEarly.channel] = -1.*W_sig[OperatorEarly.channel];
            }
        }
        
        
        
        if (lower==NULL){
            if (upper==NULL){
                
                if (OperatorEarly.flag == OperatorEarly.alpha){
                    
                    Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time ,OperatorLatter.time );
                    O_sig[OperatorEarly.alpha] = -1. * O_sig[OperatorEarly.alpha];
                    
                }
                else{
                    
                    Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time  ,OperatorLatter.time );
                }
            }
            else{
                
                if (OperatorEarly.flag == OperatorEarly.alpha){
                    
                    Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time, OperatorLatter.time );
                    O_sig[OperatorEarly.alpha] = O_sig[OperatorEarly.alpha] *(-1.);
                   
                }
                else{
                    
                    Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time, OperatorLatter.time );
                }
            }
        }
        else{
            if (lower->data.alpha == lower->data.flag){
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time ,OperatorLatter.time );

            }
            else{
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max, OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time ,OperatorLatter.time );
                O_sig[OperatorEarly.alpha]= -1. * O_sig[OperatorEarly.alpha];
            }
        }

        
    }
    else{
        
        if (lower==NULL){
            if (OperatorEarly.flag == OperatorEarly.alpha){
                
                l_sig[OperatorEarly.alpha] = -(T_max - OperatorEarly.time);
            }
            else{
                
                l_sig[OperatorEarly.alpha] = T_max - OperatorEarly.time;
            }
        }
        else{
            if (lower->data.flag == lower->data.alpha){
                
                l_sig[OperatorEarly.alpha] = T_max -  OperatorEarly.time;
            }
            else{
                
                l_sig[OperatorEarly.alpha] =    - (T_max - OperatorEarly.time );
            }
        }
        
        if (upper==NULL){
            if (OperatorLatter.flag == OperatorLatter.alpha){
                
                l_sig[OperatorLatter.alpha] = - (T_max - OperatorLatter.time);
            }
            else{
                l_sig[OperatorLatter.alpha] = T_max - OperatorLatter.time;
            }
        }
        else{
            if (upper->data.flag == upper->data.alpha){
                
                l_sig[OperatorLatter.alpha] = T_max -  OperatorLatter.time;
            }
            else{
                
                l_sig[OperatorLatter.alpha] =    -(T_max -  OperatorLatter.time) ;
            }
        }
        
        
        
        if (lower==NULL and upper==NULL){
            if (OperatorEarly.alpha == OperatorEarly.flag){
                
                W_sig[OperatorEarly.channel] = - (T_max - OperatorEarly.time) ;
            }
            else{
                
                W_sig[OperatorEarly.channel] = T_max - OperatorLatter.time;
            }
        }
        
        else if (lower != NULL and upper !=NULL){
            if (lower->data.alpha == lower->data.flag){
                
                if (OperatorEarly.time > upper->data.time){
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel]+  T_max - OperatorLatter.time;
                }
                else{
                    
                    if (AddRemove == 0){
                        W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + T_max - OperatorLatter.time;
                    }
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time,T_max );
                }
            }
            else{
                
                if (OperatorEarly.time > upper->data.time){
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + OperatorEarly.time - upper->data.time;
                    if (AddRemove==0){
                        W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    }
                    else{
                        
                        W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    }
                }
                else{
                    
                    if (AddRemove==0){
                        W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    }
                    else{
                        
                        W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    }
                    
                    double overlap_EarlyUpper =0.0;
                    Overlap_HHtilde(overlap_EarlyUpper,T_max,OperatorEarly.alpha , OperatorEarly.channel,OperatorEarly.time, upper->data.time );
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - overlap_EarlyUpper;
                }
            }
        }
        
        else if (lower==NULL){
            
            if (upper->data.alpha == upper->data.flag ){
                
                if (OperatorEarly.time > upper->data.time){
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + (T_max - OperatorLatter.time);
                }
                else{
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + (T_max - OperatorLatter.time);
                    
                    Overlap_HHtilde(W_sig[OperatorEarly.channel],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time,upper->data.time );
                }
 
            }
            else{
                
                if (OperatorEarly.time > upper->data.time){
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + OperatorEarly.time - upper->data.time;
                }
                else{
                    
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] - (T_max - upper->data.time);
                    double overlap_EarlyUpper =0.0;
                    Overlap_HHtilde(overlap_EarlyUpper,T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time,upper->data.time);
                    W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel]  - overlap_EarlyUpper;
                }
            }
        }
        else if (upper==NULL){
            
            if (lower->data.alpha == lower->data.flag ){
                
                W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + (T_max - OperatorLatter.time);

            }
            else{
                W_sig[OperatorEarly.channel] = W_sig[OperatorEarly.channel] + (OperatorEarly.time - lower->data.time)- (T_max - lower->data.time);
            }
        }
        
        
        
        if (lower==NULL){
            if (OperatorEarly.flag==OperatorEarly.alpha){
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time, T_max );
                O_sig[OperatorEarly.alpha] = -1. * O_sig[OperatorEarly.alpha] ;
            }
            else{
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time,T_max );
            }
        }
        else{
            if (lower->data.flag == lower->data.alpha){
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time , T_max);
            }
            else{
                
                Overlap_SigSigBar(O_sig[OperatorEarly.alpha],T_max,OperatorEarly.alpha , OperatorEarly.channel, OperatorEarly.time , T_max);
                O_sig[OperatorEarly.alpha] = -1.* O_sig[OperatorEarly.alpha];
            }
        }
        
        if (upper==NULL){
            if (OperatorLatter.flag==OperatorLatter.alpha){
                
                Overlap_SigSigBar(O_sig[OperatorLatter.alpha],T_max,OperatorLatter.alpha , OperatorLatter.channel, OperatorLatter.time, T_max );
                O_sig[OperatorLatter.alpha] = -1. * O_sig[OperatorLatter.alpha] ;
            }
            else{
                
                Overlap_SigSigBar(O_sig[OperatorLatter.alpha],T_max,OperatorLatter.alpha , OperatorLatter.channel, OperatorLatter.time,T_max );
            }
        }
        else{
            if (upper->data.flag == upper->data.alpha){
                
                Overlap_SigSigBar(O_sig[OperatorLatter.alpha],T_max,OperatorLatter.alpha , OperatorLatter.channel, OperatorLatter.time , T_max);
            }
            else{
                
                Overlap_SigSigBar(O_sig[OperatorLatter.alpha],T_max,OperatorLatter.alpha , OperatorLatter.channel, OperatorLatter.time , T_max);
                O_sig[OperatorLatter.alpha] =-1.* O_sig[OperatorLatter.alpha];
            }
        }
        

    }
    
}
