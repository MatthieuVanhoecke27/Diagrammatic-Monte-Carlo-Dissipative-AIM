//
//  main.cpp
//  DiagrammaticMonteCarlo Dephasing
//
//  Created by adm.matthieu.vanhoecke on 02/01/2023.
//

#include <iostream>
#include<cstdio>
using namespace std;
#include <random>
#include <complex>
typedef complex<double> Complex;
#include <armadillo>
using namespace arma;
typedef complex<double> Complex;
#include <algorithm>

#include "FonctionPrincipale.hpp"
#include "HybridizationFunction.hpp"
#include "BathOpe Determinant.hpp"

/**
 * Calculates the complex average of a vector of complex numbers.
 * @param v Vector containing the complex values.
 * @param sizeTab Reference to the number of elements to process.
 * @return The arithmetic mean as a complex<double>.
 */
complex<double> getAverage(vector<complex<double>> v, double & sizeTab){
    Complex Complexi(0.0,1.0);
    double meanRealpart=0.0;
    double meanImagPart=0.0;
 
    for (int i=0; i<sizeTab;i=i+1){
        meanRealpart = meanRealpart + real(v[i]);
        meanImagPart = meanImagPart + imag(v[i]);
    
    }
    return (meanRealpart/sizeTab + Complexi * meanImagPart/sizeTab);
}

/**
 * Computes the statistical average of the Density Matrix elements across the simulation samples.
 * @param MeanDensityMatrix Output matrix storing the averaged results for a specific time index.
 * @param TabDensityMatrix Input matrix containing all sampled density matrices.
 * @param sizeTab Number of samples (rows) in TabDensityMatrix.
 * @param index Current time-slice index (row) in the result matrix.
 */
void getAverage_DensityMatrix(Mat<Complex> & MeanDensityMatrix,Mat<Complex> & TabDensityMatrix,double & sizeTab , int & index){
    
    Complex Complexi(0.0,1.0);
    double meanRealpart=0.0;
    double meanImagPart=0.0;
    for (int k=0;k<4;k=k+1){
        meanRealpart=0.0;
        meanImagPart=0.0;
        for (int i=0; i<sizeTab;i=i+1){
            meanRealpart = meanRealpart + real(TabDensityMatrix(i,k));
            meanImagPart = meanImagPart + imag(TabDensityMatrix(i,k));
         
        }
        MeanDensityMatrix(index,k) = (meanRealpart + Complexi * meanImagPart)/sizeTab ;
    }
}

/**
 * Records the state of the density matrix based on the current worldline configuration.
 * Maps the number of operators per channel to specific basis states in the Hilbert space.
 * @param ListTime Linked-list containing the temporal operator configuration.
 * @param TabDensityMatrix Matrix to store the state components.
 * @param Trace_impurity Current trace value of the impurity.
 * @param Phase_config Complex weight/phase of the current Monte Carlo configuration.
 * @param indice Sample index for storage.
 */
void SavedensityMatrix(clist & ListTime, Mat<Complex> & TabDensityMatrix , Complex Trace_impurity, Complex Phase_config,int indice){
    for (int i=0 ; i<4 ; i=i+1){
        TabDensityMatrix(indice,i)=0.0;
    }

    if (rho_initial[0][0]== 0 and rho_initial[1][0]==0 and rho_initial[0][1]==0 and rho_initial[1][1]==0){

        if (ListTime.nb_operator[0][0]%2==0){
            if (ListTime.nb_operator[0][1]%2==0){
                TabDensityMatrix(indice,0) = Phase_config;
            }
            else{
                TabDensityMatrix(indice,2)= Phase_config;
            }
        }
         else{
             if (ListTime.nb_operator[0][1]%2==0){
                 TabDensityMatrix(indice,1) = Phase_config;
             }
             else{
                 TabDensityMatrix(indice,3)= Phase_config;
             }
         }
     }
 }

/**
 * Computes the expectation value of the number operator n_sigma for a specific channel.
 * Evaluates the operator based on the parity of vertices (diagonal vs off-diagonal).
 * @return The complex mean value of the number operator.
 */
complex<double> MeanValue_NumberOperator(clist & Listime,int & channel, Complex Trace_impurity, Complex Phase_config ){
    complex<double> Complexi(0.0,1.0);

    if (rho_initial[0][channel]==0){
        if ((Listime.nb_operator[0][channel])%2==0){
            return 0.0;
        }
        else{
            
            Trace_impurity=1.;
            return  Phase_config ;
        }
    }
    else{
        if ((Listime.nb_operator[0][channel])%2!=0){
            return 0.0;
        }
        else{
            Trace_impurity=1.;
            return  Phase_config * (real(Trace_impurity) - Complexi * imag(Trace_impurity)) / (real(Trace_impurity)* real(Trace_impurity)+ imag(Trace_impurity)* imag(Trace_impurity));
        }
    }
}

/**
 * Calculates the contribution of the impurity trace to the configuration weight.
 * Accounts for epsilon shifts, interaction U, and dephasing rates (gamma).
 * @param EffectifHamltonian Toggle for Jump operators: 0 (with jumps), 1 (no jumps/effective non-hermitian evolution).
 * @return Complex exponential factor representing the impurity trace contribution.
 */
complex<double> Trace_impurity(vector<double> & l_sig , vector<double> & W_sig,vector<double> & O_sig, double & U , double & EpsilonI , double & gammaDeph_up,double & gammaDeph_down , int EffectifHamltonian){
    if (EffectifHamltonian==1){
        W_sig[0] = 0.0;
        W_sig[1] = 0.0;
    }
    complex<double> Complexi(0.0,1.0);
    Complex Argument = -Complexi * EpsilonI * (l_sig[0] - l_sig[1]) -Complexi * U *(O_sig[0]-O_sig[1]) + gammaDeph_up* W_sig[0] + gammaDeph_down *W_sig[1] -  0.5*gammaDeph_up * l_sig[0] -  0.5*gammaDeph_down * l_sig[1];
    return exp( Argument);
}

/**
 * Computes the configuration sign resulting from the commutation of H and Htilde operators.
 * Iterates through the linked list to determine the fermionic permutation sign.
 * @return The resulting sign (+1.0 or -1.0).
 */
double Sign_HHtilde(clist & ListTime,int channel){
    double SignConfig=1.;
    node * ReadNode = ListTime.head;
    int NombreOpeH=0;
    int rho_initial[2][nbr_channel]={0,0,0,0};
    
    int NombreOpe_Htilde_UP =0;
    int NombreOpe_Htilde_DOWN =0;
    int NombreOpe_H_UP =0;
    int NombreOpe_H_DOWN =0;
    
    int NombreOpeH_channel =0;
    int NombreOpeH_channelBar =0;

    if (rho_initial[0][0]==1){
        NombreOpe_H_UP=1;
    }
    if (rho_initial[0][1]==1){
        NombreOpe_H_DOWN=1;
    }
    if (rho_initial[1][0]==1){
        NombreOpe_Htilde_UP=1;
    }
    if (rho_initial[1][1]==1){
        NombreOpe_Htilde_DOWN=1;
    }
    double signeTrace=1.0;
    while (ReadNode !=NULL){
        if ((ReadNode->data.alpha == 0)){
            NombreOpeH = NombreOpeH +1;
        }
        if ((ReadNode->data.alpha == 1)){
            if (NombreOpeH %2 !=0){
                SignConfig = SignConfig *(-1.);
            }
        }

    if (ReadNode->data.alpha == 0){
            if (ReadNode->data.channel==0){
        if (ReadNode->data.flag==0){
            NombreOpe_H_UP=0;
        }
        else{
             NombreOpe_H_UP=1;
        }
         }
        else{
        if (ReadNode->data.flag ==0){
            if ((NombreOpe_H_UP + NombreOpe_Htilde_UP)%2 !=0){
            signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_H_DOWN=0;
        }
        else{
                    if ((NombreOpe_H_UP + NombreOpe_Htilde_UP)%2 !=0){
                        signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_H_DOWN=1;
        }
         }
    }
    else{
            if (ReadNode->data.channel==0){
                if (ReadNode->data.flag ==1){
                    if (NombreOpe_H_UP ==1){
                        signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_Htilde_UP =0;
        }
                else{
            if (NombreOpe_H_UP ==1){
                        signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_Htilde_UP =1;
        }
        }
        else{
            if (ReadNode->data.flag ==1){
            if ((NombreOpe_H_UP +NombreOpe_H_DOWN+ NombreOpe_Htilde_UP)%2 !=0){
                        signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_Htilde_DOWN=0;
        }
        else{
            if ((NombreOpe_H_UP +NombreOpe_H_DOWN+ NombreOpe_Htilde_UP)%2 !=0){
                        signeTrace = signeTrace * (-1.0);
                    }
                    NombreOpe_Htilde_DOWN=1;
        }
        }
        }                                                                                                   ReadNode = ReadNode ->p_next;
    }
    SignConfig=1.0;
    return SignConfig*signeTrace;
}

/**
 * Calculates the fermionic sign associated with the specific ordering of Psi(t) and Psi(t_bar).
 * Identifies double occupancy/doublons within the worldline to adjust the weight sign.
 * @return The sign (+1.0 or -1.0) derived from time-ordering permutations.
 */
double Sign_t_tBar(clist & ListTime,int channel){
    double Signe_Between_ttBar=1.;
    int NombreDoublon=0;
    node * readNode = ListTime.head;
    int flag;
    int Nombre=0;
    if (readNode !=NULL){
        while (readNode !=NULL){
            if (readNode->data.channel ==channel){
                if (Nombre!=0){
                    if (flag==0){
                        if (readNode->data.flag == 1){
                            Signe_Between_ttBar = Signe_Between_ttBar * (-1.0);
                        }
                        else{
                            NombreDoublon = NombreDoublon+1;
                        }
                    }
                    else{
                        if (readNode->data.flag == 0){
                            Signe_Between_ttBar = Signe_Between_ttBar * 1.0;
                        }
                        else{
                            NombreDoublon = NombreDoublon+1;
                        }
                    }
                    Nombre=0;
                }
                else{
                    flag = readNode->data.flag;
                    Nombre=Nombre+1;
                }
            }
            readNode = readNode->p_next;
        }
    }

    if (NombreDoublon %2 ==0 and NombreDoublon !=0){
        Signe_Between_ttBar = Signe_Between_ttBar * pow(-1.,NombreDoublon/2 );
    }
    return Signe_Between_ttBar;
}


/**
 * @brief Performs the "Adding" update in the Diagrammatic Monte Carlo simulation.
 * * This subroutine proposes the insertion of a new vertex pair (operator and its conjugate)
 * into the configuration. It handles the transition probability calculation,
 * the impurity trace ratio, and the hybridization determinant update.
 * * @param ListTime The linked-list structure representing the worldline configuration.
 * @param ParametresMC Vector containing physical/MC parameters: [t_max, U, EpsilonI, gammaUp, gammaDown].
 * @param InvMatrix_sig Inverse hybridization matrix for channel 0.
 * @param InvMatrix_sigBar Inverse hybridization matrix for channel 1.
 * @param k_accep Counter for the number of accepted moves.
 * @param ParametresBath Parameters defining the hybridization function of the bath.
 * @param ArgumentProba The cumulative phase/sign of the Monte Carlo weight.
 * @param Traceimpurity Current cumulative value of the impurity trace.
 * @param SigneBetween_HHtilde_Channel Tracker for the fermionic sign (H/Htilde) for channel 0.
 * @param SigneBetween_HHtilde_ChannelBar Tracker for the fermionic sign (H/Htilde) for channel 1.
 * @param Signe_TtBarOld_channel Tracker for the time-ordering sign for channel 0.
 * @param Signe_TtBarOld_channelBar Tracker for the time-ordering sign for channel 1.
 * @param EffectifHamltonian Toggle for Jump operators (0: with jumps, 1: effective Hamiltonian).
 */
void MonteCarlo_Adding(clist & ListTime, vector<double> ParametresMC , Mat<Complex> & InvMatrix_sig ,Mat<Complex> & InvMatrix_sigBar ,  int & k_accep , vector<double> ParametresBath ,Complex & ArgumentProba , Complex & Traceimpurity, double & SigneBetween_HHtilde_Channel, double & SigneBetween_HHtilde_ChannelBar ,double & Signe_TtBarOld_channel,double & Signe_TtBarOld_channelBar, int EffectifHamltonian ){
    std::random_device rd;
    std::default_random_engine eng(rd());

    /* --- Local variable initialization --- */
    node * upper=NULL;
    node * lower=NULL;
    kink OperatorEarly;
    kink  OperatorLatter;
    int Index_Col = 0;
    int Index_Line= 0;
    Col<Complex> New_Col;
    Row<Complex> New_Row;
    Complex DiagElement =0.;
    complex<double> Complexi(0.0,1.0);
    
    std::vector<double> l_sig= {0. , 0.};
    std::vector<double> W_sig= {0. , 0.};
    std::vector<double> O_sig= {0. , 0.};

    if (ParametresMC[0]==0){
        cout<< "ERREUR T_max" <<endl;
    }
    
    int channel = rand()%2;
    Complex ProbabilityAccep =1.;

    /* --- Case 1: Inserting into an empty channel --- */
    if (ListTime.nb_operator[0][channel]==0 and ListTime.nb_operator[1][channel]==0){

        ListTime.propose_adding(channel,ParametresMC[0], upper,lower,OperatorEarly,OperatorLatter);
        
        vector<vector<int>> NombreOpe{{ListTime.nb_operator[0][0], ListTime.nb_operator[0][1]},{ListTime.nb_operator[1][0],ListTime.nb_operator[1][1]}};
        int move=0;
        ProbabilityAccep = ProbabilityAccep * ListTime.Probability_move(move, OperatorEarly, OperatorLatter,upper, lower, channel, ParametresMC[0], NombreOpe);

        int AddRemove =0;
        ListTime.Trace_Adding(ParametresMC[0],upper,lower,OperatorEarly,OperatorLatter, l_sig,W_sig,O_sig,AddRemove);
        Complex TrImp = Trace_impurity(l_sig ,W_sig,O_sig,ParametresMC[1] ,ParametresMC[2] ,ParametresMC[3],ParametresMC[4],EffectifHamltonian);
        ProbabilityAccep = ProbabilityAccep * TrImp;
        
        if (OperatorEarly.flag==0){
            DiagElement = Hybridization_Funct(OperatorEarly.alpha, OperatorLatter.alpha , OperatorEarly.time , OperatorLatter.time, ParametresBath);
        }
        else{
            DiagElement = Hybridization_Funct( OperatorLatter.alpha,OperatorEarly.alpha ,  OperatorLatter.time,OperatorEarly.time, ParametresBath );
        }
        ProbabilityAccep = ProbabilityAccep * DiagElement;

        std::uniform_real_distribution<double> distr_dag(0.,1.);
        double RandMC = distr_dag(eng);
        
        if (RandMC < min(1.,abs(ProbabilityAccep))){
            k_accep = k_accep+1 ;
            
            ListTime.Move_Adding(upper, lower, OperatorEarly,OperatorLatter);
            Mat<Complex> InverseMat(" 1. ");
            InverseMat[0] = (real(DiagElement) - Complexi * imag(DiagElement))/(real(DiagElement)* real(DiagElement) + imag(DiagElement) * imag(DiagElement));
            
            if (channel==0){
                InvMatrix_sig =  InverseMat;
            }
            else{
                InvMatrix_sigBar =  InverseMat;
            }

            Complex sign_CprimeC = - Complexi;

            if (OperatorEarly.alpha == OperatorLatter.alpha){
                node * EarlyNode= ListTime.alph_head[OperatorEarly.alpha][channel];
                sign_CprimeC = sign_CprimeC * ListTime.Sign_TimeOrdering(EarlyNode, EarlyNode->fwd_alpha);
            }
            else{
                sign_CprimeC = sign_CprimeC * ListTime.Sign_TimeOrdering(ListTime.alph_head[OperatorEarly.alpha][channel], ListTime.alph_head[OperatorLatter.alpha][channel]);
            }

            if (OperatorEarly.alpha != OperatorLatter.alpha){
                if (OperatorEarly.alpha==0){
                    if (OperatorEarly.flag==0){
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                    else{
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                }
                else{
                    if (OperatorEarly.flag==0){
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                    else{
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                }
            }
            
            if (channel==0){
                double SigneHHtilde_New = Sign_HHtilde(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_Channel);
                SigneBetween_HHtilde_Channel = SigneHHtilde_New;
                SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
            }
            else{
                double SigneHHtilde_New = Sign_HHtilde(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_ChannelBar);
                SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
                SigneBetween_HHtilde_Channel = SigneHHtilde_New;
            }

            if (channel==0){
                double SignTTbarNew= Sign_t_tBar(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SignTTbarNew/Signe_TtBarOld_channel);
                Signe_TtBarOld_channel = SignTTbarNew;
            }
            else{
                double SignTTbarNew= Sign_t_tBar(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SignTTbarNew/Signe_TtBarOld_channelBar);
                Signe_TtBarOld_channelBar = SignTTbarNew;
            }

            Traceimpurity = Traceimpurity * TrImp;
            ArgumentProba = ArgumentProba * sign_CprimeC * DiagElement * TrImp / abs(sign_CprimeC * DiagElement * TrImp);
        }
    }
    
    /* --- Case 2: Inserting into a non-empty channel (Updating Hybridization Matrix) --- */
    else{
        ListTime.propose_adding(channel,ParametresMC[0], upper,lower,OperatorEarly,OperatorLatter);

        vector<vector<int>> NombreOpe{{ListTime.nb_operator[0][0], ListTime.nb_operator[0][1]},{ListTime.nb_operator[1][0],ListTime.nb_operator[1][1]}};
        int move=0;
        ProbabilityAccep = ProbabilityAccep * ListTime.Probability_move(move, OperatorEarly, OperatorLatter,upper, lower, channel, ParametresMC[0], NombreOpe);
        
        int AddRemove=0;
        ListTime.Trace_Adding(ParametresMC[0],upper,lower,OperatorEarly,OperatorLatter, l_sig,W_sig,O_sig,AddRemove);
        Complex TrImp = Trace_impurity(l_sig ,W_sig,O_sig,ParametresMC[1] ,ParametresMC[2] ,ParametresMC[3],ParametresMC[4],EffectifHamltonian);
        ProbabilityAccep = ProbabilityAccep * TrImp;
        
        std::tuple<Col<Complex>,Row<Complex>> New  = Creation_ColumnRow(ListTime,OperatorEarly,OperatorLatter , Index_Col , Index_Line,ParametresBath );
         New_Col = std::get<0>(New);
         New_Row = std::get<1>(New);
        
        if (OperatorEarly.flag==0){
            DiagElement = Hybridization_Funct(OperatorEarly.alpha, OperatorLatter.alpha , OperatorEarly.time , OperatorLatter.time, ParametresBath);
        }
        else{
            DiagElement = Hybridization_Funct( OperatorLatter.alpha,OperatorEarly.alpha ,  OperatorLatter.time,OperatorEarly.time, ParametresBath );
        }
        
        Complex RationDet = 0.;
        if (channel==0){
            RationDet=  RatioDet_adding( Index_Line ,Index_Col,InvMatrix_sig , New_Col,New_Row,DiagElement);
        }
        else{
            RationDet=  RatioDet_adding(Index_Line, Index_Col ,InvMatrix_sigBar , New_Col,New_Row,DiagElement);
        }
        ProbabilityAccep =  ProbabilityAccep * RationDet;
        
        std::uniform_real_distribution<double> distr_dag(0.,1.);
        double RandMC = distr_dag(eng);
        if (RandMC < min(1.,abs(ProbabilityAccep))){
            Complex sign_CprimeC = -Complexi ;
            
            if (OperatorEarly.alpha == OperatorLatter.alpha){
                if (ListTime.nb_operator[OperatorEarly.alpha][OperatorEarly.channel]%2 !=0){
                    if (ListTime.tail[0][OperatorEarly.channel] ->data.time > ListTime.tail[1][OperatorEarly.channel] ->data.time){
                        node * dernierope = ListTime.tail[0][OperatorEarly.channel];
                        if (OperatorLatter.alpha ==1){
                            if (OperatorLatter.time > dernierope->data.time){
                                sign_CprimeC = sign_CprimeC * (1.);
                            }
                        }
                    }
                    else{
                        node * dernierope = ListTime.tail[1][OperatorEarly.channel];
                        if (OperatorLatter.alpha ==0){
                            if (OperatorLatter.time > dernierope->data.time){
                                sign_CprimeC = sign_CprimeC * (1.);
                            }
                        }
                    }
                }
            }

            k_accep = k_accep+1 ;
            ListTime.Move_Adding(upper, lower, OperatorEarly,OperatorLatter);
            
            if (channel==0){
                FastUpdate_adding(Index_Line,Index_Col,InvMatrix_sig,New_Col ,New_Row, RationDet);
            }
            else{
                FastUpdate_adding(Index_Line,Index_Col,InvMatrix_sigBar,New_Col ,New_Row, RationDet);
            }
            
            if (OperatorEarly.alpha == OperatorLatter.alpha){
                sign_CprimeC =  sign_CprimeC * ListTime.Sign_TimeOrdering(lower, lower->fwd_alpha);
            }
            else{
                if (lower==NULL and upper==NULL){
                    sign_CprimeC =  sign_CprimeC * ListTime.Sign_TimeOrdering(ListTime.alph_head[OperatorEarly.alpha][channel], ListTime.alph_head[OperatorLatter.alpha][channel]);
                }
                else if (lower==NULL){
                    sign_CprimeC =  sign_CprimeC *ListTime.Sign_TimeOrdering(ListTime.alph_head[OperatorEarly.alpha][channel], upper->fwd_alpha);
                }
                else if (upper==NULL){
                    sign_CprimeC =  sign_CprimeC * ListTime.Sign_TimeOrdering(lower->fwd_alpha, ListTime.alph_head[OperatorLatter.alpha][channel]);
                }
                else{
                    sign_CprimeC =  sign_CprimeC * ListTime.Sign_TimeOrdering(lower->fwd_alpha, upper->fwd_alpha);
                }
            }

            if (OperatorEarly.alpha != OperatorLatter.alpha){
                if (ListTime.nb_operator[OperatorEarly.alpha][OperatorEarly.channel] %2 ==0){
                    if (lower->data.time > upper->data.time){
                        if (lower->data.alpha==0){
                            sign_CprimeC = sign_CprimeC * (-Complexi);
                        }
                        else{
                            sign_CprimeC = sign_CprimeC * (-Complexi);
                        }
                    }
                    else{
                        if (upper->data.alpha==0){
                            sign_CprimeC = sign_CprimeC * (-Complexi);
                        }
                        else{
                            sign_CprimeC = sign_CprimeC * (-Complexi);
                        }
                    }
                }
                else{
                    if (OperatorEarly.alpha == 0){
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                    else{
                        sign_CprimeC = sign_CprimeC * (Complexi);
                    }
                }
            }

            if (channel==0){
                double SigneHHtilde_New = Sign_HHtilde(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_Channel);
                SigneBetween_HHtilde_Channel = SigneHHtilde_New;
                SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
            }
            else{
                double SigneHHtilde_New = Sign_HHtilde(ListTime,channel);
                sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_ChannelBar);
                SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
                SigneBetween_HHtilde_Channel = SigneHHtilde_New;
            }

            if (channel==0){
                double SignTTbarNew= Sign_t_tBar(ListTime,channel);
                sign_CprimeC = sign_CprimeC * SignTTbarNew/Signe_TtBarOld_channel;
                Signe_TtBarOld_channel = SignTTbarNew;
            }
            else{
                double SignTTbarNew= Sign_t_tBar(ListTime,channel);
                sign_CprimeC = sign_CprimeC * SignTTbarNew/Signe_TtBarOld_channelBar;
                Signe_TtBarOld_channelBar = SignTTbarNew;
            }

            Traceimpurity = Traceimpurity * TrImp;
            Complex WC = TrImp * RationDet * sign_CprimeC;
            
            ArgumentProba =  ArgumentProba * WC/abs(WC);
            if (std::isnan(real(ArgumentProba)) or std::isnan(imag(ArgumentProba))){
                cout<<"Adding Error: NaN detected" <<endl;
                cout<< "Trace Impurity =" << TrImp <<endl;
                cout<< "Signe C C prime =" << sign_CprimeC <<endl;
                cout<< "Ratio Determinant =" << RationDet <<endl;
            }
        }
    }
}






/**
 * @brief Performs the "Removing" update in the Diagrammatic Monte Carlo simulation.
 * * This subroutine proposes the removal of an existing vertex pair from the configuration.
 * It calculates the transition probability based on the impurity trace ratio and the
 * determinant ratio of the hybridization matrix, then updates the system state
 * and cumulative signs if the move is accepted.
 *
 * @param ListTime The linked-list structure representing the worldline configuration.
 * @param ParametresMC Vector of MC parameters: [t_max, U, EpsilonI, gammaUp, gammaDown].
 * @param ParametresBath Parameters defining the bath hybridization function.
 * @param InvMatrix_sig Inverse hybridization matrix for channel 0.
 * @param InvMatrix_sigBar Inverse hybridization matrix for channel 1.
 * @param k_accep Counter for the number of accepted moves.
 * @param ArgumentProba Cumulative phase/sign of the Monte Carlo configuration weight.
 * @param Traceimpurity Current cumulative value of the impurity trace.
 * @param SigneBetween_HHtilde_Channel Tracker for the H/Htilde fermionic sign (channel 0).
 * @param SigneBetween_HHtilde_ChannelBar Tracker for the H/Htilde fermionic sign (channel 1).
 * @param Signe_TtBarOld_channel Tracker for the time-ordering sign (channel 0).
 * @param Signe_TtBarOld_channelBar Tracker for the time-ordering sign (channel 1).
 * @param EffectifHamltonian Toggle for calculation mode (0: with jumps, 1: effective Hamiltonian).
 */
void MonteCarlo_Removing(clist & ListTime, vector<double> ParametresMC , vector<double> & ParametresBath ,Mat<Complex> & InvMatrix_sig ,Mat<Complex> & InvMatrix_sigBar ,  int & k_accep ,Complex & ArgumentProba , Complex & Traceimpurity,double & SigneBetween_HHtilde_Channel , double & SigneBetween_HHtilde_ChannelBar ,double & Signe_TtBarOld_channel , double & Signe_TtBarOld_channelBar , int &  EffectifHamltonian){
    
    complex<double> Complexi(0.0,1.0);
    std::random_device rd;
    std::default_random_engine eng(rd());

    /* --- Local node pointers for vertex identification --- */
    node * Ope_removeEarly=NULL;
    node * Ope_removeLatter=NULL;
    node * prev_alphaEarly = NULL;
    node * prev_alphaLatter=NULL;
    
    /* --- Hybridization weight components --- */
    std::vector<double> l_sig= {0. , 0.};
    std::vector<double> W_sig= {0. , 0.};
    std::vector<double> O_sig= {0. , 0.};
    Complex ProbabilityAccep = 1.;

    /* --- 1. Propose removal of a vertex pair --- */
    ListTime.propose_removing(Ope_removeEarly, Ope_removeLatter, prev_alphaEarly);
    
    /* --- 2. Calculate the Impurity Trace Ratio --- */
    ListTime.Trace_Removing(ParametresMC[0], Ope_removeEarly, Ope_removeLatter, prev_alphaEarly, prev_alphaLatter, l_sig, W_sig, O_sig);
    Complex TrImp = Trace_impurity(l_sig, W_sig, O_sig, ParametresMC[1], ParametresMC[2], ParametresMC[3], ParametresMC[4], EffectifHamltonian);
    ProbabilityAccep = ProbabilityAccep * TrImp;
    
    /* --- 3. Calculate Transition Probability --- */
    vector<vector<int>> NombreOpe{{ListTime.nb_operator[0][0], ListTime.nb_operator[0][1]},{ListTime.nb_operator[1][0],ListTime.nb_operator[1][1]}};
    int move = 1; // Identifying "Removing" move type
    if (Ope_removeLatter->data.alpha == Ope_removeEarly->data.alpha ){
        // Case: Removing a diagonal vertex
        ProbabilityAccep = ProbabilityAccep * ListTime.Probability_move(move, Ope_removeEarly->data, Ope_removeLatter->data, Ope_removeLatter->fwd_alpha, prev_alphaEarly, Ope_removeEarly->data.channel, ParametresMC[0], NombreOpe);
    }
    else{
        // Case: Removing an off-diagonal vertex
        ProbabilityAccep = ProbabilityAccep * ListTime.Probability_move(move, Ope_removeEarly->data, Ope_removeLatter->data, prev_alphaLatter, prev_alphaEarly, Ope_removeEarly->data.channel, ParametresMC[0], NombreOpe);
    }
    
    /* --- 4. Identify Matrix Indices and Compute Determinant Ratio --- */
    node * previousTime_Early = NULL;
    node * previousTime_Latter = NULL;
    int IndexEarly = 0;
    int IndexLatter = 0;
    ListTime.RecherchePreviousTime_and_index(Ope_removeEarly->data, IndexEarly, previousTime_Early);
    ListTime.RecherchePreviousTime_and_index(Ope_removeLatter->data, IndexLatter, previousTime_Latter);

    Complex RatioDet = 0.;
    if (Ope_removeEarly->data.channel == 0){
        if (Ope_removeEarly->data.flag == 0){
            RatioDet = RatioDet_removing(IndexEarly, IndexLatter, InvMatrix_sig);
        }
        else{
            RatioDet = RatioDet_removing(IndexLatter, IndexEarly, InvMatrix_sig);
        }
    }
    else{
        if (Ope_removeEarly->data.flag == 0){
            RatioDet = RatioDet_removing(IndexEarly, IndexLatter, InvMatrix_sigBar);
        }
        else{
            RatioDet = RatioDet_removing(IndexLatter, IndexEarly, InvMatrix_sigBar);
        }
    }
    ProbabilityAccep = ProbabilityAccep * RatioDet;

    /* --- 5. Acceptance Decision --- */
    std::uniform_real_distribution<double> distr_dag(0,1);
    double RandMC = distr_dag(eng);
    
    if (RandMC < min(1., abs(ProbabilityAccep))){
        k_accep = k_accep + 1;

        /* --- 6. Update Inverse Hybridization Matrix --- */
        if (Ope_removeEarly->data.channel == 0){
            if (Ope_removeEarly->data.flag == 0){
                FastUpdate_removing(IndexEarly, IndexLatter, InvMatrix_sig);
            }
            else{
                FastUpdate_removing(IndexLatter, IndexEarly, InvMatrix_sig);
            }
        }
        else{
            if (Ope_removeEarly->data.flag == 0){
                FastUpdate_removing(IndexEarly, IndexLatter, InvMatrix_sigBar);
            }
            else{
                FastUpdate_removing(IndexLatter, IndexEarly, InvMatrix_sigBar);
            }
        }

        /* --- 7. Global Sign and Phase Tracking --- */
        Complex sign_CprimeC = Complexi;
        sign_CprimeC = sign_CprimeC * ListTime.Sign_TimeOrdering(Ope_removeEarly, Ope_removeLatter);
        
        int ancienAlphaTail = 2; // Sentinel value
        if (Ope_removeEarly->data.alpha == Ope_removeLatter->data.alpha){
            if (ListTime.nb_operator[0][Ope_removeEarly->data.channel] % 2 != 0){
                if (ListTime.tail[0][Ope_removeEarly->data.channel]->data.time > ListTime.tail[1][Ope_removeEarly->data.channel]->data.time ){
                    ancienAlphaTail = 0;
                }
                if (ListTime.tail[1][Ope_removeEarly->data.channel]->data.time > ListTime.tail[0][Ope_removeEarly->data.channel]->data.time ){
                    ancienAlphaTail = 1;
                }
            }
        }

        int Diagonalvertex = 0;
        int alphaOperateurRemoveEarly = Ope_removeEarly->data.alpha;
        if (Ope_removeEarly->data.alpha != Ope_removeLatter->data.alpha){
            Diagonalvertex = 1;
        }
        
        int channel = Ope_removeEarly->data.channel;

        /* --- 8. Physical removal of nodes from the linked list --- */
        ListTime.Move_Removing(Ope_removeEarly, Ope_removeLatter, prev_alphaEarly, prev_alphaLatter, previousTime_Early, previousTime_Latter);
        
        /* --- 9. Sign Adjustment for Off-Diagonal Configurations --- */
        if (Diagonalvertex == 1){
            if (ListTime.nb_operator[0][channel] % 2 != 0){
                sign_CprimeC = sign_CprimeC * (Complexi);
            }
            else{
                sign_CprimeC = sign_CprimeC * (-Complexi);
            }
        }

        if (ListTime.nb_operator[0][channel] % 2 != 0){
            if (ListTime.tail[0][channel]->data.time > ListTime.tail[1][channel]->data.time ){
                if (ancienAlphaTail == 1) sign_CprimeC = sign_CprimeC * (1.);
            }
            if (ListTime.tail[1][channel]->data.time > ListTime.tail[0][channel]->data.time ){
                if (ancienAlphaTail == 0) sign_CprimeC = sign_CprimeC * (1.);
            }
        }

        /* --- 10. Evolution Operator and Time-Ordering Signs --- */
        if (channel == 0){
            double SigneHHtilde_New = Sign_HHtilde(ListTime, channel);
            sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_Channel);
            SigneBetween_HHtilde_Channel = SigneHHtilde_New;
            SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
        }
        else{
            double SigneHHtilde_New = Sign_HHtilde(ListTime, channel);
            sign_CprimeC = sign_CprimeC * (SigneHHtilde_New / SigneBetween_HHtilde_ChannelBar);
            SigneBetween_HHtilde_ChannelBar = SigneHHtilde_New;
            SigneBetween_HHtilde_Channel = SigneHHtilde_New;
        }

        if (channel == 0){
            double SignTTbarNew = Sign_t_tBar(ListTime, channel);
            sign_CprimeC = sign_CprimeC * SignTTbarNew / Signe_TtBarOld_channel;
            Signe_TtBarOld_channel = SignTTbarNew;
        }
        else{
            double SignTTbarNew = Sign_t_tBar(ListTime, channel);
            sign_CprimeC = sign_CprimeC * SignTTbarNew / Signe_TtBarOld_channelBar;
            Signe_TtBarOld_channelBar = SignTTbarNew;
        }

        /* --- 11. Weight Update and Error Checking --- */
        Traceimpurity = Traceimpurity * TrImp;
        Complex WC = TrImp * RatioDet * sign_CprimeC;
        ArgumentProba = ArgumentProba * WC / abs(WC);

        if (std::isnan(real(ArgumentProba)) || std::isnan(imag(ArgumentProba))){
            cout << "Removing Move: NaN detected" << endl;
            cout << "Trace Impurity =" << TrImp << endl;
            cout << "Signe C/C' =" << sign_CprimeC << endl;
            cout << "Ratio Determinant =" << RatioDet << endl;
        }
    }
}


/**
 * @brief Performs the "Shifting" update in the Diagrammatic Monte Carlo simulation.
 * * This move proposes changing the time of an existing operator. It requires a
 * Rank-1 update of the inverse hybridization matrix (either a row or a column
 * depending on the operator's flag) and tracks fermionic signs resulting from
 * crossing other operators in the time-ordered sequence.
 *
 * @param ListTime The worldline configuration linked-list.
 * @param ParametresMC Vector of MC parameters: [t_max, U, EpsilonI, gammaUp, gammaDown].
 * @param ParametresBath Bath hybridization function parameters.
 * @param InvMatrix_sig Inverse hybridization matrix for channel 0.
 * @param InvMatrix_sigBar Inverse hybridization matrix for channel 1.
 * @param k_accep Counter for accepted moves.
 * @param ArgumentProba Cumulative phase/sign of the configuration.
 * @param Traceimpurity Cumulative impurity trace value.
 * @param SigneBetween_HHtilde_Channel Fermionic sign tracker (channel 0).
 * @param SigneBetween_HHtilde_ChannelBar Fermionic sign tracker (channel 1).
 * @param Signe_TtBarOld_channel Time-ordering sign tracker (channel 0).
 * @param Signe_TtBarOld_channelBar Time-ordering sign tracker (channel 1).
 * @param EffectifHamltonian Effective Hamiltonian calculation toggle.
 */
void MonteCarlo_Shifting(clist & ListTime, vector<double> ParametresMC , vector<double> & ParametresBath ,Mat<Complex> & InvMatrix_sig ,Mat<Complex> & InvMatrix_sigBar ,  int & k_accep , Complex & ArgumentProba , Complex & Traceimpurity, double & SigneBetween_HHtilde_Channel , double & SigneBetween_HHtilde_ChannelBar ,double & Signe_TtBarOld_channel , double & Signe_TtBarOld_channelBar , int & EffectifHamltonian  ){
    
    std::random_device rd;
    std::default_random_engine eng(rd());
    complex<double> Complexi(0.0,1.0);

    node * NodeShift = NULL;
    node * prev_alpha = NULL;
    double New_time = 0;
    
    std::vector<double> l_sig= {0. , 0.};
    std::vector<double> W_sig= {0. , 0.};
    std::vector<double> O_sig= {0. , 0.};
    Complex ProbabilityAccep = 1.;
    
    /* --- 1. Propose the Shift move --- */
    ListTime.propose_shift(NodeShift, prev_alpha, New_time, ParametresMC[0]);
    
    /* --- 2. Calculate Impurity Trace Ratio --- */
    ListTime.Trace_Shifting(ParametresMC[0], NodeShift, prev_alpha, New_time, l_sig, W_sig, O_sig);
    Complex TrImp = Trace_impurity(l_sig, W_sig, O_sig, ParametresMC[1], ParametresMC[2], ParametresMC[3], ParametresMC[4], EffectifHamltonian);
    ProbabilityAccep = ProbabilityAccep * TrImp;
    
    /* --- 3. Count Swaps for Fermionic Signs --- */
    int NombreSwap = 0;
    int NombreOpeCHannelBar = 0;
    if (New_time < NodeShift->data.time){
        node * readNode = NodeShift;
        if (readNode->p_prev != NULL){
            while((readNode->p_prev)->data.time > New_time){
                readNode = readNode->p_prev;
                if ((readNode->data.channel == NodeShift->data.channel) && (readNode->data.flag == NodeShift->data.flag)) ++NombreSwap;
                if (readNode->data.channel != NodeShift->data.channel) ++NombreOpeCHannelBar;
                if (readNode->p_prev == NULL) break;
            }
        }
    }
    else {
        node * readNode = NodeShift;
        if (readNode->p_next != NULL){
            while((readNode->p_next)->data.time < New_time){
                readNode = readNode->p_next;
                if ((readNode->data.channel == NodeShift->data.channel) && (readNode->data.flag == NodeShift->data.flag)) ++NombreSwap;
                if (readNode->data.channel != NodeShift->data.channel) ++NombreOpeCHannelBar;
                if (readNode->p_next == NULL) break;
            }
        }
    }

    Complex RatioDeterminant = 0.0;
    
    /* --- 4. Branching Logic: Row vs Column Shift --- */
    if (NodeShift->data.flag == 0){
        /* CASE: Shift acting on a Matrix Row */
        int IndexShift = 0;
        node * previous_time = NULL;
        ListTime.RecherchePreviousTime_and_index(NodeShift->data, IndexShift, previous_time);
        
        int matrix_dim = (ListTime.nb_operator[0][NodeShift->data.channel] + ListTime.nb_operator[1][NodeShift->data.channel])/2;
        Row<Complex> NewRow(matrix_dim);
        Row<Complex> OldRow(matrix_dim);

        node * readNode = ListTime.time_head[abs(NodeShift->data.flag - 1)][NodeShift->data.channel];
        int k = 0;
        while (readNode != NULL){
            NewRow[k] = Hybridization_Funct(NodeShift->data.alpha, readNode->data.alpha, New_time, readNode->data.time, ParametresBath);
            OldRow[k] = Hybridization_Funct(NodeShift->data.alpha, readNode->data.alpha, NodeShift->data.time, readNode->data.time, ParametresBath);
            k++;
            readNode = readNode->fwd_time;
        }
        NewRow = NewRow - OldRow; // Delta row for Sherman-Morrison
        
        Mat<Complex> &targetInv = (NodeShift->data.channel == 0) ? InvMatrix_sig : InvMatrix_sigBar;
        RatioDeterminant = RatioDet_Shifting_Line(targetInv, NewRow, IndexShift) * pow(-1., NombreSwap);
        ProbabilityAccep *= RatioDeterminant;

        std::uniform_real_distribution<double> distr_dag(0,1);
        if (distr_dag(eng) < min(1., abs(ProbabilityAccep))){
            k_accep++;
            FastUpdate_shifting_line(targetInv, NewRow, IndexShift, RatioDeterminant * pow(-1., NombreSwap));
            
            // Re-order matrix to maintain time-ordering sync
            if (New_time < NodeShift->data.time)
                for (int i=0; i<NombreSwap; ++i) targetInv.swap_cols(IndexShift - i, IndexShift - i - 1);
            else
                for (int i=0; i<NombreSwap; ++i) targetInv.swap_cols(IndexShift + i, IndexShift + i + 1);

            // Update configuration and signs
            updateConfiguration(ListTime, NodeShift, prev_alpha, New_time, TrImp, RatioDeterminant, NombreOpeCHannelBar, ArgumentProba, Traceimpurity, SigneBetween_HHtilde_Channel, SigneBetween_HHtilde_ChannelBar, Signe_TtBarOld_channel, Signe_TtBarOld_channelBar);
        }
    }
    else {
        /* CASE: Shift acting on a Matrix Column */
        int IndexShift = 0;
        node * previous_time = NULL;
        ListTime.RecherchePreviousTime_and_index(NodeShift->data, IndexShift, previous_time);

        int matrix_dim = (ListTime.nb_operator[0][NodeShift->data.channel] + ListTime.nb_operator[1][NodeShift->data.channel])/2;
        Col<Complex> NewCol(matrix_dim);
        Col<Complex> OldCol(matrix_dim);
        
        int k = 0;
        node * readNode = ListTime.time_head[abs(NodeShift->data.flag - 1)][NodeShift->data.channel];
        while(readNode != NULL){
            NewCol[k] = Hybridization_Funct(readNode->data.alpha, NodeShift->data.alpha, readNode->data.time, New_time, ParametresBath);
            OldCol[k] = Hybridization_Funct(readNode->data.alpha, NodeShift->data.alpha, readNode->data.time, NodeShift->data.time, ParametresBath);
            k++;
            readNode = readNode->fwd_time;
        }
        NewCol = NewCol - OldCol;

        Mat<Complex> &targetInv = (NodeShift->data.channel == 0) ? InvMatrix_sig : InvMatrix_sigBar;
        RatioDeterminant = RatioDet_Shifting_Col(targetInv, NewCol, IndexShift) * pow(-1., NombreSwap);
        ProbabilityAccep *= RatioDeterminant;

        std::uniform_real_distribution<double> distr_dag(0,1);
        if (distr_dag(eng) < min(1., abs(ProbabilityAccep))){
            k_accep++;
            FastUpdate_shifting_Col(targetInv, NewCol, IndexShift, RatioDeterminant * pow(-1., NombreSwap));
            
            if (New_time < NodeShift->data.time)
                for (int i=0; i<NombreSwap; ++i) targetInv.swap_rows(IndexShift - i, IndexShift - i - 1);
            else
                for (int i=0; i<NombreSwap; ++i) targetInv.swap_rows(IndexShift + i, IndexShift + i + 1);

            updateConfiguration(ListTime, NodeShift, prev_alpha, New_time, TrImp, RatioDeterminant, NombreOpeCHannelBar, ArgumentProba, Traceimpurity, SigneBetween_HHtilde_Channel, SigneBetween_HHtilde_ChannelBar, Signe_TtBarOld_channel, Signe_TtBarOld_channelBar);
        }
    }
}







/**
 * @brief Calculates the sign factor when commuting all H-tilde operators to the right.
 */
double Sign_HHtilde2(clist & ListTime, int channel) {
    double SignConfig = 1.0;
    node * ReadNode = ListTime.head;
    int NombreOpeH = 0;

    while (ReadNode != NULL) {
        if (ReadNode->data.alpha == 0) {
            NombreOpeH++;
        } else if (ReadNode->data.alpha == 1) {
            // If we cross an odd number of H operators, sign flips
            if (NombreOpeH % 2 != 0) SignConfig *= -1.0;
        }
        ReadNode = ReadNode->p_next;
    }
    return SignConfig;
}

/**
 * @brief Determines the sign based on the sequence of forward (t) and backward (tBar) operators.
 */
double Sign_t_tBar2(clist & ListTime, int channel) {
    double Signe_Between_ttBar = 1.0;
    int NombreDoublon = 0;
    node * readNode = ListTime.head;

    if (readNode != NULL) {
        while (readNode != NULL && readNode->p_next != NULL) {
            if (readNode->data.flag == 0) {
                if (readNode->p_next->data.flag == 1) Signe_Between_ttBar *= -1.0;
                else NombreDoublon++;
            } else if (readNode->data.flag == 1) {
                if (readNode->p_next->data.flag == 1) NombreDoublon++;
            }
            readNode = readNode->p_next;
            if (readNode != NULL) readNode = readNode->p_next;
        }
    }

    if (NombreDoublon % 2 != 0) Signe_Between_ttBar *= -1.0;
    
    return Signe_Between_ttBar;
}
