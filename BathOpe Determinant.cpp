/**
 * @file BathOpe Determinant.cpp
 * @brief Implementation of determinant-based updates for Diagrammatic Monte Carlo (DMC).
 * @author Matthieu Vanhoecke
 * @date 2023-03-27
 */

#include "BathOpe Determinant.hpp"
#include <complex>
#include <tuple>
#include <vector>

/*----------------------------------------------------------------------------------------------------------
                                 Matrix Component Generation (Move Adding)
 ----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Generates the column and row vectors required for adding a new operator pair.
 * * This subroutine constructs the hybridization vectors by iterating through the
 *   linked list of existing operators.
 * @param List_time The linked list containing operator times and channels.
 * @param OperatorEarly The operator appearing earlier in the sequence.
 * @param OperatorLatter The operator appearing later in the sequence.
 * @param index_col Output: index of the column to be added.
 * @param index_line Output: index of the line to be added.
 * @param ParametreBath Physics parameters for the hybridization function.
 * @return std::tuple containing the new Column and Row vectors.
 */
std::tuple<arma::Col<std::complex<double>>, arma::Row<std::complex<double>>>
Creation_ColumnRow(clist & List_time, kink & OperatorEarly, kink & OperatorLatter,
                   int & index_col, int & index_line, std::vector<double> ParametreBath) {
    
    int sizeVecRow = (List_time.nb_operator[0][OperatorEarly.channel] +
                      List_time.nb_operator[1][OperatorEarly.channel]) / 2;
    
    arma::Col<std::complex<double>> ColAdd(sizeVecRow);
    arma::Row<std::complex<double>> RowAdd(sizeVecRow);

    kink Ope_t;
    kink Ope_tBar;
    
    // Assign t and tbar based on the flags (flag 0 is t, flag 1 is tbar)
    if (OperatorEarly.flag == 0) {
        Ope_t = OperatorEarly;
        Ope_tBar = OperatorLatter;
    } else {
        Ope_t = OperatorLatter;
        Ope_tBar = OperatorEarly;
    }

    // --- Row Construction (Fixing t, varying tbar) ---
    node * readNode = List_time.time_head[abs(Ope_t.flag - 1)][Ope_t.channel];
    index_col = 0;
    double tNew = Ope_t.time;
    int alpha_OpeEarl = Ope_t.alpha;
    int i = 0;

    while (readNode != nullptr) {
        if (readNode->data.time < Ope_tBar.time) {
            index_col++;
        }
        RowAdd[i] = Hybridization_Funct(alpha_OpeEarl, readNode->data.alpha, tNew, readNode->data.time, ParametreBath);
        i++;
        readNode = readNode->fwd_time;
    }

    // --- Column Construction (Fixing tbar, varying t) ---
    readNode = List_time.time_head[abs(Ope_t.flag)][Ope_t.channel];
    index_line = 0;
    tNew = Ope_tBar.time;
    alpha_OpeEarl = Ope_tBar.alpha;
    i = 0;

    while (readNode != nullptr) {
        if (readNode->data.time < Ope_t.time) {
            index_line++;
        }
        ColAdd[i] = Hybridization_Funct(readNode->data.alpha, alpha_OpeEarl, readNode->data.time, tNew, ParametreBath);
        i++;
        readNode = readNode->fwd_time;
    }

    return std::make_tuple(ColAdd, RowAdd);
}

/*----------------------------------------------------------------------------------------------------------
                                   Determinant Ratio Evaluation
 ----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Calculates the ratio of determinants when adding a row and column.
 */
std::complex<double> RatioDet_adding(int & Index, int & IndexBar, arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                     arma::Col<std::complex<double>> New_column, arma::Row<std::complex<double>> New_row,
                                     std::complex<double> DiagElem) {
    double sign = ((Index + IndexBar) % 2 == 0) ? 1.0 : -1.0;
    return sign * (DiagElem - arma::as_scalar(New_row * MatrixOld_sigma * New_column));
}

/**
 * @brief Calculates the ratio of determinants when removing a row and column.
 */
std::complex<double> RatioDet_removing(int const& IndexLine, int const& IndexCol, arma::Mat<std::complex<double>> const& MatrixOld_sigma) {
    double sign = ((IndexLine + IndexCol) % 2 == 0) ? 1.0 : -1.0;
    return sign * arma::as_scalar(MatrixOld_sigma(IndexCol, IndexLine));
}

/**
 * @brief Calculates the ratio of determinants when shifting a line.
 */
std::complex<double> RatioDet_Shifting_Line(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                           arma::Row<std::complex<double>> & Delta_Line, int & IndexLine) {
    arma::Row<std::complex<double>> Product = Delta_Line * MatrixOld_sigma;
    
    // Debug check for NaN results
    if (std::isnan(std::real(Product[IndexLine]))) {
        std::cerr << "Error: NaN detected in RatioDet_Shifting_Line" << std::endl;
    }
    
    Delta_Line = Product; // Update Delta_Line for subsequent FastUpdate
    return (1.0 + arma::as_scalar(Product[IndexLine]));
}

/**
 * @brief Calculates the ratio of determinants when shifting a column.
 */
std::complex<double> RatioDet_Shifting_Col(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                          arma::Col<std::complex<double>> & Delta_column, int & IndexCol) {
    arma::Col<std::complex<double>> Product = MatrixOld_sigma * Delta_column;
    
    if (std::isnan(std::real(Product[IndexCol]))) {
        std::cerr << "Error: NaN detected in RatioDet_Shifting_Col" << std::endl;
    }
    
    Delta_column = Product; // Update Delta_column for subsequent FastUpdate
    return (1.0 + arma::as_scalar(Product[IndexCol]));
}

/*----------------------------------------------------------------------------------------------------------
                                     Fast Matrix Updates (Sherman-Morrison)
 ----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Performs a fast rank-1 update of the inverse hybridization matrix after adding an operator.
 */
void FastUpdate_adding(int const& Index_line, int const& Index_Col, arma::Mat<std::complex<double>> & MatrixOld_sigma,
                       arma::Col<std::complex<double>> New_col, arma::Row<std::complex<double>> New_row,
                       std::complex<double> & Ratio_det) {
    
    arma::Col<std::complex<double>> L = MatrixOld_sigma * New_col;
    arma::Row<std::complex<double>> R = New_row * MatrixOld_sigma;
    
    // Expand vectors to include the new row/column contribution
    R.insert_cols(Index_line, 1);
    R[Index_line] = -1.0;
    L.insert_rows(Index_Col, 1);
    L[Index_Col] = -1.0;

    MatrixOld_sigma.insert_cols(Index_line, 1);
    MatrixOld_sigma.insert_rows(Index_Col, 1);
    
    // Calculate the update factor xi
    double sign = ((Index_line + Index_Col) % 2 == 0) ? 1.0 : -1.0;
    std::complex<double> numerator = sign * Ratio_det;
    std::complex<double> xi = std::conj(numerator) / std::norm(numerator);
    
    MatrixOld_sigma += xi * (L * R);
}

/**
 * @brief Performs a fast update after removing an operator pair.
 */
void FastUpdate_removing(int const& Index_line, int const& Index_Col, arma::Mat<std::complex<double>> & MatrixOld_sigma) {
    MatrixOld_sigma -= (MatrixOld_sigma.col(Index_line) * MatrixOld_sigma.row(Index_Col)) /
                       arma::as_scalar(MatrixOld_sigma(Index_Col, Index_line));
    
    MatrixOld_sigma.shed_row(Index_Col);
    MatrixOld_sigma.shed_col(Index_line);
}

/**
 * @brief Fast update for shifting a column.
 */
void FastUpdate_shifting_Col(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                             arma::Col<std::complex<double>> & Prod_MoldDeltaB,
                             int const& IndexCol, std::complex<double> & Ratio_det) {
    MatrixOld_sigma -= (1.0 / Ratio_det) * Prod_MoldDeltaB * MatrixOld_sigma.row(IndexCol);
}

/**
 * @brief Fast update for shifting a line.
 */
void FastUpdate_shifting_line(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                              arma::Row<std::complex<double>> & Prod_DeltaBMold,
                              int const& IndexCol, std::complex<double> & Ratio_det) {
    MatrixOld_sigma -= (1.0 / Ratio_det) * MatrixOld_sigma.col(IndexCol) * Prod_DeltaBMold;
}
