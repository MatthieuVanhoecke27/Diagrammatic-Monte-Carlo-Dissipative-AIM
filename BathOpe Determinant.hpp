/**
 * @file BathOpe Determinant.hpp
 * @brief Declarations of functions for determinant ratio evaluations and fast matrix updates.
 * * This header provides the interface for handling the hybridization matrix in
 * Diagrammatic Monte Carlo simulations using the Sherman-Morrison update formula.
 * * @author Matthieu Vanhoecke
 * @date 2023-01-02
 */

#ifndef BathOpe_Determinant_hpp
#define BathOpe_Determinant_hpp

#include <complex>
#include <iostream>
#include <vector>
#include <tuple>
#include <armadillo>

#include "FonctionPrincipale.hpp"
#include "HybridizationFunction.hpp"

/*----------------------------------------------------------------------------------------------------------
                                  Matrix Generation & Determinant Ratios
 ----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Constructs the vectors required for expanding the hybridization matrix.
 * * @param List_time Reference to the linked list of operators.
 * @param OperatorEarly The first operator of the pair.
 * @param OperatorLatter The second operator of the pair.
 * @param index_col Reference to store the calculated column index.
 * @param index_line Reference to store the calculated line index.
 * @param ParametreBath Physics parameters for the bath.
 * @return std::tuple containing the (Column, Row) vectors.
 */
std::tuple<arma::Col<std::complex<double>>, arma::Row<std::complex<double>>>
Creation_ColumnRow(clist & List_time, kink & OperatorEarly, kink & OperatorLatter,
                   int & index_col, int & index_line, std::vector<double> ParametreBath);

/**
 * @brief Computes the determinant ratio for adding a new operator pair.
 */
std::complex<double> RatioDet_adding(int & Index, int & IndexBar,
                                     arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                     arma::Col<std::complex<double>> New_column,
                                     arma::Row<std::complex<double>> New_row,
                                     std::complex<double> DiagElem);

/**
 * @brief Computes the determinant ratio for removing an operator pair.
 */
std::complex<double> RatioDet_removing(int const& IndexLine, int const& IndexCol,
                                       arma::Mat<std::complex<double>> const& MatrixOld_sigma);

/**
 * @brief Computes the determinant ratio for shifting an operator's time (Row-wise).
 */
std::complex<double> RatioDet_Shifting_Line(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                           arma::Row<std::complex<double>> & Delta_Line,
                                           int & IndexLine);

/**
 * @brief Computes the determinant ratio for shifting an operator's time (Column-wise).
 */
std::complex<double> RatioDet_Shifting_Col(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                                          arma::Col<std::complex<double>> & Delta_column,
                                          int & IndexCol);

/*----------------------------------------------------------------------------------------------------------
                                     Fast Matrix Update Operations
 ----------------------------------------------------------------------------------------------------------*/

/**
 * @brief Updates the inverse matrix after an operator addition using rank-1 update.
 */
void FastUpdate_adding(int const& Index_line, int const& Index_Col,
                       arma::Mat<std::complex<double>> & MatrixOld_sigma,
                       arma::Col<std::complex<double>> New_col,
                       arma::Row<std::complex<double>> New_row,
                       std::complex<double> & Ratio_det);

/**
 * @brief Updates the inverse matrix after an operator removal.
 */
void FastUpdate_removing(int const& Index_line, int const& Index_Col,
                         arma::Mat<std::complex<double>> & MatrixOld_sigma);

/**
 * @brief Updates the inverse matrix after a column shift.
 * @param Prod_MoldDeltaB Precomputed product from the ratio calculation.
 */
void FastUpdate_shifting_Col(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                             arma::Col<std::complex<double>> & Prod_MoldDeltaB,
                             int const& IndexCol,
                             std::complex<double> & Ratio_det);

/**
 * @brief Updates the inverse matrix after a line (row) shift.
 * @param Prod_DeltaBMold Precomputed product from the ratio calculation.
 */
void FastUpdate_shifting_line(arma::Mat<std::complex<double>> & MatrixOld_sigma,
                              arma::Row<std::complex<double>> & Prod_DeltaBMold,
                              int const& IndexCol,
                              std::complex<double> & Ratio_det);

#endif /* BathOpe_Determinant_hpp */

