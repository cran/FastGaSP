// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/FastGaSP.h"
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// Construct_W0_matern_5_2
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda);
RcppExport SEXP _FastGaSP_Construct_W0_matern_5_2(SEXP sigma2SEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_W0_matern_5_2(sigma2, lambda));
    return rcpp_result_gen;
END_RCPP
}
// Construct_W0_exp
MatrixXd Construct_W0_exp(const double sigma2, const double lambda);
RcppExport SEXP _FastGaSP_Construct_W0_exp(SEXP sigma2SEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_W0_exp(sigma2, lambda));
    return rcpp_result_gen;
END_RCPP
}
// Construct_G_matern_5_2
List Construct_G_matern_5_2(Eigen::VectorXd delta_x, double lambda);
RcppExport SEXP _FastGaSP_Construct_G_matern_5_2(SEXP delta_xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_G_matern_5_2(delta_x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// Construct_G_exp
List Construct_G_exp(Eigen::VectorXd delta_x, double lambda);
RcppExport SEXP _FastGaSP_Construct_G_exp(SEXP delta_xSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_G_exp(delta_x, lambda));
    return rcpp_result_gen;
END_RCPP
}
// Construct_W_matern_5_2
List Construct_W_matern_5_2(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0);
RcppExport SEXP _FastGaSP_Construct_W_matern_5_2(SEXP sigma2SEXP, SEXP delta_xSEXP, SEXP lambdaSEXP, SEXP W0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type W0(W0SEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_W_matern_5_2(sigma2, delta_x, lambda, W0));
    return rcpp_result_gen;
END_RCPP
}
// Construct_W_exp
List Construct_W_exp(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0);
RcppExport SEXP _FastGaSP_Construct_W_exp(SEXP sigma2SEXP, SEXP delta_xSEXP, SEXP lambdaSEXP, SEXP W0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type sigma2(sigma2SEXP);
    Rcpp::traits::input_parameter< Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< MatrixXd >::type W0(W0SEXP);
    rcpp_result_gen = Rcpp::wrap(Construct_W_exp(sigma2, delta_x, lambda, W0));
    return rcpp_result_gen;
END_RCPP
}
// Get_Q_K
List Get_Q_K(const List GG, const List W, const Eigen::MatrixXd C0, const double VV);
RcppExport SEXP _FastGaSP_Get_Q_K(SEXP GGSEXP, SEXP WSEXP, SEXP C0SEXP, SEXP VVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const List >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< const double >::type VV(VVSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_Q_K(GG, W, C0, VV));
    return rcpp_result_gen;
END_RCPP
}
// Get_log_det_S2
List Get_log_det_S2(const Eigen::VectorXd param, const bool have_noise, const Eigen::VectorXd delta_x, const Eigen::VectorXd output, const String kernel_type);
RcppExport SEXP _FastGaSP_Get_log_det_S2(SEXP paramSEXP, SEXP have_noiseSEXP, SEXP delta_xSEXP, SEXP outputSEXP, SEXP kernel_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const bool >::type have_noise(have_noiseSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_log_det_S2(param, have_noise, delta_x, output, kernel_type));
    return rcpp_result_gen;
END_RCPP
}
// Get_L_inv_y
VectorXd Get_L_inv_y(const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K, const Eigen::VectorXd output);
RcppExport SEXP _FastGaSP_Get_L_inv_y(SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_L_inv_y(GG, Q, K, output));
    return rcpp_result_gen;
END_RCPP
}
// Get_C_R_K_Q
List Get_C_R_K_Q(const VectorXi index, const List GG, const List W, const Eigen::MatrixXd C0, double VV);
RcppExport SEXP _FastGaSP_Get_C_R_K_Q(SEXP indexSEXP, SEXP GGSEXP, SEXP WSEXP, SEXP C0SEXP, SEXP VVSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const List >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< double >::type VV(VVSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_C_R_K_Q(index, GG, W, C0, VV));
    return rcpp_result_gen;
END_RCPP
}
// Get_m_a_pred
List Get_m_a_pred(const VectorXi index, const Eigen::VectorXd output_vec, const List GG, const Eigen::MatrixXd K);
RcppExport SEXP _FastGaSP_Get_m_a_pred(SEXP indexSEXP, SEXP output_vecSEXP, SEXP GGSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output_vec(output_vecSEXP);
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_m_a_pred(index, output_vec, GG, K));
    return rcpp_result_gen;
END_RCPP
}
// Get_S_KK
List Get_S_KK(const VectorXi index, const List GG, const List C, const List R);
RcppExport SEXP _FastGaSP_Get_S_KK(SEXP indexSEXP, SEXP GGSEXP, SEXP CSEXP, SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXi >::type index(indexSEXP);
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const List >::type C(CSEXP);
    Rcpp::traits::input_parameter< const List >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_S_KK(index, GG, C, R));
    return rcpp_result_gen;
END_RCPP
}
// Get_s_1st
MatrixXd Get_s_1st(const List m, const List a, const List C, const List KK);
RcppExport SEXP _FastGaSP_Get_s_1st(SEXP mSEXP, SEXP aSEXP, SEXP CSEXP, SEXP KKSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type m(mSEXP);
    Rcpp::traits::input_parameter< const List >::type a(aSEXP);
    Rcpp::traits::input_parameter< const List >::type C(CSEXP);
    Rcpp::traits::input_parameter< const List >::type KK(KKSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_s_1st(m, a, C, KK));
    return rcpp_result_gen;
END_RCPP
}
// Kalman_smoother
List Kalman_smoother(const VectorXd param, const bool have_noise, const VectorXi index_obs, const VectorXd delta_x_all, const VectorXd output, const double sigma_2_hat, const String kernel_type);
RcppExport SEXP _FastGaSP_Kalman_smoother(SEXP paramSEXP, SEXP have_noiseSEXP, SEXP index_obsSEXP, SEXP delta_x_allSEXP, SEXP outputSEXP, SEXP sigma_2_hatSEXP, SEXP kernel_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const bool >::type have_noise(have_noiseSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type index_obs(index_obsSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type delta_x_all(delta_x_allSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const double >::type sigma_2_hat(sigma_2_hatSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Kalman_smoother(param, have_noise, index_obs, delta_x_all, output, sigma_2_hat, kernel_type));
    return rcpp_result_gen;
END_RCPP
}
// Sample_KF
MatrixXd Sample_KF(const List GG, const List W, const Eigen::MatrixXd C0, const double VV, const String kernel_type, const int sample_type);
RcppExport SEXP _FastGaSP_Sample_KF(SEXP GGSEXP, SEXP WSEXP, SEXP C0SEXP, SEXP VVSEXP, SEXP kernel_typeSEXP, SEXP sample_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const List >::type W(WSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type C0(C0SEXP);
    Rcpp::traits::input_parameter< const double >::type VV(VVSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const int >::type sample_type(sample_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample_KF(GG, W, C0, VV, kernel_type, sample_type));
    return rcpp_result_gen;
END_RCPP
}
// Sample_KF_post
MatrixXd Sample_KF_post(const VectorXi index_obs, const List C_R_K_Q, const Eigen::MatrixXd W0, const List GG, const List W, const double VV, const VectorXd output, String kernel_type, const int sample_type);
RcppExport SEXP _FastGaSP_Sample_KF_post(SEXP index_obsSEXP, SEXP C_R_K_QSEXP, SEXP W0SEXP, SEXP GGSEXP, SEXP WSEXP, SEXP VVSEXP, SEXP outputSEXP, SEXP kernel_typeSEXP, SEXP sample_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXi >::type index_obs(index_obsSEXP);
    Rcpp::traits::input_parameter< const List >::type C_R_K_Q(C_R_K_QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type W0(W0SEXP);
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const List >::type W(WSEXP);
    Rcpp::traits::input_parameter< const double >::type VV(VVSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const int >::type sample_type(sample_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Sample_KF_post(index_obs, C_R_K_Q, W0, GG, W, VV, output, kernel_type, sample_type));
    return rcpp_result_gen;
END_RCPP
}
// Get_L_t_y
VectorXd Get_L_t_y(const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K, const Eigen::VectorXd output);
RcppExport SEXP _FastGaSP_Get_L_t_y(SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_L_t_y(GG, Q, K, output));
    return rcpp_result_gen;
END_RCPP
}
// Get_L_y
VectorXd Get_L_y(const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K, const Eigen::VectorXd output);
RcppExport SEXP _FastGaSP_Get_L_y(SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_L_y(GG, Q, K, output));
    return rcpp_result_gen;
END_RCPP
}
// Get_L_t_inv_y
VectorXd Get_L_t_inv_y(const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K, const Eigen::VectorXd output);
RcppExport SEXP _FastGaSP_Get_L_t_inv_y(SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_L_t_inv_y(GG, Q, K, output));
    return rcpp_result_gen;
END_RCPP
}
// Get_R_y
VectorXd Get_R_y(const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K, const Eigen::VectorXd output);
RcppExport SEXP _FastGaSP_Get_R_y(SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP, SEXP outputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type output(outputSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_R_y(GG, Q, K, output));
    return rcpp_result_gen;
END_RCPP
}
// Get_Y_minus_a_1_scaled_matrix_2d
Eigen::MatrixXd Get_Y_minus_a_1_scaled_matrix_2d(const Eigen::MatrixXd output_KF, const List GG, const Eigen::VectorXd Q, const Eigen::MatrixXd K);
RcppExport SEXP _FastGaSP_Get_Y_minus_a_1_scaled_matrix_2d(SEXP output_KFSEXP, SEXP GGSEXP, SEXP QSEXP, SEXP KSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type output_KF(output_KFSEXP);
    Rcpp::traits::input_parameter< const List >::type GG(GGSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type Q(QSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type K(KSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_Y_minus_a_1_scaled_matrix_2d(output_KF, GG, Q, K));
    return rcpp_result_gen;
END_RCPP
}
// F_Funct
double F_Funct(const Eigen::MatrixXd A_cur, const List G);
RcppExport SEXP _FastGaSP_F_Funct(SEXP A_curSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(F_Funct(A_cur, G));
    return rcpp_result_gen;
END_RCPP
}
// Get_G_log_det_cov
List Get_G_log_det_cov(const Eigen::VectorXd param, const Eigen::MatrixXd output, const Eigen::VectorXd delta_x, int d, const String kernel_type);
RcppExport SEXP _FastGaSP_Get_G_log_det_cov(SEXP paramSEXP, SEXP outputSEXP, SEXP delta_xSEXP, SEXP dSEXP, SEXP kernel_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type delta_x(delta_xSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_G_log_det_cov(param, output, delta_x, d, kernel_type));
    return rcpp_result_gen;
END_RCPP
}
// KF_cpp_eigen
List KF_cpp_eigen(double F_t, double V_t, double G_t, double W_t, const Eigen::VectorXd& y);
RcppExport SEXP _FastGaSP_KF_cpp_eigen(SEXP F_tSEXP, SEXP V_tSEXP, SEXP G_tSEXP, SEXP W_tSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type F_t(F_tSEXP);
    Rcpp::traits::input_parameter< double >::type V_t(V_tSEXP);
    Rcpp::traits::input_parameter< double >::type G_t(G_tSEXP);
    Rcpp::traits::input_parameter< double >::type W_t(W_tSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(KF_cpp_eigen(F_t, V_t, G_t, W_t, y));
    return rcpp_result_gen;
END_RCPP
}
// cubic_solver
double cubic_solver(const std::vector<double>& p);
RcppExport SEXP _FastGaSP_cubic_solver(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<double>& >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(cubic_solver(p));
    return rcpp_result_gen;
END_RCPP
}
// fmou_cpp
List fmou_cpp(const Eigen::MatrixXd& output, int d, int M, double threshold, bool est_U0, bool est_sigma0_2, bool track_iterations, bool track_neg_log_lik, Nullable<Eigen::MatrixXd> U0, Nullable<Eigen::MatrixXd> U_init, Nullable<Eigen::VectorXd> rho_init, Nullable<Eigen::VectorXd> sigma2_init, Nullable<double> sigma0_2);
RcppExport SEXP _FastGaSP_fmou_cpp(SEXP outputSEXP, SEXP dSEXP, SEXP MSEXP, SEXP thresholdSEXP, SEXP est_U0SEXP, SEXP est_sigma0_2SEXP, SEXP track_iterationsSEXP, SEXP track_neg_log_likSEXP, SEXP U0SEXP, SEXP U_initSEXP, SEXP rho_initSEXP, SEXP sigma2_initSEXP, SEXP sigma0_2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type output(outputSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type est_U0(est_U0SEXP);
    Rcpp::traits::input_parameter< bool >::type est_sigma0_2(est_sigma0_2SEXP);
    Rcpp::traits::input_parameter< bool >::type track_iterations(track_iterationsSEXP);
    Rcpp::traits::input_parameter< bool >::type track_neg_log_lik(track_neg_log_likSEXP);
    Rcpp::traits::input_parameter< Nullable<Eigen::MatrixXd> >::type U0(U0SEXP);
    Rcpp::traits::input_parameter< Nullable<Eigen::MatrixXd> >::type U_init(U_initSEXP);
    Rcpp::traits::input_parameter< Nullable<Eigen::VectorXd> >::type rho_init(rho_initSEXP);
    Rcpp::traits::input_parameter< Nullable<Eigen::VectorXd> >::type sigma2_init(sigma2_initSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type sigma0_2(sigma0_2SEXP);
    rcpp_result_gen = Rcpp::wrap(fmou_cpp(output, d, M, threshold, est_U0, est_sigma0_2, track_iterations, track_neg_log_lik, U0, U_init, rho_init, sigma2_init, sigma0_2));
    return rcpp_result_gen;
END_RCPP
}
// matern_5_2_funct
Eigen::MatrixXd matern_5_2_funct(const Eigen::Map<Eigen::MatrixXd>& d, double beta_i);
RcppExport SEXP _FastGaSP_matern_5_2_funct(SEXP dSEXP, SEXP beta_iSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd>& >::type d(dSEXP);
    Rcpp::traits::input_parameter< double >::type beta_i(beta_iSEXP);
    rcpp_result_gen = Rcpp::wrap(matern_5_2_funct(d, beta_i));
    return rcpp_result_gen;
END_RCPP
}
// rcppeigen_get_chol
MatrixXd rcppeigen_get_chol(const MatrixXd& R);
RcppExport SEXP _FastGaSP_rcppeigen_get_chol(SEXP RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const MatrixXd& >::type R(RSEXP);
    rcpp_result_gen = Rcpp::wrap(rcppeigen_get_chol(R));
    return rcpp_result_gen;
END_RCPP
}
// F_Funct_Dev_Large_k
Eigen::MatrixXd F_Funct_Dev_Large_k(const Eigen::MatrixXd A_cur, const List UD);
RcppExport SEXP _FastGaSP_F_Funct_Dev_Large_k(SEXP A_curSEXP, SEXP UDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type UD(UDSEXP);
    rcpp_result_gen = Rcpp::wrap(F_Funct_Dev_Large_k(A_cur, UD));
    return rcpp_result_gen;
END_RCPP
}
// Get_B_U_V_Large_k
List Get_B_U_V_Large_k(const Eigen::MatrixXd A_cur, const List UD);
RcppExport SEXP _FastGaSP_Get_B_U_V_Large_k(SEXP A_curSEXP, SEXP UDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type UD(UDSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_B_U_V_Large_k(A_cur, UD));
    return rcpp_result_gen;
END_RCPP
}
// Y_Funct
List Y_Funct(const Eigen::MatrixXd A_cur, const List B_U_V, double tau);
RcppExport SEXP _FastGaSP_Y_Funct(SEXP A_curSEXP, SEXP B_U_VSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type B_U_V(B_U_VSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(Y_Funct(A_cur, B_U_V, tau));
    return rcpp_result_gen;
END_RCPP
}
// F_Funct_Dev
Eigen::MatrixXd F_Funct_Dev(const Eigen::MatrixXd A_cur, const List G);
RcppExport SEXP _FastGaSP_F_Funct_Dev(SEXP A_curSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(F_Funct_Dev(A_cur, G));
    return rcpp_result_gen;
END_RCPP
}
// Get_B_U_V
List Get_B_U_V(const Eigen::MatrixXd A_cur, const List G);
RcppExport SEXP _FastGaSP_Get_B_U_V(SEXP A_curSEXP, SEXP GSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_cur(A_curSEXP);
    Rcpp::traits::input_parameter< const List >::type G(GSEXP);
    rcpp_result_gen = Rcpp::wrap(Get_B_U_V(A_cur, G));
    return rcpp_result_gen;
END_RCPP
}
// Optimization_Stiefel_Manifold
Eigen::MatrixXd Optimization_Stiefel_Manifold(const Eigen::MatrixXd A_ini, const List G, int max_iter);
RcppExport SEXP _FastGaSP_Optimization_Stiefel_Manifold(SEXP A_iniSEXP, SEXP GSEXP, SEXP max_iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd >::type A_ini(A_iniSEXP);
    Rcpp::traits::input_parameter< const List >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    rcpp_result_gen = Rcpp::wrap(Optimization_Stiefel_Manifold(A_ini, G, max_iter));
    return rcpp_result_gen;
END_RCPP
}
// A_t_times_x_particle
VectorXd A_t_times_x_particle(const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi num_neighbors_vec, const int D_y, const int N_tilde);
RcppExport SEXP _FastGaSP_A_t_times_x_particle(SEXP outputSEXP, SEXP A_all_vSEXP, SEXP num_neighbors_vecSEXP, SEXP D_ySEXP, SEXP N_tildeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v(A_all_vSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec(num_neighbors_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type D_y(D_ySEXP);
    Rcpp::traits::input_parameter< const int >::type N_tilde(N_tildeSEXP);
    rcpp_result_gen = Rcpp::wrap(A_t_times_x_particle(output, A_all_v, num_neighbors_vec, D_y, N_tilde));
    return rcpp_result_gen;
END_RCPP
}
// A_times_x_particle
VectorXd A_times_x_particle(const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi num_neighbors_vec, const int D, const int N);
RcppExport SEXP _FastGaSP_A_times_x_particle(SEXP outputSEXP, SEXP A_all_vSEXP, SEXP num_neighbors_vecSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v(A_all_vSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec(num_neighbors_vecSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(A_times_x_particle(output, A_all_v, num_neighbors_vec, D, N));
    return rcpp_result_gen;
END_RCPP
}
// IKF_CG_particle
List IKF_CG_particle(VectorXd param, const String kernel_type, const VectorXd delta_x_all, const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi sort_d_all_ix, const VectorXi num_neighbors_vec, const double tilde_nu, const int D, const int N, float tol, int maxIte);
RcppExport SEXP _FastGaSP_IKF_CG_particle(SEXP paramSEXP, SEXP kernel_typeSEXP, SEXP delta_x_allSEXP, SEXP outputSEXP, SEXP A_all_vSEXP, SEXP sort_d_all_ixSEXP, SEXP num_neighbors_vecSEXP, SEXP tilde_nuSEXP, SEXP DSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP maxIteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type delta_x_all(delta_x_allSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v(A_all_vSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type sort_d_all_ix(sort_d_all_ixSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec(num_neighbors_vecSEXP);
    Rcpp::traits::input_parameter< const double >::type tilde_nu(tilde_nuSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIte(maxIteSEXP);
    rcpp_result_gen = Rcpp::wrap(IKF_CG_particle(param, kernel_type, delta_x_all, output, A_all_v, sort_d_all_ix, num_neighbors_vec, tilde_nu, D, N, tol, maxIte));
    return rcpp_result_gen;
END_RCPP
}
// IKF_CG_particle_two_interact
List IKF_CG_particle_two_interact(VectorXd param1, VectorXd param2, const String kernel_type1, const String kernel_type2, const VectorXd delta_x_all1, const VectorXd delta_x_all2, const Eigen::VectorXd A_all_v1, const Eigen::VectorXd A_all_v2, const VectorXi sort_d_all_ix1, const VectorXi sort_d_all_ix2, const VectorXi num_neighbors_vec1, const VectorXi num_neighbors_vec2, const VectorXd output, const double tilde_nu, const int D, const int N, float tol, int maxIte);
RcppExport SEXP _FastGaSP_IKF_CG_particle_two_interact(SEXP param1SEXP, SEXP param2SEXP, SEXP kernel_type1SEXP, SEXP kernel_type2SEXP, SEXP delta_x_all1SEXP, SEXP delta_x_all2SEXP, SEXP A_all_v1SEXP, SEXP A_all_v2SEXP, SEXP sort_d_all_ix1SEXP, SEXP sort_d_all_ix2SEXP, SEXP num_neighbors_vec1SEXP, SEXP num_neighbors_vec2SEXP, SEXP outputSEXP, SEXP tilde_nuSEXP, SEXP DSEXP, SEXP NSEXP, SEXP tolSEXP, SEXP maxIteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type param1(param1SEXP);
    Rcpp::traits::input_parameter< VectorXd >::type param2(param2SEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type1(kernel_type1SEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type2(kernel_type2SEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type delta_x_all1(delta_x_all1SEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type delta_x_all2(delta_x_all2SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v1(A_all_v1SEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v2(A_all_v2SEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type sort_d_all_ix1(sort_d_all_ix1SEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type sort_d_all_ix2(sort_d_all_ix2SEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec1(num_neighbors_vec1SEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec2(num_neighbors_vec2SEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const double >::type tilde_nu(tilde_nuSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIte(maxIteSEXP);
    rcpp_result_gen = Rcpp::wrap(IKF_CG_particle_two_interact(param1, param2, kernel_type1, kernel_type2, delta_x_all1, delta_x_all2, A_all_v1, A_all_v2, sort_d_all_ix1, sort_d_all_ix2, num_neighbors_vec1, num_neighbors_vec2, output, tilde_nu, D, N, tol, maxIte));
    return rcpp_result_gen;
END_RCPP
}
// IKF_CG_particle_cell
List IKF_CG_particle_cell(VectorXd param, const String kernel_type, const VectorXd delta_x_all, const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi sort_d_all_ix, const VectorXd sigma_2_vec, const VectorXi num_neighbors_vec, const double tilde_nu, const int D, const VectorXd n_t_record, float tol, int maxIte);
RcppExport SEXP _FastGaSP_IKF_CG_particle_cell(SEXP paramSEXP, SEXP kernel_typeSEXP, SEXP delta_x_allSEXP, SEXP outputSEXP, SEXP A_all_vSEXP, SEXP sort_d_all_ixSEXP, SEXP sigma_2_vecSEXP, SEXP num_neighbors_vecSEXP, SEXP tilde_nuSEXP, SEXP DSEXP, SEXP n_t_recordSEXP, SEXP tolSEXP, SEXP maxIteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< VectorXd >::type param(paramSEXP);
    Rcpp::traits::input_parameter< const String >::type kernel_type(kernel_typeSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type delta_x_all(delta_x_allSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type output(outputSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd >::type A_all_v(A_all_vSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type sort_d_all_ix(sort_d_all_ixSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type sigma_2_vec(sigma_2_vecSEXP);
    Rcpp::traits::input_parameter< const VectorXi >::type num_neighbors_vec(num_neighbors_vecSEXP);
    Rcpp::traits::input_parameter< const double >::type tilde_nu(tilde_nuSEXP);
    Rcpp::traits::input_parameter< const int >::type D(DSEXP);
    Rcpp::traits::input_parameter< const VectorXd >::type n_t_record(n_t_recordSEXP);
    Rcpp::traits::input_parameter< float >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxIte(maxIteSEXP);
    rcpp_result_gen = Rcpp::wrap(IKF_CG_particle_cell(param, kernel_type, delta_x_all, output, A_all_v, sort_d_all_ix, sigma_2_vec, num_neighbors_vec, tilde_nu, D, n_t_record, tol, maxIte));
    return rcpp_result_gen;
END_RCPP
}
