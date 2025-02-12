#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _FastGaSP_Construct_G_exp(SEXP, SEXP);
extern SEXP _FastGaSP_Construct_G_matern_5_2(SEXP, SEXP);
extern SEXP _FastGaSP_Construct_W_exp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Construct_W_matern_5_2(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Construct_W0_exp(SEXP, SEXP);
extern SEXP _FastGaSP_Construct_W0_matern_5_2(SEXP, SEXP);
extern SEXP _FastGaSP_Get_C_R_K_Q(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_L_inv_y(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_log_det_S2(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_m_a_pred(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_Q_K(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_s_1st(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_S_KK(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Kalman_smoother(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Sample_KF(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Sample_KF_post(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_L_t_y(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_L_y(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_L_t_inv_y(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_Get_R_y(SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_fmou_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_matern_5_2_funct(SEXP,SEXP);
extern SEXP _FastGaSP_rcppeigen_get_chol(SEXP);
extern SEXP _FastGaSP_Get_G_log_det_cov(SEXP,SEXP,SEXP,SEXP, SEXP);
extern SEXP _FastGaSP_Optimization_Stiefel_Manifold(SEXP,SEXP,SEXP);
extern SEXP _FastGaSP_F_Funct(SEXP,SEXP);
extern SEXP _FastGaSP_IKF_CG_particle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_A_t_times_x_particle(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_A_times_x_particle(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_IKF_CG_particle_two_interact(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _FastGaSP_IKF_CG_particle_cell(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);





static const R_CallMethodDef CallEntries[] = {
    {"_FastGaSP_Construct_G_exp",               (DL_FUNC) &_FastGaSP_Construct_G_exp,                      2},
    {"_FastGaSP_Construct_G_matern_5_2",        (DL_FUNC) &_FastGaSP_Construct_G_matern_5_2,               2},
    {"_FastGaSP_Construct_W_exp",               (DL_FUNC) &_FastGaSP_Construct_W_exp,                      4},
    {"_FastGaSP_Construct_W_matern_5_2",        (DL_FUNC) &_FastGaSP_Construct_W_matern_5_2,               4},
    {"_FastGaSP_Construct_W0_exp",              (DL_FUNC) &_FastGaSP_Construct_W0_exp,                     2},
    {"_FastGaSP_Construct_W0_matern_5_2",       (DL_FUNC) &_FastGaSP_Construct_W0_matern_5_2,              2},
    {"_FastGaSP_Get_C_R_K_Q",                   (DL_FUNC) &_FastGaSP_Get_C_R_K_Q,                          5},
    {"_FastGaSP_Get_L_inv_y",                   (DL_FUNC) &_FastGaSP_Get_L_inv_y,                          4},
    {"_FastGaSP_Get_log_det_S2",                (DL_FUNC) &_FastGaSP_Get_log_det_S2,                       5},
    {"_FastGaSP_Get_m_a_pred",                  (DL_FUNC) &_FastGaSP_Get_m_a_pred,                         4},
    {"_FastGaSP_Get_Q_K",                       (DL_FUNC) &_FastGaSP_Get_Q_K,                              4},
    {"_FastGaSP_Get_s_1st",                     (DL_FUNC) &_FastGaSP_Get_s_1st,                            4},
    {"_FastGaSP_Get_S_KK",                      (DL_FUNC) &_FastGaSP_Get_S_KK,                             4},
    {"_FastGaSP_Kalman_smoother",               (DL_FUNC) &_FastGaSP_Kalman_smoother,                      7},
    {"_FastGaSP_Sample_KF",                     (DL_FUNC) &_FastGaSP_Sample_KF,                            6},
    {"_FastGaSP_Sample_KF_post",                (DL_FUNC) &_FastGaSP_Sample_KF_post,                       9},
    {"_FastGaSP_Get_L_t_y",                     (DL_FUNC) &_FastGaSP_Get_L_t_y,                            4},
    {"_FastGaSP_Get_L_y",                       (DL_FUNC) &_FastGaSP_Get_L_y,                              4},
    {"_FastGaSP_Get_L_t_inv_y",                 (DL_FUNC) &_FastGaSP_Get_L_t_inv_y,                        4},
    {"_FastGaSP_Get_R_y",                       (DL_FUNC) &_FastGaSP_Get_R_y,                              4},
    {"_FastGaSP_fmou_cpp", 	                (DL_FUNC) &_FastGaSP_fmou_cpp,                                 13},
    {"_FastGaSP_matern_5_2_funct", 	        (DL_FUNC) &_FastGaSP_matern_5_2_funct,                         2},
    {"_FastGaSP_rcppeigen_get_chol", 	        (DL_FUNC) &_FastGaSP_rcppeigen_get_chol,                   1},
    {"_FastGaSP_Get_G_log_det_cov", 	        (DL_FUNC) &_FastGaSP_Get_G_log_det_cov,                    5},
    {"_FastGaSP_Optimization_Stiefel_Manifold", (DL_FUNC) &_FastGaSP_Optimization_Stiefel_Manifold,        3},
    {"_FastGaSP_F_Funct",                       (DL_FUNC) &_FastGaSP_F_Funct,                              2},
    {"_FastGaSP_IKF_CG_particle",               (DL_FUNC) &_FastGaSP_IKF_CG_particle,                      12},
    {"_FastGaSP_A_t_times_x_particle",          (DL_FUNC) &_FastGaSP_A_t_times_x_particle,                 5},
    {"_FastGaSP_A_times_x_particle",            (DL_FUNC) &_FastGaSP_A_times_x_particle,                   5},
    {"_FastGaSP_IKF_CG_particle_two_interact",  (DL_FUNC) &_FastGaSP_IKF_CG_particle_two_interact,         18},
    {"_FastGaSP_IKF_CG_particle_cell",          (DL_FUNC) &_FastGaSP_IKF_CG_particle_cell,                 13},
    {NULL, NULL, 0}
};

void R_init_FastGaSP(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
