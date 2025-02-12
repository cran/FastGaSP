
// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*- 
 
// we only include RcppEigen.h which pulls Rcpp.h in for us 
#include <RcppEigen.h> 
#include <Rcpp.h> 
#include <iostream>
#include <vector>
#include <cmath> 
#include "ctools.h"
// [[Rcpp::depends(RcppEigen)]] 
//#include <chrono>
//#include <random>

using namespace Rcpp;
using namespace std;
using namespace Eigen; 
 
 

//July 16, 2018
////Construct_W0_matern_5_2 
// [[Rcpp::export]] 
MatrixXd Construct_W0_matern_5_2(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(3,3); 
  //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  W0(0,0)=sigma2; 
  W0(0,2)=W0(2,0)=-sigma2*pow(lambda,2.0)/3.0; 
  W0(1,1)=sigma2*pow(lambda,2.0)/3.0; 
  W0(2,2)=sigma2*pow(lambda,4.0); 

  return W0; 
} 

//July 16, 2018
////Construct_W0_exp 
// [[Rcpp::export]] 
MatrixXd Construct_W0_exp(const double sigma2, const double lambda){ 
  //int num_dim=sigma2.size(); 
  
  Eigen::MatrixXd W0= Eigen::MatrixXd::Zero(1,1); 

  W0(0,0)=sigma2; 
  
  return W0; 
} 


 
 ////Construct_G_matern_5_2 
 // [[Rcpp::export]] 
 List Construct_G_matern_5_2(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
   int num_obs=delta_x.size()+1; 
   //int num_dim=lambda.size();  
   List GG(num_obs);  
   GG[0]=Eigen::MatrixXd::Zero(3,3); 
   
   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
   
   // num_dim list, each is 3(num_obs)\times 3 list 
  // for(int i_GG=0;i_GG<num_dim;i_GG++){ 
  //   Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);  //the first row has all zeros  
     for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
       int j_GG_1=j_GG+1;    
       d(0,0)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)+2*lambda*delta_x[j_GG]+2; 
       d(1,0)=-pow(lambda,3.0)*pow(delta_x[j_GG],2.0); 
       d(2,0)=pow(lambda,4.0)*pow(delta_x[j_GG],2.0)-2*pow(lambda,3.0)*delta_x[j_GG]; 
       d(0,1)=2*(lambda*pow(delta_x[j_GG],2.0)+delta_x[j_GG]); 
       d(1,1)=-2*(pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-lambda*delta_x[j_GG]-1); 
       d(2,1)=2*(pow(lambda,3.0)*pow(delta_x[j_GG],2.0)-3*pow(lambda,2.0)*delta_x[j_GG]); 
       d(0,2)=pow(delta_x[j_GG],2); 
       d(1,2)=2*delta_x[j_GG]-lambda*pow(delta_x[j_GG],2.0); 
       d(2,2)=pow(lambda,2.0)*pow(delta_x[j_GG],2.0)-4*lambda*delta_x[j_GG]+2;     
       d=exp(-lambda*delta_x[j_GG])/2.0*d;
       GG[j_GG_1]=d; 
     } 
    // GG[i_GG]=d; 
   //} 
   return GG; 
} 

////Construct_G_exp
//' @keywords internal
//' @noRd
// [[Rcpp::export]] 
List Construct_G_exp(Eigen::VectorXd delta_x, double lambda){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=lambda.size();  
  List GG(num_obs);  
  GG[0]=Eigen::MatrixXd::Zero(1,1); 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_GG=0;j_GG<(num_obs-1);j_GG++){ 
    d(0,0)=exp(-delta_x[j_GG]*lambda);
    GG[j_GG+1]=d; 
  }

  return GG;
}

////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_matern_5_2(double sigma2,Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(3,3); //the first row has all zeros  
  
  
 // List Wi(num_obs); 
  //for(int i_Wi=0;i_Wi<num_dim;i_Wi++){ 
    //Eigen::MatrixXd d= Eigen::MatrixXd::Zero(num_obs,9);   
    double  lambda_delta_x;
    double exp_neg_2_lambda_delta_x;
    int  j_Wi_1;
    for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
      j_Wi_1= j_Wi+1; 
      lambda_delta_x=lambda*delta_x[j_Wi];  //close and jump then it is... 
      exp_neg_2_lambda_delta_x=exp(-2*lambda_delta_x); 
      
      d(0,0)=(exp_neg_2_lambda_delta_x*(3+6*lambda_delta_x+6*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-3 )/(-4*pow(lambda,5.0)); 
      d(1, 0)=  d(0, 1)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],4.0)/2.0; 
      d(2, 0)=  d(0, 2)=(exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)+4*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))-1 )/(4*pow(lambda,3.0)); 
      d(1, 1)= (exp_neg_2_lambda_delta_x*(1+2*lambda_delta_x+2*pow(lambda_delta_x,2.0)-4*pow(lambda_delta_x,3.0)+2*pow(lambda_delta_x,4.0))-1 )/(-4*pow(lambda,3.0)); 
      d(2, 1)=  d(1, 2)=exp_neg_2_lambda_delta_x*pow(delta_x[j_Wi],2.0)*(4-4*lambda_delta_x+pow(lambda_delta_x,2.0) )/2.0; 
      d(2, 2)=(exp_neg_2_lambda_delta_x*(-3+10*lambda_delta_x-22*pow(lambda_delta_x,2.0)+12*pow(lambda_delta_x,3.0)-2*pow(lambda_delta_x,4.0))+3 )/(4*lambda)     ;  
      d=d*(4*sigma2*pow(lambda,5.0)/3.0); 
      Wi[j_Wi_1]=d; 
      
    //} 
  } 
  return Wi; 
}


////Construct_W_matern_5_2  
// [[Rcpp::export]] 
List Construct_W_exp(double sigma2, Eigen::VectorXd delta_x, double lambda, MatrixXd W0){  //be careful about delta_x,lambda if one only wants to sample one 
  int num_obs=delta_x.size()+1; 
  //int num_dim=sigma2.size();  
  List Wi(num_obs);  
  Wi[0]=W0; 
  Eigen::MatrixXd d= Eigen::MatrixXd::Zero(1,1); 
  
  for(int j_Wi=0;j_Wi<(num_obs-1);j_Wi++){ 
    d(0,0)=1-exp(-2*delta_x[j_Wi]*lambda);
    Wi[j_Wi+1]=d;
  }
  
  return Wi;
}

////Get_Q_K  
// [[Rcpp::export]] 
List Get_Q_K(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV){ 

   int n=GG.size();
   int k=C0.rows();
   
   Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
   Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
   Eigen::MatrixXd C=C0;
    
   Eigen::MatrixXd GG_matrix;
   Eigen::MatrixXd W_matrix;
   
   Eigen::MatrixXd RR;
   
      
   // num_dim list, each is 3(num_obs)\times 3 list 
   for(int t=0;t<n;t++){ 
     GG_matrix=GG[t];
     W_matrix=W[t];
     RR=GG_matrix*C*GG_matrix.transpose()+W_matrix;
     //Q[t]=RR(0,0);
     Q[t]=RR(0,0)+VV;
     K.row(t)=RR.col(0).transpose()/Q[t];
     C=RR-RR.col(0)*RR.row(0)/Q[t];
   }

   List return_list;
   return_list.push_back(Q);
   return_list.push_back(K);
   
   return return_list;
}



////Get_log_det_S2
// [[Rcpp::export]] 
List Get_log_det_S2(const Eigen::VectorXd param,const bool have_noise,const Eigen::VectorXd delta_x,const Eigen::VectorXd output,
                    const  String kernel_type){
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  
  int n=output.rows();
  
  double gamma=1.0/exp(param[0]);
  
  double VV=0;
  if(have_noise){
    VV=exp(param[1]);
  }
  
  Eigen::MatrixXd    W0;
  List    GG;
  
  List    W;
  List    Q_K;
  
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
  
  //param.tail(k).array().exp().matrix();

  W0=Construct_W0_matern_5_2(1.0,lambda);   
  GG=Construct_G_matern_5_2(delta_x,lambda);  
  W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0);
  
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x,lambda);  
    W=Construct_W_exp(1.0,delta_x,lambda,W0);
    
  }

  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  double log_det_R=Q.array().log().sum();
  List return_vec;
  return_vec.push_back(log_det_R);
    
  //return_list.push_back(log_det_R);
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1[t]=(output[t]-a[0]);
    
   // Y_minus_a_1_scaled_vec[t]=(output[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  

  ///sigma2_S2=(FRFt+sigma_2*eta)  ####using this will lead to |R|
    
    
  double S2=(Y_minus_a_1.array()*Y_minus_a_1.array()/Q.array()).sum(); 
    
  return_vec.push_back(S2);
  
  return return_vec;
}




//  this function output the L^{-1}y given  Q matrix and K matrix
// [[Rcpp::export]] 
VectorXd Get_L_inv_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  
  //Eigen::VectorXd Q=Q_K[0];
  //Eigen::MatrixXd K=Q_K[1];
  
  //double log_det_R=Q.array().log().sum();
  //List return_vec;
  //return_vec.push_back(log_det_R);
  
  //return_list.push_back(log_det_R);
  int n=GG.size();
  //int k=C0.rows();
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(n);
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1[t]=(output[t]-a[0]);
    
    // Y_minus_a_1_scaled_vec[t]=(output[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output[t]-a[0]);
  }
  
  ///sigma2_S2=(FRFt+sigma_2*eta)  ####using this will lead to |R|
  
  Eigen::VectorXd res=(Y_minus_a_1.array()/Q.array().sqrt()).matrix();
  return res;
  //double S2=(Y_minus_a_1.array()*Y_minus_a_1.array()/Q.array()).sum(); 
  
}



/////////the following code is for prediction 
////Get_C_R_K_pred, C and R is for smoothing, K is for filtering the mean  
// let me also output Q here for determinant
//it replaces the previous C R K code
// [[Rcpp::export]] 
List Get_C_R_K_Q(const VectorXi index, const List GG,const List  W,const Eigen::MatrixXd C0, double VV){ 
  
  
  //index is a sequence where 0 means NA and 1 means observations
  int n=GG.size();//total number
  int k=C0.rows();
  
  //Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd K=Eigen::MatrixXd::Zero(n,k);
  List C(n+1);
  C[0]=C0;
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd W_matrix;
  Eigen::VectorXd Q=Eigen::VectorXd::Zero(n);
  
  Eigen::MatrixXd RR;
  
  List R(n);
  
  Eigen::MatrixXd C_cur=C[0];
  //int index_num=0;
  
  for(int t=0;t<n;t++){
    
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      W_matrix=W[t];
      RR=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      R[t]=RR;
      
      //Q[t]=RR[t](0,0);
      Q[t]=RR(0,0)+VV;
      
      
      K.row(t)=RR.col(0).transpose()/Q[t];
      C[t+1]=RR-RR.col(0)*RR.row(0)/Q[t];
      C_cur=C[t+1];
      
    }else{
      GG_matrix=GG[t];
      W_matrix=W[t];
      
      R[t]=GG_matrix*C_cur*GG_matrix.transpose()+W_matrix;
      C[t+1]=C_cur=R[t];   
      
    }
  }
  
  List return_list;
  return_list.push_back(C);
  return_list.push_back(R);
  return_list.push_back(K);
  return_list.push_back(Q);
  
  return return_list;
}



// [[Rcpp::export]] 
List Get_m_a_pred(const VectorXi index, const Eigen::VectorXd output_vec,const List GG,const Eigen::MatrixXd K){
  
  //output_KF
  int n=GG.size();//total number
  
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  int k=K.cols();
  
  List m(n);  //m should have n+1 item but the first one is (0,0,0) so it's omitted.
  List a(n);
  Eigen::VectorXd a_vec;
  Eigen::VectorXd m_cur=Eigen::VectorXd::Zero(k); //3\times 1
  
  // =Eigen::MatrixXd::Zero(k,n2);
  
  //Eigen::MatrixXd Y_minus_a_1_scaled_matrix=Eigen::MatrixXd::Zero(n1,n2); 
  //Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::MatrixXd GG_matrix;
  
  //int index_num=0;
  int index_count=0;
  for(int t=0;t<n;t++){
    if(index[t]==1){
      //index_num=index_num+1;
      GG_matrix=GG[t];
      a_vec=GG_matrix*m_cur;
      a[t]=a_vec;
      m[t]=a_vec+K.row(t).transpose()*(output_vec[index_count]-a_vec[0]);
      m_cur=m[t];
      index_count=index_count+1;
    }else{
      GG_matrix=GG[t];
      m[t]=a[t]=GG_matrix*m_cur;
      m_cur=m[t];
    }
  }
  
  List return_list;
  
  return_list.push_back(m);
  return_list.push_back(a);
  
  return return_list;
}


////Kalman_smoother, the index is 0 and 1 where 0 means missing value
// [[Rcpp::export]] 
List Get_S_KK(const VectorXi index,const List GG,const List C,const List R){
  
  int n=GG.size();
  //int k=3; //matern_5_2
  
  List S(n);
  List KK(n-1);
  
  S[n-1]= C[n];
  
  Eigen::MatrixXd GG_matrix;
  Eigen::MatrixXd C_matrix;
  Eigen::MatrixXd R_matrix;
  Eigen::MatrixXd KK_matrix;
  
  Eigen::MatrixXd S_cur=C[n];
  for(int t=n-2;t>=0;t--){
    GG_matrix=GG[t+1];
    R_matrix=R[t+1];
    C_matrix=C[t+1];  // C should have n+1 items
    
    LLT<MatrixXd> lltOfR(R_matrix);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
    
    KK_matrix= (L.transpose().triangularView<Upper>().solve(L.triangularView<Lower>().solve(GG_matrix*C_matrix))); 
    KK[t]=KK_matrix;
    
    S[t]=C_matrix-KK_matrix.transpose()*(R_matrix-S_cur)*KK_matrix;
    S_cur=S[t];
  }
  
  List return_list;
  
  return_list.push_back(S);
  return_list.push_back(KK);
  
  return return_list;
  
}

// [[Rcpp::export]] 
MatrixXd Get_s_1st(const List m,const List a,const List C,const List KK){
  
  int n=C.size()-1;
  //int k=3; //matern_5_2
  
  Eigen::VectorXd s=m[n-1];
  //int n_2=s.cols();
  VectorXd s_1st=Eigen::VectorXd::Zero(n);
  MatrixXd KK_matrix;
  s_1st[n-1]=s[0];
  VectorXd a_vec;
  VectorXd m_vec;
  for(int t=n-2;t>=0;t--){
    KK_matrix=KK[t];
    a_vec=a[t+1];
    m_vec=m[t];
    s=m_vec+KK_matrix.transpose()*(s-a_vec);
    s_1st[t]=s[0];
  }
  
  return s_1st; 
}




// [[Rcpp::export]] 
List Kalman_smoother(const VectorXd param,const bool have_noise,const VectorXi index_obs, 
                     const VectorXd delta_x_all, const VectorXd output, const double sigma_2_hat,
                     const String kernel_type){
  //int n1=output_KF.rows();
  //int n2=output_KF.cols();
  
  double gamma=1.0/exp(param[0]);
  //double lambda=sqrt(5.0)/gamma;
  
  
  double VV=0;
  if(have_noise){
    VV=exp(param[1]); 
  }
  int n=delta_x_all.size()+1;
  //param.tail(k).array().exp().matrix();
  
  
  Eigen::MatrixXd W0;
  List GG;
  List W;
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    
    //param.tail(k).array().exp().matrix();
    
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x_all,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x_all,lambda,W0);
    
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x_all,lambda);  
    W=Construct_W_exp(1.0,delta_x_all,lambda,W0);
    
  }
  
  
  
  List C_R_K_Q=Get_C_R_K_Q(index_obs, GG,W, W0, VV);
  
  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K_Q[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K_Q[0],C_R_K_Q[1]);
  
  
  MatrixXd s_1st=Get_s_1st(m_a[0],m_a[1], C_R_K_Q[0],S_KK[1]);
  
  List return_list;
  
  List S= S_KK[0];
  MatrixXd S_matrix;
  VectorXd Var=VectorXd::Zero(n);
  for(int t=0;t<n;t++){
    S_matrix=S[t];
    Var[t]=S_matrix(0,0);
  }
  Var=((Var.array()+VV)*sigma_2_hat).matrix();
  
  return_list.push_back(s_1st);
  return_list.push_back(Var);
  /*
  List return_list;
  
  return_list.push_back(W);
  */
  return return_list;
  
}

//the following code is for generating samples
//Sample the prior from Kalman Filter 
//sample_type 0 means all theta, 1 means first theta, 2 means data (contains the noise)
// [[Rcpp::export]] 
MatrixXd  Sample_KF(const List GG,const List  W,const Eigen::MatrixXd C0,const double VV,const  String kernel_type, const int sample_type){
  int n=GG.size();
  

  LLT<MatrixXd> lltOfM(C0);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfM.matrixL();   //retrieve factor L  in the decomposition
  
  Eigen::MatrixXd theta_sample;
  
  Eigen::MatrixXd W_mat;
  Eigen::MatrixXd G_mat;
  
  Eigen::MatrixXd random_norm_theta;

  if(kernel_type=="matern_5_2"){
     theta_sample=Eigen::MatrixXd::Zero(3,n);
     random_norm_theta=Eigen::MatrixXd::Zero(3,n+1);
    
     for(int t=0;t<3;t++){
        for(int i=0;i<(n+1);i++){
           random_norm_theta(t,i)= R::rnorm(0,1.0);
        }
      }

  }else if(kernel_type=="exp"){
      theta_sample=Eigen::MatrixXd::Zero(1,n);
    
      random_norm_theta=Eigen::MatrixXd::Zero(1,n+1);
      for(int i=0;i<(n+1);i++){
          random_norm_theta(0,i)= R::rnorm(0,1.0);
      }
  }
  
  Eigen::VectorXd sample_cur=L*random_norm_theta.col(0); // current sample
  
  for(int t=0;t<n;t++){
    W_mat=W[t];
    G_mat=GG[t];
    LLT<MatrixXd> lltOfM(W_mat);    
    L = lltOfM.matrixL();
    
    sample_cur=G_mat*sample_cur+L*random_norm_theta.col(t+1);
    //sample_record[t]=sample_cur(0)+sqrt(VV)*random_norm_y(t);
    theta_sample.col(t)=sample_cur;
  }
  
  Eigen::MatrixXd sample_record;
  if(sample_type==0){
    sample_record=theta_sample;
  }else if(sample_type==1){
    sample_record=theta_sample.row(0);
  }else{
    Eigen::VectorXd random_norm_y=Eigen::VectorXd::Zero(n);
    for(int i=0;i<n;i++){
      random_norm_y[i]=R::rnorm(0,1.0);
    }
    sample_record=theta_sample.row(0);
    sample_record=sample_record+sqrt(VV)*random_norm_y.transpose();
  }
  return sample_record;
}


//Sample the posterior from Kalman Filter 
//sample_type 0 means all theta, 1 means first theta, 2 means data (contains the noise)
// [[Rcpp::export]] 
MatrixXd  Sample_KF_post(const VectorXi index_obs,const List C_R_K_Q, const Eigen::MatrixXd W0,const List GG,const List W,const double VV,
                         const VectorXd output,String kernel_type, const int sample_type){
  
  int n=GG.size();
  
  List m_a=Get_m_a_pred(index_obs, output, GG,C_R_K_Q[2]);
  
  List S_KK=Get_S_KK(index_obs, GG,C_R_K_Q[0],C_R_K_Q[1]);
  
  MatrixXd S_matrix;
  List S_list=S_KK[0];
  List C=C_R_K_Q[0];
  Eigen::MatrixXd C_here=C[n]; //it looks C has n+1 terms but m only have n terms
  
  Eigen::MatrixXd theta_sample;
  Eigen::MatrixXd random_norm_theta;
  
  if(kernel_type=="matern_5_2"){
    theta_sample=Eigen::MatrixXd::Zero(3,n);
    random_norm_theta=Eigen::MatrixXd::Zero(3,n+1);
    
    for(int t=0;t<3;t++){
      for(int i=0;i<(n+1);i++){
        random_norm_theta(t,i)= R::rnorm(0,1.0);
      }
    }
  }else if(kernel_type=="exp"){
    theta_sample=Eigen::MatrixXd::Zero(1,n);
    
    random_norm_theta=Eigen::MatrixXd::Zero(1,n+1);
    for(int i=0;i<(n+1);i++){
      random_norm_theta(0,i)= R::rnorm(0,1.0);
    }
  }
  
  Eigen::VectorXd theta_minus_a; 
  List m=m_a[0];
  
  Eigen::VectorXd m_here=m[n-1];
  
  LLT<MatrixXd> lltOfC_here(C_here);    // compute the cholesky decomposition of R called lltofR
  MatrixXd L_C_here = lltOfC_here.matrixL();   //retrieve factor L  in the decomposition
  
  Eigen::VectorXd theta_sample_here;
  
  theta_sample_here=m_here+L_C_here*random_norm_theta.col(n-1);
  
  theta_sample.col(n-1)=theta_sample_here;
  
  //done for sample 1
  List KK_list=S_KK[1];
  
  List a=m_a[1];
  VectorXd a_here;
  //VectorXd m_here;
  
  VectorXd h;
  MatrixXd KK_matrix;
  MatrixXd S_matrix_plus_1;
  
  MatrixXd H_matrix;
  
  for(int t=n-2;t>=0;t--){

    KK_matrix=KK_list[t];
    a_here=a[t+1];
    m_here=m[t];
    h=m_here+KK_matrix.transpose()*(theta_sample_here-a_here);
    
    //not sure why is this
    S_matrix=S_list[t];
    S_matrix_plus_1=S_list[t+1];
    
    
    H_matrix=S_matrix-KK_matrix.transpose()*S_matrix_plus_1*KK_matrix;
    
    //LLT<MatrixXd> lltOfS(S_matrix);    // compute the cholesky decomposition of R called lltofR
    // MatrixXd L_S = lltOfS.matrixL();   //retrieve factor L  in the decomposition
    
    LLT<MatrixXd> lltOfH(H_matrix);    // compute the cholesky decomposition of R called lltofR
    MatrixXd L_H = lltOfH.matrixL();   //retrieve factor L  in the decomposition
    
    theta_sample_here=h+L_H*random_norm_theta.col(t);
    theta_sample.col(t)=theta_sample_here;
    
  }
  
  //generat the sample record based different types of tasks
  Eigen::MatrixXd sample_record;
  if(sample_type==0){
    sample_record=theta_sample;
  }else if(sample_type==1){
    sample_record=theta_sample.row(0);
  }else{
    Eigen::VectorXd random_norm_y=Eigen::VectorXd::Zero(n);
    for(int i=0;i<n;i++){
      random_norm_y[i]=R::rnorm(0,1.0);
    }
    sample_record=theta_sample.row(0);
    sample_record=sample_record+sqrt(VV)*random_norm_y.transpose();
  }
  return sample_record;
  
}





///IKF functions
// functions only using matrices
// [[Rcpp::export]] 
VectorXd Get_L_t_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  int n=GG.size();
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  
  Eigen::VectorXd res=Eigen::VectorXd::Zero(n);
  
  res(n-1)=sqrt_Q(n-1)*output(n-1);
  
  Eigen::MatrixXd GG_matrix;
  if(n>=2){
    GG_matrix=GG[n-1];
    Eigen::VectorXd u=GG_matrix*K.row(n-2).transpose();
    res(n-2)=sqrt_Q(n-2)*(u(0)*output(n-1)+output(n-2)); // this is for Matern 2.5
    if(n>=3){ 
      //Nov 2022
      Eigen::MatrixXd GG_matrix_plus_1;
      Eigen::MatrixXd g;
      Eigen::VectorXd GG_K;
      Eigen::VectorXd g_GG_K;
      g=GG_matrix.row(0)*output(n-1);
      for(int t=n-3;t>=0;t--){
        GG_matrix_plus_1=GG[t+1];
        GG_K=GG_matrix_plus_1*K.row(t).transpose();
        g_GG_K=g*GG_K;
        res(t)=sqrt_Q(t)*(g_GG_K(0)+GG_K(0)*output(t+1)+output(t));
        g=g*GG_matrix_plus_1+GG_matrix_plus_1.row(0)*output(t+1);
      }
    }
  }
  
  return res;
}

// [[Rcpp::export]] 
VectorXd Get_L_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  int n=GG.size();
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::MatrixXd m=Eigen::VectorXd::Zero(K.cols());
  
  Eigen::VectorXd a;
  
  //Eigen::VectorXd Y_minus_a_1_scaled_vec=Eigen::VectorXd::Zero(n); 
  //Eigen::VectorXd Y_minus_a_1=Eigen::VectorXd::Zero(n); 
  
  Eigen::VectorXd res=Eigen::VectorXd::Zero(n);
  Eigen::MatrixXd GG_matrix;
  
  for(int t=0;t<n;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    res[t]=a[0]+output[t]*sqrt_Q[t];
    
    //Y_minus_a_1[t]=(res[t]-a[0]);
    
    // Y_minus_a_1_scaled_vec[t]=(res[t]-a[0])/sqrt_Q[t];
    m=a+K.row(t).transpose()*(res[t]-a[0]);
  }
  //return output;
  
  return res;
}

// [[Rcpp::export]] 
VectorXd Get_L_t_inv_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  int n=GG.size();
  
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::VectorXd res=Eigen::VectorXd::Zero(n);
  
  res(n-1)=output(n-1)/sqrt_Q(n-1);
  
  Eigen::MatrixXd GG_matrix;
  if(n>=2){
    GG_matrix=GG[n-1];
    Eigen::VectorXd u=GG_matrix*K.row(n-2).transpose();
    res(n-2)=output(n-2)/sqrt_Q(n-2)-u(0)*res(n-1);//sqrt_Q(n-2)*(u(0)*output(n-1)+output(n-2)); // this is for Matern 2.5
    if(n>=3){ 
      //Aug 2023
      Eigen::MatrixXd GG_matrix_plus_1;
      Eigen::MatrixXd g;
      Eigen::VectorXd GG_K;
      Eigen::VectorXd g_GG_K;
      g=GG_matrix.row(0)*res(n-1);
      for(int t=n-3;t>=0;t--){
        GG_matrix_plus_1=GG[t+1];
        GG_K=GG_matrix_plus_1*K.row(t).transpose();
        g_GG_K=g*GG_K;
        res(t)=output(t)/sqrt_Q(t)-g_GG_K(0)-GG_K(0)*res(t+1);//sqrt_Q(t)*(g_GG_K(0)+GG_K(0)*output(t+1)+output(t));
        g=g*GG_matrix_plus_1+GG_matrix_plus_1.row(0)*res(t+1);
      }
    }
  }
  
  return res;
}


// [[Rcpp::export]] 
VectorXd Get_R_y(const List GG,const Eigen::VectorXd Q, const Eigen::MatrixXd K,  const Eigen::VectorXd output){
  
  Eigen::VectorXd tilde_z = Get_L_t_y(GG,Q,K,output);
  Eigen::VectorXd res = Get_L_y(GG,Q,K,tilde_z);
  
  return res;
}


// [[Rcpp::export]] 
Eigen::MatrixXd Get_Y_minus_a_1_scaled_matrix_2d(const Eigen::MatrixXd output_KF,const List GG,const Eigen::VectorXd Q,const Eigen::MatrixXd K){
  
  int n1=output_KF.rows();
  int n2=output_KF.cols();
  int k=K.cols();
  
  Eigen::MatrixXd m=Eigen::MatrixXd::Zero(k,n2);
  
  Eigen::MatrixXd Y_minus_a_1_scaled_matrix=Eigen::MatrixXd::Zero(n1,n2); 
  Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  
  Eigen::MatrixXd GG_matrix;
  
  Eigen::MatrixXd a;
  for(int t=0;t<n1;t++){
    GG_matrix=GG[t];
    a=GG_matrix*m;
    Y_minus_a_1_scaled_matrix.row(t)=(output_KF.row(t)-a.row(0))/sqrt_Q[t];
    m=a+K.row(t).transpose()*(output_KF.row(t)-a.row(0));
  }
  
  
  return Y_minus_a_1_scaled_matrix;
}

// [[Rcpp::export]] 
double F_Funct(const Eigen::MatrixXd A_cur,const List G){
  double return_val=0;
  int d=A_cur.cols();
  Eigen::MatrixXd G_matrix;
  for(int i=0;i<d;i++){
    G_matrix=G[i];
    return_val=return_val+A_cur.col(i).transpose()*G_matrix*A_cur.col(i);
  }
  return -return_val;
}

// [[Rcpp::export]] 
List Get_G_log_det_cov(const Eigen::VectorXd param,const Eigen::MatrixXd output,const Eigen::VectorXd delta_x,int d,
                       const String kernel_type){
  Eigen::VectorXd gamma=(1.0/param.head(d).array().exp()).matrix();
  Eigen::VectorXd tau=param.tail(d).array().exp().matrix();
  
  
  Eigen::MatrixXd    output_t=output.transpose();
  
  Eigen::VectorXd Q;
  
  
  Eigen::MatrixXd W0;
  List GG;
  List W;
  List Q_K;
  
  Eigen::VectorXd VV=1.0/tau.array();
  Eigen::VectorXd log_det_cov=Eigen::VectorXd::Zero(d);
  
  //int n=output.cols();
  
  Eigen::MatrixXd    z;
  
  List return_list;
  
  Eigen::VectorXd lambda;
  if(kernel_type=="matern_5_2"){
    // lambda=sqrt(5.0)/gamma;
    lambda=(sqrt(5.0)/gamma.array()).matrix();
    
    for(int i=0;i<d;i++){  
      W0=Construct_W0_matern_5_2(1.0,lambda[i]);  
      
      GG=Construct_G_matern_5_2(delta_x,lambda[i]);  
      
      W=Construct_W_matern_5_2(1.0,delta_x,lambda[i],W0);
      
      Q_K=Get_Q_K(GG,W,W0,VV[i]);
      
      Q=Q_K[0];
      
      log_det_cov[i]=(tau[i]*Q.array()).log().sum();
      
      z=Get_Y_minus_a_1_scaled_matrix_2d(output_t,GG,Q_K[0],Q_K[1]);
      
      return_list.push_back( (1.0/tau[i]*(z.transpose()*z).array()).matrix());
      // return_list.push_back( ((z.transpose()*z).array()).matrix());
      
    }
  }else if(kernel_type=="exp"){
    lambda=(1.0/gamma.array()).matrix();
    
    for(int i=0;i<d;i++){  
      W0=Construct_W0_exp(1.0,lambda[i]);  
      GG=Construct_G_exp(delta_x,lambda[i]);  
      W=Construct_W_exp(1.0,delta_x,lambda[i],W0);
      
      
      
      Q_K=Get_Q_K(GG,W,W0,VV[i]);
      
      Q=Q_K[0];
      
      log_det_cov[i]=(tau[i]*Q.array()).log().sum();
      
      z=Get_Y_minus_a_1_scaled_matrix_2d(output_t,GG,Q_K[0],Q_K[1]);
      
      return_list.push_back( (1.0/tau[i]*(z.transpose()*z).array()).matrix());
      // return_list.push_back( ((z.transpose()*z).array()).matrix());
      
    }
    
  }
  return_list.push_back(log_det_cov);
  
  return return_list;
}


// [[Rcpp::export]] 
List KF_cpp_eigen(double F_t, double V_t, double G_t, double W_t, const Eigen::VectorXd &y) {
  int n = y.size();
  // direct results from KF and RTS
  Eigen::VectorXd a_record(n), R_record(n), f_record(n), Q_record(n), m_record(n), C_record(n), S_record(n), s_record(n);
  // primary off-diagonal of covariance matrix: Cov(z(t),z(t+1)|Y)
  Eigen::VectorXd b_record = Eigen::VectorXd::Zero(n - 1); 
  
  a_record(0) = 0;
  R_record(0) = W_t / (1 - G_t * G_t);
  f_record(0) = F_t * a_record(0);
  Q_record(0) = F_t * R_record(0) * F_t + V_t;
  m_record(0) = a_record(0) + R_record(0) * F_t * (1 / Q_record(0)) * (y(0) - f_record(0));
  C_record(0) = R_record(0) - R_record(0) * F_t * (1 / Q_record(0)) * F_t * R_record(0);
  
  for (int i = 1; i < n; ++i) {
    a_record(i) = G_t * m_record(i - 1);
    R_record(i) = G_t * C_record(i - 1) * G_t + W_t;
    f_record(i) = F_t * a_record(i);
    Q_record(i) = F_t * R_record(i) * F_t + V_t;
    m_record(i) = a_record(i) + R_record(i) * F_t * (1 / Q_record(i)) * (y(i) - f_record(i));
    C_record(i) = R_record(i) - R_record(i) * F_t * (1 / Q_record(i)) * F_t * R_record(i);
  }
  
  s_record(n - 1) = m_record(n - 1);
  S_record(n - 1) = C_record(n - 1);
  
  for (int i = n - 2; i >= 0; --i) {
    double s_cur = m_record(i) + C_record(i) * G_t * (1 / R_record(i + 1)) * (s_record(i + 1) - a_record(i + 1));
    double S_cur = C_record(i) - C_record(i) * G_t * (1 / R_record(i + 1)) * (R_record(i + 1) - S_record(i + 1)) * (1 / R_record(i + 1)) * G_t * C_record(i);
    double b_cur = C_record(i) * G_t * (1 / R_record(i + 1)) * S_record(i + 1);
    
    s_record(i) = s_cur;  // posterior mean
    S_record(i) = S_cur;  // posterior variance
    b_record(i) = b_cur; // posterior covariance Cov(z(t),z(t+1)|Y)
  }
  
  return List::create(Named("b_vec") = b_record,
                      Named("s_t") = s_record,
                      Named("f_t") = f_record,
                      Named("Q_t") = Q_record,
                      Named("a_t") = a_record,
                      Named("S_t") = S_record,
                      Named("R_t") = R_record,
                      Named("m_t") = m_record,
                      Named("C_t") = C_record);
}

// [[Rcpp::export]] 
double cubic_solver(const std::vector<double>& p) {
  if (std::abs(p[0]) < std::pow(std::numeric_limits<double>::epsilon(), 0.95)) {
    throw std::invalid_argument("Coefficient of highest power must not be zero!");
  }
  if (p.size() != 4) {
    throw std::invalid_argument("p is not a numeric or has not 4 elements!");
  }
  
  // Normalize coefficients
  std::vector<double> a(3);
  for (size_t i = 1; i < 4; ++i) {
    a[i - 1] = p[i] / p[0];
  }
  
  double Q = (a[0] * a[0] - 3 * a[1]) / 9;
  double R = (2 * a[0] * a[0] * a[0] - 9 * a[0] * a[1] + 27 * a[2]) / 54;
  
  std::vector<std::complex<double>> x(3);
  
  // Case 1 - 3 real roots
  if (R * R < Q * Q * Q) {
    double theta = std::acos(R / std::sqrt(Q * Q * Q));
    x[0] = -2 * std::sqrt(Q) * std::cos(theta / 3.0) - a[0] / 3.0;
    x[1] = -2 * std::sqrt(Q) * std::cos((theta + 2 * M_PI) / 3.0) - a[0] / 3.0;
    x[2] = -2 * std::sqrt(Q) * std::cos((theta - 2 * M_PI) / 3.0) - a[0] / 3.0;
  } else {
    std::complex<double> A = -std::copysign(1.0, R) * std::pow(std::abs(R) + std::sqrt(R * R - Q * Q * Q), 1.0 / 3.0);
    std::complex<double> B = (A == std::complex<double>(0, 0)) ? 0.0 : Q / A;
    
    x[0] = (A + B) - a[0] / 3.0;
    x[1] = -0.5 * (A + B) - a[0] / 3.0 + std::sqrt(3.0) * std::complex<double>(0, 1) * (A - B) / 2.0;
    x[2] = -0.5 * (A + B) - a[0] / 3.0 - std::sqrt(3.0) * std::complex<double>(0, 1) * (A - B) / 2.0;
  }
  
  // Filter real roots within [-1, 1]
  std::vector<double> filteredRoots;
  for (int i = 0; i < x.size(); ++i) {
    if (std::abs(x[i].imag()) < 1e-10) { // Check if the root is real
      double realRoot = x[i].real();
      if (realRoot > -1 && realRoot < 1) { // Check if it's in range (-1, 1)
        return realRoot; // Return the valid root
      }
    }
  }
  return -1;
}

// [[Rcpp::export]] 
List fmou_cpp(const Eigen::MatrixXd & output, int d,
          int M = 50, double threshold = 1e-4,
          bool est_U0 = true,  bool est_sigma0_2=true,
          bool track_iterations=false,
          bool track_neg_log_lik=false,
          Nullable<Eigen::MatrixXd> U0 = R_NilValue,
          Nullable<Eigen::MatrixXd> U_init = R_NilValue,
          Nullable<Eigen::VectorXd> rho_init = R_NilValue,
          Nullable<Eigen::VectorXd> sigma2_init = R_NilValue,
          Nullable<double> sigma0_2 = R_NilValue) {
  
  int n = output.cols();
  int k = output.rows();
  double tr_Y_Y_t = (output.array() * output.array()).sum();
  Eigen::MatrixXd output_2 = output * output.transpose();
  std::string kernel_type = "exp";
  Eigen::VectorXd delta_x = Eigen::VectorXd::Ones(output.cols() - 1);
  double trace_output_2 = output.array().square().sum();
  
  
  // Initialization
  // (1). U_cur
  Eigen::MatrixXd U_cur;
  if (est_U0) {
    if (U_init.isNotNull()) {
      U_cur = as<Eigen::MatrixXd>(U_init);
    } else {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(output, Eigen::ComputeThinU | Eigen::ComputeThinV);
      U_cur = svd.matrixU().leftCols(d);
    }
  } else {
    U_cur = as<Eigen::MatrixXd>(U0);
  }
  
  // (2). sigma2_0
  double sigma0_2_cur;
  if(est_sigma0_2){
    sigma0_2_cur=std::log(1.5);
  }
  else{
    sigma0_2_cur = as<double>(sigma0_2);
  }
  
  // (3). rho_cur
  Eigen::VectorXd rho_cur;
  if(rho_init.isNotNull()){
    rho_cur = as<Eigen::VectorXd>(rho_init);
  }
  else{
    rho_cur = Eigen::VectorXd::LinSpaced(d, 0.8, 0.99);
  }
  
  // (4). sigma2_cur
  Eigen::VectorXd sigma2_cur;
  if(sigma2_init.isNotNull()){
    sigma2_cur = as<Eigen::VectorXd>(sigma2_init);
  }
  else{
    sigma2_cur = Eigen::VectorXd::LinSpaced(d, 0.5, 1.0);
  }
  
  // (5). Z_hat_cur, diag_Sigma, off_diag_Sigma 
  Eigen::MatrixXd output_tilde = U_cur.transpose() * output;
  Eigen::MatrixXd Z_hat_cur(d, n);
  Eigen::MatrixXd diag_Sigma(d,n);
  Eigen::MatrixXd off_diag_Sigma(d, n-1);
  for (int l = 0; l < d; ++l) {
    List kf = KF_cpp_eigen(1.0, sigma0_2_cur, rho_cur(l), sigma2_cur(l), output_tilde.row(l));
    Z_hat_cur.row(l) = as<Eigen::VectorXd>(kf["s_t"]);
    diag_Sigma.row(l) = as<Eigen::VectorXd>(kf["S_t"]);
    off_diag_Sigma.row(l) = as<Eigen::VectorXd>(kf["b_vec"]);
  }
  
  Eigen::MatrixXd pred_cur = U_cur * Z_hat_cur;
  Eigen::MatrixXd pred_pre = pred_cur + Eigen::MatrixXd::Ones(pred_cur.rows(), pred_cur.cols());
  
  // 2. EM algorithm
  Eigen::MatrixXd record_rho = Eigen::MatrixXd::Zero(d, M);
  Eigen::MatrixXd record_sigma2 = Eigen::MatrixXd::Zero(d, M);
  Eigen::VectorXd record_sigma0_2 = Eigen::VectorXd::Zero(M);
  Eigen::VectorXd record_neg_log_lik = Eigen::VectorXd::Zero(M);
  Eigen::VectorXd log_Q = Eigen::VectorXd::Zero(d);
  
  int m = 1;
  while (sqrt((pred_cur - pred_pre).squaredNorm() / (k * n)) > threshold && m <= M) {
    pred_pre = pred_cur;
    
    // Update U
    if (est_U0) {
      Eigen::MatrixXd Z_hat_Y_t = Z_hat_cur * output.transpose();
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(Z_hat_Y_t, Eigen::ComputeThinU | Eigen::ComputeThinV);
      U_cur = svd.matrixV() * svd.matrixU().transpose();
    } else {
      U_cur = as<Eigen::MatrixXd>(U0);
    }
    output_tilde = U_cur.transpose() * output;
    
    // Update sigma0_2
    if (est_sigma0_2) {
      double tr_Z_hat_Y_t_V = (Z_hat_cur.array() * output_tilde.array()).sum();
      double tr_Z_t_Z = Z_hat_cur.squaredNorm();
      sigma0_2_cur = (tr_Y_Y_t - 2 * tr_Z_hat_Y_t_V + diag_Sigma.sum() + tr_Z_t_Z) / (n * k);
    } else {
      sigma0_2_cur = as<double>(sigma0_2);
    }
    record_sigma0_2(m - 1) = sigma0_2_cur;
    
    // Update rho and sigma2
    for (int l = 0; l < d; ++l) {
      // Update rho
      double p3 = (n-1) * (Z_hat_cur.row(l).segment(1,n-2).squaredNorm() +
                   diag_Sigma.row(l).segment(1, n-2).sum());
      double p2 = (2-n) * Z_hat_cur.row(l).tail(n - 1).dot(Z_hat_cur.row(l).head(n - 1)) +
        (2-n) * off_diag_Sigma.row(l).sum();
      double p1 = -Z_hat_cur.row(l).squaredNorm() - diag_Sigma.row(l).sum() -
        n * (Z_hat_cur.row(l).segment(1,n - 2).squaredNorm() +
        diag_Sigma.row(l).segment(1, n - 2).sum());
      double p0 = n * (Z_hat_cur.row(l).tail(n - 1).dot(Z_hat_cur.row(l).head(n - 1)) +
                       off_diag_Sigma.row(l).sum());
      NumericVector coeffs = NumericVector::create(p3, p2, p1, p0);
      std::vector<double> p(coeffs.begin(), coeffs.end());
      double rho_l_cur = cubic_solver(p);
      rho_cur(l) = rho_l_cur;
      record_rho(l, m - 1) = rho_l_cur;
      
      // Update sigma2
      double sigma2_l_cur = ((1 - rho_l_cur * rho_l_cur) * Z_hat_cur(l, 0) * Z_hat_cur(l, 0) +
                             (Z_hat_cur.row(l).tail(n - 1) - rho_l_cur * Z_hat_cur.row(l).head(n - 1)).squaredNorm() +
                             diag_Sigma.row(l).sum() +
                             rho_l_cur * rho_l_cur * diag_Sigma.row(l).segment(1, n - 2).sum() -
                             2 * rho_l_cur * off_diag_Sigma.row(l).sum()) /n;
      sigma2_cur(l) = sigma2_l_cur;
      record_sigma2(l, m - 1) = sigma2_l_cur;
      
      // Update Z_hat, diag_Sigma, off_diag_Sigma
      List kf = KF_cpp_eigen(1.0, sigma0_2_cur, rho_l_cur, sigma2_l_cur, output_tilde.row(l));
      Z_hat_cur.row(l) = as<Eigen::VectorXd>(kf["s_t"]);
      diag_Sigma.row(l) = as<Eigen::VectorXd>(kf["S_t"]);
      off_diag_Sigma.row(l) = as<Eigen::VectorXd>(kf["b_vec"]);
      log_Q(l) = as<Eigen::VectorXd>(kf["Q_t"]).array().log().sum();
    }
    // (Optional) compute and record negative log likelihood
    if(track_neg_log_lik){
      double S_2 = tr_Y_Y_t - (output_tilde.array()*Z_hat_cur.array()).sum();
      double neg_log_det=-0.5*log_Q.sum()+(n*d/2.0)*std::log(sigma0_2_cur);
      record_neg_log_lik(m-1)= -(neg_log_det-(n*k)/2*std::log(S_2));
    }
    
    // Update predictive mean of observations
    pred_cur = pred_cur = U_cur * Z_hat_cur;
    m++;
  }
  // compute 95% lower and upper bounds of the predictive mean of observations
  Eigen::MatrixXd pred_mean_var = Eigen::MatrixXd::Zero(k, n);
  for (int tt = 0; tt < n; ++tt) {
    pred_mean_var.col(tt) = (U_cur.array() * (diag_Sigma.col(tt).asDiagonal() * U_cur.transpose()).transpose().array()).rowwise().sum();
  }
  
  // Compute the prediction intervals
  Eigen::MatrixXd pred_mean_95lb = pred_cur - 1.96 * pred_mean_var.array().sqrt().matrix();
  Eigen::MatrixXd pred_mean_95ub = pred_cur + 1.96 * pred_mean_var.array().sqrt().matrix();
  
  List res;
  if(track_iterations){
    res = List::create(Named("output") = output,
                       Named("d")=d,
                       Named("U")=U_cur,
                       Named("post_z_mean")=Z_hat_cur,
                       Named("post_z_var")=diag_Sigma,
                       Named("post_z_cov")=off_diag_Sigma,
                       Named("mean_obs")=pred_cur,
                       Named("mean_obs_95lb")=pred_mean_95lb,
                       Named("mean_obs_95ub")=pred_mean_95ub,
                       Named("record_sigma0_2")=record_sigma0_2.head(m-1),
                       Named("record_rho")=record_rho.leftCols(m-1),
                       Named("record_sigma2")=record_sigma2.leftCols(m-1),
                       Named("num_iterations")=m-1);
  }
  else{
    res = List::create(Named("output") = output,
                       Named("d")=d,
                       Named("U")=U_cur,
                       Named("post_z_mean")=Z_hat_cur,
                       Named("post_z_var")=diag_Sigma,
                       Named("post_z_cov")=off_diag_Sigma,
                       Named("mean_obs")=pred_cur,
                       Named("mean_obs_95lb")=pred_mean_95lb,
                       Named("mean_obs_95ub")=pred_mean_95ub,
                       Named("sigma0_2")=sigma0_2_cur,
                       Named("rho")=rho_cur,
                       Named("sigma2")=sigma2_cur,
                       Named("num_iterations")=m-1);
  }
  if(track_neg_log_lik){
    res.push_back(record_neg_log_lik.head(m-1), "record_neg_log_lik");
  }
  return res;
}

// [[Rcpp::export]]
Eigen::MatrixXd matern_5_2_funct (const Eigen::Map<Eigen::MatrixXd> &d, double beta_i){
  //inline static Mat matern_5_2_funct (const Eigen::Map<Eigen::MatrixXd> & d, double beta_i){
  const double cnst = std::sqrt(5.0);
  Eigen::MatrixXd matOnes = Eigen::MatrixXd::Ones(d.rows(),d.cols());
  Eigen::MatrixXd result = cnst*beta_i*d;
  return ((matOnes + result +
	   result.array().pow(2.0).matrix()/3.0).cwiseProduct((-result).array().exp().matrix()));
  
}

// [[Rcpp::export]]
MatrixXd rcppeigen_get_chol(const MatrixXd &R){
  
  LLT<MatrixXd> lltOfR(R);             // compute the cholesky decomposition of R called lltofR
  MatrixXd L = lltOfR.matrixL();   //retrieve factor L  in the decomposition
  return L;
}
// [[Rcpp::export]] 
Eigen::MatrixXd F_Funct_Dev_Large_k(const Eigen::MatrixXd A_cur,const List UD){
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd return_matrix=Eigen::MatrixXd::Zero(k,d); 
  Eigen::MatrixXd UD_matrix;
  Eigen::VectorXd z;
  
  for(int i=0;i<d;i++){
    UD_matrix=UD[i];
    z=(A_cur.col(i).transpose()*UD_matrix).transpose();
    //return_matrix[,i]=UD[[i]]%*%t(z_t)
    
    return_matrix.col(i)=UD_matrix*z;
  }
  return -2*return_matrix;
}

//[[Rcpp::export]] 
List Get_B_U_V_Large_k(const Eigen::MatrixXd A_cur,const List UD){
  Eigen::MatrixXd B= F_Funct_Dev_Large_k(A_cur,UD);
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd U=Eigen::MatrixXd::Zero(k,2*d); 
  Eigen::MatrixXd V=Eigen::MatrixXd::Zero(k,2*d); 
  
  U.leftCols(d)=B;
  U.rightCols(d)=A_cur;
  
  V.leftCols(d)=A_cur;
  V.rightCols(d)=-B;
  
  List return_list;
  
  return_list.push_back(B);
  return_list.push_back(U);
  return_list.push_back(V);
  
  return return_list;
  
}

//[[Rcpp::export]] 
List Y_Funct(const Eigen::MatrixXd A_cur, const List B_U_V,  double tau){
  int d=A_cur.cols();
  int k=A_cur.rows();
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  
  //I may need to consider what if this is singular
  
  Eigen::MatrixXd middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
  
  JacobiSVD<MatrixXd> svd(middle_middle_term);
  double cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
  
  //cout << tau;    
  
  while(cond>pow(10.0,16.0)){
    // tau=tau/(2*log(cond/pow(10.0,15.0)+1));
    tau=tau/2;
    
    middle_middle_term=Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U;
    
    JacobiSVD<MatrixXd>  svd(middle_middle_term);
    
    cond = svd.singularValues()(0)/svd.singularValues()(svd.singularValues().size()-1);
    
  }
  //cout << tau;    
  
  
  Eigen::MatrixXd middle_term=U*((Eigen::MatrixXd::Identity(2*d,2*d)+tau/2.0*V.transpose()*U).lu().solve(V.transpose()));
  //Eigen::MatrixXd middle_term=U*((middle_middle_term).lu().solve(V.transpose()));
  
  Eigen::MatrixXd Y_tau= (Eigen::MatrixXd::Identity(k,k)-tau*middle_term)*A_cur;
  
  // Eigen::MatrixXd Y_tau= A_cur -tau*middle_term*A_cur;
  
  Eigen::MatrixXd Y_dev_tau=-middle_term*(A_cur+Y_tau)/2;
  
  List return_list;
  return_list.push_back(Y_tau);
  return_list.push_back(Y_dev_tau);
  return_list.push_back(tau);
  
  return return_list;
  
}

// [[Rcpp::export]] 
Eigen::MatrixXd F_Funct_Dev(const Eigen::MatrixXd A_cur,const List G){
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd return_matrix=Eigen::MatrixXd::Zero(k,d); 
  Eigen::MatrixXd G_matrix;
  
  for(int i=0;i<d;i++){
    G_matrix=G[i];
    return_matrix.col(i)=G_matrix*A_cur.col(i);
  }
  return -2*return_matrix;
}

//[[Rcpp::export]] 
List Get_B_U_V(const Eigen::MatrixXd A_cur,const List G){
  Eigen::MatrixXd B= F_Funct_Dev(A_cur,G);
  int k=A_cur.rows();
  int d=A_cur.cols();
  Eigen::MatrixXd U=Eigen::MatrixXd::Zero(k,2*d); 
  Eigen::MatrixXd V=Eigen::MatrixXd::Zero(k,2*d); 
  
  U.leftCols(d)=B;
  U.rightCols(d)=A_cur;
  
  V.leftCols(d)=A_cur;
  V.rightCols(d)=-B;
  
  List return_list;
  
  return_list.push_back(B);
  return_list.push_back(U);
  return_list.push_back(V);
  
  return return_list;
  
}

//[[Rcpp::export]] 
Eigen::MatrixXd Optimization_Stiefel_Manifold(const Eigen::MatrixXd A_ini, const List G, int max_iter){
  int k=A_ini.cols();
  int d=A_ini.rows();
  
  double rho_1=pow(10.0,-4.0);
  double delta=0.2;
  double eta=0.85;
  double epsilon_1=pow(10.0,-5.0);
  double epsilon_2=pow(10.0,-10.0);
  
  double C_cur=F_Funct(A_ini,G);
  Eigen::MatrixXd  A_cur=A_ini;
  double  Q_cur=1.0;
  //double tau_cur=0.001;
  
  //List B_U_V;
  
  List B_U_V=Get_B_U_V(A_cur, G);
  Eigen::MatrixXd B=B_U_V[0];
  Eigen::MatrixXd U=B_U_V[1];
  Eigen::MatrixXd V=B_U_V[2];
  
  //double tau_cur=0.01;
  
  
  double tau_cur=1.0/((V.transpose()*U).diagonal().array().abs().sum());
  
  //   1/sum(abs(diag(t(V)%*%U)));
  
  
  Eigen::MatrixXd gradient_F_Y_tau=B-A_cur*B.transpose()*A_cur;
  
  double F_Y_tau=F_Funct(A_cur,G);
  double F_Y_0;
  
  bool find_tau;
  Eigen::MatrixXd gradient_F_A_cur;
  double norm_gradient_cur;
  
  double F_cur_val=F_Y_tau;
  double F_last_val=F_cur_val-1;
  
  List Y_tau_dev_tau_all;
  
  Eigen::MatrixXd Y_tau;
  Eigen::MatrixXd Y_tau_dev;
  
  Eigen::MatrixXd Y_0_dev;
  
  Eigen::MatrixXd B_A;
  
  Eigen::MatrixXd diff_graident;
  double F_dev_Y_0;
  
  Eigen::MatrixXd S;
  
  double SS;
  double SY_abs;
  double YY;
  
  double tau_next;  
  
  double Q_next;
  
  int count;
  for(int i_iter=0; i_iter<max_iter;i_iter++){
    find_tau=true;
    //B_U_V=B_U_V_next;
    
    gradient_F_A_cur=gradient_F_Y_tau;
    norm_gradient_cur=pow(gradient_F_A_cur.array().pow(2.0).sum(),0.5);
    
    F_cur_val=F_Y_tau;
    
    if(i_iter>1){
      if((norm_gradient_cur/(k*d))<epsilon_1 ||  (abs(F_cur_val- F_last_val))<epsilon_2 ){
        break;
      }
    }
    
    F_last_val=F_cur_val;
    
    count=0;
    while(find_tau || count>50 ){
      count++;
      Y_tau_dev_tau_all=Y_Funct(A_cur,B_U_V,tau_cur);
      
      Y_tau=Y_tau_dev_tau_all[0];
      Y_tau_dev=Y_tau_dev_tau_all[1];
      tau_cur=Y_tau_dev_tau_all[2];
      
      
      Y_0_dev=-U*V.transpose()*A_cur;
      
      F_Y_tau=F_Funct(Y_tau,G);
      
      F_Y_0=F_Funct(A_cur,G);
      
      B_A=F_Funct_Dev(A_cur,G);
      
      F_dev_Y_0=(B_A.array()*Y_0_dev.array()).sum();
      
      
      if( (F_Y_tau- (C_cur+rho_1*tau_cur*F_dev_Y_0)) <=0){
        find_tau=false;
      }else{
        tau_cur=delta*tau_cur;
      }
    }
    
    
    //t(Y_tau)%*%Y_tau
    B_U_V=Get_B_U_V(Y_tau, G);
    
    B=B_U_V[0];
    U=B_U_V[1];
    V=B_U_V[2];
    
    //compute that trace thing for tau
    S=Y_tau-A_cur;
    
    gradient_F_Y_tau=B-Y_tau*B.transpose()*Y_tau;
    
    diff_graident=gradient_F_Y_tau-gradient_F_A_cur;
    
    SS=(S.array()*S.array()).sum();
    SY_abs=abs((S.array()*diff_graident.array()).sum());
    YY=(diff_graident.array()*diff_graident.array()).sum();
    
    if(i_iter%2==0){
      tau_next=SS/SY_abs;
    }else{
      tau_next=SY_abs/YY;
    }
    
    A_cur=Y_tau;
    Q_next=eta*Q_cur+1;
    C_cur=(eta*Q_cur*C_cur+F_Y_tau)/Q_next;
    Q_cur=Q_next;
    
    //this is controversial   
    if(!isnan(tau_next)){
      tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
    }
    
    //tau_cur=max(min(tau_next,pow(10.0,20)),pow(10.0,-20));
    
  }
  
  return A_cur; 
  
}

///// particle interaction estimation
// [[Rcpp::export]] 
VectorXd A_t_times_x_particle(const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi num_neighbors_vec,  
                              const int D_y,const int N_tilde){
  
  //int N=n*T_sim*S;
  //int N=N_tilde/D;
  //int ND=N*D; //dimension of the output vector
  //int N_d=output.size();
  VectorXd vec_N_tilde = VectorXd::Zero(N_tilde);
  
  
  int j=0;
  int count_ne_here=0;
  
  
  
  for(int i=0; i<N_tilde;i++){
    vec_N_tilde[i]=(A_all_v.segment(i*D_y,D_y).transpose()*output.segment(j*D_y,D_y)).value(); //[(i-1)*D+1:D]*p[j*D+1:D]
    count_ne_here=count_ne_here+1;
    if( (count_ne_here)%num_neighbors_vec[j]==0 ){
      j=j+1;
      count_ne_here=0;
    }
  }
  return vec_N_tilde;
  
}


// [[Rcpp::export]] 
VectorXd A_times_x_particle(const VectorXd output, const Eigen::VectorXd A_all_v, const VectorXi num_neighbors_vec,  
                            const int D, const int N){
  
  //int N=n*T_sim*S;
  int n_t=N/D;
  //int ND=N*D; //dimension of the output vector
  int N_d=output.size();
  //int N_d=sort_d_all_ix.size();
  
  
  
  VectorXd vec_N_d = VectorXd::Zero(N_d);
  
  VectorXd vec_N=VectorXd::Zero(N);
  
  int sum_neighbors_start=0;
  
  VectorXd vec;
  
  int num_neighbors_vec_i_obs;
  
  for(int i_obs=0;i_obs<n_t;i_obs++){
    //index_v=sum_neighbors_start+seq(1,num_neighbors_vec[i_obs],1)-1; 
    //index_A=D*sum_neighbors_start+seq(1,num_neighbors_vec[i_obs]*D,1)-1;  
    num_neighbors_vec_i_obs=num_neighbors_vec[i_obs];
    vec=A_all_v.segment(D*sum_neighbors_start,num_neighbors_vec_i_obs*D);
    Map<MatrixXd>   A_here(vec.data(),D,num_neighbors_vec_i_obs);
    //A_here(mat.data(),D,num_neighbors_vec_i_obs);//matrix(A_all_v.segment(D*sum_neighbors_start,num_neighbors_vec[i_obs]*D),D,num_neighbors_vec)
    //for(int i_D=0;i_D<D;i_D++){
    //vec_N[i_obs*D+i_D]=(vec.segment(i_D*num_neighbors_vec_i_obs,num_neighbors_vec_i_obs).transpose()*vec_N_d.segment(sum_neighbors_start,num_neighbors_vec_i_obs)).value(); // [index_v]
    vec_N.segment(i_obs*D,D)=A_here*output.segment(sum_neighbors_start,num_neighbors_vec_i_obs); // [index_v]
    sum_neighbors_start=sum_neighbors_start+num_neighbors_vec[i_obs];
    //}
  }
  
  return vec_N;
}

// [[Rcpp::export]]
List IKF_CG_particle(VectorXd param, const  String kernel_type, const VectorXd delta_x_all,
                     const VectorXd output, const Eigen::VectorXd A_all_v, 
                     const VectorXi sort_d_all_ix,
                     const VectorXi num_neighbors_vec, const double tilde_nu,
                     const int D, const int N, 
                     float tol=1e-6, int maxIte = 1000){
  
  
  //int N=n*T_sim*S;
  int n=N/D;
  //int ND=N*D; //dimension of the output vector
  
  //int N_d=S*n*n*T_sim;
  int N_d=sort_d_all_ix.size();
  //double beta=exp(param[0]);
  //here I use tau as sigma_j_2/sigma_2_0 in the paper 
  //double  nu=exp(param[1]);
  double  tau=exp(param[1]);
  
  
  
  
  VectorXd x=VectorXd::Zero(N);
  VectorXd r = output;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1.0;
  double rs_ratio;
  
  VectorXd vec_N_d = VectorXd::Zero(N_d);
  VectorXd vec_N_d_sorted=VectorXd::Zero(N_d);
  
  VectorXd vec_N=VectorXd::Zero(N);
  VectorXd y_KF;
  
  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);
  
  VectorXd vec;
  
  int num_neighbors_vec_i_obs;
  double alpha;
  //preconstructure 
  
  double gamma=1.0/exp(param[0]);
  
  //Oct 24, 2022
  //double VV=exp(param[1]);
  //this should be tilde_nu, the stablizing par, not nu
  double VV=tilde_nu; //here tilde_nu is the tilde_sigma_2_0 to stablize the computation
  
  
  
  Eigen::MatrixXd    W0;
  List    GG;
  
  List    W;
  List    Q_K;
  
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    
    //param.tail(k).array().exp().matrix();
    
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x_all,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x_all,lambda,W0);
    
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x_all,lambda);  
    W=Construct_W_exp(1.0,delta_x_all,lambda,W0);
    
  }
  
  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  //Sep 2022
  Q = (Q.array() < 0).select(0, Q); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  
  //Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  //end 
  while((ite < maxIte) && (rs_new > tol)){
    
    int j=0;
    int count_ne_here=0;
    
    for(int i=0; i<N_d;i++){
      vec_N_d[i]=(A_all_v.segment(i*D,D).transpose()*p.segment(j*D,D)).value(); //[(i-1)*D+1:D]*p[j*D+1:D]
      count_ne_here=count_ne_here+1;
      if( (count_ne_here)%num_neighbors_vec[j]==0 ){
        j=j+1;
        count_ne_here=0;
      }
    }
    
    for(int i=0; i<N_d;i++){
      vec_N_d_sorted[i]=vec_N_d[sort_d_all_ix[i]-1];
    }
    
    y_KF=Get_R_y(GG, Q, K, vec_N_d_sorted)-tilde_nu*vec_N_d_sorted;  //remove the added  tilde_nu*z
    
    for(int i=0; i<N_d;i++){
      vec_N_d[sort_d_all_ix[i]-1]=y_KF[i];
    }
    
    
    int sum_neighbors_start=0;
    
    for(int i_obs=0;i_obs<n;i_obs++){
      num_neighbors_vec_i_obs=num_neighbors_vec[i_obs];
      vec=A_all_v.segment(D*sum_neighbors_start,num_neighbors_vec_i_obs*D);
      Map<MatrixXd>   A_here(vec.data(),D,num_neighbors_vec_i_obs);
      vec_N.segment(i_obs*D,D)=A_here*vec_N_d.segment(sum_neighbors_start,num_neighbors_vec_i_obs); // [index_v]
      sum_neighbors_start=sum_neighbors_start+num_neighbors_vec[i_obs];
      
    }
    //vec_N=vec_N+nu*p;
    vec_N=vec_N*tau+p;
    
    
    
    alpha = rs_old / (p.transpose() * vec_N).value();
    x += alpha*p;
    r -= alpha*vec_N;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;
    
  }
  
  
  List ans_list;
  
  ans_list.push_back(x);
  ans_list.push_back(resid);
  ans_list.push_back(ite);
  
  
  return ans_list;
}

// [[Rcpp::export]]
List IKF_CG_particle_two_interact(VectorXd param1, VectorXd param2, 
                                  const String kernel_type1, const String kernel_type2,
                                  const VectorXd delta_x_all1,const VectorXd delta_x_all2,
                                  const Eigen::VectorXd A_all_v1, const Eigen::VectorXd A_all_v2,
                                  const VectorXi sort_d_all_ix1,const VectorXi sort_d_all_ix2,
                                  const VectorXi num_neighbors_vec1,const VectorXi num_neighbors_vec2,
                                  const VectorXd output, const double tilde_nu,const int D, const int N, 
                                  float tol=1e-6, int maxIte = 1000){
  
  
  //int N=n*T_sim*S;
  int n=N/D;
  //int ND=N*D; //dimension of the output vector
  
  //int N_d=S*n*n*T_sim;
  int N_d1=sort_d_all_ix1.size();
  int N_d2=sort_d_all_ix2.size();
  //double beta=exp(param[0]);
  //here I use tau as sigma_j_2/sigma_2_0 in the paper 
  //double  nu=exp(param[1]);
  double  tau1=exp(param1[1]);
  double  tau2=exp(param2[1]);
  
  
  
  
  VectorXd x=VectorXd::Zero(N);
  VectorXd r = output;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1.0;
  double rs_ratio;
  
  VectorXd vec_N_d1 = VectorXd::Zero(N_d1);
  VectorXd vec_N_d2 = VectorXd::Zero(N_d2);
  VectorXd vec_N_d_sorted1=VectorXd::Zero(N_d1);
  VectorXd vec_N_d_sorted2=VectorXd::Zero(N_d2);
  
  VectorXd vec_N1=VectorXd::Zero(N);
  VectorXd vec_N2=VectorXd::Zero(N);
  VectorXd y_KF1;
  VectorXd y_KF2;
  
  VectorXd Bp = VectorXd::Zero(N);
  
  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);
  
  VectorXd vec;
  
  int num_neighbors_vec_i_obs;
  double alpha;
  //preconstructure 
  
  double gamma1=1.0/exp(param1[0]);
  double gamma2=1.0/exp(param2[0]);
  
  //Oct 24, 2022
  //double VV=exp(param[1]);
  //this should be tilde_nu, the stablizing par, not nu
  double VV=tilde_nu; //here tilde_nu is the tilde_sigma_2_0 to stablize the computation
  
  
  
  Eigen::MatrixXd    W01;
  Eigen::MatrixXd    W02;
  List    GG1;
  List    GG2;
  
  List    W1;
  List    W2;
  List    Q_K1;
  List    Q_K2;
  
  
  double lambda=0;
  if(kernel_type1=="matern_5_2"){
    lambda=sqrt(5.0)/gamma1;
    
    //param.tail(k).array().exp().matrix();
    
    W01=Construct_W0_matern_5_2(1.0,lambda);   
    GG1=Construct_G_matern_5_2(delta_x_all1,lambda);  
    W1=Construct_W_matern_5_2(1.0,delta_x_all1,lambda,W01);
    
  }else if(kernel_type1=="exp"){
    
    lambda=1.0/gamma1;
    W01=Construct_W0_exp(1.0,lambda);  
    GG1=Construct_G_exp(delta_x_all1,lambda);  
    W1=Construct_W_exp(1.0,delta_x_all1,lambda,W01);
    
  }
  Q_K1=Get_Q_K(GG1,W1,W01,VV);
  
  if(kernel_type2=="matern_5_2"){
    lambda=sqrt(5.0)/gamma2;
    
    //param.tail(k).array().exp().matrix();
    
    W02=Construct_W0_matern_5_2(1.0,lambda);   
    GG2=Construct_G_matern_5_2(delta_x_all2,lambda);  
    W2=Construct_W_matern_5_2(1.0,delta_x_all2,lambda,W02);
    
  }else if(kernel_type2=="exp"){
    
    lambda=1.0/gamma2;
    W02=Construct_W0_exp(1.0,lambda);  
    GG2=Construct_G_exp(delta_x_all2,lambda);  
    W2=Construct_W_exp(1.0,delta_x_all2,lambda,W02);
    
  }
  
  Q_K2=Get_Q_K(GG2,W2,W02,VV);
  
  
  Eigen::VectorXd Q1=Q_K1[0];
  Eigen::VectorXd Q2=Q_K2[0];
  Eigen::MatrixXd K1=Q_K1[1];
  Eigen::MatrixXd K2=Q_K2[1];
  
  //Sep 2022
  Q1 = (Q1.array() < 0).select(0, Q1); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  Q2 = (Q2.array() < 0).select(0, Q2);
  
  //Eigen::VectorXd sqrt_Q1=Q1.array().sqrt();
  //Eigen::VectorXd sqrt_Q2=Q2.array().sqrt();
  //end 
  while((ite < maxIte) && (rs_new > tol)){
    //first kernel
    int j=0;
    int count_ne_here=0;
    
    for(int i=0; i<N_d1;i++){
      vec_N_d1[i]=(A_all_v1.segment(i*D,D).transpose()*p.segment(j*D,D)).value(); //[(i-1)*D+1:D]*p[j*D+1:D]
      count_ne_here=count_ne_here+1;
      if( (count_ne_here)%num_neighbors_vec1[j]==0 ){
        j=j+1;
        count_ne_here=0;
      }
    }
    
    for(int i=0; i<N_d1;i++){
      vec_N_d_sorted1[i]=vec_N_d1[sort_d_all_ix1[i]-1];
    }
    
    
    y_KF1=Get_R_y(GG1, Q1, K1, vec_N_d_sorted1)-tilde_nu*vec_N_d_sorted1;  //remove the added  tilde_nu*z
    
    for(int i=0; i<N_d1;i++){
      vec_N_d1[sort_d_all_ix1[i]-1]=y_KF1[i];
    }
    
    
    int sum_neighbors_start=0;
    
    for(int i_obs=0;i_obs<n;i_obs++){
      num_neighbors_vec_i_obs=num_neighbors_vec1[i_obs];
      vec=A_all_v1.segment(D*sum_neighbors_start,num_neighbors_vec_i_obs*D);
      Map<MatrixXd>   A_here(vec.data(),D,num_neighbors_vec_i_obs);
      vec_N1.segment(i_obs*D,D)=A_here*vec_N_d1.segment(sum_neighbors_start,num_neighbors_vec_i_obs); // [index_v]
      sum_neighbors_start=sum_neighbors_start+num_neighbors_vec1[i_obs];
      
    }
    
    //second kernel
    j=0;
    count_ne_here=0;
    
    for(int i=0; i<N_d2;i++){
      vec_N_d2[i]=(A_all_v2.segment(i*D,D).transpose()*p.segment(j*D,D)).value(); //[(i-1)*D+1:D]*p[j*D+1:D]
      count_ne_here=count_ne_here+1;
      if( (count_ne_here)%num_neighbors_vec2[j]==0 ){
        j=j+1;
        count_ne_here=0;
      }
    }
    
    for(int i=0; i<N_d2;i++){
      vec_N_d_sorted2[i]=vec_N_d2[sort_d_all_ix2[i]-1];
    }
    
    y_KF2=Get_R_y(GG2, Q2, K2, vec_N_d_sorted2)-tilde_nu*vec_N_d_sorted2;  //remove the added  tilde_nu*z
    
    for(int i=0; i<N_d2;i++){
      vec_N_d2[sort_d_all_ix2[i]-1]=y_KF2[i];
    }
    
    
    sum_neighbors_start=0;
    
    for(int i_obs=0;i_obs<n;i_obs++){
      num_neighbors_vec_i_obs=num_neighbors_vec2[i_obs];
      vec=A_all_v2.segment(D*sum_neighbors_start,num_neighbors_vec_i_obs*D);
      Map<MatrixXd>   A_here(vec.data(),D,num_neighbors_vec_i_obs);
      vec_N2.segment(i_obs*D,D)=A_here*vec_N_d2.segment(sum_neighbors_start,num_neighbors_vec_i_obs); // [index_v]
      sum_neighbors_start=sum_neighbors_start+num_neighbors_vec2[i_obs];
      
    }
    
    
    //vec_N=vec_N+nu*p;
    Bp=vec_N1*tau1+vec_N2*tau2+p;
    
    
    alpha = rs_old / (p.transpose() * Bp).value();
    x += alpha*p;
    r -= alpha*Bp;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;
    
  }
  
  
  List ans_list;
  
  ans_list.push_back(x);
  ans_list.push_back(resid);
  ans_list.push_back(ite);
  
  
  return ans_list;
}


// [[Rcpp::export]]
List IKF_CG_particle_cell(VectorXd param, const  String kernel_type, const VectorXd delta_x_all,
                              const VectorXd output, const Eigen::VectorXd A_all_v, 
                              const VectorXi sort_d_all_ix, const VectorXd sigma_2_vec,
                              const VectorXi num_neighbors_vec, const double tilde_nu,
                              const int D, const VectorXd n_t_record, //const int N, 
                              float tol=1e-6, int maxIte = 1000){
  
  
  //int n=N/D;
  int n=n_t_record.sum();
  int N=n*D;
  int T_time=n_t_record.size();
  //VectorXd Diag=VectorXd::Zero(N);
  VectorXd Diag=VectorXd::Zero(N);
  
  for(int i=0; i<T_time;i++){
    if(i==0){
      Diag.segment(0,n_t_record(i))=VectorXd::Constant(1,n_t_record(i),sigma_2_vec(i));
    }else{
      Diag.segment(n_t_record.head(i).sum(),n_t_record(i))=VectorXd::Constant(1,n_t_record(i),sigma_2_vec(i));
    }
  }
  
  //int N_d=S*n*n*T_sim;
  int N_d=sort_d_all_ix.size();
  //double beta=exp(param[0]);
  //here I use tau as sigma_j_2/sigma_2_0 in the paper 
  //double  nu=exp(param[1]);
  double  tau=exp(param[1]);
  
  
  VectorXd x=VectorXd::Zero(N);
  VectorXd r = output;
  VectorXd p = r;
  double rs_old = (r.transpose() * r).value();
  double rs_new=1.0;
  double rs_ratio;
  
  VectorXd vec_N_d = VectorXd::Zero(N_d);
  VectorXd vec_N_d_sorted=VectorXd::Zero(N_d);
  
  VectorXd vec_N=VectorXd::Zero(N);
  VectorXd y_KF;
  
  int ite = 0;
  VectorXd resid = VectorXd::Zero(maxIte);
  
  VectorXd vec;
  
  int num_neighbors_vec_i_obs;
  double alpha;
  //preconstructure 
  
  double gamma=1.0/exp(param[0]);
  
  //Oct 24, 2022
  //double VV=exp(param[1]);
  //this should be tilde_nu, the stablizing par, not nu
  double VV=tilde_nu; //here tilde_nu is the tilde_sigma_2_0 to stablize the computation
  
  
  
  Eigen::MatrixXd    W0;
  List    GG;
  
  List    W;
  List    Q_K;
  
  
  double lambda=0;
  if(kernel_type=="matern_5_2"){
    lambda=sqrt(5.0)/gamma;
    
    //param.tail(k).array().exp().matrix();
    
    W0=Construct_W0_matern_5_2(1.0,lambda);   
    GG=Construct_G_matern_5_2(delta_x_all,lambda);  
    W=Construct_W_matern_5_2(1.0,delta_x_all,lambda,W0);
    
  }else if(kernel_type=="exp"){
    
    lambda=1.0/gamma;
    W0=Construct_W0_exp(1.0,lambda);  
    GG=Construct_G_exp(delta_x_all,lambda);  
    W=Construct_W_exp(1.0,delta_x_all,lambda,W0);
    
  }
  
  Q_K=Get_Q_K(GG,W,W0,VV);
  
  
  Eigen::VectorXd Q=Q_K[0];
  Eigen::MatrixXd K=Q_K[1];
  
  //Sep 2022
  Q = (Q.array() < 0).select(0, Q); //truncate them to be zero to avoid NA in singular case, but error could be large if it is too singular
  
  //Eigen::VectorXd sqrt_Q=Q.array().sqrt();
  //end 
  while((ite < maxIte) && (rs_new > tol)){
    
    int j=0;
    int count_ne_here=0;
    
    for(int i=0; i<N_d;i++){
      vec_N_d[i]=(A_all_v.segment(i*D,D).transpose()*p.segment(j*D,D)).value(); //[(i-1)*D+1:D]*p[j*D+1:D]
      count_ne_here=count_ne_here+1;
      if( (count_ne_here)%num_neighbors_vec[j]==0 ){
        j=j+1;
        count_ne_here=0;
      }
    }
    
    for(int i=0; i<N_d;i++){
      vec_N_d_sorted[i]=vec_N_d[sort_d_all_ix[i]-1];
    }
    
    y_KF=Get_R_y(GG, Q, K, vec_N_d_sorted)-tilde_nu*vec_N_d_sorted;  //remove the added  tilde_nu*z
    
    for(int i=0; i<N_d;i++){
      vec_N_d[sort_d_all_ix[i]-1]=y_KF[i];
    }
    
    
    int sum_neighbors_start=0;
    
    for(int i_obs=0;i_obs<n;i_obs++){
      num_neighbors_vec_i_obs=num_neighbors_vec[i_obs];
      vec=A_all_v.segment(D*sum_neighbors_start,num_neighbors_vec_i_obs*D);
      Map<MatrixXd>   A_here(vec.data(),D,num_neighbors_vec_i_obs);
      vec_N.segment(i_obs*D,D)=A_here*vec_N_d.segment(sum_neighbors_start,num_neighbors_vec_i_obs); // [index_v]
      sum_neighbors_start=sum_neighbors_start+num_neighbors_vec[i_obs];
      
    }
    //vec_N=vec_N+nu*p;
    //vec_N=vec_N*tau+p;
    vec_N=vec_N*tau+p.VectorXd::cwiseProduct(Diag);
    
    
    
    alpha = rs_old / (p.transpose() * vec_N).value();
    x += alpha*p;
    r -= alpha*vec_N;
    rs_new = (r.transpose() * r).value();
    rs_ratio = rs_new / rs_old;
    p = r + rs_ratio * p;
    rs_old = rs_new;
    resid[ite] = rs_new; //mean((res2[[1]]-z)^2)
    ite++;
    
  }
  
  
  List ans_list;
  
  ans_list.push_back(x);
  ans_list.push_back(resid);
  ans_list.push_back(ite);
  
  
  return ans_list;
}




