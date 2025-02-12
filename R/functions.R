##########################################################################
## fgasp construction function
## 
## FastGaSP Package
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, April 2013.
##
## Copyright (C) 2018-present Mengyang Gu
##							  
##    
##########################################################################


fgasp <- function(input, output, have_noise=TRUE, kernel_type='matern_5_2'){


  if( (kernel_type!='matern_5_2')&(kernel_type!='exp') ){
    stop("the current version only support the Matern covariance with smoothness 
         parameter being 2.5 or 0.5 (exponential kernel). \n")
  }
  
  ##if the input is not sorted so we will sort them
  object <- new("fgasp")
  object@num_obs=length(output)
  #object@param=param
  object@kernel_type=kernel_type
  object@have_noise=have_noise
  
  ## if the input is not sorted, then I sort it
  if(sum(which(input[2:object@num_obs]-input[1:(object@num_obs-1)]<0 ))>0){
    input_sorted_all=sort(input,index.return=TRUE)
    object@input=input_sorted_all[[1]]
    object@output=output[input_sorted_all[[2]]]
  }else{
    object@input=input
    object@output=output
  }
  object@delta_x=object@input[2:object@num_obs]-object@input[1:(object@num_obs-1)]
  
  if(length(object@delta_x)!=(length(object@output)-1)){
    stop("length(output) should be equal to length(delta_x)+1 \n")
  }
  ##the current version does not support duplicated inputs
  if(length(which(object@delta_x==0))>0){
    stop("Please delete the repeated inputs. \n")
  }
  return(object)
}

show.fgasp <- function(object) {	
  cat('Number of observations: ',object@num_obs,'\n')
  cat('Have noise: ',object@have_noise,'\n')
  cat('Type of kernel: ',object@kernel_type,'\n')
  
}


##likelihood function of GaSP using fast computation
##the Filtering algorithm is implemented using Rcpp and RcppEigen
log_lik<-function(param,object){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
       stop("Please specify both the log inverse range parameter 
            and log nugget parameter if there is noise. \n")
  }
  
  log_det_S2=Get_log_det_S2(param,object@have_noise,object@delta_x,
                            object@output, object@kernel_type);
  ##log likelihood
  -log_det_S2[[1]]/2-(object@num_obs)/2*log(log_det_S2[[2]])
}

##prediction 
predict.fgasp<-function(param,object,testing_input,var_data=TRUE,sigma_2=NULL){
  param=as.vector(param)
  if(object@have_noise==T && length(param)!=2){
    stop("Please specify both the log inverse range parameter 
          and log nugget parameter. \n")
  }
  ## if sigma_2 is not given, I will estimate it
  if(is.null(sigma_2)){
  sigma_2=Get_log_det_S2(param,object@have_noise,object@delta_x,
                             object@output,object@kernel_type)[[2]]/object@num_obs
  }
  
  predictobj<- new("predictobj.fgasp")
  
  testing_input_sorted=sort(testing_input)
  predictobj@num_testing=length(testing_input_sorted)
  predictobj@testing_input=testing_input_sorted
  predictobj@param=param
  predictobj@var_data=var_data
  ## we don't support the repeated testing input in this version
  if(length(which(testing_input_sorted[2:predictobj@num_testing]-
                  testing_input_sorted[1:(predictobj@num_testing-1)]==0))>0){
    stop("Please delete the repeated testing inputs. \n")
  }
  ##combine the input and testing input
  input_all_ori=c(object@input,testing_input_sorted)
  input_all_sort=sort(input_all_ori,index.return=TRUE)

  input_all=input_all_sort[[1]] ##all sorted input with and without the observations
  delta_x_all=input_all[2:length(input_all_ori)]-input_all[1:(length(input_all_ori)-1)]
  
  ##create a vector of inputs where the one mean there is an observation
  ##and zero means no observation
  index_all=rep(0, length(input_all))
  index_all[1:object@num_obs]=1
  index_all=index_all[input_all_sort[[2]]]
  index_all=as.integer(index_all)
  
  

  testing_loc=which(index_all==0)
  delta_x_all_here=delta_x_all
  
  if(length(which(delta_x_all==0))>0){
    #input_all=input_all[-which(delta_x_all==0)]
    index_all=index_all[-(which(delta_x_all==0)+1)]
    delta_x_all_here=delta_x_all[-which(delta_x_all==0)]
  }
  
  
  
  KF_smoother_result=Kalman_smoother(param, object@have_noise,
                                     index_all, delta_x_all_here,object@output,sigma_2,object@kernel_type)
  if(length(which(delta_x_all==0))==0){
    predictobj@mean=KF_smoother_result[[1]][testing_loc]
    if(var_data==T | object@have_noise==F){
       predictobj@var=KF_smoother_result[[2]][testing_loc]
    }else{
      predictobj@var=KF_smoother_result[[2]][testing_loc]-sigma_2*exp(param[2])
    }
  }else{##enlarge it to the orginal long sequence that may contain knots. 
    res=rep(NA,length(input_all_ori))
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[1]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    predictobj@mean=res[testing_loc]
    
    res[-(which(delta_x_all==0)+1)]=KF_smoother_result[[2]]
    res[which(delta_x_all==0)+1]= res[which(delta_x_all==0)]
    if(var_data==T | object@have_noise==F){
      predictobj@var=res[testing_loc]
    }else{
      predictobj@var=res[testing_loc]-sigma_2*exp(param[2])
    }
    
  }
  
  return(predictobj)
  
}

###simulate
#simulate.fgasp <- function (param,object,...){
  
#}

## fmou class
fmou <- function(output, d=NULL, est_d=FALSE,est_U0=TRUE, est_sigma0_2=TRUE, 
                 U0 = NULL, sigma0_2=NULL){
  object <- new("fmou", output, d, est_d, est_U0, est_sigma0_2, 
                U0, sigma0_2)
  return(object)
}

#### fmou k-step-ahead prediction
predict_fmou <- function(param_est, step=1, interval=FALSE, interval_data=TRUE){
  n = ncol(param_est$output)
  U0 = param_est$U
  rho = param_est$rho
  z_n_hat = matrix(param_est$post_z_mean[,n], ncol=1)
  z_n_var = matrix(param_est$post_z_var[,n], ncol=1)
  sigma2 = param_est$sigma2
  sigma0_2 = param_est$sigma0_2
  
  pred_mean_obs = matrix(NA, nrow=nrow(U0), ncol=length(step))
  for(t in step){
    pred_mean_latent = diag(rho^{t}) %*% z_n_hat
    pred_mean_obs[,t] = U0 %*% pred_mean_latent
  }
  
  if(interval){
    pred_z_var = matrix(NA, nrow=length(rho), ncol=length(step))
    pre_var = z_n_var
    for(t in step){
      pred_z_var[,t] = rho^2 * pre_var + sigma2
      pre_var = pred_z_var[,t]
    }
    
    if(interval_data){
      data_var = matrix(NA, nrow=nrow(U0), ncol=length(step))
      for(t in step){
        data_var[,t] = rowSums(U0*t(matrix(pred_z_var[,t],dim(U0)[2],dim(U0)[1]))*U0) + sigma0_2
      }
      pred_interval_95lb = pred_mean_obs - 1.96 * sqrt(data_var)
      pred_interval_95ub = pred_mean_obs + 1.96 * sqrt(data_var)
    }
    else{
      mean_var = matrix(NA, nrow=nrow(U0), ncol=length(step))
      for(t in step){
        mean_var[,t] = rowSums(U0*t(matrix(pred_z_var[,t],dim(U0)[2],dim(U0)[1]))*U0)
      }
      pred_interval_95lb = pred_mean_obs - 1.96 * sqrt(mean_var)
      pred_interval_95ub = pred_mean_obs + 1.96 * sqrt(mean_var)
    }
  }
  
  res = list(pred_mean = pred_mean_obs)
  if(interval){
    res$pred_interval_95lb = pred_interval_95lb
    res$pred_interval_95ub = pred_interval_95ub
  }
  return(res)
}

## gppca class
neg_log_lik_shared_cov_FFBS<-function(param,kernel_type, output, delta_x,d){
  # compute negative likelihood when latent processes share parameters
  
  n = ncol(output)
  k = nrow(output)
  output_2=output%*%t(output)
  trace_output_2=sum(output^2)
  G_log_det_cov=Get_G_log_det_cov(param, output, delta_x,d=1,kernel_type=kernel_type)
  G=output_2-G_log_det_cov[[1]]
  
  eigen_G=eigen(G)
  
  res = -(-(sum(G_log_det_cov[[2]]))*d/2-(n*k)/2*log(trace_output_2-sum(eigen_G$values[1:d]) ))
  return(res)
}

neg_log_lik_diff_cov_FFBS<-function(param,A_ini,kernel_type='matern_5_2',output,delta_x){
  # compute negative likelihood when latent processes have distinct parameters
  
  d = ncol(A_ini)
  n = ncol(output)
  k = nrow(output)
  G_log_det_cov=Get_G_log_det_cov(param,output, delta_x,d,kernel_type)
  G=as.list(1:d)
  output_2=output%*%t(output)
  trace_output_2=sum(output^2)
  
  for(i in 1:d){
    G[[i]]=output_2-G_log_det_cov[[i]]
  }
  
  A_hat_here=Optimization_Stiefel_Manifold(A_ini, G=G,max_iter=200)
  
  S_2=trace_output_2+F_Funct(A_hat_here,G)
  neg_log_lik=-(-(sum(G_log_det_cov[[d+1]]))/2-(n*k)/2*log(S_2))
  return(neg_log_lik)
}

pred_FFBS_FastGaSP<-function(param,A_hat,input,testing_input, output_here,d,var_data=F,kernel_type){
  # predict posterior mean of observations
  
  beta=param[1]
  sigma_2=param[2]
  sigma_2_0=param[3]
  k = nrow(output_here)
  
  num_testing=length(testing_input)
  num_obs=length(input)
  
  predict_all_dlm=matrix(0, num_testing,d)
  var_all_dlm=matrix(0, num_testing,d)
  output_t_A=t(output_here)%*%A_hat
  var_est=0
  
  for(i_d in 1:d){
    
    m.model=fgasp(input,output_t_A[,i_d],have_noise=TRUE,kernel_type=kernel_type)
    m.pred=predict(param=c(log(beta),log(sigma_2_0/sigma_2)),object=m.model,
                   testing_input=testing_input,var_data=FALSE,sigma_2=sigma_2)
    predict_all_dlm[,i_d]=m.pred@mean
    
    var_all_dlm[,i_d]= m.pred@var  
    
    var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
  }
  
  return.list=as.list(1:2)
  return.list[[1]]=A_hat%*%t(predict_all_dlm)
  if(var_data==F){
    return.list[[2]]=var_est
  }else{
    return.list[[2]]=var_est+rep(sigma_2_0,k)
  }
  
  return(return.list)
  
}

pred_FFBS_FastGaSP_diff_cov<-function(param,A_hat,input,testing_input, output_here,d,var_data=F,kernel_type='matern_5_2'){
  
  
  beta=param[1:d]
  sigma_2=param[(d+1):(2*d)]
  sigma_2_0=param[2*d+1]
  num_testing=length(testing_input)
  num_obs=length(input)
  k = nrow(output_here)
  
  predict_all_dlm=matrix(0, num_testing,d)
  var_all_dlm=matrix(0, num_testing,d)
  output_t_A=t(output_here)%*%A_hat
  var_est=0
  
  for(i_d in 1:d){
    
    m.model=fgasp(input,output_t_A[,i_d],have_noise=TRUE,kernel_type=kernel_type)
    m.pred=predict(param=c(log(beta[i_d]),log(sigma_2_0/sigma_2[i_d])),object=m.model,
                   testing_input=testing_input,var_data=FALSE,sigma_2=sigma_2[i_d])
    predict_all_dlm[,i_d]=m.pred@mean
    
    var_all_dlm[,i_d]= m.pred@var  
    
    var_est= var_est+(A_hat[,i_d]^2)%*%t(as.matrix(var_all_dlm[,i_d],num_testing,1))
  }
  
  return.list=as.list(1:2)
  return.list[[1]]=A_hat%*%t(predict_all_dlm)
  if(var_data==F){
    return.list[[2]]=var_est
  }else{
    return.list[[2]]=var_est+rep(sigma_2_0,k)
  }
  
  return.list
  
}

get_chol<-function(x,beta){
  # get cholesky decomposition
  R0_00=abs(outer(x,x,'-'))
  R=matern_5_2_funct(R0_00,beta)
  rcppeigen_get_chol(R)
}

gppca <- function(input, output, d,est_d=FALSE, shared_params = TRUE, kernel_type="matern_5_2"){
  object <- new("gppca", input, output, d, est_d, shared_params, kernel_type)
  return(object)
}

gppca_fit <- function(input, output, d, est_d, shared_params, kernel_type, sigma0_2, d_ub){
  input=sort(input)
  delta_x=input[2:length(input)]-input[1:(length(input)-1)]
  
  k = nrow(output)
  n = ncol(output)
  output_2=output%*%t(output)
  trace_output_2=sum(output^2)
  svd_output=svd(output)
  
  if(est_d){
    if(is.null(d_ub)){
      d_ub = k
    }
    if(is.null(sigma0_2)){
      # variance of noise if unknown, use IC to estimate d
      IC_vals = rep(NA, d_ub)
      for(i_d in 1:d_ub){
        criteria_val = log(mean((output - svd_output$u[,1:i_d]%*%t(svd_output$u[,1:i_d])%*%output)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
        IC_vals[i_d] = criteria_val
      }
      d = which.min(IC_vals)
      print(d)
      
    }
    else{
      # variance is known, use binary search
      left = 1
      right = d_ub
      
      while(left <= right){
        mid = left + floor((right-left)/2)
        if(shared_params){
          param_ini=c(log(.1),log(10))
          m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type,output=output,
                      delta_x=delta_x, d=mid,method="L-BFGS-B"),silent=T)
          
          while(is.character(m)){
            param_ini=param_ini+runif(2)
            m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS,kernel_type=kernel_type, output=output,
                        delta_x=delta_x,d=mid, method="L-BFGS-B"),silent=T)
          }
          G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d=1,kernel_type=kernel_type)
          G_hat=output_2- G_log_det_cov_hat[[1]]
          eigen_G_hat=eigen(G_hat)
          sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:mid]))/(n*k)
        }
        else{
          param_ini=c(rep(log(.1),mid),rep(log(10),mid))
          m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,output=output,delta_x=delta_x,
                      lower=c(rep(log(10^{-10}),mid),rep(log(10^{-10}),mid)),
                      upper=c(rep(log(10^10),mid),rep(log(10^10),mid)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
          
          while(is.character(m)){
            param_ini=param_ini+runif(2*mid)
            #print(is.character(m))
            m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,output=output,delta_x=delta_x,
                        lower=c(rep(log(10^{-10}),mid),rep(log(10^{-10}),mid)),
                        upper=c(rep(log(10^10),mid),rep(log(10^10),mid)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
          }
          
          tau_hat=exp(m$par[(mid+1):(2*mid)])
          G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d,kernel_type = kernel_type)
          for(i in 1:mid){
            G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
          }
          A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat,max_iter=100)
          sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(n*k)
        }
        if(sigma_2_0_hat < sigma0_2){
          # overfit, reduce d
          right = mid-1
        }
        else{
          # underfit, increase d
          left = mid + 1
        }
      }
      d = mid
    }
    
  }
  
  ## start to implement GPPCA, shared_params=TRUE
  A_ini=svd_output$u[,1:d] 
  if(shared_params){
    param_ini=c(log(.1),log(10))
    m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS, kernel_type=kernel_type,output=output,
                delta_x=delta_x, d=d, method="L-BFGS-B"),silent=T)
    
    while(is.character(m)){
      #print(is.character(m))
      param_ini=param_ini+runif(2)
      m=try(optim(param_ini,neg_log_lik_shared_cov_FFBS,kernel_type=kernel_type, output=output,
                  delta_x=delta_x,d=d, method="L-BFGS-B"),silent=T)
    } 
    
    beta_hat=exp(m$par[1])
    tau_hat=exp(m$par[2])
    
    G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d=1,kernel_type=kernel_type)
    G_hat=output_2- G_log_det_cov_hat[[1]]
    eigen_G_hat=eigen(G_hat)
    sigma_2_0_hat=(trace_output_2-sum(eigen_G_hat$values[1:d]))/(n*k)
    sigma2_hat = tau_hat*(sigma_2_0_hat)
    est_param = c(beta_hat,sigma2_hat, sigma_2_0_hat)
    
    
    AZ_posterior=pred_FFBS_FastGaSP(param=est_param,A_hat=eigen_G_hat$vectors[,1:d],input=input,testing_input=input, output_here=output,d=d,kernel_type=kernel_type)
    Y_hat=AZ_posterior[[1]]
    Y_var=AZ_posterior[[2]]
    Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
    Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
    
    res = list(est_A = eigen_G_hat$vectors[,1:d],
               est_beta = beta_hat,
               est_sigma0_2 = sigma_2_0_hat,
               est_sigma2 = sigma2_hat,
               mean_obs = Y_hat,
               mean_obs_95lb = Y_95_lower,
               mean_obs_95ub = Y_95_upper)
  }
  else{
    param_ini=c(rep(log(.1),d),rep(log(10),d))
    m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,output=output,delta_x=delta_x,
                lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
    
    while(is.character(m)){
      param_ini=param_ini+runif(2*d)
      #print(is.character(m))
      m=try(optim(param_ini,neg_log_lik_diff_cov_FFBS,A_ini=A_ini,output=output,delta_x=delta_x,
                  lower=c(rep(log(10^{-10}),d),rep(log(10^{-10}),d)),
                  upper=c(rep(log(10^10),d),rep(log(10^10),d)),control = list(maxit = 100), method="L-BFGS-B"),silent=T)
    }
    
    beta_hat=exp(m$par[1:d])
    tau_hat=exp(m$par[(d+1):(2*d)])
    G_log_det_cov_hat=Get_G_log_det_cov(m$par, output, delta_x,d,kernel_type = kernel_type)
    for(i in 1:d){
      G_hat[[i]]=output_2-G_log_det_cov_hat[[i]]
    }
    A_hat=Optimization_Stiefel_Manifold(A_ini, G=G_hat,max_iter=100)
    sigma_2_0_hat=(trace_output_2+F_Funct(A_hat,G_hat))/(n*k)
    sigma2_hat = tau_hat*(sigma_2_0_hat)
    
    est_param=c(beta_hat, sigma2_hat, sigma_2_0_hat)
    pred_GPPCA=pred_FFBS_FastGaSP_diff_cov(param=est_param,A_hat=A_hat,input=input,
                                  testing_input=input, output_here=output,
                                  d=d,var_data=F,kernel_type)
    
    Y_hat=pred_GPPCA[[1]]
    Y_var=pred_GPPCA[[2]]
    Y_95_lower=Y_hat+sqrt(Y_var)*qnorm(0.025)
    Y_95_upper=Y_hat+sqrt(Y_var)*qnorm(0.975)
    
    res = list(est_A = A_hat,
               est_beta = beta_hat,
               est_sigma0_2 = sigma_2_0_hat,
               est_sigma2 = sigma2_hat,
               mean_obs = Y_hat,
               mean_obs_95lb = Y_95_lower,
               mean_obs_95ub = Y_95_upper)
  }
  
  return(res)
}

predict_gppca <- function(param, A_hat, input, step=1, output,d, kernel_type, shared_param, interval=FALSE, interval_data=TRUE){
  n = ncol(output)
  testing_input = n + step
  if(shared_param){
    pred_gppca = pred_FFBS_FastGaSP(param,A_hat,input,testing_input=testing_input, output_here=output,d=d,kernel_type=kernel_type,var_data=interval_data)
  }
  else{
    pred_GPPCA=pred_FFBS_FastGaSP_diff_cov(param,A_hat,input,testing_input=testing_input, output_here=output,d,var_data=interval_data,kernel_type=kernel_type)
  }
  res = list(pred_mean = pred_gppca[[1]])
  Y_hat = pred_gppca[[1]]
  Y_var=pred_gppca[[2]]
  if(interval){
    res$pred_interval_95lb = Y_hat+sqrt(Y_var)*qnorm(0.025)
    res$pred_interval_95ub = Y_hat+sqrt(Y_var)*qnorm(0.975)
  }
  return(res)
}
