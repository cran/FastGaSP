IKF = function(beta, tilde_nu, delta_x, output, kernel_type='matern_5_2'){
  if(kernel_type=='matern_5_2'){
    lambda=sqrt(5)*beta
    
    W0=Construct_W0_matern_5_2(1,lambda)  
    GG=Construct_G_matern_5_2(delta_x,lambda) 
    W=Construct_W_matern_5_2(1.0,delta_x,lambda,W0)
  }else if(kernel_type=='exp'){
    lambda=beta
    W0=Construct_W0_exp(1,lambda)
    GG=Construct_G_exp(delta_x,lambda) 
    W=Construct_W_exp(1,delta_x,lambda,W0)
  }
  Q_K=Get_Q_K(GG,W,W0,tilde_nu)
  res=Get_R_y(GG,pmax(Q_K[[1]],0),Q_K[[2]],output)#-tilde_nu*output
  
  return(res)
}


initialization_Vicsek = function(n_t,v_abs){ ##2D initial position, uniform 
  p0=sqrt(n_t)*runif(n_t*2) ##uniform initial positions, enlarge the box size proportional to the number of particles 
  theta0=2*pi*runif(n_t)#-pi ###uniform direction 
  v0=rep(NA,2*n_t)
  for(i in 1:n_t){
    v0[(i-1)*2+1:2]=c(v_abs*cos(theta0[i]),v_abs*sin(theta0[i]))
  }
  return(list(p0 = p0, v0 = v0, theta0=theta0))
}


#
Vicsek = function(p0,v0,theta0,v_abs,n_t,T_sim,h,cut_r,sigma_0,noise_type='Gaussian'){ ##position, velocity, particle number, time, time interval, sigma is the noise
  pos=matrix(NA,2*n_t,T_sim+1)
  v=matrix(NA,2*n_t,T_sim+1)
  theta=matrix(NA,n_t,T_sim+1)
  
  v[,1]=v0
  pos[,1]=p0
  theta[,1]=theta0
  
  for(t in 1:T_sim){
    input_here=matrix(pos[,t],2,n_t)
    v_vec_here=matrix(v[,t],2,n_t)
    v_here=sqrt(colSums(v_vec_here^2))
    
    for(i in 1:n_t){
      input_i=as.vector(pos[(i-1)*2+1:2,t])
      d_vec_here_all=(input_i-input_here)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      
      index_neighbor=which(d_here<cut_r&v_here>0)
      theta_mean = atan2(mean(v_vec_here[2,index_neighbor]),mean(v_vec_here[1,index_neighbor]))
      
      if(noise_type=='Gaussian'){
        theta_i=theta_mean+sigma_0*rnorm(1)  ##only learn this function
      }else if(noise_type=='Uniform'){
        theta_i=theta_mean+sigma_0*(runif(1)-0.5)  ##only learn this function
      }
      
      v[(i-1)*2+1:2,t+1]=c(v_abs*cos(theta_i),v_abs*sin(theta_i))
      pos[(i-1)*2+1:2,t+1]=pos[(i-1)*2+1:2,t]+v[(i-1)*2+1:2,t]*h
      theta[i,t+1]=theta_i
      
    }
  }
  
  
  return(list(pos = pos, v=v, theta=theta))
}


simulate_Vicsek = function(v_abs, n_t, T_sim, h, cut_r, sigma_0, noise_type = 'Gaussian') {
  # Initialize the system
  initial_all = initialization_Vicsek(n_t = n_t, v_abs = v_abs)
  p0 = initial_all$p0 # Initial location
  v0 = initial_all$v0 # Initial speed
  theta0 = initial_all$theta0 # Initial velocity angle
  
  # Adjust theta0 to be within [-pi, pi]
  theta0[which(theta0 > pi)] = theta0[which(theta0 > pi)] - 2 * pi
  
  # Simulate the trajectory
  m_Vicsek = Vicsek(
    p0 = p0, v0 = v0, theta0 = theta0, 
    v_abs = v_abs, n_t = n_t, 
    T_sim = T_sim, h = h, cut_r = cut_r, 
    sigma_0 = sigma_0, noise_type = noise_type
  )
  
  
  
  # Split position matrix into x and y lists
  px_list = split(m_Vicsek$pos[seq(1, nrow(m_Vicsek$pos), 2), ], 
                  col(m_Vicsek$pos[seq(1, nrow(m_Vicsek$pos), 2), ]))
  py_list = split(m_Vicsek$pos[seq(2, nrow(m_Vicsek$pos), 2), ], 
                  col(m_Vicsek$pos[seq(2, nrow(m_Vicsek$pos), 2), ]))
  
  # Split velocity matrix into x and y lists
  vx_list = split(m_Vicsek$v[seq(1, nrow(m_Vicsek$v), 2), ], 
                  col(m_Vicsek$v[seq(1, nrow(m_Vicsek$v), 2), ]))
  vy_list = split(m_Vicsek$v[seq(2, nrow(m_Vicsek$v), 2), ], 
                  col(m_Vicsek$v[seq(2, nrow(m_Vicsek$v), 2), ]))
  
  # Split theta matrix into list
  theta_list = split(m_Vicsek$theta, col(m_Vicsek$theta))
  
  new("particle.data",
      px_list = px_list,
      py_list = py_list,
      vx_list = vx_list,
      vy_list = vy_list,
      theta_list = theta_list,
      data_type = "simulation",
      n_particles = n_t,
      T_time = T_sim,
      sigma_0 = sigma_0,
      radius = cut_r,
      model = "Vicsek",
      D_y = 1)
  
  
}


get_boundary_grid=function(px_min,px_max,py_min,py_max,nx,ny){
  grid_boundary_mat=matrix(NA,nx*ny,4)
  colnames(grid_boundary_mat)=c('pos_x_min','pos_x_max','pos_y_min','pos_y_max')
  
  ##get slightly larger grid 
  len_x_ori=px_max-px_min
  len_y_ori=py_max-py_min
  delta_x=len_x_ori/nx*0.1
  delta_y=len_y_ori/ny*0.1
  
  px_seq=seq(px_min-delta_x,px_max+delta_x,(px_max-px_min+2*delta_x)/(nx))
  py_seq=seq(py_min-delta_y,py_max+delta_y,(py_max-py_min+2*delta_y)/(ny))
  grid_boundary_mat[,1]=rep(px_seq[1:nx],ny)
  grid_boundary_mat[,2]=rep(px_seq[2:(nx+1)],ny)
  grid_boundary_mat[,3]=as.numeric(t(matrix(py_seq[1:ny],ny,nx)))
  grid_boundary_mat[,4]=as.numeric(t(matrix(py_seq[2:(ny+1)],ny,nx)))
  
  my_grid=list()
  
  
  Lx_min = min(grid_boundary_mat[,1:2]);
  Lx_max = max(grid_boundary_mat[,1:2]);
  Ly_min = min(grid_boundary_mat[,3:4]);
  Ly_max = max(grid_boundary_mat[,3:4]);
  
  len_x=  (max(grid_boundary_mat[,1:2])-min(grid_boundary_mat[,1:2]))/nx
  len_y=  (max(grid_boundary_mat[,3:4])-min(grid_boundary_mat[,3:4]))/ny
  
  
  grid_boundary_info=list()
  grid_boundary_info$grid_boundary_mat=grid_boundary_mat
  grid_boundary_info$grid_info=as.matrix(c(Lx_min,Lx_max,Ly_min,Ly_max,nx,ny,len_x,len_y),8,1)
  rownames(  grid_boundary_info$grid_info)=c('Lx_min','Lx_max','Ly_min','Ly_max','nx','ny','len_x','len_y')
  
  return(grid_boundary_info)
  
}



initiate_grid=function(grid_boundary_info){
  #include_theta = !is.null(theta) #if provided theta, then update theta to grid
  
  # Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  # Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  # Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  # Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  # len_x=unname(grid_boundary_info$grid_info['len_x',])
  # len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  ##form neighboring particle
  ##from x first
  neighbor_index_list=as.list(1:(nx*ny))
  ##exterior
  for(i in 1:(nx*ny)){
    i_x=(i%%nx)
    if(i_x==0){
      i_x=nx
    }
    i_y=ceiling(i/nx)
    
    
    if((i_x-1)>0&(i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x-1))
    }
    
    if((i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x))
    }
    
    if((i_x+1)<=nx&(i_y-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-2)*nx+(i_x+1))
    }
    
    if((i_x-1)>0){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x-1))
    }
    
    if((i_x+1)<=nx){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y-1)*nx+(i_x+1))
    }
    
    if((i_x-1)>0&(i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x-1))
    }
    
    if((i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x))
    }
    if((i_x+1)<=nx&(i_y+1)<=ny){
      neighbor_index_list[[i]]=c( neighbor_index_list[[i]],(i_y)*nx+(i_x+1))
    }
  }
  
  
  
  
  #return(list(m_grid = m_grid, neighbor_index_list = neighbor_index_list))
  return(neighbor_index_list)
}

create_particle_grid=function(grid_boundary_info, pos_x,pos_y,vel_x,vel_y,neighbor_index_list){ #,return_theta=FALSE
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  m_grid=vector(mode = 'list', nx*ny) #as.list(rep(NA,nx*ny))
  n_t=length(pos_x)
  
  for(i in 1:(nx*ny)){
    m_grid[[i]]=list(particle_pos = NULL, neighbor_pos=NULL, particle_vel=NULL, neighbor_vel=NULL)
    # if(return_theta){
    #   m_grid[[i]]=list(particle_pos = NULL, neighbor_pos=NULL, particle_vel=NULL, neighbor_vel=NULL,
    #                    particle_theta=NULL, neighbor_theta=NULL)
    # }else{
    #   
    # }
  }
  
  
  for(i in 1:n_t){
    i_x=ceiling((pos_x[i]-Lx_min)/len_x)
    i_y=ceiling((pos_y[i]-Ly_min)/len_y)
    
    index_grid=(i_y-1)*nx+i_x
    
    ##update pos and vel
    m_grid[[index_grid]]$particle_pos=cbind((m_grid[[index_grid]]$particle_pos),c(pos_x[i],pos_y[i]))
    m_grid[[index_grid]]$particle_vel=cbind((m_grid[[index_grid]]$particle_vel),c(vel_x[i],vel_y[i]))
    
    #if(return_theta) m_grid[[index_grid]]$particle_theta = c(m_grid[[index_grid]]$particle_theta,atan2(vel_y[i],vel_x[i]))
    
  }
  
  for(i in 1:(nx*ny)){
    #print(i)
    if(!is.null(m_grid[[i]]$particle_pos)){# only care about the grid that has particles
      neighbor_index = neighbor_index_list[[i]]
      for(idx in neighbor_index){
        if(!is.null(m_grid[[idx]]$particle_pos)){
          m_grid[[i]]$neighbor_pos=cbind(m_grid[[i]]$neighbor_pos,m_grid[[idx]]$particle_pos)
          m_grid[[i]]$neighbor_vel=cbind(m_grid[[i]]$neighbor_vel,m_grid[[idx]]$particle_vel)
          
          #if(return_theta) m_grid[[i]]$neighbor_theta = c(m_grid[[i]]$neighbor_theta, m_grid[[idx]]$particle_theta)
        }
      }
    }
    
  }
  
  return(m_grid)
  
}

find_grid_neighbors=function(pos_x_list,pos_y_list, vel_x_list,vel_y_list, time_range, #n_t,T_time,
                             grid_boundary_info){ #
  
  
  # Lx_min=grid_boundary_info$grid_info[1]
  # Lx_max=grid_boundary_info$grid_info[2]
  # Ly_min=grid_boundary_info$grid_info[3]
  # Ly_max=grid_boundary_info$grid_info[4]
  # nx=grid_boundary_info$grid_info[5]
  # ny=grid_boundary_info$grid_info[6]
  # len_x=grid_boundary_info$grid_info[7]
  # len_y=grid_boundary_info$grid_info[8]
  
  if(min(time_range)<1 | max(time_range)>length(pos_x_list)){
    stop('invalid time_range')
  }
  
  
  neighbor_index_list=initiate_grid(grid_boundary_info=grid_boundary_info)
  
  neighbors_info = vector(length(time_range), mode = "list")
  names(neighbors_info) <- paste0("time", time_range)
  for(t in time_range){
    #print(t)
    
    m_grid_here=create_particle_grid(grid_boundary_info=grid_boundary_info,
                                     pos_x=pos_x_list[[t]],pos_y=pos_y_list[[t]],
                                     vel_x=vel_x_list[[t]],vel_y=vel_y_list[[t]],
                                     neighbor_index_list=neighbor_index_list)
    
    neighbors_info[[t]] = m_grid_here
    
    
  }
  
  
  return(neighbors_info)
}


#### unnormalized Vicsek #####
unnormalized_Vicsek <- function(p0,v0,n_t,T_sim,h,cut_r,sigma_0,noise_type='Gaussian'){ ##position, velocity, particle number, time, time interval, sigma is the noise
  pos=matrix(NA,2*n_t,T_sim+1)
  v=matrix(NA,2*n_t,T_sim+1)
  
  v[,1]=v0
  pos[,1]=p0
  
  for(t in 1:T_sim){
    input_here=matrix(pos[,t],2,n_t)
    v_vec_here=matrix(v[,t],2,n_t)
    #v_here=sqrt(colSums(v_vec_here^2))
    
    for(i in 1:n_t){
      input_i=as.vector(pos[(i-1)*2+1:2,t])
      d_vec_here_all=(input_i-input_here) 
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      
      index_neighbor=which(d_here<cut_r) #&v_here>0
      vy_neighbor=v_vec_here[2,index_neighbor]
      vx_neighbor=v_vec_here[1,index_neighbor]
      #theta_neighbor=atan2(v_vec_here[2,index_neighbor],v_vec_here[1,index_neighbor])
      
      v[(i-1)*2+1,t+1]=mean(vx_neighbor)
      v[(i-1)*2+2,t+1]=mean(vy_neighbor)
      
      ##add noise
      if(noise_type=='Gaussian'){
        v[(i-1)*2+1,t+1]=  v[(i-1)*2+1,t+1]+sigma_0*rnorm(1)  ##
        v[(i-1)*2+2,t+1]= v[(i-1)*2+2,t+1]+sigma_0*rnorm(1)  ##
      }else if(noise_type=='Uniform'){
        v[(i-1)*2+1,t+1]=  v[(i-1)*2+1,t+1]+sigma_0*(runif(1)-0.5)  ##
        v[(i-1)*2+2,t+1]= v[(i-1)*2+2,t+1]+sigma_0*(runif(1)-0.5)
      }
      
      pos[(i-1)*2+1:2,t+1]=pos[(i-1)*2+1:2,t]+v[(i-1)*2+1:2,t]*h
      
    }
  }
  
  
  return(list(pos = pos, v=v))
}

simulate_unnormalized_Vicsek = function(v_abs, n_t, T_sim, h, cut_r, sigma_0, noise_type = 'Gaussian') {
  # Initialize the system
  initial_all = initialization_Vicsek(n_t = n_t, v_abs = v_abs)
  p0 = initial_all$p0 # Initial location
  v0 = initial_all$v0 # Initial speed
  
  
  # Simulate the trajectory
  m_unnormalized_Vicsek = unnormalized_Vicsek(p0=p0,v0=v0,n_t=n_t,T_sim=T_sim,
                                              h=h,cut_r=cut_r,sigma_0=sigma_0,noise_type=noise_type)
  
  # input_pos_all = m_Vicsek$pos
  # v_all = m_Vicsek$v
  # theta_all = m_Vicsek$theta
  
  # Split position matrix into x and y lists
  px_list = split(m_unnormalized_Vicsek$pos[seq(1, nrow(m_unnormalized_Vicsek$pos), 2), ], 
                  col(m_unnormalized_Vicsek$pos[seq(1, nrow(m_unnormalized_Vicsek$pos), 2), ]))
  py_list = split(m_unnormalized_Vicsek$pos[seq(2, nrow(m_unnormalized_Vicsek$pos), 2), ], 
                  col(m_unnormalized_Vicsek$pos[seq(2, nrow(m_unnormalized_Vicsek$pos), 2), ]))
  
  # Split velocity matrix into x and y lists
  vx_list = split(m_unnormalized_Vicsek$v[seq(1, nrow(m_unnormalized_Vicsek$v), 2), ], 
                  col(m_unnormalized_Vicsek$v[seq(1, nrow(m_unnormalized_Vicsek$v), 2), ]))
  vy_list = split(m_unnormalized_Vicsek$v[seq(2, nrow(m_unnormalized_Vicsek$v), 2), ], 
                  col(m_unnormalized_Vicsek$v[seq(2, nrow(m_unnormalized_Vicsek$v), 2), ]))
  
  new("particle.data",
      px_list = px_list,
      py_list = py_list,
      vx_list = vx_list,
      vy_list = vy_list,
      data_type = "simulation",
      n_particles = n_t,
      T_time = T_sim,
      sigma_0 = sigma_0,
      radius = cut_r,
      model = "unnormalized_Vicsek",
      D_y = 2)
}



form_neighbors_unnormalized_Vicsek_with_r=function(threshold_r,pos_x_list,pos_y_list,vel_x_list,vel_y_list,
                                                   time_range,grid_boundary_info,neighbors_info,D_y){
  
  T_total = length(time_range)
  n_t = length(pos_x_list[[1]])
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  A_vec=rep(NA,n_t*T_total*15*D_y)
  v_vec=rep(NA,n_t*T_total*15*D_y)
  num_neighbors_vec=rep(NA,n_t*T_total) ##here N neighbor each particle (obs) has
  
  count=0
  
  for(i_t in 1:length(time_range)){
    t = time_range[i_t]
    pos_x_t=pos_x_list[[t]]
    pos_y_t=pos_y_list[[t]]
    vel_x_t=vel_x_list[[t]]
    vel_y_t=vel_y_list[[t]]
    m_grid_here = neighbors_info[[t]]
    
    for(i in 1:n_t){
      input_pos_i=as.vector(c(pos_x_t[i],pos_y_t[i]))
      
      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_pos_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      
      index_neighbor=which(d_here<threshold_r)
      # if(apolar_vicsek==F){
      #   index_neighbor=which(d_here<threshold_r)
      # }else{
      #   index_neighbor=which(d_here<threshold_r)
      #   index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>=0)
      #   index_neighbor=intersect(index_neighbor,index_same_v_direction)
      # }
      n_neighbor=length(index_neighbor)
      num_neighbors_vec[n_t*(i_t-1)+i]=n_neighbor
      
      v_vec[count*D_y+(1:(n_neighbor*D_y))] = c(m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor],
                                                m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor])
      A_vec[2*count*D_y+(1:(2*n_neighbor*D_y))]=c(rep(c(1/n_neighbor,0),n_neighbor),rep(c(0,1/n_neighbor),n_neighbor))
      
      
      
      
      count=count+n_neighbor
      
      
    }
  }
  
  
  A_vec=A_vec[1:(2*count*D_y)]
  v_vec=v_vec[1:(count*D_y)]
  
  return(list(A_vec=A_vec, v_vec=v_vec, num_neighbors_vec=num_neighbors_vec))
}



pred_ho_output_unnormalized_Vicsek_log_RMSE=function(param, kernel_type, neighbors_info,grid_boundary_info,
                                                     pos_x_list, pos_y_list, vel_x_list, vel_y_list,
                                                     T_index_train, T_index_ho, output, ho_output, 
                                                     D_y, cut_r_max, tilde_nu, tol, maxIte){
  beta=exp(param[1])
  tau=exp(param[2])  ###sigma_2/sigma_2_0
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  T_train=length(T_index_train) ##length of training data
  T_ho=length(T_index_ho) ##hold out prediction
  n_t = length(pos_x_list[[1]])
  
  
  ans_neighbors_train=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                                                vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_train,
                                                                grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  ans_neighbors_ho=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                                             vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_ho,
                                                             grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  N_train=n_t*T_train*D_y 
  N_ho=n_t*T_ho*D_y
  
  ##sort d_train
  sort_d_train=sort(ans_neighbors_train$v_vec,index.return=T) ###sorted d
  N_tilde=length(ans_neighbors_train$v_vec) ###N_tilde, distances  
  delta_x_train=sort_d_train$x[-1]-sort_d_train$x[-N_tilde]
  sort_d_train_ix=sort_d_train$ix
  
  ##form augmented samples for cross-validation
  N_ho_tilde=length(ans_neighbors_ho$v_vec)  #number of hold out distance
  d_aug=c(ans_neighbors_ho$v_vec,(ans_neighbors_train$v_vec))
  d_aug_sort=sort(d_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
  d_aug_sort_x=d_aug_sort$x
  d_aug_sort_rev_ix=sort(d_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort 
  
  delta_x_aug=d_aug_sort_x[2:length(d_aug_sort_x)]-d_aug_sort_x[1:(length(d_aug_sort_x)-1)]
  
  ###finish construction, now start to predict
  
  m_CG=IKF_CG_particle( param,  kernel_type,   delta_x_train,   output, ans_neighbors_train$A_vec, # [[1]], 
                        sort_d_train_ix,  ans_neighbors_train$num_neighbors_vec*2, tilde_nu,
                        D_y,   N_train,   tol=tol,  maxIte = maxIte)
  
  ans_CG=m_CG[[1]] 
  
  ##this change this back to original parameterization
  ans_CG_tilde=ans_CG*tau ##this gives R_inv_y 
  
  ###z=A_t_sparse_times_x, 
  w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=ans_neighbors_train$A_vec,  num_neighbors_vec=2*ans_neighbors_train$num_neighbors_vec,  
                            D_y=D_y, N_tilde=N_tilde)
  w_aug=c(rep(0,N_ho_tilde),w_CG)
  
  # param_here=log(c(beta,tilde_nu)) ##tilde nu is one to stablize the computation
  # pred_mean_aug=R_times_z(param_here, have_noise=T, delta_x=delta_x_aug, z=w_aug[d_aug_sort$ix],
  #                         kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
  pred_mean_aug = IKF(beta=beta, tilde_nu=tilde_nu, 
                      delta_x=delta_x_aug, output=w_aug[d_aug_sort$ix], kernel_type=kernel_type) - tilde_nu*w_aug[d_aug_sort$ix]
  
  pred_mean_fast=pred_mean_aug[d_aug_sort_rev_ix][1:N_ho_tilde]
  
  ##can only observes output so cross-validation on output
  pred_mean_ho_output=A_times_x_particle(output= pred_mean_fast, A_all_v= ans_neighbors_ho$A_vec,  num_neighbors_vec = 2*ans_neighbors_ho$num_neighbors_vec,  
                                         D_y,  N_ho)
  
  log_RMSE_ho=1/2*log(mean( (ho_output-pred_mean_ho_output)^2)) ##many pars should work as it contains noises
  
  #print(c(beta,tau,threshold_r,log_RMSE_ho))
  
  return(log_RMSE_ho)
  
  
}



particle_interaction_est_unnormalized_Vicsek = function(data_obj, param, cut_r_max, est_param=TRUE, nx=NULL, ny=NULL,
                                                        kernel_type='matern_5_2', tilde_nu=0.1, tol=1e-6, maxIte=1000, 
                                                        output=NULL, ho_output=NULL, testing_input=NULL, compute_CI = TRUE){ 
  
  px_list = data_obj@px_list
  py_list = data_obj@py_list
  vx_list = data_obj@vx_list
  vy_list = data_obj@vy_list
  n_t = data_obj@n_particles
  T_sim = data_obj@T_time
  D_y = data_obj@D_y
  
  N=n_t*T_sim*D_y 
  
  
  T_index_time = 1:T_sim
  
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  #print(grid_boundary_info)
  
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  
  
  
  if(est_param){
    ## split train and hold-out validation
    T_index_ho=seq(5,T_sim,5) ##every 5 use the last one as holdout
    T_index_train=(1:T_sim)[-T_index_ho]
    
    T_train=length(T_index_train) ##length of training data
    T_ho=length(T_index_ho) ##hold out prediction
    
    if(is.null(output)){
      output=as.vector(rbind(
        unlist(vx_list[T_index_train + 1]),
        unlist(vy_list[T_index_train + 1])
      ))#as.vector(v_all[,1+T_index_train])  
    } 
    if(is.null(ho_output)){
      ho_output=as.vector(rbind(
        unlist(vx_list[T_index_ho + 1]),
        unlist(vy_list[T_index_ho + 1])
      ))# as.vector(v_all[,1+T_index_ho])  
    }
    
    
    m_IKF=optim(param,pred_ho_output_unnormalized_Vicsek_log_RMSE, control=list(maxit=200),
                #lower=c(-8,-8,-8), upper=c(5,1,2),
                kernel_type=kernel_type, neighbors_info=neighbors_info,grid_boundary_info=grid_boundary_info,
                pos_x_list=px_list, pos_y_list=py_list, vel_x_list=vx_list, vel_y_list=vy_list,
                T_index_train=T_index_train, T_index_ho=T_index_ho, output=output, ho_output=ho_output, 
                D_y=D_y, cut_r_max=cut_r_max, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte)
    
    while(m_IKF$par[2]>log(10^6)){
      param[2]=param[2]-log(10)+runif(1)
      m_IKF=optim(param,pred_ho_output_unnormalized_Vicsek_log_RMSE, control=list(maxit=200),
                  #lower=c(-8,-8,-8), upper=c(5,1,2),
                  kernel_type=kernel_type, neighbors_info=neighbors_info,grid_boundary_info=grid_boundary_info,
                  pos_x_list=px_list, pos_y_list=py_list, vel_x_list=vx_list, vel_y_list=vy_list,
                  T_index_train=T_index_train, T_index_ho=T_index_ho, output=output, ho_output=ho_output, 
                  D_y=D_y, cut_r_max=cut_r_max, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte)
      
    }
    
    param=m_IKF$par
  }
  
  beta=exp(param[1])
  tau=exp(param[2])  ###sigma_2/sigma_2_0
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  
  parameters = c(beta, tau, threshold_r)
  names(parameters) = c('beta', 'tau', 'radius')
  
  
  # # prediction
  output_all=as.vector(rbind(
    unlist(vx_list[1+T_index_time]),
    unlist(vy_list[1+T_index_time])
  ))
  
  ans_neighbors_all=form_neighbors_unnormalized_Vicsek_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                                              vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                                              grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  A_all_vec=ans_neighbors_all$A_vec
  v_all_vec=ans_neighbors_all$v_vec
  num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
  sort_v_all=sort(v_all_vec,index.return=T)
  N_tilde_all=length(v_all_vec) ###this is N_j in the paper
  
  
  delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
  
  
  m_CG=IKF_CG_particle(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                       A_all_v = A_all_vec, sort_d_all_ix=sort_v_all$ix,  num_neighbors_vec=2*num_neighbors_all_vec, tilde_nu=tilde_nu,
                       D=D_y,  N=N,   tol=tol,  maxIte = maxIte)
  
  ans_CG=m_CG[[1]] 
  
  ans_CG_tilde=ans_CG*tau
  
  sigma_2_0_est = output_all%*%ans_CG/length(output_all) ##sometimes negative? solved
  
  ###z=A_t_sparse_times_x, get the weight; maybe write this in C++
  w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=A_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,  
                            D_y=D_y, N_tilde=N_tilde_all)
  
  
  if(!is.null(testing_input)){
    testing_n = length(testing_input)
    
    
    sigma_2_est=sigma_2_0_est*tau
    
    param=log(c(beta,tau))
    
    d_aug=c(testing_input,(v_all_vec))
    d_aug_sort=sort(d_aug,index.return=T)
    d_aug_sort_x=d_aug_sort$x
    d_aug_sort_rev_ix=sort(d_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort 
    
    delta_x_aug=d_aug_sort_x[2:length(d_aug_sort_x)]-d_aug_sort_x[1:(length(d_aug_sort_x)-1)]
    
    
    w_aug=c(rep(0,testing_n),w_CG)
    
    ###this should go back to nu
    # param_tilde=log(c(beta,tilde_nu)) 
    # pred_mean_aug=R_times_z(param_tilde, have_noise=T, delta_x=delta_x_aug, z=w_aug[d_aug_sort$ix],
    #                         kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
    pred_mean_aug = IKF(beta=beta, tilde_nu=tilde_nu, 
                        delta_x=delta_x_aug, output=w_aug[d_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_aug[d_aug_sort$ix]
    # if(kernel_type=='matern_5_2'){
    #   lambda=sqrt(5)*beta
    #   
    #   W0=Construct_W0_matern_5_2(1,lambda)  
    #   GG=Construct_G_matern_5_2(delta_x_aug,lambda) 
    #   W=Construct_W_matern_5_2(1.0,delta_x_aug,lambda,W0)
    # }else if(kernel_type=='exp'){
    #   lambda=beta
    #   W0=Construct_W0_exp(1,lambda)
    #   GG=Construct_G_exp(delta_x_aug,lambda) 
    #   W=Construct_W_exp(1,delta_x_aug,lambda,W0)
    # }
    # Q_K=Get_Q_K(GG,W,W0,tilde_nu)
    # pred_mean_aug=Get_R_y(GG,Q_K[[1]],Q_K[[2]],w_aug[d_aug_sort$ix])-tilde_nu*w_aug[d_aug_sort$ix]
    
    
    pred_mean_fast=pred_mean_aug[d_aug_sort_rev_ix][1:testing_n]
    
    
    #NRMSE = mean( (pred_mean_fast-testing_output)^2)/sd(testing_output)
    if(compute_CI){
      #predictive variance
      c_star=rep(NA,testing_n)
      r0=abs(outer(testing_input,(v_all_vec),'-'))
      if(kernel_type=='exp'){
        r = exp(-beta*r0)
      }else if(kernel_type=='matern_5_2'){
        r = matern_5_2_funct(r0, beta)
      }
      
      #system.time(
      print("Computing the predictive variance ...")
      for(i in 1:testing_n ){
        #print(i)
        
        A_r_i=A_times_x_particle(output=r[i,], A_all_v=A_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,
                                 D=D_y, N)
        
        #tol=sd(A_r_i)^2*0.01*N_tilde ##can make it smaller
        tol_interval=tol*10^{-4}
        R_inv_r_all=IKF_CG_particle( param=param,  kernel_type=kernel_type,   delta_x_all=delta_x_all,   output=A_r_i,
                                     A_all_v=A_all_vec, sort_v_all$ix,  num_neighbors_vec=2*num_neighbors_all_vec, tilde_nu,
                                     D=D_y,   N=N,   tol=tol_interval,  maxIte = maxIte)
        R_inv_r=R_inv_r_all[[1]]*tau
        r_R_inv_r=A_r_i%*%R_inv_r
        
        c_star[i]=1-r_R_inv_r
        
      }
      c_star = abs(c_star)
      ##95 intervals
      LB95=    pred_mean_fast+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.025)
      UB95=    pred_mean_fast+sqrt(as.numeric(sigma_2_est)*c_star)*qnorm(0.975)
    }
    
    # coverage95 = length(which(testing_output<UB95 & testing_output>LB95))/length(testing_output)
    # length95=mean(UB95-LB95)
    
    #est_par_val=c(exp(m_IKF$par), exp(m_IKF$value))
    
    #return(list(NRMSE=NRMSE,coverage95=coverage95,length95=length95,est_par_val=est_par_val,pred_mean_fast=pred_mean_fast,LB95=LB95,UB95=UB95))
    #return(list(est_par_val=est_par_val,pred_mean_fast=pred_mean_fast,LB95=LB95,UB95=UB95))
    
  }
  
  
  
  # return(list(input_pos_all=input_pos_all,v_all = v_all,theta_all = theta_all,
  #             parameters = parameters, 
  #             pred_mean_fast = pred_mean_fast, LB95 = LB95, UB95 = UB95))
  new("particle.est",
      D_y=D_y,
      parameters = parameters,  # This contains the estimated parameters
      sigma_2_0_est = sigma_2_0_est[1,1],
      predictions = if(!is.null(testing_input)) {
        if(compute_CI) {
          list(mean = pred_mean_fast, lower95 = LB95, upper95 = UB95)
        } else {
          list(mean = pred_mean_fast)
        }
      } else {
        NULL
      },
      training_data = list(
        training_velocity = v_all_vec,
        A_v = A_all_vec,
        num_neighbors = num_neighbors_all_vec
        
      ),
      gp_weights = matrix(w_CG)  # Reshape weights if necessary
  )
}



#### Vicsek variation #####

f_Vicsek_variation=function(r,a=0.02,b=1,r_min=0.01,r_max=0.8, beta=20){
  fr=-a/(r+r_min)-b*(r-r_max)+a/r_max
  return(beta*fr)
}
##here D has to be 2
Vicsek_variation = function(p0,v0,n_t,D=2,T_sim,h,cut_r,sigma_0,
                            noise_type='Gaussian'){ ##position, velocity, particle number, time, time interval, sigma is the noise
  pos=matrix(NA,D*n_t,T_sim+1)
  v=matrix(NA,D*n_t,T_sim+1)
  #theta=matrix(NA,n_t,T_sim+1)
  
  v[,1]=v0
  pos[,1]=p0
  #ans_list=as.list(1:2)
  #set.seed(0)
  for(t in 1:T_sim){
    input_here=matrix(pos[,t],D,n_t)
    v_vec_here=matrix(v[,t],D,n_t)
    
    for(i in 1:n_t){
      input_i=input_here[,i]#as.vector(pos[(i-1)*D+1:D,t])
      d_vec_here_all=(input_here-input_i) ###same as in the paper
      
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      v_here=sqrt(colSums(v_vec_here^2))
      ##need to delete itself
      index_neighbor=which(d_here<cut_r&d_here>0) #v_here>0
      
      
      if(length(index_neighbor)>0){
        #n_neighbor=length(index_neighbor)
        vx_sum_neighbor=v_vec_here[1,i]+sum(v_vec_here[1,index_neighbor])
        vy_sum_neighbor=v_vec_here[2,i]+sum(v_vec_here[2,index_neighbor])
        
        f_neighbor=f_Vicsek_variation(d_here[index_neighbor],r_max=cut_r)
        d_norm_x=d_vec_here_all[1,index_neighbor]/d_here[index_neighbor] ##normalized version
        d_norm_y=d_vec_here_all[2,index_neighbor]/d_here[index_neighbor] ##normalized version
        
        fx_sum_neighbor=sum(f_neighbor*(d_norm_x))
        fy_sum_neighbor=sum(f_neighbor*(d_norm_y))
        
        v_norm= (length(index_neighbor))
        
        ##D is 2
        v[(i-1)*D+1,t+1]= (vx_sum_neighbor/(v_norm+1)+fx_sum_neighbor/(v_norm))
        v[(i-1)*D+2,t+1]=(vy_sum_neighbor/(v_norm+1)+fy_sum_neighbor/(v_norm))
        
      }else{
        
        v[(i-1)*D+1,t+1]=v_vec_here[1,i]
        v[(i-1)*D+2,t+1]=v_vec_here[2,i]
        
      }
      
      #v_norm= sqrt(v[(i-1)*D+1,t+1]^2+ v[(i-1)*D+2,t+1]^2)
      #v_norm=1
      
      ##add noise
      if(noise_type=='Gaussian'){
        v[(i-1)*D+1,t+1]=  v[(i-1)*D+1,t+1]+sigma_0*rnorm(1)  ##
        v[(i-1)*D+2,t+1]= v[(i-1)*D+2,t+1]+sigma_0*rnorm(1)  ##
      }else if(noise_type=='Uniform'){
        v[(i-1)*D+1,t+1]=  v[(i-1)*D+1,t+1]+sigma_0*(runif(1)-0.5)  ##
        v[(i-1)*D+2,t+1]= v[(i-1)*D+2,t+1]+sigma_0*(runif(1)-0.5)
      }
      
      pos[(i-1)*D+1:D,t+1]=pos[(i-1)*D+1:D,t]+v[(i-1)*D+1:D,t]*h
      
      #theta[i,t+1]=atan2(v[(i-1)*D+1,t+1],v[(i-1)*D+2,t+1]) ##
      
      #c(v_abs*cos(theta_i),v_abs*sin(theta_i))
      
    }
  }
  ans_list = list(pos=pos, v=v) #, theta=theta
  
  return(ans_list)
  
}



form_neighbors_Vicsek_variation_with_r=function(threshold_r,pos_x_list,pos_y_list,vel_x_list,vel_y_list,
                                                time_range,grid_boundary_info,neighbors_info,D_y){
  T_total = length(time_range)
  n_t = length(pos_x_list[[1]])
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  A_v_vec=rep(NA,n_t*T_total*15*D_y)
  A_f_vec=rep(NA,n_t*T_total*15*D_y)
  d_pos_vec=rep(NA,n_t*T_total*15)
  v_vec=rep(NA,n_t*T_total*15*D_y)
  num_neighbors_vec=rep(NA,n_t*T_total) 
  
  count=0
  
  #count_neighbor_start=0
  
  for(i_t in 1:length(time_range)){
    t = time_range[i_t]
    pos_x_t=pos_x_list[[t]]
    pos_y_t=pos_y_list[[t]]
    vel_x_t=vel_x_list[[t]]
    vel_y_t=vel_y_list[[t]]
    m_grid_here = neighbors_info[[t]]
    
    for(i in 1:n_t){
      input_pos_i=as.vector(c(pos_x_t[i],pos_y_t[i]))
      input_vel_i=as.vector(c(vel_x_t[i],vel_y_t[i]))
      
      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=as.matrix(m_grid_here[[index_grid]]$neighbor_pos)-input_pos_i
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      
      #index_neighbor=which(d_here<cut_r_max)
      index_neighbor=which(d_here<threshold_r&d_here>0) # add i itself (d_here == 0) at the end later
      index_itself = which(d_here==0)
      
      n_neighbor=length(index_neighbor)+1
      num_neighbors_vec[n_t*(i_t-1)+i]=n_neighbor
      
      v_vec[count*D_y+(1:(n_neighbor*D_y))] = c(m_grid_here[[index_grid]]$neighbor_vel[1,c(index_neighbor,index_itself)],
                                                m_grid_here[[index_grid]]$neighbor_vel[2,c(index_neighbor,index_itself)])
      A_v_vec[2*count*D_y+(1:(2*n_neighbor*D_y))]=c(rep(c(1/n_neighbor,0),n_neighbor),rep(c(0,1/n_neighbor),n_neighbor))
      
      if(n_neighbor>1){
        e_vec_here=as.vector(d_vec_here_all[,index_neighbor]/rep(d_here[index_neighbor],each=D_y))
        A_f_vec[count*D_y+(1:(n_neighbor*D_y))]=c(e_vec_here/(n_neighbor-1),0,0)
      }else{
        #e_vec_here=c(0,0)
        A_f_vec[count*D_y+(1:D_y)]=c(0,0)
      }
      
      d_pos_vec[count+(1:n_neighbor)]=c(d_here[index_neighbor],0)
      
      count=count+n_neighbor
      
      
    }
  }
  A_v_vec=A_v_vec[1:(2*count*D_y)]
  A_f_vec=A_f_vec[1:(count*D_y)]
  d_pos_vec=d_pos_vec[1:count]
  v_vec=v_vec[1:(count*D_y)]
  
  ans_list = list(A_v_vec=A_v_vec,A_f_vec=A_f_vec,
                  v_vec=v_vec,d_pos_vec=d_pos_vec,
                  num_neighbors_vec=num_neighbors_vec)
  return(ans_list)
}

pred_ho_output_Vicsek_variation_log_RMSE=function(param, kernel_type, neighbors_info,grid_boundary_info,
                                                  pos_x_list, pos_y_list, vel_x_list, vel_y_list,
                                                  T_index_train, T_index_ho, output, ho_output, 
                                                  D_y, cut_r_max, tilde_nu, tol, maxIte){
  beta_v=exp(param[1])
  beta_f=exp(param[2])
  tau_v=exp(param[3])  
  tau_f=exp(param[4])  
  threshold_r=exp(param[5])/(1+exp(param[5])) * cut_r_max
  
  T_train=length(T_index_train) ##length of training data
  T_ho=length(T_index_ho) ##hold out prediction
  n_t = length(pos_x_list[[1]])
  
  ##form neighbor of train
  # ans_neighbors_train=form_neighbors_Vicsek_variation_with_max_cut_r(T_time=T_train,n_t=n_t,D_y=D_y,threshold_r=threshold_r,
  #                                                                    d_pos_vec_max=max_neighbors_train$d_pos_vec,
  #                                                                    v_vec_max=max_neighbors_train$v_vec,
  #                                                                    e_vec_max=max_neighbors_train$e_vec,
  #                                                                    num_neighbors_vec_max=max_neighbors_train$num_neighbors_vec)
  
  ans_neighbors_train=form_neighbors_Vicsek_variation_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                                             vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_train,
                                                             grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  ##form neighbor of hold-out cross-validation
  # ans_neighbors_ho=form_neighbors_Vicsek_variation_with_max_cut_r(T_time=T_ho,n_t=n_t,D_y=D_y,threshold_r=threshold_r,
  #                                                                 d_pos_vec_max=max_neighbors_ho$d_pos_vec,
  #                                                                 v_vec_max=max_neighbors_ho$v_vec,
  #                                                                 e_vec_max=max_neighbors_ho$e_vec,
  #                                                                 num_neighbors_vec_max=max_neighbors_ho$num_neighbors_vec)
  
  ans_neighbors_ho=form_neighbors_Vicsek_variation_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                                          vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_ho,
                                                          grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  
  
  N_train=n_t*T_train*D_y ##this one is N in training, the output y dimension
  N_ho=n_t*T_ho*D_y ##this one is N in testing
  
  ######### process for vx ###########
  ##sort vx_train
  sort_v_train=sort(ans_neighbors_train$v_vec,index.return=T) ###sorted vx
  N_v_tilde=length(ans_neighbors_train$v_vec) ###N_tilde, distances
  delta_v_train=sort_v_train$x[-1]-sort_v_train$x[-N_v_tilde]
  sort_v_train_ix=sort_v_train$ix
  
  ##form augmented samples for cross-validation
  N_v_ho_tilde=length(ans_neighbors_ho$v_vec)  #number of hold out distance
  v_aug=c(ans_neighbors_ho$v_vec,(ans_neighbors_train$v_vec))
  v_aug_sort=sort(v_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
  v_aug_sort_x=v_aug_sort$x
  v_aug_sort_rev_ix=sort(v_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort
  
  delta_v_aug=v_aug_sort_x[2:length(v_aug_sort_x)]-v_aug_sort_x[1:(length(v_aug_sort_x)-1)]
  
  ######### process for f ###########
  ##sort f_train
  sort_f_train=sort(ans_neighbors_train$d_pos_vec,index.return=T) ###sorted f
  N_f_tilde=length(ans_neighbors_train$d_pos_vec) ###N_tilde, distances
  delta_f_train=sort_f_train$x[-1]-sort_f_train$x[-N_f_tilde]
  sort_f_train_ix=sort_f_train$ix
  
  ##form augmented samples for cross-validation
  N_f_ho_tilde=length(ans_neighbors_ho$d_pos_vec)  #number of hold out distance
  f_aug=c(ans_neighbors_ho$d_pos_vec,(ans_neighbors_train$d_pos_vec))
  f_aug_sort=sort(f_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
  f_aug_sort_x=f_aug_sort$x
  f_aug_sort_rev_ix=sort(f_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort
  
  delta_f_aug=f_aug_sort_x[2:length(f_aug_sort_x)]-f_aug_sort_x[1:(length(f_aug_sort_x)-1)]
  
  ###finish construction, now start to predict
  m_CG=IKF_CG_particle_two_interact(param1=log(c(beta_v,tau_v)), param2=log(c(beta_f,tau_f)), 
                                    kernel_type1=kernel_type, kernel_type2=kernel_type, 
                                    delta_x_all1=delta_v_train, delta_x_all2=delta_f_train, 
                                    A_all_v1=ans_neighbors_train$A_v_vec, A_all_v2=ans_neighbors_train$A_f_vec, 
                                    sort_d_all_ix1=sort_v_train_ix, sort_d_all_ix2=sort_f_train_ix, 
                                    num_neighbors_vec1=2*ans_neighbors_train$num_neighbors_vec, num_neighbors_vec2=ans_neighbors_train$num_neighbors_vec, 
                                    output=output, tilde_nu=tilde_nu, 
                                    D=D_y, N=N_train, tol = tol, maxIte = maxIte)
  
  
  ans_CG=m_CG[[1]]
  
  ##this change this back to original parameterization
  ans_CG_v_tilde=ans_CG*tau_v ##this gives R_inv_y
  
  ###z=A_t_sparse_times_x,
  w_CG_v=A_t_times_x_particle(output=ans_CG_v_tilde, A_all_v=ans_neighbors_train$A_v_vec,  num_neighbors_vec=2*ans_neighbors_train$num_neighbors_vec,
                              D_y=D_y, N_tilde=N_v_tilde)
  w_v_aug=c(rep(0,N_v_ho_tilde),w_CG_v)
  
  # param_here=log(c(beta_v,tilde_nu)) ##tilde nu is one to stablize the computation
  # pred_mean_v_aug=R_times_z(param_here, have_noise=T, delta_x=delta_v_aug, z=w_v_aug[v_aug_sort$ix],
  #                           kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
  
  pred_mean_v_aug = IKF(beta=beta_v, tilde_nu=tilde_nu, 
                        delta_x=delta_v_aug, output=w_v_aug[v_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
  
  pred_mean_v_fast=pred_mean_v_aug[v_aug_sort_rev_ix][1:N_v_ho_tilde]
  
  ##can only observes output so cross-validation on output
  pred_mean_v_ho_output=A_times_x_particle( pred_mean_v_fast,  ans_neighbors_ho$A_v_vec,  2*ans_neighbors_ho$num_neighbors_vec,
                                            D_y,  N_ho)
  
  ##this change this back to original parameterization
  ans_CG_f_tilde=ans_CG*tau_f ##this gives R_inv_y
  
  ###z=A_t_sparse_times_x,
  w_CG_f=A_t_times_x_particle(output=ans_CG_f_tilde, A_all_v=ans_neighbors_train$A_f_vec,  num_neighbors_vec=ans_neighbors_train$num_neighbors_vec,
                              D_y=D_y, N_tilde=N_f_tilde)
  w_f_aug=c(rep(0,N_f_ho_tilde),w_CG_f)
  
  # param_here=log(c(beta_f,tilde_nu)) ##tilde nu is one to stablize the computation
  # pred_mean_f_aug=R_times_z(param_here, have_noise=T, delta_x=delta_f_aug, z=w_f_aug[f_aug_sort$ix],
  #                           kernel_type=kernel_type)-tilde_nu*w_f_aug[f_aug_sort$ix]
  pred_mean_f_aug = IKF(beta=beta_f, tilde_nu=tilde_nu, 
                        delta_x=delta_f_aug, output=w_f_aug[f_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_f_aug[f_aug_sort$ix]
  
  pred_mean_f_fast=pred_mean_f_aug[f_aug_sort_rev_ix][1:N_f_ho_tilde]
  
  ##can only observes output so cross-validation on output
  pred_mean_f_ho_output=A_times_x_particle( pred_mean_f_fast,  ans_neighbors_ho$A_f_vec,  ans_neighbors_ho$num_neighbors_vec,
                                            D_y,  N_ho)
  
  
  pred_mean_ho_output=pred_mean_v_ho_output+pred_mean_f_ho_output
  
  log_RMSE_ho=1/2*log(mean( (ho_output-pred_mean_ho_output)^2)) ##many pars should work as it contains noises
  
  #print(c(beta_v,beta_f,tau_v,tau_f,threshold_r,log_RMSE_ho))
  
  return(log_RMSE_ho)
  
  
}


simulate_Vicsek_variation = function(v_abs, n_t, T_sim, h, cut_r, sigma_0, noise_type = 'Gaussian'){
  initial_all=initialization_Vicsek(n_t=n_t,v_abs=v_abs)
  p0=initial_all$p0 # initial location
  v0=initial_all$v0 # initial speed
  
  ##simulate the trajectory
  m_Vicsek_variation=Vicsek_variation(
    p0=p0,v0=v0,n_t=n_t, D = 2,
    T_sim=T_sim,h=h,cut_r=cut_r,
    sigma_0=sigma_0,noise_type=noise_type
  )
  # input_pos_all = m_Vicsek_variation$pos 
  # v_all = m_Vicsek_variation$v
  
  # Split position matrix into x and y lists
  px_list = split(m_Vicsek_variation$pos[seq(1, nrow(m_Vicsek_variation$pos), 2), ], 
                  col(m_Vicsek_variation$pos[seq(1, nrow(m_Vicsek_variation$pos), 2), ]))
  py_list = split(m_Vicsek_variation$pos[seq(2, nrow(m_Vicsek_variation$pos), 2), ], 
                  col(m_Vicsek_variation$pos[seq(2, nrow(m_Vicsek_variation$pos), 2), ]))
  
  # Split velocity matrix into x and y lists
  vx_list = split(m_Vicsek_variation$v[seq(1, nrow(m_Vicsek_variation$v), 2), ], 
                  col(m_Vicsek_variation$v[seq(1, nrow(m_Vicsek_variation$v), 2), ]))
  vy_list = split(m_Vicsek_variation$v[seq(2, nrow(m_Vicsek_variation$v), 2), ], 
                  col(m_Vicsek_variation$v[seq(2, nrow(m_Vicsek_variation$v), 2), ]))
  
  new("particle.data",
      px_list = px_list,
      py_list = py_list,
      vx_list = vx_list,
      vy_list = vy_list,
      data_type = "simulation",
      n_particles = n_t,
      T_time = T_sim,
      sigma_0 = sigma_0,
      radius = cut_r,
      model = "two_interactions_Vicsek",
      D_y = 2)
}

particle_interaction_est_Vicsek_variation=function(data_obj, param, cut_r_max, est_param=TRUE, nx=NULL, ny=NULL,#px_list, py_list, vx_list, vy_list, n_t, T_sim, D_y=2, 
                                                   kernel_type='matern_5_2',tilde_nu=0.1, tol=1e-6, maxIte=1000,
                                                   output=NULL, ho_output=NULL, testing_v_input=NULL, testing_d_input=NULL, compute_CI = TRUE){
  
  px_list = data_obj@px_list
  py_list = data_obj@py_list
  vx_list = data_obj@vx_list
  vy_list = data_obj@vy_list
  n_t = data_obj@n_particles
  T_sim = data_obj@T_time
  D_y = data_obj@D_y
  
  N=n_t*T_sim*D_y 
  
  T_index_time = 1:T_sim
  
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  
  
  if(est_param){
    T_index_ho=seq(5,T_sim,5) ##every 5 use the last one as holdout
    T_index_train=(1:T_sim)[-T_index_ho]
    
    T_train=length(T_index_train) ##length of training data
    T_ho=length(T_index_ho) ##hold out prediction
    
    if(is.null(output)){
      output=as.vector(rbind(
        unlist(vx_list[T_index_train + 1]),
        unlist(vy_list[T_index_train + 1])
      ))#as.vector(v_all[,1+T_index_train])  
    } 
    if(is.null(ho_output)){
      ho_output=as.vector(rbind(
        unlist(vx_list[T_index_ho + 1]),
        unlist(vy_list[T_index_ho + 1])
      ))# as.vector(v_all[,1+T_index_ho])  
    }
    
    m_IKF=optim(param,pred_ho_output_Vicsek_variation_log_RMSE, control=list(maxit=200),
                #lower=c(-8,-8,-8), upper=c(5,1,2),
                kernel_type=kernel_type, neighbors_info=neighbors_info,grid_boundary_info=grid_boundary_info,
                pos_x_list=px_list, pos_y_list=py_list, vel_x_list=vx_list, vel_y_list=vy_list,
                T_index_train=T_index_train, T_index_ho=T_index_ho, output=output, ho_output=ho_output, 
                D_y=D_y, cut_r_max=cut_r_max, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte)
    
    #while(m_IKF$par[3]>log(10^6) | m_IKF$par[4]>log(10^6)){
    ##10^8 for the second one
    while(m_IKF$par[3]>log(10^6) | m_IKF$par[4]>log(10^6)){
      if(m_IKF$par[3]>log(10^6)){
        param[3]=param[3]+runif(1)
      }
      if(m_IKF$par[4]>log(10^6)){
        param[4]=param[4]+runif(1)
      }
      
      m_IKF=optim(param,pred_ho_output_Vicsek_variation_log_RMSE, control=list(maxit=200),
                  #lower=c(-8,-8,-8), upper=c(5,1,2),
                  kernel_type=kernel_type, neighbors_info=neighbors_info,grid_boundary_info=grid_boundary_info,
                  pos_x_list=px_list, pos_y_list=py_list, vel_x_list=vx_list, vel_y_list=vy_list,
                  T_index_train=T_index_train, T_index_ho=T_index_ho, output=output, ho_output=ho_output, 
                  D_y=D_y, cut_r_max=cut_r_max, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte)
    }
    
    
    param = m_IKF$par
  }
  
  beta_v=exp(param[1])
  beta_f=exp(param[2])
  tau_v=exp(param[3])  
  tau_f=exp(param[4])  
  threshold_r=exp(param[5])/(1+exp(param[5])) * cut_r_max
  
  parameters = c(beta_v, beta_f, tau_v, tau_f, threshold_r)
  names(parameters) = c('beta_v', 'beta_f', 'tau_v', 'tau_f', 'radius')
  
  
  
  ###for predicting
  # ans_neighbors_all=find_neighbors_Vicsek_variation_fast_grid(pos_x_list=px_list[T_index_time],pos_y_list=py_list[T_index_time], 
  #                                                             vel_x_list=vx_list[T_index_time],vel_y_list=vy_list[T_index_time], 
  #                                                             n_t=n_t,T_time=T_sim,D_y=D_y,grid_boundary_info=grid_boundary_info,cut_r=threshold_r)
  
  ans_neighbors_all=form_neighbors_Vicsek_variation_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                                           vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                                           grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,D_y=D_y)
  
  A_v_all_vec=ans_neighbors_all$A_v_vec
  A_f_all_vec=ans_neighbors_all$A_f_vec
  v_all_vec=ans_neighbors_all$v_vec
  d_pos_all_vec=ans_neighbors_all$d_pos_vec
  num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
  
  sort_v_all=sort(v_all_vec,index.return=T)
  N_tilde_v_all=length(v_all_vec) ###this is N_j in the paper
  delta_x_v_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_v_all]
  
  
  sort_f_all=sort(d_pos_all_vec,index.return=T)
  N_tilde_f_all=length(d_pos_all_vec) ###this is N_j in the paper
  delta_x_f_all=sort_f_all$x[-1]-sort_f_all$x[-N_tilde_f_all]
  
  output_all=as.vector(rbind(
    unlist(vx_list[2:(T_sim+1)]),
    unlist(vy_list[2:(T_sim+1)])
  ))#as.vector(v_all[,2:(T_sim+1)]) 
  
  m_CG=IKF_CG_particle_two_interact(param1=log(c(beta_v,tau_v)), param2=log(c(beta_f,tau_f)), 
                                    kernel_type1=kernel_type, kernel_type2=kernel_type, 
                                    delta_x_all1=delta_x_v_all, delta_x_all2=delta_x_f_all, 
                                    A_all_v1=A_v_all_vec, A_all_v2=A_f_all_vec, 
                                    sort_d_all_ix1=sort_v_all$ix, sort_d_all_ix2=sort_f_all$ix, 
                                    num_neighbors_vec1=2*num_neighbors_all_vec, num_neighbors_vec2=num_neighbors_all_vec, 
                                    output=output_all, tilde_nu=tilde_nu, 
                                    D=D_y, N=N, tol = tol, maxIte = maxIte)
  ans_CG=m_CG[[1]]
  sigma_2_0_est = output_all%*%ans_CG/length(output_all) ##sometimes negative? solved
  
  #first interaction
  ans_CG_v_tilde=ans_CG*tau_v ##this gives R_inv_y
  w_CG_v=A_t_times_x_particle(output=ans_CG_v_tilde, A_all_v=A_v_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,
                              D_y=D_y, N_tilde=N_tilde_v_all)
  #second interaction
  ans_CG_f_tilde=ans_CG*tau_f ##this gives R_inv_y
  w_CG_f=A_t_times_x_particle(output=ans_CG_f_tilde, A_all_v=A_f_all_vec,  num_neighbors_vec=num_neighbors_all_vec,
                              D_y=D_y, N_tilde=N_tilde_f_all)
  
  
  if(!is.null(testing_v_input) & !is.null(testing_d_input)){
    testing_n = length(testing_v_input)
    
    
    sigma_2_v_est=sigma_2_0_est*tau_v
    sigma_2_f_est=sigma_2_0_est*tau_f
    
    
    
    v_aug=c(testing_v_input,v_all_vec)
    v_aug_sort=sort(v_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
    v_aug_sort_x=v_aug_sort$x
    v_aug_sort_rev_ix=sort(v_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort
    
    delta_v_aug=v_aug_sort_x[2:length(v_aug_sort_x)]-v_aug_sort_x[1:(length(v_aug_sort_x)-1)]
    
    
    w_v_aug=c(rep(0,testing_n),w_CG_v)
    
    # param_tilde=log(c(beta_v,tilde_nu)) ##tilde nu is one to stablize the computation
    # pred_mean_v_aug=R_times_z(param_tilde, have_noise=T, delta_x=delta_v_aug, z=w_v_aug[v_aug_sort$ix],
    #                           kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
    
    pred_mean_v_aug = IKF(beta=beta_v, tilde_nu=tilde_nu, 
                          delta_x=delta_v_aug, output=w_v_aug[v_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
    
    pred_mean_v_fast=pred_mean_v_aug[v_aug_sort_rev_ix][1:testing_n]
    
    #NRMSE_v=sqrt(mean((pred_mean_v_fast-testing_v_output)^2))/sd(testing_v_output)
    
    
    
    f_aug=c(testing_d_input,d_pos_all_vec)
    f_aug_sort=sort(f_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
    f_aug_sort_x=f_aug_sort$x
    f_aug_sort_rev_ix=sort(f_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort
    
    delta_f_aug=f_aug_sort_x[2:length(f_aug_sort_x)]-f_aug_sort_x[1:(length(f_aug_sort_x)-1)]
    
    
    
    w_f_aug=c(rep(0,testing_n),w_CG_f)
    
    # param_tilde=log(c(beta_f,tilde_nu)) ##tilde nu is one to stablize the computation
    # pred_mean_f_aug=R_times_z(param_tilde, have_noise=T, delta_x=delta_f_aug, z=w_f_aug[f_aug_sort$ix],
    #                           kernel_type=kernel_type)-tilde_nu*w_f_aug[f_aug_sort$ix]
    
    pred_mean_f_aug = IKF(beta=beta_f, tilde_nu=tilde_nu, 
                          delta_x=delta_f_aug, output=w_f_aug[f_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_f_aug[f_aug_sort$ix]
    
    pred_mean_f_fast=pred_mean_f_aug[f_aug_sort_rev_ix][1:testing_n]
    
    
    if(compute_CI){
      #######variance
      #first kernel
      c_v_star=rep(NA,testing_n)
      r0_v=abs(outer(testing_v_input,(v_all_vec),'-'))
      if(kernel_type=='exp'){
        r_v = exp(-beta*r0_v)
      }else if(kernel_type=='matern_5_2'){
        r_v = matern_5_2_funct(r0_v, beta_v)
      }
      
      print("Computing the predictive variance for the first interaction ...")
      #system.time(
      for(i in 1:testing_n){
        #print(i)
        #A_r_i=A_times_x_particle(output= r[i,],  A_all_v=A_all_v,  num_neighbors_vec=num_neighbors_vec,  
        #                       D=D, N_tilde=N_tilde,  T_sim=T_sim,  S=S)
        
        A_r_v_i=A_times_x_particle(output=r_v[i,], A_all_v=A_v_all_vec,  num_neighbors_vec=2*num_neighbors_all_vec,  
                                   D=D_y, N)
        
        #tol=sd(A_r_i)^2*0.01*N_tilde ##can make it smaller
        tol_interval=tol*10^{-4}
        R_inv_r_v_all=IKF_CG_particle_two_interact(param1=log(c(beta_v,tau_v)), param2=log(c(beta_f,tau_f)), 
                                                   kernel_type1=kernel_type, kernel_type2=kernel_type, 
                                                   delta_x_all1=delta_x_v_all, delta_x_all2=delta_x_f_all, 
                                                   A_all_v1=A_v_all_vec, A_all_v2=A_f_all_vec, 
                                                   sort_d_all_ix1=sort_v_all$ix, sort_d_all_ix2=sort_f_all$ix, 
                                                   num_neighbors_vec1=2*num_neighbors_all_vec, num_neighbors_vec2=num_neighbors_all_vec, 
                                                   output=A_r_v_i, tilde_nu=tilde_nu, 
                                                   D=D_y, N=N, tol = tol_interval, maxIte = maxIte)
        
        # R_inv_r_all=fast_pred_sparse_CG( param=param,  kernel_type=kernel_type,   delta_x_all=delta_x_all,   output=A_r_i, 
        #                                  A_all_v=A_all_vec, sort_d_all$ix,  num_neighbors_vec=num_neighbors_all_vec, tilde_nu,
        #                                  D=D_y,   N=N,   tol=tol_interval,  maxIte = maxIte)
        R_inv_r_v=R_inv_r_v_all[[1]]*tau_v
        r_R_inv_r_v=A_r_v_i%*%R_inv_r_v
        
        c_v_star[i]=1-r_R_inv_r_v
        ##Ur=sparse_A_times_z( A,   P=permutation_A, r[i,])
        
      }
      c_v_star=abs(c_v_star)
      
      LB95_v=    pred_mean_v_fast+sqrt(as.numeric(sigma_2_v_est)*c_v_star)*qnorm(0.025)
      UB95_v=    pred_mean_v_fast+sqrt(as.numeric(sigma_2_v_est)*c_v_star)*qnorm(0.975)
      
      
      #second kernel
      c_f_star=rep(NA,testing_n)
      r0_f=abs(outer(testing_d_input,(d_pos_all_vec),'-'))
      if(kernel_type=='exp'){
        r_f = exp(-beta*r0_f)
      }else if(kernel_type=='matern_5_2'){
        r_f = matern_5_2_funct(r0_f, beta_f)
      }
      
      print("Computing the predictive variance for the second interaction ...")
      #system.time(
      for(i in 1:testing_n){
        #print(i)
        #A_r_i=A_times_x_particle(output= r[i,],  A_all_v=A_all_v,  num_neighbors_vec=num_neighbors_vec,  
        #                       D=D, N_tilde=N_tilde,  T_sim=T_sim,  S=S)
        
        A_r_f_i=A_times_x_particle(output=r_f[i,], A_all_v=A_f_all_vec,  num_neighbors_vec=num_neighbors_all_vec,  
                                   D=D_y, N)
        
        #tol=sd(A_r_i)^2*0.01*N_tilde ##can make it smaller
        tol_interval=tol*10^{-4}
        R_inv_r_f_all=IKF_CG_particle_two_interact(param1=log(c(beta_v,tau_v)), param2=log(c(beta_f,tau_f)), 
                                                   kernel_type1=kernel_type, kernel_type2=kernel_type, 
                                                   delta_x_all1=delta_x_v_all, delta_x_all2=delta_x_f_all, 
                                                   A_all_v1=A_v_all_vec, A_all_v2=A_f_all_vec, 
                                                   sort_d_all_ix1=sort_v_all$ix, sort_d_all_ix2=sort_f_all$ix, 
                                                   num_neighbors_vec1=2*num_neighbors_all_vec, num_neighbors_vec2=num_neighbors_all_vec, 
                                                   output=A_r_f_i, tilde_nu=tilde_nu, 
                                                   D=D_y, N=N, tol = tol_interval, maxIte = maxIte)
        
        # R_inv_r_all=fast_pred_sparse_CG( param=param,  kernel_type=kernel_type,   delta_x_all=delta_x_all,   output=A_r_i, 
        #                                  A_all_v=A_all_vec, sort_d_all$ix,  num_neighbors_vec=num_neighbors_all_vec, tilde_nu,
        #                                  D=D_y,   N=N,   tol=tol_interval,  maxIte = maxIte)
        R_inv_r_f=R_inv_r_f_all[[1]]*tau_f
        r_R_inv_r_f=A_r_f_i%*%R_inv_r_f
        
        c_f_star[i]=1-r_R_inv_r_f
        ##Ur=sparse_A_times_z( A,   P=permutation_A, r[i,])
        
      }
      
      c_f_star=abs(c_f_star)
      
      LB95_f=    pred_mean_f_fast+sqrt(as.numeric(sigma_2_f_est)*c_f_star)*qnorm(0.025)
      UB95_f=    pred_mean_f_fast+sqrt(as.numeric(sigma_2_f_est)*c_f_star)*qnorm(0.975)
      
    }
    
    #est_par_val=c(exp(m_IKF$par), exp(m_IKF$value))
    
  }
  
  
  
  
  
  
  # return(list(input_pos_all=input_pos_all,v_all = v_all,
  #             parameters = parameters, 
  #             pred_mean_v_fast = pred_mean_v_fast, pred_mean_f_fast = pred_mean_f_fast,
  #             LB95_v = LB95_v, UB95_v = UB95_v, LB95_f=LB95_f, UB95_f = UB95_f))
  new("particle.est",
      D_y = D_y,
      # data = list(
      #   positions = input_pos_all,
      #   velocities = v_all
      # ),
      parameters = parameters,  # This contains the estimated parameters
      sigma_2_0_est = sigma_2_0_est[1,1], 
      predictions = if(!is.null(testing_v_input) & !is.null(testing_d_input)){
        if(compute_CI) {
          list(mean_v = pred_mean_v_fast,lower95_v = LB95_v,upper95_v = UB95_v,
               mean_f = pred_mean_f_fast,lower95_f = LB95_f,upper95_f = UB95_f)
        } else {
          list(mean_v = pred_mean_v_fast,
               mean_f = pred_mean_f_fast)
        }
      }else{
        NULL
      },
      training_data = list(
        training_velocity = v_all_vec,
        training_distance = d_pos_all_vec,
        A_v = A_v_all_vec,
        A_f = A_f_all_vec,
        num_neighbors =num_neighbors_all_vec
      ),
      gp_weights = cbind(w_v=w_CG_v,w_f=w_CG_f) #matrix(w_CG, ncol = D_y)
  )
  
}


### general functions

simulate_particle = function(v_abs, n_t=100, T_sim=5, h=0.1, 
                             cut_r=0.5, sigma_0=0.1, noise_type = 'Gaussian', model = "Vicsek"){
  
  if (!model %in% c("Vicsek", "unnormalized_Vicsek", "two_interactions_Vicsek")){
    stop("Invalid model specified. Model must be either 'Vicsek', 'unnormalized_Vicsek', or 'two_interactions_Vicsek'")
  }
  if(model == "Vicsek"){
    sim = simulate_Vicsek(v_abs = v_abs, n_t = n_t, T_sim = T_sim, h = h,
                          cut_r = cut_r, sigma_0 = sigma_0, noise_type = noise_type)
  }else if(model == "unnormalized_Vicsek"){
    sim = simulate_unnormalized_Vicsek(v_abs = v_abs, n_t = n_t, T_sim = T_sim, h = h,
                                       cut_r = cut_r, sigma_0 = sigma_0, noise_type = noise_type)
  }else if(model == "two_interactions_Vicsek"){
    sim = simulate_Vicsek_variation(v_abs = v_abs, n_t = n_t, T_sim = T_sim, h = h, 
                                    cut_r = cut_r, sigma_0 = sigma_0, noise_type = noise_type)
  }
  return(sim)
}



### cell
trajectory_data <- function(particle_data) {
  # Get time range
  T_start <- min(particle_data$time)
  T_end <- max(particle_data$time)
  T_time <- T_end - T_start 
  
  # Initialize lists
  px_list <- vector("list", T_time + 1)
  py_list <- vector("list", T_time + 1)
  vx_list <- vector("list", T_time + 1)
  vy_list <- vector("list", T_time + 1)
  theta_list <- vector("list", T_time + 1)
  n_record <- numeric(T_time)
  particle_tracking <- vector("list", T_time)
  
  # Pre-compute time indices for faster lookup
  time_indices <- split(seq_len(nrow(particle_data)), particle_data$time)
  
  # First, store all particle data for each time point
  for(t in T_start:T_end) {
    current_idx <- t - T_start + 1
    index_current <- time_indices[[as.character(t)]]
    
    current_data <- particle_data[index_current, ]
    px_list[[current_idx]] <- current_data$px
    py_list[[current_idx]] <- current_data$py
    vx_list[[current_idx]] <- current_data$vx
    vy_list[[current_idx]] <- current_data$vy
    theta_list[[current_idx]] <- atan2(current_data$vy, current_data$vx)
    if(t < T_end) {
      n_record[current_idx] <- length(index_current)
    }
    
  }
  
  # Then, create tracking information between consecutive frames
  for(t in T_start:(T_end-1)) {
    current_idx <- t - T_start + 1
    
    # Get indices for current and next frame
    index_current <- time_indices[[as.character(t)]]
    index_next <- time_indices[[as.character(t + 1)]]
    
    
    # Find shared particles
    particles_current <- particle_data$particleID[index_current]
    particles_next <- particle_data$particleID[index_next]
    shared_particleID <- intersect(particles_current, particles_next)
    
    if(length(shared_particleID) == 0) next
    
    # Get indices of shared particles in both frames
    match_current <- match(shared_particleID, particles_current)
    match_next <- match(shared_particleID, particles_next)
    
    # Store tracking information
    particle_tracking[[current_idx]] <- data.frame(
      t_idx = current_idx,
      particle_ids = shared_particleID,
      current_indices = match_current,
      next_indices = match_next
    )
  }
  
  # Create and return particle.data object
  new("particle.data",
      px_list = px_list,
      py_list = py_list,
      vx_list = vx_list,
      vy_list = vy_list,
      theta_list = theta_list,
      particle_tracking = particle_tracking,
      data_type = "experiment",
      n_particles = n_record,
      T_time = T_time,
      D_y = 1)
}

# # Utility function to get paired data from consecutive frames
get_consecutive_data <- function(data_obj, variable = c("vx", "vy", "px", "py", "theta")) {
  variable <- match.arg(variable)
  
  # Get the appropriate list based on the variable
  data_list <- switch(variable,
                      "vx" = data_obj@vx_list,
                      "vy" = data_obj@vy_list,
                      "px" = data_obj@px_list,
                      "py" = data_obj@py_list,
                      "theta" = data_obj@theta_list)
  
  if(is.null(data_list)) {
    stop(paste(variable, "data not available"))
  }
  
  T_time <- data_obj@T_time
  start_list <- vector("list", T_time)
  end_list <- vector("list", T_time)
  
  # For simulation data, it's straightforward
  if(data_obj@data_type == "simulation") {
    for(t in 1:T_time) {
      start_list[[t]] <- data_list[[t]]
      end_list[[t]] <- data_list[[t + 1]]
    }
  } else {
    # For experimental data, use particle tracking
    for(t in 1:T_time) {
      tracking <- data_obj@particle_tracking[[t]]
      if(!is.null(tracking)) {
        # Use tracking indices to get corresponding data
        current_frame_data <- data_list[[t]]
        next_frame_data <- data_list[[t + 1]]
        
        # Extract data in the correct order using tracking indices
        start_list[[t]] <- current_frame_data[tracking$current_indices]
        end_list[[t]] <- next_frame_data[tracking$next_indices]
      }
    }
  }
  
  list(start = start_list, end = end_list)
}


form_neighbors_cell_with_r = function(threshold_r,pos_x_list,pos_y_list,vel_x_list,vel_y_list,
                                      time_range,grid_boundary_info,neighbors_info,
                                      direction,apolar_vicsek = FALSE){
  
  n_record = sapply(pos_x_list, length)
  
  Lx_min=unname(grid_boundary_info$grid_info['Lx_min',])
  Lx_max=unname(grid_boundary_info$grid_info['Lx_max',])
  Ly_min=unname(grid_boundary_info$grid_info['Ly_min',])
  Ly_max=unname(grid_boundary_info$grid_info['Ly_max',])
  nx=unname(grid_boundary_info$grid_info['nx',])
  ny=unname(grid_boundary_info$grid_info['ny',])
  len_x=unname(grid_boundary_info$grid_info['len_x',])
  len_y=unname(grid_boundary_info$grid_info['len_y',])
  
  
  
  A_v_neighbor_record=rep(NA,sum(n_record)*15)
  v_neighbor_record=rep(NA,sum(n_record)*15)
  num_neighbors_vec=rep(NA,sum(n_record))  
  
  count=0
  
  for(i_t in 1:length(time_range)){
    t = time_range[i_t]
    pos_x_t=pos_x_list[[t]]
    pos_y_t=pos_y_list[[t]]
    vel_x_t=vel_x_list[[t]]
    vel_y_t=vel_y_list[[t]]
    m_grid_here = neighbors_info[[t]]
    n_t = n_record[t]
    
    index_start = ifelse(i_t==1, 0, sum(n_record[time_range][1:(i_t-1)]))
    
    for(i in 1:n_t){
      input_pos_i=as.vector(c(pos_x_t[i],pos_y_t[i]))
      input_vel_i=as.vector(c(vel_x_t[i],vel_y_t[i]))
      
      i_x=ceiling((input_pos_i[1]-Lx_min)/len_x)
      i_y=ceiling((input_pos_i[2]-Ly_min)/len_y)
      
      index_grid=(i_y-1)*nx+i_x
      d_vec_here_all=input_pos_i-as.matrix(m_grid_here[[index_grid]]$neighbor_pos)
      
      d_here=sqrt(colSums(d_vec_here_all^2))
      if(apolar_vicsek){
        index_neighbor=which(d_here<threshold_r)
        #index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>0)
        #if(length(index_same_v_direction)==0){  ##if v is (0,0), the above can be integer(0)
        #  index_same_v_direction=1:(length(m_grid_here[[index_grid]]$neighbor_vel)/2) ##then all may be counted as same direction
        #}
        index_same_v_direction=which(colSums(m_grid_here[[index_grid]]$neighbor_vel*input_vel_i)>=0)
        
        index_neighbor=intersect(index_neighbor,index_same_v_direction)
      }else{
        index_neighbor=which(d_here<threshold_r)
      }
      
      n_neighbor = length(index_neighbor)
      num_neighbors_vec[index_start+i]=n_neighbor
      
      
      if(direction == "x"){
        v_neighbor_record[count+(1:n_neighbor)]=m_grid_here[[index_grid]]$neighbor_vel[1,index_neighbor]
      }else if(direction == "y"){
        v_neighbor_record[count+(1:n_neighbor)]=m_grid_here[[index_grid]]$neighbor_vel[2,index_neighbor]
      }
      
      A_v_neighbor_record[count+(1:n_neighbor)]=rep(1/n_neighbor,n_neighbor)
      
      
      
      count=count+n_neighbor
      
      
      
    }
    
  }
  
  v_neighbor_record=v_neighbor_record[1:count]
  A_v_neighbor_record=A_v_neighbor_record[1:count]
  
  ans_list = list(A_v_neighbor_record=A_v_neighbor_record,v_neighbor_record=v_neighbor_record,
                  num_neighbors_vec=num_neighbors_vec)
  return(ans_list)
}



pred_ho_output_cell_log_RMSE = function(param, kernel_type, neighbors_info,grid_boundary_info,
                                        pos_x_list, pos_y_list, vel_x_list, vel_y_list,sigma_2_record,
                                        T_index_train, T_index_ho, output, ho_output, direction, apolar_vicsek,
                                        D_y, cut_r_max, tilde_nu, tol, maxIte){
  beta=exp(param[1]) 
  tau=exp(param[2])  
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  n_record = sapply(pos_x_list, length)
  
  # ##form neighbor of train
  # ans_neighbors_train=form_neighbors_cell_with_max_neighbor(v_max_neighbor_record=v_max_neighbor_record_train,
  #                                                           d_pos_max_vec=max_neighbors_apolar_train$d_pos_vec,
  #                                                           num_neighbors_max_vec=max_neighbors_apolar_train$num_neighbors_vec,
  #                                                           n_t_record=n_record[T_index_train],T_time=T_train,cut_r=threshold_r)
  
  ans_neighbors_train=form_neighbors_cell_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                                 vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_train,
                                                 grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,
                                                 direction=direction,apolar_vicsek = apolar_vicsek)
  
  ##form neighbor of hold-out 
  # ans_neighbors_ho=form_neighbors_cell_with_max_neighbor(v_max_neighbor_record=v_max_neighbor_record_ho,
  #                                                        d_pos_max_vec=max_neighbors_apolar_ho$d_pos_vec,
  #                                                        num_neighbors_max_vec=max_neighbors_apolar_ho$num_neighbors_vec,
  #                                                        n_t_record=n_record[T_index_ho],T_time=T_ho,cut_r=threshold_r)
  
  ans_neighbors_ho=form_neighbors_cell_with_r(threshold_r=threshold_r,pos_x_list=pos_x_list,pos_y_list=pos_y_list,
                                              vel_x_list=vel_x_list,vel_y_list=vel_y_list,time_range=T_index_ho,
                                              grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,
                                              direction=direction,apolar_vicsek = apolar_vicsek)
  
  N_train=D_y*sum(n_record[T_index_train])
  N_ho=D_y*sum(n_record[T_index_ho]) ##this one is N in testing
  
  ##sort vx_train
  sort_v_train=sort(ans_neighbors_train$v_neighbor_record,index.return=T) ###sorted vx
  N_v_tilde=length(ans_neighbors_train$v_neighbor_record) ###N_tilde, distances
  delta_v_train=sort_v_train$x[-1]-sort_v_train$x[-N_v_tilde]
  sort_v_train_ix=sort_v_train$ix
  
  ##form augmented samples for cross-validation
  N_v_ho_tilde=length(ans_neighbors_ho$v_neighbor_record)  #number of hold out distance
  v_aug=c(ans_neighbors_ho$v_neighbor_record,(ans_neighbors_train$v_neighbor_record))
  v_aug_sort=sort(v_aug,index.return=T) ##sort augmented samples, this will have N_aug_tilde log(N_aug_tilde) order?
  v_aug_sort_x=v_aug_sort$x
  v_aug_sort_rev_ix=sort(v_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort
  
  delta_v_aug=v_aug_sort_x[2:length(v_aug_sort_x)]-v_aug_sort_x[1:(length(v_aug_sort_x)-1)]
  
  ###finish construction, now start to predict
  
  m_CG=IKF_CG_particle_cell(param=param, kernel_type=kernel_type, delta_x_all=delta_v_train, output=output, 
                            A_all_v=ans_neighbors_train$A_v_neighbor_record, sort_d_all_ix=sort_v_train_ix, 
                            sigma_2_vec=sigma_2_record[T_index_train], num_neighbors_vec=ans_neighbors_train$num_neighbors_vec, 
                            tilde_nu=tilde_nu, 
                            D=D_y, n_t_record=n_record[T_index_train], tol = tol, maxIte = maxIte)
  
  
  
  ans_CG=m_CG[[1]]
  
  ##this change this back to original parameterization
  ans_CG_v_tilde=ans_CG*tau ##this gives R_inv_y
  
  ###z=A_t_sparse_times_x,
  w_CG_v=A_t_times_x_particle(output=ans_CG_v_tilde, A_all_v=ans_neighbors_train$A_v_neighbor_record,  num_neighbors_vec=ans_neighbors_train$num_neighbors_vec,
                              D_y=D_y, N_tilde=N_v_tilde)
  w_v_aug=c(rep(0,N_v_ho_tilde),w_CG_v)
  
  # param_here=log(c(beta,tilde_nu)) ##tilde nu is one to stablize the computation
  # pred_mean_v_aug=R_times_z(param_here, have_noise=T, delta_x=delta_v_aug, z=w_v_aug[v_aug_sort$ix],
  #                           kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
  pred_mean_v_aug = IKF(beta=beta, tilde_nu=tilde_nu, 
                        delta_x=delta_v_aug, output=w_v_aug[v_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_v_aug[v_aug_sort$ix]
  
  
  pred_mean_v_fast=pred_mean_v_aug[v_aug_sort_rev_ix][1:N_v_ho_tilde]
  
  
  ##can only observes output so cross-validation on output
  pred_mean_v_ho_output=A_times_x_particle( pred_mean_v_fast,  ans_neighbors_ho$A_v_neighbor_record,  ans_neighbors_ho$num_neighbors_vec,
                                            D_y,  N_ho)
  
  
  log_RMSE_ho=1/2*log(mean( (ho_output-pred_mean_v_ho_output)^2)) ##many pars should work as it contains noises
  
  #print(c(beta,tau,threshold_r,log_RMSE_ho))
  
  return(log_RMSE_ho)
  
}


particle_interaction_est_cell = function(data_obj, param, cut_r_max, est_param=TRUE, nx=NULL, ny=NULL, direction, 
                                         kernel_type='matern_5_2',tilde_nu=0.1, tol=1e-6, maxIte=1000, 
                                         output=NULL, ho_output = NULL, testing_input=NULL, compute_CI = TRUE, apolar_vicsek=FALSE){
  
  px_list = get_consecutive_data(data_obj, "px")$start
  py_list = get_consecutive_data(data_obj, "py")$start
  vx_pairs = get_consecutive_data(data_obj, "vx")
  vx_list = vx_pairs$start
  vx_end_list = vx_pairs$end
  vy_pairs = get_consecutive_data(data_obj, "vy")
  vy_list = vy_pairs$start
  vy_end_list = vy_pairs$end
  D_y = data_obj@D_y
  
  
  T_time = data_obj@T_time
  n_record = sapply(px_list,length)
  
  
  T_index_time = 1:T_time
  
  
  
  px_min=min(unlist(px_list))
  px_max=max(unlist(px_list))
  py_min=min(unlist(py_list))
  py_max=max(unlist(py_list))
  
  
  if(is.null(nx)){
    nx=floor((px_max-px_min)/cut_r_max)
  }else{
    if(cut_r_max>(px_max-px_min)/nx) nx=floor((px_max-px_min)/cut_r_max)
  }
  
  if(is.null(ny)){
    ny=floor((py_max-py_min)/cut_r_max)
  }else{
    if(cut_r_max>(py_max-py_min)/ny) ny=floor((py_max-py_min)/cut_r_max)
  }
  
  grid_boundary_info=get_boundary_grid(px_min=px_min,px_max=px_max,
                                       py_min=py_min,py_max=py_max,nx=nx,ny=ny)
  
  neighbors_info = find_grid_neighbors(pos_x_list=px_list,pos_y_list=py_list,
                                       vel_x_list=vx_list,vel_y_list=vy_list, 
                                       time_range=T_index_time, grid_boundary_info=grid_boundary_info)
  
  
  #sigma_2_record = rep(NA,T_time)
  
  
  if(est_param){
    T_index_ho=seq(5,T_time,5) ##every 5 use the last one as holdout
    T_index_train=(1:T_time)[-T_index_ho]
    
    T_train=length(T_index_train) ##length of training data
    T_ho=length(T_index_ho) ##hold out prediction
    
    if(direction == "x"){
      if(is.null(output)) output=unlist(vx_end_list[T_index_train])
      if(is.null(ho_output)) ho_output=unlist(vx_end_list[T_index_ho])
      sigma_2_record=sapply(vx_list, var)
      #testing_input=seq(min(unlist(vx_list[T_index_time])),max(unlist(vx_list[T_index_time])),length.out=testing_n)
    }else if(direction == "y"){
      if(is.null(output)) output=unlist(vy_end_list[T_index_train])
      if(is.null(ho_output)) ho_output=unlist(vy_end_list[T_index_ho])
      sigma_2_record=sapply(vy_list, var)
      #testing_input=seq(min(unlist(vy_list[T_index_time])),max(unlist(vy_list[T_index_time])),length.out=testing_n)
    }
    
    #param=log(c(1,50,50)) 
    m_IKF=optim(param,pred_ho_output_cell_log_RMSE, control=list(maxit=200),
                #lower=c(-8,-8,-8), upper=c(5,1,2),
                kernel_type=kernel_type, neighbors_info=neighbors_info,grid_boundary_info=grid_boundary_info,
                pos_x_list=px_list, pos_y_list=py_list, vel_x_list=vx_list, vel_y_list=vy_list,sigma_2_record=sigma_2_record,
                T_index_train=T_index_train, T_index_ho=T_index_ho, output=output, ho_output=ho_output, direction=direction,
                apolar_vicsek=apolar_vicsek, D_y=D_y, cut_r_max=cut_r_max, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte)
    
    param = m_IKF$par
  }
  
  
  beta=exp(param[1])
  tau=exp(param[2]) 
  threshold_r=exp(param[3])/(1+exp(param[3])) * cut_r_max
  
  parameters = c(beta, tau, threshold_r)
  names(parameters) = c('beta', 'tau', 'radius')
  
  # # prediction
  # ans_neighbors_all=form_neighbors_cell_fast_grid(pos_x_list=px_list[T_index_time],pos_y_list=py_list[T_index_time], 
  #                                                 vel_x_list=vx_list[T_index_time],vel_y_list=vy_list[T_index_time],direction,
  #                                                 n_record=n_record[T_index_time],T_time,grid_boundary_info,cut_r=threshold_r,apolar_vicsek=apolar_vicsek)
  
  if(direction == "x"){
    output_all=unlist(vx_end_list[T_index_time])
    sigma_2_record=sapply(vx_list, var)
    #testing_input=seq(min(unlist(vx_list[T_index_time])),max(unlist(vx_list[T_index_time])),length.out=testing_n)
  }else if(direction == "y"){
    output_all=unlist(vy_end_list[T_index_time])
    sigma_2_record=sapply(vy_list, var)
    #testing_input=seq(min(unlist(vy_list[T_index_time])),max(unlist(vy_list[T_index_time])),length.out=testing_n)
  }
  
  ans_neighbors_all=form_neighbors_cell_with_r(threshold_r=threshold_r,pos_x_list=px_list,pos_y_list=py_list,
                                               vel_x_list=vx_list,vel_y_list=vy_list,time_range=T_index_time,
                                               grid_boundary_info=grid_boundary_info,neighbors_info=neighbors_info,
                                               direction=direction,apolar_vicsek = apolar_vicsek)
  
  A_v_all_vec=ans_neighbors_all$A_v_neighbor_record
  v_all_vec=ans_neighbors_all$v_neighbor_record
  #d_pos_vec=ans_neighbors_all$d_pos_vec
  num_neighbors_all_vec=ans_neighbors_all$num_neighbors_vec
  #mean(num_neighbors_all_vec)
  
  N=D_y*sum(n_record[T_index_time])
  sort_v_all=sort(v_all_vec,index.return=T)
  N_tilde_all=length(v_all_vec) ###this is N_j in the paper
  
  delta_x_all=sort_v_all$x[-1]-sort_v_all$x[-N_tilde_all]
  
  m_CG=IKF_CG_particle_cell(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=output_all, 
                            A_all_v=A_v_all_vec, sort_d_all_ix=sort_v_all$ix, 
                            sigma_2_vec=sigma_2_record[T_index_time], num_neighbors_vec=num_neighbors_all_vec, tilde_nu=tilde_nu, 
                            D=D_y, n_t_record=n_record[T_index_time], tol = tol, maxIte = maxIte)
  ans_CG=m_CG[[1]]
  ans_CG_tilde=ans_CG*tau
  
  sigma_2_0_prop_est = output_all%*%m_CG[[1]]/length(output_all)
  
  w_CG=A_t_times_x_particle(output=ans_CG_tilde, A_all_v=A_v_all_vec,  num_neighbors_vec=num_neighbors_all_vec,  
                            D_y=D_y, N_tilde=N_tilde_all)
  
  
  if(!is.null(testing_input)){
    testing_n = length(testing_input)
    
    sigma_2_0_est=(as.numeric(sigma_2_0_prop_est)*sigma_2_record[T_index_time])
    sigma_2_est=sigma_2_0_prop_est*tau
    
    v_aug=c(testing_input,(v_all_vec))
    v_aug_sort=sort(v_aug,index.return=T)
    v_aug_sort_x=v_aug_sort$x
    v_aug_sort_rev_ix=sort(v_aug_sort$ix,index.return=T)$ix ###this is to reverse the previous sort 
    
    delta_x_aug=v_aug_sort_x[2:length(v_aug_sort_x)]-v_aug_sort_x[1:(length(v_aug_sort_x)-1)]
    
    
    w_aug=c(rep(0,testing_n),w_CG)
    
    ###this should go back to nu
    # param_tilde=log(c(beta,tilde_nu)) 
    # pred_mean_aug=R_times_z(param_tilde, have_noise=T, delta_x=delta_x_aug, z=w_aug[v_aug_sort$ix],
    #                         kernel_type=kernel_type)-tilde_nu*w_aug[v_aug_sort$ix]
    pred_mean_aug = IKF(beta=beta, tilde_nu=tilde_nu, 
                        delta_x=delta_x_aug, output=w_aug[v_aug_sort$ix], kernel_type=kernel_type)-tilde_nu*w_aug[v_aug_sort$ix]
    
    
    
    pred_mean_fast=pred_mean_aug[v_aug_sort_rev_ix][1:testing_n]
    #plot(testing_input, pred_mean_fast, type='l', main = paste(round(beta,2), round(tau,2), round(threshold_r,2)))
    
    
    if(compute_CI){
      c_v_star=rep(NA,testing_n)
      r0_v=abs(outer(testing_input,(v_all_vec),'-'))
      if(kernel_type=='exp'){
        r_v = exp(-beta*r0_v)
      }else if(kernel_type=='matern_5_2'){
        r_v = matern_5_2_funct(r0_v, beta)
      }
      
      N=length(num_neighbors_all_vec)
      A_r_v_rec=matrix(NA,N,testing_n)
      R_inv_A_r_v_rec=matrix(NA,N,testing_n)
      
      print("Computing the predictive variance ...")
      # #system.time(
      for(i in 1:testing_n){
        #print(i)
        
        A_r_v_i=A_times_x_particle(output=r_v[i,], A_all_v=A_v_all_vec,  num_neighbors_vec=num_neighbors_all_vec,
                                   D=D_y, N)
        A_r_v_rec[,i]=A_r_v_i
        #tol=sd(A_r_i)^2*0.01*N_tilde ##can make it smaller
        tol_interval=tol*10^{-14}
        R_inv_r_v_all=IKF_CG_particle_cell(param=log(c(beta,tau)), kernel_type=kernel_type, delta_x_all=delta_x_all, output=A_r_v_i,
                                           A_all_v=A_v_all_vec, sort_d_all_ix=sort_v_all$ix,
                                           sigma_2_vec=sigma_2_record, num_neighbors_vec=num_neighbors_all_vec, tilde_nu=tilde_nu,
                                           D=D_y, n_t_record=n_record[T_index_time], tol = tol_interval, maxIte = maxIte)
        R_inv_A_r_v_rec[,i]=R_inv_r_v_all[[1]]
        
        
        R_inv_r_v=R_inv_r_v_all[[1]]*tau
        r_R_inv_r_v=A_r_v_i%*%R_inv_r_v
        
        c_v_star[i]=1-r_R_inv_r_v
        
      }
      c_v_star=abs(c_v_star)
      
      #credible interval of the mean
      LB95 = pred_mean_fast+sqrt(as.numeric(sigma_2_est)*c_v_star)*qnorm(0.025)
      UB95 = pred_mean_fast+sqrt(as.numeric(sigma_2_est)*c_v_star)*qnorm(0.975)
    }
    
    
  }
  
  
  
  new("particle.est",
      D_y=D_y,
      # data = list(
      #   positions = input_pos_all,
      #   velocities = v_all,
      #   angles = theta_all
      # ),
      parameters = parameters,  # This contains the estimated parameters
      sigma_2_0_est = sigma_2_0_prop_est[1,1], 
      predictions = if(!is.null(testing_input)) {
        if(compute_CI) {
          list(mean = pred_mean_fast, lower95 = LB95, upper95 = UB95)
        } else {
          list(mean = pred_mean_fast)
        }
      } else {
        NULL
      },
      training_data =list(
        training_velocity = v_all_vec,
        A_v = A_v_all_vec,
        num_neighbors =num_neighbors_all_vec
      ),
      gp_weights = matrix(w_CG)
  )
  
}


extract_time_window <- function(data_obj, first_frame, last_frame) {
  # Input validation
  if (!inherits(data_obj, "particle.data")) {
    stop("Input must be a particle.data object")
  }
  
  # We need one more frame than the number of pairs
  if (first_frame < 1 || (last_frame + 1) > length(data_obj@px_list)) {
    stop("Time indices out of bounds")
  }
  
  if (first_frame >= last_frame) {
    stop("first_frame must be less than last_frame")
  }
  
  # Extract time window for each list
  # We need frames from first_frame to last_frame + 1 to have complete pairs
  px_list_new <- data_obj@px_list[first_frame:(last_frame + 1)]
  py_list_new <- data_obj@py_list[first_frame:(last_frame + 1)]
  vx_list_new <- data_obj@vx_list[first_frame:(last_frame + 1)]
  vy_list_new <- data_obj@vy_list[first_frame:(last_frame + 1)]
  
  # Handle optional theta_list
  theta_list_new <- if (!is.null(data_obj@theta_list)) {
    data_obj@theta_list[first_frame:(last_frame + 1)]
  } else {
    NULL
  }
  
  # Handle particle tracking
  particle_tracking_new <- if (!is.null(data_obj@particle_tracking)) {
    data_obj@particle_tracking[first_frame:last_frame]
  } else {
    NULL
  }
  
  # Extract relevant n_particles
  if (length(data_obj@n_particles) == 1) {
    # For simulation data
    n_particles_new <- data_obj@n_particles
  } else {
    # For experimental data
    n_particles_new <- data_obj@n_particles[first_frame:last_frame]  # Don't include last frame in n_particles
  }
  
  # Create new particle.data object
  new("particle.data",
      px_list = px_list_new,
      py_list = py_list_new,
      vx_list = vx_list_new,
      vy_list = vy_list_new,
      theta_list = theta_list_new,
      particle_tracking = particle_tracking_new,
      data_type = data_obj@data_type,
      n_particles = n_particles_new,
      T_time = last_frame - first_frame + 1,  # Number of pairs
      model = data_obj@model,
      sigma_0 = data_obj@sigma_0,
      radius = data_obj@radius,
      D_y = data_obj@D_y)
}



fit.particle.data = function(data, param, cut_r_max=1, est_param = TRUE, nx = NULL, ny = NULL,
                             # px_list = data@px_list, py_list = data@py_list, 
                             # vx_list = data@vx_list, vy_list = data@vy_list, theta_list = data@theta_list,
                             # n_t = data@n_particles,T_time = data@T_time, D_y=data@D_y, 
                             kernel_type = "matern_5_2", tilde_nu=0.1, tol=1e-6, maxIte=1000,
                             output=NULL, ho_output = NULL, 
                             testing_inputs, compute_CI = TRUE, num_interaction = (length(param)-1)/2,
                             data_type = data@data_type, model = data@model,  apolar_vicsek=FALSE, direction = NULL){
  if(data_type == "simulation"){
    if (!model %in% c("unnormalized_Vicsek", "two_interactions_Vicsek")){
      stop("Invalid model specified. Model must be either 'unnormalized_Vicsek' or 'two_interactions_Vicsek'")
    }
  }
  
  if(data_type == "experiment"){
    if (is.null(direction)){
      stop("Please specify the modeling direction ('x' or 'y')")
    }
  }
  
  
  if (!is.null(testing_inputs)) testing_inputs = matrix(testing_inputs, ncol = num_interaction)
  
  if(data_type == "simulation"){
    if(model == "unnormalized_Vicsek"){
      
      if (!is.null(testing_inputs)){
        testing_input = testing_inputs[,1]
      } else {
        testing_input = NULL
      }
      est = particle_interaction_est_unnormalized_Vicsek(data_obj=data, param=param, cut_r_max=cut_r_max, est_param=est_param, nx = nx, ny = ny,
                                                         kernel_type=kernel_type, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte, 
                                                         output=output, ho_output = ho_output, testing_input = testing_input, compute_CI = compute_CI)
      est@model = model
    }else if(model == "two_interactions_Vicsek"){
      if (!is.null(testing_inputs)){
        testing_v_input = testing_inputs[,1]
        testing_d_input = testing_inputs[,2]
      } else{
        testing_v_input = NULL
        testing_d_input = NULL
      }
      est = particle_interaction_est_Vicsek_variation(data_obj=data, param=param, cut_r_max=cut_r_max, est_param=est_param, nx = nx, ny = ny,
                                                      kernel_type=kernel_type,tilde_nu=tilde_nu, tol=tol, maxIte=maxIte,
                                                      output=output, ho_output = ho_output, testing_v_input=testing_v_input, testing_d_input=testing_d_input, compute_CI = compute_CI)
      est@model = model
    }
  }else if(data_type == "experiment"){
    if (!is.null(testing_inputs)) {
      testing_input = testing_inputs[,1]
    }else{
      testing_input = NULL
    }
    est = particle_interaction_est_cell(data_obj=data, param=param, cut_r_max=cut_r_max, est_param=est_param, nx = nx, ny = ny, direction=direction, 
                                        kernel_type=kernel_type,tilde_nu=tilde_nu, tol=tol, maxIte=maxIte, 
                                        output=output, ho_output = ho_output, testing_input=testing_input, compute_CI = compute_CI, apolar_vicsek=apolar_vicsek)
  }
  
  est@data_type = data_type
  est@num_interaction = num_interaction
  
  return(est)
}


