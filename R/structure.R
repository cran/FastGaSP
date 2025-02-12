## fgasp Class
setClass("fgasp", 		
         representation( 
           num_obs = "integer",          ## observations number
           have_noise="logical",            ## whether the data contain noise
           kernel_type="character",           ## information of the kernel
           ## data
           input="vector",              ## location of the sort input
           delta_x="vector",            ## distance between each sorted input
           output = "vector"           ## the observations, size nx1
         ),
)

## show
if(!isGeneric("show")) {
  setGeneric(name = "show",
             def = function(object) standardGeneric("show")
  )
}

# if (!isGeneric("fit")) {
#   setGeneric("fit", function(object, ...) standardGeneric("fit"))
# }
setGeneric("fit",
           def = function(object, ...) {
             standardGeneric("fit")
           }
)

setMethod("show", "fgasp",
          function(object){
            show.fgasp(object)
          }
)



## pred.fgasp Class
setClass("predictobj.fgasp", representation(
  num_testing="vector",     ##the number of testing input
  testing_input="vector",   ##sorted testing input
  param="vector",            ##param
  mean = "vector",          ##predictive mean 
  var="vector",              ##predictive variance
  var_data="logical"        ##whether to calculate the predictive variance of the data.
),
)

if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object,...) standardGeneric("predict")
  )
}

setMethod("predict", "fgasp",
          definition=function(param=NA,object, testing_input, var_data=TRUE,sigma_2=NULL){
            predict.fgasp(param=param,object=object,testing_input=testing_input,var_data=var_data,
                          sigma_2=sigma_2)
          }
)


# fmou class
setClassUnion("numericOrNULL", c("numeric", "NULL"))
setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))

setClass(
  "fmou",
  representation  = list(
    output = "matrix",          # observation data matrix
    d = "numericOrNULL",        # number of latent factor
    est_d = "logical",          # whether to estimate d
    est_U0 = "logical",         # whether to estimate factor loading matrix 
    est_sigma0_2 = "logical",   # whether to estimate variance of noise
    U0 = "matrixOrNULL",        # fixed factor loading matrix
    sigma0_2 = "numericOrNULL" # fixed variance of noise
  )
)

setMethod(
  "initialize", "fmou",
  function(.Object, output, d = NULL, est_d = FALSE,
           est_U0 = TRUE, est_sigma0_2 = TRUE, U0 = NULL, 
           sigma0_2 = NULL) {
    
    .Object@output <- output
    .Object@d <- d
    .Object@est_d <- est_d
    .Object@est_U0 <- est_U0
    .Object@est_sigma0_2 <- est_sigma0_2
    .Object@U0 <- U0
    .Object@sigma0_2 <- sigma0_2
    
    # Add validation checks here (from the original function)
    if (.Object@est_d == FALSE & is.null(.Object@d)) {
      stop("d should be provided if est_d=FALSE")
    }
    
    if (.Object@est_U0 == FALSE & is.null(.Object@U0)) {
      stop("U0 must be provided if it is fixed.")
    }
    
    if (!is.null(.Object@U0)) {
      if (!all(abs(t(.Object@U0) %*% .Object@U0 - diag(ncol(.Object@U0))) < 1e-10)) {
        stop("U0 must be an orthogonal matrix.")
      }
    }
    
    if(.Object@est_sigma0_2==F & is.null(.Object@sigma0_2)){
      stop("sigma0_2 must be provided if it is fixed.")
    }
    if(.Object@est_sigma0_2==F & (!is.null(.Object@sigma0_2))){
      if(.Object@sigma0_2<0){
        stop("sigma0_2 must be nonnegative.")
      }
    }
    
    if(!is.null(.Object@d)){
      if(!is.null(.Object@U0)){
        if(nrow(.Object@U0)!=nrow(.Object@output) | ncol(.Object@U0)!=.Object@d){
          stop("U0 must be a k*d matrix.")
        }
      }
    }
    return(.Object)
  }
)


setGeneric("fit.fmou", function(object,...) standardGeneric("fit.fmou"))

setMethod(
  "fit.fmou", "fmou",
  function(object,M=50, threshold=1e-4,track_iterations=FALSE,track_neg_log_lik=FALSE,
           U_init=NULL, rho_init=NULL, sigma2_init=NULL, d_ub=NULL) {
    k <- nrow(object@output)
    n <- ncol(object@output)
    
    if(!is.null(U_init) & !is.null(object@d)){
      if(nrow(U_init)!=k | ncol(U_init)!=object@d){
        stop("U_init must be a k*d matrix.")
      }
    }
    
    if(!is.null(rho_init) & !is.null(object@d)){
      if(length(rho_init)!=object@d){
        stop("rho_init is a d-dimensional vector.")
      }
    }
    if(!is.null(sigma2_init) & !is.null(object@d)){
      if(length(sigma2_init)!=object@d){
        stop("sigma2_init is a d-dimensional vector.")
      }
      if(any(sigma2_init<0)){
        stop("sigma2_init is a non-negative vector.")
      }
    }
    
    # estimate d
    if(object@est_d & is.null(object@d)){ 
      if(object@est_sigma0_2){
        # noise level is unknown, perform IC
        if(is.null(d_ub)){
          d_ub = k
        }
        IC_vals = rep(NA, d_ub)
        U_IC = svd(object@output)$u
        
        for(i_d in 1:d_ub){
          criteria_val = log(mean((object@output - U_IC[,1:i_d]%*%t(U_IC[,1:i_d])%*%object@output)^2)) + i_d*(k+n)/(k*n)*log(k*n/(k+n))
          IC_vals[i_d] = criteria_val
        }
        object@d = which.min(IC_vals)
      }
      else{
        # noise level is known, perform variance matching (VM) with binary search
        if(is.null(d_ub)){
          d_ub=k
        }
        if(!is.null(U_init)){
          d_ub = min(d_ub, ncol(U_init))
        }
        if(!object@est_U0){
          d_ub = min(d_ub, ncol(object@U0))
        }
        if(!is.null(rho_init)){
          d_ub = min(d_ub, length(rho_init))
        }
        if(!is.null(sigma2_init)){
          d_ub = min(d_ub, length(sigma2_init))
        }
        left = 1
        right = d_ub
        
        while(left <= right){
          mid = left + floor((right-left)/2)
          em_fit = fmou_cpp(object@output, d=mid, M=M, threshold=threshold,
                            est_U0=object@est_U0, est_sigma0_2=T, 
                            track_iterations=F, 
                            track_neg_log_lik=F,
                            U0 = object@U0[,1:mid], U_init=U_init[,1:mid], 
                            rho_init=rho_init[1:mid],
                            sigma2_init=sigma2_init[1:mid],sigma0_2=NULL)
          est_noise = em_fit$sigma0_2
          if(est_noise < object@sigma0_2){
            # overfit, reduce d
            right = mid-1
          }
          else{
            # underfit, increase d
            left = mid + 1
          }
        }
        object@d = mid
      }
    }
    # estimate parameters in fmou
    fmou_fit <- fmou_cpp(
      object@output, 
      d = object@d, 
      M = M, 
      threshold = threshold,
      est_U0 = object@est_U0, 
      est_sigma0_2 = object@est_sigma0_2, 
      track_iterations = track_iterations, 
      track_neg_log_lik = track_neg_log_lik,
      U0 = object@U0[,1:object@d], 
      U_init = U_init[,1:object@d], 
      rho_init = rho_init[1:object@d],
      sigma2_init = sigma2_init[1:object@d], 
      sigma0_2 = object@sigma0_2
    )
    return(fmou_fit)
  }
)

setGeneric("predict.fmou", function(object,...) standardGeneric("predict.fmou"))

setMethod(
  "predict.fmou", "fmou",
  function(object, param_est, step=1, interval=FALSE, interval_data=TRUE) {
    d = object@d
    res = predict_fmou(param_est, step=1, interval=interval, interval_data=interval_data)
    return(res)
  }
)


#### GPPCA


setClass(
  "gppca",
  representation  = list(
    input = "numeric",           # input
    output = "matrix",          # observation matrix
    d = "numericOrNULL",        # number of latent factor
    est_d = "logical",          # whether to estimate d
    shared_params = "logical",  # whether latent processes share parameters
    kernel_type = "character"  # kernel type of latent processes
  )
)

setMethod(
  "initialize", "gppca",
  function(.Object, input, output, d,est_d=FALSE, shared_params = TRUE, 
           kernel_type="matern_5_2") {
    
    .Object@input <- input
    .Object@output <- output
    .Object@d <- d
    .Object@est_d <- est_d
    .Object@shared_params <- shared_params
    .Object@kernel_type <- kernel_type
    return(.Object)
  }
)

setGeneric("fit.gppca", function(object,...) standardGeneric("fit.gppca"))

setMethod(
  "fit.gppca", "gppca",
  function(object, sigma0_2=NULL, d_ub=NULL) {
    res <- gppca_fit(object@input, object@output, object@d, object@est_d,
                     object@shared_params, object@kernel_type, sigma0_2=NULL,
                     d_ub=NULL)
    return(res)
  }
)


setGeneric("predict.gppca", function(object,...) standardGeneric("predict.gppca"))

setMethod(
  "predict.gppca", "gppca",
  function(object, param, A_hat, step=1, interval=FALSE, interval_data=TRUE) {
    res <- predict_gppca(param=param, A_hat=A_hat, input=object@input, step=step,
                         output=object@output, d=object@d, kernel_type=object@kernel_type,
                         shared_param=object@shared_params,interval=interval, interval_data=interval_data)
    return(res)
  }
)

### particle simulations

# Define particle simulation class
# setClass("particle.data",
#          slots = list(
#            positions = "matrix",  # Matrix of particle positions
#            velocities = "matrix",  # Matrix of particle velocities
#            angles = "matrixOrNULL", # Matrix of angles
#            n_particles = "numeric",  # Number of particles
#            T_time = "numeric", # Time steps
#            sigma_0 = "numeric",  # Noise variance sigma_0
#            radius = "numeric", # Interaction radius
#            model = "character",  # Model type
#            D_y = "numeric" # Dimension of output
#          ))



# setMethod("initialize", "particle.data",
#           function(.Object, 
#                    positions = matrix(),velocities = matrix(),angles = NULL, # only needed for Vicsek
#                    n_particles = numeric(),T_time = numeric(),
#                    sigma_0 = numeric(),radius = numeric(),model = "Vicsek",  # Default model
#                    D_y = numeric(), 
#                    ...) {
#             # Assign values to slots
#             .Object@positions <- positions
#             .Object@velocities <- velocities
#             .Object@angles <- angles
#             .Object@n_particles <- n_particles
#             .Object@T_time <- T_time
#             .Object@sigma_0 <- sigma_0
#             .Object@radius <- radius
#             .Object@model <- model
#             .Object@D_y <- D_y
#             
#             # Validation
#             if (!identical(dim(positions), dim(velocities))) {
#               stop("Dimensions of positions and velocities must match")
#             }
#             
#             validObject(.Object)
#             return(.Object)
#           })
# 
# 
# setMethod("show", "particle.data",
#           function(object) {
#             cat("Model:", object@model, "\n")
#             cat("Number of particles:", object@n_particles, "\n")
#             cat("Number of time frames:", object@T_time, "\n")
#             cat("sigma_0 (noise variance):", object@sigma_0, "\n")
#             cat("Interaction radius:", object@radius, "\n")
#           })


setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
setClass("particle.data",
         slots = list(
           # Core data storage
           px_list = "list",  # List of x positions at each time
           py_list = "list",  # List of y positions at each time
           vx_list = "list",     # List of x velocities at each time
           vy_list = "list",     # List of y velocities at each time
           theta_list = "listOrNULL",  # Optional list of angles at each time
           
           # Tracking information (primarily for experimental data)
           particle_tracking = "listOrNULL",  # List of data frames for particle mappings
           
           # Metadata
           data_type = "character",    # "simulation" or "experiment"
           n_particles = "numeric",    # Number of particles (constant for simulation)
           T_time = "numeric",         # Total time steps
           D_y = "numeric",           # Dimension of output
           
           # Simulation-specific parameters (NULL for experimental)
           model = "characterOrNULL",      # Model type (e.g., "Vicsek")
           sigma_0 = "numericOrNULL",      # Noise variance
           radius = "numericOrNULL"        # Interaction radius
         ))

# Initialize method
setMethod("initialize", "particle.data",
          function(.Object,
                   # For both types
                   px_list = list(), py_list = list(),
                   vx_list = list(), vy_list = list(),
                   theta_list = NULL,
                   n_particles = numeric(),
                   T_time = numeric(),
                   D_y = numeric(),  # Default dimension
                   
                   # For experimental data
                   particle_tracking = NULL,
                   data_type = c("simulation", "experiment"),
                   
                   # For simulation data
                   model = NULL,
                   sigma_0 = NULL,
                   radius = NULL,
                   ...) {
            
            # Match data_type argument
            data_type <- match.arg(data_type)
            
            # Assign core data
            .Object@px_list <- px_list
            .Object@py_list <- py_list
            .Object@vx_list <- vx_list
            .Object@vy_list <- vy_list
            .Object@theta_list <- theta_list
            
            # Assign metadata
            .Object@data_type <- data_type
            .Object@n_particles <- n_particles
            .Object@T_time <- T_time
            .Object@D_y <- D_y
            
            # Assign tracking information for experimental data
            .Object@particle_tracking <- particle_tracking
            
            # Assign simulation-specific parameters
            .Object@model <- model
            .Object@sigma_0 <- sigma_0
            .Object@radius <- radius
            
            # Validation
            if (length(px_list) != length(py_list) ||
                length(px_list) != length(vx_list) ||
                length(px_list) != length(vy_list)) {
              stop("All list lengths must match")
            }
            
            # if (data_type == "simulation") {
            #   # For simulation, check constant particle numbers
            #   n_particles_per_frame <- vapply(px_list, length, numeric(1))
            #   if (!all(n_particles_per_frame == n_particles)) {
            #     stop("Simulation must have constant number of particles")
            #   }
            #   
            #   # Check simulation parameters
            #   if (is.null(model) || is.null(sigma_0) || is.null(radius)) {
            #     stop("Simulation data requires model, sigma_0, and radius parameters")
            #   }
            # } else {
            #   # For experimental, check particle tracking if provided
            #   if (!is.null(particle_tracking) && length(particle_tracking) != (length(px_list) - 1)) {
            #     stop("particle_tracking length must be T_time - 1")
            #   }
            # }
            
            # Additional validation for theta if provided
            if (!is.null(theta_list) && length(theta_list) != length(px_list)) {
              stop("theta_list length must match other lists")
            }
            
            validObject(.Object)
            return(.Object)
          })

# Show method
setMethod("show", "particle.data",
          function(object) {
            cat("Particle Data Object\n")
            cat("-------------------\n")
            cat("Data type:", object@data_type, "\n")
            if (object@data_type == "simulation") {
              cat("Model:", object@model, "\n")
              cat("Time steps:", object@T_time, "\n")
              cat("Number of particles:", object@n_particles, "\n")
              cat("Sigma_0 (noise variance):", object@sigma_0, "\n")
              cat("Interaction radius:", object@radius, "\n")
            } else {
              cat("Time steps:", object@T_time, "\n")
              cat("Average number of particles:", 
                  round(mean(vapply(object@px_list, length, numeric(1))),2), "\n")
            }
            # cat("Dimension:", object@D_y, "\n")
            # if (!is.null(object@theta_list)) {
            #   cat("Theta information: Available\n")
            # }
            # if (!is.null(object@particle_tracking)) {
            #   cat("Particle tracking: Available\n")
            # }
          })


# Define particle estimation class
setClass("particle.est",
         slots = list(
           data_type = "character",
           model = "characterOrNULL",  # Model type
           D_y = "numeric", # Dimension of output
           num_interaction = "numeric", # Number of interactions
           # data = "list",  # Input data list
           parameters = "numeric",   # Estimated parameters (beta, tau, radius)
           sigma_2_0_est = "numeric",
           predictions = "listOrNULL",   # Predictions and intervals
           training_data = "list",  # Training data used in GP model
           gp_weights = "matrix"    # GP model weights for prediction
         ))

# Initialize method
setMethod("initialize", "particle.est",
          function(.Object,
                   data_type = c("simulation", "experiment"),
                   model = NULL,
                   D_y = numeric(),
                   num_interaction = numeric(),
                   parameters = numeric(),
                   sigma_2_0_est = numeric(),
                   predictions = NULL,
                   training_data = list(),
                   gp_weights = matrix(),
                   ...) {
            
            # Match data_type argument
            data_type <- match.arg(data_type)
            
            # Assign slots
            .Object@data_type <- data_type
            .Object@model <- model
            .Object@D_y <- D_y
            .Object@num_interaction <- num_interaction
            .Object@parameters <- parameters
            .Object@sigma_2_0_est <- sigma_2_0_est
            .Object@predictions <- predictions
            .Object@training_data <- training_data
            .Object@gp_weights <- gp_weights
            
            # # Validation
            # if (data_type == "simulation" && 
            #     !model %in% c("Vicsek", "two_interactions_Vicsek")) {
            #   stop("For simulation data, model must be either 'Vicsek' or 'two_interactions_Vicsek'")
            # }
            # 
            # if (data_type == "experiment" && !is.null(model)) {
            #   stop("For experimental data, model should be NULL")
            # }
            
            # if (length(parameters) != 2*num_interaction + 1) {
            #   stop("Length of parameters vector must be 2*num_interaction + 1")
            # }
            
            validObject(.Object)
            return(.Object)
          })

# Show method
setMethod("show", "particle.est",
          function(object) {
            cat("Particle Model Estimation\n")
            cat("-------------------------\n")
            cat("Data type:", object@data_type, "\n")
            
            if (object@data_type == "simulation") {
              cat("Model type:", object@model, "\n")
            }
            
            cat("Output dimension:", object@num_interaction, "\n")
            cat("\nEstimated parameters:\n")
            cat("  beta (inverse range):", object@parameters[1:object@num_interaction], "\n")
            cat("  tau (variance-noise ratio):", 
                object@parameters[(object@num_interaction + 1):(2*object@num_interaction)], "\n")
            cat("  radius:", object@parameters[2*object@num_interaction + 1], "\n")
          })



setMethod("fit", "particle.data",
          definition=function(object, param, cut_r_max=1, est_param = TRUE, nx=NULL, ny=NULL, 
                              # px_list = object@px_list, py_list = object@py_list, 
                              # vx_list = object@vx_list, vy_list = object@vy_list, theta_list = object@theta_list,
                              # n_t = object@n_particles,T_time = object@T_time, D_y=object@D_y, 
                              kernel_type = "matern_5_2", tilde_nu=0.1, tol=1e-6, maxIte=1000,
                              output=NULL, ho_output = NULL, 
                              testing_inputs=NULL, compute_CI = TRUE, num_interaction = (length(param)-1)/2, 
                              data_type = object@data_type, model= object@model,apolar_vicsek=FALSE, direction = NULL){
            
            fit.particle.data(data=object, param = param, cut_r_max = cut_r_max, est_param = est_param, nx = nx, ny = ny,
                              # px_list = px_list,py_list = py_list, 
                              # vx_list = vx_list, vy_list = vy_list, theta_list = theta_list,
                              # n_t = n_t,T_time = T_time, D_y=D_y, 
                              kernel_type = kernel_type, tilde_nu=tilde_nu, tol=tol, maxIte=maxIte,
                              output=output, ho_output = ho_output, 
                              testing_inputs = testing_inputs, compute_CI=compute_CI, num_interaction = num_interaction,
                              data_type = data_type, model = model, apolar_vicsek=apolar_vicsek, direction = direction)
          }
)

