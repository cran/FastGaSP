#' Multiply two numbers
#'
#' @param param a vector of parameters. The first parameter is the natural logarithm of the inverse range parameter in the kernel function. If the data contain noise, the second parameter is the logarithm of the nugget-variance ratio parameter.
#' @param have_noise a bool value. If it is true, it means the model contains a noise.
#' @param delta_x a vector with dimension (num_obs-1) x 1 for the differences between the sorted input locations.
#' @param output a vector with dimension num_obs x 1 for the observations at the sorted input locations.
#' @param kernel_type A character specifying the type of kernel.
#' @return A list where the first value is the natural logarithm of the determinant of the correlation matrix and the second value is the estimated sum of squares.
#' @keywords internal
#' @export

Get_log_det_S2_r <- function(param,have_noise,delta_x,output,kernel_type){
  return(Get_log_det_S2(param,have_noise,delta_x,output,kernel_type))
}