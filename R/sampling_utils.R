#' Sampler to sample parameters from a uniform distribution
#'
#' Funtion that will sample parameter values from a uniform distribution based
#' on latin hypercube sampling.
#'
#' @param n numeric for the number of samples. Default is set to 100.
#' @param paramdf dataframe for parameters names and min and max values to use
#'    in the uniform distribution to sample parameters.
#'
#' @return dataframe for sampled parameters.
#'
#' @details This function will sample from a uniform distribution 3 different
#'  parameter values of interest. Maximum and minimum values to sample from
#'  uniform distribution should be provided.
#' @export
sampler <- function(n = 100, paramdf){

  ## number of dimension (the number of parameters that we will estimate)
  d = length(names(paramdf))
  params <- randomLHS(n, d)
  params[,1] <- qunif(params[,1], min = paramdf[1,1], max = paramdf[2,1])
  params[,2] <- qunif(params[,2], min = paramdf[1,2], max = paramdf[2,2])
  params[,3] <- qunif(params[,3], min = paramdf[1,3], max = paramdf[2,3])



  params <- data.frame(params)
  colnames(params) <- names(paramdf)


  return (params)

}


#' Sampler to sample parameters from a uniform distribution
#'
#' Funtion that will sample parameter values from a uniform distribution based
#' on latin hypercube sampling.
#'
#' @param n numeric for the number of samples. Default is set to 100.
#' @param paramdf dataframe for parameters names and values
#'
#' @return dataframe for sampled parameters.
#'
#' @details This function will sample from a uniform distribution 2 different
#'  parameter values of interest. Maximum and minimum values to sample from
#'  uniform distribution should be provided.
#' @export
sampler2 <- function(n = 100, paramdf){

  d = length(names(paramdf)) # number of dimension (the number of parameters that we will estimate)
  params <- randomLHS(n, d)
  params[,1] <- qunif(params[,1], min = paramdf[1,1], max = paramdf[2,1])
  params[,2] <- qunif(params[,2], min = paramdf[1,2], max = paramdf[2,2])

  params <- data.frame(params)
  colnames(params) <- names(paramdf)


  return (params)

}
