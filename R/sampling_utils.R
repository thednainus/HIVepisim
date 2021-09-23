#' Sampler to sample parameters from a uniform distribution
#'
#' Funtion that will sample parameter values from a uniform distribution based
#' on latin hypercube sampling.
#'
#' @param n numeric for the number of samples. Default is set to 100.
#' @param paramdf <- dataframe for parameters names and values
#'
#' @return dataframe for sampled parameters.
#'
#' @details This function will sample from a uniform distribution parameter
#'  values for initial population size for infected individuals (init.pop.param),
#'  act rate parameter (act.rate.param), and infection probability parameters
#'  (inf.prob.paramter)
#'  Boundaries for each parameter:
#'  init.pop.param: min = 100, max = 10000
#'  act.rate.param: min = , max =
#'  inf.prob.parameter: min = 0.01, max = 1
#' @export
#'
#' @examples
#' #the example below will sample 100 parameter values for each b, initS, import,
#' and st from a uniform ditribution
#' parameters <- sampler()
sampler <- function(n = 100, paramdf){

  d = length(names(paramdf)) # number of dimension (the number of parameters that we will estimate)
  params <- randomLHS(n, d)
  params[,1] <- qunif(params[,1], min = paramdf$init_pop1_param[1],
                      max = paramdf$init_pop1_param[2]) # init_pop_param
  params[,2] <- qunif(params[,2], min = paramdf$act.rate_param[1],
                      max = paramdf$act.rate_param[2]) # act.rate_param
  #params[,3] <- qunif(params[,3], min = 0.01, max = 1) # inf.prob.param
  params[,3] <- qunif(params[,3], min = paramdf$inf.prob.param[1],
                      max = paramdf$inf.prob.param[2]) # inf.prob.param



  params <- data.frame(params)
  colnames(params) <- names(paramdf)


  return (params)

}
