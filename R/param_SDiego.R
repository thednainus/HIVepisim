
# MSM -----------------------------------------------------------------

#' @title Epidemic Model Parameters
#'
#' @description Sets the epidemic parameters for stochastic network models
#'              simulated with \code{\link{netsim}} for EpiModelHIV
#'
#' @param netstats Target statistics and related network initialization data from
#'        the standard ARTnet workflow.
#'
#' @param time.unit Unit of time relative to one day.
#'
#' @param hiv.test.rate Mean probability of HIV testing per week for
#'        MSM.
#' @param test.window.int Length of the HIV test window period in weeks.
#' @param stage_prog_rate1 Rate per time unit spent in CD4+ compartment for stage 1 of HIV
#'        infection (Cori et al. 2015).
#' @param stage_prog_rate2 Rate per time unit spent in CD4+ compartment for stage 2 of HIV
#'        infection (Cori et al. 2015).
#' @param stage_prog_rate3 Rate per time unit spent in CD4+ compartment for stage 3 of HIV
#'        infection (Cori et al. 2015).
#'
#' @param f1 Probability of an individual entering stage 1 of HIV infection after
#'        initial HIV infection (acute and early HIV infection).
#' @param f2 Probability of an individual entering stage 2 of HIV infection after
#'        initial HIV infection (acute and early HIV infection).
#' @param f3 Probability of an individual entering stage 3 of HIV infection after
#'        initial HIV infection (acute and early HIV infection).
#' @param f4 Probability of an individual entering stage 4 of HIV infection after
#'        initial HIV infection (acute and early HIV infection).
#'
#' @param trans.r Baseline for transmission rate per act.
#'
#' @param ws0 Risk ratio of an individual in acute and early HIV infection (stage 0)
#'        will transmit to another individual.
#' @param ws1 Risk ratio of an individual in stage 1 will transmit to an HIV negative
#'        individual.
#' @param ws2 Risk ratio of an individual in stage 2 will transmit to an HIV negative
#'        individual.
#' @param ws3 Risk ratio of an individual in stage 3 will transmit to an HIV negative
#'        individual.
#' @param ws4 Risk ratio of an individual in stage 4 will transmit to an HIV negative
#'        individual.
#' @param wc1 Risk ratio of an undiagnosed HIV + individual will transmit to an HIV
#'        negative individual.
#' @param wc2 Risk ratio of an diagnosed and untreated HIV + individual will
#'        transmit to an HIV negative individual.
#' @param wc3 Risk ratio of an diagnosed and treated HIV + individual will
#'        transmit to an HIV negative individual.
#' @param wr1 Risk ratio of individuals in risk group 1 will transmit to an HIV
#'        negative individual.
#' @param wr2 Risk ratio of individuals in risk group 2 will transmit to an HIV
#'        negative individual.
#'
#'
#' @param tx.init.prob Probability per time step that an MSM who has
#'        tested positive will initiate treatment.
#' @param tx.halt.prob Probability per time step that an
#'        MSM who have started treatment will stop treatment.
#' @param tx.reinit.prob Probability per time step that an
#'        MSM who has stopped treatment will restart treatment.
#' @param aids.mr Mortality rate of persons in the AIDS stage who are currently
#'        off ART.
#'
#' @param a.rate Rate at which MSM enter the population.
#' @param arrival.age Age (in years) of new arrivals.
#'
#'
#' @param ... Additional arguments passed to the function.
#'
#' @return
#' A list object of class \code{param_msm}, which can be passed to
#' EpiModel function \code{netsim}.
#'
#' @keywords msm
#'
#' @references
#' Cori A, Pickles M, van Sighem A, et al.
#' CD4+ cell dynamics in untreated HIV-1 infection: overall rates, and effects
#' of age, viral load, sex and calendar time. AIDS. 2015;29(18):2435-2446.
#' doi:10.1097/QAD.0000000000000854
#'
#' @export
#'
param_msm_sdiego <- function(netstats,
                             time.unit = 7,

                             # Clinical
                             stage_prog_rate1 = 1/((3.32 * 365) / 7),
                             stage_prog_rate2 = 1/((2.7 * 365) / 7),
                             stage_prog_rate3 = 1/((5.5 * 365) / 7),
                             f1 = 0.76,
                             f2 = 0.19,
                             f3 = 0.05,
                             f4 = 0,
                             hiv.test.rate = 0.01325,
                             test.window.int = 21/7,
                             tx.init.prob = 0.092,
                             tx.halt.prob = 0.0102,
                             tx.reinit.prob = 0.00066,

                             # Baseline for transmission rate
                             trans.r = 0.0001,

                             # Transmission risk ratio
                             # by HIV stage of infection
                             ws0 = 1,
                             ws1 = 0.1,
                             ws2 = 0.1,
                             ws3 = 0.1,
                             ws4 = 0.3,
                             # Transmission risk ratio
                             # by care status (undiagnosed; diagnosed and untreated; and diagnosed an treated)
                             wc1 = 1,
                             wc2 = 0.5,
                             wc3 = 0.05,
                             # Transmission risk ratio
                             # by risk group
                             wr1 = 1,
                             wr2 = 10,

                             # HIV natural history
                             aids.mr = 1/((5.06 * 365) / 7),

                             # Demographic
                             a.rate = 0.00052,
                             arrival.age = 18,

                             ...) {

  p <- get_args(formal.args = formals(sys.function()),
                dot.args = list(...))

  class(p) <- "param.net"
  return(p)
}
