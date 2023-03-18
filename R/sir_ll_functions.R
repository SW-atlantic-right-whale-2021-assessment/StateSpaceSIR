#' Compute analytic estimates of q, the scaling parameter between indices and
#' absolute population size
#'
#' @param rel.abundance Relative abundance index
#' @param add_CV Coefficient of variation
#' @param Pred_N Predicted population
#' @param start_yr Initial year
#' @param num.IA Index of abundance
#'
#' @return A numeric estimator for $q$.
#' @export
#'
CALC.ANALYTIC.Q <- function(rel.abundance, Pred_N, start_yr,
                            add_CV = 0, num.IA) {
    ## Vector to store the q values
    analytic.Q <- rep(NA, num.IA)

    for (i in 1:num.IA) {
        ## Subseting across each index of abundance
        IA <- rel.abundance[rel.abundance$Index == i,]
        ## Years for which IAs are available
        IA.yrs <- IA$Year-start_yr + 1
        ## Numerator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
        qNumerator <- sum((log(IA$IA.obs / Pred_N[IA.yrs])) /
                              (IA$Sigma * IA$Sigma + add_CV * add_CV))
        ## Denominator of the analytic q estimator (Zerbini et al., 2011 - eq. (3))
        qDenominator <- sum(1 / (IA$Sigma * IA$Sigma + add_CV * add_CV))
        ## Estimate of q
        analytic.Q[i] <- exp(qNumerator / qDenominator)
    }
    analytic.Q
}


#' Compute analytic estimates of q, the scaling parameter between indices and
#' absolute population size
#'
#' @param rel.abundance Relative abundance index
#' @param rel.var.covar Variance covariance (narrow)
#' @param rel.hess Hessian matrix of relative abundance (narrow)
#' @param add_CV Coefficient of variation
#' @param Pred_N Predicted population
#' @param start_yr Initial year
#' @param beta Density dependent catchability coefficient I = q*N^(1+beta)
#' @param num.IA Index of abundance
#'
#' @return A numeric estimator for $q$.
#' @export
#'
CALC.ANALYTIC.Q.MVLNORM <- function(rel.abundance, rel.var.covar, rel.hess, Pred_N, start_yr,
                                    add_CV = 0, beta, num.IA) {
    ## Vector to store the q values
    analytic.Q <- rep(NA, num.IA)

    for (i in 1:num.IA) {
        ## Subseting across each index of abundance
        IA <- rel.abundance[rel.abundance$Index == i,]
        HESS <- rel.hess[rel.abundance$Index == i,]
        HESS <- as.matrix(HESS[,1:nrow(HESS)])

        ## Years for which IAs are available
        IA.yrs <- IA$Year-start_yr + 1

        ## Estimate of q
        analytic.Q[i] <-  exp(sum(HESS %*% log(IA$IA.obs / (Pred_N[IA.yrs]^(1+beta[i])))) / sum(HESS))
    }
    analytic.Q
}

#' Compute the log-likelihood of indices of abundance assuming univariate lognormal
#'
#' @param rel.abundance Relative abundance
#' @param Pred_N Predicted population size
#' @param start_yr Initial year
#' @param q.values Scaling parameter
#' @param add.CV Coefficient of variation
#' @param log Boolean, return log likelihood (default TRUE) or
#'   likelihood.
#'
#' @return List of likelihood based on Zerbini et al. (2011) eq. 5 or using `dnorm`
#' @export
#'
LNLIKE.IAs <- function(rel.abundance, Pred_N, start_yr,
                       q.values, add.CV, log = TRUE) {
    loglike.IA1 <- 0
    IA.yrs <- rel.abundance$Year-start_yr + 1
    loglike.IA1 <- -sum(
        dlnorm( # NOTE: can be changed to dlnorm_zerb
            x = rel.abundance$IA.obs,
            meanlog = log( q.values[rel.abundance$Index] * Pred_N[IA.yrs] ),
            sdlog = rel.abundance$Sigma + add.CV,
            log))

    loglike.IA1
}


#' Compute the log-likelihood of indices of abundance assuming multivariate lognormal
#'
#' @param rel.abundance Relative abundance indices
#' @param rel.var.covar Variance covariance matrix
#' @param Pred_N Predicted population size
#' @param start_yr Initial year
#' @param q.sample.IA1 Scaling parameter
#' @param q.sample.IA2 Scaling parameter 2
#' @param add.CV Coefficient of variation [UNUSED]
#' @param log Boolean, return log likelihood (default TRUE) or
#'   likelihood.
#'
#' @return List of likelihood based on multivariate lognormal distribution
#' @export
#'
LNLIKE.MVLNORM.IAs <- function(rel.abundance, rel.var.covar, Pred_N, start_yr,
                               q.sample.IA1, q.sample.IA2, add.CV, log = TRUE) {

    loglike.IA1 <- 0
    IA.yrs <- rel.abundance$Year-start_yr + 1 # Starts at start year

    loglike.IA1 <- -sum(
        mvtnorm::dmvnorm(
            x = log(rel.abundance$IA.obs),
            mean = log( q.sample.IA1[rel.abundance$Index] * (Pred_N[IA.yrs] ^ (1 +  q.sample.IA2[rel.abundance$Index])) ) - 0.5 * diag(rel.var.covar), # Lognormal bias correction
            sigma = rel.var.covar,
            log))

    loglike.IA1
}


#' Predict indices of abundance
#'
#' @param SIR SIR object
#'
#' @return List of predicted indices based on Zerbini et al. (2011) eq. 5 or using `dnorm`
#' @export
#'
PREDICT.IAs <- function(rel.abundance, Pred_N, start_yr,
                        q.values, add.CV, log = TRUE) {
    loglike.IA1 <- 0
    IA.yrs <- rel.abundance$Year-start_yr + 1
    loglike.IA1 <- -sum(
        rlnorm( # NOTE: can be changed to dlnorm_zerb
            x = rel.abundance$IA.obs,
            meanlog = log( q.values[rel.abundance$Index] * Pred_N[IA.yrs] ),
            sdlog = rel.abundance$Sigma + add.CV,
            log))

    loglike.IA1
}

#' LOG LIKELIHOOD OF ABSOLUTE ABUNDANCE
#'
#' This function computes two estimates of the log-likelihood of the estimated
#' absolute abundance using the equation from Zerbini et al. 2011 (eq. 4) and a
#' lognormal distribution from \code{\link{CALC.LNLIKE}}.
#'
#' @param Obs.N Observed absoluted abundance in numbers as a data.frame
#'   containing year, estimate of absolute abundance, and standard deviation
#' @param Pred_N Predicted absolute abundance in numbers from
#'   \code{\link{GENERALIZED_LOGISTIC}}.
#' @param start_yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param add_CV Additional CV to add to variance of lognormal distribution
#'   sampled from \code{priors$add_CV}.
#' @param log Return the log of the likelihood (TRUE/FALSE)
#'
#' @return A list of two numeric scalars of estimates of log-likelihood.
#'
#' @export
#' @examples
#' Obs.N  <-  data.frame(Year = 2005, Sigma = 5, Obs.N = 1000)
#' Pred_N  <-  1234
#' start_yr  <-  2005
#' LNLIKE.Ns(Obs.N, Pred_N, start_yr, add_cv = 0, log=TRUE)
LNLIKE.Ns <- function(Obs.N, Pred_N, start_yr, add_cv, log = TRUE) {
    N.yrs <- Obs.N$Year-start_yr+1
    nll_n <- -sum(
        dlnorm( # NOTE: can be changed to dlnorm_zerb
            x = Obs.N$N.obs,
            meanlog = log( Pred_N[N.yrs] ),
            sdlog = Obs.N$Sigma + add_cv,
            log))  ## Years for which Ns are available
    nll_n
}

#' Calculate the log-likelihood of the growth rate
#'
#' Calculates the log-likelihood of the estimated growth rate given the observed
#' growth rate and the standard deviation of the observed growth rate.
#'
#' @param Obs.GR Observed growth rate
#' @param Pred.GR Predicted growth rate
#' @param GR.SD.Obs Standard error of the observed growth rate
#'
#' @return A \code{list} containing \code{loglike.GR1} and \code{loglike.GR2}
#' @export
#' @examples
#' LNLIKE.GR(0.1, 0.1, 0.1)
LNLIKE.GR <- function(Obs.GR, Pred.GR, GR.SD.Obs, log = T) {
    ## This can be converted the likelihood from Zerbini et al. 2011 (eq. 6)
    -sum( dnorm( x = Obs.GR , mean = Pred.GR , sd = GR.SD.Obs, log = log ) )
}

#' Calculate the log-likelihood of the beach whale data
#'
#' @param beached_data Data on whale strandingsas a data.frame
#'   containing year, estimate of absolute abundance, and standard deviation
#' @param Pred_N Predicted absolute abundance in numbers from
#'   \code{\link{GENERALIZED_LOGISTIC}}.
#' @param start_yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param add_CV Additional CV to add to variance of lognormal distribution
#'   sampled from \code{priors$add_CV}.
#' @param q_anthro is the proportion of the population that will be killed each year from anthropogenic mortality
#' @param d_anthro is the detection probability of carcasses on the beach (including the probability of washing up on shore)
#' @param p_anthro is the proportion of the range of the species covered by monitoring
#' @param log Return the log of the likelihood (TRUE/FALSE)
#'
#' @export
#' @return the negative log-likelihood
LNLIKE.BEACHED <- function(beached_data,
                           Pred_N,
                           start_yr,
                           q_anthro,
                           d_anthro,
                           p_anthro,
                           log=TRUE){

    N.yrs <- beached_data$Year-start_yr+1
    nll_n <- -sum(
        dpois( # NOTE: can be changed to dlnorm_zerb
            x = beached_data$N.obs,
            lambda = q_anthro * d_anthro * p_anthro * Pred_N[N.yrs],
            log))  ## Years for which Ns are available
    nll_n
}

#' Function to calculate the log-likelihood using a lognormal distribution
#'
#' @param Obs.N Time series of observed abundance
#' @param Pred_N Time series of estimated abundance
#' @param CV coefficient of variation
#' @param log whether to export as log-likelihood
#'
#' @return returns a scalar of the likelihood
#' @export
#' @examples
#' Obs.N <- 2000
#' Pred_N <- 2340
#' CV <- 4
#' CALC.LNLIKE(Obs.N, Pred_N, CV)
CALC.LNLIKE <- function(Obs.N, Pred_N, CV, log = FALSE) {
    sum(dnorm(x = log(Obs.N), mean = log(Pred_N), sd = CV, log = log))
}
