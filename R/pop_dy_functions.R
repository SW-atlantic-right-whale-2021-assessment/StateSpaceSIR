
#' PREDICTED GROWTH RATE
#'
#' \code{PRED.GROWTH.RATE} computes the predicted growth rate if such
#' information is available from an independent estimate rather than being
#' estimated from data. Growth rate is calculated as: $$r_{t_0 - t_{fin}}^{pred}
#' = \frac{ \sum_{t = t_0} ^{t_{fin - 1}} ln \left( \frac{N_{t+1}^{pred}} {
#' N_t^{pred}} \right) } { t_{fin} - t_0 } = \frac{ ln \left( N_{fin}^{pred}
#' \right) - ln \left( N_{0}^{pred} \right)} { t_{fin} - t_0 }$$ where
#' $N^{pred}$ is the model predicted population size, in numbers, at time $t$ or
#' $t+1$ in years, $t_0$ is the start year of the equation (1995 in Zerbini et
#' al. 2011), and $t_{fin}$ is the last year of the equation (1998 in Zerbini et
#' al. 2011).
#'
#' @param growth.rate.Yrs The years to be used for growth rate computation. 1995 - 1998 are used in Zerbini et al. 2011.
#' @param Pred_N Time series of predicted abundance, in numbers, from \code{\link{GENERALIZED_LOGISTIC}}.
#' @param start_yr The first year of the projection (assumed to be the first year in the catch series).
#'
#' @return A numeric scalar representing predicted growth rate.
#'
#' @export
#' @examples
#' growth.rate.Yrs  <-  c(1995:1998)
#' Pred_N <- c(1000, 1500, 1500, 2000)
#' start_yr  <-  1995
#' PRED.GROWTH.RATE(growth.rate.Yrs, Pred_N, start_yr=start_yr)
PRED.GROWTH.RATE <- function(growth.rate.Yrs, Pred_N, start_yr = start_yr) {
    ## Computing the growth rate years
    GR.Yrs <- growth.rate.Yrs - start_yr + 1
    Pred_N.GR <- Pred_N[GR.Yrs]

    ## FIXME Just return this line?
    Pred.GR <- (log(Pred_N.GR[length(Pred_N.GR)]) -
                    log(Pred_N.GR[1])) / (length(Pred_N.GR) - 1)

    Pred.GR
}

#' Computes the predicted rate of increase for a set of specified years for
#' comparison with trends estimated separately with any of the indices of
#' abundance or count data
#'
#' @param data Count data or relative abundance index to use
#' @param Pred_N Number of individuals predicted
#' @param start_yr Initial year
#'
#' @return Vector of rates of increase, one per index
#' @export
#'
COMPUTING.ROI <- function(data = data, Pred_N = Pred_N, start_yr = NULL) {
    num.indices <- max(data$Index)
    Pred.ROI <- rep(NA, num.indices)

    for (i in 1:num.indices) {
        index.ini.year <- (head(subset(data, Index == i)$Year, 1) - start_yr)
        index.final.year <- (tail(subset(data, Index == i)$Year, 1) - start_yr)
        elapsed.years <- index.final.year - index.ini.year

        Pred.ROI[i] <- exp((log(Pred_N[index.final.year]) -
                                log(Pred_N[index.ini.year])) /
                               (elapsed.years)) - 1
    }
    Pred.ROI
}

#' Calculate a target K for the bisection method
#'
#' @param allee_model Switch to determine the depensation model. 0 = no Allee effect; 1 = Hilborn et al 2014 P50 Allee Effect; 2 = Logistic Allee effect; 3 = Lin and Li 2002; 4 = Haider et al 2017.
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param K Pre-expoitation population size in numbers or biomass
#'   (depending on input).
#' @param N1 Population size in numbers or biomass at year 1 (generally
#'   assumed to be K).
#' @param z Generalized logistic shape parameter, determines population
#'   size where productivity is masimum (assumed to be 2.39 by the ISC
#'   SC).
#' @param P50 Parameter that determines the level of depensation.
#' @param num_Yrs The number of projection years. Set as the last year
#'   in the catchor abundance series whichever is most recent, minus the
#'   start year.
#' @param start_yr First year of the projection (assumed to be the first
#'   year in the catch series).
#' @param target.Pop Target population size.
#' @param catches Catch time series. Cannot include NAs,
#' @param proc_error The time series of lognormal process errors. Currently
#' does not handle NAs or 0s
#' @param MVP Minimum Viable Population Size; `4 * num.haplotypes`
#'
#' @return Vector of differences between predicted population and target
#'   population.
#' @export
#'
#' @examples
#' TARGET.K(r_max, K, N1, z, start_yr=start_yr, num_Yrs=bisection.Yrs,
#'          target.Pop=target.Pop, catches=catches, MVP=MVP)
TARGET.K <- function(allee_model, r_max, K, N1, z, P50,
                     num_Yrs, start_yr,
                     target.Pop, catches, proc_error,
                     MVP = 0) {

    Pred_N <- GENERALIZED_LOGISTIC(allee_model = allee_model,
                                   r_max = r_max,
                                   K = K,
                                   N1 = K,
                                   z = z,
                                   P50 = P50,
                                   start_yr = start_yr,
                                   num_Yrs = num_Yrs,
                                   catches = catches,
                                   proc_error = proc_error,
                                   MVP = MVP)
    Pred_N$Pred_N[num_Yrs] - target.Pop
}

#' LOGISTIC BISECTION
#'
#' Method of Butterworth and Punt (1995) where the prior distribution of the
#' current absolute abundance $N_{2005}$ and maximum net recruitment rate
#' \code{r_max} are sampled and then used to determine the unique value of the
#' population abundance $N$ in \code{start_yr} (assumed to correspond to
#' carrying capacity $K$). Requires \code{\link{TARGET.K}} and subsequent
#' dependencies.
#'
#' @param K.low Lower bound for $K$ when preforming the bisection method of Punt
#'   and Butterworth (1995). Default is 1.
#' @param K.high Upper bound for $K$ when preforming the bisection method of
#'   Punt and Butterworth (1995). Default is 500,000.
#' @param allee_model Switch to determine the depensation model. 0 = no Allee effect; 1 = Hilborn et al 2014 P50 Allee Effect; 2 = Logistic Allee effect; 3 = Lin and Li 2002; 4 = Haider et al 2017.
#' @param r_max The maximum net recruitment rate ($r_{max}$).
#' @param z The parameter that determines the population size where productivity
#'   is maximum (e.g. assumed to be 2.39 by the IWC SC for logistic).
#' @param P50 Parameter that determines the level of depensation.
#' @param num_Yrs The number of projection years. Set as the last year in the
#'   catch or abundance series, whichever is most recent, minus the
#'   \code{start_yr}.
#' @param start_yr The first year of the projection (assumed to be the first
#'   year in the catch series).
#' @param target.Pop A sample of the prior on population abundance $N$, in
#'   numbers, set as \code{sample.N.obs} sampled from \code{priors$N.obs}
#' @param catches The time series of catch in numbers or biomass. Currently does
#'   not handle NAs and zeros will have to input a priori for years in which
#'   there were no catches.
#' @param proc_error The time series of lognormal process errors. Currently
#' does not handle NAs or 0s
#' @param MVP The minimum viable population size in numbers or biomass. Computed
#'   as 3 * \code{\link{num.haplotypes}} to compute minimum viable population
#'   (from Jackson et al., 2006 and IWC, 2007).
#' @param tol The desired accuracy (convergence tolerance) of
#'   \code{\link{stats::uniroot}}.
#'
#' @return A numeric scalar of an estimate of  carrying capacity $K$.
#'
#' @export
#' @examples
#' LOGISTIC.BISECTION.K(K.low = 1, K.high = 100000, r_max = r_max, z = z,
#'                      num_Yrs = bisection.Yrs, start_yr = start_yr,
#'                      target.Pop = target.Pop, catches = catches, proc_error = proc_error, MVP = MVP,
#'                      tol = 0.001)
LOGISTIC.BISECTION.K <- function(K.low,
                                 K.high,
                                 allee_model,
                                 r_max,
                                 z,
                                 P50,
                                 num_Yrs,
                                 start_yr,
                                 target.Pop,
                                 catches,
                                 proc_error,
                                 MVP,
                                 tol = 0.001) {
    Kmin <- uniroot(TARGET.K,
                    tol = tol,
                    c(K.low,
                      K.high),
                    allee_model = allee_model,
                    r_max = r_max,
                    z = z,
                    P50 = P50,
                    num_Yrs = num_Yrs,
                    start_yr = start_yr,
                    target.Pop = target.Pop,
                    catches = catches,
                    proc_error = proc_error,
                    MVP = MVP)
    Kmin$root
}
