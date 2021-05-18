#' Function to calculate the accumulated number of whales using previously estimated regression coefficients
#'
#' @param a Parameter from log-link GLM representing the intercept
#' @param b Parameter from log-link GLM representing the categorical effect of year
#' @param c Parameter from log-link GLM representing the effect of day
#' @param d Parameter from log-link GLM representing the effect of day^2
#' @param mu Mean of normal distribution representing time whale remains in area
#' @param sigma Standard deviation of normal distribution representing time whale remains in area
#' @param year Year for calculating accumulated number of whales
#' @param x Day of year to calculated accumulated number of whales: Assuming 320 because after 320 no new whales come in
#'
#' @export
#'
#' @examples
#' accum_fun(x = 320)
accum_fun <- function(a = -77.41,
                      b =  0.03234,
                      c = 0.1605,
                      d = -0.0003323,
                      mu = 60,
                      sigma = sqrt(60),
                      x = 320,
                      ReportList = FALSE){

    t <- 1:365 # Julian day
    W_t <- c(0, exp(a + b + c * t + d * t^2)) # Estimated number of whales on day t: W_0 = 0
    W_t[1:100] <- 0 # No Whales until April
    P_t <- dnorm(t, mu, sigma, FALSE) # Probability whale remains in area
    dW_t <- c(W_t[2:366] - W_t[1:365]) # delta w - 2:366 in the first term because W_t goes from day 0 to day 365
    A_x <- rep(0, 366) # Accumulated number of whales A_0 = 0

    # Loop through K
    # -- Should be as follows:
    # -- A_0 = 0
    # -- A_1 = dW_1 + p1 * A_0
    # -- A_2 = dW_2 + p1 * A_1 + p2 * A_0
    # -- A_3 = dW_3 + p1 * A_2 + p2 * A_1 + p3 * A_0
    for(k in 1:x){
        A_x[k+1] <- sum(dW_t[1:k]) + sum(P_t[1:k] * A_x[k:1])
    }

    if(ReportList){
        return(list(A_x = A_x, W_t = W_t, dW_t = dW_t))
    } else{
        return(A_x[x+1])
    }
}
