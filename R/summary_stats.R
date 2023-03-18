
#' OUTPUT FUNCTION
#'
#' Function that provides a summary of SIR outputs including: mean, median, 95%
#' credible interval, 90% predicitive interval, max, and sample size.
#'
#' @param x A data.frame of mcmc samples.
#' @param object Name of the model object as specified by the user.
#' @param file_name name of a file to identified the files exported by the
#'   function. If NULL, will not save .csv file.
#'
#' @return Returns a data.frame with summary of SIR outputs
#' @export
#' @examples
#' x  <-  rnorm(1000, 5, 7)
#' y  <-  rnorm(1000, 6, 9)
#' df <- data.frame(x = x, y = y)
#' summary_sir( df , object = "example_summary")
summary_sir <- function(x, object = "USERDEFINED", file_name = "NULL") {

    # Change name if not supplied
    if(object == "USERDEFINED"){
        if(TRUE %in% grepl(pattern = "N_", names(x), ignore.case = FALSE)){object == "trajectory_summary"}
        if(TRUE %in% grepl(pattern = "r_max", names(x), ignore.case = FALSE)){object == "parameter_summary"}
    }

    row_names <- c("mean", "median",
                   "2.5%PI", "97.5%PI",
                   "5%PI", "95%PI",
                   "min", "max", "n")

    output_summary <- matrix(nrow = length(row_names), ncol = dim(x)[2])
    output_summary[1, ] <- sapply(x, mean)
    output_summary[2:6, ] <- sapply(x, quantile, probs= c(0.5, 0.025, 0.975, 0.05, 0.95))
    output_summary[7, ] <- sapply(x, min)
    output_summary[8, ] <- sapply(x, max)
    output_summary[9, ] <- sapply(x, length)

    output_summary <- data.frame(output_summary)
    names(output_summary) <- names(x)
    row.names(output_summary) <- row_names
    noquote(format(output_summary, digits = 3, scientific = FALSE))

    if(!is.null(file_name)){
        write.csv(output_summary,
                  paste0(file_name, "_", object, ".csv"))
    }

    list(object = object, date=Sys.time(), output_summary = output_summary)
}

#' Function to write tables of logistic model parameter and derived quantities for StateSpaceSIR similar to Table 5 from Zerbini et al (2011).
#'
#' @param SIR Resample summary from StateSpaceSIR
#' @param file_name Desired filename to where csv file will be saved. If NULL, will not save.
summary_table <- function( SIR, file_name = NULL){

  num.IA <- max(SIR$inputs$rel.abundance$Index)

  # Vars of interest
  years <- sort(c( SIR$inputs$target.Yr, SIR$inputs$output.Years))

  if(SIR$inputs$allee_model == 0){
  vars <- c("r_max", "K", "z", "Pmsy", "var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")
  vars_latex <- c("$r_{max}$", "$K$", "$z$", "$Pmsy$","$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years), paste0("$q_{flt", 1:num.IA, "}$"), paste0("$\beta_{q_{flt", 1:num.IA,"}}$"), "$sigma_q$")
  } else{

      vars <- c("r_max", "K", "z", "Pmsy", "P50","var_N", "Nmin", paste0("N", years), "Max_Dep", paste0("status", years), paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")
      vars_latex <- c("$r_{max}$", "$K$", "$z$", "$Pmsy$", "$P_{50}$","$sigma$", "$N_{min}$", paste0("$N_{", years, "}$"), "Max depletion", paste0("Depletion in ", years), paste0("$q_{flt", 1:num.IA, "}$"), paste0("$\beta_{q_{flt", 1:num.IA,"}}$"), "$sigma_q$")
  }


  pop_vars <- c("K", "Nmin", paste0("N", years))
  depletion_vars <- c("Max_Dep", paste0("status", years), paste0("q_IA1", 1:num.IA), paste0("q_IA2", 1:num.IA), "add_VAR_IA")

  results <- data.frame(matrix(NA, nrow = length(vars), ncol = 8))
  colnames(results) <- c("Parameter","Mean", "Median", "2.5% CI", "25% CI", "75% CI", "97.5% CI", "Unique")

  x <- SIR$resamples_output[,vars]
  x$var_N <- sqrt(x$var_N)

  # Get posterior of q
  q_posteriors <- list() # Each layer is an index

  # -- Determining the number of Indices of Abundance available
  rel.abundance <- SIR$inputs$rel.abundance
  indices <- unique(rel.abundance$Index)
  IA.yrs <- rel.abundance$Year
  N_hat <- SIR$resamples_trajectories[, paste0("N_", IA.yrs)] # Estimates of N within IOA years

  # -- Q2 for exponent
  q1_cols <- grep("q_IA1", colnames(SIR$resamples_output)) # Columns of resample Q estimates
  q1_est <- SIR$resamples_output[, q1_cols]
  q1_est <- as.matrix(q1_est, ncol = length(q1_cols))

  q2_cols <- grep("q_IA2", colnames(SIR$resamples_output)) # Columns of resample Q estimates
  q2_est <- SIR$resamples_output[, q2_cols]
  q2_est <- as.matrix(q2_est, ncol = length(q2_cols))

  # -- Make var-covar into wide and tall with cov = 0 for different indices
  rel.var.covar.tall <-  subset(rel.abundance, select = -c(Index,Year,IA.obs,IndYear))
  rel.var.covar.wide <- rel.var.covar.tall[which(rel.abundance$Index == 1),]
  rel.var.covar.wide <- rel.var.covar.wide[1:nrow(rel.var.covar.wide),1:nrow(rel.var.covar.wide)]

  rel.hess.wide <- solve(rel.var.covar.wide[1:nrow(rel.var.covar.wide), 1: nrow(rel.var.covar.wide)])

  if(num.IA>1){
    for(i in 2:length(unique(rel.abundance$Index))){
      var.cov.tmp <- as.matrix(rel.var.covar.tall[which(rel.abundance$Index == i),])
      var.cov.tmp <- var.cov.tmp[1:nrow(var.cov.tmp), 1:nrow(var.cov.tmp)]
      colnames(var.cov.tmp) <- NULL
      rownames(var.cov.tmp) <- NULL
      rel.var.covar.wide <- Matrix::bdiag(as.matrix(rel.var.covar.wide), var.cov.tmp)
      rel.hess.tall <- plyr::rbind.fill.matrix(rel.hess.tall, solve(var.cov.tmp))
    }
  }
  rel.var.covar.wide <- as.matrix(rel.var.covar.wide)

  # -- Loop through posterior draws
  for(j in 1:nrow(SIR$resamples_trajectories)){
    # -- Sample q
    q_posteriors_tmp <- exp(MASS::mvrnorm(
      n = 5,
      mu = as.numeric(log(rel.abundance$IA.obs/(N_hat[j,] ^ (1+q2_est[j,rel.abundance$Index]))) - diag(rel.var.covar.wide)/2),
      Sigma = rel.var.covar.wide))

    # q_est <- exp(sum(rel.hess.wide %*% as.numeric(log(rel.abundance$IA.obs/N_hat[j,] ^ (q2_est[j,rel.abundance$Index] + 1))))/(sum(rel.hess.wide))) # q_i

    # -- Assign to list
    for(i in indices){
      if(j == 1){
        q_posteriors[[i]] <- c(q_posteriors_tmp[,which(rel.abundance$Index == i)])
      } else {
        q_posteriors[[i]] <- c(q_posteriors[[i]], c(q_posteriors_tmp[,which(rel.abundance$Index == i)]))
      }
    }
  }


  # Get summary statistics
  results[,1] <- vars_latex
  results[,2] <- sapply(x, mean)
  results[,3:7] <- t(sapply(x, quantile, probs= c(0.5, 0.025, 0.25, 0.75, 0.975)))
  results[,8] <- sapply(x, function(x) length(unique(x)))

  # Update q for posteriors
  posterior_q_results <- data.frame(matrix(NA, nrow = num.IA, ncol = 8))
  posterior_q_results[,1] <- paste0("$p(q)_{flt", 1:num.IA, "}$")
  posterior_q_results[,2] <- round(sapply(q_posteriors, mean),3)
  posterior_q_results[,3:7] <- round(t(sapply(q_posteriors, quantile, probs= c(0.5, 0.025, 0.25, 0.75, 0.975))), 3)
  colnames(posterior_q_results) <- c("Parameter","Mean", "Median", "2.5% CI", "25% CI", "75% CI", "97.5% CI", "Unique")

  # Format things
  results[c(1,3:4),2:7] <- round(results[c(1,3:4),2:7], 3)
  results[5,2:7] <- round(results[5,2:7], 4)
  results[which(vars %in% depletion_vars),2:7] <- round(results[which(vars %in% depletion_vars),2:7], 3)
  results[which(vars %in% pop_vars),2:7] <- format(round(results[which(vars %in% pop_vars),2:7], 0),big.mark=",",scientific=FALSE)
  results[,8] <- format(round(results[,8], 0),big.mark=",",scientific=FALSE)
  results <- rbind(results, posterior_q_results)

  if(!is.null(file_name)){
    write.csv(results, file = paste0(file_name, "_summary_table.csv"))
  }
  return(results)
}


#' Bayes factor
#'
#' @param SIR SIR Fit model or list of SIR fit models
#' @param prior_probs prior probabilities of each model
#'
#' @return vector of bayes factors
bayes_factor <- function( SIR , prior_probs = NULL){
  # If it is a single SIR, make into a list
  if(class(SIR) == "SIR"){
    stop("Error: only one SIR model provided")
  }

  if(is.null(prior_probs)){
    prior_probs <- rep(1/length(SIR), length(SIR)) # Make uniform prior probs
  }

  # Get average likelihoods
  data_probs <- sapply(SIR, function(x) sum(x$resamples_output$Likelihood)/ length(x$resamples_output$Likelihood) )

  harmonic_mean <- sapply(SIR, function(x) (sum(x$resamples_output$Likelihood ^ -1)/ length(x$resamples_output$Likelihood)) ^-1 ) # Should not be used, facvors parameter rich models

  post_model <- data_probs * prior_probs

  bayes_factor <- post_model/sum(post_model)
  return(bayes_factor)
}




#' Weighted SIR model
#'
#' Function to create a weighted model using bayes factors
#'
#' @param SIR List of SIR fit models
#' @param bayes_factor vector of bayes_factors from \code{\link{bayes_factor}}
#'
#' @return an object of class "SIR"
#' @export
weight_model <- function(SIR, bayes_factor){

  weighted_SIR <- SIR[[1]]

  names_output <- colnames(weighted_SIR$resamples_output)

  # Size of vector
  n_samples <- nrow( weighted_SIR$resamples_trajectories )
  subs_samples <- rmultinom(n = 1, size = n_samples, prob = bayes_factor) # How many values to take from each SIR

  sample_ind <- 1
  for(i in 1:length(SIR)){

    names_resample_sub <- names(SIR[[i]]$resamples_output)
    random_rows <- sample(1:n_samples, size = subs_samples[i], replace = FALSE)

    weighted_SIR$resamples_output[sample_ind:(sample_ind + subs_samples[i] - 1), names_output[names_output %in% names_resample_sub] ] <- SIR[[i]]$resamples_output[random_rows,names_output[names_output %in% names_resample_sub]]
    weighted_SIR$resamples_trajectories[sample_ind:(sample_ind + subs_samples[i] - 1),] <- SIR[[i]]$resamples_trajectories[random_rows, ]
    weighted_SIR$catch_trajectories[sample_ind:(sample_ind + subs_samples[i] - 1),] <- SIR[[i]]$catch_trajectories[random_rows, ]
    sample_ind <- sample_ind + subs_samples[i]
  }

  return(weighted_SIR)
}