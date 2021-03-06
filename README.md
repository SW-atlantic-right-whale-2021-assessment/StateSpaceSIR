# StateSpaceSIR
This is a package to implement a state-space generalized logistic model in a Bayesian framework using a Sampling-Importance-Resampling (SIR) algorithm adapted following Romero et al. (2022). The code has been adapted from from Zerbini et al. (2011 & 2019). The SIR algorithm follows that implemented by McAllister et al. (1994).

Code and data to run the analysis in Romero et al. (2022) can be found at: https://github.com/SW-atlantic-right-whale-2021-assessment

# Population dynamics
The population dynamics are modeled as follows:

N_(t+1)= (N_t + N_t * r_max * [1 - (N_t / K) ^z] - C_t * SLR_p(t)) * e_t

where Nt is the estimated population abundance in year t, K is the estimated population carrying capacity, z is the assumed shape parameter corresponding to the percentage of K at which maximum sustainable yield production is achieved, r_max is the estimated maximum population growth rate, C_t is the annual catch, SLR_p(t) is a correction factor for the period of years that includes year t to account for whales that were struck and lost, and e_t is the lognormally distributed mean 0 process error. The estimable parameters of this model are K, rmax, θ, and $\sigma^2$, where θ determines the true catch for the pre-modern era given uncertainty in landings numbers and $\sigma^2$ is the variance of the lognormal process error:

C_t = C_(t,min) + θ * (C_(t,max) - C_(t,min) )

where C_(t,min) is the minimum estimate of catch in year t, C_(t,max) is the maximum estimated catch in year t. The parameter is K is not assigned a prior. Rather, abundance was projected using a “backwards” approach, which avoids explicitly defining a prior for K by instead assigning a prior to a recent abundance, Nrecent and back-calculating the abundance trajectory.  Priors are therefore defined for a recent estimate of absolute abundance, rmax, SLRp(t), θ, and $\sigma^2$. Model likelihoods were constructed for the absolute abundance and relative indices of relative abundance data assuming log-normal distributions. Catchability coefficients for the indices of relative abundance are analytically integrated out to produce a marginal likelihoods.

More information on the model can be found in Romero et al. (2022).

# Installation and usage
To install the model the following code can be ran:
```{r}
library(devtools)
devtools::install_github(repo = "SW-atlantic-right-whale-2021-assessment/StateSpaceSIR")
```

# An example of how the base model from Romero et al. (2022) can be run is:
```{r}
################################################################################
# Load packages
################################################################################
library(StateSpaceSIR)
library(EnvStats)


################################################################################
# Read in data
################################################################################
# Downlowd R project from https://github.com/SW-atlantic-right-whale-2021-assessment
# Load R project

# -- Catch
sw_right_data<-read.delim("Data/datosModeloBallenasmiles2020Miles1648to2019.csv", sep=";",header=FALSE)   
names(sw_right_data)<- c("Year","CatchMin","CatchMax","Nt")

# Four periods of struck and loss rates (SLRs)
# - Period 1: 1648-1770: SLR = 1
# - Period 2: 1771-1850: SLR ~ N(1.6, 0.04)
# - Period 3: 1851-1973: SLR ~ N(1.09, 0.04)
# - Period 4: 1974-Present: SLR = 1
catch_list <- list(sw_right_data[which(sw_right_data$Year < 1771),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1771 & sw_right_data$Year <= 1850),1:3],
                   sw_right_data[which(sw_right_data$Year >= 1851 & sw_right_data$Year <= 1973),1:3],
                   sw_right_data[which(sw_right_data$Year > 1973),1:3])

# -- Absolute abundance
Abs.Abundance.2010 <- data.frame(Year = 2010, N.obs = 4245, CV.obs = 245/4245) # 2010: 4245 (SE: 245, 95% CI 3,765, 4,725).

# -- Relative abundance
sw_right_rel_abundance<-read.csv("Data/Accumulated_n_whales_1999_to_2019.csv") 

Rel.Abundance.SWRight <- data.frame(Index = rep(1, nrow(sw_right_rel_abundance)), 
                                    Year = sw_right_rel_abundance$Year, 
                                    IA.obs = sw_right_rel_abundance$A_xy_mu_sim) #Using 0.2 as a proxy
Rel.Abundance.SWRight = cbind(Rel.Abundance.SWRight, sw_right_rel_abundance[,paste0("X",1:17)]) # This is binding the MVNORM Var-Covar Matrix for the index


################################################################################
# Base model
################################################################################
dir.create("Model runs/Base")
file_name <- "Model runs/Base/Base"

sir_base <- list()
for(i in 1:2){ # i = 1, fit model; i = 2; run prior predictive check
  sir_base[[i]] <-  StateSpaceSIR(
    file_name = NULL,
    n_resamples = 20000,
    priors = make_prior_list(r_max =  make_prior(runif, 0, 0.11),
                             N_obs = make_prior(runif, 100, 10000),
                             var_N = make_prior(runif, 6.506055e-05, 6.506055e-05 * 10),
                             z = make_prior(use = FALSE),
                             Pmsy = make_prior(runif, 0.5, 0.8)),
    catch_multipliers = make_multiplier_list(
      make_prior(1),
      make_prior(rnorm, 1.60 , 0.04), 
      make_prior(rnorm, 1.09, 0.04),
      make_prior(1)),
    target.Yr = 2019,
    num.haplotypes = 24,
    output.Yrs = c(2021, 2030),
    abs.abundance = Abs.Abundance.2010,
    abs.abundance.key = TRUE,
    rel.abundance = Rel.Abundance.SWRight,
    rel.abundance.key = TRUE, # Indices of abundance
    count.data = Count.Data, # Not used
    count.data.key = FALSE, # Don't use count data
    growth.rate.obs = c(0.074, 0.033, FALSE), # Do not include growth rate
    growth.rate.Yrs = c(1995, 1996, 1997, 1998), # Not used
    catch.data = catch_list,
    control = sir_control(threshold = 1e-5, progress_bar = TRUE),
    realized_prior = ifelse(i == 1, "FALSE", "TRUE"))
}

# Output table of posterior samples of N_yr and parameters
resample_summary_reference <- summary_sir(sir_base[[1]]$resamples_output, object = "Resample_Summary", file_name = file_name)
trajectory_summary_reference <- summary_sir(sir_base[[1]]$resamples_trajectories, object = "Trajectory_Summary", file_name = file_name)

# Plot fitted trajectory
plot_trajectory(sir_base[[1]],  file_name = file_name)

# Plot prior predictive trajectory
plot_trajectory(sir_base[[2]],  file_name = paste0(file_name, "prior"))

# Plot densities
plot_density(SIR = list(sir_base[[1]]),  file_name = file_name,   priors = list(sir_base[[2]]), inc_reference = FALSE)

# Plot index fit
plot_ioa(sir_base[[1]],  file_name = file_name, ioa_names = NULL )

# Get estimated quantities
summary_table(sir_base[[1]],  file_name = file_name)
```

# 2022 Example figures
Estimated population trajectory and time series of catches of southern right whales from Romero et al. (2022). The solid blue line represents the median estimated trajectory of the population abundance (Ny), while the shaded areas correspond to the 50% and 95% credible intervals. The solid red line represents the average number of whaling catches as estimated by the catch parameter (π), while the red shaded areas correspond to the 50% and 95% credible intervals. The grey and black dots represent the estimated and observed, respectively, absolute abundance in 2010.
![alt text](https://github.com/SW-atlantic-right-whale-2021-assessment/RightwhaleRuns/blob/main/Model%20runs/Base/Base_trajectory_summary.png "Population trajectory")

# References
McAllister, M. K., Pikitch, E. K., Punt, A. E., Hilborn, R. 1994. A Bayesian approach to stock assessment and harvest decisions using the sampling/importance resampling algorithm. Canadian Journal of Fisheries and Aquatic Sciences. 12, 2673-2687. 

Romero, M.A., Coscarella, M.A., Adams, G.D. et al. Historical reconstruction of the population dynamics of southern right whales in the southwestern Atlantic Ocean. Sci Rep 12, 3324 (2022). https://doi.org/10.1038/s41598-022-07370-6

Zerbini, A. N., Ward, E., Engel, M., Andriolo, A., Kinas, P. G. 2011. A Bayesian assessment of the conservation status of humpback whales (Megaptera novaeangliae) in the western South Atlantic Ocean (Breeding Stock A). J. Cetacean Res. Manage. (special issue 3). 131-144. 

Zerbini, A. N., Adams, G. D., Best, J. K., Clapham, P. J., Jackson, J. A., Punt, A. E. In Review. Assessing the Recovery of an Antarctic Predator from Historical Exploitation. Royal Society Open Science: Proceedings B.
