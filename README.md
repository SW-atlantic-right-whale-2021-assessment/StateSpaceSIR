# StateSpaceSIR
This is a package to implement a state-space generalized logistic model in a Bayesian framework using a Sampling-Importance-Resampling (SIR) algorithm adapted following Romero et al. (2022). The code has been adapted from from Zerbini et al. (2011) and Zerbini et al (2019). The SIR algorithm follows that implemented by McAllister et al. (1994).

Code and data to run the analysis in Romero et al. (2022) can be found at: https://github.com/SW-atlantic-right-whale-2021-assessment

# Population dynamics
The population dynamics are modeled as follows:

N_(t+1)= (N_t + N_t * r_max * [1 - (N_t / K) ^z] - C_t * SLR_p(t)) * e_t

where Nt is the estimated population abundance in year t, K is the estimated population carrying capacity, z is the assumed shape parameter corresponding to the percentage of K at which maximum sustainable yield production is achieved, r_max is the estimated maximum population growth rate, C_t is the annual catch, SLR_p(t) is a correction factor for the period of years that includes year t to account for whales that were struck and lost, and e_t is the lognormally distributed mean 0 process error. The estimable parameters of this model are K, rmax, θ, and $\sigma^2$, where θ determines the true catch for the pre-modern era given uncertainty in landings numbers and $\sigma^2$ is the variance of the lognormal process error:

C_t = C_(t,min) + θ * (C_(t,max) - C_(t,min) )

where C_(t,min) is the minimum estimate of catch in year t, C_(t,max) is the maximum estimated catch in year t. The parameter is K is not assigned a prior. Rather, abundance was projected using a “backwards” approach, which avoids explicitly defining a prior for K by instead assigning a prior to a recent abundance, Nrecent and back-calculating the abundance trajectory.  Priors are therefore defined for a recent estimate of absolute abundance, rmax, SLRp(t), θ, and $\sigma^2$. Model likelihoods were constructed for the absolute abundance and relative indices of relative abundance data assuming log-normal distributions. Catchability coefficients for the indices of relative abundance are analytically integrated out to produce a marginal likelihoods.

More information on the model can be found in Zerbini et al. (2011), Zerbini et al (2019), and Romero et al (In Revision).

# Installation and usage
To install the model the following code can be ran:
```{r}
library(devtools)
devtools::install_github(repo = "SW-atlantic-right-whale-2021-assessment/StateSpaceSIR")
```

An example of how the model can be run is:
```{r}
library(StateSpaceSIR)

# Load data from 2011 assessment
data("Abs.Abundance.2005") # Absolute abundance data
data("Catch.data") # Catch data
data("Count.Data") # Count
data("Rel.Abundance") # Relative abundance

# Set up priors
prior_list <- make_prior_list(
    r_max = make_prior(runif, 0, 0.106), # Population growth rate
    K = make_prior(use = FALSE), # Carrying capacity, no prior because backwars method
    var_N = make_prior(rinvgamma(4,0.1)), # Variance of process error
    N_obs = make_prior(runif, 500, 20000), # Prior on a recent abundance estimate
    premodern_catch_sample = make_prior(runif, 0, 1), # Samples between the minimum and maximum premodern catch values
    z = make_prior(2.39)) # Shape parameter for generalized logistic 
    
# Set up catch multiplier for struck and loss rate
slr_prior <- make_multiplier_list(c_mult_1 = make_prior(rnorm, 1.2, 0.05)) # Normal prior with mean = 1.2 and sd = 0.05 

# Run SIR
sirMod <- HUMPBACK.SIR(file_name = NULL, # File name to save
    n_resamples = 1000,
    priors = prior_list,
    catch_multipliers = slr_prior,
    target.Yr = 2005, # Year of the target population estimate for the backwards method.
    num.haplotypes = 0,
    output.Yrs = c(2005, 2006), # Years to output
    abs.abundance = Abs.Abundance.2005,
    rel.abundance = Rel.Abundance,
    rel.abundance.key = TRUE,
    count.data = Count.Data,
    count.data.key = FALSE, # Do not include in likelihood
    growth.rate.obs = c(0.074, 0.033, TRUE), # Prior on growth rat
    growth.rate.Yrs = c(1995, 1996, 1997, 1998),
    catch.data = Catch.data,
    control = sir_control(threshold = 1e-20, progress_bar = TRUE))
    
```
The model can be tuned to acheive the desired resampling rate by changing the threshold.

Once the model has run, summary statistics and model fits can be exported as follows:
```{r}
# Summarize parameters including mean, median, probability intervals of estimated parameters   
resample_summary_base <- summary_sir(sirMod$resamples_output, 
                                    object = "Resample_Summary", file_name = NULL)

# Summary of trajectory including mean, median, probability intervals of annual abundance                                 
trajectory_summary_base <- summary_sir(sirMod$resamples_trajectories, 
                            object = "Trajectory_Summary", file_name = NULL)
                            
# Plot population trajector
plot_trajectory(sirMod,  file_name = NULL)

# Plot posterior parameter densities
plot_density(SIR = sirMod,  file_name = NULL, inc_reference = FALSE)

# Plot fits to indices
plot_ioa(SIR = sirMod,  file_name = NULL, ioa_names = c("Feeding ground", "Breeding ground") )

# Get summary table
summary_trable( sirMod, file_name = NULL )
```

# 2019 Example figures
Estimated population trajectory and time series of catches of WSA humpback whales from Zerbini et al. (2019). The gray line represents the model averaged median trajectory, and the dark and light shaded areas correspond, respectively, to the 50% and 95% confidence intervals. The dashed black line represents the median trajectory of the reference case scenario, and the red line represents the catches, with shaded areas corresponding to uncertainty in the pre-modern whaling catches. The model is fit to the absolute abundance estimates in 2005 (black dots with confidence interval) and the model predicted abundance estimates in the same year (offset gray dots with confidence interval).
![alt text](https://github.com/SW-atlantic-right-whale-2021-assessment/StateSpaceSIR/blob/master/data-raw/Example/2019_Zerbini/Reference_state_space_trajectory_summary.png "Population trajectory")

# References
McAllister, M. K., Pikitch, E. K., Punt, A. E., Hilborn, R. 1994. A Bayesian approach to stock assessment and harvest decisions using the sampling/importance resampling algorithm. Canadian Journal of Fisheries and Aquatic Sciences. 12, 2673-2687. 

Romero, M.A., Coscarella, M.A., Adams, G.D. et al. Historical reconstruction of the population dynamics of southern right whales in the southwestern Atlantic Ocean. Sci Rep 12, 3324 (2022). https://doi.org/10.1038/s41598-022-07370-6

Zerbini, A. N., Ward, E., Engel, M., Andriolo, A., Kinas, P. G. 2011. A Bayesian assessment of the conservation status of humpback whales (Megaptera novaeangliae) in the western South Atlantic Ocean (Breeding Stock A). J. Cetacean Res. Manage. (special issue 3). 131-144. 

Zerbini, A. N., Adams, G. D., Best, J. K., Clapham, P. J., Jackson, J. A., Punt, A. E. In Review. Assessing the Recovery of an Antarctic Predator from Historical Exploitation. Royal Society Open Science: Proceedings B.
