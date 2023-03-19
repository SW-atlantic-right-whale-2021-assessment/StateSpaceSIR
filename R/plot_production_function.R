# Define models
# - Theta logistic model
surplus <- function(r = rmax, k = Kinit, n = 50, z = Zinit){ # Theta logistic
    r*n*(1-(n/k)^z)
}

# - Haider et al 2017 (improved upon Courchamp)
haider <- function(k = Kinit, nmin = pmin * Kinit, n = 50){
    nmin/(k-nmin) * (n/nmin - 1)
}

# - Logistic allee
logistic <- function(nmin = pmin * Kinit, n = 50){
    1/(1+exp(-n/nmin)) * 2 - 1
}

# - Lin and Li 2002
linli <- function(k = Kinit, nmin = pmin * Kinit, n = 50){
    n/k - nmin/k
}

# - Hilborn et al 2014
hilborn <- function(nmin = pmin * Kinit, n = 50){
    1-exp(log(0.5)*n/nmin)
}


plot_suplus_prod <- function(SIRlist, coolors = c("#941B0C", "#104F55"), file_name = NULL){

    # Set up objects
    Nvec <- seq(from = 0, to = 100, length.out = 100)
    surplus_prod_trajectories <- list()
    output_summary <- list()

    # Loop through SIR models
    for(i in 1:length(SIRlist)){

        # - No Allee
        if(SIRlist[[i]]$inputs$allee_model == 0 | is.null(SIRlist[[i]]$inputs$allee_model)){
            surplus_prod_trajectories[[i]] <- sapply(Nvec, function(x) surplus(r = SIRlist[[i]]$resamples_output$r_max, k = 100, n = x, z = SIRlist[[i]]$resamples_output$z))
        }

        # - Hilborn Allee
        if(SIRlist[[i]]$inputs$allee_model == 1){
            surplus_prod_trajectories[[i]] <- sapply(Nvec, function(x) surplus(r = SIRlist[[i]]$resamples_output$r_max, k = 100, n = x, z = SIRlist[[i]]$resamples_output$z) * hilborn(nmin = SIRlist[[i]]$resamples_output$P50 * 100, n = x))
        }

        # - Logistic Allee
        if(SIRlist[[i]]$inputs$allee_model == 1){
            surplus_prod_trajectories[[i]] <- sapply(Nvec, function(x) surplus(r = SIRlist[[i]]$resamples_output$r_max, k = 100, n = x, z = SIRlist[[i]]$resamples_output$z) * logistic(nmin = SIRlist[[i]]$resamples_output$P50 * 100, n = x))
        }

        # - Lin Li Allee
        if(SIRlist[[i]]$inputs$allee_model == 1){
            surplus_prod_trajectories[[i]] <- sapply(Nvec, function(x) surplus(r = SIRlist[[i]]$resamples_output$r_max, k = 100, n = x, z = SIRlist[[i]]$resamples_output$z) * linli(k = 100, nmin = SIRlist[[i]]$resamples_output$P50 * 100, n = x))
        }

        # - Haider Allee
        if(SIRlist[[i]]$inputs$allee_model == 1){
            surplus_prod_trajectories[[i]] <- sapply(Nvec, function(x) surplus(r = SIRlist[[i]]$resamples_output$r_max, k = 100, n = x, z = SIRlist[[i]]$resamples_output$z) * haider(k = 100, nmin = SIRlist[[i]]$resamples_output$P50 * 100, n = x))
        }


        # Extact summary statistics
        row_names <- c("mean", "median",
                       "2.5%PI", "97.5%PI",
                       "25%PI", "75%PI",
                       "min", "max")

        output_summary[[i]] <- matrix(nrow = length(row_names), ncol = length(Nvec))
        output_summary[[i]][1, ] <- colMeans(surplus_prod_trajectories[[i]])
        output_summary[[i]][2:6, ] <- apply(surplus_prod_trajectories[[i]], 2, function(x) quantile(x, probs= c(0.5,  0.025, 0.975, 0.25, 0.75)))
        output_summary[[i]][7, ] <- apply(surplus_prod_trajectories[[i]], 2, min)
        output_summary[[i]][8, ] <- apply(surplus_prod_trajectories[[i]], 2, max)

        output_summary[[i]] <- as.data.frame(output_summary[[i]])
        rownames(output_summary[[i]]) <- row_names
        colnames(output_summary[[i]]) <- paste0("P",Nvec)
    }


    # Get plotting dimensions
    minprod <- min(sapply(output_summary, function(x) min(x[3,])))
    maxprod <- max(sapply(output_summary, function(x) max(x[4,])))

    # Set up plot save
    # Plot trajectory
    for(j in 1: (1 + as.numeric(!is.null(file_name)) * 2)){

        # PNG
        if(j == 2){
            filename <- paste0(file_name, "_surplus_production_function", ".png")
            png( file = filename , width=7.5, height = 100 / 25.4, family = "serif", units = "in", res = 300)
        }

        # Plot suplus production function
        plot(y = NA, x = NA,
             ylim = c(minprod, maxprod),
             xlim = c(0, 100),
             xlab = NA, ylab = NA, font = 2)

        abline(h = 0, col = "grey60")
        mtext(side = 2, "Surplus production (N)", line = 1.6, font = 2, adj = 0.6)
        mtext(side = 1, "Depletion (%)", line = 1.6, font = 2)

        # - Loop across models
        for(i in 1:length(output_summary)){

            # - Credible interval
            polygon(
                x = c(Nvec,rev(Nvec)),
                y = c(output_summary[[i]][3, ],rev(output_summary[[i]][4, ])),
                col = adjustcolor(coolors[i], alpha = 0.2), border = NA) # 95% CI

            polygon(
                x = c(Nvec,rev(Nvec)),
                y = c(output_summary[[i]][5, ], rev(output_summary[[i]][6, ])),
                col = adjustcolor(coolors[i], alpha = 0.5), border = NA) # 50% CI

            # Median
            lines( x = Nvec, y = output_summary[[i]][2, ], lwd = 3, col = coolors[i]) # Median
        }

        legend("topleft", model_names,
               lty = rep(1, length(model_names)),
               lwd = rep(3, length(model_names)),
               col = coolors,
               bty = "n")

        if(j == 2){
            dev.off()
        }
    }

}



