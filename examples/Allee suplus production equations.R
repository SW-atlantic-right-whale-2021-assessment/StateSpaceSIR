# Set up parameters
rmax = 0.2
pmin = 0.10 # Depletion threshold for allee effect
Kinit = 100

# Derive of theta logistic to find Pmsy
NmsyKz <- function(z,Pmsy) { 1-(z+1)*Pmsy^z }
Zinit <- uniroot(NmsyKz,Pmsy=0.6,lower=1,upper=100)$root

# Set up models
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

pmsy_z_logistic <- function(z,Pmsy, k = Kinit, r = rmax, q = pmin) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*(1/(1+exp(-n/(k*q))) * 2 - 1), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}

Zinit_log <- uniroot(pmsy_z_logistic,Pmsy=0.6,lower=0.1,upper=100)$root

# - Lin and Li 2002
linli <- function(k = Kinit, nmin = pmin * Kinit, n = 50){
  n/k - nmin/k
}

# - Courchamp et al 1999
courchamp <- function(k = Kinit, nmin = pmin * Kinit, n = 50){
  n/nmin - 1
}

# - Hilborn et al 2014
hilborn <- function(nmin = pmin * Kinit, n = 50){
 1-exp(log(0.5)*n/nmin)
}

pmsy_z_hilborn <- function(z,Pmsy, k = Kinit, r = rmax, q = pmin) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*(1-exp(log(0.5)*n/(k*q))), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}

Zinit_hill <- uniroot(pmsy_z_hilborn,Pmsy=0.6,lower=1,upper=100)$root

# - Threshold sensu Quinn and Collie 1990
threshold <- function(r = rmax, k = Kinit, n = 50, nmin = pmin * Kinit, z = Zinit){
  r*(n-nmin)*(1-((n-nmin)/(k-nmin))^z)
}

# - Dennis 1989 - Two parameter
dennis <- function(r = rmax, k = Kinit, n = 50, theta = pmin * Kinit, # Population size where fitness is halved
                   alpha = .1, # Strenth of Allee effect
                   z = Zinit){
  n*(r*(1-(n/k)^z)-(alpha*theta/(theta+n)))
}


# Plot
curve(surplus(n = x, z = Zinit), from = 0, to = 100, ylab = "Surplus production", xlab = "Depletion", lty = 1)
# curve(surplus(n = x) * haider(n = x), from = 0, to = 100, add = TRUE, lty = 2)
curve(surplus(n = x, z = Zinit_log) * logistic(n = x), from = 0, to = 100, add = TRUE, lty = 3)
# curve(surplus(n = x) * linli(n = x), from = 0, to = 100, add = TRUE, lty = 4)
# curve(surplus(n = x) * courchamp(n = x), from = 0, to = 100, add = TRUE, lty = 5)
curve(surplus(n = x, z = Zinit_hill) * hilborn(n = x), from = 0, to = 100, add = TRUE, lty = 6)
# curve(threshold(n = x), from = 0, to = 100, add = TRUE, lty = 7)
# curve(dennis(n = x), from = 0, to = 100, add = TRUE, lty = 8)
abline(h = 0, col = "grey")
abline(v = 60, col = "grey")

legend("topleft", bty = "n", legend = c("No Allee", "Allee", "Logistic Allee","Lin & Li", "Courchamp", "Hilborn", "Threshold", "Dennis"), lty = 1:8)


# Bibtex
#' @article{Dennis1989,
#'   author = {Brian Dennis},
#'   doi = {10.1111/j.1939-7445.1989.tb00119.x},
#'   issn = {08908575},
#'   issue = {4},
#'   journal = {Natural Resource Modeling},
#'   month = {9},
#'   pages = {481-538},
#'   title = {ALLEE EFFECTS: POPULATION GROWTH, CRITICAL DENSITY, AND THE CHANCE OF EXTINCTION},
#'   volume = {3},
#'   url = {https://onlinelibrary.wiley.com/doi/10.1111/j.1939-7445.1989.tb00119.x},
#'   year = {1989},
#' }
#' @article{Courchamp1999,
#'   author = {Franck Courchamp and Tim Clutton-Brock and Bryan Grenfell},
#'   doi = {10.1016/S0169-5347(99)01683-3},
#'   issn = {01695347},
#'   issue = {10},
#'   journal = {Trends in Ecology & Evolution},
#'   month = {10},
#'   pages = {405-410},
#'   title = {Inverse density dependence and the Allee effect},
#'   volume = {14},
#'   url = {https://linkinghub.elsevier.com/retrieve/pii/S0169534799016833},
#'   year = {1999},
#' }
#' @article{Haider2017,
#'   abstract = {Potential biological removal (PBR) is an approach used to calculate sustainable harvest and “take” limits for populations. PBR was originally derived assuming logistic growth while ignoring the effects of small population size (i.e., an Allee effect). We derived a version of PBR that includes an Allee effect (i.e., small population size or densities limiting population growth rates). We found that PBR becomes less conservative when it fails to consider an Allee effect. Specifically, sustainable harvest and take levels based upon PBR with an Allee effect were between approximately 51% and 66% of levels based upon PBR without an Allee effect. Managers and biologists using PBR may need to consider the limitations if an Allee effect may be present in the species being modeled. Considerations for Management. Based upon our finding, management considerations may include: acknowledging that populations under stress may also be subject to Allee effects; recognizing limitations of approaches such as PBR when applying them to small populations; and broadly considering the Allee effects when using population models for natural resource management.},
#'   author = {Humza S. Haider and Sarah C. Oldfield and Tiffany Tu and Rosa K. Moreno and Jay E. Diffendorfer and Eric A. Eager and Richard A. Erickson},
#'   doi = {10.1111/nrm.12133},
#'   issn = {19397445},
#'   issue = {3},
#'   journal = {Natural Resource Modeling},
#'   keywords = {Endangered Species Act,Indiana bat (Myotis sodalis),Marine Mammal Protection Act,cetacean,logistic growth,population assessment},
#'   month = {8},
#'   publisher = {Rocky Mountain Mathematics Consortium},
#'   title = {Incorporating Allee effects into the potential biological removal level},
#'   volume = {30},
#'   year = {2017},
#' }
#' @article{Hilborn2014,
# abstract = {Previous meta-analysis of spawner-recruit relationships suggested that depensatory behaviour is uncommon, and stocks pushed to low abundance are unlikely to suffer decreases in recruitment more severe than would be expected based on the decline in spawning stock. Using an updated database that has over 100 stocks that were depleted to less than 20% of their maximum observed stock size, we tested for depensatory behaviour in both total surplus production and recruitment and we also examined the probability of stock increase as a function of stock size and fishing pressure. The number of stocks that showed a significant improvement with depensatory models was less than that expected by chance. Hierarchical meta-analysis showed that the majority of the evidence was for no depensatory behaviour but could not rule out depensation at very low stock sizes. Stocks that are depleted to low abundance are expected to rebuild when fishing pressure is reduced if the environment has not changed but there is considerable evidence that the majority of fish stocks are impacted by changes in productivity regimes. Nevertheless, if stocks are very heavily depleted and fishing pressure is not reduced to quite low levels, the expected recovery time is both uncertain and long. Very low abundance should clearly be avoided for many reasons and the range of abundance where depensation cannot be ruled out is well below commonly adopted limit reference points.},
# author = {Ray Hilborn and Daniel J. Hively and Olaf P. Jensen and Trevor A. Branch},
# doi = {10.1093/icesjms/fsu035},
# issn = {10959289},
# issue = {8},
# journal = {ICES Journal of Marine Science},
# keywords = {Depensation,low density dynamic,rebuilding,recovery,regime changes},
# month = {10},
# pages = {2141-2151},
# publisher = {Oxford University Press},
# title = {The dynamics of fish populations at low abundance and prospects for rebuilding and recovery},
# volume = {71},
# year = {2014},
# }
#' @report{Lin2002,
#'   abstract = {In this study an attempt is made to investigate comprehensively the maximum sustainable yield (MSY) of Allee population dynamic system. The results show that: (1) because the net production rate of Allee system is a function of the shape parameter d, all Allee systems which possess different B (B=A/K, here K is the carrying capacity, and A the parameter of Allee effect) are adapted themselves to the environments by controlling their production rate. Allee effect does not change the interaction between the system and environment, but it will reduce the MSY; (2) the MSY of Allee system is approximately one-nineteenth of the height of the allometry curve of body size for all A/K50.01. However when the net production rate R takes the values from 3 to 15, the MSY of Allee system is also approximately 1/19 of the height of the allometry curve of body size for all A/K 50.001.},
#'   author = {Zhen-Shan Lin and Bai-Lian Li},
#'   journal = {Ecological Modelling},
#'   keywords = {Allee-effect,MSY wwwelseviercom/locate/ecolmodel,Population dynamics},
#'   pages = {1-7},
#'   title = {The maximum sustainable yield of Allee dynamic system},
#'   volume = {154},
#'   url = {www.elsevier.com/locate/ecolmodel},
#'   year = {2002},
#' }
#' @report{Quinn1990,
#'   abstract = {Under a threshold management policy, harvesting occurs at a constant rate but ceases when a population drops below a threshold. A simulation model of an age-structured population with stochastic recruitment was constructed with such a harvest policy with several threshold levels. Other factors were fishing mortality, recruitment, and initial biomass. The objective function was a weighted function of average yield and standard deviation over a planning horizon. First, we determined the optimal threshold given fishing mortality. Secondly, we determined optimal threshold and fishing mortality, simultaneously. In application to eastern Bering Sea pollock, a threshold management policy always increased average yield over a non-threshold policy. For the first problem, optimal threshold levels ranged from 20 to 30% of pristine biomass. For the second problem, each scenario had a unique threshold and fishing mortality, with fishing mortality slightly above the maximum sustainable yield (MSY) level and a threshold range of 25-58%. These results were robust in regard to other factors. Benefits of the threstpld policy were greater with a Ricker spawner-recruit model and with higher fishing mortality. t h e success of the threshold management policy is due to the relatively rapid rebuilding of a population to levels producing MSY. Une pditique de gestion bake sur un seuil permet de maintenir un taux de capture constant; tsutefois, il faut cesser la pecke 602s que Iteffectif se trsuve ssus ce seuil. On a mis au point un rnodde de simulation d'une population dont on connaissait la structure par 2ge et dont le recrutement ktait aleatsire; le rnsdele etait bas6 sur la technique de gestion mentionnee ci-dessus et pr6vcsyait plusieurs seuils. On a tenu compte, entre autres, de facteurs comme la mortalit4 due 2 la peche, te recruternent et la bismasse initiale. La fsnctisn objective etait une fonction ponder6e d'un rendernent moyen et d'un 6cart type au cours d'uwe p6riode 4onn4e de planifisation. On a d'abord determine le seuil optimal compte tenu de la mortalit6 due 3 la peche. Ensuite, on a determine? simultan4rnent le seuil optimal et la mortalit4 due 2 la p6che. bfapplicatiow de cette technique de gestion pour une population de goberges dans I'est de la mer de Bering a invariablement Bait accroitre le rendement rnoyen par rapport 2 une rn6thode de gestisn ne prevoyant aucun seuil. Dans le premier cas, ies seuils optimaux variaient de 20 2 30 % de la bismasse intouchke. Dans le deuxi6me cas, chaque sc6nario avait son propre seuil et son propre taux de mortalit6 par peche : la mortalit6 due 3 la peche depassait kggrement le rendement 6quilibr6 maximal (MSY) et la plage des seuils variait de 25 50 %. Ces r6sultats etaient rsbustes par rapport 21 d'autres facteurs. Les avantages d'une telle rnethode de gestisn 6taient plus nornbreux lorsqu'on employait un modgle geniteur-recrutement de Ricker et un taux 68ev6 de mortalit6 due 3 la peche, Le succPs de eette rn6thode de gestion est dO a la reconstitution assez rapide d'une population de rnaniere a ce qukelle produise un rendement kquilibrk maximal.},
#'   author = {Terrance J Quinn and Robert Fagen and Jie Zheng},
#'   journal = {Can. j. Fish. Aquat. Sci},
#'   pages = {2016-2829},
#'   title = {Threshold management policies for exploited populations},
#'   volume = {47},
#'   year = {1990},
#' }
#'

