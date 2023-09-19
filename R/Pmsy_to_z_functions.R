

#' Function to convert between Pmsy and Z for Hilborn et al 2014 depensation function
#'
#' @param z
#' @param Pmsy
#' @param k
#' @param r
#' @param q
#'
#'
#' @description
#' article{Hilborn2014,
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
#'
#' @return
#' @export
#'
pmsy_z_hilborn <- function(z,Pmsy, k = 100, r = .1, q = .02) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*(1-exp(log(0.5)*n/(k*q))), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}


#' Function to convert between Pmsy and Z for Lin and Li 2002 depensation function
#'
#' @param z
#' @param Pmsy
#' @param k
#' @param r
#' @param q
#'
#'
#' @description
#' report{Lin2002,
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
#'
#' @return
#' @export
#'
pmsy_z_linli <- function(z,Pmsy, k = 100, r = .1, q = .02) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*(n/k - q), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}

#' Function to convert between Pmsy and Z for logistic depensation function
#'
#' @param z
#' @param Pmsy
#' @param k
#' @param r
#' @param q
#'
#'
#' @return
#' @export
#'
pmsy_z_logistic <- function(z,Pmsy, k = 100, r = .1, q = .02) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*(1/(1+exp(-n/(q*k))) * 2 - 1), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}


#' Function to convert between Pmsy and Z for Haider et al 2017 depensation function
#'
#' @param z
#' @param Pmsy
#' @param k
#' @param r
#' @param q
#'
#'
#' @description
#' article{Haider2017,
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
#'
#' @return
#' @export
#'
pmsy_z_haider <- function(z,Pmsy, k = 100, r = .1, q = .02) {
    n = Pmsy * k
    derivallee <- deriv(~r*n*(1-(n/k)^z)*((q*k)/(k-(q*k)) * (n/(q*k) - 1)), c("n","z","k","r","q"), function.arg = TRUE)
    check <- derivallee(n, z, k, r, q)
    with(attributes(check), gradient)[1]
}

uniroot(pmsy_z_haider, Pmsy=.7, lower=.1,upper=100)$root

