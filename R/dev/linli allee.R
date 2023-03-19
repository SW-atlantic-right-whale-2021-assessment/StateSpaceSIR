pmsy_z_linli <- function(z,Pmsy, q = .02) {
    n = Pmsy * 100
    derivallee <- deriv(~0.1*n*(1-(n/100)^z)*(n/100 - q), c("n","z","q"), function.arg = TRUE)
    check <- derivallee(n, z, q)
    with(attributes(check), gradient)[1]
}

uniroot(pmsy_z_linli,Pmsy=.65, q = .00002, lower= .0001,upper=100)$root

pmsy_vec <- seq(from = 0.5, to = 0.8, length.out = 100)
q_vec <- seq(0, to = 0.02, length.out = 100)

mat <- matrix(NA, 100,100)

for(i in 1:100){
    for(j in 1:100){
        mat[i,j] = try(uniroot(pmsy_z_linli,Pmsy=pmsy_vec[i], q = q_vec[j], lower= .01,upper=100)$root, silent = TRUE)
    }
}

mat[mat == "Error in uniroot(pmsy_z_linli, Pmsy = pmsy_vec[i], q = q_vec[j], lower = 0.01,  : \n  f() values at end points not of opposite sign\n"] <- NA
contour(y = q_vec, x = pmsy_vec, z = mat)

curve(0.1*x*(1-(x/100)^.1)*(x/100 - .2), from = 0, to = 100)
