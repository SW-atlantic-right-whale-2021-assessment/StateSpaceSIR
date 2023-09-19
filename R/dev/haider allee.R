pmsy_z_haider <- function(z,Pmsy, q = .02) {
    n = Pmsy * 100
    haider_fn <- function(n,z,q){0.1*n*(1-(n/100)^z)*((q*100)/(100-(q*100)) * (n/(q*100) - 1))}
    haider <- formula(~n*(1-(n/100)^z)*((q*100)/(100-(q*100)) * (n/(q*100) - 1)),  c("n","z","q"))
    derivallee <- deriv(haider, c("n","z","q"), function.arg = TRUE, hessian = TRUE)
    check <- derivallee(n, z, q)
    with(attributes(check), gradient)[1]
}

sample.z <- NA
sample.z <- uniroot(pmsy_z_haider, Pmsy = .7, q = .2, lower= .4,upper=1)$root


z_vec <- seq(from = 0.1, to = 15, length.out = 100)
pmsy_vec <- seq(from = 0.6, to = 0.8, length.out = 100)
q_vec <- seq(0, to = 0.2, length.out = 100)



mat2 <- matrix(NA, 100,100)
mat <- matrix(NA, 100,100)

for(i in 1:100){
    for(j in 1:100){
        check <- derivallee(70, z_vec[i], q_vec[j])
        matrix(with(attributes(check), hessian),3,3)[3,1]
        mat2[i,j] <- pmsy_z_haider(z = z_vec[i], Pmsy = 0.7, q = q_vec[j])
        mat[i,j] = try(uniroot(pmsy_z_haider,Pmsy=pmsy_vec[i], q = q_vec[j], lower= .01,upper=100)$root, silent = TRUE)
    }
}

library(plotly)
mat2 <- mat2 >= 0
mat[mat == "Error in uniroot(pmsy_z_haider, Pmsy = pmsy_vec[i], q = q_vec[j], lower = 0.01,  : \n  f() values at end points not of opposite sign\n"] <- -999
mat[mat == "Error in uniroot(pmsy_z_haider, Pmsy = pmsy_vec[i], q = q_vec[j], lower = 0.01,  : \n  f.lower = f(lower) is NA\n"] <- -999

plot_ly(y = z_vec, x = q_vec, z = mat2, type = "contour")
plot_ly(y = pmsy_vec, x = q_vec, z = mat, type = "contour")



# Second deriv of q given Pmsy
Pmsy = 0.7
haider <- formula(~0.1*n*(1-(n/100)^z)*((q*100)/(100-(q*100)) * (n/(q*100) - 1)),  c("n","z","q"))
derivallee <- deriv(haider, c("n","z","q"), function.arg = TRUE, hessian = TRUE)
derivallee
