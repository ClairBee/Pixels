
require(reshape)

# bivariate normal density function
bvn <- function(x, y, x0, y0, sig.x, sig.y, rho) {
    const <- 1 / (2 * pi * sig.x * sig.y * sqrt(1 - rho^2))
    rr <- 1 / (2 * (1 - rho^2))^-1
    
    z.x <- (x - x0) / sig.x
    z.y <- (y - y0) / sig.y
    
    const * exp(- rr * (z.x^2 + z.y^2 - (2 * rho * z.x * z.y)))  
}

# create some toy data
x <- seq(0, 500, 25); y <- seq(0, 500, 25)
z <- outer(x, y, bvn, x0 = 250, y0 = 250, sig.x = 100, sig.y = 100, rho = 0)

fdat <- setNames(melt(z), nm = c("x", "y", "z"))

gauss.2d.ls <- function(param, x, y, z) {
    x0 <- param["x0"]; y0 <- param["y0"]
    sig.x <- param["sig.x"]; sig.y <- param["sig.y"]; rho <- param["rho"]
    
    est <- bvn(x, y, x0, y0, sig.x, sig.y, rho)
    
    sum((est - z)^2, na.rm = T)
}

# run using starting values close to true values
g2d <- optim(c(x0 = 240, y0 = 240, sig.x = 80, sig.y = 80, rho = 0),
             gauss.2d.ls, x = fdat$x, y = fdat$y, z = fdat$z, method = "L-BFGS-B",
             lower = c(-Inf, 1, 1, 0, 0, 0), upper = c(Inf, 250, 250, Inf, Inf, 1))
g2d$par

#           x0           y0        sig.x        sig.y          rho 
# 2.400000e+02 2.400000e+02 8.000000e+01 8.000000e+01 1.985233e-20