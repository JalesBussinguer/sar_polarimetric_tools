library(maxLik)
 
samples <- read.csv("D:/thesis_data/VEG_INDICES/samples/stratified/campestre/20m/FC_20170112_20m_patches.csv", sep=",")

GI0.Estimator.m1m2 <- function(z, L) {

    m1 <- mean(z)
    m2 <- mean(z**2)
    m212 <- m2 / m1**2

    a <- -2 - (L + 1) / (L * m212)
    g <- m1 * (2 + (L + 1) / (L * m212))

    return(list("alpha" = a, "gamma" = g))
}

GI0.Estimator.m1m2(samples, 4)