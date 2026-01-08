rm(list=ls())
library(multinomineq)
A <- matrix(c(-1, 1, 0,
              0, -1, 1),
            nrow = 2, byrow = TRUE)
b <- c(0, 0)
res <- bf_binom(k = c(16, 4, 2), n = c(40, 36, 15),
                A = A, b = b, M = 100000)
res

