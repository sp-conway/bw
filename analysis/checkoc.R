rm(list=ls())
ctd <- .26
tcd <- .24
tdc <- .2
cdt <- .19
dtc <- 0
dct <- 0
sum(c(ctd,tcd,tdc,cdt,dtc,dct))

bc <- ctd+cdt
bt <- tdc+tcd
wc <- dtc+tdc
wt <- cdt+dct
print(bc)
print(bt)
print(wc)
print(wt)
bc>bt
wc<wt
