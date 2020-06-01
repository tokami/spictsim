library(TMB)

compile("../src/thorson.cpp")
dyn.load(dynlib("thorson"))

data1 <- read.table("bh.dat", header=TRUE)
data2 <- read.table("bh0.dat",header=TRUE)


## no priors
data <- as.list(data1)
data$priors <- c(1.8085,-12.183)
data$priorsCov <- rbind(c(0.01619256, 0.03828887),
                        c(0.03828887, 0.10517698))
data$priorFlag = 0

param <- list()
param$loga <- 0
param$logb <- 0
param$logsigma <-0

obj <- MakeADFun(data, param, DLL="bh")
opt <- nlminb(obj$par, obj$fn, obj$gr)
summary(sdreport(obj))

obj$fn(opt$par)

## with priors
data <- as.list(data1)
data$priors <- c(1.8085,-12.183)
data$priorsCov <- rbind(c(0.01619256, 0.03828887),
                        c(0.03828887, 0.10517698))
data$priorFlag = 1

param <- list()
param$loga <- 0
param$logb <- 0
param$logsigma <-0

obj <- MakeADFun(data, param, DLL="bh")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- summary(sdreport(obj))

obj$fn(opt$par)


fitmcmc2 <- tmbstan(obj, chains=1,
                    iter=10000, init=list(opt$par),
                    lower=c(-50,-50,-50), upper=c(10,10,10))


## MCMC
require(tmbstan)

mc <- extract(fitmcmc2, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)

npar<-dim(mc)[3]
layout(rbind(rep(1,npar),c(2:(npar+1))))
matplot(mc[,1,], type="l")
for(i in 1:npar){
    pn <- attr(mc,"dimnames")$parameters[i]
    plot(density(mc[,1,i]), main=pn, col=i, lwd=3)
    abline(v=sdr[i,1], lwd=3)
    xx <- pretty(range(mc[,1,i]),100)
    lines(xx,dnorm(xx,sdr[i,1],sdr[i,2]),
          col=gray(.5), lwd=3)
}

mcP <- mc
denP <- kde2d(mcP[,1,1],mcP[,1,2], n=1000)
image(denP, col = terrain.colors(100), axes = FALSE)
contour(denP, levels = seq(0.1, 20, by = 1),
        add = TRUE, col = "black")



## integrated
data <- list()
data$logR <- c(data1$logR,data2$logR)
data$ssb <- c(data1$ssb,data2$ssb)
data$priors <- c(1.8085,-12.183)
data$priorsCov <- rbind(c(0.01619256, 0.03828887),
                        c(0.03828887, 0.10517698))
data$priorFlag = 0

param <- list()
param$loga <- 0
param$logb <- 0
param$logsigma <-0

obj <- MakeADFun(data, param, DLL="bh")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- summary(sdreport(obj))

obj$fn(opt$par)



## MCMC
require(tmbstan)

fitmcmc2 <- tmbstan(obj, chains=1,
                    iter=10000, init=list(opt$par),
                    lower=c(-50,-50,-50), upper=c(10,10,10))


mc <- extract(fitmcmc2, pars=names(obj$par),
              inc_warmup=TRUE, permuted=FALSE)

dev.new()
npar<-dim(mc)[3]
layout(rbind(rep(1,npar),c(2:(npar+1))))
matplot(mc[,1,], type="l")
for(i in 1:npar){
    pn <- attr(mc,"dimnames")$parameters[i]
    plot(density(mc[,1,i]), main=pn, col=i, lwd=3)
    abline(v=sdr[i,1], lwd=3)
    xx <- pretty(range(mc[,1,i]),100)
    lines(xx,dnorm(xx,sdr[i,1],sdr[i,2]),
          col=gray(.5), lwd=3)
}


mcI <- mc


plot(mc[,,1], mc[,,2], ty='p', xlim = c(1.4, 2.2),
     col="grey70",pch=16,
     ylim = c(-14, -11),xlab = "loga", ylab = "logb")
den <- kde2d(mc[,,1], mc[,,2])
contour(den, add=TRUE, levels  =  c(0.1, 0.5, 0.95,2,4) )


denI <- kde2d(mcI[,1,1],mcI[,1,2], n=1000)
image(denI, col = terrain.colors(100), axes = FALSE)
contour(denI, levels = seq(0.1, 20, by = 1),
        add = TRUE, col = "black")



par(mfrow = c(1,2))
image(denP, col = terrain.colors(100), axes = TRUE,ylim=c(-13,-11.5),xlim=c(1.56,2.2))
contour(denP, levels = seq(2.1, 10, by = 1),
        add = TRUE, col = "black")
image(denI, col = terrain.colors(100), axes = TRUE,ylim=c(-13,-11.5),xlim=c(1.56,2.2))
contour(denI, levels = seq(2.1, 10, by = 1),
        add = TRUE, col = "black")
