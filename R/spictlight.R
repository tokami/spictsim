

library(TMB)


## surp prod func
predict.b <- function(gamma,m,K,n,B,F,sdb,dt,type = "lamperti"){
    if(type == "lamperti"){
        exp(log(B) + (gamma * m/K - gamma * m/K * (B/K)^(n-1) -
                      F - 0.5*sdb^2  ) *dt) ## + rnorm(1, 0, sqrt(dt) * sdb)) * dt)
    }else if(type == "naive"){
        max(B + (gamma * m/K * B - gamma * m * (B/K)^n - F * B) * dt, 1e-3)
    }
}
## parameters
K <- 4000
m <- 200
n <- 2
gamma <- n^(n/(n-1))/(n-1)
q <- 0.05
B1K <- 1
sdb <- 0.01
sdi <- 0.1
sdf <- 0.1
sdc <- 0.1
fstart <- 0.001
fmax <- 0.3
dt <- 1/16
nyears <- 140
C <- rep(NA,nyears*1/dt)
B <- rep(NA,nyears*1/dt)
B[1] <- B1K * K
C[1] <- B[1] * F[1] * dt
##
F <- exp(log(c(seq(fstart,fmax,length.out=50*1/dt),
               rep(fmax,20*1/dt),
               seq(fmax,fstart,length.out=50*1/dt),
               rep(fstart, 20* 1/dt))) + rnorm(nyears,0,sqrt(dt) * sdf))
for(i in 2:(nyears*1/dt)){
    B[i] <- exp(log(predict.b(gamma,m,K,n,B[i-1],F[i-1],sdb,dt,"lamperti")) + rnorm(1, 0, sqrt(dt) * sdb))
    C[i] <- B[i] * F[i] * dt
}
##
Cobs <- rep(0,nyears)
Iobs <- rep(NA,nyears)
newYears <- seq(1,nyears*1/dt, 1/dt)
years <- rep(1:nyears, each = 1/dt)
for(i in 1:(nyears*1/dt)){
    if(i %in% newYears)
        Iobs[years[i]] <- exp(log(q) + log(B[i]) + rnorm(1, 0, sdi))
    Cobs[years[i]] <- Cobs[years[i]] + ifelse(is.na(C[i]),0,C[i])
    if(i %in% newYears-1) Cobs[years[i]] <- exp(log(Cobs[years[i]]) +
                                                rnorm(1, 0, sdc))
}
logobsI <- log(Iobs)
logobsC <- log(Cobs)


## compile
compile("../src/spictlight.cpp")
dyn.load(dynlib("../src/spictlight"))


##
data <- list(logobsC = logobsC,
             logobsI = logobsI,
             dt = 1/16,
             ii = seq(1, 16*length(logobsI), 16),
             ic = seq(1, 16*length(logobsI), 16),
             nc = rep(16, length(logobsC)),
             priorlogn = c(0,2,0.5))


param <- list()
param$logsdb <- 0
param$logsdi <- 0
param$logsdf <- 0
param$logsdc <- 0
param$logn <- log(2)
param$logm <- log(max(Iobs, na.rm=TRUE)/3)
param$logK <- log(max(Iobs, na.rm=TRUE))
param$logB <- rep(0, length(B))
param$logF <- rep(0, length(F))
param$logq <- log(0.01)


obj <- MakeADFun(data, param, DLL="spictlight", random=c("logB","logF"))
##                 map = list(logn = factor(NA)))

opt <- nlminb(obj$par, obj$fn, obj$grp)
sdr <- sdreport(obj)
summary(sdr)


pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")
vals <- names(sdr$value)
Cpred <- sdr$value[vals == "Cpredsub"]


## Plot
par(mfrow=c(2,2))
## biomass
plot(B, ty='l', lwd=2, ylim=c(0,1.2)*range(B),
     xlab = "Time")
lines(exp(pl$logB), col="darkorange", lwd=2)
mtext("Biomass",font=2,line=0.5)
## fish mort
plot(F, ty='l', lwd=2, ylim=c(0,1.2)*range(F),
     xlab = "Time")
lines(exp(pl$logF), col="darkorange", lwd=2)
mtext("Fishing mortality",font=2,line=0.5)
## catch
plot(C, ty='l', lwd=2, ylim=c(0,1.2)*range(C,na.rm=TRUE),
     xlab = "Time")
lines(Cpred, col="darkorange", lwd=2)
mtext("Catch",font=2,line=0.5)
## sur prod
sp <- diff(B) + C[-1]
plot(B[-1], sp, ty='p', pch=1,
     ylim=c(1.2,1.2)*range(sp,na.rm=TRUE),
     xlab = "B", ylab = "SP")
points(exp(pl$logB)[-1], diff(exp(pl$logB)) + Cpred[-1],
       col="darkorange", pch=3)
mtext("Production curve",font=2,line=0.5)




##
library(ramlegacy)

dat <- load_ramlegacy(tables = c("tb.data","tc.data"))

## common species
speciesAll <- colnames(dat$tb.data)
speciesAll <- speciesAll[speciesAll %in% colnames(dat$tc.data)]
nPerSpec <- sapply(speciesAll, function(x)
    length(na.omit(dat$tb.data[,colnames(dat$tb.data) == x])))
species <- speciesAll[which.max(nPerSpec)]
## species <- "ESOLEPCOAST"

species <- speciesAll[head(tail(nPerSpec,11),1)]
## common years
years.tb <- rownames(dat$tb.data)
years.tc <- rownames(dat$tc.data)
years <- years.tb[which(years.tb %in% years.tc)]
##
Cobs <- dat$tc.data[rownames(dat$tc.data) %in% years,
                    colnames(dat$tc.data) == species]
SBobs <- dat$tb.data[rownames(dat$tb.data) %in% years,
                     colnames(dat$tb.data) == species]
ind <- !(is.na(Cobs) | is.na(SBobs))
Cobs <- Cobs[ind]
SBobs <- SBobs[ind]

##
data <- list(logobsC = log(Cobs),
             logobsI = log(SBobs),
             dt = 1/16,
             ii = seq(1, 16*length(Cobs), 16),
             ic = seq(1, 16*length(SBobs), 16),
             nc = rep(16, length(Cobs)),
             priorlogn = c(0,2,0.5))


param <- list()
param$logsdb <- 0
param$logsdi <- 0
param$logsdf <- 0
param$logsdc <- 0
param$logn <- log(2)
param$logm <- log(max(SBobs, na.rm=TRUE)/3)
param$logK <- log(max(SBobs, na.rm=TRUE))
param$logB <- rep(0, length(Cobs) * 1/dt)
param$logF <- rep(0, length(Cobs) * 1/dt)
param$logq <- log(0.01)




obj <- MakeADFun(data, param, DLL="spictlight", random=c("logB","logF"))
##                 map = list(logn = factor(NA)))

opt <- nlminb(obj$par, obj$fn, obj$grp)
sdr <- sdreport(obj)
summary(sdr)


pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")
vals <- names(sdr$value)
Cpred <- sdr$value[vals == "Cpredsub"]


## Plot
par(mfrow=c(2,2))
## biomass
plot(B, ty='l', lwd=2, ylim=c(0,1.2)*range(B),
     xlab = "Time")
lines(exp(pl$logB), col="darkorange", lwd=2)
mtext("Biomass",font=2,line=0.5)
## fish mort
plot(F, ty='l', lwd=2, ylim=c(0,1.2)*range(F),
     xlab = "Time")
lines(exp(pl$logF), col="darkorange", lwd=2)
mtext("Fishing mortality",font=2,line=0.5)
## catch
plot(C, ty='l', lwd=2, ylim=c(0,1.2)*range(C,na.rm=TRUE),
     xlab = "Time")
lines(Cpred, col="darkorange", lwd=2)
mtext("Catch",font=2,line=0.5)
## sur prod
sp <- diff(B) + C[-1]
plot(B[-1], sp, ty='p', pch=1,
     ylim=c(1.2,1.2)*range(sp,na.rm=TRUE),
     xlab = "B", ylab = "SP")
points(exp(pl$logB)[-1], diff(exp(pl$logB)) + Cpred[-1],
       col="darkorange", pch=3)
mtext("Production curve",font=2,line=0.5)
