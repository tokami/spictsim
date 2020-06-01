
library(TMB)
## remotes::install_github("ropensci/ramlegacy")
library(ramlegacy)

## surp prod func
sp <- function(gamma,y,SB,SB0,phi){
    if(phi != 1){
        gamma * y * SB - gamma * y * SB0 * (SB/SB0)^phi
    }else{
        exp(1) * y * SB * log(SB/SB0)
    }
}
## parameters
SB0 <- 4000
y <- 0.3
phi <- 2
gamma <- phi^(phi/(phi-1))/(phi-1)
SB1SB0 <- 0.9
sigma_p <- 0.05
sigma_m <- 0.3

Cobs <- rep(NA,100)
SB <- rep(NA,100)
SB[1] <- SB1SB0 * SB0
SBobs <- rep(NA,100)
SBobs[1] <- exp(log(SB[1]) + rnorm(1, 0, sigma_p))
for(i in 2:100){
    Cobs[i-1] <- SB[i-1] * c(seq(0.1,0.9,length.out=50),rep(0.9,25),
                             seq(0.9,0.3,length.out=25))[i-1]   ## harvest rate
    SB[i] <- exp(log(SB[i-1] + sp(gamma,y,SB[i-1], SB0, phi) - Cobs[i-1]) + rnorm(1, 0, sigma_p))
    SBobs[i] <- exp(log(SB[i]) + rnorm(1, 0, sigma_m))
}

plot(SB)

SBobs <- SBobs[-100]
SB <- SB[-100]
Cobs <- Cobs[-100]



## compile
compile("../src/thorson.cpp")
dyn.load(dynlib("../src/thorson"))


## example data
## download_ramlegacy()
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
data <- list(Cobs = Cobs,
             SBobs = SBobs)
SPobs <- diff(SBobs) + Cobs[-length(Cobs)]
plot(SBobs[-1], SPobs)


param <- list()
param$logPhi <- log(2)
param$logSB0 <- log(max(SBobs, na.rm=TRUE))
## param$logSB1SB0 <- log(0.9) ## log(SBobs[1]/exp(param$logSB0))
param$logY <- log(0.3) ## log(mean(SBobs, na.rm=TRUE)/ exp(param$logSB0))
param$logSB <- log(SBobs)
param$logSigma_p <- log(0.1)
param$logSigma_m <- log(0.1)

obj <- MakeADFun(data, param, DLL="thorson", random="logSB")
##                 map = list(logPhi = factor(NA)))

opt <- nlminb(obj$par, obj$fn, obj$grp)
sdr <- sdreport(obj)
summary(sdr)


pl <- as.list(sdr, "Est")
plsd <- as.list(sdr, "Std")


vals <- names(sdr$value)
sp <- sdr$value[vals == "S"]
sb <- exp(sdr$par.random)
plot(sb, ty='l')
points(SBobs)

plot(sb, sp, ty="p",ylim = c(-3000,8000))
points(SBobs[-1], SPobs,col=4)

plot(sp)

exp(pl$logepsilon)
exp(pl$logtau)

obj$fn(opt$par)


obj$par
opt
