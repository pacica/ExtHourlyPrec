rm(list=ls())

# libraries
library(parallel)

setwd("~/Downloads/ExtHourlyPrec-main")

initvalue <- "nogiven" # if initvalue=="given", it skips the procedure for the computation of starting values

if (initvalue=="given"){
  ####################################################################################################################################
  ################ this skips the procedure for the computation of starting values for the estimation of the FESD-GPD model
  ####################################################################################################################################
  # previous data----
  load("output/01_WIC_GPD_static.RData")
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_WIC <- tvth[[q]]$dt
  dt_WIC$ID <- "WIC"
  
  load("output/01_KAN_GPD_static.RData")
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_KAN <- tvth[[q]]$dt
  dt_KAN$ID <- "KAN"
  
  load("output/01_STL_GPD_static.RData")
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_STL <- tvth[[q]]$dt
  dt_STL$ID <- "STL"
  
  load("output/01_EVA_GPD_static.RData")
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_EVA <- tvth[[q]]$dt
  dt_EVA$ID <- "EVA"
  
  load("output/01_LOU_GPD_static.RData")
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_LOU <- tvth[[q]]$dt
  dt_LOU$ID <- "LOU"
  
  # create full panel dataset
  panel <- rbind(dt_WIC,dt_KAN,dt_STL,dt_EVA,dt_LOU)
  panel$exc <- panel$HourlyPrecipitation-panel$tvth_mv
  units <- unique(panel$ID)
  nunits <- length(units)
  
  par_ini <- c(0.0704771468,  0.0716293027,  0.1224237459,  0.1665797461,  0.1668300412, # intercept of sigma for each location
               0.0029396486,  0.0008774869, -0.0031856875, -0.0049889719, -0.0050614738, # intercept of xi for each location
               -1.0498221245, -3.2188758249,  0.6190392084,  4.4987990588) # common parameters (A^sigma, A^xi, B^sigma, B^xi)
  
  # load functions
  source("functions/GPD.stat.fit.R")
  source("functions/param.gradrep.R")
  source("functions/GPD.GAS.nllk.fe_rep.R")
  source("functions/VCOV.R")
  source("functions/GPD.scaled.score.R")
  source("functions/GPD.SD.VaR_tvth.R")
  source("functions/GPD.par.dyn.R")
}else{
  ############################################################################################################################
  ################ this computes the starting values for the estimation of the FESD-GPD model
  ############################################################################################################################
  # load function for static GPD estimation
  source("functions/GPD.stat.fit.R")
  
  # previous data----
  load("output/01_WIC_GPD_static.RData")
  set.seed(1)
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_WIC <- tvth[[q]]$dt
  dt_WIC$ID <- "WIC"
  th <- dt_WIC[,which(colnames(dt_WIC) %in% "tvth_mv")]
  long_term_mean_WIC <- stat.GPD(dt_WIC$HourlyPrecipitation-th,0)$par
  
  load("output/01_KAN_GPD_static.RData")
  set.seed(1)
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_KAN <- tvth[[q]]$dt
  dt_KAN$ID <- "KAN"
  th <- dt_KAN[,which(colnames(dt_KAN) %in% "tvth_mv")]
  long_term_mean_KAN <- stat.GPD(dt_KAN$HourlyPrecipitation-th,0)$par
  
  load("output/01_STL_GPD_static.RData")
  set.seed(1)
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_STL <- tvth[[q]]$dt
  dt_STL$ID <- "STL"
  th <- dt_STL[,which(colnames(dt_STL) %in% "tvth_mv")]
  long_term_mean_STL <- stat.GPD(dt_STL$HourlyPrecipitation-th,0)$par
  
  load("output/01_EVA_GPD_static.RData")
  set.seed(1)
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_EVA <- tvth[[q]]$dt
  dt_EVA$ID <- "EVA"
  th <- dt_EVA[,which(colnames(dt_EVA) %in% "tvth_mv")]
  long_term_mean_EVA <- stat.GPD(dt_EVA$HourlyPrecipitation-th,0)$par
  
  load("output/01_LOU_GPD_static.RData")
  set.seed(1)
  q=1 # change to q=3 for threshold at the 97th percentile
  dt_LOU <- tvth[[q]]$dt
  dt_LOU$ID <- "LOU"
  th <- dt_LOU[,which(colnames(dt_LOU) %in% "tvth_mv")]
  long_term_mean_LOU <- stat.GPD(dt_LOU$HourlyPrecipitation-th,0)$par
  
  panel <- rbind(dt_WIC,dt_KAN,dt_STL,dt_EVA,dt_LOU)
  panel$exc <- panel$HourlyPrecipitation-panel$tvth_mv
  units <- unique(panel$ID)
  nunits <- length(units)
  
  # set grid for static common parameters
  guess_A_sigma <- c(0.1,0.2,0.3,0.35,0.4,0.45,0.5)
  guess_A_xi <- c(0.01,0.02,0.03,0.04,0.045,0.05,0.06)
  guess_B_sigma <- c(0.6,0.65,0.7,0.75,0.8,0.85,0.9)
  guess_B_xi <- c(0.96,0.965,0.97,0.975,0.98,0.985,0.989)
  AB <- as.data.frame(expand.grid(guess_A_sigma,guess_A_xi,guess_B_sigma,guess_B_xi))
  names(AB) <- c("A_sigma","A_xi","B_sigma","B_xi")
  
  # load functions
  source("functions/param.gradrep.R")
  source("functions/GPD.GAS.nllk.fe_rep.R")
  source("functions/VCOV.R")
  source("functions/GPD.scaled.score.R")
  source("functions/GPD.SD.VaR_tvth.R")
  source("functions/GPD.par.dyn.R")
  
  # initialization
  guess_omega_WIC <- guess_omega_KAN <- guess_omega_STL <- guess_omega_EVA <- guess_omega_LOU <- matrix(NA,2,nrow(AB))
  initpar <- initpar_1 <- matrix(NA,nrow(AB),length(units)*2+4)
  nloglik <- numeric()
  
  for (l in 1:nrow(AB)){
    guess_omega_WIC[1,l] <- log(long_term_mean_WIC[1])*(1-AB[l,3])
    guess_omega_WIC[2,l] <- log(long_term_mean_WIC[2]/(1-long_term_mean_WIC[2]))*(1-AB[l,4])
    guess_omega_KAN[1,l] <- log(long_term_mean_KAN[1])*(1-AB[l,3])
    guess_omega_KAN[2,l] <- log(long_term_mean_KAN[2]/(1-long_term_mean_KAN[2]))*(1-AB[l,4])
    guess_omega_STL[1,l] <- log(long_term_mean_STL[1])*(1-AB[l,3])
    guess_omega_STL[2,l] <- log(long_term_mean_STL[2]/(1-long_term_mean_STL[2]))*(1-AB[l,4])
    guess_omega_EVA[1,l] <- log(long_term_mean_EVA[1])*(1-AB[l,3])
    guess_omega_EVA[2,l] <- log(long_term_mean_EVA[2]/(1-long_term_mean_EVA[2]))*(1-AB[l,4])
    guess_omega_LOU[1,l] <- log(long_term_mean_LOU[1])*(1-AB[l,3])
    guess_omega_LOU[2,l] <- log(long_term_mean_LOU[2]/(1-long_term_mean_LOU[2]))*(1-AB[l,4])
    
    
    initpar[l,] <- c(guess_omega_WIC[1,l],guess_omega_KAN[1,l],guess_omega_STL[1,l],guess_omega_EVA[1,l],guess_omega_LOU[1,l],
                     guess_omega_WIC[2,l],guess_omega_KAN[2,l],guess_omega_STL[2,l],guess_omega_EVA[2,l],guess_omega_LOU[2,l],
                     AB[l,1],AB[l,2],AB[l,3],AB[l,4])
    
    initpar_1[l,] <- c(initpar[l,1:(nunits*2)],log(initpar[l,nunits*2+1]),log(initpar[l,nunits*2+2]),log(initpar[l,nunits*2+3]/(1-initpar[l,nunits*2+3])),log(initpar[l,nunits*2+4]/(1-initpar[l,nunits*2+4])))
    
    nloglik[l] <- GPD.GAS.nllk.fe(initpar_1[l,],panel,"heavy",units)
  }
  par_ini <- initpar_1[which(nloglik==min(nloglik,na.rm=T)),]
}

# SD GPD estimation----
tails_type <- "heavy"
GPD_SD_par <- 6
MLE_par_GPD_fe <- numeric(GPD_SD_par) # n.SD.par
llk_GPD_fe <- numeric()
dyn_par_GPD_fe <- quant_GPD_SD_fe <- vector(mode="list")

# SD estimation
print("SD estimation")
MLE_estimates <- optim(par = par_ini,
                       fn = GPD.GAS.nllk.fe,
                       control=list(maxit=5000000),
                       # Custom Inputs
                       x = panel,
                       tails = tails_type,
                       groups = units
)
cat("convergence==",MLE_estimates$convergence)
print("--")

MLE_par_opt_fe <- MLE_estimates$par
MLE_par_GPD_fe <- c(MLE_par_opt_fe[1:(nunits*2)],exp(MLE_par_opt_fe[nunits*2+1]),exp(MLE_par_opt_fe[nunits*2+2]),1/(1+exp(-MLE_par_opt_fe[nunits*2+3])),1/(1+exp(-MLE_par_opt_fe[nunits*2+4])))
llk_GPD_fe <- -GPD.GAS.nllk.fe(MLE_estimates$par, panel, tails_type, units)

print("quantile estimation")
qlevel <- 0.99
for (i in 1:length(units)){
  # parameter dynamics and quantiles
  dyn_par_GPD_fe <- append(dyn_par_GPD_fe,list(par.dyn.GPD(MLE_par_GPD_fe[c(i,i+length(units),(length(units)*2+1):length(MLE_par_GPD_fe))], (panel$HourlyPrecipitation-panel$tvth_mv)[panel$ID==units[i]], tails_type)))
  panel_1<-panel[panel$ID==units[i],]
  quant_GPD_SD_fe <- append(quant_GPD_SD_fe,list(GPD.SD.VaR(dyn_par_GPD_fe[[i]][,1],dyn_par_GPD_fe[[i]][,2],panel_1$tvth_mv,qlevel,panel_1$HourlyPrecipitation)))
}
names(dyn_par_GPD_fe) <- names(quant_GPD_SD_fe) <- units

# compute SE
print("SE estimation")
nllk <- GPD.GAS.nllk.fe
MLE_par <- MLE_par_opt_fe
obs <- panel
nobs <- length(panel$HourlyPrecipitation[panel$HourlyPrecipitation>panel$tvth_mv])

VCOV <- vcovdiff.fe(nllk,MLE_par,obs,"heavy",units)
SE <- sqrt(diag(VCOV$SNDW)/nobs)
SE_rep <- sqrt(diag(VCOV$SNDW_rep)/nobs)

save.image(paste0("output/02_FE_GPD_dynamic_q",q,"_rep.RData"))

