##  Loading requried R package
library(mrgsolve)    ## R-package for Loading mrgsolve code into r via mcode from the 'mrgsolve' pckage
library(magrittr)    ## R-package for the pipe, %>% , comes from the magrittr package by Stefan Milton Bache
library(dplyr)       ## R-package for transform and summarize tabular data with rows and columns; %>%
library(tidyverse)   ## R-package for transform and summarize tabular data with rows and columns; %>%
library(ggplot2)     ## R-package for GGplot
library(gridExtra)   ## R-package for grid
library(lattice)     ## for plotting the figure
library(FME)         ## R-package for MCMC simulation and model fitting
library(minpack.lm)  ## R-package for model fitting
library(reshape)     ## Package for melt function to reshape the table
library(truncnorm)   ## Package for the truncated normal distribution function   
library(EnvStats)    ## Package for Environmental Statistics, Including US EPA Guidance
library(invgamma)    ## Package for inverse gamma distribution function
library(foreach)     ## Package for parallel computing
library(doParallel)  ## Package for parallel computing
library(bayesplot)   ## Package for MCMC traceplot
library(patchwork)   ## Package for merge images
library(readxl)      ## Package for load excel
library(ggExtra)     ## R-package for extend the ggplot2

## Input mrgsolve-based PBPK Model
source (file = "mice_mod.R")

## Build mrgsolve-based PBPK Model
Gmod_R    <- mcode ("GMicePBPK.code", GMicePBPK.code) # bulid the gestational model

## Loading datasets for model calibration
Data_0    <- read.csv(file = "DATA_Mice.csv")

# Model calibration for gestational PBPK model based on the data of TK study                                                                                        #
# Dose regimen: 0.08 mg/kg, exposure in GD13
# Abreviation: maternal plasma (MP), maternal liver (ML), fetal brain (FB)
# A1. : Pregnant C57BL/6J mice oral dose to 0.08 mg/kg, matrix: MP, Sampling time: GD13.08, 13.17, 13.33 
# 13.50, 14.00, 14.50, 15.00, 16.00, 17.00
# A2. : Pregnant C57BL/6J mice oral dose to 0.08 mg/kg, matrix: ML, Sampling time: GD13.08, 13.17, 13.33
# 13.50, 14.00, 15.00, 17.00                                          
# A3. : Pregnant C57BL/6J mice oral dose to 0.08 mg/kg, matrix: FB, Sampling time: GD13.08, 13.17, 13.33
# 13.50, 14.00, 15.00, 17.00                               
# B1. : Pregnant C57BL/6J mice IV dose to 0.08 mg/kg, matrix: MP, Sampling time: GD13.02, 13.08, 13.17, 13.33 
# 13.50, 14.00, 14.50, 15.00, 16.00, 17.00
# B2. : Pregnant C57BL/6J mice IV dose to 0.08 mg/kg, matrix: ML, Sampling time: GD13.08, 13.17, 13.33
# 13.50, 14.00, 15.00, 17.00                                          
# B3. : Pregnant C57BL/6J mice IV dose to 0.08 mg/kg, matrix: FB, Sampling time: GD13.08, 13.17, 13.33
# 13.50, 14.00, 15.00, 17.00                              
#===================================================================================================

## Read these datasets and later used in model calibration

OBS.A1  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.A2  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.A3  <- Data_0 %>% filter(Study == 1 & Sample == "FB" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CB_Fet= "Conc")
OBS.B1  <- Data_0 %>% filter(Study == 2 & Sample == "MP" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CPlas = "Conc")
OBS.B2  <- Data_0 %>% filter(Study == 2 & Sample == "ML" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CL    = "Conc")
OBS.B3  <- Data_0 %>% filter(Study == 2 & Sample == "FB" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CB_Fet= "Conc")

## Define the prediction function based on study design of literature of TK study
## Exposure scenario: Dosing in GD13

pred.A <- function(Gpars) { ## Gpars: input parameters
  
  ## Get out of log domain
  Gpars <- exp(Gpars)                   ## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 0.025                  ## Body weight during gestation (measrument data if available); Default value adopted from Chou and Lin (2019)
  tinterval    = 144                    ## Time interval; 
  GTDOSE       = 1                      ## Total dosing/Dose times
  DOSE         = 0.08                   ## Input oral dose  
  GDOSEoral    = DOSE*GBW               ## Amount of oral dose
  
  ## To create a data set of 1 subject receiving DOSEoral for 1 total doses
  Gex.oral.A<- ev (ID = 1,               ## One individual
                   time = 24*13,         ## Dossed start time (GD13)
                   amt  = GDOSEoral,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = GTDOSE - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  Gex.oral.B<- ev (ID = 1,               ## One individual
                   time = 24*13,         ## Dossed start time (GD13)
                   amt  = GDOSEoral,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = GTDOSE - 1,    ## Addtional doseing 
                   cmt  = "APlas_free",  ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  
  Gtsamp  = tgrid(24*13, tinterval*(GTDOSE - 1) + 24*19, 1) ## set up the output time; dtout = 1 hours; start = 0h, end = 24*50h
  
  ## Simulation of exposure scenaior (oral dose to 80 μg/kg)
  Gout.A <- 
    Gmod_R %>% ## Gestational PBPK model
    param (Gpars) %>% ## Update the parameter list with Gpars
    update(atol = 1E-6, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
    mrgsim_d (data = Gex.oral.A, tgrid = Gtsamp)   
  
  ## Extract the "Time", "CPlas", "CL" ,"CB_Fet" from Gout; 
  Goutdf.A = cbind.data.frame(Time   = Gout.A$time, 
                              CPlas  = Gout.A$Plasma,      # CPlas: F-53B concentration in maternal plasma
                              CL     = Gout.A$Liver,       # CL: F-53B concentration in maternal liver
                              CB_Fet = Gout.A$Brain_Fet,   # CB_Fet: F-53B concentration in fetal brain 
                              Mbal   = Gout.A$Bal,         # Mbal: maternal mass balance
                              MbalF  = Gout.A$Bal_Fet)     # MbalF: fetal mass balance
  
  Gout.B <- 
    Gmod_R %>% ## Gestational PBPK model
    param (Gpars) %>% ## Update the parameter list with Gpars
    update(atol = 1E-6, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
    mrgsim_d (data = Gex.oral.B, tgrid = Gtsamp)   
  
  ## Extract the "Time", "CPlas", "CL" , "CPla" ,"CAM" , "CB_Fet" from Gout; 
  Goutdf.B = cbind.data.frame(Time   = Gout.B$time, 
                              CPlas  = Gout.B$Plasma,      # CPlas: F-53B concentration in maternal plasma
                              CL     = Gout.B$Liver,       # CL: F-53B concentration in maternal liver
                              CB_Fet = Gout.B$Brain_Fet)   # CB_Fet: F-53B concentration in fetal brain
  
  
  Goutdf.A <- Goutdf.A %>% filter (Time > 0) # filter the value at time = 0
  Goutdf.B <- Goutdf.B %>% filter (Time > 0)                                 
  
  return (list("Goutdf.A"=Goutdf.A,
               "Goutdf.B"=Goutdf.B)) # Return Goutdf
}

result <- pred.A(Gtheta.int)

plot(result$Goutdf.A$Time,result$Goutdf.A$CPlas,type="l",lwd=2,xlab="Time",ylab="F-53B concentration in plasma")
plot(result$Goutdf.B$Time,result$Goutdf.B$CPlas,type="l",lwd=2,xlab="Time",ylab="F-53B concentration in plasma")
plot(result$Goutdf.A$Time,result$Goutdf.A$Mbal,type="l",lwd=2,xlab="Time",ylab="Mass")
plot(result$Goutdf.A$Time,result$Goutdf.A$MbalF,type="l",lwd=2,xlab="Time",ylab="Mass_Fetal")

## Create a cost fuction and later used in model calibration 
## Estimate the model residual with experimental data by modCost function (from FME package) ,err = "SD"
Gcost<-function (pars){
  
  outdf <- pred.A (Gpars = pars)
  
  cost<- modCost  (model = outdf$Goutdf.A, obs = OBS.A1, x ="Time" ,weight = "mean")
  cost<- modCost  (model = outdf$Goutdf.A, obs = OBS.A2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$Goutdf.A, obs = OBS.A3, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$Goutdf.B, obs = OBS.B1, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$Goutdf.B, obs = OBS.B2, x ="Time" ,weight = "mean",cost = cost)
  cost<- modCost  (model = outdf$Goutdf.B, obs = OBS.B3, x ="Time" ,weight = "mean",cost = cost)
  
  return(cost)
}

## Local sensitivity analysis
## Choose the senstiive parameters in the model
## initial parmaeters
Gtheta.int <- log(c(
  Vmax_baso_invitro              = 393.45,                      
  Km_baso                        = 27.2,                        
  Vmax_apical_invitro            = 4185,                        
  Km_apical                      = 52.3,                        
  RAFbaso                        = 3.99,                        
  RAFapi                         = 2.81,                       
  Kdif                           = 4.6e-5,                       
  KeffluxC                       = 5.6,                        
  KbileC                         = 0.00039,                       
  KurineC                        = 0.015,                        
  Free                           = 0.009,                        
  PL                             = 2.10,                        
  PK                             = 0.80,                         
  PM                             = 0.16,
  PF                             = 0.13,
  PRest                          = 0.22,                         
  PPla                           = 0.59,
  K0C                            = 1,                           
  Kabsc                          = 2.43,                        
  KunabsC                        = 5.40e-4,                     
  Ktrans1C                       = 1.27,
  Ktrans2C                       = 1,
  Ktrans3C                       = 0.23,
  Ktrans4C                       = 0.001,
  Free_Fet                       = 0.009,
  PL_Fet                         = 2.10,
  PB_Fet                         = 1.55,
  PRest_Fet                      = 0.22
))

## Senstivity function (FME) 
## Check the senstive parameters in the model
SnsPlasma <- sensFun(func = Gcost, parms = Gtheta.int, varscale = 1)
Sen       <- summary(SnsPlasma)
plot(Sen)

sen1 <- Sen %>% filter(Mean != 0)
SnsPlasma1 <- SnsPlasma[,row.names(sen1)]

## Selected senstive parameters
Gtheta <- Gtheta.int[abs(Sen$Mean) > 1.2*mean(abs(Sen$Mean))]
Gtheta 

## Selected parameters
Gtheta.int <- log(c(
  #Vmax_baso_invitro              = 393.45,                      
  Km_baso                        = 27.2,                        
  Vmax_apical_invitro            = 4185,                        
  Km_apical                      = 52.3,                        
  #RAFbaso                        = 3.99,                        
  #RAFapi                         = 2.81,                       
  Kdif                           = 4.6e-5,                       
  #KeffluxC                       = 5.6,                        
  KbileC                         = 3.9e-4,                       
  KurineC                        = 0.015,                        
  Free                           = 0.009,                        
  #PL                             = 2.10,                        
  #PK                             = 0.80,                         
  #PM                             = 0.16,
  PF                             = 0.13,
  PRest                          = 0.22,                         
  PPla                           = 0.59,
  K0C                            = 1,                           
  #Kabsc                          = 2.43,                        
  #KunabsC                        = 5.40e-4,                     
  Ktrans1C                       = 1.27,
  Ktrans2C                       = 1,
  #Ktrans3C                       = 0.23,
  #Ktrans4C                       = 0.001,
  Free_Fet                       = 0.009
  #PL_Fet                         = 2.10
  #PB_Fet                         = 1.55
  #PRest_Fet                      = 0.22
))
## PBPK model calibration
## Least squres fit using method "Nelder-Mead" algorithm

lower <- c(-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-Inf,-2.97,-Inf,-Inf,-Inf,-Inf)
upper <- c(Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf)

GFit <- modFit(f = Gcost, p = Gtheta.int, method = "Nelder-Mead",
               control = nls.lm.control(nprint = 1),upper = upper,lower = lower)

## Summary of fitting results
summary(GFit)
exp(GFit$par)  ## Get the arithmetic value out of the log domain
 
Gcost(GFit$par)

Sim.fit.A = pred.A(GFit$par)
df.sim.A1 <- cbind.data.frame(Time = Sim.fit.A$Goutdf.A$Time,
                              CPlas = Sim.fit.A$Goutdf.A$CPlas,
                              CL = Sim.fit.A$Goutdf.A$CL,
                              CB_Fet = Sim.fit.A$Goutdf.A$CB_Fet)

df.sim.B1 <- cbind.data.frame(Time = Sim.fit.A$Goutdf.B$Time,
                              CPlas = Sim.fit.A$Goutdf.B$CPlas,
                              CL = Sim.fit.A$Goutdf.B$CL,
                              CB_Fet = Sim.fit.A$Goutdf.B$CB_Fet)

OBS.A1.1  <- Data_0 %>% filter(Study == 1 & Sample == "MP" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.A2.1  <- Data_0 %>% filter(Study == 1 & Sample == "ML" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = "SD")
OBS.A3.1  <- Data_0 %>% filter(Study == 1 & Sample == "FB" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CB_Fet= "Conc", SD = "SD")
OBS.B1.1  <- Data_0 %>% filter(Study == 2 & Sample == "MP" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CPlas = "Conc", SD = "SD")
OBS.B2.1  <- Data_0 %>% filter(Study == 2 & Sample == "ML" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CL    = "Conc", SD = "SD")
OBS.B3.1  <- Data_0 %>% filter(Study == 2 & Sample == "FB" & Dose == 0.08 & Time != 0) %>% select(Time = "Time", CB_Fet= "Conc", SD = "SD")

plot.A1 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CPlas), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A1, aes(Time, CPlas), col = "#e3716e", size = 2.5) +
  geom_errorbar(data = OBS.A1.1, aes(x = Time, ymin = CPlas - SD, ymax = CPlas + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/mL)") + 
  xlab("Time (h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 0.301), expand = c(0, 0))+
  theme_classic()+
  annotate("text", x = 380, y = 0.29, label = "Plasma", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.3), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.3, xend = 462, yend = 0.3), 
               color = "black", size = 0.5)

plot.A2 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CL), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A2, aes(Time, CL), col = "#e3716e", size = 2.5) + 
  geom_errorbar(data = OBS.A2.1, aes(x = Time, ymin = CL - SD, ymax = CL + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/g)") + 
  xlab("Time (h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.501), expand = c(0, 0)) + 
  theme_classic()+
  annotate("text", x = 380, y = 0.48, label = "Liver", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.5), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.5, xend = 462, yend = 0.5), 
               color = "black", size = 0.5)

plot.A3 <- ggplot() +
  geom_line(data = df.sim.A1, aes(Time, CB_Fet), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.A3, aes(Time, CB_Fet), col = "#e3716e", size = 2.5) + 
  geom_errorbar(data = OBS.A3.1, aes(x = Time, ymin = CB_Fet - SD, ymax = CB_Fet + SD), col = "#e3716e", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/g)") + 
  xlab("Time (h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  theme_classic()+
  annotate("text", x = 380, y = 0.48, label = "Fetal brain", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.5), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.5, xend = 462, yend = 0.5), 
               color = "black", size = 0.5)

plot.B1 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CPlas), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B1, aes(Time, CPlas), col = "#7ac7e2", size = 2.5, shape = 17, fill = "#7ac7e2") + 
  geom_errorbar(data = OBS.B1.1, aes(x = Time, ymin = CPlas - SD, ymax = CPlas + SD), col = "#7ac7e2", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/mL)") + 
  xlab("Time (h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.0, 0.601), expand = c(0, 0))+
  theme_classic()+
  annotate("text", x = 380, y = 0.58, label = "Plasma", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.6), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.6, xend = 462, yend = 0.6), 
               color = "black", size = 0.5)


plot.B2 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CL), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B2, aes(Time, CL), col = "#7ac7e2", size = 2.5, shape = 17, fill = "#7ac7e2") + 
  geom_errorbar(data = OBS.B2.1, aes(x = Time, ymin = CL - SD, ymax = CL + SD), col = "#7ac7e2", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/g)") + 
  xlab("Time(h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  theme_classic()+
  annotate("text", x = 380, y = 0.48, label = "Liver", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.5), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.5, xend = 462, yend = 0.5), 
               color = "black", size = 0.5)

plot.B3 <- ggplot() +
  geom_line(data = df.sim.B1, aes(Time, CB_Fet), col = "#5A5A5A", lwd = 0.8) +
  geom_point(data = OBS.B3, aes(Time, CB_Fet), col = "#7ac7e2", size = 2.5, shape = 17, fill = "#7ac7e2") + 
  geom_errorbar(data = OBS.B3.1, aes(x = Time, ymin = CB_Fet - SD, ymax = CB_Fet + SD), col = "#7ac7e2", width = 2, size = 0.8) +
  ylab("Concentration of F-53B (μg/g)") +
  xlab("Time(h)") +
  scale_x_continuous(breaks = seq(312, 462, 30), limits = c(312, 462), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  theme_classic()+
  annotate("text", x = 380, y = 0.48, label = "Fetal brain", size = 6, color = "black")+
  geom_segment(aes(x = 462, y = 0.0, xend = 462, yend = 0.5), 
               color = "black", size = 0.5) +
  geom_segment(aes(x = 312, y = 0.5, xend = 462, yend = 0.5), 
               color = "black", size = 0.5)

grid.arrange(plot.A1, plot.A2, plot.A3, plot.B1, plot.B2, plot.B3, ncol = 3)

## Save the fitting results to RDS files
saveRDS(GFit, file = "GFit_R.rds") 

GFit_R <- readRDS(file = "GFit_R.rds")

RatCost   <- Gcost(GFit_R$par)

GPDat <- cbind.data.frame (name= RatCost$residuals$name,
                           OBS = RatCost$residuals$obs, 
                           PRE = RatCost$residuals$mod)

## Transformed the predicted and obseved values using log10-sacle to do the plot
GPDat %<>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Rat")

## Estimating the R-squared and goodness-of-fit using linear regression model
fit <- lm(Log.OBS ~Log.PRE, data = GPDat)
summary(fit)

R2 <- summary(fit)$adj.r.squared

GPDat %<>% mutate(res = residuals(fit), 
                  prediction = predict(fit), 
                  OPR = PRE/OBS,
                  logOPR =  log(OPR,10)) ## OPR: the ratio of prediction value and observed data

p <- ggplot(GPDat, aes(Log.OBS, Log.PRE)) + 
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(3, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.33, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(2, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_abline(intercept = log(0.5, 10), slope = 1, color = "grey", linetype = "dashed", linewidth = 1, alpha = 0.8) +
  geom_point(colour = "steelblue4", size = 4) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2, 1), labels = scales::math_format(10^.x)) +
  scale_x_continuous(limits = c(-2, 1), labels = scales::math_format(10^.x))

## Set up your theme and font
windowsFonts("Times" = windowsFont("Times New Roman"))

p1 <- p + 
  theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position         ='none') +
    labs (x = "Observed values (μg/mL)",  y = "Predicted values (μg/mL)") +
    annotate("text", x = Inf, y = -1.8, label = paste0("R² = ", round(R2, 2)), size = 5, fontface = "bold", hjust = 1.1, vjust = 1.1)

p2 <-
  ggplot(GPDat, aes(Log.PRE, logOPR)) +
  geom_hline(yintercept = log10(2),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(0.5),linetype = 3, color   = "black", size =1) +
  geom_hline(yintercept = log10(3),linetype = 3, color   = "grey", size =1) +
  geom_hline(yintercept = log10(0.33),linetype = 3, color   = "grey", size =1) +
  geom_point(color   = "steelblue4", 
             aes(shape= as.factor(Species)),size = 3) +
  annotation_logticks() +
  scale_y_continuous(limits = c(-2,1), labels = scales::math_format(10^.x))+
  scale_x_continuous(limits = c(-2,1),labels = scales::math_format(10^.x))

p2 <- p2 +
  theme (
    plot.background         = element_rect (fill ="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border            = element_rect (colour = "black", fill=NA, linewidth =2),
    panel.background        = element_rect (fill ="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(), 
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes 
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
    legend.position='none') +
  labs (x = "Predicted values (μg/mL)" ,y = "Predicted/Observed")

combined_plot <- p1 + p2 
combined_plot

GPDat$ID <- 1:nrow(GPDat)
p3 <- ggplot(GPDat, aes(ID, logOPR)) + geom_point() + theme_bw()+
  scale_y_continuous(limits = c(-2, 1), labels = scales::math_format(10^.x)) +
  annotation_logticks()

p3 <- ggMarginal(p3, type = "histogram", margins = "y", size = 4,
                 color="black", alpha = 0.4, bins = 20,
                 fill = "#00AFBB", position="identity")
print(p3)
### Sensitivity analysis
mice.theta.G <- log(c(
  # Physiological parameters
  BW                             = 0.025,
  Htc                            = 0.48,
  QLC                            = 0.161,
  QKC                            = 0.091,
  QMC                            = 0.002,
  QFC                            = 0.07,
  VLC                            = 0.055,
  VKC                            = 0.017,
  VMC                            = 0.01,
  VFC                            = 0.068,
  VPlasC                         = 0.049,
  VFilC                          = 0.0017,
  VPTCC                          = 1.35e-4,
  protein                        = 2.0e-6,
  GFRC                           = 59,
  GEC                            = 0.54,
  QLC_Fet                        = 0.161,
  VPlasC_Fet                     = 0.049,
  QBC_Fet                        = 0.1055,
  # Chemical-specific parameters (final mean values)
  Vmax_baso_invitro              = 393.45,     
  Km_baso                        = 16.58,    ### fitting parameters
  Vmax_apical_invitro            = 4168.90,  ### fitting parameters
  Km_apical                      = 72.85,   ### fitting parameters                  
  RAFbaso                        = 3.99,                       
  RAFapi                         = 2.81,    
  KeffluxC                       = 5.60,   
  KbileC                         = 9.18e-6,  ### fitting parameters
  KurineC                        = 3.33e-2,  ### fitting parameters                     
  Free                           = 0.067,    ### fitting parameters
  PL                             = 2.1,     
  PK                             = 0.8,
  PM                             = 0.16,
  PF                             = 0.27,     ### fitting parameters
  PRest                          = 0.43,     ### fitting parameters
  PPla                           = 0.07,     ### fitting parameters
  K0C                            = 5.64,     ### fitting parameters
  Kabsc                          = 2.43,
  Kdif                           = 9.94e-6,  ### fitting parameters                      
  KunabsC                        = 5.4e-4,    
  Ktrans1C                       = 1.42,     ### fitting parameters
  Ktrans2C                       = 0.59,     ### fitting parameters
  Ktrans3C                       = 0.23, 
  Ktrans4C                       = 0.001,
  Free_Fet                       = 0.070,    ### fitting parameters
  PL_Fet                         = 2.10,    
  PB_Fet                         = 1.55,    
  PRest_Fet                      = 0.22,  
  # Growth equatons
  SA_VPlas            = 1,  #: Unitless, SA for maternal Volume of plasma (VPlas)
  SA_VAmX             = 1,  #: Unitless, SA for maternal volume of amniotic fluid volume (VAmX)
  SA_VPla             = 1,  #: Unitless, SA for maternal volume of placenta (VPla)
  SA_VF_P             = 1,  #: Unitless, SA for maternal volume of fat during pregnant (VF_P)
  SA_VM_P             = 1,  #: Unitless, SA for maternal volume of mammary during pregnant (VM_P)
  SA_QC_P             = 1,  #: Unitless, SA for maternal Cardiac output during pregnancy (QC)
  SA_QF_P             = 1,  #: Unitless, SA for maternal Blood flows of fat during prengnacy (QF_P)
  SA_QM_P             = 1,  #: Unitless, SA for maternal Blood flows of mammary gland during prengnacy (QM_P)
  SA_QPla             = 1,  #: Unitless, SA for maternal Blood flows of placenta (QPla)
  SA_VFet             = 1,  #: Unitless, SA for volume of fetus (VFet)
  SA_VL_Fet           = 1,  #: Unitless, SA for volume of fetal liver (VL_Fet)
  SA_VB_Fet           = 1,  #: Unitless, SA for volume of fetal brain (VB_Fet)
  SA_QC_Fet           = 1   #: Unitless, SA for fetal Cardiac output (QC_Fet)
))

pred.mice <- function(Gpars) {
  
  ## Get out of log domain
  Gpars <- exp(Gpars)                   ## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 0.025                  ## Body weight during gestation (measurement data if available); Default value adopted from Chou and Lin (2019)
  tinterval    = 144                    ## Time interval; 
  GTDOSE       = 1                      ## Total dosing/Dose times
  DOSE         = 0.08                   ## Input oral dose  
  GDOSEoral    = DOSE*GBW               ## Amount of oral dose
  
  # To create a exposure scenario
  Gex.oral <- ev (ID   = 1,             ## One individual
                  amt  = GDOSEoral,     ## Amount of dose 
                  ii   = tinterval,     ## Time interval
                  addl = GTDOSE - 1,    ## Addtional doseing 
                  cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                  replicate = FALSE)    ## No replicate
  
  Gtsamp  = tgrid(24*13, tinterval*(GTDOSE - 1) + 24*19, 1) ## Simulation time 
  
  ## Simulation of exposure scenaior (oral dose to 80 μg/kg)
  Gout <- 
    Gmod_R %>%
    param (Gpars) %>%
    update(atol = 1E-6, maxsteps = 500000) %>%          
    mrgsim_d (data = Gex.oral, tgrid = Gtsamp)
  
  Goutdf = cbind.data.frame(Time   = Gout$time, 
                            CPlas  = Gout$Plasma, 
                            CPlas_Fet = Gout$Plasma_Fet,
                            AUC_CPlas = Gout$AUC_CPlas,
                            AUC_CPlas_Fet = Gout$AUC_CPlas_Fet,
                            AUC_CL = Gout$AUC_CL,
                            AUC_CPla = Gout$AUC_CPla,
                            AUC_CL_Fet=Gout$AUC_CL_Fet,
                            AUC_CB_Fet=Gout$AUC_CB_Fet)
  
  
  Goutdf <- Goutdf %>% filter (Time == 19*24) ## The model output on GD19
  
  
  return (list("G" = Goutdf))
}

result  <-  pred.mice(GFit$par)

par.G <- GFit_R$par
mice.theta.G[names(par.G)]  <- as.numeric(par.G) 

NSC_func <- function (Gpars, Pred) {
  nG <- length(Gpars)
  NSC_GCA     = matrix(NA, nrow = length(Gpars) , ncol = 6)
  
  for (i in 1:nG) {
    Gpars.new      <- Gpars %>% replace(i, log(exp((Gpars[i]))*1.01))
    Mnew.G         <- Pred(Gpars.new)
    M.G            <- Pred(Gpars)
    delta.Gpars    <- exp(Gpars[i])/(exp(Gpars[i])*0.01)
    
    ## Estimated the AUC
    AUC.GPlas.new       =  Mnew.G$G %>% select (AUC_CPlas)
    AUC.GPlas.ori       =  M.G$G    %>% select (AUC_CPlas)
    AUC.GPlas_Fet.new   =  Mnew.G$G %>% select (AUC_CPlas_Fet)
    AUC.GPlas_Fet.ori   =  M.G$G    %>% select (AUC_CPlas_Fet)
    AUC.GL.new          =  Mnew.G$G %>% select (AUC_CL)
    AUC.GL.ori          =  M.G$G    %>% select (AUC_CL)
    AUC.GPla.new        =  Mnew.G$G %>% select (AUC_CPla)
    AUC.GPla.ori        =  M.G$G    %>% select (AUC_CPla)
    AUC.GL_Fet.new      =  Mnew.G$G %>% select (AUC_CL_Fet)
    AUC.GL_Fet.ori      =  M.G$G    %>% select (AUC_CL_Fet)
    AUC.GB_Fet.new      =  Mnew.G$G %>% select (AUC_CB_Fet)
    AUC.GB_Fet.ori      =  M.G$G    %>% select (AUC_CB_Fet)
    
    delta.AUC.GPlas     =  AUC.GPlas.new - AUC.GPlas.ori
    delta.AUC.GPlas_Fet =  AUC.GPlas_Fet.new - AUC.GPlas_Fet.ori
    delta.AUC.GL        =  AUC.GL.new -  AUC.GL.ori
    delta.AUC.GPla      =  AUC.GPla.new - AUC.GPla.ori
    delta.AUC.GL_Fet    =  AUC.GL_Fet.new - AUC.GL_Fet.ori
    delta.AUC.GB_Fet    =  AUC.GB_Fet.new - AUC.GB_Fet.ori
    
    NSC_GCA     [i, 1]   <- as.numeric((delta.AUC.GPlas/AUC.GPlas.ori) * delta.Gpars)
    NSC_GCA     [i, 2]   <- as.numeric((delta.AUC.GPlas_Fet /AUC.GPlas_Fet.ori) * delta.Gpars)
    NSC_GCA     [i, 3]   <- as.numeric((delta.AUC.GL /AUC.GL.ori) * delta.Gpars)
    NSC_GCA     [i, 4]   <- as.numeric((delta.AUC.GPla /AUC.GPla.ori) * delta.Gpars)
    NSC_GCA     [i, 5]   <- as.numeric((delta.AUC.GL_Fet /AUC.GL_Fet.ori) * delta.Gpars)
    NSC_GCA     [i, 6]   <- as.numeric((delta.AUC.GB_Fet /AUC.GB_Fet.ori) * delta.Gpars)
  }
  return (NSC_GCA)
}

# Model results
A <- NSC_func (mice.theta.G, pred.mice)

rownames (A)  = names(mice.theta.G)
rownames(A) <- gsub("SA_", "", rownames(A))
colnames (A)  = c("NSC_CPlas", "NSC_CPlas_Fet", "NSC_CL", "NSC_CPla", "NSC_CL_Fet", "NSC_CB_Fet")

NSC_GCA_M <- data.frame(A)

##NSC_GCA_M <- NSC_GCA_M %>%
#mutate_all(~replace(., . == 0, "<1e-5"))
#write.csv(NSC_GCA_M, file = 'Table_3_M_Full.csv')

##################################### Circle barplot function ###############################################
## plot modifed from "R graph gallery: https://www.r-graph-gallery.com/297-circular-barplot-with-groups/ "  #
#############################################################################################################
melt.maternal.Plas        = melt(NSC_GCA_M[,1]) 
melt.maternal.Plas$group  = c("Plasma") 
melt.maternal.Liver       = melt(NSC_GCA_M[,3])
melt.maternal.Liver$group = c("Liver")
melt.maternal.Placenta    = melt(NSC_GCA_M[,4])
melt.maternal.Placenta$group = c("Placenta")

melt.data.maternal         = rbind (melt.maternal.Plas,melt.maternal.Liver,melt.maternal.Placenta)
melt.data.maternal$par     = rep(rownames(NSC_GCA_M),3) 
data.maternal              = melt.data.maternal%>%filter(abs(value)>=0.3)
data.maternal$group <- factor(data.maternal$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data.maternal$group), ncol(data.maternal)))
colnames(to_add) <- colnames(data.maternal)
to_add$group <- rep(levels(data.maternal$group), each=empty_bar)
data.maternal <- rbind(data.maternal, to_add)
data.maternal <- data.maternal%>% arrange(group)
data.maternal$id <- seq(1, nrow(data.maternal))

# Get the name and the y position of each label
label_data <- data.maternal
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data.maternal %>% 
  group_by(group) %>% 
  summarize(start=min(id)-0.2, end=max(id) - empty_bar+0.2) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

windowsFonts(Times=windowsFont("Times New Roman"))

# Make the plot
p1 <- ggplot(data.maternal, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 90, xend = start, yend = 90), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
 
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data.maternal$id),3), y = c(30, 60, 90), label = c("30%", "60%", "90%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text= element_text (family = "Times"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id,  y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(0.5,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

melt.fetal.Plas        = melt(NSC_GCA_M[,2]) 
melt.fetal.Plas$group  = c("Plasma") 
melt.fetal.Liver       = melt(NSC_GCA_M[,5])
melt.fetal.Liver$group = c("Liver")
melt.fetal.Brain       = melt(NSC_GCA_M[,6])
melt.fetal.Brain$group = c("Brain")

melt.data.fetal        = rbind (melt.fetal.Plas,melt.fetal.Liver,melt.fetal.Brain)
melt.data.fetal$par = rep(rownames(NSC_GCA_M),3) 
data.fetal             = melt.data.fetal%>%filter(abs(value)>=0.3)
data.fetal$group <- factor(data.fetal$group)

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 4
to_add <- data.frame( matrix(NA, empty_bar*nlevels(data.fetal$group), ncol(data.fetal)))
colnames(to_add) <- colnames(data.fetal)
to_add$group <- rep(levels(data.fetal$group), each=empty_bar)
data.fetal <- rbind(data.fetal, to_add)
data.fetal <- data.fetal%>% arrange(group)
data.fetal$id <- seq(1, nrow(data.fetal))

# Get the name and the y position of each label
label_data <- data.fetal
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

#prepare a data frame for base lines
base_data <- data.fetal %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

# prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

windowsFonts(Times=windowsFont("Times New Roman"))

# Make the plot
p2 <- ggplot(data.fetal, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  
  # Add a val=90/60/30 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 90, xend = start, yend = 90), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 30, xend = start, yend = 30), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data.fetal$id),3), y = c(30, 60, 90), label = c("30%", "60%", "90%")  , color="red", size=4, angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=abs(value*100), fill=group), stat="identity", alpha=0.5) +
  ylim(-100,150) +
  theme_minimal() +
  theme(
    legend.position = "none",
    text= element_text (family = "Times"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id,  y=abs(value*100)+10, label=par, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
  # Add base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -18, label=group), hjust=c(0.5,1,0), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)

combined_plot <- p1 + p2
combined_plot

## Loading dataset for model evaluation
evaluation <- read.csv(file = "Mice Evaluation.csv")

OBS.1_MF <- evaluation %>% filter(Study==1) %>% filter( Sample == "MF" )%>% select(Time = "Time", Study = "Study", CF = "Conc", SD = "SD")
OBS.1_MR <- evaluation %>% filter(Study==1) %>% filter( Sample == "MR" )%>% select(Time = "Time", Study = "Study", CR = "Conc", SD = "SD")
OBS.1_AM <- evaluation %>% filter(Study==1) %>% filter( Sample == "AM" )%>% select(Time = "Time", Study = "Study", CAm = "Conc", SD = "SD")
OBS.1_FL <- evaluation %>% filter(Study==1) %>% filter( Sample == "FL" )%>% select(Time = "Time", Study = "Study", CL_Fet = "Conc", SD = "SD")
OBS.2_MF <- evaluation %>% filter(Study==2) %>% filter( Sample == "MF" )%>% select(Time = "Time", Study = "Study", CF = "Conc", SD = "SD")
OBS.2_MR <- evaluation %>% filter(Study==2) %>% filter( Sample == "MR" )%>% select(Time = "Time", Study = "Study", CR = "Conc", SD = "SD")
OBS.2_AM <- evaluation %>% filter(Study==2) %>% filter( Sample == "AM" )%>% select(Time = "Time", Study = "Study", CAm = "Conc", SD = "SD")
OBS.2_FL <- evaluation %>% filter(Study==2) %>% filter( Sample == "FL" )%>% select(Time = "Time", Study = "Study", CL_Fet = "Conc", SD = "SD")

pred.A <- function(Gpars) { ## Gpars: input parameters
  
  ## Get out of log domain
  Gpars <- exp(Gpars)                   ## Return a list of exp (parametrs for gestational model) from log scale
  
  ## Exposure scenario for gestational exposure
  GBW          = 0.025                  ## Body weight during gestation (measrument data if available); Default value adopted from Chou and Lin (2019)
  tinterval    = 144                    ## Time interval; 
  GTDOSE       = 1                      ## Total dosing/Dose times
  DOSE         = 0.08                   ## Input oral dose  
  GDOSEoral    = DOSE*GBW               ## Amount of oral dose
  
  ## To create a data set of 1 subject receiving DOSEoral for 1 total doses
  Gex.oral.A<- ev (ID = 1,               ## One individual
                   time = 24*13,         ## Dossed start time (GD13)
                   amt  = GDOSEoral,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = GTDOSE - 1,    ## Addtional doseing 
                   cmt  = "AST",         ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  Gex.oral.B<- ev (ID = 1,               ## One individual
                   time = 24*13,         ## Dossed start time (GD13)
                   amt  = GDOSEoral,     ## Amount of dose 
                   ii   = tinterval,     ## Time interval
                   addl = GTDOSE - 1,    ## Addtional doseing 
                   cmt  = "APlas_free",  ## The dosing comaprtment: AST Stomach  
                   replicate = FALSE)    ## No replicate
  
  
  Gtsamp  = tgrid(24*13, tinterval*(GTDOSE - 1) + 24*19, 1) ## set up the output time; dtout = 1 hours; start = 0h, end = 24*50h
  
  ## Simulation of exposure scenaior (oral dose to 0.08 mg/kg)
  Gout.A <- 
    Gmod_R %>% ## Gestational PBPK model
    param (Gpars) %>% ## Update the parameter list with Gpars
    Req(Liver_Fet, Fat, Rest, AmnioticFluid)%>%                        # select model output
    update(atol = 1E-6, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
    mrgsim_d (data = Gex.oral.A, tgrid = Gtsamp)   
  
   Goutdf.A = cbind.data.frame(Time   = Gout.A$time, 
                              CL_Fet = Gout.A$Liver_Fet,   # CL_Fet: F-53B concentration in fetal liver
                              CF     = Gout.A$Fat,         # CF: F-53B concentration in maternal Fat
                              CR     = Gout.A$Rest,        # CR: F-53B concentration in Rest of body
                              CAm    = Gout.A$AmnioticFluid) # CAm: F-53B concentration in Amniotic Fluid)
  
  Gout.B <- 
    Gmod_R %>% ## Gestational PBPK model
    param (Gpars) %>% ## Update the parameter list with Gpars
    Req(Liver_Fet, Fat, Rest,  AmnioticFluid,)%>%   # select model output
    update(atol = 1E-6, maxsteps = 500000) %>%  ## Atol: Absolute tolerance parameter; maxsteps:maximum number of steps; mindt: simulation output time below which there model will assume to have not advanced          
    mrgsim_d (data = Gex.oral.B, tgrid = Gtsamp)   
  
  Goutdf.B = cbind.data.frame(Time   = Gout.B$time, 
                              CL_Fet = Gout.B$Liver_Fet,   # CL_Fet: F-53B concentration in fetal liver
                              CF     = Gout.B$Fat,         # CF: F-53B concentration in maternal Fat
                              CR     = Gout.B$Rest,        # CR: F-53B concentration in Rest of body
                              CAm    = Gout.B$AmnioticFluid)# CAm: F-53B concentration in Amniotic Fluid)
   
  Goutdf.A <- Goutdf.A %>% filter (Time > 0) # filter the value at time = 0
  Goutdf.B <- Goutdf.B %>% filter (Time > 0) 
  
  return (list("Goutdf.A"=Goutdf.A,
               "Goutdf.B"=Goutdf.B)) # Return Goutdf
}

out_1 <- pred.A(GFit_R$par)[[1]] %>% mutate(Study = 1, SD = 0)%>% filter(Time > 0)
out_2 <- pred.A(GFit_R$par)[[2]] %>% mutate(Study = 2, SD = 0)%>% filter(Time > 0)
out_M_E <- rbind.data.frame (out_1, out_2)

out_M_E_408 <- out_M_E %>% filter(Time == 408)

out_R_subset <- out_M_E_408 %>% select(CR, SD, Study) %>% mutate(Matrix = c("Pre.Rest"))

OBS.1_MR_subset <- OBS.1_MR %>% select(CR, SD, Study) %>% mutate(Matrix = c("Obs.Rest"))
OBS.2_MR_subset <- OBS.2_MR %>% select(CR, SD, Study) %>% mutate(Matrix = c("Obs.Rest"))

combined_data <- rbind(out_R_subset, OBS.1_MR_subset,OBS.2_MR_subset)
combined_data$Study <- ifelse(combined_data$Study == 1, "Oral", "IV")

ggplot(combined_data, aes(x = Matrix, y = CR, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge", width = 1)+
  geom_errorbar(aes(ymin=CR, ymax = CR + SD), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Study) +
  scale_fill_manual(values = c("Pre.Rest.Oral" = "grey", "Pre.Rest.IV" = "#808080",
                               "Obs.Rest.Oral" = "#e3716e", "Obs.Rest.IV" = "#7ac7e2")) +
  labs(x = "Rest of body", y = "Concentration of F-53B (μg/g)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(colour = "black"),  
        legend.position = "none")+
  scale_y_continuous(limits = c(0, 0.08), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(-0.5, nrow(combined_data) + 0.5))

out_A_subset <- out_M_E_408 %>% select(CAm, SD, Study) %>% mutate(Matrix = c("Pre.AM"))

OBS.1_AM_subset <- OBS.1_AM %>% filter(Time == 408) %>% select(CAm, SD, Study) %>% mutate(Matrix = c("Obs.AM"))
OBS.2_AM_subset <- OBS.2_AM %>% filter(Time == 408) %>% select(CAm, SD, Study) %>% mutate(Matrix = c("Obs.AM"))

combined_data <- rbind(out_A_subset, OBS.1_AM_subset,OBS.2_AM_subset)
combined_data$Study <- ifelse(combined_data$Study == 1, "Oral", "IV")

ggplot(combined_data, aes(x = Matrix, y = CAm, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge", width = 1)+
  geom_errorbar(aes(ymin=CAm, ymax = CAm + SD), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Study) +
  scale_fill_manual(values = c("Pre.AM.Oral" = "grey", "Pre.AM.IV" = "#808080",
                               "Obs.AM.Oral" = "#e3716e", "Obs.AM.IV" = "#7ac7e2")) +
  labs(x = "Amniotic Fluid", y = "Concentration of F-53B (μg/mL)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(colour = "black"),  
        legend.position = "none")+
  scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(-0.5, nrow(combined_data) + 0.5))

out_A_subset <- out_M_E_408 %>% select(CF, SD, Study) %>% mutate(Matrix = c("Pre.Fat"))

OBS.1_MF_subset <- OBS.1_MF %>% filter(Time == 408) %>% select(CF, SD, Study) %>% mutate(Matrix = c("Obs.Fat"))
OBS.2_MF_subset <- OBS.2_MF %>% filter(Time == 408) %>% select(CF, SD, Study) %>% mutate(Matrix = c("Obs.Fat"))

combined_data <- rbind(out_A_subset, OBS.1_MF_subset,OBS.2_MF_subset)
combined_data$Study <- ifelse(combined_data$Study == 1, "Oral", "IV")

ggplot(combined_data, aes(x = Matrix, y = CF, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge", width = 1)+
  geom_errorbar(aes(ymin=CF, ymax = CF + SD), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Study) +
  scale_fill_manual(values = c("Pre.Fat.Oral" = "grey", "Pre.Fat.IV" = "#808080",
                               "Obs.Fat.Oral" = "#e3716e", "Obs.Fat.IV" = "#7ac7e2")) +
  labs(x = "Fat", y = "Concentration of F-53B (μg/g)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(colour = "black"),  
        legend.position = "none")+
  scale_y_continuous(limits = c(0, 0.08), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(-0.5, nrow(combined_data) + 0.5))

out_FL_subset <- out_M_E %>%filter(Time %in% c(314, 316, 320, 324, 336, 360, 408)) %>%  
  select(CL_Fet, SD, Study, Time) %>%
  mutate(Matrix = "Pre") 

OBS.1_FL_subset <- OBS.1_FL  %>% select(CL_Fet, SD, Study, Time) %>% mutate(Matrix = c("Obs"))
OBS.2_FL_subset <- OBS.2_FL  %>% select(CL_Fet, SD, Study, Time) %>% mutate(Matrix = c("Obs"))

combined_data <- rbind(out_FL_subset, OBS.1_FL_subset,OBS.2_FL_subset)
combined_data$Study <- ifelse(combined_data$Study == 1, "Oral", "IV")

oral_plot <- ggplot(combined_data %>% filter(Study == "Oral"), 
                    aes(x = Matrix, y = CL_Fet, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = CL_Fet, ymax = CL_Fet + SD), width = .2, position = position_dodge(.9)) +
  facet_grid(Study ~ Time, scales = "free", space = "free") +
  scale_fill_manual(values = c("Pre.Oral" = "grey", "Obs.Oral" = "#e3716e")) +
  labs(x = "Fetal Liver", y = "Concentration of F-53B (μg/g)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))+
  scale_y_continuous(limits = c(0, 0.601), expand = c(0, 0)) 

iv_plot <- ggplot(combined_data %>% filter(Study == "IV"), 
                  aes(x = Matrix, y = CL_Fet, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_errorbar(aes(ymin = CL_Fet, ymax = CL_Fet + SD), width = .2, position = position_dodge(.9)) +
  facet_grid(Study ~ Time, scales = "free", space = "free") +
  scale_fill_manual(values = c("Pre.IV" = "#808080", "Obs.IV" = "#7ac7e2")) +
  labs(x = "Fetal Liver", y = "Concentration of F-53B (μg/g)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        strip.text = element_blank(),
        axis.line = element_line(color = "black", size = 0.5))+
  scale_y_continuous(limits = c(0, 0.601), expand = c(0, 0)) 

combined_plot <- oral_plot / iv_plot
combined_plot

combined_data <- combined_data  %>% filter(Time == 408)

ggplot(combined_data, aes(x = Matrix, y = CL_Fet, fill = interaction(Matrix, Study))) +
  geom_bar(stat = "identity", position = "dodge", width = 1)+
  geom_errorbar(aes(ymin=CL_Fet, ymax = CL_Fet + SD), width = 0.2, position = position_dodge(0.9)) +
  facet_wrap(~Study) +
  scale_fill_manual(values = c("Pre.Oral" = "grey", "Pre.IV" = "#808080",
                               "Obs.Oral" = "#e3716e", "Obs.IV" = "#7ac7e2")) +
  labs(x = "Fetal Liver", y = "Concentration of F-53B (μg/g)") +
  theme_classic() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), 
        axis.line = element_line(colour = "black"),  
        legend.position = "none")+
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0)) + 
  coord_cartesian(xlim = c(-0.5, nrow(combined_data) + 0.5))

######Fig 2.A ######
data_oral <- read_excel("Accumulation.xlsx", sheet = "ORAL")
data_iv <- read_excel("Accumulation.xlsx", sheet = "IV")

data_oral <- data_oral %>%
  mutate(ConcPercent = Conc / sum(Conc)*100)%>%
  mutate(ypos = cumsum(ConcPercent)- 0.5*ConcPercent)%>% 
  mutate(label = paste0(Sample," (",round(ConcPercent,2),"%",")"))

data_iv <- data_iv %>%
  mutate(ConcPercent = Conc / sum(Conc)*100)%>%
  mutate(ypos = cumsum(ConcPercent)- 0.5*ConcPercent)%>% 
  mutate(label = paste0(Sample," (",round(ConcPercent,2),"%",")"))

sample_colors <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", 
                   "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
                   "#FFFF99", "#B15928", "#B3B3B3")

plot_oral <- ggplot(data_oral, aes(x = "", y = ConcPercent, fill = Sample)) +
  geom_bar(stat = "identity", color='white', width = 1) +
  scale_fill_manual(values = sample_colors, 
                    labels=sort(data_oral$label))+ 
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") + 
  geom_text(aes(x=-0.5, y=ypos, label=''))

plot_iv <- ggplot(data_iv, aes(x = "", y = ConcPercent, fill = Sample)) +
  geom_bar(stat = "identity", color='white', width = 1) +
  scale_fill_manual(values = sample_colors, 
                    labels=sort(data_iv$label))+ 
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.position = "right") + 
  geom_text(aes(x=-0.5, y=ypos, label=''))


combined_plot <- plot_oral+ plot_iv

print(combined_plot)
