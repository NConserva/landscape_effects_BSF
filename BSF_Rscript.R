#----------------------------#
#------ BSF in the AF ----------#
#----------------------------#


#packages:                                                                                      #### 
{
  library(corrplot) #
  library(gridExtra) #
  library(ggplot2) #
  library(glmmTMB) #
  library(MASS) #
  library(bbmle) #
  library(lctools) #
  library(tidyr) #
  library(car) #
  library(effects) #
  library(lattice) #
  library(segmented) #
  library(plotly) #
  library(terra) #
  library(sf) #
  library(spdep) #
  library(rgdal) #
  library(DHARMa) #
  library(performance)
  library(forcats)
}

#data:                                                                                          ####
{
  data <- read.csv("BSF_LULC_2001_2019.csv",header=T)
  
  #numeric ordinal and character variables as factor:
  data$CD_NUM <- as.factor(data$CD_NUM)
  data$Year <- as.factor(data$Year)
  
  #NA.omit:                                                                           
  na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
  na_count <- data.frame(na_count)
  
  library(tidyr)
  data <- data %>% drop_na(ED)  
  data <- data %>% drop_na(PD)
  #when no forest cover in the municipality, these variables display NAs (see FRAGSTATS softaware:)
  
  # scaling:
  data[,c(3:12)] <- scale(data[,c(3:12)])

  
}

# LM BSF cases                                                                                  ####
{ # linear regressions to investigate significance of the increase in cases with the years
  #LM data, from DataSUS
  {
  dt <- as.data.frame(cbind(
    c(
    2007,
    2008,
    2009,
    2010,
    2011,
    2012,
    2013,
    2014,
    2015,
    2016,
    2017,
    2018,
    2019),
    c(104,
      94,
      128,
      122,
      153,
      148,
      134,
      182,
      183,
      151,
      189,
      262,
      282),
    c(20,
      22,
      32,
      32,
      59,
      58,
      46,
      73,
      77,
      54,
      69,
      95,
      82)))
  colnames(dt) <- c("Year","BSF_cases","dead")
  }

  plot(dt$Year,dt$BSF_cases)
  LM <- lm(BSF_cases~Year,data=dt)
  summary(LM)
  LM <- lm(dead~Year,data=dt)
  summary(LM)
  dt$fatality <- dt$dead/dt$BSF_cases
  LM <- lm(fatality~Year,data=dt)
  summary(LM)
  # to assess significant increases in BSF cases

} 

#autochtone vs resident cases correlation                                                       ####
{ # average proportion of autochthon cases per year and correlation between autochton and residential 
  # cases (municipality of residence of the patient)
  #LM data, from DataSUS
  {
    dt <- as.data.frame(cbind(
      c(
        2007,
        2008,
        2009,
        2010,
        2011,
        2012,
        2013,
        2014,
        2015,
        2016,
        2017,
        2018,
        2019),
      c(104,
        94,
        128,
        122,
        153,
        148,
        134,
        182,
        183,
        151,
        189,
        262,
        282),
      c(20,
        22,
        32,
        32,
        59,
        58,
        46,
        73,
        77,
        54,
        69,
        95,
        82)))
    colnames(dt) <- c("Year","BSF_cases","dead")
  }
  auto <- read.csv("autoctone cases 2007-2019.csv",header=T,fileEncoding="UTF-8-BOM")
  data <- merge.data.frame(data,auto,by.x=c("CD_NUM","Year"),by.y =c("CD_NUM","Year"),all.y = F,all.x=T)
  library(dplyr)
  data <- mutate_at(data, c(dim(data)[2]), ~replace(., is.na(.), 0))
  d <- rbind(data[data$Year==2007,],data[data$Year==2008,],data[data$Year==2009,],data[data$Year==2010,]
             ,data[data$Year==2011,],data[data$Year==2012,],data[data$Year==2013,],data[data$Year==2014,]
             ,data[data$Year==2015,],data[data$Year==2016,],data[data$Year==2017,],data[data$Year==2018,]
             ,data[data$Year==2019,])
  
  cor(d$BSF_cases,d$autoch,method="spearman") # correlation 
  auto <- as.data.frame(aggregate(data$autoch,by=list(data$Year),sum))
  mean(auto[7:19,2]/dt$BSF_cases) # average proportion of autochthon cases per year
  
}

#collinearity - spearman coefficient:                                                           ####
{
  
library(corrplot)

  # To avoid colinearity issues in models, allowing variables together in a same model when
  # spearman coefficient is < 0.7, spearman coefficient was used because variables did not display
  # normal distributions
Corr <- cor(data[,c(3:12)], method = "spearman", use="pairwise.complete.obs")
is.na(Corr) <- abs(Corr) < 0.7
Corr[is.na(Corr)] <- 0
corrplot(Corr, method = "number", type="upper",pch.cex = 1, number.cex =0.5 ,tl.cex = 1, tl.col="black" 
         ,col=colorRampPalette(c("#BB4444", "#EE9988", "white", "#77AADD", "#4477AA"))(200),)

}

#MODELS:                                                                                        ####

library(glmmTMB) #Generalized linear mixed models glmmTMB package

# Null model (intercepts only):

system.time(mod_null <- glmmTMB(BSF_cases ~ 1 
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data))
save(mod_null,file="mod_null.Rdata")


# landscape composition and climate models:                                                                              ####
{
  # Models
  {
  lc_1 <- glmmTMB(BSF_cases ~ PLAND + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_1,file="lc_1.Rdata")
  
  lc_2 <- glmmTMB(BSF_cases ~ Secondary_Forest + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_2,file="lc_2.Rdata")
  
  lc_3 <- glmmTMB(BSF_cases ~ Riparian_Forest + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_3,file="lc_3.Rdata")
  
  lc_4 <- glmmTMB(BSF_cases ~ PLAND_Pasture + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_4,file="lc_4.Rdata")
  
  lc_5 <- glmmTMB(BSF_cases ~ PLAND_Agriculture + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_5,file="lc_5.Rdata")
  
  lc_6 <- glmmTMB(BSF_cases ~ PLAND_Mosaic_of_Uses + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_6,file="lc_6.Rdata")
  
  lc_7 <- glmmTMB(BSF_cases ~ PLAND + PLAND_Pasture + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_7,file="lc_7.Rdata")
  
  lc_8 <- glmmTMB(BSF_cases ~ PLAND + PLAND_Agriculture + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_8,file="lc_8.Rdata")
  
  lc_9 <- glmmTMB(BSF_cases ~ PLAND + PLAND_Mosaic_of_Uses + temp + prec
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_9,file="lc_9.Rdata")
  
  lc_10 <- glmmTMB(BSF_cases ~ Riparian_Forest + PLAND_Pasture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_10,file="lc_10.Rdata")
  
  lc_11 <- glmmTMB(BSF_cases ~  Riparian_Forest + PLAND_Agriculture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_11,file="lc_11.Rdata")
  
  lc_12 <- glmmTMB(BSF_cases ~  Riparian_Forest + PLAND_Mosaic_of_Uses + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(lc_12,file="lc_12.Rdata")
  
  
  # landscape structure (composition and configuration) and climate models:  ####
  
  
  ls_1 <- glmmTMB(BSF_cases ~ PLAND + PD + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_1,file="ls_1.Rdata")
  
  ls_2 <- glmmTMB(BSF_cases ~ PLAND*PD + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_2,file="ls_2.Rdata")
  
  ls_3 <- glmmTMB(BSF_cases ~ PLAND*PD + PLAND_Pasture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_3,file="ls_3.Rdata")
  
  ls_4 <- glmmTMB(BSF_cases ~ PLAND*PD + PLAND_Agriculture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_4,file="ls_4.Rdata")
  
  ls_5 <- glmmTMB(BSF_cases ~ PLAND*PD + PLAND_Mosaic_of_Uses + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_5,file="ls_5.Rdata")
  
  ls_6 <- glmmTMB(BSF_cases ~ Secondary_Forest*PD + PLAND_Pasture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_6,file="ls_6.Rdata")
  
  ls_7 <- glmmTMB(BSF_cases ~ Secondary_Forest*PD + PLAND_Agriculture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_7,file="ls_7.Rdata")
  
  ls_8 <- glmmTMB(BSF_cases ~ Secondary_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_8,file="ls_8.Rdata")
  
  ls_9 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Pasture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_9,file="ls_9.Rdata")
  
  ls_10 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Agriculture + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_10,file="ls_10.Rdata")
  
  ls_11 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_11,file="ls_11.Rdata")
  
  ls_12 <- glmmTMB(BSF_cases ~ ED + PLAND_Pasture + temp + prec
                   + offset(log(Population)) + (1|Year) + (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_12,file="ls_12.Rdata")
  
  ls_13 <- glmmTMB(BSF_cases ~ ED + PLAND_Agriculture + temp + prec
                   + offset(log(Population)) + (1|Year) + (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_13,file="ls_13.Rdata")
  
  ls_14 <- glmmTMB(BSF_cases ~ ED + PLAND_Mosaic_of_Uses + temp + prec
                   + offset(log(Population)) + (1|Year) + (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(ls_14,file="ls_14.Rdata")
  
  
  # a first AICc selection, to select the best predictors:
  library(bbmle)
  {
    AICctab(
      mod_null,
      
      lc_1,
      lc_2,
      lc_3,
      lc_4,
      lc_5,
      lc_6,
      lc_7,
      lc_8,
      lc_9,
      lc_10,
      lc_11,
      lc_12,
      
      ls_1,
      ls_2,
      ls_3,
      ls_4,
      ls_5,
      ls_6,
      ls_7,
      ls_8,
      ls_9,
      ls_10,
      ls_11,
      ls_12,
      ls_13,
      ls_14,
      
      base=T,delta=T,sort=T,weights=T,nobs=59104)
  }
  
  
  # Trying zero-inflated negative binomiale models with the best AIC selected model from the ones above:
  
  
  ls_11zi <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                   + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~1,family=nbinom1(link="log"), data=data)
  save(ls_11zi,file="ls_11zi.Rdata")
  
  ls_11zi2 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                         + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                         + offset(log(Population)),family=nbinom1(link="log"), data=data ) #note: no random intercept in ziformula
  save(ls_11zi2,file="ls_11zi2.Rdata")
  
  ls_11zi3 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data ) #CONVERGENCE ISSUES when adding 
                                                                                                                                    # random intercept in the ziformula
  save(ls_11zi3,file="ls_11zi3.Rdata")
  
  
  #switching to different optimizer to see if convergence can be achieved:
  
  time <- system.time(
    ls_11zi4 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                        + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                        + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data,
                                        control = glmmTMBControl(optimizer = optim, optArgs = list(method = "L-BFGS-B")))) #25min, CONVERGENCE ISSUES
  save(ls_11zi4,file="ls_11zi4.Rdata")
  
  library(nloptr)
  time <- system.time(
    ls_11zi5 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data,
                                          control = glmmTMBControl(optimizer = bobyqa, optArgs = list(method = "NLOPT_LN_BOBYQA")))) #17min, CONVERGENCE ISSUES
  save(ls_11zi5,file="ls_11zi5.Rdata")
  
  time <- system.time(
    ls_11zi6 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                                          + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data,
                                          control = glmmTMBControl(optimizer = optim, optArgs = list(method = "CG")))) # few days
  
  
  #ziformula with only climate variables that explain structural zero (more ecologically relevant)
  
  ls_11zi6 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                      + offset(log(Population)) + (1|Year)+ (1|CD_NUM), ziformula=~ temp + prec
                      + offset(log(Population)),family=nbinom1(link="log"), data=data ) #note: only climate for probability of structural zero
  save(ls_11zi6,file="ls_11zi6.Rdata")
  
  
  #Exploring different personalized dispersion functions in models:
  
  
  ls_11dis <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                   + offset(log(Population))  + (1|Year)+ (1|CD_NUM),dispformula = ~Year,family=nbinom1(link="log"), data=data)
  save(ls_11dis,file="ls_11dis.Rdata")
  
  
  ls_11dis2 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                    + offset(log(Population))  + (1|Year)+ (1|CD_NUM),dispformula = ~PD,family=nbinom1(link="log"), data=data)
  save(ls_11dis2,file="ls_11dis2.Rdata")
  
  ls_11dis3 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                       + offset(log(Population))  + (1|Year)+ (1|CD_NUM),dispformula = ~ temp+ prec,family=nbinom1(link="log"), data=data)
  save(ls_11dis3,file="ls_11dis3.Rdata")
  
  ls_11diszi3 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                    + offset(log(Population))  + (1|Year)+ (1|CD_NUM),ziformula=~Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                    + offset(log(Population)),dispformula = ~Year,family=nbinom1(link="log"), data=data)
  save(ls_11diszi3,file="ls_11diszi3.Rdata")
  

  ls_11po <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp 
                   + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=poisson(link="log"), data=data)
  save(ls_11po,file="ls_11po.Rdata") #with poisson distribution

  
  ls_11gen <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp + prec
                   + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=genpois(link="log"), data=data)
  save(ls_11gen,file="ls_11gen.Rdata") #with general poisson distribution

  }
  
  # Loading all models previously fitted 
  {
    load("mod_null.Rdata")
    
    load("lc_1.Rdata")
    load("lc_2.Rdata")
    load("lc_3.Rdata")
    load("lc_4.Rdata")
    load("lc_5.Rdata")
    load("lc_6.Rdata")
    load("lc_7.Rdata")
    load("lc_8.Rdata")
    load("lc_9.Rdata")
    load("lc_10.Rdata")
    load("lc_11.Rdata")
    load("lc_12.Rdata")
    
    load("ls_1.Rdata")
    load("ls_2.Rdata")
    load("ls_3.Rdata")
    load("ls_4.Rdata")
    load("ls_5.Rdata")
    load("ls_6.Rdata")
    load("ls_7.Rdata")
    load("ls_8.Rdata")
    load("ls_9.Rdata")
    load("ls_10.Rdata")
    load("ls_11.Rdata")
    load("ls_12.Rdata")
    load("ls_13.Rdata")
    load("ls_14.Rdata")
    
    load("ls_11zi.Rdata")
    load("ls_11zi2.Rdata")
    load("ls_11zi3.Rdata")

    load("ls_11zi5.Rdata")
    load("ls_11zi6.Rdata")
    
    load("ls_11dis.Rdata")
    load("ls_11dis2.Rdata")
    load("ls_11diszi3.Rdata")
    load("ls_11po.Rdata")
    load("ls_11gen.Rdata")

  }
  
  # VIF
  library(performance)
  { #check the VIF for collinearity
  
  check_collinearity(lc_1)
  check_collinearity(lc_2)
  check_collinearity(lc_3)
  check_collinearity(lc_4)
  check_collinearity(lc_5)
  check_collinearity(lc_6)
  check_collinearity(lc_7)
  check_collinearity(lc_8)
  check_collinearity(lc_9)
  check_collinearity(lc_10)
  check_collinearity(lc_11)
  check_collinearity(lc_12)
  
  check_collinearity(ls_1)
  check_collinearity(ls_2)
  check_collinearity(ls_3)
  check_collinearity(ls_4)
  check_collinearity(ls_5)
  check_collinearity(ls_6)
  check_collinearity(ls_7)
  check_collinearity(ls_8)
  check_collinearity(ls_9)
  check_collinearity(ls_10)
  check_collinearity(ls_11)
  check_collinearity(ls_12)
  check_collinearity(ls_13)
  check_collinearity(ls_14)

  check_collinearity(ls_11zi)
  check_collinearity(ls_11zi2) #collinearity
  check_collinearity(ls_11zi6)
  
  check_collinearity(ls_11dis)
  check_collinearity(ls_11dis2)
  check_collinearity(ls_11diszi3) #collinearity

}
  
  # Full AICc model selection
  library(bbmle)
  {
    AICctab(
      mod_null,
      
      lc_1,
      lc_2,
      lc_3,
      lc_4,
      lc_5,
      lc_6,
      lc_7,
      lc_8,
      lc_9,
      lc_10,
      lc_11,
      lc_12,
      
      ls_1,
      ls_2,
      ls_3,
      ls_4,
      ls_5,
      ls_6,
      ls_7,
      ls_8,
      ls_9,
      ls_10,
      ls_11,
      ls_12,
      ls_13,
      ls_14,
      
      ls_11zi,
      ls_11zi2,
      ls_11zi6,
      
      ls_11dis,
      ls_11dis2,
      ls_11diszi3,
      
      base=T,delta=T,sort=T,weights=T,logLik=T,nobs=59126)
    
  }


  # Models assumptions and 'performance'  ####
  library(DHARMa)
  
  #residuals from the variant models of ls_11 (take a few secondes each)
  sim11 <- simulateResiduals(fittedModel = ls_11)
  sim11zi <- simulateResiduals(fittedModel = ls_11zi)
  sim11zi2 <- simulateResiduals(fittedModel = ls_11zi2)
  sim11zi6 <- simulateResiduals(fittedModel = ls_11zi6)
  sim11dis <- simulateResiduals(fittedModel = ls_11dis)
  sim11dis2 <- simulateResiduals(fittedModel = ls_11dis2)
  sim11diszi3 <- simulateResiduals(fittedModel = ls_11diszi3)
  sim11gen <- simulateResiduals(fittedModel = ls_11gen)
  
  # inspection of the default ls_11: Despite slightly significative, DHARMA tests looks ok. Also, tests used here are more sensible when dealing with large dataset 
  {
  sim <- sim11
  plot(sim,form=data$Year) # (slight deviation make the test significant with large dataset, should not be consider as a problem)
  plot(sim) # direct interpretation  of the general QQ plot is difficult because there is many data point
  testDispersion(sim,alternative = "less") # a bit underdispersed d=0.04, p=0.016, this test however is not reliable for underdispersion (DHARMa 'help')
  testZeroInflation(sim) #no zero inflation detected (note: test not always reliable with GLMM)
  
  #plot deviance residuals:
  data$fitted <- predict(ls_11,type='response',se=F)
  m<-ls_11
  {
  plot(data$fitted,residuals(m,type='deviance'))
  abline(a=2,b=0,col="red")
  abline(a=-2,b=0,col="red")
  prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
  title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
  }
  }
  
  # ls_11zi: DHARMa test less good than ls_11
  {
    sim <- sim11zi
    plot(sim,form=data$Year) # seems bit worse than ls_11
    plot(sim)
    testDispersion(sim,alternative = "less") # bit more underdispersed than ls_11, d=0.02, p=0.016
    testZeroInflation(sim) 
    
    #plot deviance residuals:
    data$fitted <- predict(ls_11zi,type='response',se=F)
    m<-ls_11zi
    {
      plot(data$fitted,residuals(m,type='deviance'))
      abline(a=2,b=0,col="red")
      abline(a=-2,b=0,col="red")
      prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
      title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
    }
  }

  # ls_11zi2: Dharma tests quite good, non significant, Howerver VIF detect collinearity
  {
  sim <- sim11zi2
  plot(sim,form=data$Year) # seems better
  plot(sim) # direct interpretation  of the general QQ plot is difficult because there is many data point
  testDispersion(sim,alternative = "less") # less underdispersed d=0.094, p=0.084
  testZeroInflation(sim) 
  
  #plot deviance residuals:
  data$fitted <- predict(ls_11zi2,type='response',se=F)
  m<-ls_11zi2
  {
    plot(data$fitted,residuals(m,type='deviance'))
    abline(a=2,b=0,col="red")
    abline(a=-2,b=0,col="red")
    prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
    title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
  }
  }
  
  # ls_11zi6: Dharma tests quite good, bit better than ls_11, 
  {
    sim <- sim11zi6
    plot(sim,form=data$Year) # seems better
    plot(sim) 
    testDispersion(sim,alternative = "less") # bit more underdispersed d=0.056, p=0.028
    testZeroInflation(sim) 
    
    #plot deviance residuals:
    data$fitted <- predict(ls_11zi6,type='response',se=F)
    m<-ls_11zi6
    {
      plot(data$fitted,residuals(m,type='deviance'))
      abline(a=2,b=0,col="red")
      abline(a=-2,b=0,col="red")
      prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
      title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
    }
  }
  
  # ls_11dis: Tests here are much more significative, moreover,deviance residulas could indicate overfitting
  {
  sim <- sim11dis
  plot(sim,form=data$Year) 
  plot(sim) # direct interpretation  of the general QQ plot is difficult because there is many data point
  testDispersion(sim,alternative = "less") # much more underdispersed d=0.0094, p=2.2e-16
  testZeroInflation(sim) 
  testOutliers(sim)
  
  #plot deviance residuals:
  data$fitted <- predict(ls_11dis,type='response',se=F)
  m<-ls_11dis
  {
    plot(data$fitted,residuals(m,type='deviance'))
    abline(a=0.2,b=0,col="red")
    abline(a=-0.2,b=0,col="red")
    prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -0.2) + sum(residuals(m,type='deviance')> 0.2)))/length(residuals(m,type='deviance'))
    title(c("Porportion of Deviance residuals within -0.2 and 0.2 values:", round(prop,5)))
  } #much smaller deviance residuals, overfitting model?
  }
  
  # ls_11dis2: Dharma tests look all good
  {
    sim <- sim11dis2
    plot(sim,form=data$Year) # everything looks ok
    plot(sim) # direct interpretation  of the general QQ plot is difficult because there is many data point
    testDispersion(sim,alternative = "less") # less underdispersed d=0.016, p=0.04
    testZeroInflation(sim) 
    testOutliers(sim)
    
    #plot deviance residuals:
    data$fitted <- predict(ls_11dis2,type='response',se=F)
    m<-ls_11dis2
    {
      plot(data$fitted,residuals(m,type='deviance'))
      abline(a=2,b=0,col="red")
      abline(a=-2,b=0,col="red")
      prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
      title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
    } #seems not overfitting
  }
  
  # ls_11gen: non uniformity of residuals in year 2002
  {
    sim <- sim11gen
    plot(sim,form=data$Year) # 2002 year not uniform residual distrtibution, but again this is a year with little information
    plot(sim) # direct interpretation  of the general QQ plot is difficult because there is many data point
    testDispersion(sim,alternative = "less") # less underdispersed d=0.044, p=0.032
    testZeroInflation(sim) 
    testOutliers(sim)
  }
  
  # ls_11diszi3:  test worst than ls_11, almost significant outliers, might be overfitting regarding the deviance residuals plot
  {
    sim <- sim11diszi3
    plot(sim,form=data$Year)
    plot(sim) 
    testDispersion(sim,alternative = "less") # more underdispersed than ls_11 d=0.018, p=0.004
    testZeroInflation(sim) 
    testOutliers(sim) #Outliers
    
    #plot deviance residuals:
    data$fitted <- predict(ls_11diszi3,type='response',se=F)
    m<-ls_11diszi3
    {
      plot(data$fitted,residuals(m,type='deviance'))
      abline(a=2,b=0,col="red")
      abline(a=-2,b=0,col="red")
      prop <- (length(residuals(m,type='deviance')) -(sum(residuals(m,type='deviance')< -2) + sum(residuals(m,type='deviance')> 2)))/length(residuals(m,type='deviance'))
      title(c("Porportion of Deviance residuals within -2 and 2 values:", round(prop,5)))
    } #seems not overfitting
  }
  
  
  #Cross-validation with RMSE:
  library(cv)
  time <- system.time(cv_ls11 <- cv(ls_11,data=data, k=10, criterion=rmse)) 
  time <- system.time(cv_ls11zi <-cv(ls_11zi,data=data, k=10, criterion=rmse))
  time <- system.time(cv_ls11zi2 <-cv(ls_11zi2,data=data, k=10, criterion=rmse)) 
  time <- system.time(cv_ls11zi6 <-cv(ls_11zi6,data=data, k=10, criterion=rmse)) 
  time <- system.time(cv_ls11dis <-cv(ls_11dis,data=data, k=10, criterion=rmse)) 
  time <- system.time(cv_ls11dis2 <-cv(ls_11dis2,data=data, k=10, criterion=rmse))
  time <- system.time(cv_ls11diszi3 <-cv(ls_11diszi3,data=data, k=10, criterion=rmse))
  time <- system.time(cv_ls14 <-cv(ls_14,data=data, k=10, criterion=rmse))
  time <- system.time(cv_null <-cv(mod_null,data=data, k=10, criterion=rmse))
  
  save(cv_ls11,file="cv_ls11.RData")
  save(cv_ls11zi,file="cv_ls11zi.RData")
  save(cv_ls11zi2,file="cv_ls11zi2.RData")
  save(cv_ls11zi6,file="cv_ls11zi6.RData")
  save(cv_ls11dis,file="cv_ls11dis.RData")
  save(cv_ls11dis2,file="cv_ls11dis2.RData")
  save(cv_ls11diszi3,file="cv_ls11diszi3.RData")
  save(cv_ls14,file="cv_ls14.RData")
  save(cv_null,file="cv_null.RData")
  
  cv_rmse <- as.data.frame(
    rbind(
      c("cv_ls11",round(cv_ls11$`CV crit`,4),round(cv_ls11$`full crit`,4)),
      c("cv_ls11zi",round(cv_ls11zi$`CV crit`,4),round(cv_ls11zi$`full crit`,4)),
      c("cv_ls11zi2",round(cv_ls11zi2$`CV crit`,4),round(cv_ls11zi2$`full crit`,4)),
      c("cv_ls11zi6",round(cv_ls11zi6$`CV crit`,4),round(cv_ls11zi6$`full crit`,4)),
      c("cv_ls11dis",round(cv_ls11dis$`CV crit`,4),round(cv_ls11dis$`full crit`,4)),
      c("cv_ls11dis2",round(cv_ls11dis2$`CV crit`,4),round(cv_ls11dis2$`full crit`,4)),
      c("cv_ls11diszi3",round(cv_ls11diszi3$`CV crit`,4),round(cv_ls11diszi3$`full crit`,4))
      
      
      )
    )
  colnames(cv_rmse) <- c("cv_model_label","cv_rmse","sample_rmse")
  cv_rmse
  save(cv_rmse,file="cv_rmse.RData")
  load("cv_rmse.RData")
  

  
  
  #rmse doesn't improve much between model, rmse also stay the similar when performing external cv with 2020 and 2021 (done eslewhere)
  
  r2_nakagawa(ls_11zi6) # cond: 0.261 ; marginal : 0.020
  r2_nakagawa(ls_11zi2) # cond: 0.228 ; marginal : 0.024
  r2_nakagawa(ls_11) # cond: 0.270 ; marginal : 0.020
  
  # Pseudo R2 provided for informative purpose, not my religion tho


  
  
  # plot coefficients    ####
  {
    
    m<-ls_11
    summary(m)
    
    data_plot <- as.data.frame(summary(m)$coefficients$cond)[-c(1:2),1]
    data_plot <- cbind(Variables=rownames(summary(m)$coefficients$cond)[-c(1:2)],Coefficients=data_plot)
    data_plot <- as.data.frame(data_plot) 
    colnames(data_plot)[2] <- "Coefficients"
    data_plot$Coefficients <- as.numeric(data_plot$Coefficients)
    conf <- cbind(confint(m,full=F,level=0.95)[,1],confint(ls_11,full=F,level=0.95)[,2])
    conf <- conf[-c(1,2,8,9),]
    conf85 <- cbind(confint(m,full=F,level=0.85)[,1],confint(ls_11,full=F,level=0.85)[,2])
    conf85 <- conf85[-c(1,2,8,9),]
    data_plot <- cbind(data_plot,conf1=conf[,1],conf2=conf[,2],conf85.1=conf85[,1],conf85.2=conf85[,2])
    
    library(forcats)
    data_plot$Variables <- as.factor(c("PD **","PLAND Mosaic of Uses ***","temp ***","prec","Riparian Forest : PD ***"))
    data_plot$Variables<- fct_relevel(data_plot$Variables,c("PD **","Riparian Forest : PD ***","PLAND Mosaic of Uses ***","temp ***","prec"))
    data_plot$Variables<- fct_rev(data_plot$Variables)
    
    plot <- ggplot(data_plot) +
      geom_hline(yintercept = 0,linetype = "dashed") +
      geom_errorbar( aes(x=Variables, ymin=conf1, ymax=conf2,colour=Variables),alpha=0.9, linewidth=1.2, width=0) +
      geom_errorbar( aes(x=Variables, ymin=conf85.1, ymax=conf85.2,colour=Variables), linewidth=2.5, alpha=0.9, width=0) +
      geom_errorbar( aes(x=Variables, ymin=Coefficients, ymax=Coefficients), linewidth=1.25, alpha=0.9, width=0.12) +
      geom_point( aes(x=Variables, y=Coefficients), stat="identity", fill="purple", alpha=0,size=0) +
      scale_alpha_manual(values = c("0.3"=0.3, "0.9"=0.9), guide='none')+
      scale_colour_manual(values = c("temp ***"="#E57373","prec"="#0D47A1","PLAND Mosaic of Uses ***"="#FFD54F","PD **"="#AED581","Riparian Forest : PD ***"="#90CAF9"), guide='none')+
      coord_flip() + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                                        axis.title=element_text(size=15,face="bold"))
    plot
    
    save(plot,file="plot.Rdata")
    
    library(effects)
    plot(allEffects(ls_11)) # overview of the marginal effects 
    
    
  }
  
  # plot predicted values of the observed data against predictors and identification of ecological tresholds by segmented models
  {
    
    # re-import data without scaling :                                                                                          ####
    {
      data <- read.csv("BSF_LULC_2001_2019.csv",header=T)
      
      #numeric ordinal and character variables as factor:
      data$CD_NUM <- as.factor(data$CD_NUM)
      data$Year <- as.factor(data$Year)
      
      #NA.omit:                                                                           
      na_count <-sapply(data, function(y) sum(length(which(is.na(y)))))
      na_count <- data.frame(na_count)
      
      library(tidyr)
      data <- data %>% drop_na(ED)
      data <- data %>% drop_na(PD)
      data <- data %>% drop_na(ENN_MN)
      
    }
    
    load("ls_11.RData")
    pred_ls_11 <- predict(ls_11,type="response",se.fit=T) 
    save(pred_ls_11,file="pred_ls_11.Rdata")
    load("pred_ls_11.Rdata")
    pred <- pred_ls_11

    data$pred <- (pred$fit/data$Population)*100000 #in BSF incidence per 100.000 inhabitants
    data$lci <- data$pred - (((pred$se.fit/data$Population)*100000))
    data$uci <- data$pred + (((pred$se.fit/data$Population)*100000))
  
    # PD:
    plot1 <- ggplot(data,aes(x=PD)) +
      geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred),method="loess",se=F,size=1,color="#AED581")+
      geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      xlab("PD") + ylab("Pred BSF incidence")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot1
    
    # segmented model of PD:
    {
      lo <- loess(pred~PD,data=data)
      data$lo <- lo$fitted
      data_ord <- data[order(data$PD, decreasing = F),]
      
      library(segmented)
      
      model <- lm(lo ~ PD, data = data_ord)
      segmented_model <- segmented(model, seg.Z = ~PD,npsi=5,control = seg.control(n.boot=20,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$PD, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "PD",
        ylab = "Pred BSF incidence",
        col = "lightgray"
      )
      lines(data_ord$PD, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients --> tresholds
      summary(segmented_model) #segmented model need to fit right the loess curve, 
      # adaptation of breakpoints number if needed (npsi) 
      confint.segmented(segmented_model)
      
    }
    
    # PD and Riparian Forest:
    
    lo <- loess(data$pred ~ cbind(data$Riparian_Forest,data$PD))
    lolci <- loess(data$lci ~ cbind(data$Riparian_Forest,data$PD))
    louci <- loess(data$uci ~ cbind(data$Riparian_Forest,data$PD))
    save(lo,file="lo")
    save(lo,file="lolci")
    save(lo,file="louci")
    load("lo")
    load("lolci")
    load("louci")
    data$lo <- lo$fitted
    data$lolci <- lolci$fitted
    data$louci <- louci$fitted
    
    library(rgl) # 3D representation
    
    myColorRamp <- function(col, val) {
      v <- (val - min(val))/diff(range(val))
      x <- colorRamp(col)(v)
      rgb(x[,1], x[,2], x[,3], maxColorValue = 255)
    }
    
    cols <- myColorRamp(c("lightgrey","#64B5F6","#0D47A1","#020C75"),val=data$lo)
    plot3d(y=data$PD, x=data$Riparian_Forest, z=data$lo,col=cols,axes=F,ylab="",zlab="",xlab="",
           xlim=c(0,100),ylim=c(0,8),zlim=c(0,1.2))
    grid3d(c("x+", "y+","z-"), col = "gray", lty = 1, lwd = 1)
    axis3d('x-', tick = T,at = seq(0,100,20))
    axis3d('y-', tick = T,at = seq(0,8,2))
    axis3d('z+', tick = T,at = seq(0,1.2,0.2))
    mtext3d("Riparian Forest", edge="x-", line=4)
    mtext3d("PD", edge="y-", line=4)
    mtext3d("BSF incidence", edge="z+", line=6)


    
    plot2 <- ggplot(data, aes(PD, Riparian_Forest, col = lo)) + 
      geom_point() +
      coord_fixed()+
      theme_bw() +
      coord_flip()+
      ylab("Riparian Forest")+
      labs(color = "BSF\nincidence")+
      #scale_colour_gradient(low = "#D1C4E9",high = "#311B92")+
      scale_colour_gradientn(colours = c("lightgrey","#64B5F6","#0D47A1","#020C75"),values=c(0,0.2,0.6,1.2))+
      theme(legend.position = c(0.9,0.75),axis.text=element_text(size=15,face="bold"),
            axis.title=element_text(size=15,face="bold"))
    
    plot2 # in 2 dimensions
    
    # segmented model of Ripairan Forest
    {
      lo <- loess(pred~Riparian_Forest,data=data)
      data$lo <- lo$fitted
      data_ord <- data[order(data$Riparian_Forest, decreasing = F),]
      
      library(segmented)
      model <- lm(lo ~ Riparian_Forest, data = data_ord)
      segmented_model <- segmented(model, seg.Z = ~Riparian_Forest,npsi=3,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$Riparian_Forest, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "Riparian Forest",
        ylab = "Pred BSF incidence",
        col = "lightgray"
      )
      lines(data_ord$Riparian_Forest, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
      confint.segmented(segmented_model)
      
    }
    
    # Mosaic of Uses:
    
    plot3 <- ggplot(data,aes(x=PLAND_Mosaic_of_Uses)) +
      geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred),method="loess",se=F,size=1,color="#FFD54F")+
      geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
      xlab("PLAND Mosaic of Uses") + ylab("Pred BSF incidence")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot3
    
    # segmented model of Mosaic of Uses:
    {
      
      lo <- loess(pred~PLAND_Mosaic_of_Uses,data=data)
      data$lo <- lo$fitted
      data_ord <- data[order(data$PLAND_Mosaic_of_Uses, decreasing = F),]
      
      library(segmented)
      
      model <- lm(lo ~ PLAND_Mosaic_of_Uses, data = data_ord)
      segmented_model <- segmented(model, seg.Z = ~PLAND_Mosaic_of_Uses,npsi=4,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$PLAND_Mosaic_of_Uses, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "PLAND Mosaic of Uses",
        ylab = "Pred BSF incidence",
        col = "lightgray"
      )
      lines(data_ord$PLAND_Mosaic_of_Uses, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
      confint.segmented(segmented_model)
    }

    
    plot4 <- ggplot(data,aes(x=temp)) +
      #geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred),method="loess",se=F,size=1,color="#E57373")+
      #geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      xlab("temp") + 
      ylab(" ")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot4
    
    # segmented model of temp:
    {
      
      lo <- loess(pred~temp,data=data)
      data$lo <- lo$fitted
      data_ord <- data[order(data$temp, decreasing = F),]
      
      library(segmented)
      model <- lm(lo ~ temp, data = data_ord)
      segmented_model <- segmented(model, seg.Z = ~temp,npsi=2,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$temp, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "temp",
        ylab = "Pred BSF incidence",
        col = "lightgray"
      )
      lines(data_ord$temp, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
      confint.segmented(segmented_model)
    }
    
    
    # prec:
    plot5 <- ggplot(data,aes(x=prec)) +
      geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred),method="loess",se=F,size=1,color="#0D47A1")+
      geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      xlab("prec") + ylab("Pred BSF incidence")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot5
    
    # segmented model of prec:
    {
      
      lo <- loess(pred~prec,data=data)
      data$lo <- lo$fitted
      data_ord <- data[order(data$prec, decreasing = F),]
      
      library(segmented)
      
      model <- lm(lo ~ prec, data = data_ord)
      segmented_model <- segmented(model, seg.Z = ~prec,npsi=4,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$prec, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "prec",
        ylab = "Pred BSF incidence",
        col = "lightgray"
      )
      lines(data_ord$prec, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
      confint.segmented(segmented_model)
    }
    
  }
  
} 


# Moran'I test (spatial autocorrelation):                                                       ####
{
  library(terra) 
  library(joyn)
  library(sf) 
  library(spdep) 
  library(rgdal) 
  library(dplyr)
  
  load("ls_11.Rdata")
  Coord <- read.csv("Coordinates_Centroids_Munic.csv",header=T) # Centroid coordinates of municipality,
  data_coord <- data 
  resid <- simulateResiduals(ls_11) # residuals

  data_coord$resid <- residuals(resid)
  data_coord <- merge(data_coord, Coord, by=c("CD_NUM"),all.x=TRUE)               
  
  shp <- readOGR(dsn = "Zone_studied.shp", layer = "Zone_studied")       # repeat that code for every year
  shp@data <- left_join(shp@data,data_coord[data_coord$Year==2019,],by="CD_NUM")
  
  
  wm_q <- poly2nb(shp, queen = TRUE) # queen neighbourhood relation between the municipalities
  save(wm_q,file="wm_q.Rdata")
  load("wm_q.Rdata")
  
  rswm_q <- nb2listw(wm_q, style = "W", zero.policy = TRUE) #list of neighbours
  
  moran.test(as.numeric(shp$resid), listw = rswm_q, zero.policy = TRUE, na.action = na.omit)
  
  
  plot <- ggplot(data_coord[data_coord$Year==2013,],aes(x=X,y=Y,color=resid)) +
    geom_point()+ labs(title="Model residuals across the Atlantic Forest")+
    theme_bw() + theme(axis.text=element_text(size=15,face="bold"),
                       axis.title=element_text(size=15,face="bold")) +scale_color_gradientn(colors=c("blue","grey","red")) 
  
  plot


  #The Moran's I statistic ranges from -1 to 1. Values in the interval (-1, 0) indicate negative spatial 
  #autocorrelation (low values tend to have neighbours with high values and vice versa), values near 0 
  #indicate no spatial autocorrelation (no spatial pattern - random spatial distribution) and values in 
  #the interval (0,1) indicate positive spatial autocorrelation (spatial clusters of similarly low or high 
  #values between neighbour municipalities should be expected.)
  
  
}


