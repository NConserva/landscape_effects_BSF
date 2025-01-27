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
  data <- data %>% drop_na(ENN_MN)
  #when no forest cover in the municipality, these variables display NAs
  
  # scaling:
  data[,c(3:22)] <- scale(data[,c(3:22)])

  
}

# LM BSF cases                                                                                  ####
{
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
{
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
Corr <- cor(data[,c(3:22)], method = "spearman", use="pairwise.complete.obs")
is.na(Corr) <- abs(Corr) < 0.7
Corr[is.na(Corr)] <- 0
corrplot(Corr, method = "number", type="upper",pch.cex = 1, number.cex =0.5 ,tl.cex = 1, tl.col="black" 
         ,col=colorRampPalette(c("#BB4444", "#EE9988", "white", "#77AADD", "#4477AA"))(200))

}

#MODELS:                                                                                        ####

library(glmmTMB)

mod_null <- glmmTMB(BSF_cases ~ 1 
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
save(mod_null,file="mod_null.Rdata")

# choice of the most relevant climatic variable to include in the LULC models:                  ####
{
  Corr <- cor(data[,c(19:22)], method = "spearman", use="pairwise.complete.obs")
  is.na(Corr) <- abs(Corr) < 0.7
  Corr[is.na(Corr)] <- 0
  corrplot(Corr, method = "number", type="upper",pch.cex = 1, number.cex = ,tl.cex = 1, tl.col="black" 
           ,col=colorRampPalette(c("#BB4444", "#EE9988", "gray", "#77AADD", "#4477AA"))(200))
  
  clima_1 <- glmmTMB(BSF_cases ~ temp 
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_1,file="clima_1.Rdata")
  
  clima_2 <- glmmTMB(BSF_cases ~ prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_2,file="clima_2.Rdata")
  
  clima_3 <- glmmTMB(BSF_cases ~ temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_3,file="clima_3.Rdata")
  
  clima_4 <- glmmTMB(BSF_cases ~ cold_temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_4,file="clima_4.Rdata")
  
  clima_5 <- glmmTMB(BSF_cases ~ hot_temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_5,file="clima_5.Rdata")
  
  clima_6 <- glmmTMB(BSF_cases ~ cold_temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_6,file="clima_6.Rdata")
  
  clima_7 <- glmmTMB(BSF_cases ~ hot_temp + prec
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(clima_7,file="clima_7.Rdata")
  
  load("mod_null.Rdata")
  load("clima_1.Rdata")
  load("clima_2.Rdata")
  load("clima_3.Rdata")
  load("clima_4.Rdata")
  load("clima_5.Rdata")
  load("clima_6.Rdata")
  load("clima_7.Rdata")
  
  library(bbmle)
  
  tab <- AICctab(mod_null,clima_1,clima_2,clima_3,clima_4,clima_5,clima_6,clima_7
                 ,base=T,delta=T,sort=T,weights=T,nobs=59104)
  
  tab <- as.data.frame(tab)
  tab
  
  summary(modclima_1)
  summary(modclima_4)
  summary(modclima_3)
  
}
  

# landscape model:                                                                              ####
{
  # Models
  {
  mod_null <- glmmTMB(BSF_cases ~ 1 
                      + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(mod_null,file="mod_null.Rdata")
  
  land_1 <- glmmTMB(BSF_cases ~ PLAND + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_1,file="land_1.Rdata")
  
  land_2 <- glmmTMB(BSF_cases ~ Primary_Forest + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_2,file="land_2.Rdata")
  
  land_3 <- glmmTMB(BSF_cases ~ Secondary_Forest + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_3,file="land_3.Rdata")
  
  land_4 <- glmmTMB(BSF_cases ~ PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_4,file="land_4.Rdata")
  
  land_5 <- glmmTMB(BSF_cases ~ Primary_Forest*PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_5,file="land_5.Rdata")
  
  land_6 <- glmmTMB(BSF_cases ~ Secondary_Forest*PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_6,file="land_6.Rdata")
  
  land_7 <- glmmTMB(BSF_cases ~ PLAND*PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_7,file="land_7.Rdata")
  
  land_8 <- glmmTMB(BSF_cases ~ Primary_Forest + PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_8,file="land_8.Rdata")
  
  land_9 <- glmmTMB(BSF_cases ~ Secondary_Forest + PD + temp
                    + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_9,file="land_9.Rdata")
  
  land_10 <- glmmTMB(BSF_cases ~ PLAND + PD + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_10,file="land_10.Rdata")
  
  land_11 <- glmmTMB(BSF_cases ~  ED + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_11,file="land_11.Rdata")
  
  land_12 <- glmmTMB(BSF_cases ~  ENN_MN + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_12,file="land_12.Rdata")
  
  land_13 <- glmmTMB(BSF_cases ~  PD + ENN_MN + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_13,file="land_13.Rdata")
  
  land_14 <- glmmTMB(BSF_cases ~  PD + ED + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_14,file="land_14.Rdata")
  
  land_15 <- glmmTMB(BSF_cases ~  Riparian_Forest + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_15,file="land_15.Rdata")
  
  land_16 <- glmmTMB(BSF_cases ~ PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_16,file="land_16.Rdata")
  
  land_17 <- glmmTMB(BSF_cases ~  PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_17,file="land_17.Rdata")
  
  land_18 <- glmmTMB(BSF_cases ~ PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_18,file="land_18.Rdata")
  
  land_19 <- glmmTMB(BSF_cases ~ Secondary_Forest + PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_19,file="land_19.Rdata")
  
  land_20 <- glmmTMB(BSF_cases ~ Secondary_Forest + PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_20,file="land_20.Rdata")
  
  land_21 <- glmmTMB(BSF_cases ~ Secondary_Forest + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_21,file="land_21.Rdata")
  
  land_22 <- glmmTMB(BSF_cases ~ Primary_Forest + PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_22,file="land_22.Rdata")
  
  land_23 <- glmmTMB(BSF_cases ~ Primary_Forest + PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_23,file="land_23.Rdata")
  
  land_24 <- glmmTMB(BSF_cases ~ Primary_Forest + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_24,file="land_24.Rdata")
  
  land_25 <- glmmTMB(BSF_cases ~ Riparian_Forest + PD + PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_25,file="land_25.Rdata")
  
  land_26 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_26,file="land_26.Rdata")
  
  land_27 <- glmmTMB(BSF_cases ~ Riparian_Forest + PLAND_Agriculture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_27,file="land_27.Rdata")
  
  land_28 <- glmmTMB(BSF_cases ~ Riparian_Forest + PD + PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_28,file="land_28.Rdata")
  
  land_29 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_29,file="land_29.Rdata")
  
  land_30 <- glmmTMB(BSF_cases ~ Riparian_Forest + PLAND_Pasture + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_30,file="land_30.Rdata")
  
  land_31 <- glmmTMB(BSF_cases ~ Riparian_Forest + PD + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_31,file="land_31.Rdata")
  
  land_32 <- glmmTMB(BSF_cases ~ Riparian_Forest*PD + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_32,file="land_32.Rdata")
  
  land_33 <- glmmTMB(BSF_cases ~ Riparian_Forest + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_33,file="land_33.Rdata")
  
  land_34 <- glmmTMB(BSF_cases ~ Primary_Forest*PD + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_34,file="land_34.Rdata")
  
  land_35 <- glmmTMB(BSF_cases ~ Secondary_Forest*PD + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_35,file="land_35.Rdata")
  
  land_36 <- glmmTMB(BSF_cases ~ PLAND*PD + PLAND_Mosaic_of_Uses + temp
                     + offset(log(Population)) + (1|Year)+ (1|CD_NUM),family=nbinom1(link="log"), data=data)
  save(land_36,file="land_36.Rdata")
  
  }
  
  # Loading models
  {
    load("mod_null.Rdata")
    load("land_1.Rdata")
    load("land_2.Rdata")
    load("land_3.Rdata")
    load("land_4.Rdata")
    load("land_5.Rdata")
    load("land_6.Rdata")
    load("land_7.Rdata")
    load("land_8.Rdata")
    load("land_9.Rdata")
    load("land_10.Rdata")
    load("land_11.Rdata")
    load("land_12.Rdata")
    load("land_13.Rdata")
    load("land_14.Rdata")
    load("land_15.Rdata")
    load("land_16.Rdata")
    load("land_17.Rdata")
    load("land_18.Rdata")
    load("land_19.Rdata")
    load("land_20.Rdata")
    load("land_21.Rdata")
    load("land_22.Rdata")
    load("land_23.Rdata")
    load("land_24.Rdata")
    load("land_25.Rdata")
    load("land_26.Rdata")
    load("land_27.Rdata")
    load("land_28.Rdata")
    load("land_29.Rdata")
    load("land_30.Rdata")
    load("land_31.Rdata")
    load("land_32.Rdata")
    load("land_33.Rdata")
    load("land_34.Rdata")
    load("land_35.Rdata")
    load("land_36.Rdata")
    
  }
  
  # VIF
  {
  library(performance)
  
  check_collinearity(land_1)
  check_collinearity(land_2)
  check_collinearity(land_3)
  check_collinearity(land_4)
  check_collinearity(land_5)
  check_collinearity(land_6)
  check_collinearity(land_7)
  check_collinearity(land_8)
  check_collinearity(land_9)
  check_collinearity(land_10)
  check_collinearity(land_11)
  check_collinearity(land_12)
  check_collinearity(land_13)
  check_collinearity(land_14)
  check_collinearity(land_15)
  check_collinearity(land_16)
  check_collinearity(land_17)
  check_collinearity(land_18)
  check_collinearity(land_19)
  check_collinearity(land_20)
  check_collinearity(land_21)
  check_collinearity(land_22)
  check_collinearity(land_23)
  check_collinearity(land_24)
  check_collinearity(land_25)
  check_collinearity(land_26)
  check_collinearity(land_27)
  check_collinearity(land_28)
  check_collinearity(land_29)
  check_collinearity(land_30)
  check_collinearity(land_31)
  check_collinearity(land_32)
  check_collinearity(land_33)
  check_collinearity(land_34)
  check_collinearity(land_35)
  check_collinearity(land_36)
}
  
  # AICc selection
  {
    AICctab(
      mod_null,
      land_1,
      land_2,
      land_3,
      land_4,
      land_5,
      land_6,
      land_7,
      land_8,
      land_9,
      land_10,
      land_11,
      land_12,
      land_13,
      land_14,
      land_15,
      land_16,
      land_17,
      land_18,
      land_19,
      land_20,
      land_21,
      land_22,
      land_23,
      land_24,
      land_25,
      land_26,
      land_27,
      land_28,
      land_29,
      land_30,
      land_31,
      land_32,
      land_33,
      land_34,
      land_35,
      land_36,
      base=T,delta=T,sort=T,weights=T,nobs=59104)
    
  }

  summary(land_32)
  
  # plot coefficients 
  {
    
    data_plot <- as.data.frame(summary(land_32)$coefficients$cond)[-c(1:2),1]
    data_plot <- cbind(Variables=rownames(summary(land_32)$coefficients$cond)[-c(1:2)],Coefficients=data_plot)
    data_plot <- as.data.frame(data_plot) 
    colnames(data_plot)[2] <- "Coefficients"
    data_plot$Coefficients <- as.numeric(data_plot$Coefficients)
    conf <- cbind(confint(land_32,full=F,level=0.95)[,1],confint(land_32,full=F,level=0.95)[,2])
    conf <- conf[-c(1,2,7,8),]
    conf85 <- cbind(confint(land_32,full=F,level=0.85)[,1],confint(land_32,full=F,level=0.85)[,2])
    conf85 <- conf85[-c(1,2,7,8),]
    data_plot <- cbind(data_plot,conf1=conf[,1],conf2=conf[,2],conf85.1=conf85[,1],conf85.2=conf85[,2])
    
    data_plot$Variables <- as.factor(c("PD **","PLAND Mosaic of Uses ***","temp ***","Riparian Forest : PD ***"))
    data_plot$Variables<- fct_relevel(data_plot$Variables,c("PD **","Riparian Forest : PD ***","PLAND Mosaic of Uses ***","temp ***"))
    data_plot$Variables<- fct_rev(data_plot$Variables)
    
    plot <- ggplot(data_plot) +
      geom_hline(yintercept = 0,linetype = "dashed") +
      geom_errorbar( aes(x=Variables, ymin=conf1, ymax=conf2,colour=Variables),alpha=0.9, linewidth=1.2, width=0) +
      geom_errorbar( aes(x=Variables, ymin=conf85.1, ymax=conf85.2,colour=Variables), linewidth=2.5, alpha=0.9, width=0) +
      geom_errorbar( aes(x=Variables, ymin=Coefficients, ymax=Coefficients), linewidth=1.25, alpha=0.9, width=0.12) +
      geom_point( aes(x=Variables, y=Coefficients), stat="identity", fill="purple", alpha=0,size=0) +
      scale_alpha_manual(values = c("0.3"=0.3, "0.9"=0.9), guide='none')+
      scale_colour_manual(values = c("temp ***"="#E57373","PLAND Mosaic of Uses ***"="#FFD54F","PD **"="#AED581","Riparian Forest : PD ***"="#90CAF9"), guide='none')+
      coord_flip() + theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                                        axis.title=element_text(size=15,face="bold"))
    plot
    
    save(plot,file="plot.Rdata")
    
    library(effects)
    plot(allEffects(land_32)) # overview of the marginal effects 
    
    # Contour plot (better representation of the interaction effect Riparian_Forest:PD)
    {
      library(lattice)
  
      rip <- rep(seq(0,100,1),each=101)
      PD<-rep(seq(0,8,0.08),101) 
      predT<-exp(fixef(land_34)$cond[1]+fixef(land_34)$cond[2]*scale(rip)[,1]+fixef(land_34)$cond[3]*scale(PD)[,1] 
                 +fixef(land_34)$cond[6]*scale(rip*PD)[,1])

      d <- as.data.frame(cbind(PD,rip,predT))
      d$predT <- d$predT*(10^8)
      
      #color palette for contour plot
      {
        cvi_colours = list(
          cvi_purples = c("#381532", "#4b1b42", "#5d2252", "#702963",
                          "#833074", "#953784", "#a83e95"),
          my_favourite_colours = c("#702963", "#637029",    "#296370"),
          cvi_PD_rip = c("#1B5E20","#AED581","#64B5F6","#0D47A1"),
          pur_PD_rip = c("lightgrey","#953784","#5d2252","#381532"),
          red_PD_rip = c("#fedccd","#fcb398","#fc8666","#f6573e","#dd2a25","#b31218","#67000d")
        )
        
        cvi_palettes = function(name, n, all_palettes = cvi_colours, type = c("discrete", "continuous")) {
          palette = all_palettes[[name]]
          if (missing(n)) {
            n = length(palette)
          }
          type = match.arg(type)
          out = switch(type,
                       continuous = grDevices::colorRampPalette(palette)(n),
                       discrete = palette[1:n]
          )
          structure(out, name = name, class = "palette")
        }
        # save in ggplot library:
        cvi_palettes("my_favourite_colours", type = "discrete")
        
        scale_colour_cvi_d = function(name) {
          ggplot2::scale_colour_manual(values = cvi_palettes(name,type = "discrete"))
        }
        scale_fill_cvi_d = function(name) {
          ggplot2::scale_fill_manual(values = cvi_palettes(name,
                                                           type = "discrete"))
        }
        
        
        scale_colour_cvi_c = function(name) {
          ggplot2::scale_colour_gradientn(colours = cvi_palettes(name = name,
                                                                 type = "continuous"))
        }
        scale_fill_cvi_c = function(name) {
          ggplot2::scale_fill_gradientn(colours = cvi_palettes(name = name,
                                                               type = "continuous"))
        }
        scale_color_cvi_d = scale_colour_cvi_d
        scale_color_cvi_c = scale_colour_cvi_c
        
      }
      
      plot(cvi_palettes("red_PD_rip", type = "continuous", n = 10))
      col <- c(cvi_palettes("red_PD_rip", type = "continuous", n =2560))
      
      conplot <- ggplot(d,aes(x=rip,y=PD,z=predT)) +
        ylim(0,NA)+ xlab("Riparian Forest") + ylab("PD")+ xlim(0,NA)+
        theme_bw() + theme(axis.text=element_text(size=15,face="bold"),
                           axis.title=element_text(size=15,face="bold")) +
        geom_contour_filled(bins = 10,fill=col,col='black') +
        guides( fill=guide_legend(title="BSF incidence"))
      conplot 

    }
    
  }
  
  # plot predicted values and identification of ecological tresholds by segmented models
  {
    
    #data without scaling :                                                                                          ####
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
    
    load("land_32.RData")
    pred <- predict(land_32,type="response",se.fit=T) 
    save(pred,file="pred.Rdata")
    load("pred.Rdata")
    
    #cor(data$BSF_cases, pred$fit)
    #cor((data$BSF_cases/data$Population)*100000,(pred$fit/data$Population)*100000)

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
      data <- data[order(data$PD, decreasing = F),]
      
      library(segmented)
      
      model <- lm(lo ~ PD, data = data)
      segmented_model <- segmented(model, seg.Z = ~PD,npsi=5,control = seg.control(n.boot=20,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$PD, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "Independent Variable (x)",
        ylab = "Dependent Variable (y)",
        col = "blue"
      )
      lines(data$PD, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients --> tresholds
      summary(segmented_model) #segmented model need to fit right the loess curve, 
      # adaptation of breakpoints number if needed (npsi) 
    }
    
    # PD and Riparian Forest:
    library(plotly)
    
    lo <- loess(data$pred ~ cbind(data$Riparian_Forest,data$PD))
    save(lo,file="lo")
    load("lo.Rdata")
    data$lo <- lo$fitted
    
    plot_ly(y=data$PD, x=data$Riparian_Forest, z=data$lo,col = 'blue') # 3D representation
    
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
      data <- data[order(data$Riparian_Forest, decreasing = F),]
      library(segmented)
      model <- lm(lo ~ Riparian_Forest, data = data)
      segmented_model <- segmented(model, seg.Z = ~Riparian_Forest,npsi=3,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$Riparian_Forest, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "Independent Variable (x)",
        ylab = "Dependent Variable (y)",
        col = "blue"
      )
      lines(data$Riparian_Forest, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
    }
    
    # Mosaic of Uses:
    
    plot3 <- ggplot(data,aes(x=PLAND_Mosaic_of_Uses)) +
      geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred,weight=w),method="loess",se=F,size=1,color="#FFD54F")+
      geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      scale_y_continuous(breaks = seq(0, 1, by = 0.1))+
      xlab("PLAND Mosaic of Uses") + ylab("Pred BSF incidence")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot3
    
    # segmented model of Mosaic of Uses:
    {
      library(segmented)
      
      model <- lm(lo ~ PLAND_Mosaic_of_Uses, data = data)
      segmented_model <- segmented(model, seg.Z = ~PLAND_Mosaic_of_Uses,npsi=4,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$PLAND_Mosaic_of_Uses, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "Independent Variable (x)",
        ylab = "Dependent Variable (y)",
        col = "blue"
      )
      lines(data$PLAND_Mosaic_of_Uses, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
    }

    
    plot4 <- ggplot(data,aes(x=temp)) +
      geom_smooth(aes(y=lci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      geom_smooth(aes(y=pred),method="loess",se=F,size=1,color="#E57373")+
      geom_smooth(aes(y=uci),method="loess",se=F,size=1,linetype = "dashed",color="grey")+
      coord_cartesian(xlim = c(NA, NA), ylim = c(0, NA))+
      xlab("temp") + 
      ylab(" ")+
      theme_bw() + theme(legend.position = "none",axis.text=element_text(size=15,face="bold"),
                         axis.title=element_text(size=15,face="bold"))
    
    plot4
    
    # segmented model of temp:
    {
      library(segmented)
      model <- lm(lo ~ temp, data = data)
      segmented_model <- segmented(model, seg.Z = ~temp,npsi=2,control = seg.control(n.boot=10,quant=T))
      seg_preds <- predict(segmented_model)
      seg_res <- data$lo - seg_preds
      # Plot the original data with the fitted model
      plot(
        data$temp, data$lo,
        main = "Piecewise Regression Fit",
        xlab = "Independent Variable (x)",
        ylab = "Dependent Variable (y)",
        col = "blue"
      )
      lines(data$temp, seg_preds,col = "red", lwd = 2)
      # View breakpoints and coefficients
      summary(segmented_model)
    }
    
    
  }
  
} 


# Moran'I test (spatial autocorrelation):                                                       ####
{
  library(terra) 
  library(sf) 
  library(spdep) 
  library(rgdal) 
  library(DHARMa) 
  
  load("land_32.Rdata")
  Coord <- read.csv("Coordinates_Centroids_Munic.csv",header=T) # Centroid coordinates of municipality,
  data_coord <- data 
  resid <- simulateResiduals(land_32) # residuals

  data_coord$resid <- residuals(resid)
  data_coord <- merge(data_coord, Coord, by=c("CD_NUM"),all.x=TRUE)
  data_coord <- data_coord[data_coord$Year==2019,]                      # repeat that code for every year
  
  shp <- readOGR(dsn = "Zone_studied.shp", layer = "Zone_studied")
  
  shp@data <- left_join(shp@data,data_coord,by="CD_NUM")
  
  wm_q <- poly2nb(shp, queen = TRUE) # queen neighbourhood relation between the municipalities
  save(wm_q,file="wm_q.Rdata")
  load("wm_q.Rdata")
  
  rswm_q <- nb2listw(wm_q, style = "W", zero.policy = TRUE) #list of neighbours
  
  moran.test(as.numeric(shp$resid), listw = rswm_q, zero.policy = TRUE, na.action = na.omit)


  #The Moran's I statistic ranges from -1 to 1. Values in the interval (-1, 0) indicate negative spatial 
  #autocorrelation (low values tend to have neighbours with high values and vice versa), values near 0 
  #indicate no spatial autocorrelation (no spatial pattern - random spatial distribution) and values in 
  #the interval (0,1) indicate positive spatial autocorrelation (spatial clusters of similarly low or high 
  #values between neighbour municipalities should be expected.)
  
}

  

