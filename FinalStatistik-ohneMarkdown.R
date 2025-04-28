### ---------  Final statistics for dependence of metabolic rates on condition
### --- Jan 2025
### Authors: S. Bauer, S Hahn, 

remove.packages("rlang")
install.packages("rlang")
install.packages("brms")
install.packages("Rcpp")
install.packages("RcppParallel")

install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(lme4)
library(arm)
library(brm)
library(brms)


# ---- START: Data preparations ------------------------------------------------
##
## - Reading data
## - Assigning variable types


    s.test <- data.frame(read.csv("/Users/bauersil/Dropbox/MigMalaria_NF/subprojects/Ventotene/Stats_Data/MR_condition/MR_Parasites_Condition_19012025.csv", h=T, sep=";"))
    s <- s.test[c(1:243),c(1:33)]
    str(s)
    
    s$sex <- factor(s$sex, levels=c("1","2"), labels=c("male","female"))
    s$infection <- factor(s$infection, levels=c("0","1"), labels=c("non-inf","infected"))
    
    s.PhoPho <- s[s$species == "PhoPho",]
    s.SaxRub <- s[s$species == "SaxRub",]
    
    d <- s
    d.inf <- d[d$infection == "infected",]
    d.noninf <- d[d$infection == "non-inf",]
    d.lowfuel <- d[d$fuel < 2,]
    
    ### ---- Sample sizes per species, sex and infection status
    table(d$species, d$sex)
    table(d$species, d$infection)
    table(d$species, d$fuel)

# ---- END: Data preparations --------------------------------------------------



#------ START: Parasitemia vs fuel ---------------------------------------------

  m1.parasitemia <- lm(log10(parasitemia) ~ fuel, data = d)
  summary(m1.parasitemia)
  plot(jitter(d$fuel, amount=0.2), log(d$parasitemia), pch=16)

#------ END: Parasitemia vs fuel --------------------------------------


#------START: Analysis Max Metabolic Rates vs fuel score ----------------------

## - Test whether Max Metabolic rate "origMassSpec.PMR" can be fitted as a 
## non-linear function of "fuel"
## - fit the function: MMR = Asymp/(1+ exp(-(fuel - fuel_crit)/Slope)
## - assume the parameters "Asymp" and "cond_crit" are similar for infected/non-infected individuals
## - allow parameter "slope" to vary between infected/non-infected
## - use Bayesian approach


rm(MMR_cond)
MMR_cond <- brm(
  bf(origMassSpec.PMR ~ Asymp / (1 + exp(-(fuel - condcrit)/slope)),
     Asymp ~ 1, condcrit ~ 1, slope ~ infection, 
     nl = TRUE),
  data = d, family = gaussian(), 
  iter = 2000,
  prior = c(
    prior(normal(0.3, 0.1), nlpar = "Asymp"),
    prior(normal(-0.1, 0.2), nlpar = "condcrit"),
    prior(normal(0.9, 0.1), nlpar = "slope")
  ),
  control = list(adapt_delta = 0.9)
)

summary(MMR_cond)
plot(MMR_cond)
summary(MMR_cond)$fixed
conditional_effects(MMR_cond)

## ---- Plotting data, model predictions for infected/uninfected and credibility intervals 

  # simulate new data based on model
  # create newdat: contains the predictor values

    rm(newdat3, newdat4, fitmat3, fitmat4)
    newdat3 <- data.frame(fuel = seq(0,7,l=100), infection = c(rep('non-inf',100)))
    fitmat3 <- posterior_epred(MMR_cond, newdata=newdat3, re.form = NA)
    newdat3$fit <- apply(fitmat3,2,mean)  # the estimate for each row in newdat
    newdat3$lwr <- apply(fitmat3,2,quantile,probs=0.025)  # estimate for lower and
    newdat3$upr <- apply(fitmat3,2,quantile,probs=0.975)  # upper bound of 95% UI
    
    newdat4 <- data.frame(fuel = seq(0,7,l=100), infection = c(rep('infected',100)))
    fitmat4 <- posterior_epred(MMR_cond, newdata=newdat4, re.form = NA)
    newdat4$fit <- apply(fitmat4,2,mean)  # the estimate for each row in newdat
    newdat4$lwr <- apply(fitmat4,2,quantile,probs=0.025)  # estimate for lower and
    newdat4$upr <- apply(fitmat4,2,quantile,probs=0.975)  # upper bound of 95% UI
    
  # plot simulated data based on model together with measured data
    plot.new()
    plot(jitter(d$fuel, amount= 0.22), d$origMassSpec.PMR, pch=20, ylim=c(0.09,0.55),
         xlim=c(0,7),
         xlab="fuel", ylab="origMassSpec.PMR", las=1, col=d$infection)
    lines(newdat3$fuel, newdat3$fit, lwd=2)
    lines(newdat3$fuel, newdat3$lwr, lty=3, lwd=2)
    lines(newdat3$fuel, newdat3$upr, lty=3, lwd=2)
    
    lines(newdat4$fuel, newdat4$fit, lwd=2, col="red")
    lines(newdat4$fuel, newdat4$lwr, lty=3, lwd=2, col="red")
    lines(newdat4$fuel, newdat4$upr, lty=3, lwd=2, col="red")

#------END: Analysis Max Metabolic Rates vs fuel score ----------------------    


    
#-----START: Linear model Max Metabolic rates vs hemoglobin only for low fuel 

## - Calculate linear model for dependence of MMR and hemoglobin, with
##.  accounting for infected and uninfected 
## - plot original data and model-predicted results
    
    
      m.MMR <- lm(origMassSpec.PMR ~ hemo + infection, data = d.lowfuel); summary(m.MMR)
      plot(d.lowfuel$hemo, d.lowfuel$origMassSpec.PMR, col=d.lowfuel$infection, pch=20)
      
      d.lowfuel.infected <- d.lowfuel[d.lowfuel$infection=="infected",]
      m.MMR.infected <-  lm(origMassSpec.PMR ~ hemo, data = d.lowfuel.infected); summary(m.MMR.infected)
      abline(lm(origMassSpec.PMR ~ hemo, data = d.lowfuel.infected), col="red")
      
      
      newx <- data.frame(hemo= seq(min(d.lowfuel.infected$hemo, na.rm=T)-10, max(d.lowfuel.infected$hemo+10, na.rm=T),by = 0.05))
      
      conf_intervals <- predict(m.MMR.infected,newdata=data.frame(hemo=newx), 
                                interval="confidence", level=0.95)
      
      lines(newx$hemo, conf_intervals[,2], col="red", lty="dotted", lwd=1) 
      lines(newx$hemo, conf_intervals[,3], col="red", lty="dotted", lwd=1) 
      
      polygon(
        c(newx$hemo, rev(newx$hemo)),
        c(
          conf_intervals[, 2], 
          rev(conf_intervals[, 3])
        ), 
        col = adjustcolor( "red", alpha.f = 0.3), alpha= 0.4, border = NA)
      
      
      d.lowfuel.noninfected <- d.lowfuel[d.lowfuel$infection=="non-inf",]
      m.MMR.noninfected <-  lm(origMassSpec.PMR ~ hemo, data = d.lowfuel.noninfected); summary(m.MMR.noninfected)
      abline(lm(origMassSpec.PMR ~ hemo, data = d.lowfuel.noninfected), col="grey", lty="longdash")

#-----END: Linear model Max Metabolic rates vs hemoglobin only for low fuel 
      


#-------START: Hemoglobin vs fuel score ----------------------------------------
## - using same rationale as above for MMR
## - using different priors
      
      
    rm(Hemo_cond)
    Hemo_cond <- brm(
      bf(hemo ~ Asymp / (1 + exp(-(fuel - condcrit)/slope)),
         Asymp ~ 1, condcrit ~ 1, slope ~ infection, 
         nl = TRUE),
      data = d, family = gaussian(), 
      iter = 2000,
      prior = c(
        prior(normal(160, 5), nlpar = "Asymp"),
        prior(normal(-0.5, 0.1), nlpar = "condcrit"),
        prior(normal(0.5, 0.1), nlpar = "slope")
      ),
      control = list(adapt_delta = 0.9)
    )
    prior_summary(Hemo_cond)
    
    summary(Hemo_cond)
    plot(Hemo_cond)
    summary(Hemo_cond)$fixed
    conditional_effects(Hemo_cond)

## --- Plotting: model with credibility intervals on top of empirical values

          # simulate new data based on model
          # create newdat: contains the predictor values you need for the plotting
          
          # (if you have >1 predictor, all need to be specified, e.g. set to their
          #  mean, or to the baseline level for factors; for factors, don't forget
          #  to provide the same levels and the same level order as in the dataframe
          #  used for model fitting; if you have transformed predictors, apply the
          #  same transformations in newdat, too)
          
          newdat1 <- data.frame(fuel = seq(0,7,l=100), infection = c(rep('non-inf',100)))
          fitmat1 <- posterior_epred(Hemo_cond, newdata=newdat1, re.form = NA)
          newdat1$fit <- apply(fitmat1,2,mean)  # the estimate for each row in newdat
          newdat1$lwr <- apply(fitmat1,2,quantile,probs=0.025)  # estimate for lower and
          newdat1$upr <- apply(fitmat1,2,quantile,probs=0.975)  # upper bound of 95% UI
          
          newdat2 <- data.frame(fuel = seq(0,7,l=100), infection = c(rep('infected',100)))
          fitmat2 <- posterior_epred(Hemo_cond, newdata=newdat2, re.form = NA)
          newdat2$fit <- apply(fitmat2,2,mean)  # the estimate for each row in newdat
          newdat2$lwr <- apply(fitmat2,2,quantile,probs=0.025)  # estimate for lower and
          newdat2$upr <- apply(fitmat2,2,quantile,probs=0.975)  # upper bound of 95% UI
          
          
          # plot simulated data based on model together with measured data
          plot.new()
          plot(jitter(d$fuel, amount= 0.1), d$hemo, col=d$infection, pch=16, xlab="fuel", ylab="hemoglobin", las=1)
          lines(newdat1$fuel, newdat1$fit, lwd=2)
          lines(newdat1$fuel, newdat1$lwr, lty=3, lwd=2)
          lines(newdat1$fuel, newdat1$upr, lty=3, lwd=2)
          
          lines(newdat2$fuel, newdat2$fit, lwd=2, col="red")
          lines(newdat2$fuel, newdat2$lwr, lty=3, lwd=2, col="red")
          lines(newdat2$fuel, newdat2$upr, lty=3, lwd=2, col="red")

                  
#-------END: Hemoglobin vs fuel score ----------------------------------------
          

# -------START: RMR vs fuel scores -------------------------------------------


    m1.RMR <- lm(origMassSpec.BMR ~ fuel + infection, data = d)
    summary(m1.RMR)
          
   ## - plotting
   plot(jitter(d$fuel, amount= 0.22), d$origMassSpec.BMR, pch=20, ylim=c(0.03,0.09),
         xlim=c(-0.2,7.2), xaxp = c(0, 7, 7),
         xlab="fuel", ylab="BMR", las=1, col=d$infection)
    