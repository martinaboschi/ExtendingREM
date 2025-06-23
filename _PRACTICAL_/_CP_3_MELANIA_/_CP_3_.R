## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE, show = FALSE-------------------------------
if (!require("mgcv", quietly = TRUE)) {
  install.packages("mgcv")
  library("mgcv")
} else {
  if (!require("mgcViz", quietly = TRUE)) {
    install.packages("mgcViz")
    library("mgcViz")
  } else {
    if (!require("ggplot2", quietly = TRUE)) {
      install.packages("ggplot2")
      library("ggplot2")
    } else {
      if (!require("survival", quietly = TRUE)) {
        install.packages("survival")
        library("survival")
      } else {
        if (!require("RColorBrewer", quietly = TRUE)){
          install.packages("RColorBrewer")
        } else {
          if (!require("dplyr", quietly = TRUE)){
          install.packages("dplyr")
          } else {
            library("mgcv")
            library("mgcViz")
            library("ggplot2")
            library("survival")
            library("RColorBrewer")
            library(dplyr)
          }
        }
      }
    }
  }
}

## -----------------------------------------------------------------------------
load("_INPUT_CP_3_/dat_gam_1.RData")
head(dat_gam_1)

## -----------------------------------------------------------------------------
gam_fit_tve <- gam(y ~ 
                diff_female
               #  + female
               #individual_activity
               + dyadic_activity
               + closure
               + s(EVENT_INTERVAL, by = individual_activity)
               #+ s(EVENT_INTERVAL, by = dyadic_activity)
               #+ s(EVENT_INTERVAL, by = closure)
                - 1
               , data = dat_gam_1, 
               family="binomial")
summary(gam_fit_tve)

plot(gam_fit_tve)

## -----------------------------------------------------------------------------
t.mat <- cbind(rep(1, nrow(dat_gam_1)),rep(-1, nrow(dat_gam_1)))
gam_fit_nle <- gam(y ~ 
                diff_female
               #  + female
              #individual_activity
               + dyadic_activity
               +closure
               +  s(cbind(individual.activity_ev,individual.activity_nv), 
                    by = t.mat) 
              #+ s(cbind(dyadic.activity_ev, dyadic.activity_nv), by = t.mat)
              #+ s(cbind(closure_ev,closure_nv),by = t.mat)
                - 1
               , data = dat_gam_1, 
               family="binomial")
summary(gam_fit_nle)

plot(gam_fit_nle)

## -----------------------------------------------------------------------------
t.mat <- cbind(rep(1, nrow(dat_gam_1)),rep(-1, nrow(dat_gam_1)))
gam_fit_nle_log <- gam(y ~ 
                diff_female
               #  + female
              #individual_activity
               + dyadic_activity
               +closure
               +  s(cbind(log(individual.activity_ev+1),
                          log(individual.activity_nv+1)), by = t.mat) 
              #+ s(cbind(dyadic.activity_ev, dyadic.activity_nv), by = t.mat)
              #+ s(cbind(closure_ev,closure_nv),by = t.mat)
                - 1
               , data = dat_gam_1, 
               family="binomial")
summary(gam_fit_nle)
plot(gam_fit_nle_log)


ind_act<- sort(dat_gam_1$individual.activity_ev)
pred_nle<-predict(gam_fit_nle_log,
                  newdata = list(individual.activity_ev = log(ind_act + 1), 
                                 individual.activity_nv = rep(0,length(ind_act)),
                                 t.mat = cbind(rep(1, length(ind_act)),
                                               rep(0, length(ind_act))),
                                 closure =rep(0,length(ind_act)), 
                                 dyadic_activity=rep(0, length(ind_act)),
                                 diff_female = rep(0,length(ind_act))
                                 ), se.fit =T, 
                  type = "terms")

par(mar=c(5,6,4,1)+1, mgp=c(5,2,0))
plot(ind_act, pred_nle$fit[,4], t='l',xlab ="Individual activity", 
     ylab = "Non-linear effect on characters' co-occurence", ylim = c(-5,0)
     ,cex.lab=1, cex.axis=1.2, cex.main=3,col = "black",lwd = 2) 
lines(ind_act,pred_nle$fit[,4] + 1.96*pred_nle$se.fit[,4], 
      lty = 2, col = "black", lwd = 2)
lines(ind_act,pred_nle$fit[,4] - 1.96*pred_nle$se.fit[,4], 
      lty = 2, col = "black", lwd = 2)
legend("bottomright",c("Estimate", "Conf. interval"),
       lty = c(1,2),lwd = c(2,2), col = c("black","black"), cex=1.1, bty = "n")

## -----------------------------------------------------------------------------
source_factor <- factor(c(dat_gam_1$SOURCE_ev,dat_gam_1$SOURCE_nv))
dim(source_factor) <- c(nrow(dat_gam_1),2)

gam_fit_re <- gam(y ~ s(source_factor, by = t.mat, bs = "re") +
                    diff_female
                    #  + female
                    # + individual_activity
                    + dyadic_activity
                    #+ triadic_activity
                    + closure
                    #+ s(cbind(log(individual.activity_ev+1),
                    #+ log(individual.activity_nv+1)), by = t.mat) 
                    + s(EVENT_INTERVAL, by = individual_activity)
                    - 1
               , data = dat_gam_1, 
               family="binomial", method = "REML")

summary(gam_fit_re)
plot(gam_fit_re,select = 2)
gam.vcomp(gam_fit_re)

re_coeff<-coefficients(gam_fit_re)[4:412]
names(re_coeff)<-levels(source_factor)
print(sort(re_coeff, decreasing = TRUE)[1:5])
hist(re_coeff, breaks = 50)

## -----------------------------------------------------------------------------
save(dat_gam_1, 
     gam_fit_tve,
     gam_fit_nle,
     gam_fit_nle_log,
     gam_fit_re, file="_OUTPUT_CP_3_/nonlinear_models.RData")

