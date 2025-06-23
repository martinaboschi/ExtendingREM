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
source("_INPUT_CP_4_/_FUNCTIONS_.R")

## -----------------------------------------------------------------------------
load("_INPUT_CP_4_/gam_fit.RData")
load("_INPUT_CP_4_/dat_gam_1.RData")

## -----------------------------------------------------------------------------
class(gam_fit)
summary(gam_fit)

## -----------------------------------------------------------------------------
gam_fit <- gam(formula = y ~ diff_female + individual_activity + dyadic_activity + closure - 
    1, family = "binomial", data = dat_gam_1)
summary(gam_fit)

## -----------------------------------------------------------------------------
load("_INPUT_CP_4_/nonlinear_models.RData")

## -----------------------------------------------------------------------------
head(dat_gam_1)

## -----------------------------------------------------------------------------
gam_fit$formula

## -----------------------------------------------------------------------------
gam_fit_nle$formula
gam_fit_nle_log$formula
gam_fit_tve$formula
gam_fit_re$formula

## -----------------------------------------------------------------------------
set.seed(1234)
BB9 <- BB.single(dim.k = 9)
BB10 <- BB.single(dim.k = 10)


## -----------------------------------------------------------------------------
# linear model 
AIC(gam_fit)
# model including time-varying effect for individual activity
AIC(gam_fit_tve)
# model including non-linear effect for individual activity
AIC(gam_fit_nle) 
# model including non-linear effect for log-individual activity
AIC(gam_fit_nle_log) 
# model including random effect for source actor
AIC(gam_fit_re)

## -----------------------------------------------------------------------------
GOF_linear_DF = GOF_univariate(gam.fit = gam_fit, 
                            index = 1)
GOF_linear_DF[[1]]

## -----------------------------------------------------------------------------
GOF_linear_IA = GOF_univariate(gam.fit = gam_fit, 
                            index = 2)
GOF_linear_IA[[1]]

## -----------------------------------------------------------------------------
GOF_linear_DA = GOF_univariate(gam.fit = gam_fit, 
                            index = 3)
GOF_linear_DA[[1]]

## -----------------------------------------------------------------------------
GOF_linear_C = GOF_univariate(gam.fit = gam_fit, 
                            index = 4)
GOF_linear_C[[1]]

## -----------------------------------------------------------------------------
coefficients(gam_fit_nle)[1]
coefficients(gam_fit_nle)[2]
coefficients(gam_fit_nle)[3]
coefficients(gam_fit_nle)[4:12]

## -----------------------------------------------------------------------------
GOF_nle_DF = GOF_univariate(gam.fit = gam_fit_nle, 
                            index = 1)
GOF_nle_DF[[1]]

## -----------------------------------------------------------------------------
GOF_nle_DA = GOF_univariate(gam.fit = gam_fit_nle, 
                            index = 2)
GOF_nle_DA[[1]]

## -----------------------------------------------------------------------------
GOF_nle_C = GOF_univariate(gam.fit = gam_fit_nle, 
                            index = 3)
GOF_nle_C[[1]]

## -----------------------------------------------------------------------------
GOF_nle_IA = GOF_multivariate(gam.fit = gam_fit_nle, 
                  index = 4:12, BB.stat = BB9[[2]])
GOF_nle_IA[[1]]

## -----------------------------------------------------------------------------
GOF_nle_log_DF = GOF_univariate(gam.fit = gam_fit_nle_log, 
                            index = 1)
GOF_nle_log_DF[[1]]

## -----------------------------------------------------------------------------
GOF_nle_log_DA = GOF_univariate(gam.fit = gam_fit_nle_log, 
                            index = 2)
GOF_nle_log_DA[[1]]

## -----------------------------------------------------------------------------
GOF_nle_log_C = GOF_univariate(gam.fit = gam_fit_nle_log, 
                            index = 3)
GOF_nle_log_C[[1]]

## -----------------------------------------------------------------------------
GOF_nle_log_IA = GOF_multivariate(gam.fit = gam_fit_nle_log, 
                  index = 4:12, BB.stat = BB9[[2]])
GOF_nle_log_IA[[1]]

## -----------------------------------------------------------------------------
coefficients(gam_fit_tve)[1]
coefficients(gam_fit_tve)[2]
coefficients(gam_fit_tve)[3]
coefficients(gam_fit_tve)[4:13]

## -----------------------------------------------------------------------------
GOF_tve_DF = GOF_univariate(gam.fit = gam_fit_tve, 
                            index = 1)
GOF_tve_DF[[1]]

## -----------------------------------------------------------------------------
GOF_tve_DA = GOF_univariate(gam.fit = gam_fit_tve, 
                            index = 2)
GOF_tve_DA[[1]]

## -----------------------------------------------------------------------------
GOF_tve_C = GOF_univariate(gam.fit = gam_fit_tve, 
                            index = 3)
GOF_tve_C[[1]]

## -----------------------------------------------------------------------------
GOF_tve_IA = GOF_multivariate(gam.fit = gam_fit_tve, 
                  index = 4:13, BB.stat = BB10[[2]])
GOF_tve_IA[[1]]

## -----------------------------------------------------------------------------
plot(
  dat_gam_1$EVENT_INTERVAL,
  apply(GOF_tve_IA[[2]], 1, function(x) crossprod(x,x)),
  type="l", 
  lwd=2,
  col="black", 
  xlab="Time Points", 
  ylab="Squared Norm of the Martingale-Residual Process", 
  main="Goodness of Fit of Individual Activity with TVE",
  ylim=c(0,max(BB10[[1]]))
)
for (iter in 1:100){
  lines(
    dat_gam_1$EVENT_INTERVAL,
    BB10[[1]][seq(1, 2000, 
                  length.out=length(dat_gam_1$EVENT_INTERVAL)),iter], 
    col="darkgray", lwd=1, 
    lty=3)
}
lines(
  dat_gam_1$EVENT_INTERVAL,
  apply(GOF_tve_IA[[2]], 1, function(x) crossprod(x,x)),
  lwd=2,
  ylim= c(0,10), 
  col="black")
legend("topright",c("Observed", "BBridge"),
       lty = c(1,3),lwd = c(2,1), col = c("black","darkgray"))

