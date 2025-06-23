## ----setup, include=FALSE-----------------------------------------------------
knitr::knit_hooks$set(purl = knitr::hook_purl)
knitr::opts_chunk$set(echo = TRUE)

## ----message=FALSE, warning=FALSE---------------------------------------------
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
data_original <- read.csv("_INPUT_CP_2_/jean_events_EVENTS.csv")

## -----------------------------------------------------------------------------
head(data_original)

## -----------------------------------------------------------------------------
data <- data_original[,setdiff(colnames(data_original), 
                               c("TARGET", "TYPE", "EVENT_INTERVAL_ID", "EVENT", 
                                 "INTEGER_TIME", "TIME_POINT", 
                                 "TIME_UNIT", "num.actors"))]
rm(data_original)
head(data)

## -----------------------------------------------------------------------------
m = 20

## -----------------------------------------------------------------------------
data_ev <- data %>% filter(IS_OBSERVED == 1)
data_nv <- data %>% filter(IS_OBSERVED == 0)

## -----------------------------------------------------------------------------
nrow(data_ev)

## -----------------------------------------------------------------------------
nrow(data_nv)
nrow(data_nv) == nrow(data_ev) * m

## -----------------------------------------------------------------------------
coxph_fit <- coxph(Surv(time = rep(1,nrow(data)), 
                        event = data$IS_OBSERVED) ~ 
                         + diff.female
                         # + female
                         + individual.activity 
                         + dyadic.activity 
                         + closure
                         + strata(EVENT_INTERVAL)
                         , data = data)
summary(coxph_fit)

## -----------------------------------------------------------------------------
clogit_fit <- clogit(IS_OBSERVED ~ 
                         + diff.female
                         # + female
                         +  individual.activity 
                         + dyadic.activity 
                         + closure
                         + strata(EVENT_INTERVAL)
                         , data = data)
summary(clogit_fit)

## -----------------------------------------------------------------------------
coefficients(coxph_fit)
coefficients(clogit_fit)

## -----------------------------------------------------------------------------
merge_id_cols <- c("EVENT_INTERVAL")
# Consider 4 different values of m
ms <- c(1, 5, 10, 15)

## -----------------------------------------------------------------------------
# Add a new column to label each row in 'data_ev' and 'data_nv' 
# as event data and non-event data respectively
data_ev_tagged <- data_ev %>%
  mutate(.row_type = "ev")
data_nv_tagged <- data_nv %>%
  mutate(.row_type = "nv")
head(data_ev_tagged)
# head(data_nv_tagged)
rm(data_ev)
rm(data_nv)

## -----------------------------------------------------------------------------
# List to store the sampled datasets for each sampling size
dat_sampled <- list()

# Loop through each sampling size
for (mm in ms){
  
  set.seed(mm) 

  # Sample 'm' rows (if available) from each group in 'data_nv'
  data_nv_m_sampled <- data_nv_tagged %>%
    group_by(across(all_of(merge_id_cols))) %>%
    slice_sample(n = mm, replace = FALSE) %>%
    ungroup()

  # Combine the event data and the sampled non-event data
  dat_sampled_m <- bind_rows(data_ev_tagged, data_nv_m_sampled) %>%
    
    # Sort the combined data by the grouping columns and by row type 
    # (ev first, then nv)
    arrange(across(all_of(merge_id_cols)), .row_type)

  # Store the resulting data frame in the list
  dat_sampled[[match(mm, ms)]] <- dat_sampled_m
}

# Assign names to each list element
names(dat_sampled) <- ms

# Remove intermediate variables to clean up the workspace
rm(data_nv_m_sampled, dat_sampled_m)

## -----------------------------------------------------------------------------
# let's check
nrow(dat_sampled[[1]]) == nrow(data_ev_tagged)*(ms[1]+1)
nrow(dat_sampled[[2]]) == nrow(data_ev_tagged)*(ms[2]+1)
nrow(dat_sampled[[3]]) == nrow(data_ev_tagged)*(ms[3]+1)
nrow(dat_sampled[[4]]) == nrow(data_ev_tagged)*(ms[4]+1)

## -----------------------------------------------------------------------------
head(dat_sampled[[2]], 10)

## -----------------------------------------------------------------------------
clogit_fits <- list()

## -----------------------------------------------------------------------------
for (mm in ms){
  print(mm)
  clogit_fits[[match(mm, ms)]] <- clogit(IS_OBSERVED ~ 
                         + diff.female
                         # + female
                         + individual.activity 
                         + dyadic.activity 
                         + closure
                         + strata(EVENT_INTERVAL)
                         , data = dat_sampled[[match(mm, ms)]])
}

## -----------------------------------------------------------------------------
# summary(clogit_fits[[1]])

## -----------------------------------------------------------------------------
summary(clogit_fits[[2]])

## -----------------------------------------------------------------------------
summary(clogit_fits[[3]])

## -----------------------------------------------------------------------------
summary(clogit_fits[[4]])

## -----------------------------------------------------------------------------
set.seed(10)

# For each group defined by the merge_id_cols, take one random non-event
data_nv_1_sampled <- data_nv_tagged %>%
  group_by(across(all_of(merge_id_cols))) %>%
  slice_sample(n = 1) %>%
  ungroup()

# Perform a left join: 
# for each row in data_ev, attach the corresponding non-event row (if available)
# based on matching values in the merge_id_cols
# The suffixes "_ev" and "_nv" will be added to columns from data_ev and data_nv 
dat_gam_1 <- data_ev_tagged %>%
  left_join(data_nv_1_sampled, 
            by = merge_id_cols, suffix = c("_ev", "_nv"))

# All responses are set equal to 1
dat_gam_1$y <- 1

rm(data_nv_1_sampled)

## -----------------------------------------------------------------------------
# covariates defined as difference between 
# covariate values for event and corresponding non-event
dat_gam_1$female <- 
  dat_gam_1$female_ev - dat_gam_1$female_nv
dat_gam_1$diff_female <- 
  dat_gam_1$diff.female_ev - dat_gam_1$diff.female_nv
dat_gam_1$individual_activity <- 
  dat_gam_1$individual.activity_ev - dat_gam_1$individual.activity_nv
dat_gam_1$dyadic_activity <- 
  dat_gam_1$dyadic.activity_ev - dat_gam_1$dyadic.activity_nv
dat_gam_1$closure <- 
  dat_gam_1$closure_ev - dat_gam_1$closure_nv

## -----------------------------------------------------------------------------
gam_fit <- glm(y ~ 
                + diff_female
                # + female
                + individual_activity 
                + dyadic_activity 
                + closure 
                - 1 # no intercept
               , data = dat_gam_1, 
               family="binomial")
summary(gam_fit)

## -----------------------------------------------------------------------------
save(dat_gam_1, file="_OUTPUT_CP_2/dat_gam_1.RData")
save(gam_fit, file="_OUTPUT_CP_2/gam_fit.RData")

## -----------------------------------------------------------------------------
# Store the case-control datasets for different m
dat_gam <- list()

# Loop through each sampling size defined in 'ms'
for (mm in ms){
  
  set.seed(mm)

  # Sample randomly 'm' rows from each group in 'data_nv'
  data_nv_m_sampled <- data_nv_tagged %>%
    group_by(across(all_of(merge_id_cols))) %>%
    slice_sample(n = mm, replace = FALSE) %>%
    ungroup()

  # Left join the sampled non-event data to the event data using the merge keys
  # this combines each row of event data with 'm' non-event rows (if available)
  dat_gam_m <- data_ev_tagged %>%
    left_join(data_nv_m_sampled, by = merge_id_cols, suffix = c("_ev", "_nv"))

  # Compute covariate difference
  dat_gam_m$female <- 
    dat_gam_1$female_ev - dat_gam_1$female_nv
  dat_gam_m$diff_female <- 
    dat_gam_1$diff.female_ev - dat_gam_1$diff.female_nv
  dat_gam_m$individual_activity <- 
    dat_gam_m$individual.activity_ev - dat_gam_m$individual.activity_nv
  dat_gam_m$dyadic_activity <- 
    dat_gam_m$dyadic.activity_ev - dat_gam_m$dyadic.activity_nv
  dat_gam_m$closure <- 
    dat_gam_m$closure_ev - dat_gam_m$closure_nv
  
  # All responses are set equal to 1
  dat_gam_m$y <- 1
  
  # Store the resulting data frame in the list
  dat_gam[[match(mm, ms)]] <- dat_gam_m
}

# Assign names to each list element
names(dat_gam) <- ms

# Remove intermediate variables to clean up the workspace
rm(data_nv_m_sampled, dat_gam_m)

## -----------------------------------------------------------------------------
glm_fits <- list()

## ----message=FALSE, warning=FALSE---------------------------------------------
for (mm in ms){
  glm_fits[[match(mm, ms)]] <- glm(y ~ 
                                    + diff_female
                                    # + female
                                    + individual_activity 
                                    + dyadic_activity 
                                    + closure 
                                    - 1 ,# no intercept
                         family = "binomial",
                         data = dat_gam[[match(mm, ms)]])
}

## -----------------------------------------------------------------------------
summary(glm_fits[[1]])

## -----------------------------------------------------------------------------
summary(glm_fits[[2]])

## -----------------------------------------------------------------------------
summary(glm_fits[[3]])

## -----------------------------------------------------------------------------
summary(glm_fits[[4]])

