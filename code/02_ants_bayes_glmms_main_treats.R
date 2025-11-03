### wog-wog-postfire-ants

### 02_ants_bayes_glmms_main_treats.R

### Matt Bitters
### matthew.bitters@colorado.edu








# ============================================================
#  0. Setup
# ============================================================

### Install and load required packages
packages <- c("here", "dplyr", "tidyr", "brms", "rstan")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(dplyr)
library(tidyr)
library(brms)
library(rstan)


# Stan options
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



# ============================================================
#  1. Read in data and factor
# ============================================================


# Read CSV
ant_mod_data <- read.csv(here("data", "derived", "ants_model_ready_10-31-2025.csv"))

# Examine data
str(ant_mod_data)
head(ant_mod_data)

# Factor and set levels
# ant_mod_data$bio_year <- factor(ant_mod_data$bio_year)            # 1-35 (with some missing years)
ant_mod_data$treat <- factor(ant_mod_data$treat, levels=c(2,1,3))   # 2=controls, 1=fragments, 3=matrix
ant_mod_data$size <- factor(ant_mod_data$size,levels=c(4,1,2,3))    # 4=controls, 1=small, 2=medium, 3=large, NA=matrix
ant_mod_data$edge <- factor(ant_mod_data$edge, levels=c(3,1,2))     # 3=controls, 1=core, 2=edge, NA=matrix
ant_mod_data$site <- factor(ant_mod_data$site)                      # 1-72=fragments, 73-120=controls, 121-144=fragments, 145-188=matrix
ant_mod_data$rep <- factor(ant_mod_data$rep)                        # 1-6
ant_mod_data$patch <- factor(ant_mod_data$patch)                    # 1-30
ant_mod_data$topo <- factor(ant_mod_data$topo)                      # 1=drain, 2=slope
ant_mod_data$time_period <- factor(ant_mod_data$time_period, levels=c("post-frag", "post-fire"))

# Scale yrs1_2_abund for modelling
ant_mod_data$yrs1_2_abund_log_s <- scale(ant_mod_data$yrs1_2_abund_log)
ant_mod_data$yrs1_2_abund_log_s <- as.numeric(ant_mod_data$yrs1_2_abund_log_s)

# Confirm factoring
str(ant_mod_data)




# Filter to species
aphan_data <- ant_mod_data %>%
  filter(species == "aphaenogaster_longiceps")

lepto_data <- ant_mod_data %>%
  filter(species == "leptomyrmex_erythrocephalus")



# ============================================================
#  2. Priors
# ============================================================

priors_null <- c(
  prior(normal(0, 2), class = "Intercept"),  # intercept
  prior(student_t(3, 0, 2), class = "sd")   # random effect SDs
)

priors <- c(
  prior(normal(0, 1), class = "b"),         # slopes
  prior(normal(0, 2), class = "Intercept"), # intercept
  prior(student_t(3, 0, 2), class = "sd")   # random effect SDs
)





# ============================================================
#  3. Aphan treatment models
# ============================================================


m0_null_aphan <- brm(
  presence ~ 1 + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors_null,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m1_aphan <- brm(
  presence ~ treat + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m2_aphan <- brm(
  presence ~ treat + time_period + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)


m3_aphan <- brm(
  presence ~ treat + time_period + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m4_aphan <- brm(
  presence ~ treat*time_period + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)



### Compare initial models

# Look at summaries
summary(m0_null_aphan)
summary(m1_aphan)
summary(m2_aphan)
summary(m3_aphan)
summary(m4_aphan)

# Compare loo
loo_compare(loo(m0_null_aphan), 
            loo(m1_aphan), 
            loo(m2_aphan), 
            loo(m3_aphan),
            loo(m4_aphan)
)

# Compare marginal R^2s
bayes_R2(m0_null_aphan)
bayes_R2(m1_aphan)
bayes_R2(m2_aphan)
bayes_R2(m3_aphan)
bayes_R2(m4_aphan)




# Bio_year*treatment interaction for time series plot
m5_years_aphan <- brm(
  presence ~ treat*bio_year + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = aphan_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(m5_years_aphan)
bayes_R2(m5_years_aphan)



### Save  models for model comparison table in paper

saveRDS(m0_null_aphan, here("model_outputs", "m0_null_aphan.rds"))
saveRDS(m1_aphan, here("model_outputs", "m1_aphan.rds"))
saveRDS(m2_aphan, here("model_outputs", "m2_aphan.rds"))
saveRDS(m3_aphan, here("model_outputs", "m3_aphan.rds"))
saveRDS(m4_aphan, here("model_outputs", "m4_aphan.rds"))
saveRDS(m5_years_aphan, here("model_outputs", "m5_years_aphan.rds"))



# ============================================================
#  4. Lepto treatment models
# ============================================================



m0_null_lepto <- brm(
  presence ~ 1 + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors_null,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m1_lepto <- brm(
  presence ~ treat + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m2_lepto <- brm(
  presence ~ treat + time_period + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)


m3_lepto <- brm(
  presence ~ treat + time_period + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 6000,
  warmup = 3000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

m4_lepto <- brm(
  presence ~ treat*time_period + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)



### Compare initial models

# Look at summaries
summary(m0_null_lepto)
summary(m1_lepto)
summary(m2_lepto)
summary(m3_lepto)
summary(m4_lepto)

# Compare loo
loo_compare(loo(m0_null_lepto), 
            loo(m1_lepto), 
            loo(m2_lepto), 
            loo(m3_lepto),
            loo(m4_lepto)
)

# Compare marginal R^2s
bayes_R2(m0_null_lepto)
bayes_R2(m1_lepto)
bayes_R2(m2_lepto)
bayes_R2(m3_lepto)
bayes_R2(m4_lepto)




# Bio_year*treatment interaction for time series plot
m5_years_lepto <- brm(
  presence ~ treat*bio_year + yrs1_2_abund_log_s + (1 | patch) + (1 | rep),
  data = lepto_data,
  family = bernoulli(link = "logit"),
  prior = priors,
  cores = 4, chains = 4,
  iter = 10000,
  warmup = 5000,
  control = list(adapt_delta = 0.99, max_treedepth = 15)
)

summary(m5_years_lepto)
bayes_R2(m5_years_lepto)



### Save  models for model comparison table in paper

saveRDS(m0_null_lepto, here("model_outputs", "m0_null_lepto.rds"))
saveRDS(m1_lepto, here("model_outputs", "m1_lepto.rds"))
saveRDS(m2_lepto, here("model_outputs", "m2_lepto.rds"))
saveRDS(m3_lepto, here("model_outputs", "m3_lepto.rds"))
saveRDS(m4_lepto, here("model_outputs", "m4_lepto.rds"))
saveRDS(m5_years_lepto, here("model_outputs", "m5_years_lepto.rds"))










































#For models without matrix, we want 7 categories (C, SI, SO, MI, MO, LI, LO)
c  <- ifelse(               Ants$treat==2,              1, 0 )
si <- ifelse( Ants$size ==1 & Ants$treat==1 & Ants$edge==1, 2, 0 )
so <- ifelse( Ants$size==1 & Ants$treat==1 & Ants$edge==2, 3, 0 )
mi <- ifelse( Ants$size==2 & Ants$treat==1 & Ants$edge==1, 4, 0 )
mo <- ifelse( Ants$size==2 & Ants$treat==1 & Ants$edge==2, 5, 0 )
li <- ifelse( Ants$size==3 & Ants$treat==1 & Ants$edge==1, 6, 0 )
lo <- ifelse( Ants$size==3 & Ants$treat==1 & Ants$edge==2, 7, 0 )
Ants$treat_size_edge <- c + li + lo + mi + mo + si + so
rm(c,li,lo,mi,mo,si,so)
Ants$treat_size_edge <- factor(Ants$treat_size_edge)



