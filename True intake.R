###                               This code creates the true intake from the prescription data,                              ###
### then generates event times using the permutational algorithm with a prespecified hazard ratio for every exposure metric, ###
###                  and evaluates the discrepancies when using the observed exposure instead of the truth.                  ###
###                         This code is specific for Scenario AAAA, and the same for AAAB and AAAC                          ###
###                         Modifications for other scenarios are written as comments throughout                             ###
###                       This code is specifically for a sample size of 500 and null hazard ratio                           ###
###               It is also for a single replication, as this was parallelized to split the 500 replications                ###

rm(list = ls())
gc()

library(survival)

source("path/perm_algo_function.R")

exposure.dat = read.table("path/exposure_dat_Scenario_AAAA.txt", header = T)

n = 500

set.seed(seed)  

### SAMPLE INDIVIDUALS ###
  
sample.ids = sample(unique(exposure.dat$exposure_id), n, replace = T)

# For scenario BAAA, individuals are resampled such that the sample has, on average, smaller adherence
# mean.MPR_by.ID <- aggregate( exposure_MPR ~ exposure_id, exposure.dat[which(exposure.dat$exposure_last.RX == 0),], mean)  
# 
# sample.ids = sample(x = unique(mean.MPR_by.ID$exposure_id), size = n, replace = T, prob = 1.01 - mean.MPR_by.ID$exposure_MPR) 
  
sample.rows = c()
  
  for (i in c(1:length(sample.ids))) {
    id.rows = which(exposure.dat$exposure_id == sample.ids[i])
    sample.rows = c(sample.rows,id.rows)
  }
  
sample.dat = exposure.dat[sample.rows,]
  
sample.dat$exposure_id[which(as.numeric(row.names(sample.dat)) != floor(as.numeric(row.names(sample.dat))))] = 
    sample.dat$exposure_id[which(as.numeric(row.names(sample.dat)) != floor(as.numeric(row.names(sample.dat))))] + 
    (as.numeric(row.names(sample.dat)) %% 1)[which(as.numeric(row.names(sample.dat)) != floor(as.numeric(row.names(sample.dat))))]
  
sample.dat$exposure_id = as.numeric(factor(sample.dat$exposure_id))

### GENERATE DAY-BY-DAY OBSERVED EXPOSURE ###
  
exposure.dat_long = sample.dat[rep(seq_len(nrow(sample.dat)), sample.dat$exposure_lenght), ]
exposure.dat_long$f.up_day <- ave(1:nrow(exposure.dat_long), exposure.dat_long$exposure_id, FUN = seq_along)
exposure.dat_long$RX_day <- ave(1:nrow(exposure.dat_long), exposure.dat_long$exposure_id, exposure.dat_long$exposure_RX.episode, FUN = seq_along)
exposure.dat_long$last_RX_day = as.numeric(c(tail(exposure.dat_long$RX_day, -1), 1) == 1)
  
### GENERATE DAILY TRUE EXPOSURE ###
  
remaining_pills = 1
indices  = 1:nrow(exposure.dat_long)
  
  while (remaining_pills > 0 ) {
    
    remaining_pills = 0
    remaining_pills_id = c()
    remaining_pills_RX.episode = c()
    
    for (i in indices) {
      
      if(exposure.dat_long$RX_day[i] == 1){
        exposure.dat_long$true.prob[i] = ifelse(exposure.dat_long$exposure_MPR[i] < 0.9, 0.9, exposure.dat_long$exposure_MPR[i])
        exposure.dat_long$true.exposure[i] = rbinom(1,1,exposure.dat_long$true.prob[i])
        exposure.dat_long$cum_true.exposure[i] = exposure.dat_long$true.exposure[i]
      }
      else {
        exposure.dat_long$true.prob[i] = ifelse(exposure.dat_long$true.exposure[i-1] == 1,
                                                exposure.dat_long$exposure_MPR[i] + 0.5*(1-exposure.dat_long$exposure_MPR[i]), ### For scenarios ABAA and ACAA, 0.5 was replaced by 0.75 and 0.25, respectively.
                                                exposure.dat_long$exposure_MPR[i])
        exposure.dat_long$true.exposure[i] = ifelse(exposure.dat_long$cum_true.exposure[i-1] == exposure.dat_long$exposure_max.pills[i-1], 0,
                                                    rbinom(1,1,exposure.dat_long$true.prob[i]))
        exposure.dat_long$cum_true.exposure[i] = exposure.dat_long$true.exposure[i] + exposure.dat_long$cum_true.exposure[i-1]
      }
      
      if (exposure.dat_long$last_RX_day[i] == 1 & exposure.dat_long$cum_true.exposure[i] < exposure.dat_long$exposure_max.pills[i]){
        remaining_pills = remaining_pills + 1
        remaining_pills_id = c(remaining_pills_id, exposure.dat_long$exposure_id[i])
        remaining_pills_RX.episode = c(remaining_pills_RX.episode, exposure.dat_long$exposure_RX.episode[i])
      }
      
    }
    
    indices.rows = match(do.call("paste", exposure.dat_long[,c("exposure_id", "exposure_RX.episode")]),
                         do.call("paste", as.data.frame(cbind(remaining_pills_id, remaining_pills_RX.episode))))!="NA"
    indices = which(indices.rows==TRUE)
    
  }

# For Scenario AABA, the algorithm to generate true exposure allows individuals to stock medication before renewals
#
# indices  = 1:nrow(exposure.dat_long)
# 
# for (i in indices) {
#   
#   if(exposure.dat_long$RX_day[i] == 1){
#     exposure.dat_long$true.prob[i] = ifelse(exposure.dat_long$exposure_MPR[i] < 0.9, 0.9, exposure.dat_long$exposure_MPR[i])
#     exposure.dat_long$true.exposure[i] = rbinom(1,1,exposure.dat_long$true.prob[i])
#     exposure.dat_long$cum_true.exposure[i] = exposure.dat_long$true.exposure[i]
#   }
#   else {
#     exposure.dat_long$true.prob[i] = ifelse(exposure.dat_long$true.exposure[i-1] == 1, 
#                                             exposure.dat_long$exposure_MPR[i] + 0.5*(1-(exposure.dat_long$exposure_MPR[i])), 
#                                             exposure.dat_long$exposure_MPR[i])
#     exposure.dat_long$true.exposure[i] = ifelse(exposure.dat_long$cum_true.exposure[i-1] == exposure.dat_long$exposure_max.pills[i-1], 0,
#                                                 rbinom(1,1,exposure.dat_long$true.prob[i]))
#     exposure.dat_long$cum_true.exposure[i] = exposure.dat_long$true.exposure[i] + exposure.dat_long$cum_true.exposure[i-1]
#   }
#   
#   if (exposure.dat_long$last_RX_day[i] == 1 & exposure.dat_long$exposure_last.RX[i] != 1 & 
#       exposure.dat_long$cum_true.exposure[i] < exposure.dat_long$exposure_max.pills[i]){
#     
#     indice_id = exposure.dat_long$exposure_id[i]
#     indice_RX = exposure.dat_long$exposure_RX.episode[i+1]
#     indices.rows = match(do.call("paste", exposure.dat_long[,c("exposure_id", "exposure_RX.episode")]),
#                          do.call("paste", as.data.frame(cbind(indice_id, indice_RX))))!="NA"
#     indices = which(indices.rows==TRUE)
#     exposure.dat_long$exposure_max.pills[indices] = exposure.dat_long$exposure_max.pills[i+1] +  
#       (exposure.dat_long$exposure_max.pills[i] - exposure.dat_long$cum_true.exposure[i])
#   }
#   
# }

### CREATE DIFFERENT EXPOSURE METRICS ###
  
exposure.dat_long$cum_exposure = ave(exposure.dat_long$exposure, exposure.dat_long$exposure_id, FUN = cumsum) #From beginning
exposure.dat_long$true.cum_exposure = ave(exposure.dat_long$true.exposure, exposure.dat_long$exposure_id, FUN = cumsum) #From beginning
  
exposure.dat_long$cum14_exposure = c(rep(NA, 14), diff(exposure.dat_long$cum_exposure, lag=14))
exposure.dat_long$cum14_exposure[exposure.dat_long$f.up_day <= 14] <- 
    exposure.dat_long$cum_exposure[exposure.dat_long$f.up_day <= 14]
  
exposure.dat_long$true.cum14_exposure = c(rep(NA, 14), diff(exposure.dat_long$true.cum_exposure, lag=14))
exposure.dat_long$true.cum14_exposure[exposure.dat_long$f.up_day <= 14] <- 
    exposure.dat_long$true.cum_exposure[exposure.dat_long$f.up_day <= 14]
  
exposure.dat_long$any_exposure_14 = ifelse(exposure.dat_long$cum14_exposure == 0, 0, 1)
exposure.dat_long$true.any_exposure_14 = ifelse(exposure.dat_long$true.cum14_exposure == 0, 0, 1)
  
exposure.dat_long$cum30_exposure = c(rep(NA, 30), diff(exposure.dat_long$cum_exposure, lag=30))
exposure.dat_long$cum30_exposure[exposure.dat_long$f.up_day <= 30] <- 
    exposure.dat_long$cum_exposure[exposure.dat_long$f.up_day <= 30]
  
exposure.dat_long$true.cum30_exposure = c(rep(NA, 30), diff(exposure.dat_long$true.cum_exposure, lag=30))
exposure.dat_long$true.cum30_exposure[exposure.dat_long$f.up_day <= 30] <- 
    exposure.dat_long$true.cum_exposure[exposure.dat_long$f.up_day <= 30]
  
exposure.dat_long$cum60_exposure = c(rep(NA, 60), diff(exposure.dat_long$cum_exposure, lag=60))
exposure.dat_long$cum60_exposure[exposure.dat_long$f.up_day <= 60] <- 
    exposure.dat_long$cum_exposure[exposure.dat_long$f.up_day <= 60]
  
exposure.dat_long$true.cum60_exposure = c(rep(NA, 60), diff(exposure.dat_long$true.cum_exposure, lag=60))
exposure.dat_long$true.cum60_exposure[exposure.dat_long$f.up_day <= 60] <- 
    exposure.dat_long$true.cum_exposure[exposure.dat_long$f.up_day <= 60]
  
### PERMUTATIONAL ALGORITHM ###
  
true.beta = log(1)
true.beta_cum_30 = log(1)
true.beta_cum_60 = log(1)
  
death.times = t(read.table("path/death_times.txt"))
event.times = sort(ceiling(sample(death.times,n, replace = T)))
max_f.up = tapply(exposure.dat_long$f.up_day, exposure.dat_long$exposure_id, max)
  
final_data = perm_algo(seed, n, event.times, max_f.up, exposure.dat_long, true.beta)
n.events = sum(final_data$death.status)
  
final_data_any_14 = perm_algo_any_14(seed, n, event.times, max_f.up, exposure.dat_long, true.beta)
n.events_any_14 = sum(final_data_any_14$death.status)
  
final_data_cum_30 = perm_algo_cum_30(seed, n, event.times, max_f.up, exposure.dat_long, true.beta_cum_30)
n.events_cum_30 = sum(final_data_cum_30$death.status)
  
final_data_cum_60 = perm_algo_cum_60(seed, n, event.times, max_f.up, exposure.dat_long, true.beta_cum_60)
n.events_cum_60 = sum(final_data_cum_60$death.status)

### RESULTS ###
  
cox.reg.true = coxph(Surv(start, f.up_day, death.status) ~ true.exposure, data = final_data)
coef.true = exp(cox.reg.true$coefficients)
se.true = sqrt(cox.reg.true$var)
coverage.true = confint(cox.reg.true)[1,1] <= true.beta & true.beta <= confint(cox.reg.true)[1,2]
power.true = confint(cox.reg.true)[1,1] <= 0 & 0 <= confint(cox.reg.true)[1,2]

cox.reg.observed = coxph(Surv(start, f.up_day, death.status) ~ exposure, data = final_data)
coef.observed = exp(cox.reg.observed$coefficients)
se.observed = sqrt(cox.reg.observed$var)
coverage.observed =  confint(cox.reg.observed)[1,1] <= true.beta & true.beta <= confint(cox.reg.observed)[1,2]
power.observed = confint(cox.reg.observed)[1,1] <= 0 & 0 <= confint(cox.reg.observed)[1,2]

cox.reg.true_any_14 = coxph(Surv(start, f.up_day, death.status) ~ true.any_exposure_14, data = final_data_any_14)
coef.true_any_14 = exp(cox.reg.true_any_14$coefficients)
se.true_any_14 = sqrt(cox.reg.true_any_14$var)
coverage.true_any_14 = confint(cox.reg.true_any_14)[1,1] <= true.beta & true.beta <= confint(cox.reg.true_any_14)[1,2]
power.true_any_14 = confint(cox.reg.true_any_14)[1,1] <= 0 & 0 <= confint(cox.reg.true_any_14)[1,2]

cox.reg.observed_any_14 = coxph(Surv(start, f.up_day, death.status) ~ any_exposure_14, data = final_data_any_14)
coef.observed_any_14 = exp(cox.reg.observed_any_14$coefficients)
se.observed_any_14 = sqrt(cox.reg.observed_any_14$var)
coverage.observed_any_14 =  confint(cox.reg.observed_any_14)[1,1] <= true.beta & true.beta <= confint(cox.reg.observed_any_14)[1,2]
power.observed_any_14 = confint(cox.reg.observed_any_14)[1,1] <= 0 & 0 <= confint(cox.reg.observed_any_14)[1,2]

cox.reg.true_cum_30 = coxph(Surv(start, f.up_day, death.status) ~ true.cum30_exposure, final_data_cum_30)
coef.true_cum_30 = exp(cox.reg.true_cum_30$coefficients)^30
se.true_cum_30 = sqrt(cox.reg.true_cum_30$var)
coverage.true_cum_30 = confint(cox.reg.true_cum_30)[1,1] <= true.beta_cum_30 & true.beta_cum_30 <= confint(cox.reg.true_cum_30)[1,2]
power.true_cum_30 = confint(cox.reg.true_cum_30)[1,1] <= 0 & 0 <= confint(cox.reg.true_cum_30)[1,2]

cox.reg.observed_cum_30 = coxph(Surv(start, f.up_day, death.status) ~ cum30_exposure, final_data_cum_30)
coef.observed_cum_30 = exp(cox.reg.observed_cum_30$coefficients)^30
se.observed_cum_30 = sqrt(cox.reg.observed_cum_30$var)
coverage.observed_cum_30 = confint(cox.reg.observed_cum_30)[1,1] <= true.beta_cum_30 & true.beta_cum_30 <= confint(cox.reg.observed_cum_30)[1,2]
power.observed_cum_30 = confint(cox.reg.observed_cum_30)[1,1] <= 0 & 0 <= confint(cox.reg.observed_cum_30)[1,2]

cox.reg.true_cum_60 = coxph(Surv(start, f.up_day, death.status) ~ true.cum60_exposure, final_data_cum_60)
coef.true_cum_60 = exp(cox.reg.true_cum_60$coefficients)^60
se.true_cum_60 = sqrt(cox.reg.true_cum_60$var)
coverage.true_cum_60 = confint(cox.reg.true_cum_60)[1,1] <= true.beta_cum_60 & true.beta_cum_60 <= confint(cox.reg.true_cum_60)[1,2]
power.true_cum_60 = confint(cox.reg.true_cum_60)[1,1] <= 0 & 0 <= confint(cox.reg.true_cum_60)[1,2]

cox.reg.observed_cum_60 = coxph(Surv(start, f.up_day, death.status) ~ cum60_exposure, final_data_cum_60)
coef.observed_cum_60 = exp(cox.reg.observed_cum_60$coefficients)^60
se.observed_cum_60 = sqrt(cox.reg.observed_cum_60$var)
coverage.observed_cum_60 = confint(cox.reg.observed_cum_60)[1,1] <= true.beta_cum_60 & true.beta_cum_60 <= confint(cox.reg.observed_cum_60)[1,2]
power.observed_cum_60 = confint(cox.reg.observed_cum_60)[1,1] <= 0 & 0 <= confint(cox.reg.observed_cum_60)[1,2]
