### This code implements the permutational algorithm to generate event times corresponding to a prespecified hazard ratio ###
###         The code contains four functions, one for each exposure metric investigated in the accompanying paper         ###
### The functions are not optimized and require a dataset with columns named "exposure_id", "f_up.day", "true.exposure"  ###

perm_algo <- function(seed, n, event_times, censoring_times, exposure_data, true_beta){
  
  set.seed(seed)  
  
  # CREATING A DATASET CONTAINING THE EXPOSURE VALUE FOR EACH INDIVIDUAL AT EACH EVENT TIME
  # IF THE INDIVIDUAL IS CENSORED BEFORE THE EVENT TIME THEN THE EXPOSURE IS NA
  
  event_times_by.id = as.data.frame(cbind(rep(unique(exposure_data$exposure_id), each = n), rep(event_times, times = n)))
  colnames(event_times_by.id) = c("exposure_id","f.up_day")
  exposure.at.event_times = exposure_data[which(exposure_data$f.up_day %in% event_times), c("exposure_id", "f.up_day","true.exposure")]
  
  exposure_all.event_times = merge(event_times_by.id, exposure.at.event_times, all = T)

  exposure.mat = matrix(exposure_all.event_times$true.exposure,n,n)
  colnames(exposure.mat) = names(censoring_times)
  rownames(exposure.mat) = event_times
  
  # CREATING A DATASET CONTAINING THE CENSORING TIMES TO BE REPLACED BY EVENT TIMES IF SAMPLED
  
  max.cens.event_by.id = as.data.frame(cbind(as.numeric(names(censoring_times)),censoring_times, rep(0,n)), stringsAsFactors = F)
  colnames(max.cens.event_by.id) = c("exposure_id", "censoring_times", "death")

  for (i in c(1:n)) {
    
    if(all(is.na(exposure.mat[i,]))){
      break()
    }
    else{
      # PARTIAL LIKELIHOOD
      p = exp(exposure.mat[i,]*true_beta)/sum(exp(exposure.mat[i,]*true_beta), na.rm = T)
      event_id = sample(colnames(exposure.mat)[!is.na(p)], size = 1, replace = FALSE, prob = p[!is.na(p)])
      
      #REPLACING THE CENSORING TIME BY THE EVENT TIME
      max.cens.event_by.id[which(max.cens.event_by.id$exposure_id == event_id),] = as.numeric(c(event_id, row.names(exposure.mat)[i], 1))
      
      #SETTING ALL FOLLOWING TIMES TO NA FOR THE INDIVIDUAL THAT HAD THE EVENT
      exposure.mat[c(i:n),event_id] = NA
    }
    
  }
  
  exposure_data_death = merge(exposure_data, max.cens.event_by.id, by = "exposure_id")
  exposure_data_death = exposure_data_death[which(exposure_data_death$f.up_day <= as.numeric(exposure_data_death$censoring_times)),]

  exposure_data_death$death.status = ifelse(exposure_data_death$death == 1 & c(tail(exposure_data_death$f.up_day, -1), 1) == 1, 1, 0)

  exposure_data_death$start = exposure_data_death$f.up_day - 1
  
  final_data = exposure_data_death[,c("exposure_id", "start", "f.up_day", "true.exposure", "exposure", "death.status")]
  
  return(final_data)
  
}

perm_algo_any_14 <- function(seed, n, event_times, censoring_times, exposure_data, true_beta){
  
  set.seed(seed)  
  
  # CREATING A DATASET CONTAINING THE EXPOSURE VALUE FOR EACH INDIVIDUAL AT EACH EVENT TIME
  # IF THE INDIVIDUAL IS CENSORED BEFORE THE EVENT TIME THEN THE EXPOSURE IS NA
  
  event_times_by.id = as.data.frame(cbind(rep(unique(exposure_data$exposure_id), each = n), rep(event_times, times = n)))
  colnames(event_times_by.id) = c("exposure_id","f.up_day")
  exposure.at.event_times = exposure_data[which(exposure_data$f.up_day %in% event_times), c("exposure_id", "f.up_day","true.any_exposure_14")]
  
  exposure_all.event_times = merge(event_times_by.id, exposure.at.event_times, all = T)

  exposure.mat = matrix(exposure_all.event_times$true.any_exposure_14,n,n)
  colnames(exposure.mat) = names(censoring_times)
  rownames(exposure.mat) = event_times
  
  # CREATING A DATASET CONTAINING THE CENSORING TIMES TO BE REPLACED BY EVENT TIMES IF SAMPLED
  
  max.cens.event_by.id = as.data.frame(cbind(as.numeric(names(censoring_times)),censoring_times, rep(0,n)), stringsAsFactors = F)
  colnames(max.cens.event_by.id) = c("exposure_id", "censoring_times", "death")

  for (i in c(1:n)) {
    
    if(all(is.na(exposure.mat[i,]))){
      break()
    }
    else{
      # PARTIAL LIKELIHOOD
      p = exp(exposure.mat[i,]*true_beta)/sum(exp(exposure.mat[i,]*true_beta), na.rm = T)
      event_id = sample(colnames(exposure.mat)[!is.na(p)], size = 1, replace = FALSE, prob = p[!is.na(p)])
      
      #REPLACING THE CENSORING TIME BY THE EVENT TIME
      max.cens.event_by.id[which(max.cens.event_by.id$exposure_id == event_id),] = as.numeric(c(event_id, row.names(exposure.mat)[i], 1))
      
      #SETTING ALL FOLLOWING TIMES TO NA FOR THE INDIVIDUAL THAT HAD THE EVENT
      exposure.mat[c(i:n),event_id] = NA
    }
    
  }
  
  exposure_data_death = merge(exposure_data, max.cens.event_by.id, by = "exposure_id")
  exposure_data_death = exposure_data_death[which(exposure_data_death$f.up_day <= as.numeric(exposure_data_death$censoring_times)),]

  exposure_data_death$death.status = ifelse(exposure_data_death$death == 1 & c(tail(exposure_data_death$f.up_day, -1), 1) == 1, 1, 0)

  exposure_data_death$start = exposure_data_death$f.up_day - 1
  
  final_data = exposure_data_death[,c("exposure_id", "start", "f.up_day", "true.any_exposure_14", "any_exposure_14", "death.status")]
  
  return(final_data)
  
}

perm_algo_cum_60 <- function(seed, n, event_times, censoring_times, exposure_data, true_beta){
  
  set.seed(seed)  
  
  # CREATING A DATASET CONTAINING THE EXPOSURE VALUE FOR EACH INDIVIDUAL AT EACH EVENT TIME
  # IF THE INDIVIDUAL IS CENSORED BEFORE THE EVENT TIME THEN THE EXPOSURE IS NA
  
  event_times_by.id = as.data.frame(cbind(rep(unique(exposure_data$exposure_id), each = n), rep(event_times, times = n)))
  colnames(event_times_by.id) = c("exposure_id","f.up_day")
  exposure.at.event_times = exposure_data[which(exposure_data$f.up_day %in% event_times), c("exposure_id", "f.up_day","true.cum60_exposure")]
  
  exposure_all.event_times = merge(event_times_by.id, exposure.at.event_times, all = T)

  exposure.mat = matrix(exposure_all.event_times$true.cum60_exposure,n,n)
  colnames(exposure.mat) = names(censoring_times)
  rownames(exposure.mat) = event_times
  
  # CREATING A DATASET CONTAINING THE CENSORING TIMES TO BE REPLACED BY EVENT TIMES IF SAMPLED
  
  max.cens.event_by.id = as.data.frame(cbind(as.numeric(names(censoring_times)),censoring_times, rep(0,n)), stringsAsFactors = F)
  colnames(max.cens.event_by.id) = c("exposure_id", "censoring_times", "death")

  for (i in c(1:n)) {
    
    if(all(is.na(exposure.mat[i,]))){
      break()
    }
    else{
      # PARTIAL LIKELIHOOD
      p = exp(exposure.mat[i,]*true_beta)/sum(exp(exposure.mat[i,]*true_beta), na.rm = T)
      event_id = sample(colnames(exposure.mat)[!is.na(p)], size = 1, replace = FALSE, prob = p[!is.na(p)])
      
      #REPLACING THE CENSORING TIME BY THE EVENT TIME
      max.cens.event_by.id[which(max.cens.event_by.id$exposure_id == event_id),] = as.numeric(c(event_id, row.names(exposure.mat)[i], 1))
      
      #SETTING ALL FOLLOWING TIMES TO NA FOR THE INDIVIDUAL THAT HAD THE EVENT
      exposure.mat[c(i:n),event_id] = NA
    }
    
  }
  
  exposure_data_death = merge(exposure_data, max.cens.event_by.id, by = "exposure_id")
  exposure_data_death = exposure_data_death[which(exposure_data_death$f.up_day <= as.numeric(exposure_data_death$censoring_times)),]

  exposure_data_death$death.status = ifelse(exposure_data_death$death == 1 & c(tail(exposure_data_death$f.up_day, -1), 1) == 1, 1, 0)

  exposure_data_death$start = exposure_data_death$f.up_day - 1
  
  final_data = exposure_data_death[,c("exposure_id", "start", "f.up_day", "true.cum60_exposure", "cum60_exposure", "death.status")]
  
  return(final_data)
  
}

perm_algo_cum_30 <- function(seed, n, event_times, censoring_times, exposure_data, true_beta){
  
  set.seed(seed)  
  
  # CREATING A DATASET CONTAINING THE EXPOSURE VALUE FOR EACH INDIVIDUAL AT EACH EVENT TIME
  # IF THE INDIVIDUAL IS CENSORED BEFORE THE EVENT TIME THEN THE EXPOSURE IS NA
  
  event_times_by.id = as.data.frame(cbind(rep(unique(exposure_data$exposure_id), each = n), rep(event_times, times = n)))
  colnames(event_times_by.id) = c("exposure_id","f.up_day")
  exposure.at.event_times = exposure_data[which(exposure_data$f.up_day %in% event_times), c("exposure_id", "f.up_day","true.cum30_exposure")]
  
  exposure_all.event_times = merge(event_times_by.id, exposure.at.event_times, all = T)

  exposure.mat = matrix(exposure_all.event_times$true.cum30_exposure,n,n)
  colnames(exposure.mat) = names(censoring_times)
  rownames(exposure.mat) = event_times
  
  # CREATING A DATASET CONTAINING THE CENSORING TIMES TO BE REPLACED BY EVENT TIMES IF SAMPLED
  
  max.cens.event_by.id = as.data.frame(cbind(as.numeric(names(censoring_times)),censoring_times, rep(0,n)), stringsAsFactors = F)
  colnames(max.cens.event_by.id) = c("exposure_id", "censoring_times", "death")

  for (i in c(1:n)) {
    
    if(all(is.na(exposure.mat[i,]))){
      break()
    }
    else{
      # PARTIAL LIKELIHOOD
      p = exp(exposure.mat[i,]*true_beta)/sum(exp(exposure.mat[i,]*true_beta), na.rm = T)
      event_id = sample(colnames(exposure.mat)[!is.na(p)], size = 1, replace = FALSE, prob = p[!is.na(p)])
      
      #REPLACING THE CENSORING TIME BY THE EVENT TIME
      max.cens.event_by.id[which(max.cens.event_by.id$exposure_id == event_id),] = as.numeric(c(event_id, row.names(exposure.mat)[i], 1))
      
      #SETTING ALL FOLLOWING TIMES TO NA FOR THE INDIVIDUAL THAT HAD THE EVENT
      exposure.mat[c(i:n),event_id] = NA
    }
    
  }
  
  exposure_data_death = merge(exposure_data, max.cens.event_by.id, by = "exposure_id")
  exposure_data_death = exposure_data_death[which(exposure_data_death$f.up_day <= as.numeric(exposure_data_death$censoring_times)),]

  exposure_data_death$death.status = ifelse(exposure_data_death$death == 1 & c(tail(exposure_data_death$f.up_day, -1), 1) == 1, 1, 0)

  exposure_data_death$start = exposure_data_death$f.up_day - 1
  
  final_data = exposure_data_death[,c("exposure_id", "start", "f.up_day", "true.cum30_exposure", "cum30_exposure", "death.status")]
  
  return(final_data)
  
}