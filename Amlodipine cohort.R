### THIS CODE EXTRACTS AMLODIPINE PRESCRIPTIONS FROM THE CNODES PRESCRIPTION DATASET ###
###  IT ALSO GENERATES THE FINAL COHORT FROM WHICH THE INDIVIDUALS WILL BE RESAMPLED ###

rm(list = ls())
gc()

library(haven)
CNODES_drug_data = read_sas("path/sd_table_drug.sas7bdat")

# Restrict to Amlodipine (C08CA01)
drug_data = CNODES_drug_data[drug_data$atc == "C08CA01",]

#Create start and stop date
drug_data$date = as.Date(drug_data$date)
drug_data$stop.date = as.Date(drug_data$date + drug_data$dur - 1)

#### DATA CLEANING ####

#Removing incoherent drug episodes

drug_data$lag.start.date = c(as.Date(NA), head(drug_data$date, - 1)) ; drug_data$lag.start.date[!duplicated(drug_data$id)] = as.Date(NA)                   
drug_data$new.episode = ifelse(drug_data$date < drug_data$lag.start.date, 1, 0)
drug_data$new.episode[!duplicated(drug_data$id)] = 0
drug_data$RX_in_new.episode = ave(drug_data$new.episode, drug_data$id, FUN=cumsum)
drug_data_single.episode = drug_data[drug_data$RX_in_new.episode==0,]

#Shifting the gap distribution by 30 days

drug_data_single.episode$RX.count <- ave(1:nrow(drug_data_single.episode), drug_data_single.episode$id, FUN = seq_along) - 1
drug_data_single.episode$start.date = as.Date(drug_data_single.episode$date - 31*drug_data_single.episode$RX.count)

##################### FINAL DATASET #####################

C08CA01_shift30_single.episode = drug_data_single.episode[,c("id", "atc", "start.date", "dur")]
write.table(C08CA01_shift30_single.episode, "path/drug_data_C08CA01_shift30_single.episode.txt", row.names = F)


################### WITH INCLUSION CRITERIA

rm(list = ls())
gc()

drug_data = read.table("path/drug_data_C08CA01_shift30_single.episode.txt", header = T)

library(haven)
demo.death_data = as.data.frame(read_sas("path/sd_table_death.sas7bdat"))
demo_18 = demo.death_data[demo.death_data$age >= 18,]

drug_data_18 = merge(drug_data, demo_18[,c("id", "end", "death")], by = "id")
colnames(drug_data_18)[5] <- "f.up_end"

drug_data_18 = drug_data_18[order(drug_data_18$id,drug_data_18$start.date),]

write.table(drug_data_18, "path/drug_data_C08CA01_shift30_single.episode_18plus.txt", row.names = F)

