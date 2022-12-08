### This code extracts actual event times in the CNODES data, to be matched with the true drug exposure using the permutational algorithm ###

rm(list = ls())
gc()

drug_data = read.table("path/drug_data_C08CA01_shift30_single.episode.txt", header = T)

library(haven)
demo.death_data = as.data.frame(read_sas("path/sd_table_death.sas7bdat"))
demo_18 = demo.death_data[demo.death_data$age >= 18,]

drug_data_18 = merge(drug_data, demo_18[,c("id", "end", "death")], by = "id")
colnames(drug_data_18)[5] <- "f.up_end"

drug_data_18 = drug_data_18[order(drug_data_18$id,drug_data_18$start.date),]

## IDENTIFY DEATH DISTRIBUTION ##

drug_data_18_unique = drug_data_18[!duplicated(drug_data_18$id),]

drug_data_18_unique$death.time = as.numeric(as.Date(drug_data_18_unique$f.up_end) - as.Date(drug_data_18_unique$start.date))

death.times = round(drug_data_18_unique$death.time[drug_data_18_unique$death == 1],0)

write.table(death.times, "path/death_times.txt", row.names = F, col.names = F)

