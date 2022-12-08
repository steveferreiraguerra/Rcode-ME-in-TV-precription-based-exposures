### This code creates the observed prescription data from the previously generated amlodipine dataset ###
###              The resulting data is used in scenarios AAAA, ABAA, ACAA, AABA, and BAAA             ###
###     For scenarios AAAB and AAAC, grace periods of 30 and 60 days are introduced, respectively     ###

rm(list = ls())
gc()

# LOAD THE DATASET ALREADY CONTAINING THE INCLUSION CRITERIA AND THE CORRECTED START DATES

drug_data = read.table("path/drug_data_C08CA01_shift30_single.episode_18plus.txt", header = T)
drug_data$start.date = as.Date(drug_data$start.date)
drug_data$f.up_end = as.Date(drug_data$f.up_end)

# SUM DURATION OF DUPLICATED START DATES

drug_data_no.dup = aggregate(dur ~ id + start.date + f.up_end, drug_data, sum)

# CREATE RX STOP DATE

drug_data_no.dup$stop.date = as.Date(drug_data_no.dup$start.date + drug_data_no.dup$dur - 1)

# ORDER AND RENAME ROWS

drug_data_ordered = drug_data_no.dup[order(drug_data_no.dup$id,drug_data_no.dup$start.date),]
row.names(drug_data_ordered) = 1:nrow(drug_data_ordered)

# DETERMINING GAPS AND OVERLAPS IN BETWEEN CONSECUTIVE RXs -- NEGATIVE GAPS INDICATE OVERLAPS

drug_data_ordered$lag_stop.date = c(as.Date(NA), drug_data_ordered$stop.date[1:(length(drug_data_ordered$stop.date)-1)])
drug_data_ordered$lag_stop.date[!duplicated(drug_data_ordered$id)] = as.Date(NA)
drug_data_ordered$gap =  as.numeric(drug_data_ordered$start.date - drug_data_ordered$lag_stop.date) - 1

# REMOVING OVERLAPS -- END OVERLAPPED RX IF HALF THE DURATION OR LESS WAS TAKEN 
# OTHERWISE ADD NUMBER OF OVERLAPPING DAYS TO THE END OF OVERLAPPING RX, OR SUBSEQUENT POSITIVE GAP

drug_data_ordered$new.lag_stop.date = ifelse(drug_data_ordered$gap <= -drug_data_ordered$dur/2,
                                             drug_data_ordered$start.date - 1, drug_data_ordered$lag_stop.date)
drug_data_ordered$new.lag_stop.date = as.Date(drug_data_ordered$new.lag_stop.date, origin = "1970-01-01")

drug_data_ordered$new.stop.date = c(tail(drug_data_ordered$new.lag_stop.date, -1), NA)
drug_data_ordered$new.stop.date[is.na(drug_data_ordered$new.stop.date)] = drug_data_ordered$stop.date[is.na(drug_data_ordered$new.stop.date)]
drug_data_ordered$new.gap =  as.numeric(drug_data_ordered$start.date - drug_data_ordered$new.lag_stop.date) - 1

drug_data_ordered$new.dur = as.numeric(drug_data_ordered$new.stop.date - drug_data_ordered$start.date + 1)

overlap.correction = drug_data_ordered
overlap.correction$new.start.date = overlap.correction$start.date

any.overlap = 1
while (any.overlap > 0) {

  overlap.correction$new.start.date[which(overlap.correction$new.gap < 0)] = overlap.correction$new.lag_stop.date[which(overlap.correction$new.gap < 0)] +  1

  overlap.correction$new.stop.date = as.Date(overlap.correction$new.start.date + overlap.correction$dur - 1)

  overlap.correction$new.lag_stop.date = c(as.Date(NA), head(overlap.correction$new.stop.date,-1))
  overlap.correction$new.lag_stop.date[!duplicated(overlap.correction$id)] = as.Date(NA)

  overlap.correction$new.gap =  as.numeric(overlap.correction$new.start.date - overlap.correction$new.lag_stop.date) - 1

  any.overlap = sum(na.omit(overlap.correction$new.gap) < 0)

}

drug_data_no.overlaps = overlap.correction[,c("id","f.up_end","new.start.date","new.dur","new.stop.date","new.lag_stop.date", "new.gap")]

# REMOVE START DATES THAT ARE AFTER THE FOLLOW-UP END DATE

drug_data_no.overlaps_cens = drug_data_no.overlaps[which(as.Date(drug_data_no.overlaps$new.start.date) < as.Date(drug_data_no.overlaps$f.up_end)),]

# CREATING INFORMATION ON THE DRUG EPISODE

drug_data_no.overlaps_cens$RX.cens = ifelse(drug_data_no.overlaps_cens$new.stop.date > drug_data_no.overlaps_cens$f.up_end, 1, 0 )
drug_data_no.overlaps_cens$Rx.episode = ave(1:nrow(drug_data_no.overlaps_cens), drug_data_no.overlaps_cens$id, FUN = seq_along)
drug_data_no.overlaps_cens$next_RX = c(tail(drug_data_no.overlaps_cens$new.gap, -1), NA)
drug_data_no.overlaps_cens$last.RX = as.numeric(is.na(drug_data_no.overlaps_cens$next_RX))

# CORRECTING STOP DATES AND DURATION DUE TO ADMINISTRATIVE CENSORING

drug_data_no.overlaps_cens$new.stop.date[drug_data_no.overlaps_cens$RX.cens == 1] = drug_data_no.overlaps_cens$f.up_end[drug_data_no.overlaps_cens$RX.cens == 1]
drug_data_no.overlaps_cens$new.dur[drug_data_no.overlaps_cens$RX.cens == 1] = as.numeric(drug_data_no.overlaps_cens$new.stop.date[drug_data_no.overlaps_cens$RX.cens == 1] -
                                                                                           drug_data_no.overlaps_cens$new.start.date[drug_data_no.overlaps_cens$RX.cens == 1] + 1)

# EXTENDING THE LAST FOLLOW-UP TIME UNTIL THE FOLLOW-UP END DATE

drug_data_no.overlaps_cens$next_RX[drug_data_no.overlaps_cens$last.RX == 1] = drug_data_no.overlaps_cens$f.up_end[drug_data_no.overlaps_cens$last.RX == 1] -
  drug_data_no.overlaps_cens$new.stop.date[drug_data_no.overlaps_cens$last.RX == 1]

# CALCULATING A MEASURE OF ADHERENCE IN BETWEEN TWO RXs -- TO BE USED IN ASSIGNING THE TRUE DRUG EXPOSURE

drug_data_no.overlaps_cens$MPR = drug_data_no.overlaps_cens$new.dur/(drug_data_no.overlaps_cens$new.dur + drug_data_no.overlaps_cens$next_RX)

# FOR SCENARIOS AAAB AND AAAC ONLY -- INTRODUCE GRACE PERIODS
# IN THE FINAL DATA CREATION, THESE NEW VARIABLES REPLACE new.dur and next_RX

# drug_data_no.overlaps_cens$new.dur_gap30 = ifelse(drug_data_no.overlaps_cens$next_RX < 30, drug_data_no.overlaps_cens$new.dur + drug_data_no.overlaps_cens$next_RX,
#                                                   drug_data_no.overlaps_cens$new.dur + 30)
# drug_data_no.overlaps_cens$next_RX_gap30 = ifelse(drug_data_no.overlaps_cens$next_RX < 30, 0, drug_data_no.overlaps_cens$next_RX - 30)
# 
# drug_data_no.overlaps_cens$new.dur_gap60 = ifelse(drug_data_no.overlaps_cens$next_RX < 60, drug_data_no.overlaps_cens$new.dur + drug_data_no.overlaps_cens$next_RX,
#                                                   drug_data_no.overlaps_cens$new.dur + 60)
# drug_data_no.overlaps_cens$next_RX_gap60 = ifelse(drug_data_no.overlaps_cens$next_RX < 60, 0, drug_data_no.overlaps_cens$next_RX - 60)

#### CREATE FINAL EXPOSURE DATA ####

exposure_lenght = c(rbind(drug_data_no.overlaps_cens$new.dur, drug_data_no.overlaps_cens$next_RX))
exposure = rep(c(1,0), dim(drug_data_no.overlaps_cens)[1])
exposure_id = rep(as.numeric(row.names(table(drug_data_no.overlaps_cens$id))), table(drug_data_no.overlaps_cens$id)*2)
exposure_MPR = c(rbind(drug_data_no.overlaps_cens$MPR, drug_data_no.overlaps_cens$MPR))
exposure_RX.episode = c(rbind(drug_data_no.overlaps_cens$Rx.episode, drug_data_no.overlaps_cens$Rx.episode))
exposure_RX.cens = c(rbind(drug_data_no.overlaps_cens$RX.cens, drug_data_no.overlaps_cens$RX.cens))
exposure_last.RX = c(rbind(drug_data_no.overlaps_cens$last.RX,drug_data_no.overlaps_cens$last.RX))
exposure_max.pills = c(rbind(drug_data_no.overlaps_cens$new.dur,drug_data_no.overlaps_cens$new.dur))

exposure.dat = as.data.frame(cbind(exposure_id, exposure, exposure_lenght, exposure_MPR, exposure_RX.episode,
                                   exposure_RX.cens, exposure_last.RX, exposure_max.pills )) 

write.table(exposure.dat, "path/exposure_dat_Scenario_AAAA.txt", row.names = F)

