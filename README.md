# Rcode-ME-in-TV-precription-based-exposures
R code to recreate the simulation in the paper "Impact of measurement error in time-varying prescription-based exposures on estimated hazard ratios from the Cox model" (in preparation for submission to Pharmacoepidemiology and Drug Safety). The code is based on the CNODES Simulated Dataset (https://www.cnodes.ca/simulated-data-set/).

The way the code is structured is that one first creates the amlodipine cohort, from which individuals will be resampled in the simulation, using 
"Amlodipine cohort.R".

The information needed to generate the day-by-day observed prescription-based exposure is contained in the file "Observed prescription data.R".

One then generates the day-by-day observed and true exposure using "True intake.R". This file also contains the code to extract results from a single sample (repetition). Each repetition was ran independently on servers provided by Calcul Québec (https://www.calculquebec.ca/en/) and Compute Canada (www.computecanada.ca). Note that this file uses the Permutational Algorithm (Sylvestre MP, Abrahamowicz M. Comparison of algorithms to generate event times conditional on time-dependent covariates. Statistics in Medicine. 2008;27(14):2618-2634. doi:10.1002/sim.3092) to assign an event time (generated in "Event time extraction.R") to each true exposure history. The permutational algorithm is here coded manually in "perm_algo_function.R", given limitations regarding the data type with the current R package.
