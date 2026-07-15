# ---
# title: Wrangle ESMs using R package 'hotrstuff'
# output: ESM projection for Chl-a 
# ---

#Workflow adapted from: 
# https://github.com/SnBuenafe/hotrstuff.
# Citation: Buenafe K, Schoeman D, Everett J (2025). hotrstuff: Facilitates the rapid download, wrangling and processing of Earth System Model (ESM) output from the Coupled Model Intercomparison Project (CMIP). R package version 0.0.2, https://github.com/SnBuenafe/hotrstuff.

# #Install package
# devtools::install_github("SnBuenafe/hotrstuff")

#List of ESMs used in generating ensembles
## IPSL-CM6A-LR
## CMCC-ESM2
## CNRM-ESM2-1
## GFDL-ESM4
## UKESM1-0-LL

# load packages
  library(hotrstuff)
  library(tidyverse)

# define variable
  var <- "chlos"
  base_dir <- "/home/bcalp/UQ"

# to download ESM outputs using the following functions (alternatively, manually download from ESGF)
  #prior step: save wget files in the path directory of 'indir' attribute of 'htr_download_ESM'
  #check wget list
    wget_list <-   list.files(file.path(base_dir, "data", "wget"))  

  #download
    htr_download_ESM(
      hpc = NA,
      indir = file.path(base_dir, "data", "wget"),
      outdir = file.path(base_dir, "data", "raw", var))

  #detect experiment
   #identify ssp scenario
   ensemble_list <- list.files(file.path(base_dir, "data", "raw", var))  
   
   # A. function to detect scenario and assign inputs for later functions
   detect_esm <- function(ensemble_list){
     exp_list <- data.frame(id = numeric(length = length(ensemble_list)),
                            source_id = character(length = length(ensemble_list)),
                            experiment = character(length = length(ensemble_list)))
     
     for(i in 1:length(ensemble_list)){
       ssp_scenario <- str_extract(basename(ensemble_list[i]), "ssp\\d{3}|historical")
       exp_list$id[i] <- i
       exp_list$source_id <- basename(ensemble_list[i])
       exp_list$experiment[i] <- ssp_scenario
     }
     
     print(exp_list %>% group_by(experiment) %>% summarise(count = n()))
     
     exp_lvl_count <- nlevels(factor(exp_list$experiment))
     
     if(exp_lvl_count == 1){
       scenario <- levels(factor(exp_list$experiment))
       print(paste0("Single experiment detected: ", scenario))
       
       assign("ssp_scenario", scenario, envir = .GlobalEnv)
       if(scenario == "historical"){
         year_start = "1850"
         year_end = "2014"
         assign("year_start", year_start, envir = .GlobalEnv)
         assign("year_end", year_end, envir = .GlobalEnv)
       }else if(scenario %in% c("ssp126","ssp245","ssp370","ssp585")){
         year_start = "2015"
         year_end = "2100"
         assign("year_start", year_start, envir = .GlobalEnv)
         assign("year_end", year_end, envir = .GlobalEnv)
       }else{
         print("Experiment cannot be determined")
       }
     }else{
       print("Multiple experiments detected")
     }
   }
   detect_esm(ensemble_list)
   
   # B. type in values:
   # year_start <- "1850"
   # year_end <- "2014"
   # scenario <- "historical"

# to merge files
  htr_merge_files(
    indir = file.path(base_dir, "data", "raw", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "merged", var), # output directory
    year_start = year_start, # earliest year across all the scenarios considered (e.g., historical, ssp126, ssp245, ssp585)
    year_end = year_end # latest year across all the scenarios considered
  )

# to adjust and reframe time periods to 2015-2100 for ScenarioMIPs and 1850-2014 for historical data of ESMs
  htr_slice_period(
    indir = file.path(base_dir, "data", "proc", "merged", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "sliced", var), # output directory
    freq = "Omon", # ocean, daily
    scenario = ssp_scenario,
    year_start = year_start,
    year_end = year_end,
    overwrite = FALSE
  )

# to fix calendar periods (if needed in case the covered time period has a leap year)
  htr_fix_calendar(indir = file.path(base_dir, "data", "proc", "sliced", var)) # will be rewritten

# to change frequency of climate data to yearly
  htr_change_freq(
    freq = "yearly",
    indir = file.path(base_dir, "data", "proc", "sliced", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "yearly", var)
  )

# to regrid into 1.0 x 1.0 cell size (make it consistent for all projections)
  htr_regrid_esm(
    indir = file.path(base_dir, "data", "proc", "yearly", var),
    outdir = file.path(base_dir, "data", "proc", "regridded", "yearly", var),
    cell_res = 1.0,
    layer = "annual"
  )
  
#to generate ensemble
#ensemble by mean
  htr_create_ensemble(
    indir = file.path(base_dir, "data", "proc", "regridded", "yearly", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "ensemble", "mean", var), # output directory
    model_list = c("ACCESS-ESM1-5","CanESM5-1","CESM2-WACCM", "CMCC-ESM2","IPSL-CM6A-LR","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM"), # list of models for ensemble
    variable = var, # variable name
    freq = "Omon", # original frequency of data
    scenario = ssp_scenario, # scenario
    mean = TRUE # if false, takes the median
  )
  # #to plot
  # ensemble_model <- list.files(file.path(base_dir, "data", "proc", "ensemble", "mean", var), full.names = TRUE)
  # ensemble <- rast(ensemble_model)
  # plot(ensemble$chlos_86)

#ensemble by median
  htr_create_ensemble(
    indir = file.path(base_dir, "data", "proc", "regridded", "yearly", var), # input directory
    outdir = file.path(base_dir, "data", "proc", "ensemble", "median", var), # output directory
    model_list = c("ACCESS-ESM1-5","CanESM5-1","CESM2-WACCM", "CMCC-ESM2","IPSL-CM6A-LR","MPI-ESM1-2-HR","MPI-ESM1-2-LR","NorESM2-LM","NorESM2-MM"), # list of models for ensemble
    variable = var, # variable name
    freq = "Omon", # original frequency of data
    scenario = ssp_scenario, # scenario
    mean = FALSE # if false, takes the median
  )
  # #to plot
  # ensemble_model <- list.files(file.path(base_dir, "data", "proc", "ensemble", "median", var), full.names = TRUE)
  # ensemble <- rast(ensemble_model)
  # plot(ensemble$chlos_86)

##End
#05062026 - added historical scenario
#09062026 - removed CanESM5; kept CanESM5-1
#12061015 - updated ensemble composition following Petrik et al (2022)  