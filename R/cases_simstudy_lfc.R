###
#   Author: 
#   Max Westphal
###
#   Date:
#   2022-08-21
###
#   Project:
#   cases_simstudy: Comparison of multiple  
###
#   Task: 
#   cases_simstudy_lfc.R: Comparison of multiple comparison procedures 
#   for diagnostic accuracy studies implemented in cases package
#   (LFC: least-favorable parameter configurations setting)
###



# PACKAGES ------------------------------------------------------------------------------------
library(cases)
library(batchtools)
library(data.table)
library(parallel)
library(dplyr)
library(readr)
library(tictoc)


# PREPARATION ---------------------------------------------------------------------------------
# rm(list = ls())                           

## simulation name:
sim <- "cases_SIM_LFC"                            
nsim <- 100                         
cr <- 0.8                                 
req.packages <- c("cases") 

## setup folders: 
main.dir <- file.path("E:/cases_SIM")              
reg.dir <- file.path(main.dir, sim)  
out.dir <- file.path(main.dir, "cases_SIM_results")
dir.create(main.dir)
dir.create(out.dir)

## setup problem parameters: 
prob <- list()
prob[["generate_instance_lfc"]] <- 
  CJ(nrep = 100,
     n = c(100, 200, 400, 800, 4000), 
     prev = c(0.1, 0.25), 
     random = FALSE,
     m = c(5, 10), 
     se = c(0.8, 0.9),
     sp = 0.8,
     L = 1,
     rhose = c(0.25, 0.5),  
     cortype = c("equi")) %>% 
  mutate(rhosp = rhose) %>% 
  filter(! (m != 5 & cortype != "equi"))

## setup algo parameters: 
algo <- list()
algo[["process_instance"]] <- 
  CJ(contrast = "cases::define_contrast('raw', NA)",
     benchmark = 0.5,
     alpha = 0.025,
     alternative = "greater",
     adjustment = c("none", "bonferroni", "maxt", "bootstrap", "mbeta"),
     transformation = c("none", "logit"),
     analysis = c("co-primary"),
     regu = c("c(1,1/2,1/4)"), 
     pars = c("list()", "list(type='wild')")) %>% 
  filter(! (adjustment!="bootstrap" & pars != "list()")) %>% 
  filter(! (adjustment!="maxt" & transformation == "logit")) 
  #filter(! ( !(adjustment %in% c("maxt", "mbeta")) & regu != "c(1,1/2,1/4)")) 

message("Number of problem parameters: ", nrow(prob[[1]]))
message("Number of algorithm arguments: ", nrow(algo[[1]]))


# REGISTRY SETUP ------------------------------------------------------------------------------
## create/load registry:
reg <- makeExperimentRegistry(file.dir = reg.dir, seed = 1)
reg <- loadRegistry(reg.dir, writeable = TRUE)

reg$packages <- req.packages
reg$cluster.functions <- makeClusterFunctionsSocket(round(cr*detectCores()))

## add problems/algorithms:
addProblem(name="generate_instance_lfc", data=NULL, 
           fun=cases::generate_instance_lfc, reg=reg) 
addAlgorithm(name="process_instance", 
             fun=cases::process_instance, reg=reg) 

## add experiments:
addExperiments(prob.designs = prob,
               algo.designs = algo,
               repls=nsim,
               reg=reg)

summarizeExperiments(reg=reg)


# Test jobs -----------------------------------------------------------------------------------
# j <- findExperiments(algo.pars = adjustment == "mbeta")
j <- 1
tic()
r <- testJob(j, reg=reg)
toc()
r
unwrap(getJobPars(j, reg=reg))


# Execute jobs --------------------------------------------------------------------------------
tic()
# jj <- findExperiments(prob.pars = m <=10) %>% ijoin(findNotDone())
jj <- findNotDone(reg=reg)
submitJobs(jj, reg=reg)
toc()


# Optional steps ------------------------------------------------------------------------------
getStatus(reg=reg)
# findExpired(reg=reg)
# removeExperiments(..., reg=reg)
# waitForJobs(reg=reg)
# getErrorMessages()
# killJobs(reg=reg)
# clearRegistry(reg=reg)


# Save results --------------------------------------------------------------------------------
# options(batchtools.progress = TRUE)
JP <- getJobPars(findDone(), reg=reg) %>% unwrap()

tic()
R <- reduceResultsList(findDone(), reg=reg, fun = as.data.table) %>% rbindlist(fill=TRUE)
toc()

D <- merge(JP, R)
head(D)
dim(D)

readr::write_csv2(D, file.path(out.dir, paste0(sim, ".csv")))
