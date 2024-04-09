###
#   Author: 
#   Max Westphal
###
#   Date:
#   2022-08-21
###
#   Project:
#   cases_simstudy
###
#   Task: 
#   cases_simstudy_lfc.R: Comparison of multiple comparison procedures 
#   for diagnostic accuracy studies implemented in cases package
#   (ROC: realistic parameter configurations from generative 'biomarker' model)
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
sim <- "cases_SIM_REV"                            
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
prob[["generate_instance_roc"]] <- 
  CJ(nrep = 100,
     n = c(100, 200, 400, 800, 4000),
     prev = c(0.25),  
     random = FALSE,
     m = c(5), 
     auc = c("seq(0.85, 0.90, length.out = 5)"),   
     rhose = c(0.5), 
     dist = c("normal"),
     e = c(100), 
     k = c(100),
     delta = 0) %>%
  mutate(rhosp = rhose)  

prob[["generate_instance_lfc"]] <- 
  CJ(nrep = 100,
     n = c(100, 200, 400, 800, 4000), 
     prev = c(0.25), 
     random = FALSE,
     m = c(5), 
     se = c(0.8),
     sp = 0.8,
     L = 1,
     rhose = c(0.5),  
     cortype = c("equi")) %>% 
  mutate(rhosp = rhose)


## setup algo parameters: 
algo <- list()
algo[["process_instance"]] <- 
  CJ(contrast = "define_contrast()",
     benchmark = 0.5,
     alpha = 0.025,
     alternative = "greater",
     adjustment = c("bonferroni", "mbeta"),
     transformation = "none",
     analysis = "co-primary",
     regu = "c(1,1/2,1/4)",
     pars = c("list()",
              "list(lfc_type='full')",
              "list(lfc_type='full', lfc_pr=0.5)",
              "list(lfc_type='full', lfc_pr=0.0625)" # = 0.5^(m-1) for m=5
              )
  ) %>% 
  filter(! (adjustment!="mbeta" & pars != "list()" ))

message("Number of problem parameters: ", nrow(prob[[1]]))
message("Number of problem parameters: ", nrow(prob[[2]]))
message("Number of algorithm arguments: ", nrow(algo[[1]]))


# REGISTRY SETUP ------------------------------------------------------------------------------
## create/load registry:
reg <- makeExperimentRegistry(file.dir = reg.dir, seed = 1)
reg <- loadRegistry(reg.dir, writeable = TRUE)

reg$packages <- req.packages
reg$cluster.functions <- makeClusterFunctionsSocket(round(cr*detectCores()))

## add problems/algorithms:
addProblem(name="generate_instance_roc", data=NULL, 
           fun=cases::generate_instance_roc, reg=reg) 
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
j <- 5
r <- testJob(j, reg=reg)
r
unwrap(getJobPars(j, reg=reg))


# Execute jobs --------------------------------------------------------------------------------
tic()
# jj <- findExperiments(prob.pars = m == 10) %>% ijoin(findNotDone())
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

## split up results:
ids1 <- findExperiments(prob.name = "generate_instance_roc")
JP1 <- getJobPars(ids1) %>% unwrap()

tic()
R1 <- reduceResultsList(ids1, reg=reg, fun = as.data.table) %>% rbindlist(fill=TRUE)
toc()

D1 <- merge(JP1, R1)
head(D1)
dim(D1)

readr::write_csv2(D1, file.path(out.dir, paste0(sim, "_roc", ".csv")))


ids2 <- findExperiments(prob.name = "generate_instance_lfc")
JP2 <- getJobPars(ids2) %>% unwrap()

tic()
R2 <- reduceResultsList(ids2, reg=reg, fun = as.data.table) %>% rbindlist(fill=TRUE)
toc()

D2 <- merge(JP2, R2)
head(D2)
dim(D2)

readr::write_csv2(D2, file.path(out.dir, paste0(sim, "_lfc", ".csv")))  
