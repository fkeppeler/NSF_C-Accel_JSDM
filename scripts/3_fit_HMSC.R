#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, Hmsc, snow)

#Load models
ms <- readRDS ("unfitted_models/m_list_LS.rds")

#Set parameters
samples = 250
thin = 500
nChains = 5
nParallel = 5

# Fit each model individually. 
#CAUTION. 
# Each model takes a long time to run (hours to days).
# For this project, we run models in parallel using computer clusters at 
# Uw Madison and Rutgers

for (i in 1:length (ms)){
  model_i <- ms [[i]]
  
  start <- Sys.time()
  m = sampleMcmc(model_i,  samples = samples, thin=thin,
                 transient = (samples*thin)*0.2,
                 nChains = nChains,nParallel = nParallel,
                 updater=list(GammaEta=FALSE))
  proc_time <- Sys.time() - start
  
  file_name <- paste0 ("model_",i,"_thin",thin,"_LS.RData")
  save(m, proc_time, file="fitted_models/model_1_thin500_LS.RData")
}




