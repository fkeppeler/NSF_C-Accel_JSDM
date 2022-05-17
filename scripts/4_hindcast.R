load ("data/test/test_data_LS_2012-4.RData")

#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, abind, dplyr)

#Options for parallel processing
cores= detectCores()-1
cl <- makeCluster(cores) 
registerDoParallel(cl)

#Retrieve name of the model
Sp_Spec_occ <- foreach(
  i= c(1:520),
  .combine=rbind,
  .packages = c("dplyr","Hmsc", "abind") 
) %dopar% {
  
  model_number <- i
  
  file_name <- paste0("fitted_models/model_",model_number,"_thin500_LS.RData")
  
  #Load model
  load (file_name)
  
  #Define test data
  if (model_number == 1){
    Y_test <- extrap_datasets [[1]]
    X_test <- extrap_datasets [[3]]
  } else if (model_number == 2){
    Y_test <- extrap_datasets [[2]]
    X_test <- extrap_datasets [[3]]
  } else if (model_number > 2 & model_number < 261){
    Y_test <- extrap_datasets [[1]][,colnames(m$Y)]
    X_test <- extrap_datasets [[3]]
  } else if (model_number == 261){
    Y_test <- interp_datasets [[1]]
    X_test <- interp_datasets [[3]]
  } else if (model_number == 262){
    Y_test <- interp_datasets [[2]]
    X_test <- interp_datasets [[3]]
  } else {
    Y_test <- interp_datasets [[1]][,colnames(m$Y)]
    X_test <- interp_datasets [[3]]
  }
  
  LatLong <- X_test[,c("lon","lat")] %>% as.matrix ()
  colnames(LatLong) <- c("Long","Lat")
  
  Year <- X_test[,c("year")] %>% as.matrix ()
  colnames(Year) <- "Year"
  
  #Prediction 
  Gradient <- prepareGradient (hM=m, XDataNew=X_test, 
                               sDataNew= list(LatLong_LS=LatLong,
                                              Year=Year))
  set.seed(666)
  predY_prob <- predict (object=m, Gradient= Gradient, predictEtaMean=TRUE,expected=TRUE,
                         post = poolMcmcChains(postList=m$postList,start=226))
  
  set.seed(666)
  predY_real <- predict (object=m, Gradient= Gradient, predictEtaMean=TRUE,expected=FALSE,
                         post = poolMcmcChains(postList=m$postList,start=226))
  
  # Summarize info  
  if (model_number == 1 | model_number == 2 | model_number == 261 | model_number == 262){
    EpredY = Reduce("+", predY_prob)/length(predY_prob)
  } else {
    EpredY = apply(abind(predY_prob,along = 3), c(1,2), mean)
  }
  
  #Save info
  file_name <- paste0("eval_models/Fit_eval_model_",model_number,"_thin500_LS.RData")
  
  save (EpredY, predY_real, predY_prob, file = file_name)
  
}

parallel::stopCluster(cl = cl)
