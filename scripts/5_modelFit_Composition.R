#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, dplyr, ggplot2, 
               ggExtra, cowplot, ggExtra, vegan, 
               sdmCom, tidyr, cowplot)
# sdmCom is available for download at https://github.com/aminorberg/SDM-comparison

#Functions 
#This R function takes a list object filled with dataframes (same number of named columns) and converts them into an array with 
#the same number of rows, the same number of columns and length(list) number of sheets.
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}

#Load test data
load ("data/test/test_data_LS_2012-4.RData")

#Type of models
pred_type <- c("Extrapolation", "Interpolation")
model_type <- c('JSDM_All','JSDM_Abun','SDM_All','SDM_Abun')

#Define name of abundant species
load ("eval_models/Fit_eval_model_2_thin500_LS.RData")
abun_sp_names_extra <- colnames(predY_real[[1]]) 

load ("eval_models/Fit_eval_model_262_thin500_LS.RData")
abun_sp_names_int <- colnames(predY_real[[1]]) 

# Calculate performance metrics for each species in each model

n_spp <- 258
n_models <- 520

So_Comm <- data.frame (model_type=NA, pred_type=NA,
                       accuracy_Rich=NA, discrimination_Rich=NA, calibration_Rich=NA, precision_Rich=NA,
                       accuracy_Comm_sim=NA, discrimination_Comm_sim=NA, calibration_Comm_sim=NA, precision_Comm_sim=NA,
                       accuracy_Comm_sne=NA, discrimination_Comm_sne=NA, calibration_Comm_sne=NA, precision_Comm_sne=NA,
                       accuracy_Comm_sor=NA, discrimination_Comm_sor=NA, calibration_Comm_sor=NA, precision_Comm_sor=NA)

ct_final <- 1
for (i in 1:length (pred_type)){
  for (k in 1:length (model_type)){
    
    
    #Define test data & load model predictions
    
    if (pred_type [i] == "Extrapolation" & model_type [k] == "JSDM_All"){
      Y_test <- extrap_datasets [[1]]
      
      load ("eval_models/Fit_eval_model_1_thin500_LS.RData")
      
      EpredY_samples <- sapply (predY_real, FUN=rowSums)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Extrapolation" & model_type [k] == "JSDM_Abun"){
      Y_test <- extrap_datasets [[2]]
      load ("eval_models/Fit_eval_model_2_thin500_LS.RData")
      
      EpredY_samples <- sapply (predY_real, FUN=rowSums)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Extrapolation" & model_type [k] == "SDM_All"){
      Y_test <- extrap_datasets [[1]]
      
      Group_Y <- list ()
      Group_Y2 <- list ()
      ct <- 1
      
      for (j in 3:260){
        file_name <- paste0 ("eval_models/Fit_eval_model_",j,"_thin500_LS.RData")
        load (file_name)
        
        #Richness
        n <- length(predY_prob[[1]])
        DF <- structure(predY_real, row.names = c(NA, -n), class = "data.frame")
        
        features <- c(sprintf("Post_%02d", seq(1,dim(DF)[2])))
        features2 <- c(sprintf("Site_%02d", seq(1,dim(DF)[1])))
        colnames (DF) <- features 
        rownames (DF) <- features2 
        
        M <- as.matrix(DF)
        
        Group_Y [[ct]] <- M
        
        #Composition
        if (length(Group_Y2) ==0){
          Group_Y2 <- predY_real
        } else {
          Group_Y2 <- Map(cbind, Group_Y2, predY_real)
        }
        
        ct <- ct + 1
      }
      
      EpredY_samples <- Reduce("+", Group_Y)
      EpredY <- rowMeans(EpredY_samples)
      
      
    } else if (pred_type [i] == "Extrapolation" & model_type [k] == "SDM_Abun"){
      Y_test <- extrap_datasets [[2]]
      
      Group_Y <- list ()
      Group_Y2 <- list ()
      
      ct <- 1
      
      for (j in 3:260){
        file_name <- paste0 ("eval_models/Fit_eval_model_",j,"_thin500_LS.RData")
        load (file_name)
        
        if (! colnames(predY_real[[1]]) %in% abun_sp_names_extra){
          next
        }
        #Richness
        n <- length(predY_prob[[1]])
        DF <- structure(predY_real, row.names = c(NA, -n), class = "data.frame")
        
        features <- c(sprintf("Post_%02d", seq(1,dim(DF)[2])))
        features2 <- c(sprintf("Site_%02d", seq(1,dim(DF)[1])))
        colnames (DF) <- features 
        rownames (DF) <- features2 
        
        M <- as.matrix(DF)
        
        Group_Y [[ct]] <- M
        
        #Composition
        if (length(Group_Y2) ==0){
          Group_Y2 <- predY_real
        } else {
          Group_Y2 <- Map(cbind, Group_Y2, predY_real)
        }
        
        ct <- ct + 1
      }
      
      EpredY_samples <- Reduce("+", Group_Y)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Interpolation" & model_type [k] == "JSDM_All"){
      Y_test <- interp_datasets [[1]]
      load ("eval_models/Fit_eval_model_261_thin500_LS.RData")
      
      EpredY_samples <- sapply (predY_real, FUN=rowSums)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Interpolation" & model_type [k] == "JSDM_Abun"){
      Y_test <- interp_datasets [[2]]
      load ("eval_models/Fit_eval_model_262_thin500_LS.RData")
      
      EpredY_samples <- sapply (predY_real, FUN=rowSums)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Interpolation" & model_type [k] == "SDM_All"){
      Y_test <- interp_datasets [[1]]
      
      Group_Y <- list ()
      Group_Y2 <- list ()
      ct <- 1
      
      for (j in 263:520){
        file_name <- paste0 ("eval_models/Fit_eval_model_",j,"_thin500_LS.RData")
        load (file_name)
        
        #Richness
        n <- length(predY_prob[[1]])
        DF <- structure(predY_real, row.names = c(NA, -n), class = "data.frame")
        
        features <- c(sprintf("Post_%02d", seq(1,dim(DF)[2])))
        features2 <- c(sprintf("Site_%02d", seq(1,dim(DF)[1])))
        colnames (DF) <- features 
        rownames (DF) <- features2 
        
        M <- as.matrix(DF)
        
        Group_Y [[ct]] <- M
        
        #Composition
        if (length(Group_Y2) ==0){
          Group_Y2 <- predY_real
        } else {
          Group_Y2 <- Map(cbind, Group_Y2, predY_real)
        }
        
        ct <- ct + 1
      }
      
      EpredY_samples <- Reduce("+", Group_Y)
      EpredY <- rowMeans(EpredY_samples)
      
    } else if (pred_type [i] == "Interpolation" & model_type [k] == "SDM_Abun"){
      Y_test <- interp_datasets [[2]]
      
      Group_Y <- list ()
      Group_Y2 <- list ()
      ct <- 1
      
      for (j in 263:520){
        file_name <- paste0 ("eval_models/Fit_eval_model_",j,"_thin500_LS.RData")
        load (file_name)
        
        if (! colnames(predY_real[[1]]) %in% abun_sp_names_int){
          next
        }
        #Richness
        n <- length(predY_prob[[1]])
        DF <- structure(predY_real, row.names = c(NA, -n), class = "data.frame")
        
        features <- c(sprintf("Post_%02d", seq(1,dim(DF)[2])))
        features2 <- c(sprintf("Site_%02d", seq(1,dim(DF)[1])))
        colnames (DF) <- features 
        rownames (DF) <- features2 
        
        M <- as.matrix(DF)
        
        Group_Y [[ct]] <- M
        
        #Composition
        if (length(Group_Y2) ==0){
          Group_Y2 <- predY_real
        } else {
          Group_Y2 <- Map(cbind, Group_Y2, predY_real)
        }
        
        ct <- ct + 1
      }
      
      EpredY_samples <- Reduce("+", Group_Y)
      EpredY <- rowMeans(EpredY_samples)
    }
    
    
    ### Summarize results
    
    #Create community object
    
    if (k ==3 | k ==4){
      pred_ary <- list2ary(Group_Y2)
    } else {
      pred_ary <- list2ary(predY_real)
    }
    
    predictions <- list()
    predictions$predictions[[1]][[1]]$predictions <- pred_ary
    Y_test_composition <- Y_test
    
    community_composition <- sdmCom:::calc_betapart(predictions=predictions,
                                                    yvalid=Y_test_composition, 
                                                    n_of_site_pairs = 300, seed = 666)
    
    #Create richness object
    Y_test_richness <- rowSums (Y_test)
    EpredY <- EpredY %>% as.matrix()

    species_richness = list ()
    species_richness$predicted$DT1$temp_data <- EpredY
    species_richness$validation <- Y_test_richness %>% as.matrix()
    
    #Accuracy
    accuracy <- sdmCom:::pm_accuracy_rmse (community_composition = community_composition, 
                                           species_richness = species_richness,
                                           as_array=FALSE)
    So_Comm [ct_final,"accuracy_Rich"] <- accuracy$sp_richness_rmse [[1]][[1]]
    
    So_Comm [ct_final,"accuracy_Comm_sim"] <- accuracy$community_composition_rmse$beta_sim_rmse [[1]][[1]]
    So_Comm [ct_final,"accuracy_Comm_sne"] <- accuracy$community_composition_rmse$beta_sne_rmse [[1]][[1]]
    So_Comm [ct_final,"accuracy_Comm_sor"] <- accuracy$community_composition_rmse$beta_sor_rmse [[1]][[1]]
    
    ##Discrimination
    discrimination <- sdmCom:::pm_discrimination_spearm (community_composition = community_composition, 
                                                         species_richness = species_richness,
                                                         as_array=FALSE)
    So_Comm [ct_final,"discrimination_Rich"] <- discrimination$sp_richness_spearm [[1]][[1]]
    
    So_Comm [ct_final,"discrimination_Comm_sim"] <- discrimination$community_composition_spearm$beta_sim_spear [[1]][[1]]
    So_Comm [ct_final,"discrimination_Comm_sne"] <- discrimination$community_composition_spearm$beta_sne_spear [[1]][[1]]
    So_Comm [ct_final,"discrimination_Comm_sor"] <- discrimination$community_composition_spearm$beta_sor_spear [[1]][[1]]
    
    
    #Re-Create richness object
    species_richness = list ()
    species_richness$predicted$DT1$temp_data <- EpredY_samples
    species_richness$validation <- Y_test_richness %>% as.matrix()
    
    
    ##Calibration
    
    calibration <- sdmCom:::pm_calibration_pred_interval (community_composition = community_composition,
                                                          species_richness = species_richness,
                                                          as_array=FALSE)
    So_Comm [ct_final,"calibration_Rich"] <- calibration$sp_richness_predint [[1]][[1]]
    
    So_Comm [ct_final,"calibration_Comm_sim"] <- calibration$community_composition_predint$beta_sim_predint [[1]][[1]]
    So_Comm [ct_final,"calibration_Comm_sne"] <- calibration$community_composition_predint$beta_sne_predint [[1]][[1]]
    So_Comm [ct_final,"calibration_Comm_sor"] <- calibration$community_composition_predint$beta_sor_predint [[1]][[1]]
    
    #Precision
    precision <- sdmCom:::pm_precision_sd (community_composition = community_composition, 
                                           species_richness = species_richness,
                                           as_array=FALSE)
    So_Comm [ct_final,"precision_Rich"] <- precision$sp_richness_sd [[1]][[1]]
    
    So_Comm [ct_final,"precision_Comm_sim"] <- precision$community_composition_sd$beta_sim_sd [[1]][[1]]
    So_Comm [ct_final,"precision_Comm_sne"] <- precision$community_composition_sd$beta_sne_sd [[1]][[1]]
    So_Comm [ct_final,"precision_Comm_sor"] <- precision$community_composition_sd$beta_sor_sd [[1]][[1]]
    
    
    So_Comm [ct_final,"model_type"] <- model_type [k]
    So_Comm [ct_final,"pred_type"] <- pred_type [i]
    
    print (ct_final)
    ct_final <- ct_final + 1
  }
}

#Save results
write.table (So_Comm,file="results/perf_metrics_com.csv",sep=",",row.names=FALSE)

#Split "model type" col into two new cols
So_Comm <- So_Comm %>%
  separate (model_type, c("Model","Data"), sep = "_")

So_Comm$Data <- rep (c('All species','Only common species'),4)
So_Comm$Model <- factor (So_Comm$Model,levels = c("SDM","JSDM"))

#Richness
P1 <- ggplot (So_Comm,aes (x=Model,y=accuracy_Rich,fill=pred_type))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~Data)+theme_minimal()+
  ylab ('Accuracy')+
  theme (legend.position="none")

P2 <- ggplot (So_Comm,aes (x=Model,y=discrimination_Rich,fill=pred_type))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~Data)+theme_minimal()+
  ylab ('Discrimination')+
  theme (legend.position="none")

P3 <- ggplot (So_Comm,aes (x=Model,y=calibration_Rich,fill=pred_type))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~Data)+theme_minimal()+
  ylab ('Calibration')+
  theme (legend.position="none")

P4 <- ggplot (So_Comm,aes (x=Model,y=precision_Rich,fill=pred_type))+
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~Data)+theme_minimal()+
  ylab ('Precision')+
  theme (legend.position="none")

tiff(file="Results/Richness_all.tiff",width=2000,height=2000,units = "px",
     compression="lzw",res=300)
plot_grid (P1,P2,P3,P4)
dev.off()

ys <- c('Rich','Comm_sim','Comm_sne','Comm_sor')

for (y in ys){
  #Composition sim
  P1 <- ggplot (So_Comm,aes_string (x="Model",y=paste0("accuracy_",y),fill="pred_type"))+
    geom_bar(stat="identity", position=position_dodge())+
    facet_wrap(~Data)+theme_minimal()+
    ylab ('Accuracy')+
    theme (legend.position="none")
  
  P2 <- ggplot (So_Comm,aes_string (x="Model",y=paste0("discrimination_",y),fill="pred_type"))+
    geom_bar(stat="identity", position=position_dodge())+
    facet_wrap(~Data)+theme_minimal()+
    ylab ('Discrimination')+
    theme (legend.position="none")
  
  P3 <- ggplot (So_Comm,aes_string (x="Model",y=paste0("calibration_",y),fill="pred_type"))+
    geom_bar(stat="identity", position=position_dodge())+
    facet_wrap(~Data)+theme_minimal()+
    ylab ('Calibration')+
    theme (legend.position="none")
  
  P4 <- ggplot (So_Comm,aes_string (x="Model",y=paste0("precision_",y),fill="pred_type"))+
    geom_bar(stat="identity", position=position_dodge())+
    facet_wrap(~Data)+theme_minimal()+
    ylab ('Precision')+
    theme (legend.position="none")
  
  tiff(file=paste0("Results/",y,".tiff"),width=2000,height=2000,units = "px",
       compression="lzw",res=300)
  plot_grid (P1,P2,P3,P4) %>% print()
  dev.off()
}






