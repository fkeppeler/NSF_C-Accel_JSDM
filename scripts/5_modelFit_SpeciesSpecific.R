#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(foreach, doParallel, dplyr, ggplot2, 
               ggExtra, cowplot, ggExtra, vegan, sdmCom)
# sdmCom is available for download at https://github.com/aminorberg/SDM-comparison

#Options for parallel processing
cores= detectCores()-1
cl <- makeCluster(cores) 
registerDoParallel(cl)

#Load test data
load ("data/test/test_data_LS_2012-4.RData")

# Calculate performance metrics for each species in each model

n_predictions <- 2 #Extrapolation and Interpolation 
n_spp <- 258
n_models <- 520
n_JSDM <- 2 # One for all species and another with only abundant species


model_type <- rep (
  c( 
    rep('JSDM',2), #Include models with all species and with only abundant species
    rep ('SDM',n_spp)
  ),
  n_predictions)
  
pred_type <- c(
  rep ("Extrapolation", n_spp + n_JSDM),
  rep ("Interpolation", n_spp + n_JSDM)
                  )

Sp_Spec_occ <- foreach(
  i= 1:n_models,
  .combine=rbind,
  .packages = c("dplyr","Hmsc","sdmCom") 
) %dopar%{
  Sp_Spec_occ <- data.frame (model_type=NA, pred_type=NA, type_JSDM = NA,
                             species =NA, Accuracy=NA, Discrimination=NA, Calibration=NA, Precision=NA)
  ct <- 1
  
  file_name <- paste0 ("eval_models/Fit_eval_model_",i,"_thin500_LS.RData")
  load (file_name)
  #Define test data
  if (i == 1){
    Y_test <- extrap_datasets [[1]]
  } else if (i == 2){
    Y_test <- extrap_datasets [[2]]
  } else if (i > 2 & i < 261){
    Y_test <- extrap_datasets [[1]][,colnames(EpredY)]
  } else if (i == 261){
    Y_test <- interp_datasets [[1]]
  } else if (i == 262){
    Y_test <- interp_datasets [[2]]
  } else {
    Y_test <- interp_datasets [[1]][,colnames(EpredY)]
  }
  
  for (k in 1:ncol (EpredY)){
    #Data_required
    occurrence_probs = list ()
    species_name <-colnames (EpredY)[k]
    occurrence_probs$DT1$temp_species <- as.matrix(EpredY[,k])
    
    yvalid <-  as.matrix (Y_test)[,k, drop=FALSE]
    
    #Model type & Species
    Sp_Spec_occ [ct,"model_type"] <- model_type [i]
    Sp_Spec_occ [ct , "pred_type"] <- pred_type [i]
    
    if (i == 1 | i == 261){
      Sp_Spec_occ [ct , "type_JSDM"] <- "All"
    }
    
    if (i == 2 | i == 262){
      Sp_Spec_occ [ct , "type_JSDM"] <- "Abun_species"
    }
    
    Sp_Spec_occ [ct , "species"] <- species_name
    
    #Accuracy 
    Sp_Spec_occ [ct, "Accuracy"] <- sdmCom:::pm_accuracy_expected_me(occurrence_probs = occurrence_probs, 
                                                                     yvalid = yvalid, as_array = TRUE) [3]
    
    #Discrimination
    Sp_Spec_occ [ct, "Discrimination"] <- sdmCom:::pm_discrimination_auc(occurrence_probs = occurrence_probs, 
                                                                         yvalid = yvalid, as_array = TRUE) [3]
    
    #Calibration
    Sp_Spec_occ [ct, "Calibration"] <- sdmCom:::pm_calibration_prob_bin_mse(occurrence_probs = occurrence_probs, 
                                                                            yvalid = yvalid, as_array = TRUE) [3]
    
    #Precision
    Sp_Spec_occ [ct, "Precision"] <- sdmCom:::pm_precision_sqrt_probs(occurrence_probs = occurrence_probs, 
                                                                      as_array = TRUE) [3]
    ct <- ct + 1
  }
  Sp_Spec_occ
}

parallel::stopCluster(cl = cl)


########
#Save data
write.table (Sp_Spec_occ, file="results/perf_metrics.csv",row.names=FALSE,sep=',')
Sp_Spec_occ <- read.csv ("results/perf_metrics.csv")

#Explore results
############
### INCLUDING RARE SPECIES
############

Sp_Spec_occ <- Sp_Spec_occ %>%
  mutate (
    Accuracy = as.numeric (Accuracy),
    Discrimination = as.numeric (Discrimination),
    Calibration = as.numeric (Calibration),
    Precision = as.numeric (Precision),
    model_type = as.factor(model_type),
    pred_type = as.factor (pred_type),
    type_JSDM = as.factor (type_JSDM),
    species = as.factor (species)
  )

JSDM_all <- Sp_Spec_occ %>%
  filter (type_JSDM == "All") %>%
  arrange(species,pred_type)
  
SDM <-  Sp_Spec_occ %>%
  filter (is.na (type_JSDM)) %>%
  arrange(species,pred_type)

#Scatter plots
P1 <- ggplot(data = data.frame(y = SDM$Accuracy, 
                               x = JSDM_all$Accuracy,
                               pred_type = JSDM_all$pred_type),
       aes(x = x, y = y, colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Accuracy at species level", x= "JSDM Accuracy at species level")+
  theme_classic()+theme (legend.position="none")
P1 <- ggMarginal (P1,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P2 <- ggplot(data = data.frame(y = SDM$Discrimination, 
                               x = JSDM_all$Discrimination,
                               pred_type = JSDM_all$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Discrimination at species level", x= "JSDM Discrimination at species level")+
  theme_classic()+theme (legend.position="none")
P2 <- ggMarginal (P2,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P3 <- ggplot(data = data.frame(y = SDM$Calibration, 
                               x = JSDM_all$Calibration,
                               pred_type = JSDM_all$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Calibration at species level", x= "JSDM Calibration at species level")+
  theme_classic()+theme (legend.position="none")
P3 <- ggMarginal (P3,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P4 <- ggplot(data = data.frame(y = SDM$Precision, 
                               x = JSDM_all$Precision,
                               pred_type = JSDM_all$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Precision at species level", x= "JSDM Precision at species level")+
  theme_classic()+theme (legend.position="none")
P4 <- ggMarginal (P4,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

tiff(file="results/Species_all.tiff",width=2000,height=2000,units = "px",
     compression="lzw",res=300)
plot_grid(P1,
          P2,
          P3,
          P4,
          ncol=2,nrow=2)
dev.off()


#density plot

ct <- 1

pred_type <- c("Extrapolation","Interpolation")
param <- c("Accuracy", "Discrimination", "Calibration", "Precision")

results <- list ()
for (j in pred_type){
  for (k in param){
    y <- subset (JSDM_all, pred_type==j)[,k]
    x <- subset (SDM, pred_type==j)[,k]
    
    residuals_1to1 <- resid(lm(y-x ~ 0))
    
    results [[ct]] <- data.frame (residuals_1to1=residuals_1to1,
                                     pred_type=j,
                                     param=k)
    ct <- ct + 1
  }

}
df_results <- do.call ("rbind",results)


Upp_Low_Int <- df_results %>%
  group_by(pred_type,param) %>%
  summarise(
    Lower = quantile (residuals_1to1, probs= (0.025)),
    Medium = quantile (residuals_1to1, probs= (0.5)),
    Upper = quantile (residuals_1to1, probs= (0.975)))

ct <- 1
param_graph <- list ()
for (i in param){
  df <- subset (df_results, param==i)
  Upp_Low_temp <- subset (Upp_Low_Int, param==i)
  
  param_graph [[ct]]<- ggplot(data=df, aes(x=residuals_1to1,y= ..scaled..,
                      fill=pred_type)) +
    geom_density(color=NA)+
    facet_wrap(~pred_type,scales = "free")+
    geom_vline(xintercept = 0)+
    geom_vline(data = Upp_Low_temp, aes(xintercept = Lower,color=pred_type),linetype="dashed")+
    geom_vline(data = Upp_Low_temp, aes(xintercept = Upper,color=pred_type),linetype="dashed")+
    theme_classic()+theme (legend.position="none",plot.title = element_text(hjust = 0.5))+
    labs (y='Scaled density',x='Residuals - 1 to 1 line')+
    ggtitle (i)
    
    ct <- ct + 1
}

tiff(file="results/1to1residuals_all.tiff",width=3000,height=2000,units = "px",
     compression="lzw",res=300)
plot_grid (param_graph [[1]],param_graph [[2]],
           param_graph [[3]],param_graph [[4]])
dev.off()

############
### EXCLUDING RARE SPECIES
############

JSDM_abun <- Sp_Spec_occ %>%
  filter (type_JSDM == "Abun_species") %>%
  arrange(species,pred_type)

SDM <-  Sp_Spec_occ %>%
  filter (species %in% (JSDM_abun$species %>% unique ()) &
            is.na (type_JSDM)) %>%
  arrange(species,pred_type)


#Scatter plots
P1 <- ggplot(data = data.frame(y = SDM$Accuracy, 
                               x = JSDM_abun$Accuracy,
                               pred_type = JSDM_abun$pred_type),
             aes(x = x, y = y, colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Accuracy at species level", x= "JSDM Accuracy at species level")+
  theme_classic()+theme (legend.position="none")
P1 <- ggMarginal (P1,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P2 <- ggplot(data = data.frame(y = SDM$Discrimination, 
                               x = JSDM_abun$Discrimination,
                               pred_type = JSDM_abun$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Discrimination at species level", x= "JSDM Discrimination at species level")+
  theme_classic()+theme (legend.position="none")
P2 <- ggMarginal (P2,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P3 <- ggplot(data = data.frame(y = SDM$Calibration, 
                               x = JSDM_abun$Calibration,
                               pred_type = JSDM_abun$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Calibration at species level", x= "JSDM Calibration at species level")+
  theme_classic()+theme (legend.position="none")
P3 <- ggMarginal (P3,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)

P4 <- ggplot(data = data.frame(y = SDM$Precision, 
                               x = JSDM_abun$Precision,
                               pred_type = JSDM_abun$pred_type), 
             aes(x = x, y = y,colour=pred_type)) +
  geom_point ()+
  geom_abline(intercept=0, slope=1)+
  labs (y= "SDM Precision at species level", x= "JSDM Precision at species level")+
  theme_classic()+theme (legend.position="none")

P4 <- ggMarginal (P4,
                  type = 'hist', groupColour = TRUE, groupFill = TRUE,alpha=0.5)


tiff(file="results/Species_abund.tiff",width=2000,height=2000,units = "px",
     compression="lzw",res=300)
plot_grid(P1,
          P2,
          P3,
          P4,
          ncol=2,nrow=2)
dev.off()


#density plot

ct <- 1

pred_type <- c("Extrapolation","Interpolation")
param <- c("Accuracy", "Discrimination", "Calibration", "Precision")

results <- list ()
for (j in pred_type){
  for (k in param){
    x <- subset (JSDM_abun, pred_type==j)[,k]
    y <- subset (SDM, pred_type==j)[,k]
    
    residuals_1to1 <- resid(lm(y-x ~ 0))
    
    results [[ct]] <- data.frame (residuals_1to1=residuals_1to1,
                                  pred_type=j,
                                  param=k)
    ct <- ct + 1
  }
  
}
df_results <- do.call ("rbind",results)


Upp_Low_Int <- df_results %>%
  group_by(pred_type,param) %>%
  summarise(
    Lower = quantile (residuals_1to1, probs= (0.025)),
    Medium = quantile (residuals_1to1, probs= (0.5)),
    Upper = quantile (residuals_1to1, probs= (0.975)))

ct <- 1
param_graph <- list ()
for (i in param){
  df <- subset (df_results, param==i)
  Upp_Low_temp <- subset (Upp_Low_Int, param==i)
  
  param_graph [[ct]]<- ggplot(data=df, aes(x=residuals_1to1,  y=..scaled..,
                                           fill=pred_type)) +
    geom_density(color=NA)+
    facet_wrap(~pred_type,scales = "free")+
    geom_vline(xintercept = 0)+
    geom_vline(data = Upp_Low_temp, aes(xintercept = Lower,color=pred_type),linetype="dashed")+
    geom_vline(data = Upp_Low_temp, aes(xintercept = Upper,color=pred_type),linetype="dashed")+
    theme_classic()+theme (legend.position="none",plot.title = element_text(hjust = 0.5))+
    labs (y='Scaled density',x='Residuals - 1 to 1 line')+
    ggtitle (i)
  
  ct <- ct + 1
}

tiff(file="results/1to1residuals_abun.tiff",width=3000,height=2000,units = "px",
     compression="lzw",res=300)
plot_grid (param_graph [[1]],param_graph [[2]],
           param_graph [[3]],param_graph [[4]])
dev.off()


