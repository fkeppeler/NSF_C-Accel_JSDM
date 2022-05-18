#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(dplyr, Hmsc, tidyverse,maps,sf,rgeos)

#Load data
load ('data/processed/db_for_HMSC.RData')

col_to_retain <- c("haul_id", "lon", "lat", "year", "latgrid","longrid",
                   "complex_cl", "seabed_forms", "GRAINSIZE",
                   "monthly_mean_Po2_bott",
                   "monthly_mean_salt_bott", "monthly_mean_salt_surf",
                   "monthly_mean_temp_surf","monthly_mean_temp_bott", 
                   "annual_max_temp_surf", "annual_min_temp_surf",
                   "annual_max_temp_bott","annual_min_temp_bott")

hauls_filtered  <- transform(haul_env,
                             latgrid = lat  %>% round (digits=1),
                             longrid = lon %>% round (digits=1)) %>%
  dplyr::select (one_of(col_to_retain)) %>%
  filter (seabed_forms < 50) %>% #Remove points that felt out of the landmask
  mutate (complex_cl = as.factor (complex_cl),
          seabed_forms = as.factor (seabed_forms))
  
####################

Lat_long <- hauls_filtered[,c("lat","lon")]

Lat_long_sf = st_as_sf(Lat_long, coords = c("lon", "lat"))

# Create grid
AtlGrid <- st_make_grid(x = Lat_long_sf, what = "polygons", cellsize = 1)
AtlGrid  <- st_sf(idcell = 1:length(AtlGrid), geom = AtlGrid ) %>% 
  st_cast("POLYGON")

intersection <- st_intersects(AtlGrid, Lat_long_sf)

group_sites <- data.frame ()
for (i in 1:length(intersection)){
  temp_df <- data.frame (su_id = intersection [[i]], grid = rep (i,intersection [[i]] %>% length()))
  group_sites <- rbind (group_sites, temp_df)
  
}

group_sites <- group_sites[!duplicated (group_sites$su_id), ]

hauls_filtered$su_id <- 1:nrow(hauls_filtered)
hauls_filtered <- base::merge (hauls_filtered, group_sites,
                               by = "su_id", )


ctrs <- lapply( unique( hauls_filtered$grid) , function(x) gCentroid( SpatialPoints( hauls_filtered[ hauls_filtered$grid == x , 
                                                                                                     c('lon','lat') ] ) ) )
ctrs <- setNames( ctrs , unique(hauls_filtered$grid))

ctrs_df <- data.frame()
for (i in 1:length(ctrs)){
  ctrs_df[i,"Long_LS"] <- ctrs[[i]]$x
  ctrs_df[i,"Lat_LS"] <- ctrs[[i]]$y
  ctrs_df[i,"grid"] <- names(ctrs)[i] %>% as.numeric()
}

ggplot() +
  geom_point(data=hauls_filtered, aes (x=lon,y=lat,col=hauls_filtered$grid %>% factor()))+
  geom_point(data=ctrs_df, mapping= aes (x=Long_LS,y=Lat_LS))+
  geom_sf(data=AtlGrid,aes(),fill=NA)+
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

hauls_filtered <- base::merge (hauls_filtered, ctrs_df, by="grid")

rownames(hauls_filtered) <- hauls_filtered$haul_id

#####################
#Order matrices
hauls_filtered <- hauls_filtered [order(hauls_filtered %>% rownames()),]
wide_dat <- wide_dat [rownames (wide_dat) %in% hauls_filtered$haul_id,]
wide_dat<- wide_dat [order(wide_dat %>% rownames()),]

#Threshold to be considered a rare species
threshold_common_sp <- (nrow (wide_dat) * 0.05) %>% floor()

#Create new variable associated with vessel type
hauls_filtered <- hauls_filtered %>%
  mutate (Vessel = ifelse (year < 2009, "Albatross", "Bigelow") %>% factor())

#Define train datasets
dt <- list ()
dt_names <- c('extrap','interp')

dt [[1]] <- hauls_filtered %>% subset (year <= 2011) %>% rownames() #extrapolation

set.seed(666)
dt [[2]] <- sample(rownames(wide_dat), dt [[1]] %>% length()) #interpolation

#Define response matrices 
Y <- as.matrix(wide_dat)
Ypa <- 1*(Y>0) #presence-absence
Ypa_cs <- Ypa[,colSums(Ypa) >= threshold_common_sp]

#Define test datasets
extrap_datasets <- list ()
extrap_datasets[[1]] <- Ypa[!(row.names(Ypa) %in% dt[[1]]),]
extrap_datasets[[2]] <- Ypa_cs[!(row.names(Ypa_cs) %in% dt[[1]]),]
extrap_datasets[[3]] <- hauls_filtered[!(row.names(hauls_filtered) %in% dt[[1]]),]
names(extrap_datasets) <- c('Ypa_test','Ypa_cs_test','hauls_test')

interp_datasets <- list ()
interp_datasets[[1]] <- Ypa[!(row.names(Ypa) %in% dt[[2]]),]
interp_datasets[[2]] <- Ypa_cs[!(row.names(Ypa_cs) %in% dt[[2]]),]
interp_datasets[[3]] <- hauls_filtered[!(row.names(hauls_filtered) %in% dt[[2]]),]
names(interp_datasets) <- c('Ypa_test','Ypa_cs_test','hauls_test')

save (extrap_datasets, interp_datasets,
      file='data/test/test_data_LS.RData')

#For loop to create models 
m_list <- list ()
k <- 1
for (i in 1:length(dt)){

  #Defining training datasets 
  Ypa_train <-Ypa[dt[[i]],]
  Ypa_cs_train <-Ypa_cs[dt[[i]],]
  hauls_train <- hauls_filtered [dt[[i]],]
  
  #Env Data 
  XData <- as.data.frame(hauls_train)
  XData$haul_id <- as.factor (XData$haul_id)
  

  XFormula <- ~  poly(monthly_mean_Po2_bott, degree = 2, raw = TRUE)+ 
    poly(monthly_mean_salt_bott, degree = 2, raw = TRUE)+
    poly(monthly_mean_salt_surf, degree = 2, raw = TRUE)+
    poly(monthly_mean_temp_surf, degree = 2, raw = TRUE)+
    poly(monthly_mean_temp_bott, degree = 2, raw = TRUE)+
    poly(annual_max_temp_surf, degree = 2, raw = TRUE)+
    poly(annual_min_temp_surf, degree = 2, raw = TRUE)+
    poly(annual_max_temp_bott, degree = 2, raw = TRUE)+
    poly(annual_min_temp_bott, degree = 2, raw = TRUE)+
    complex_cl +
    seabed_forms +
    poly(GRAINSIZE, degree = 2, raw = TRUE) + 
    Vessel
  
  
  #Study design
  hauls_train  <- transform(hauls_train,
                            LatLong_id = paste0 (latgrid,longrid) %>% factor() %>% as.numeric)
  
  studyDesign <- data.frame (Year = factor(hauls_train$year),
                             LatLong_LS = factor(hauls_train$grid))
  
  #Random component
  Coord_coarse <- cbind.data.frame (studyDesign$LatLong_LS,hauls_train[,c("Long_LS","Lat_LS")]) %>%
    unique ()
  row.names(Coord_coarse) <- Coord_coarse [,1]
  Coord_coarse [1] <- NULL
  colnames(Coord_coarse) <- c("Long_LS","Lat_LS")
  rL.spatial_coarse <- HmscRandomLevel(sData = Coord_coarse, longlat = FALSE)
  
  Time <- cbind.data.frame (studyDesign$Year,hauls_train[,c("year")]) %>% 
    unique ()
  row.names(Time) <- Time [,1]
  Time [1] <- NULL
  colnames(Time) <- c("Year") 
  rL.temporal <- HmscRandomLevel(sData = Time)
  
  
  #Models
  m_list [[k]] <- Hmsc(Y = Ypa_train,
                       XData = XData,
                       XFormula = XFormula,
                       distr="probit",
                       studyDesign=studyDesign,
                       ranLevels=list(Year = rL.temporal,
                                      LatLong_LS = rL.spatial_coarse
                       ))
  
  names (m_list)[k] <- paste0 ('m.Ypa_',dt_names[i])
  
  
  m_list [[k+1]] <- Hmsc(Y = Ypa_cs_train,
                         XData = XData,
                         XFormula = XFormula,
                         distr="probit",
                         studyDesign=studyDesign,
                         ranLevels=list(Year = rL.temporal,
                                        LatLong_LS = rL.spatial_coarse
                         ))
  names (m_list)[k+1] <- paste0 ('m.Ypa.cs_',dt_names[i])
  
  k <- k + 2
  
  for (j in colnames(Ypa_train)){
    m_list [[k]] <- Hmsc(Y = Ypa_train[,j,drop = FALSE],
                         XData = XData,
                         XFormula = XFormula,
                         distr="probit",
                         studyDesign=studyDesign,
                         ranLevels=list(Year = rL.temporal,
                                        LatLong_LS = rL.spatial_coarse
                                        ))
    names (m_list)[k] <- paste (j,dt_names[i],sep='_')
    k <- k + 1
    
  }
}


saveRDS (m_list,file="unfitted_models/m_list_LS.rds")




