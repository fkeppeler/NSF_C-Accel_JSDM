#Packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyr, tidyverse, dplyr, leaflet, lubridate, corrplot, foreach,
               doParallel, cowplot, vegan, PerformanceAnalytics, sf, stars)

# Define number of cores - parallel processing
n.cores <- 8

#Increase memory limit
memory.limit(size = 100000)

# Hauls Data - Morley
load ("data/trawl/hauls_catch_Dec2019.RData")

cols_hauls <- c('haulid','rugosity',
                'GRAINSIZE','SST.seasonal.mean','SBT.seasonal',
                'SBT.min','SBT.max',
                'SST.min','SST.max',
                'depth')

hauls <- hauls %>% 
  select (one_of(cols_hauls))

#Climate data - Samantha
haul_info <- read.csv ("data/trawl/NEFSC_trawl_data_1980_2014.csv",header=T,sep=",")[,c(1,4,5)]

haul_info <- left_join(haul_info, hauls, 
                        by = c("haul_id" = "haulid"))

var_names <- c("temp_surf", "temp_bott", "salt_bott", "salt_surf", "Po2_bott")

drop_columns <- c("X1980","X1981","X1982","X1983","X1984","X1985","X1986","X1987",
                  "X1988","X1989","X1990","X1991", "X1992", "X1993", "X1994", "X1995", 
                  "X1996", "X1997", "X1998", "X1999" , "X2000", "X2001", "X2002", "X2003", 
                  "X2004", "X2005", "X2006", "X2007", "X2008", "X2009", "X2010", "X2011", 
                  "X2012", "X2013","X2014", "X")

my.cluster <- makeCluster(
  n.cores, 
  type = "PSOCK"
)

registerDoParallel(cl = my.cluster)

m_clim_var <- foreach(i=1:length(var_names), 
                      .combine='rbind', .packages = c("dplyr","tidyr")) %dopar% {
                        file_name <- paste0 ("data/climate/NWA_",var_names[i],".csv")
                        df <- read.csv (file_name,header=T,sep=",")
                        df <- cbind (haul_info, df) 
                        df <- df %>% 
                          select(-all_of(drop_columns)) %>% # Remove columns associated with year-level values
                          drop_na () %>% # remove nas
                          pivot_longer(cols =X1980FEB:X2014DEC,
                                       values_to = "Climate_values") %>% # Transform df - long to wide
                          transform(years = substr(name, 2, 5), 
                                    months = substr(name, 6, 8),
                                    Climate_variable = var_names[i]) %>% # Create year and month columns from 'date id' column (name)
                          select(-name)  # Remove 'name' column
  return (df)
}

parallel::stopCluster(cl = my.cluster)

#filtering
m_clim_var2 <- m_clim_var %>%
  filter (years == year) 

m_clim_var3 <- m_clim_var2 %>% 
  group_by(haul_id, month, lon, lat, Climate_variable,
           year,rugosity, GRAINSIZE, SST.seasonal.mean, SBT.seasonal, 
           SBT.min, SBT.max, SST.min, SST.max, depth) %>%
  summarize(annual_mean = mean(Climate_values, na.rm=TRUE),#
            annual_min = min(Climate_values, na.rm=TRUE),#
            annual_max = max(Climate_values, na.rm=TRUE))#

m_clim_var4 <- m_clim_var2 %>%
  filter (month == months) %>%
  select(haul_id, Climate_values, Climate_variable) %>%
  rename (monthly_mean = Climate_values)


haul_env <- left_join(m_clim_var3, m_clim_var4, 
                      by = c("haul_id","Climate_variable"))

haul_env <- haul_env %>%
  pivot_wider(names_from = Climate_variable, 
              values_from = c(annual_mean, annual_min, annual_max, monthly_mean)) %>%
  arrange(desc(haul_id)) %>%
  drop_na ()

#calculating metrics 
OceanModel<- c("monthly_mean_temp_surf","annual_min_temp_surf", "annual_max_temp_surf",
               "monthly_mean_temp_bott","annual_min_temp_bott", "annual_max_temp_bott")
Trawl <- c("SST.seasonal.mean", "SST.min","SST.max",
           "SBT.seasonal", "SBT.min", "SBT.max")

labels_OceanModel <- c("SST.seasonal.mean.OceanModel","SST.annual.min.OceanModel","SST.annual.max.OceanModel",
                       "SBT.seasonal.mean.OceanModel","SBT.annual.min.OceanModel","SBT.annual.max.OceanModel")
labels_Trawl <- c("SST.seasonal.mean.Trawl", "SST.annual.min.Trawl", "SST.annual.max.Trawl",
                  "SBT.seasonal.mean.Trawl", "SBT.annual.min.Trawl", "SBT.annual.max.Trawl")

gg_plots <- list ()

for (i in 1:length (OceanModel)){
  cor_value <- cor (haul_env[,Trawl[i]], haul_env[,OceanModel[i]])[1] %>%
    round (digits=2)
  
  gg_plots [[i]] <- ggplot (haul_env, aes_string(x=Trawl[i],y=OceanModel[i]))+
    geom_point ()+
    labs (y=labels_OceanModel[i],x=labels_Trawl[i])+
    ggtitle(paste0("Correlation = ",cor_value))+
    theme_minimal()
}

##############################
#Filter biological data
drop_cols_dat <- c('presfit', 'logwtcpue', 'Freq', 'region')

wide_dat <- dat %>%
  separate(sppocean, sep = "_", into=c('species','region')) %>%
  select (-one_of(drop_cols_dat)) %>%
  filter (haulid %in% haul_env$haul_id) %>% 
  pivot_wider (names_from = species, values_from =wtcpue,
               values_fill=0) %>%
  replace(is.na(.), 0) %>%
  arrange(desc(haulid)) %>%
  column_to_rownames(var = "haulid") %>%
  rename_with( ~ gsub (" ",".",.x))

#table (rownames (wide_dat) == haul_env$haul_id)


# Benthic data - Marta
files <- c("bathy_100m",
           "seabed_forms",
           "bathy_classes",
           "bpi_broad_classes",
           "bpi_fine_classes",
           "complex_cl",
           "sedim_class")

my.cluster <- makeCluster(
  n.cores, 
  type = "PSOCK"
)
registerDoParallel(cl = my.cluster)

ben_data <- foreach(i=1:length(files), 
                    .combine='c', .packages = c("dplyr","stars","sf")) %dopar% {
                      file_name <- paste0 ("data/benthic/",files[i],".tif")
                      ri <- read_stars (file_name)
                      
                      coords <- haul_env [, c("lon", "lat")] %>% as.matrix()
                      
                      output <- stars::st_extract (ri,
                                            at=coords)
                      
                      return (output)
                    }

parallel::stopCluster(cl = my.cluster)

ben_data <- do.call ("rbind",ben_data) %>% t()

colnames (ben_data) <- files

#merge Morley data and Benthic data/Marta
haul_env <- cbind.data.frame (haul_env, ben_data)



#Var selection
var_int <- haul_env[,c("rugosity", "GRAINSIZE", "depth",
                       "monthly_mean_temp_bott", "monthly_mean_temp_surf",
                       "monthly_mean_salt_bott", "monthly_mean_salt_surf",
                       "monthly_mean_Po2_bott",
                       "annual_max_temp_surf", "annual_mean_temp_surf", "annual_min_temp_surf",
                       "annual_max_temp_bott", "annual_mean_temp_bott", "annual_min_temp_bott",
                       "annual_max_salt_surf", "annual_mean_salt_surf", "annual_min_salt_surf",
                       "annual_max_salt_bott", "annual_mean_salt_bott", "annual_min_salt_bott",
                       "annual_max_Po2_bott", "annual_mean_Po2_bott", "annual_min_Po2_bott")]

#save 
save (haul_env,wide_dat, 
      file = "data/processed/db_for_HMSC.RData")

