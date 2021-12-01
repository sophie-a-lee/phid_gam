#######################################################
####                                               ####
#### Load data and packages for PHID mgcv tutorial ####
####                                               ####
#######################################################


#### Load packages ####
install.packages("pacman")

# p_load installs packages and loads them in the same step
pacman::p_load(tidyverse, data.table, mgcv, sf, mgcViz, INLA,
               spdep, geobr, cowplot, pROC)



#### Load Brazilian shapefile ####
## Due to time constraints, select municipality-level data from Rio de Janeiro
# Load data from geobr package, taken from IBGE
shp <- read_municipality(code_muni = "RJ") %>% 
  rename(municip_code_ibge = code_muni) %>% 
  arrange(municip_code_ibge) %>%
  # Add municipaliity index for INLA random effects
  mutate(municip_index = 1:nrow(.)) 



#### Load epidemiological and socioeconomic data ####
## Use data from RJ 2010 - 2020
df <- fread("data/dengue_rj.csv") %>% 
  filter(year >= 2010) %>% 
  # Convert year to a factor to include as a fixed effect
  mutate(fyear = factor(year),
         # Create a time index value for INLA random effects
         time = year - 2009,
         # Create a categorical level of influence variables with labels
         regic18 = factor(level18_num, levels = 1:5,
                          labels = c("Metropolis",
                                     "Regional capital",
                                     "Sub-regional centre",
                                     "Zone centre",
                                     "Local centre")))
