library(ncdf4)
library(tidyverse)
library(RColorBrewer)
library(extrafont)
loadfonts(device="win") 

# Set working directory
workdir="C:/Users/Ina Nilsen/Atlantis/Oil_spill/" #Finner frem til mappen med filene

# Select forlders 
folder0 = "output_00/"
folder1 = "output_10/"
folder2 = "output_11/"
folder3 = "output_12/"

# Import data from folders
nc_0=nc_open(paste(workdir,folder0,"ina_results_001.nc",sep=""))
nc_1=nc_open(paste(workdir,folder1,"ina_results_001.nc",sep=""))
nc_2=nc_open(paste(workdir,folder2,"ina_results_001.nc",sep=""))
nc_3=nc_open(paste(workdir,folder3,"ina_results_001.nc",sep=""))

# Create timeframe
time=ncvar_get(nc_0,"t")

# Choose species
sel.Fish=c("North_atl_cod","Haddock","Norwegian_ssh") 

############    Make dataframes of runs   ##################################################
# This is done 4 times for the 4 runs selected
# The base run (Biomass_00) is added as a column in the dataframes of the modified runs
# This was done because these runs (MX10, MX50, MX90) will be compared to the base run

# BASE (Run 00)
Biomass_00 <- lapply(sel.Fish, function(f) {          # All age-structured groups have
  if (f == "Capelin") {                               # 10 age classes, except:
    ycl = c(1:5)                                      # Capelin, snow crab and sperm whale
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){ 
    nums = ncvar_get(nc_0, paste(f, y, "_Nums", sep = ""))      # Extract numbers from nc-file
    struct = ncvar_get(nc_0,paste(f, y, "_StructN", sep = ""))  # Extract structural weight
    res = ncvar_get(nc_0,paste(f, y, "_ResN", sep = ""))        # Extract reserved weight
    biom=(struct+res)*nums*5.7*20/1e9                           # Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))                              #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>%                                        
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*73, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=61) %>%      # Filter out years before mortality event
  rename(Biomass00=Biomass) # Rename biomass columne to add to other data frames

# MX10 (run with 10 % mortality added)
Biomass_10 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_1, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_1,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_1,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*73, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=61) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

# MX50 (run with 50 % mortality added) 
Biomass_50 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_2, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_2,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_2,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*73, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=61) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

# MX90 (run with 90 % mortality added)
Biomass_90 <- lapply(sel.Fish, function(f) {
  if (f == "Capelin") {
    ycl = c(1:5)
  } else if (f == "Snow_crab") {
    ycl = c(1:6)
  } else if (f == "Sperm_whale") {
    ycl = c(1:8)
  } else{
    ycl = c(1:10)
  }
  sapply(ycl, function(y){
    nums = ncvar_get(nc_3, paste(f, y, "_Nums", sep = ""))
    struct = ncvar_get(nc_3,paste(f, y, "_StructN", sep = ""))
    res = ncvar_get(nc_3,paste(f, y, "_ResN", sep = ""))
    biom=(struct+res)*nums*5.7*20/1e9 #Convert from Nmg to tonnes
    biom_ts=colSums(colSums(biom))    #Sums layers and polygons
    biom_ts
    
  }) %>% as.tibble() %>% 
    mutate(Species = f) %>% 
    rowid_to_column("time") %>% 
    mutate(time1=time-1) %>% 
    mutate(day = time1*73, Year = day/365) %>% 
    select(-time, -time1)
  
}) %>% bind_rows() %>% 
  gather(key = Age_class, value = Biomass, -Species, -Year) %>% 
  mutate(Age_class = gsub("V","", Age_class), 
         Age_class = as.integer(Age_class)) %>% 
  filter(Year>=61) %>%  # Filter out years before mortality event
  left_join(Biomass_00) # Add column from base run

############    Compare runs with base run  #################################################
# Make combined badaframe of all runs
Biom_all <- bind_rows(`MX10` = Biomass_10,   # (run with 10 % mortality added)
                      `MX50` = Biomass_50,   # (run with 50 % mortality added)
                      `MX90` = Biomass_90,   # (run with 90 % mortality added)
                      .id="Run")%>%
            mutate(Season = Year %% 1) %>%
            filter(Season==0) %>%      # As data is printed 5 times per year we choose the first timestep
            mutate(Diff=(Biomass-Biomass00)/Biomass00*100) #Difference (%) between base and modified runs

#Add black line representing total change of all age classes
Biom_tot <- Biom_all %>% 
  group_by(Run, Species, Year) %>% 
  summarise(Biomass00_sum=sum(Biomass00),
            Biomass_sum=sum(Biomass)) 

Bar_plot <- left_join(Biom_tot,Biom_all) %>% 
  mutate(Change=Diff*(Biomass00/Biomass00_sum)) 

Line_plot <-Biom_tot %>% 
  mutate(Tot_change=(Biomass_sum-Biomass00_sum)/Biomass00_sum*100)


############    Make plot   #################################################################

Plot <- Bar_plot %>% filter(Species %in% c( "North_atl_cod")) %>%    # Choose species
        filter( Age_class==1 | Age_class==2 | Age_class==3 | Age_class==4 | Age_class==5 | Age_class==6 |Age_class==7| Age_class==8 | Age_class==9 | Age_class==10) %>% 
        ggplot(aes(x=Year, y=Change, fill=as.factor(Age_class))) +
          scale_fill_brewer(palette = "Paired")+  
          geom_col() +
          geom_hline(yintercept = 0, size=1) + # Add base line at 0
          geom_line(aes(x=Year, y=Tot_change),Line_plot %>%   # Add line to plot
             filter(Species %in% c( "North_atl_cod")), inherit.aes = FALSE , size =1) +
          labs(x="Time (Year post MME)", y="Change in biomass relative \n to baseline (%)", title="Northeast Arctic Cod", fill="Age class") +
          theme_minimal()+
          theme(text=element_text(family="Calibri",size=25), plot.title = element_text(size=30,face="bold")) +
          facet_wrap(~Run)    # Plot all runs seperatly  (same scale on axis)      
       Plot 

ggsave("Figure_3.tiff", units="in", height=6, width=16, dpi=300 )

