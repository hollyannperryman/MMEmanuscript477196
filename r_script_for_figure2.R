
library(tidyverse)
library(RColorBrewer)
library(extrafont)
loadfonts(device="win") 

# Set working directory
workdir="C:/Users/Ina Nilsen/Atlantis/Oil_spill/" 

# List of species codes linked to full species name 
nordic <- read.table("C:/Users/Ina Nilsen/Atlantis/Oil_spill/nordic_groups_Oil.csv", dec=".", sep=",", header=TRUE) %>% 
  select(Code,Species_name) 


# Select folders with output data
folder0  = "output_00"
folder1  = "output_01"
folder2  = "output_02"
folder3  = "output_03"
folder4  = "output_04"
folder5  = "output_05"
folder6  = "output_06"
folder7  = "output_07"
folder8  = "output_08"
folder9  = "output_09"
folder10 = "output_10"
folder11 = "output_11"
folder12 = "output_12"

# Import data from folders
run_00 <- read.table(paste(workdir, folder0, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_01 <- read.table(paste(workdir, folder1, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_02 <- read.table(paste(workdir, folder2, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_03 <- read.table(paste(workdir, folder3, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_04 <- read.table(paste(workdir, folder4, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_05 <- read.table(paste(workdir, folder5, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_06 <- read.table(paste(workdir, folder6, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_07 <- read.table(paste(workdir, folder7, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_08 <- read.table(paste(workdir, folder8, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_09 <- read.table(paste(workdir, folder9, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_10 <- read.table(paste(workdir, folder10, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_11 <- read.table(paste(workdir, folder11, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)
run_12 <- read.table(paste(workdir, folder12, "/ina_results_001BiomIndx.txt", sep=""), header=TRUE)


# Make dataframe
biom.df <-bind_rows(`BASE` = run_00,
                    `CO10` = run_01,
                    `CO50` = run_02,
                    `CO90` = run_03,
                    `HE10` = run_04,
                    `HE50` = run_05,
                    `HE90` = run_06,
                    `HA10` = run_07,
                    `HA50` = run_08,
                    `HA90` = run_09,
                    `MX10` = run_10,
                    `MX50` = run_11,
                    `MX90` = run_12,
                  .id="Runs")%>%
  gather(key=Species_code, value=Biomass, -Time, -Runs) %>%
  rename(Day =Time) %>% 
  mutate(Year=Day/365) %>% 
  mutate(Season = Year %% 1) %>%
  filter(Season==0)  %>%
  filter(Year>58) %>% 
  filter(Species_code %in% c("HAD", "SSH", "NCO")) %>%   # Choose species
  left_join(nordic, by=c("Species_code"="Code"))         # Add complete name of species


# Set order of runs, species and colours for plot
run.order <- c("CO10","CO50","CO90","HA10","HA50","HA90","HE10","HE50","HE90","MX10","MX50","MX90","BASE")
sp.order <- c("Northeast Arctic Cod","Northeast Arctic Haddock","Norwegian Spring Spawning Herring")
colours <- c("#a6cee3","#80b1d3","#1f78b4","#b3de69","#33a02c","#006837","#fdb462","#ff7f00","#b15928","#cab2d6","#bc80bd","#6a3d9a","#000000")


# Make plot
plot <- biom.df %>% 
  filter(Species_code %in% c("SSH","HAD","NCO")) %>% 
    mutate(Runs = factor(Runs, levels = run.order)) %>%
    arrange(Runs) %>% 
    mutate(Species_name = factor(Species_name, levels = sp.order)) %>%
    arrange(Species_name)  %>% 
    mutate(Biomass=Biomass/250000) %>%                                   # Create index of biomass
    ggplot(aes(x=Year, y=Biomass, color=Runs)) +  
       geom_line(lwd=1.5) +  
       scale_colour_manual(values=colours)+  
       labs(x="Year", y="Biomass index", title="", colour="Scenarios") +
       theme_minimal()+
       theme(text=element_text(family="Calibri",size=25), plot.title = element_text(size=24)) +
       facet_wrap(~Species_name, scales = "free_y") # Plot species seperatly (different scale on y-axis)
plot


ggsave("Figure_2.tiff", units="in", height=5.5, width=17, dpi=300 )



