
#
# The following script generates Figures S2 and S3 for 
#
# Ecological effects and ecosystem shifts caused by mass mortality events on early life stages of fish
# Erik Olsen*, Cecilie H. Eide, Ina Nilsen, Holly A. Perryman, Frode Vikebø
# Submitted to Journal:Frontiers in Marine Science
# Specialty Section:Global Change and the Future Ocean
# Article type:Original Research Article
#

# This script was developed by Holly Perryman

#
#
# Create violin plots of Rn and Sn
#
#

# Load libraries
# not all libraries may be necessary, I copied list from another scritp file
library(sm)
library(ncdf4)
library(beepr)
library(atlantistools)
library(devtools)
library(rbgm)
library(shiny)
library(plyr)
library(dplyr)
library(DT)
library(ggplot2);library(grid);library(gridExtra);library(gtable)
library(stringr)
library(tidyr)
library(squash)
library(reshape2)
library(colorspace)
library(shiny)
library(shinyjs)
library(fmsb)
library(RColorBrewer)

# load self-made functions
# below is a function I mad for extracting information from the NC file
# it depends on functions from atlantistools
ExtractFromNCatltools <- function(fids, # path to fid csv file
                                  init, # path to initialization nc file 
                                  prm_biol, # path to biological prm file
                                  plgn, # path to atlantis input for bgm map
                                  nc_out, # path to output nc file
                                  prm_run, # path to run prm file
                                  fids_age, # selected age-structure fids
                                  fids_pool){ # selected biomass pool fids
  # #
  # # DEBUGGING
  # #
  # fids <- "nordic_groups_v04.csv"
  # init <- "nordic_biol_v23.nc"
  # prm_biol <- "nordic_biol_incl_harv_v01.prm"
  # plgn <- "Nordic02.bgm"
  # nc_out <- paste(outputfolders[1],"/ina_results_001.nc",sep="")
  # prm_run <- "nordic_run_v01.prm"
  # fids_age <- grps
  # fids_pool <- c("Small_phytop")
  print("!!!!!")
  print("Please double check the bat file to make sure all of the provided paths correspond to those used in the simulation.")
  print("!!!!!")
  # get benthic groups from input nc file
  bps <- load_bps(fgs = file.path(fids), init = file.path(init))
  # get bio_conv from biological prm file
  bio_conv <- get_conv_mgnbiot(prm_biol = file.path(prm_biol))
  # get boundary boxes from bgm file
  bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(plgn)))
  # 1 - Get nc data to calc biomass for age structured groups:
  # Get Nums from nc file
  Nums <- load_nc(nc = file.path(nc_out),
                  bps = bps,
                  fgs = file.path(fids),
                  prm_run = file.path(prm_run),
                  bboxes = bboxes,
                  select_groups = fids_age,
                  select_variable = "Nums")
  # get StructN from nc file
  sn  <- load_nc(nc = file.path(nc_out),
                 bps = bps,
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run),
                 bboxes = bboxes,
                 select_groups = fids_age,
                 select_variable = "StructN")
  # get ResN from nc file
  rn  <- load_nc(nc = file.path(nc_out),
                 bps = bps,
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run),
                 bboxes = bboxes,
                 select_groups = fids_age,
                 select_variable = "ResN")
  # 2 - Get nc data to calc biomass for biomass pooled groups:
  # get N from nc file
  n   <- load_nc(nc = file.path(nc_out), 
                 bps = bps, 
                 fgs = file.path(fids),
                 prm_run = file.path(prm_run), 
                 bboxes = bboxes,
                 select_groups = fids_pool, 
                 select_variable = "N")
    # get physics, specifically vol and dz, from nc file
  vol <- load_nc_physics(nc = file.path(nc_out),
                         prm_run = file.path(prm_run), 
                         bboxes = bboxes, 
                         aggregate_layers = FALSE,
                         select_physics = c("volume", "dz"))
  # CALC THE BIOMASS OF FIDS ACROSS SPACE
  print("Calculating group(s) biomass using calculate_biomass_spatial from atlantistools")
  B <- calculate_biomass_spatial(nums = Nums, 
                                 sn = sn, 
                                 rn = rn, 
                                 n = n, 
                                 vol_dz = vol, 
                                 bio_conv = bio_conv, 
                                 bps = bps)
  print("Calculating complete.")
  print("Merging output into one dataframe. This may take a couple minutes...")
  #Return
  list(B,rn,sn,Nums)
}
#

#
# Initialize 
#
# time steps between save events (see run prm file)
toutinc <- 73 
# mortality event (yr)
eventyr <- 62
# functional group name and info
FunctionalGroupInfo <- read.table("nordic_groups_v04.csv",header=TRUE,sep=",") 
FunctionalGroupInfo <- FunctionalGroupInfo[which(FunctionalGroupInfo$IsTurnedOn == 1),] 
FunctionalGroupInfo$OutID <- c(1:nrow(FunctionalGroupInfo)) 
grps <- c("North_atl_cod", "Haddock", "Norwegian_ssh")
grpids <- FunctionalGroupInfo[which(FunctionalGroupInfo$Name %in% grps),'OutID']
grps <- as.character(FunctionalGroupInfo[which(FunctionalGroupInfo$OutID %in% grpids),'Name']) # MAKE SURE grps IS IN THE CORRECT ORDER
grplgnms <- as.character(FunctionalGroupInfo[which(FunctionalGroupInfo$OutID %in% grpids),'Long.Name']) 
grpshnms <- as.character(FunctionalGroupInfo[which(FunctionalGroupInfo$OutID %in% grpids),'Code']) 
# upload bgm file
File_bgm <- "Nordic02.bgm" 
# define active boxes within modeling domain 
bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(File_bgm)))
#
# begin making plots
# define locations with simulation outputs 
outputfolders <- list.files()[which(substr(list.files(),1,6) == "output")]
# ID outfiles for processing
runs <- c(2,3,4,5,6,7,8,9,10,11,12,13,14,1)
# set corhort cut off
ageclcutoff <- 1
# define locations for nc data
outlist <- list()
outlistNC <- list()
# Loop
# for each file in runs
# 1) extract data from the nc file
# 2) sum across polygon and depth
# 3) apply ageclcutoff
# 4) remove data from burnin
# 5) loop through species in grps to get data:
#    RN and SN; annual value, annual min, annual max
for (i in 1:length(runs)){
  # Get info from NC file:
  nc <- ExtractFromNCatltools("nordic_groups_v04.csv",
                              "nordic_biol_v23.nc",
                              "nordic_biol_incl_harv_v01.prm",
                              "Nordic02.bgm",
                              paste(outputfolders[runs[i]],"/ina_results_001.nc",sep=""),
                              "nordic_run_v01.prm",
                              grps, 
                              c("Small_phytop")) # NEED TO DECLARE AT LEAST ONE B POOL GROUP
  # https://davetang.org/muse/2016/10/13/using-dplyr-aggregate-r/
  # sum across polygon and depth: nc[[2]]
  temp <- as.data.frame(group_by(nc[[2]], species, agecl, time) %>% summarise(RN = sum(atoutput)))
  outlist[[i]] <- temp
  # apply ageclcutoff
  temp2 <- temp[which(temp$agecl == ageclcutoff),] # ONLY PULL AGECL "agecutoff" DATA FOR NOW
  # remove data from burnin 
  temp2 <- temp2[which(temp2$time >= eventyr),]
  # remove temp and repeat for nc[[3]]
  remove(temp)
  # sum across polygon and depth: nc[[3]]
  temp <- as.data.frame(group_by(nc[[3]], species, agecl, time) %>% summarise(SN = sum(atoutput)))
  # apply ageclcutoff
  temp3 <- temp[which(temp$agecl == ageclcutoff),] # ONLY PULL AGECL "agecutoff" DATA FOR NOW
  # remove data from burnin 
  temp3 <- temp3[which(temp3$time >= eventyr),]
  # discard excess
  remove(temp,nc)
  # prepare for grps loop
  # get save intervals  
  printxintervals <- seq(1,length(unique(temp2$time)),365/toutinc) 
  # create save df
  tempdf <- data.frame(species=character(),
                       time=double(), 
                       Amax=double(),
                       Amin=double(),
                       Ax=double(),
                       metric=character())
  # Loop
  # go through species in grps to get data:
  #    RN and SN; annual value, annual min, annual max 
  for (j in 1:length(grps)) {
    # pull data
    temp2a <- temp2[which(temp2$species == grplgnms[j]),] 
    temp3a <- temp3[which(temp3$species == grplgnms[j]),] 
    # get maxs for temp2a
    annualmaxstemp <- lapply(1:(length(printxintervals)-1), function(x){
      subtemp <- temp2a[which(temp2a$time %in% (((c(1:(365/toutinc))+((365/toutinc)*(x-1)))/(365/toutinc))+eventyr)),]
      max(subtemp$RN / as.numeric(temp2a[which(temp2a$time == eventyr),"RN"]))
    })
    annualmax2 <- do.call("rbind", annualmaxstemp)
    # get maxs for temp3a
    annualmaxstemp <- lapply(1:(length(printxintervals)-1), function(x){
      subtemp <- temp3a[which(temp3a$time %in% (((c(1:(365/toutinc))+((365/toutinc)*(x-1)))/(365/toutinc))+eventyr)),]
      max(subtemp$SN / as.numeric(temp3a[which(temp3a$time == eventyr),"SN"]))
    })
    annualmax3 <- do.call("rbind", annualmaxstemp)
    # discard excess
    remove(annualmaxstemp)
    # get mins for temp2a
    annualminstemp <- lapply(1:(length(printxintervals)-1), function(x){
      subtemp <- temp2a[which(temp2a$time %in% (((c(1:(365/toutinc))+((365/toutinc)*(x-1)))/(365/toutinc))+eventyr)),]
      min(subtemp$RN / as.numeric(temp2a[which(temp2a$time == eventyr),"RN"]))
    })
    annualmin2 <- do.call("rbind", annualminstemp)
    # get mins for temp3a
    annualminstemp <- lapply(1:(length(printxintervals)-1), function(x){
      subtemp <- temp3a[which(temp3a$time %in% (((c(1:(365/toutinc))+((365/toutinc)*(x-1)))/(365/toutinc))+eventyr)),]
      min(subtemp$SN / as.numeric(temp3a[which(temp3a$time == eventyr),"SN"]))
    })
    annualmin3 <- do.call("rbind", annualminstemp)
    # discard excess
    remove(annualminstemp)
    # get annual values for the RN/R0 ratio
    annual2 <- temp2a[printxintervals,"RN"] / temp2a[which(temp2a$time == eventyr),"RN"]
    annual3 <- temp3a[printxintervals,"SN"] / temp3a[which(temp3a$time == eventyr),"SN"]
    # save RN
    tempdf2 <- data.frame(species = temp2a[printxintervals,"species"],
                         time = temp2a[printxintervals,"time"],
                         Amax = c(1,annualmax2),
                         Amin = c(1,annualmin2),
                         Ax = annual2,
                         metric = rep("RN",length(annual2)))
    tempdf <- rbind(tempdf,tempdf2)
    # save SN
    tempdf3 <- data.frame(species = temp3a[printxintervals,"species"],
                          time = temp3a[printxintervals,"time"],
                          Amax = c(1,annualmax3),
                          Amin = c(1,annualmin3),
                          Ax = annual3,
                          metric = rep("SN",length(annual3)))
    tempdf <- rbind(tempdf,tempdf3)
    # discard excess
    remove(tempdf3, tempdf2, annual3, annual2, annualmin3, annualmin2, annualmax3, annualmax2, temp3a, temp2a)
  }; remove(j)
  # SAVE 
  outlistNC[[i]] <- tempdf
  # discard excess
  remove(tempdf, printxintervals, temp2, temp3)
}; remove(i)
#
# Plot 
# make violin plots with last 10 years of simulation data
# prep data
# make flag identifying group data being plotted
i <- 3
# make vector of names for plots
plotnames <- c("Haddock", "Herring", "Cod")
# T/F
# make a binomial tag for identifying if we are plotting SN (True) or RN (False)
plotRN <- F
if (plotRN) {
  plotmetric <- "RN"
  ploty <-" RNt/RN0"
} else {
  plotmetric <- "SN"
  ploty <-" SNt/SN0"
}
# get data
forplot2 <- outlistNC[[1]][which(outlistNC[[1]]$species == grplgnms[i] & outlistNC[[1]]$metric == plotmetric),]
forplot3 <- outlistNC[[2]][which(outlistNC[[2]]$species == grplgnms[i] & outlistNC[[2]]$metric == plotmetric),]
forplot4 <- outlistNC[[3]][which(outlistNC[[3]]$species == grplgnms[i] & outlistNC[[3]]$metric == plotmetric),]
forplot5 <- outlistNC[[4]][which(outlistNC[[4]]$species == grplgnms[i] & outlistNC[[4]]$metric == plotmetric),]
forplot6 <- outlistNC[[5]][which(outlistNC[[5]]$species == grplgnms[i] & outlistNC[[5]]$metric == plotmetric),]
forplot7 <- outlistNC[[6]][which(outlistNC[[6]]$species == grplgnms[i] & outlistNC[[6]]$metric == plotmetric),]
forplot8 <- outlistNC[[7]][which(outlistNC[[7]]$species == grplgnms[i] & outlistNC[[7]]$metric == plotmetric),]
forplot9 <- outlistNC[[8]][which(outlistNC[[8]]$species == grplgnms[i] & outlistNC[[8]]$metric == plotmetric),]
forplot10 <- outlistNC[[9]][which(outlistNC[[9]]$species == grplgnms[i] & outlistNC[[9]]$metric == plotmetric),]
forplot11 <- outlistNC[[10]][which(outlistNC[[10]]$species == grplgnms[i] & outlistNC[[10]]$metric == plotmetric),]
forplot12 <- outlistNC[[11]][which(outlistNC[[11]]$species == grplgnms[i] & outlistNC[[11]]$metric == plotmetric),]
forplot13 <- outlistNC[[12]][which(outlistNC[[12]]$species == grplgnms[i] & outlistNC[[12]]$metric == plotmetric),]
forplot17 <- outlistNC[[13]][which(outlistNC[[13]]$species == grplgnms[i] & outlistNC[[13]]$metric == plotmetric),]
# merge data
forplotmerge <- Reduce(function(x, y) merge(x, y, by = c("species", "time"), all=TRUE), 
                       list(forplot2[,c("species", "time", "Ax")],forplot3[,c("species", "time", "Ax")],forplot4[,c("species", "time", "Ax")],
                            forplot5[,c("species", "time", "Ax")],forplot6[,c("species", "time", "Ax")],forplot7[,c("species", "time", "Ax")],
                            forplot8[,c("species", "time", "Ax")],forplot9[,c("species", "time", "Ax")],forplot10[,c("species", "time", "Ax")],
                            forplot11[,c("species", "time", "Ax")],forplot12[,c("species", "time", "Ax")], forplot13[,c("species", "time", "Ax")], 
                            forplot17[,c("species", "time", "Ax")]))#,forplot0[,c("species", "time", "Ax")]))
# changer header
names(forplotmerge) <- c("species", "time", "CO10", "CO50", "CO90", "HE10", "HE50", "HE90", "HA10", "HA50", "HA90", "MX10", "MX50", "MX90", "BASE")
# identifier saying last "stopyr" years of data will be used for plots 
stopyr <- 10
# T/F 
# is the data coming from the end of the data (T) 
# or the beginning of the data (F)
dataend <- T
# get data for plot
if(dataend){
  #last stopyr years of the data:
  cutoff <- tail(sort(forplotmerge$time),10)
  dataforVplot <- forplotmerge[which(forplotmerge$time %in% cutoff),]
} else {
  # first stopyr years of the data:
  cutoff <- seq(63,(63+(stopyr - 1)),1)
  dataforVplot <- forplotmerge[which(forplotmerge$time %in% cutoff),]
}
# convert from wide to long for plot
dataforVplot <- gather(dataforVplot, scenario, ratio, CO10:BASE, factor_key=TRUE)
# make plot
# 1, save plot
if(dataend){
  #last stopyr years of the data:
  tiff(paste(plotmetric,grpshnms[i],stopyr,"end.tiff",sep=""), units="in", width=8, height=2.5, res=300)
  #pdf(paste(plotmetric,grpshnms[i],stopyr,"end.pdf",sep=""),width=8, height=2.5)
} else {
  # first stopyr years of the data:
  tiff(paste(plotmetric,grpshnms[i],stopyr,".tiff",sep=""), units="in", width=8, height=2.5, res=300)
  #pdf(paste(plotmetric,grpshnms[i],stopyr,".pdf",sep=""),width=8, height=2.5)
}
# 2, generate plot
ggplot(dataforVplot, aes(x=scenario, y=ratio, fill=scenario)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  labs(title=" ",x="Model scenario", y = paste(plotnames[i],ploty,sep="")) +
  #theme_minimal() + 
  theme_classic() + 
  #scale_fill_brewer(palette="Blues") +
  scale_fill_manual(values=c(cCO10,cCO50,cCO90,cHA10,cHA50,cHA90,cHE10,cHE50,cHE90,cMX10,cMX50,cMX90,cBASE)) +
  theme(legend.position="none") # Remove legend
# 3, close and save plot
dev.off()
