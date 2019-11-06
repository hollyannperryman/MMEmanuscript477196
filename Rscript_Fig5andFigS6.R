
#
# The following script generates Figure 5 for 
#
# Olsen, E., Eide, C.H., Nilsen, I., Perryman, H.A. and 
# Vikebø, F., 2019. Ecological effects and ecosystem 
# shifts caused by mass mortality events on early life stages
# of fish. Frontiers in Marine Science, 6, p.669.

# This script was developed by Holly Perryman based on 
# discussions with Cecilie Hansen, Ina Nelson and Erik Olsen
# and inspired by indicators listed in:
# Olsen, E., Kaplan, I.C., Ainsworth, C., Fay, G., Gaichas, S., Gamble, R., 
# Girardin, R., Eide, C.H., Ihde, T.F., Morzaria-Luna, H.N. and Johnson,
# K.F., 2018. Ocean futures under ocean acidification, marine protection, 
# and changing fishing pressures explored using a worldwide suite of ecosystem 
# models. Frontiers in Marine Science, 5, p.64.

#
#
# Create spider plots of ecological and fisheries indicators
#
#

# Upload libraries
# not all of the following is required for this script, 
# I simply used all the libraries called from another script file
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
library(ReactiveAtlantis)

# upload self made function 
# this function was made by Holly based on discussions with Ina and Cecilie
get_PrimaryProduction <- function(fids, 
                                  PRODnc, 
                                  BGMfile){
  # FOR DEBUGGING:
  # PRODnc<-paste(outputfolders[1],"/ina_results_001PROD.nc",sep="")
  # fids <- sp
  # BGMfile <- File_bgm
  #
  # get info from PRODnc file
  nc <- nc_open(PRODnc)
  # get info from BGMfile
  bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(BGMfile)))
  boxdetails <- bgmfile(BGMfile)
  # create matrix for saving output  
  forout <- matrix(NA, 
                   nrow = dim(ncvar_get(nc,paste(as.character(FunctionalGroupInfo$Name[fids[1]]),"Prodn", sep ="")))[2], 
                   ncol = length(fids))
  # loop through groups (ie fids)
  for (spp in c(1:length(fids))) {
    # get data
    temp <- ncvar_get(nc,paste(as.character(FunctionalGroupInfo$Name[fids[spp]]),"Prodn", sep =""))
    # convert from (mg N m-3 d-1) to (tons d-1)  
    temp=temp*5.7*20/(1000*1000*1000) # Convert from (mg N m-3 d-1) to (tons m-3 d-1) 
    prod.tot <- matrix(NA,nrow=dim(temp)[1],ncol=dim(temp)[2])
    for(i in c(1:dim(temp)[1])){ # Convert from (tons m-3 d-1) to (tons d-1) 
      prod.tot[i,]=temp[i,] * as.numeric(boxdetails$boxes[i,4]) * -as.numeric(boxdetails$boxes[i,3])
    } 
    # sum across polygons 
    # !!!! you will need to remove boundary polygons (bboxes)
    # !!!! you need to correct to match with polygon output location 
    temp <- temp[!(c(1:dim(temp)[1]) %in% (bboxes + 1)),] 
    forout[,spp] <- colSums(temp)
  }
  PPout <- rowSums(forout)
  return(PPout)
}
#

# Prep for the loop calculating indicators 
FunctionalGroupInfo <- read.table("nordic_groups_v04.csv",header=TRUE,sep=",") # UPLOAD GROUP INFORMATION USED IN MODEL RUN
FunctionalGroupInfo <- FunctionalGroupInfo[which(FunctionalGroupInfo$IsTurnedOn == 1),] # REMOVE FIDS THAT ARE TURNED OFF
FunctionalGroupInfo$OutID <- c(1:nrow(FunctionalGroupInfo)) # SAVE OUTPUT LOCATION
File_bgm <- "Nordic02.bgm" 
bboxes <- get_boundary(boxinfo = load_box(bgm = file.path(File_bgm)))
toutinc <- 73 # DECLARING LENGTH BETWEEN TIMESTEPS (SEE THE RUN PRM FILE)
#
# Output folders
outputfolders <- list.files()[which(substr(list.files(),1,6) == "output")]
tier <- outputfolders[-1] 
# subset the outputfolders if desired
#tier <- outputfolders[c(11,12,13,14)]
#
# titles of scenarios (i.e., output files in tier)
tier.t=c("CO10","CO50","CO90",
         "HE10","HE50","HE90",
         "HA10","HA50","HA90",
         "MX10","MX50","MX90","BASE")
#tier.t=c("MX10","MX50","MX90","BASE")
#
# Basic info of functional groups for computing indicators 
basic=read.table("NOBA_BasicInfo_HAP.csv",header=T,sep=';',na.strings = "NA",stringsAsFactors = F)
#
# create list of the indicators for calculations 
indicators_e <- c("PelB.PP", # ECOLOGICAL INDICATORS
                  "Bio.PP",
                  "MTL.B",
                  "PropPred",
                  "DemT.PelT",
                  "Dem.Pel",
                  "Dem.PP",
                  
                  "TotPC", # FISHERIES INDICATORS 
                  "TotC",
                  "MTL.C",
                  "FishExRt",
                  "ExRt",
                  "Val",
                  "TotFC",
                  "TotDC")
# create array for collecting indicator computations for radar plot
forplot_radarc=array(0,dim=c((length(tier)+2),length(indicators_e)))



#
#
# RADAR PLOTS 
# AVERAGES I) 5 YEARS, II) 10 YEARS, AND III) 20 YEARS POST MME
#
#

# Define the temporal span of data going into calculations
t.span.out <- 20

# seasonality flag (T/F):
# (T) - you are collecting all of the seasonal data reported in the out files 
# (F) - you are collecting only the data reported at the end of the year (i.e., every 365 time steps) 
t.season <- F # !!! I DO NOT HAVE IT PROGRAMMED FOR t.season <- T YET !!! 

# The time step just before the infliction of the mortality event 
eventyr <- 61
#
# Loop through each model scenario out to collect data for radar plot
# (this method works but is inefficient)
for(i in 1:length(tier)){ 
  # get output data
  biom=read.table(paste(tier[i],'/ina_results_001BiomIndx.txt',sep=''),header=T) # you will need to correct this if file in different directory (wd, see above)
  catch=read.table(paste(tier[i],'/ina_results_001Catch.txt',sep=''),header=T) # you will need to correct this if file in different directory (wd, see above)
  # Subset temporal span of data going into indicator computation 
  # create markers indicating the start and end of the temporal range
  last.yr <- floor(biom$Time[length(biom$Time)]/365) 
  # sometimes a model stops before it reaches the next year 
  # to make sure you are collecting consistant points in the output
  # you need to correct for this by grabing the last annual time step
  # and not some time step througout the year
  if(t.season){
    # !!! NOT PROGRAMED - DID YOU DO !!! 
  } else {
    # !!!! the above line is assuming you are collecting data from the last (t.span.out years) of the simulation 
    data.grab <- c((eventyr + 1):(eventyr + t.span.out))
    # !!!! the above line is assuming you are collecting data immediately following the mass mortality event 
  }
  # subset to get data fro caluculations 
  b.data.for.calc <- biom[which(biom$Time %in% (data.grab*365)),] # you need to correct data.grab so it is off the same scale as the txt out files
  c.data.for.calc <- catch[which(catch$Time %in% (data.grab*365)),] # you need to correct data.grab so it is off the same scale as the txt out files
  # create index for generating dataframe for radar plot
  jcnt=1 
  # loop through identifiers in var.sp to collect info for indicator computations 
  for(j in 1:length(indicators_e)){
    # loop through the indicators listed in indicators_e and compute means for radar plot
    if(indicators_e[j] == "PelB.PP"){
      # get total pelagic B
      guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
      # get primary production PP
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                     paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      # temporal subset of PPout
      # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / tempPP
      # purge old
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "Bio.PP"){ 
      # get total B
      biom_var1=rowSums(b.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
      # get primary production PP
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                     paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      # temporal subset of PPout
      # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / tempPP
      # purge old
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "MTL.B"){ 
      # create place holder
      b.data.for.calc$MTL <- NA
      # need to loop through each yr and fid
      for (tempyr in c(1:nrow(b.data.for.calc))) {
        TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
        biom.var <- b.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
        numtr <- biom.var * TLs
        b.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
      }
      # set tempvalue for calc avg of indicator 
      tempvalue <- b.data.for.calc$MTL
      # purge old
      remove(TLs, biom.var, numtr)
    } else if(indicators_e[j] == "PropPred"){ 
      # get predatory fish B
      guild=which(basic$Atlantis.species.code=='IsPredatoryFish') # row number in basic that corresponds to IsPredatoryFish 
      sp=which(basic[guild,]==1)                                  # functional groups (columns) that meet IsPredatoryFish
      sp=sp-1                                                     # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                     # total biomass
      # get fish B
      guild=which(basic$Atlantis.species.code=='IsFish')     # row number in basic that corresponds to IsFish 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
      sp=sp-1                                                # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / biom_var2
      # purge old
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "DemT.PelT"){ 
      # get DemT B
      guild=which(basic$Atlantis.species.code=='IsDemersalTeleost') # row number in basic that corresponds to IsDemersalTeleost 
      sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
      sp=sp-1                                                       # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
      # get PelT B
      guild=which(basic$Atlantis.species.code=='IsPelagicTeleost') # row number in basic that corresponds to IsPelagicTeleost 
      sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
      sp=sp-1                                                      # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / biom_var2
      # remove 
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "Dem.Pel"){ 
      # get Dem B
      guild=which(basic$Atlantis.species.code=='IsDemersal') # row number in basic that corresponds to IsDemersalTeleost 
      sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
      sp=sp-1                                                       # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
      # get Pel B
      guild=which(basic$Atlantis.species.code=='IsPelagic') # row number in basic that corresponds to IsPelagicTeleost 
      sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
      sp=sp-1                                                      # correct functional group identifer 
      biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / biom_var2
      # remove 
      remove(guild,sp,biom_var1,biom_var2)
    } else if(indicators_e[j] == "Dem.PP"){ 
      # get total demersal B   IsDemersal
      guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
      # get primary production PP
      guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
      sp=sp-1                                                # correct functional group identifer 
      # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
      PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                     paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                     File_bgm)
      # temporal subset of PPout
      # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
      tempPP <- PPout[(biom$Time / 365) %in% (data.grab)]
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / tempPP
      # remove
      remove(biom_var1,tempPP,PPout,sp,guild)
    } else if(indicators_e[j] == "TotPC"){ 
      # get total pelagic B
      guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1
      # remove
      remove(biom_var1,sp,guild)
    } else if(indicators_e[j] == "TotC"){ 
      # get total B
      biom_var1=rowSums(c.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1
      # remove
      remove(biom_var1)
    } else if(indicators_e[j] == "MTL.C"){ 
      # create place holder
      c.data.for.calc$MTL <- NA
      # need to loop through each yr and fid
      for (tempyr in c(1:nrow(c.data.for.calc))) {
        TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
        biom.var <- c.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
        numtr <- biom.var * TLs
        c.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
      }
      # set tempvalue for calc avg of indicator 
      tempvalue <- c.data.for.calc$MTL
      # remove
      remove(TLs, biom.var, numtr)
    } else if(indicators_e[j] == "FishExRt"){ 
      # get targeted fish groups
      guild1=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
      guild2=which(basic$Atlantis.species.code=='IsFish')   # row number in basic that corresponds to IsTarget 
      sp1=which(basic[guild1,]==1)                          # functional groups (columns) 
      sp2=which(basic[guild2,]==1)                          # functional groups (columns) 
      sp <- sp1[sp1 %in% sp2]                               # targeted groups that are also fish
      sp=sp-1                                               # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])               # total catch
      biom_var2=rowSums(b.data.for.calc[,sp])               # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / biom_var2
      # remove
      remove(guild1, guild2, sp1, sp2, biom_var1,  biom_var2)
    } else if(indicators_e[j] == "ExRt"){ 
      # get targeted groups
      guild=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
      sp=which(basic[guild,]==1)                             # functional groups (columns) 
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total catch
      biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1 / biom_var2
      # remove
      remove(guild, sp, biom_var1, biom_var2)
    } else if(indicators_e[j] == "Val"){ 
      # get total fish C
      guild=which(basic$Atlantis.species.code=='USDollarsPerTon') # row number in basic that corresponds to USDollarsPerTon 
      sp=which(is.na(basic[guild,]) == F)                         # functional groups (columns) that that have USDollarsPerTon
      biom_prc <- basic[guild,sp[-1]]                            # prices
      sp=sp[-1]-1                                                 # correct functional group identifer 
      biom_totc=c.data.for.calc[,sp]                              # catch (ton) of fids 
      valprod <- t(t(as.matrix(biom_totc)) * as.numeric(biom_prc))
      biom_var1=rowSums(valprod)                    
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1
      # remove
      remove(guild, sp, biom_prc, biom_totc, valprod, biom_var1)
    } else if(indicators_e[j] == "TotFC"){ 
      # get total fish C
      guild=which(basic$Atlantis.species.code=='IsFish')  # row number in basic that corresponds to IsFish 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1
      # remove
      remove(biom_var1,sp,guild)
    } else if(indicators_e[j] == "TotDC"){ 
      # get total demersal C
      guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
      sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
      sp=sp-1                                                # correct functional group identifer 
      biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
      # set tempvalue for calc avg of indicator 
      tempvalue <- biom_var1
      # remove
      remove(biom_var1,sp,guild)
    } 
    # Save computed indicator 
    forplot_radarc[(i+2),j]=mean(tempvalue)
  } # end of j loop
}
#
# Prep data for plot
# compute relative values 
forplot <- t(t(forplot_radarc)/forplot_radarc[dim(forplot_radarc)[1],]) 
# create names
# !!! THIS NEEDS TO BE BASED ON THE ORDER IN forplot
colnames(forplot) <- c('Pel.bio/PP',
                       'Bio/PP',
                       'MTL.bio',
                       'Predfish.Prop',
                       'Dem/pel.fish',
                       'Dem/Pel',
                       'Dem.bio/PP',
                       
                       "Pel.C",
                       "Tot.C",
                       "MTL.C",
                       "Fish.exp.rt",
                       "Exp.rt",
                       "Val",
                       "Fish.C",
                       "Dem.C")
rownames(forplot)=c('max','min',tier.t)
# Create data frame for the radar plots 
forplot1 <- forplot[,c(1:7)]
forplot2 <- forplot[,c(8:15)]
# add min/max values for the first two rows
forplot1[1,] <- 1.1 #max(forplot,na.rm = T)
forplot1[2,] <- 0.9 #min(forplot[c(1,3:(length(tier)+2)),] ,na.rm = T)
forplot2[1,] <- 1.1 #max(forplot,na.rm = T)
forplot2[2,] <- 0.7 #min(forplot[c(1,3:(length(tier)+2)),] ,na.rm = T)
#
seqforplot1 <- seq(0.9,1.1,0.05)
seqforplot2 <- seq(0.7,1.1,0.1)
#
forplot1=data.frame(forplot1)
forplot2=data.frame(forplot2)
#
# Colors for plots
#http://research.stowers.org/mcm/efg/R/Color/Chart/ColorChart.pdf
cCO10<-"#BEC1D4" ;cCO50<-"#7D87B9" ;cCO90<-"#023FA5" 
cHA10<-"#F0F2D5" ;cHA50<-"#B7D0A0" ;cHA90<-"#57945C" 
cHE10<-"#F1DE81" ;cHE50<-"#EEAB65" ;cHE90<-"#B9534C" 
cMX10<-"#D8BFD8" ;cMX50<-"#BA55D3" ;cMX90<-"#800080" 
#
# MAKE PLOTS
tiff(paste("Radar",t.span.out,".tiff",sep = ""), width = 12, height = 4, units = "in", res = 300)
par(mfrow=c(1,2),xpd=TRUE)
par(#oma=c(0,0,0,0),
  mar=c(0,0,0,0),xpd=TRUE)
# Create radar plot of ecological indicators 
radarchart(forplot1,
           plwd=rep(2,13),
           plty=c(rep(1,9),rep(2,4)),
           pcol=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
           axistype=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seqforplot1, cglwd=0.8,
           vlcex=0.8)
par(mar=c(0,0,0,6),xpd=TRUE)
# Create radar plot of fisheries indicators 
radarchart(forplot2,
           plwd=rep(2,13),
           plty=c(rep(1,9),rep(2,4)),
           pcol=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
           axistype=1,
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seqforplot2, cglwd=0.8,
           vlcex=0.8)
legend("right",
       lty=c(rep(1,9),rep(2,4)),
       lwd = rep(2,13),
       col=c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey","black"),
       legend=tier.t,
       bty = "n",
       inset=c(-0.22,0))
dev.off()




#
#
# TIME SERIES PLOTS FOR SUPPLEMENTARY MATERIALS (FIG. S6)
# 
#

# We want to avoid plotting seasonal variability - that is not the focus here
# so first I need to define the vector of ts that we want to plot
data.grab <- (seq(eventyr,(eventyr + 50),1) * 365)

# We need to loop through output file and create dataframes for plotting
# rather than edit the loop for the previous radar plots
# I prefer just makeing a new loop
# That way I dont mess up the cod for the radar plots 

# define data frames for data
# ecological indicators:
data.for.plot.PelBPP <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.BioPP <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.MTLbio <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.PredFishrop <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.DemPelFish <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.DemPel <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.DemBioPP <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
# fisheries indicators:
data.for.plot.PelC <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.TotC <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.MTLC <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.FishExpRT <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.ExpRT <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.Val <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.FishC <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)
data.for.plot.DemC <- data.frame(Time = ((data.grab / 365) - eventyr),temp = 1)

# Loop through output folders
# select data
# calculatin indicators
# safe into defined data frames
for(i in 1:length(tier)){ 
  # get output data
  remove(biom,catch)
  biom=read.table(paste(tier[i],'/ina_results_001BiomIndx.txt',sep=''),header=T) # you will need to correct this if file in different directory (wd, see above)
  catch=read.table(paste(tier[i],'/ina_results_001Catch.txt',sep=''),header=T) # you will need to correct this if file in different directory (wd, see above)
  # drop data prior to MME (MME should be t 0 in plots)
  remove(b.data.for.calc,c.data.for.calc)
  b.data.for.calc <- biom[which(biom$Time %in% data.grab),] # you need to correct data.grab so it is off the same scale as the txt out files
  c.data.for.calc <- catch[which(catch$Time %in% data.grab),] # you need to correct data.grab so it is off the same scale as the txt out files
  # 
  #
  # Get Time Series data for PelB.PP
  #
  #
  # get total pelagic B
  guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
  # get primary production PP
  guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
  sp=sp-1                                                # correct functional group identifer 
  # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
  PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                 paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                 File_bgm)
  # temporal subset of PPout
  # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
  tempPP <- PPout[(biom$Time %in% data.grab)]
  # save data into data frame 
  data.for.plot.PelBPP$PelB.PP <- biom_var1 / tempPP
  names(data.for.plot.PelBPP)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,tempPP,PPout,sp,guild)
  # 
  #
  # Get Time Series data for Bio.PP
  #
  #
  # get total B
  biom_var1=rowSums(b.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
  # get primary production PP
  guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
  sp=sp-1                                                # correct functional group identifer 
  # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
  PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                 paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                 File_bgm)
  # temporal subset of PPout
  # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
  tempPP <- PPout[(biom$Time %in% data.grab)]
  # save data into data frame 
  data.for.plot.BioPP$data <- biom_var1 / tempPP
  names(data.for.plot.BioPP)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,tempPP,PPout,sp,guild)
  # 
  #
  # Get Time Series data for MTL.bio
  #
  #
  # create place holder
  b.data.for.calc$MTL <- NA
  # need to loop through each yr and fid
  for (tempyr in c(1:nrow(b.data.for.calc))) {
    TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
    biom.var <- b.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
    numtr <- biom.var * TLs
    b.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
  }
  # save data into data frame 
  data.for.plot.MTLbio$data <- b.data.for.calc$MTL
  names(data.for.plot.MTLbio)[i+2] <- tier[i]
  # purge old
  remove(TLs, biom.var, numtr)
  # 
  #
  # Get Time Series data for Predfish.Prop
  #
  #
  # get predatory fish B
  guild=which(basic$Atlantis.species.code=='IsPredatoryFish') # row number in basic that corresponds to IsPredatoryFish 
  sp=which(basic[guild,]==1)                                  # functional groups (columns) that meet IsPredatoryFish
  sp=sp-1                                                     # correct functional group identifer 
  biom_var1=rowSums(b.data.for.calc[,sp])                     # total biomass
  # get fish B
  guild=which(basic$Atlantis.species.code=='IsFish')     # row number in basic that corresponds to IsFish 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
  sp=sp-1                                                # correct functional group identifer 
  biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
  # save data into data frame 
  data.for.plot.PredFishrop$data <- biom_var1 / biom_var2
  names(data.for.plot.PredFishrop)[i+2] <- tier[i]
  # purge old
  remove(guild,sp,biom_var1,biom_var2)
  # 
  #
  # Get Time Series data for Dem.pel.fish
  #
  #
  # get DemT B
  guild=which(basic$Atlantis.species.code=='IsDemersalTeleost') # row number in basic that corresponds to IsDemersalTeleost 
  sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
  sp=sp-1                                                       # correct functional group identifer 
  biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
  # get PelT B
  guild=which(basic$Atlantis.species.code=='IsPelagicTeleost') # row number in basic that corresponds to IsPelagicTeleost 
  sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
  sp=sp-1                                                      # correct functional group identifer 
  biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
  # save data into data frame 
  data.for.plot.DemPelFish$data <- biom_var1 / biom_var2
  names(data.for.plot.DemPelFish)[i+2] <- tier[i]
  # purge old
  remove(guild,sp,biom_var1,biom_var2)
  # 
  #
  # Get Time Series data for Dem.Pel
  #
  #
  # get Dem B
  guild=which(basic$Atlantis.species.code=='IsDemersal') # row number in basic that corresponds to IsDemersalTeleost 
  sp=which(basic[guild,]==1)                                    # functional groups (columns) that meet IsDemersalTeleost
  sp=sp-1                                                       # correct functional group identifer 
  biom_var1=rowSums(b.data.for.calc[,sp])                       # total biomass
  # get Pel B
  guild=which(basic$Atlantis.species.code=='IsPelagic') # row number in basic that corresponds to IsPelagicTeleost 
  sp=which(basic[guild,]==1)                                   # functional groups (columns) that meet IsPelagicTeleost
  sp=sp-1                                                      # correct functional group identifer 
  biom_var2=rowSums(b.data.for.calc[,sp])                      # total biomass
  # save data into data frame 
  data.for.plot.DemPel$data <- biom_var1 / biom_var2
  names(data.for.plot.DemPel)[i+2] <- tier[i]
  # purge old
  remove(guild,sp,biom_var1,biom_var2)
  # 
  #
  # Get Time Series data for Dem.bio.PP
  #
  #
  # get total demersal B   IsDemersal
  guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(b.data.for.calc[,sp])                 # total biomass
  # get primary production PP
  guild=which(basic$Atlantis.species.code=='IsPrimaryProducer')  # row number in basic that corresponds to IsPrimaryProducer 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPrimaryProducer
  sp=sp-1                                                # correct functional group identifer 
  # the following function get Prodn for each fid, sums it across fids, and converts to the appropriate units (B) for computing the ratio
  PPout <- get_PrimaryProduction(sp, # NOTE !!!! this function was developed under the assumption that FunctionalGroupInfo is in the Environment
                                 paste(tier[i],"/ina_results_001PROD.nc",sep=""), 
                                 File_bgm)
  # temporal subset of PPout
  # !!!! you need to determine the save locations corresponding to the appropriate annual time steps of concern
  tempPP <- PPout[(biom$Time) %in% (data.grab)]
  # save data into data frame 
  data.for.plot.DemBioPP$data <- biom_var1 / tempPP
  names(data.for.plot.DemBioPP)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,tempPP,PPout,sp,guild)
  # 
  #
  # Get Time Series data for Pel.C
  #
  #
  # get total pelagic B
  guild=which(basic$Atlantis.species.code=='IsPelagic')  # row number in basic that corresponds to IsPelagic 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsPelagic
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
  # save data into data frame 
  data.for.plot.PelC$data <- biom_var1
  names(data.for.plot.PelC)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,sp,guild)
  # 
  #
  # Get Time Series data for Tot.C
  #
  #
  # get total B
  biom_var1=rowSums(c.data.for.calc[,c(2:dim(basic)[2])]) # total biomass
  # save data into data frame 
  data.for.plot.TotC$data <- biom_var1
  names(data.for.plot.TotC)[i+2] <- tier[i]
  # purge old
  remove(biom_var1)
  # 
  #
  # Get Time Series data for MTL.C
  #
  #
  # create place holder
  c.data.for.calc$MTL <- NA
  # need to loop through each yr and fid
  for (tempyr in c(1:nrow(c.data.for.calc))) {
    TLs <- as.numeric(as.character(basic[which(basic$Atlantis.species.code=='TrophicLevel'),-1]))
    biom.var <- c.data.for.calc[tempyr,c(2:(ncol(basic)-1))]
    numtr <- biom.var * TLs
    c.data.for.calc[tempyr,"MTL"] <- sum(numtr) / sum(biom.var)
  }
  # save data into data frame 
  data.for.plot.MTLC$data <- c.data.for.calc$MTL
  names(data.for.plot.MTLC)[i+2] <- tier[i]
  # purge old
  remove(TLs, biom.var, numtr)
  # 
  #
  # Get Time Series data for Fish.exp.rt
  #
  #
  # get targeted fish groups
  guild1=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
  guild2=which(basic$Atlantis.species.code=='IsFish')   # row number in basic that corresponds to IsTarget 
  sp1=which(basic[guild1,]==1)                          # functional groups (columns) 
  sp2=which(basic[guild2,]==1)                          # functional groups (columns) 
  sp <- sp1[sp1 %in% sp2]                               # targeted groups that are also fish
  sp=sp-1                                               # correct functional group identifer 
  biom_var1=rowSums(c.data.for.calc[,sp])               # total catch
  biom_var2=rowSums(b.data.for.calc[,sp])               # total biomass
  # save data into data frame 
  data.for.plot.FishExpRT$data <- biom_var1 / biom_var2
  names(data.for.plot.FishExpRT)[i+2] <- tier[i]
  # purge old
  remove(guild1, guild2, sp1, sp2, biom_var1,  biom_var2)
  # 
  #
  # Get Time Series data for Exp.rt
  #
  #
  # get targeted groups
  guild=which(basic$Atlantis.species.code=='IsTarget') # row number in basic that corresponds to IsTarget 
  sp=which(basic[guild,]==1)                             # functional groups (columns) 
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(c.data.for.calc[,sp])                 # total catch
  biom_var2=rowSums(b.data.for.calc[,sp])                 # total biomass
  # save data into data frame 
  data.for.plot.ExpRT$data <- biom_var1 / biom_var2
  names(data.for.plot.ExpRT)[i+2] <- tier[i]
  # purge old
  remove(guild, sp, biom_var1, biom_var2)
  # 
  #
  # Get Time Series data for Val
  #
  #
  # get total fish C
  guild=which(basic$Atlantis.species.code=='USDollarsPerTon') # row number in basic that corresponds to USDollarsPerTon 
  sp=which(is.na(basic[guild,]) == F)                         # functional groups (columns) that that have USDollarsPerTon
  biom_prc <- basic[guild,sp[-1]]                            # prices
  sp=sp[-1]-1                                                 # correct functional group identifer 
  biom_totc=c.data.for.calc[,sp]                              # catch (ton) of fids 
  valprod <- t(t(as.matrix(biom_totc)) * as.numeric(biom_prc))
  biom_var1=rowSums(valprod)                    
  # save data into data frame 
  data.for.plot.Val$data <- biom_var1
  names(data.for.plot.Val)[i+2] <- tier[i]
  # purge old
  remove(guild, sp, biom_prc, biom_totc, valprod, biom_var1)
  # 
  #
  # Get Time Series data for Fish.C
  #
  #
  # get total fish C
  guild=which(basic$Atlantis.species.code=='IsFish')  # row number in basic that corresponds to IsFish 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsFish
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
  # save data into data frame 
  data.for.plot.FishC$data <- biom_var1
  names(data.for.plot.FishC)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,sp,guild)
  # 
  #
  # Get Time Series data for Dem.C
  #
  #
  # get total demersal C
  guild=which(basic$Atlantis.species.code=='IsDemersal')  # row number in basic that corresponds to IsDemersal 
  sp=which(basic[guild,]==1)                             # functional groups (columns) that meet IsDemersal
  sp=sp-1                                                # correct functional group identifer 
  biom_var1=rowSums(c.data.for.calc[,sp])                 # total biomass
  # save data into data frame 
  data.for.plot.DemC$data <- biom_var1
  names(data.for.plot.DemC)[i+2] <- tier[i]
  # purge old
  remove(biom_var1,sp,guild)
  #
  #
  # Finished
  #
  #
}

#
#
# MAKE PLOTS
#
#
  
# # convert data from wide to long
library(tidyr)

#
#
# Get Data For Plot
#
#

temp1 <- data.for.plot.PelBPP; temp1$temp <- NULL
temp1[,-1] <- temp1[,-1] / temp1$output_17_26.12.19
temp1_long <- gather(temp1, model, PelBioPP, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp2 <- data.for.plot.BioPP; temp2$temp <- NULL
temp2[,-1] <- temp2[,-1] / temp2$output_17_26.12.19
temp2_long <- gather(temp2, model, BioPP, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp3 <- data.for.plot.MTLbio; temp3$temp <- NULL
temp3[,-1] <- temp3[,-1] / temp3$output_17_26.12.19
temp3_long <- gather(temp3, model, MTLBio, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp4 <- data.for.plot.PredFishrop; temp4$temp <- NULL
temp4[,-1] <- temp4[,-1] / temp4$output_17_26.12.19
temp4_long <- gather(temp4, model, PredFishProp, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp5 <- data.for.plot.DemPelFish; temp5$temp <- NULL
temp5[,-1] <- temp5[,-1] / temp5$output_17_26.12.19
temp5_long <- gather(temp5, model, DemPelFish, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp6 <- data.for.plot.DemPel; temp6$temp <- NULL
temp6[,-1] <- temp6[,-1] / temp6$output_17_26.12.19
temp6_long <- gather(temp6, model, DemPel, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp7 <- data.for.plot.DemBioPP; temp7$temp <- NULL
temp7[,-1] <- temp7[,-1] / temp7$output_17_26.12.19
temp7_long <- gather(temp7, model, DemBioPP, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp8 <- data.for.plot.PelC; temp8$temp <- NULL
temp8[,-1] <- temp8[,-1] / temp8$output_17_26.12.19
temp8_long <- gather(temp8, model, PelC, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp9 <- data.for.plot.TotC; temp9$temp <- NULL
temp9[,-1] <- temp9[,-1] / temp9$output_17_26.12.19
temp9_long <- gather(temp9, model, TotC, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp10 <- data.for.plot.MTLC; temp10$temp <- NULL
temp10[,-1] <- temp10[,-1] / temp10$output_17_26.12.19
temp10_long <- gather(temp10, model, MTLC, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp11 <- data.for.plot.FishExpRT; temp11$temp <- NULL
temp11[,-1] <- temp11[,-1] / temp11$output_17_26.12.19
temp11_long <- gather(temp11, model, FishExpRT, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp12 <- data.for.plot.ExpRT; temp12$temp <- NULL
temp12[,-1] <- temp12[,-1] / temp12$output_17_26.12.19
temp12_long <- gather(temp12, model, ExpRT, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp13 <- data.for.plot.Val; temp13$temp <- NULL
temp13[,-1] <- temp13[,-1] / temp13$output_17_26.12.19
temp13_long <- gather(temp13, model, Val, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp14 <- data.for.plot.FishC; temp14$temp <- NULL
temp14[,-1] <- temp14[,-1] / temp14$output_17_26.12.19
temp14_long <- gather(temp14, model, FishC, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)

temp15 <- data.for.plot.DemC; temp15$temp <- NULL
temp15[,-1] <- temp15[,-1] / temp15$output_17_26.12.19
temp15_long <- gather(temp15, model, DemC, output_01_26.01.19:output_17_26.12.19, factor_key=TRUE)


# The data needs to be split between a normal scale and log scal for the plots

# First make a plot for the results on a normal scale

# Merge data for plot:
tempmerge <- Reduce(function(x, y) merge(x, y, by=c("Time","model"), all=TRUE), list(temp3_long,
                                                                                     temp4_long, temp5_long, temp6_long,
                                                                                     temp8_long, temp9_long,
                                                                                     temp10_long, temp11_long, temp12_long,
                                                                                     temp13_long, temp14_long, temp15_long))
# convert data for plot
tempmerge_long <- gather(tempmerge, indicator, value, MTLBio:DemC, factor_key=TRUE)

# plot
png("TimeSeriesPlot_indicators.png", width = 12, height = 8, units = "in", res = 300)
ggplot(tempmerge_long, aes(x = Time, y = value, color = model, linetype = model)) +
  geom_line() +
  theme_bw() +
  scale_linetype_manual("", 
                        values = c(rep("solid",12), "dashed"),
                        labels = tier.t) +
  scale_color_manual("", 
                     values = c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey"),
                     labels = tier.t) +
  facet_wrap(~ indicator, scales = "free", ncol = 4) +
  #theme(axis.text.x = element_text(size = 8)) +
  ylab("Value (relative to the baseline)") +
  xlab("")
dev.off()

# now make plots for data on log-scale

# merge data
tempmergelog <- Reduce(function(x, y) merge(x, y, by=c("Time","model"), all=TRUE), list(temp1_long, 
                                                                                     temp2_long, 
                                                                                     temp7_long))

# convert data
tempmergelog_long <- gather(tempmergelog, indicator, value, PelBioPP:DemBioPP, factor_key=TRUE)

# make plot
png("TimeSeriesPlot_indicatorslog.png", width = 11.25, height = 2.5, units = "in", res = 300)
ggplot(tempmergelog_long, aes(x = Time, y = value, color = model, linetype = model)) +
  geom_line() +
  theme_bw() +
  scale_linetype_manual("", 
                        values = c(rep("solid",12), "dashed"),
                        labels = tier.t,
                        guide = guide_legend(override.aes = list(color = "white"))) + # make legend go away while keeping plot dimensions
  scale_color_manual("", 
                     values = c(cCO10,cCO50,cCO90,cHE10,cHE50,cHE90,cHA10,cHA50,cHA90,cMX10,cMX50,cMX90,"grey"),
                     labels = tier.t,
                     guide = "none") +
  theme(legend.text = element_text(color = "white"), # make legend go away while keeping plot dimensions
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white")) +
  scale_y_continuous(breaks=c(0.0004,0.04,1,400,40000), trans='log10') +
  facet_wrap(~ indicator, scales = "free") +
  ylab("(log scale)") +
  xlab("Time (post MME)") 
dev.off()
