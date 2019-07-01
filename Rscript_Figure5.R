
#
# The following script generates Figure 5 for 
#
# Ecological effects and ecosystem shifts caused by mass mortality events on early life stages of fish
# Erik Olsen*, Cecilie H. Eide, Ina Nilsen, Holly A. Perryman, Frode Vikebø
# Submitted to Journal:Frontiers in Marine Science
# Specialty Section:Global Change and the Future Ocean
# Article type:Original Research Article
#

# This script was developed by Holly Perryman
# based on discussions with Cecilie Hansen, Ina Nelson and Erik Olsen
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
# Define the temporal span of data going into calculations
t.span.out <- 20
# seasonality flag (T/F):
# (T) - you are collecting all of the seasonal data reported in the out files 
# (F) - you are collecting only the data reported at the end of the year (i.e., every 365 time steps) 
t.season <- F
# !!! I DO NOT HAVE IT PROGRAMMED FOR t.season <- T YET !!! 
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
    # !!! NEED TO PROGRAM !!! 
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
