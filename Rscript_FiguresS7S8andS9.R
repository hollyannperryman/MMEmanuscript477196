
#
# The following script generates Figures S7-9 
#
# Olsen, E., Eide, C.H., Nilsen, I., Perryman, H.A. and 
# Vikebø, F., 2019. Ecological effects and ecosystem 
# shifts caused by mass mortality events on early life stages
# of fish. Frontiers in Marine Science, 6, p.669.
#

# This script was developed by Holly Perryman based on 
# discussions with Cecilie Hansen, Ina Nelson and Erik Olsen

# initialize misc. parmeters

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

# define locations with simulation outputs 

outputfolders <- list.files()[which(substr(list.files(),1,6) == "output")]


# get data
#  this script focuses on the extreme, single spp runs

Temp17 <- read.csv(paste(outputfolders[14],"/ina_results_001DietCheck.txt",sep=""),
                   sep="",
                   stringsAsFactors=FALSE)
Temp03 <- read.csv(paste(outputfolders[4],"/ina_results_001DietCheck.txt",sep=""),
                   sep="",
                   stringsAsFactors=FALSE)
Temp06 <- read.csv(paste(outputfolders[7],"/ina_results_001DietCheck.txt",sep=""),
                   sep="",
                   stringsAsFactors=FALSE)
Temp09 <- read.csv(paste(outputfolders[10],"/ina_results_001DietCheck.txt",sep=""),
                   sep="",
                   stringsAsFactors=FALSE)

# define time steps of focus
#  in this script you are grabbing all the ts that relate to the years of focus
#  so you can sum across the whole year 
#tofcon <- c(c((eventyr - 1):(eventyr + 4))*365,(eventyr + 10)*365,(eventyr + 25)*365,(eventyr + 45)*365)
tofcon <- round(c(seq(eventyr - 1.8, eventyr - 1, 0.2)*365,
            seq(eventyr - 0.8, eventyr - 0, 0.2)*365,
            seq(eventyr + 0.2, eventyr + 1, 0.2)*365,
            seq(eventyr + 1.2, eventyr + 2, 0.2)*365,
            #seq(eventyr + 2.2, eventyr + 3, 0.2)*365,
            #seq(eventyr + 3.2, eventyr + 4, 0.2)*365,
            #seq(eventyr + 9.2, eventyr + 10, 0.2)*365,
            seq(eventyr + 24.2, eventyr + 25, 0.2)*365))
            #seq(eventyr + 44.2, eventyr + 45, 0.2)*365))
# # check
# tofcon %in% unique(Temp17$Time) 

# subset data
#  remove t that is not of interest: 
Temp17sub <- Temp17[which(Temp17$Time %in% tofcon),] 
Temp03sub <- Temp03[which(Temp03$Time %in% tofcon),] 
Temp06sub <- Temp06[which(Temp06$Time %in% tofcon),] 
Temp09sub <- Temp09[which(Temp09$Time %in% tofcon),] 


#
#
# First plot 
# radar plot - change in predation on prey groups 
#
#


# select prey group 
#  this can be set to 1, 2, or 3
grpg <- 3

# get predators
# keep in mind predators can change through the diff scenarios
temp1<-Temp17sub[which(Temp17sub[,grpshnms[grpg]] > 0),]
temp1<-aggregate(temp1[,grpshnms[grpg]],by=list(Predator = temp1$Predator),FUN=sum)
#temp1$stuff <- round(temp1$x / sum(temp1$x),3);temp1
temp.pred17 <- temp1[which(temp1$x / sum(temp1$x) > 0.05),"Predator"]
remove(temp1)
if (grpg == 1) { #haddock
  temp1<-Temp09sub[which(Temp09sub[,grpshnms[grpg]] > 0),]
  temp1<-aggregate(temp1[,grpshnms[grpg]],by=list(Predator = temp1$Predator),FUN=sum)
  temp.scenario <- temp1[which(temp1$x / sum(temp1$x) > 0.05),"Predator"]
} else if (grpg == 2) { #herring
  temp1<-Temp06sub[which(Temp06sub[,grpshnms[grpg]] > 0),]
  temp1<-aggregate(temp1[,grpshnms[grpg]],by=list(Predator = temp1$Predator),FUN=sum)
  temp.scenario <- temp1[which(temp1$x / sum(temp1$x) > 0.05),"Predator"]
} else if (grpg == 3) { #cod
  temp1<-Temp03sub[which(Temp03sub[,grpshnms[grpg]] > 0),]
  temp1<-aggregate(temp1[,grpshnms[grpg]],by=list(Predator = temp1$Predator),FUN=sum)
  temp.scenario <- temp1[which(temp1$x / sum(temp1$x) > 0.05),"Predator"]
}
remove(temp1)
#
my.predators <- unique(c(temp.pred17,temp.scenario))
remove(temp.pred17,temp.scenario)
# we decided to specify the key predators of focus for cod (grpg = 3) but the others we can search for
if(grpg == 3){my.predators <- c("DEL","DEO","HAD","GRH","PEL")}
# check
my.predators

# subset data - BASE
# subset by my.predators
# and prep for creating data for plot
Temp17sub2 <- Temp17sub[which(Temp17sub$Predator %in% my.predators),] # subset by predators
Temp17sub2 <- Temp17sub2[,which(names(Temp17sub2) %in% c(names(Temp17sub2)[1:5],grpshnms[grpg]))] # drop prey items that are not of concern
Temp17sub2$Time <- round((Temp17sub2$Time / 365)-(eventyr-1),1) # change scale of time
Temp17sub2$Year <- ceiling(as.numeric(Temp17sub2$Time))-1 # create new temporal variable indicating the year each timestep corresponds to
Temp17.predators <- aggregate(Temp17sub2[,grpshnms[grpg]], by=list(Year=Temp17sub2$Year, 
                                                                   Predator = as.character(Temp17sub2$Predator)), 
                              FUN=sum) # aggregate by Year and sum
Temp17.predators <- spread(Temp17.predators, Predator, x) # convert from long to wide 

# subset data - extreme scenarios
# subset by my.predators
# and prep for creating data for plot
if(grpg==3){ # cod
  # subset data - CO90
  Temp03sub2 <- Temp03sub[which(Temp03sub$Predator %in% my.predators),]
  Temp03sub2 <- Temp03sub2[,which(names(Temp03sub2) %in% c(names(Temp03sub2)[1:5],grpshnms[grpg]))]
  
  Temp03sub2$Time <- round((Temp03sub2$Time / 365)-(eventyr-1),1)
  Temp03sub2$Year <- ceiling(as.numeric(Temp03sub2$Time))-1
  Temp03.predators <- aggregate(Temp03sub2[,grpshnms[grpg]], by=list(Year=Temp03sub2$Year, 
                                                                     Predator = as.character(Temp03sub2$Predator)), 
                                FUN=sum)
  Temp03.predators <- spread(Temp03.predators, Predator, x)
  
  runcompare <- Temp03.predators
  
} else if (grpg==1){ #haddock
  # subset data - HA90
  Temp09sub2 <- Temp09sub[which(Temp09sub$Predator %in% my.predators),]
  Temp09sub2 <- Temp09sub2[,which(names(Temp09sub2) %in% c(names(Temp09sub2)[1:5],grpshnms[grpg]))]
  Temp09sub2$Time <- round((Temp09sub2$Time / 365)-(eventyr-1),1) 
  Temp09sub2$Year <- ceiling(as.numeric(Temp09sub2$Time))-1
  Temp09.predators <- aggregate(Temp09sub2[,grpshnms[grpg]], by=list(Year=Temp09sub2$Year, 
                                                                     Predator = as.character(Temp09sub2$Predator)), 
                                FUN=sum)
  Temp09.predators <- spread(Temp09.predators, Predator, x)
  
  runcompare <- Temp09.predators

}else if (grpg==2){ #herring
  # subset data - HE90
  Temp06sub2 <- Temp06sub[which(Temp06sub$Predator %in% my.predators),]
  Temp06sub2 <- Temp06sub2[,which(names(Temp06sub2) %in% c(names(Temp06sub2)[1:5],grpshnms[grpg]))] 
  Temp06sub2$Time <- round((Temp06sub2$Time / 365)-(eventyr-1),1) 
  Temp06sub2$Year <- ceiling(as.numeric(Temp06sub2$Time))-1
  Temp06.predators <- aggregate(Temp06sub2[,grpshnms[grpg]], by=list(Year=Temp06sub2$Year, 
                                                                     Predator = as.character(Temp06sub2$Predator)), 
                                FUN=sum)
  Temp06.predators <- spread(Temp06.predators, Predator, x)
  
  runcompare <- Temp06.predators
}

#
# create data for plot
#

dataforplot <- rbind(rep(1,length(my.predators)),
                     rep(-1, length(my.predators)),
                     (runcompare[,my.predators]/Temp17.predators[,my.predators])-1)
#
# make radar plot
#
library(fmsb)
#https://stackoverflow.com/questions/42029645/how-to-change-the-default-label-font-type-to-italic-in-radarchart-fmsb-packag
#https://rdrr.io/cran/fmsb/man/radarchart.html

# make texts for plot
if (grpg==1) {
  plotname <- "NE Arctic Haddock"
  txtforplot <- c("Large\ndemersals", "Other\ndemersals", "NE Arctic\ncod")
  modelrun <- "[HA90]"
} else if (grpg==2) {
  plotname <- "Spring Spawning Herring"
  txtforplot <- c("Killer\nwhales", "Other\nredfish", "Boreal\nseabirds")
  modelrun <- "[HE90]"
} else if (grpg==3) {
  plotname <- "NE Arctic Cod"
  txtforplot <- c("Large\ndemersals", "Other\ndemersals", "Haddock", "Greenland\nhalibut", "Large\npelagics")
  modelrun <- "[CO90]"
}

# make plot
png(file=paste("CHGpredON_",grpshnms[grpg],"90.png",sep = ""),
    width = 7, height = 6, units = "in",
    res = 300)
radarchart(dataforplot,
           axistype=1,caxislabels=seq(-1,1,0.5),axislabcol="grey40", calcex = 1,vlcex=1, 
           cglcol="grey",cglty=1,cglwd=0.8,
           vlabels = txtforplot, 
           pcol = c("black", "#F9B282", "#ED7C97", "#BC5AA9", "#419F44"), 
           plty = c(2,1,1,1,3),pty = rep(32,5), plwd = rep(2,5),
           #pcol = c("black", "#F3E79A", "#F9B282", "#ED7C97", "#BC5AA9", "#704D9E","#B4DCAB","#419F44","#004616"),
           title = paste("Change in Predation upon ",plotname,"\n",modelrun,sep="")) 
legend(1,1.3,
       title = "Time since MME",
       legend=c("1","2","3","25"),
       col=c("#F9B282", "#ED7C97", "#BC5AA9", "#419F44"),
       lty=c(1,1,1,3),
       bty = "n")
dev.off()





#
#
# Second plot 
# radar plot - change in prey consuption amongst key predators
#  these plots are to be plotted for each key predator of each focal prey group  
#
#




# select predator from list of my.predators
# for prey group 1 there are 3 main predators
# for prey group 2 there are 3 main predators 
# for prey group 3 there are 5 main predators 
predindex <- 3
remove(plotnotes)

# prep BASE data for plot
Temp17.predsub <- Temp17sub[which(Temp17sub$Predator == my.predators[predindex]),] # subset by pred
Temp17.predsub <- Temp17.predsub[, colSums(Temp17.predsub != 0) > 0] # remove groups that are not prey items
library(dplyr)
Temp17.predprey <- Temp17.predsub[!(names(Temp17.predsub) %in% c("Predator","Cohort","Updated"))] %>% group_by(Time) %>% summarise_each(list(sum)) # aggregate by cohort and updated 
Temp17.predprey$Time <- round((Temp17.predprey$Time / 365)-(eventyr-1),1) # change scale of time
Temp17.predprey$Year <- ceiling(as.numeric(Temp17.predprey$Time))-1 # make new temporal variable Year to id the year each ts corresponds to
Temp17.prey <- aggregate(Temp17.predprey[,!(names(Temp17.predprey) %in% c("Time","Year"))], by=list(Year=Temp17.predprey$Year), 
                         FUN=sum) # aggregate by Year

# prep scenario data for plot
if (grpg == 1) {
  Temp09.predsub <- Temp09sub[which(Temp09sub$Predator == my.predators[predindex]),] 
  Temp09.predsub <- Temp09.predsub[, colSums(Temp09.predsub != 0) > 0] 
  Temp09.predprey <- Temp09.predsub[!(names(Temp09.predsub) %in% c("Predator","Cohort","Updated"))] %>% group_by(Time) %>% summarise_each(list(sum))
  Temp09.predprey$Time <- round((Temp09.predprey$Time / 365)-(eventyr-1),1) 
  Temp09.predprey$Year <- ceiling(as.numeric(Temp09.predprey$Time))-1
  Temp09.prey <- aggregate(Temp09.predprey[,!(names(Temp09.predprey) %in% c("Time","Year"))], by=list(Year=Temp09.predprey$Year), 
                           FUN=sum)
  scenario.data <- Temp09.prey
} else if (grpg == 2) {
  Temp06.predsub <- Temp06sub[which(Temp06sub$Predator == my.predators[predindex]),] 
  Temp06.predsub <- Temp06.predsub[, colSums(Temp06.predsub != 0) > 0] 
  Temp06.predprey <- Temp06.predsub[!(names(Temp06.predsub) %in% c("Predator","Cohort","Updated"))] %>% group_by(Time) %>% summarise_each(list(sum))
  Temp06.predprey$Time <- round((Temp06.predprey$Time / 365)-(eventyr-1),1) 
  Temp06.predprey$Year <- ceiling(as.numeric(Temp06.predprey$Time))-1
  Temp06.prey <- aggregate(Temp06.predprey[,!(names(Temp06.predprey) %in% c("Time","Year"))], by=list(Year=Temp06.predprey$Year), 
                           FUN=sum)
  scenario.data <- Temp06.prey
} else if(grpg == 3){
  Temp03.predsub <- Temp03sub[which(Temp03sub$Predator == my.predators[predindex]),] 
  Temp03.predsub <- Temp03.predsub[, colSums(Temp03.predsub != 0) > 0] 
  Temp03.predprey <- Temp03.predsub[!(names(Temp03.predsub) %in% c("Predator","Cohort","Updated"))] %>% group_by(Time) %>% summarise_each(list(sum))
  Temp03.predprey$Time <- round((Temp03.predprey$Time / 365)-(eventyr-1),1) 
  Temp03.predprey$Year <- ceiling(as.numeric(Temp03.predprey$Time))-1
  Temp03.prey <- aggregate(Temp03.predprey[,!(names(Temp03.predprey) %in% c("Time","Year"))], by=list(Year=Temp03.predprey$Year), 
                           FUN=sum)
  scenario.data <- Temp03.prey
}

# select all prey of predator
# my.prey <- unique(c(names(Temp17.prey)[(!names(Temp17.prey) %in% c("Time","Year"))],
#                     names(scenario.data)[(!names(scenario.data) %in% c("Time","Year"))]))
#
# the question arose - how important is the focal prey to the predator?
# base
temp1 <- Temp17.predsub[!(names(Temp17.predsub) %in% c("Time","Cohort","Updated"))] %>% group_by(Predator) %>% summarise_each(list(sum)) # aggregate 
temp1 <- temp1[,which(temp1[,-1]/sum(temp1[,-1]) > 0.05)+1] # select prey, need to add 1 to account for skipping over predator column
keyprey1 <- names(temp1)
remove(temp1)
# senarios
if (grpg == 1){
  temp1 <- Temp09.predsub[!(names(Temp09.predsub) %in% c("Time","Cohort","Updated"))] %>% group_by(Predator) %>% summarise_each(list(sum)) # aggregate 
  temp1 <- temp1[,which(temp1[,-1]/sum(temp1[,-1]) > 0.05)+1]
  keyprey2 <- names(temp1)
  remove(temp1)
} else if (grpg == 2){
  temp1 <- Temp06.predsub[!(names(Temp06.predsub) %in% c("Time","Cohort","Updated"))] %>% group_by(Predator) %>% summarise_each(list(sum)) # aggregate 
  temp1 <- temp1[,which(temp1[,-1]/sum(temp1[,-1]) > 0.05)+1]
  keyprey2 <- names(temp1)
  remove(temp1)
} else if (grpg == 3){
  temp1 <- Temp03.predsub[!(names(Temp03.predsub) %in% c("Time","Cohort","Updated"))] %>% group_by(Predator) %>% summarise_each(list(sum)) # aggregate 
  temp1 <- temp1[,which(temp1[,-1]/sum(temp1[,-1]) > 0.05)+1]
  keyprey2 <- names(temp1)
  remove(temp1)
}
# get key prey
my.prey <- unique(c(keyprey1,keyprey2))
remove(keyprey1,keyprey2)
# is focal prey a key prey, if not add it for the sake of the graph
if (!(grpshnms[grpg] %in% my.prey)){
  my.prey <- c(my.prey,grpshnms[grpg])
  plotnotes <- "focal prey not key prey group"
}
plotnotes
my.prey

# We had to specify prey for KWH as not enough were selected based on the criteria
if (grpg == 2 & predindex == 1){my.prey <- c("SSH","CEP","NCO","MAC")}

# make data for plot
dataforplot <- rbind(rep(1,length(my.prey)),
                     rep(-1, length(my.prey)),
                     (scenario.data[,my.prey]/Temp17.prey[,my.prey])-1)

# make texts for plot
if (grpg==1) {
  modelrun <- "[HA90]"
  if(predindex == 1){
    plotname <- "Large Demersal";
    txtforplot <- c("Squid","Large\nzooplankton","Medium zooplankton","Haddock*")
  } else if(predindex == 2){
    plotname <- "Other Demersal";
    txtforplot <- c("Prawns","Large\nzooplankton","Medium zooplankton","Detritus","Haddock*")
  } else if(predindex == 3){
    plotname <- "NE Arctic Cod";
    txtforplot <- c("Capelin","Prawns","Squid","Large\nzooplankton","Haddock*")
  }
} else if (grpg==2) {
  modelrun <- "[HE90]"
  if(predindex == 1){
    plotname <- "Killer Whale";
    txtforplot <- c("S. S. herring","Squid","NE Arctic cod","Mackerel")
  } else if(predindex == 2){
    plotname <- "Other Redfish";
    txtforplot <- c("Blue whiting","S. S.\nherring","Large zooplankton","Medium\nzooplankton")
  } else if(predindex == 3){
    plotname <- "Boreal Seabird";
    txtforplot <- c("S. S. herring","Prawns","Large zooplankton","Medium\nzooplankton")
  }
} else if (grpg==3) {
  modelrun <- "[CO90]"
  if(predindex == 1){
    plotname <- "Large Demersal"
    txtforplot <- c("Mesopelagics","Large\nzooplankton","Medium zooplankton","NE Arctic\ncod*")
  } else if(predindex == 2){
    plotname <- "Other Demersal"
    txtforplot <- c("Prawns","Large\nzooplankton","Medium\nzooplankton","Detritus","NE Arctic\ncod*")
  } else if(predindex == 3){
    plotname <- "NE Arctic Haddock"
    txtforplot <- c("Prawns","Large\nzooplankton","Medium zooplankton","NE Arctic\ncod*")
  } else if(predindex == 4){
    plotname <- "Greenland Halibut"
    txtforplot <- c("Large zooplankton","Medium\nzooplankton","NE Arctic\ncod*")
  } else if(predindex == 5){
    plotname <- "Large Pelagic"
    txtforplot <- c("Mesopelagics","Squid","Large\nzooplankton","Medium\nzooplankton","NE Arctic\ncod*")
  }
}

# make plot

png(file=paste("CHGpredOF_",my.predators[predindex],"_",grpshnms[grpg],"90.png",sep = ""),
    width = 7, height = 5, units = "in",
    res = 300)
radarchart(dataforplot,
           axistype=1,caxislabels=seq(-1,1,0.5),axislabcol="grey40", calcex = 1,vlcex=1, 
           cglcol="grey",cglty=1,cglwd=0.8,
           vlabels = txtforplot, 
           pcol = c("black", "#F9B282", "#ED7C97", "#BC5AA9", "#419F44"), 
           plty = c(2,1,1,1,3),pty = rep(32,5), plwd = rep(2,5),
           title = paste("Change in ",plotname," Predation\n",modelrun,sep="")) 
# legend(1.2,1.3,
#        title = "Time since MME",
#        legend=c("0","1","2","25"),
#        col=c("#F9B282", "#ED7C97", "#BC5AA9", "#419F44"),
#        lty=c(1,1,1,3),
#        bty = "n")
dev.off()
plotnotes
