#!/usr/bin/env Rscript

library(ncdf4)

model     <- commandArgs(trailingOnly = TRUE)[1]
exp       <- commandArgs(trailingOnly = TRUE)[2]
startyear <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
endyear   <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
inpath    <- commandArgs(trailingOnly = TRUE)[5]
outpath   <- commandArgs(trailingOnly = TRUE)[6]

years <- c(startyear:endyear)

## 1) load files 
global_file  <- nc_open(paste0(inpath,"/nbp_global.nc"))
north_file   <- nc_open(paste0(inpath,"/nbp_north.nc"))
tropics_file <- nc_open(paste0(inpath,"/nbp_tropics.nc"))
south_file   <- nc_open(paste0(inpath,"/nbp_south.nc"))

global  <- ncvar_get(global_file,"nbp")
north   <- ncvar_get(north_file,"nbp")
tropics <- ncvar_get(tropics_file,"nbp")
south   <- ncvar_get(south_file,"nbp")


## 2) calculate annual NBP per zone per year
nbp_global  <- apply(global,3,sum,na.rm=T)
nbp_north   <- apply(north,3,sum,na.rm=T)
nbp_tropics <- apply(tropics,3,sum,na.rm=T)
nbp_south   <- apply(south,3,sum,na.rm=T)

## 3) write ascii txt files
outdat <- data.frame(years,nbp_global,nbp_north,nbp_tropics,nbp_south)
colnames(outdat) <- c("Year","Global","North","Tropics","South")   # see protocol

write.table(outdat,file=paste0(outpath,"/",model,"_zonalNBP_",exp,".txt"),quote=FALSE,row.names=FALSE)
