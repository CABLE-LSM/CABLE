#!/usr/bin/env Rscript

library(ncdf4)

exp       <- commandArgs(trailingOnly = TRUE)[1]
vars      <- commandArgs(trailingOnly = TRUE)[2]
startyear <- as.numeric(commandArgs(trailingOnly = TRUE)[3])
endyear   <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
inpath    <- commandArgs(trailingOnly = TRUE)[5]
outpath   <- commandArgs(trailingOnly = TRUE)[6]

sectoday <- 86400
kgtog    <- 1000

vars  <- unlist(strsplit(vars," "))
years <- c(startyear:endyear)

res_annual  <- data.frame(matrix(NA,ncol=0,nrow=length(years),dimnames=list(years,NULL)))
res_monthly <- data.frame(matrix(NA,ncol=0,nrow=length(years)*12,dimnames=list(paste(sort(rep(years,12)),c(1:12),sep="."),NULL)))


## 1) load files and get variables
for (var in vars){

  varfiles <- list.files(path=inpath,pattern=paste0(var,"_"))

  for (varfile in varfiles){

    # naming
    lon1 <- unlist(strsplit(unlist(strsplit(varfile,"_"))[2],"\\."))[1]
    lon2 <- unlist(strsplit(unlist(strsplit(varfile,"_"))[3],"\\."))[1]
    lat1 <- unlist(strsplit(unlist(strsplit(varfile,"_"))[4],"\\."))[1]
    lat2 <- unlist(strsplit(unlist(strsplit(varfile,"_"))[5],"\\."))[1]
    
    if (lon1 == "-180" & lon2 == "180"){   # latitudinal bands?
      varname  <- paste(var,lat1,lat2,sep="_")
      if ((lat1 == "-60" | lat1 == "-90") & lat2 == "90" ){  # global?
        varname <- paste(var,"global",sep="_")
      }
    } else { # regions
      varname <- paste(var,lon1,lon2,lat1,lat2,sep="_")
    }
    
    
    # get files and put them into the result matrices
    ncfile  <- nc_open(varfile)
    vals    <- ncvar_get(ncfile,var)



    # unit conversions if necessary
    if (var %in% c("gpp","nbp","ra","npp","rh","fLUC")){
      vals <- vals * kgtog * sectoday # --> gC m-2 d-1
    } else if (var %in% c("pr","mrro","evapotrans","tran","evspsblveg","evspsblsoi")){
      vals <- vals * sectoday   # --> mm d-1
    }

    if (grepl("annual",varfile)){
      if (var == vars[1] & length(vals) > length(years)){
        nmax <- length(years) 
      } else {
        nmax <- length(vals)
      }
      res_annual[,varname] <- vals[1:nmax]
    } else if (grepl("monthly",varfile)){
      res_monthly[,varname] <- vals[1:(nmax*12)]
    }

  }
  file.remove(paste0(inpath,"/",varfiles))   #remove varfiles
}


## 3) write ascii txt files
write.table(res_annual,file=paste0(outpath,"/",exp,"_annual.txt"),quote=FALSE,row.names=T)
write.table(res_monthly,file=paste0(outpath,"/",exp,"_monthly.txt"),quote=FALSE,row.names=T)
