#!/usr/bin/env Rscript

## Splits up global land mask into n number of land masks to be used in serial runs
## Makes use of raster objects rather than gridlists to speed things up.
## Expects land mask as input with 1=land and 0 or NA = non-land

library("ncdf4")
library("raster")

landmask_file <- commandArgs(trailingOnly=TRUE)[1]
nmasks        <- as.numeric(commandArgs(trailingOnly=TRUE)[2])
outdir        <- commandArgs(trailingOnly=TRUE)[3] 
bbox          <- commandArgs(trailingOnly=TRUE)[4]

gl   <- raster(landmask_file)  # gl = global landmask
if (bbox != "global"){
  bounds    <- extent(as.numeric(unlist(strsplit(bbox,","))))  # spatial extent object
  lmask     <- crop(gl,bounds)
  gl[]      <- NA
  gl        <- merge(gl,lmask)
  names(gl) <- "land"
}

vals <- getValues(gl)

ntot     <- sum(vals,na.rm=T)    # total number of land grid cells
npermask <- ceiling(ntot/nmasks) # land grid cells per mask

vals[!is.na(vals)] <- vals[!is.na(vals)] + sort(rep(1:nmasks-1,npermask))[1:ntot]

for (i in 1:nmasks){
  mask <- vals
  mask[mask != i] <- NA_integer_
  mask[mask == i] <- 1
  mask <- setValues(gl,mask)
  
  outdir_full <- paste0(outdir,"/run",i,"/landmask")
  dir.create(outdir_full,recursive=TRUE,showWarnings=FALSE)
  writeRaster(mask,paste0(outdir_full,"/landmask",i,".nc"),format="CDF",datatype="INT1S",overwrite=TRUE,
              NAflag=0,varunit="0:no land, 1:land",longname="land")
}


