#!/bin/bash

wdir=./

file_ext=Cumberland_POP_CNP



cp $wdir/restart_files/${file_ext}_pop_rst_transient.nc  ./restart_files/${file_ext}_pop_rst.nc 
cp  $wdir/restart_files/${file_ext}_climate_rst_transient.nc ./restart_files/${file_ext}_climate_rst.nc 
cp  $wdir/restart_files/${file_ext}_casa_rst_transient.nc ./restart_files/${file_ext}_casa_rst.nc 
cp  $wdir/restart_files/${file_ext}_cable_rst_transient.nc ./restart_files/${file_ext}_cable_rst.nc.
