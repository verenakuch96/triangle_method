# ----- Skript Triangle Method -----
# Script handels LST and NDVI images to calculte the traingle method.
# Author: verena Kuch following the script of Dr. Verena Huber Garcia (LMU)
#Last update: 2.8.2021
#Copyright (c) Verena Kuch, 2021
#Email: kuch.verena@gmail.com

#------------------------------------------------

#packages needed

rm(list=ls())
library(raster)
library(rgdal)
#install.packages('tidyverse')

library(tidyverse)
library(fields)

#install.packages ("rCurl")
library(RCurl)
#install.packages("devtools")
library(devtools)
#install_url("https://github.com/Terradue/rLandsat8/releases/download/v0.1-SNAPSHOT/rLandsat8_0.1.0.tar.gz")
library(rLandsat8)
library(tidyr)
library(dplyr)
#install.packages("scales")
#install.packages("rlang")




# Path for Data

LST_path='C:/Users/../LST'
NDVI_path='C:/Users/.../NDVI'



#Path for plots

plot_path='C:/Users/...'


# additional information

z=600 #height above sea level
LST_name='LST_date.tif'
NDVI_name='NDVI_date.tif'
ec_data='C:/Users/...'
mask1='C:/Users/mask'
image_date=as.Date(c("YYYY-MM-DD"))


#load the LST Data - important same size as NDVI is needed

setwd(LST_path)
LST = raster(LST_name)   
LST_C = LST - 273.15 #convert to °C

#load NDVI Data - important same size as LST is needed
#exclude all NDVI out of range (-1 to 1)

setwd(NDVI_path)
NDVI=raster(NDVI_name)
NDVI[][which(NDVI[]>=1)] = NA
NDVI[][which(NDVI[]<=0)] = NA

# Mask if NDVI and LST don't have the same size

#mask_schech=readOGR(mask1)
#LST=mask(x=LST,mask=mask_schech)
#NDVI=mask(x=NDVI, mask=mask_schech)
#K='klein'

# EC Data with temperature and net radiation

obs2 = read.csv(ec_data,header = TRUE, stringsAsFactors = FALSE,sep=';')  
obs2$date = as.Date(obs2$date,format="%d.%m.%Y")

# find line for overfly date and time
EC_overfly_hour_2 = obs2[which(obs2$date==image_date & obs2$time=='12:00:00'),] #here the overfly time is needed

# reclassify NDVI Data for the NDVI interval

a = seq((-1),0.95,by=0.05)
b = seq((-0.95),1,by=0.05)
c = seq((-0.95),1,by=0.05)
abc = cbind(a,b,c)
class.ndvi = reclassify(NDVI,abc)

#--------------calculations for wet and dry edge-----------

y1=sort(LST[], decreasing = T)[1]   # highest temperature, 
x1=0 # minimal NDVI
y2= min(LST[which(class.ndvi[]==sort(class.ndvi[], decreasing = T)[1])],na.rm = T )#min LST of highest NDVI
x2=sort(NDVI[], decreasing = T)[1]#highest NDVI
x2=1
m=(y2-y1)/(x2-x1)  # calculation of the slope
t=y1-m*x1 #calculation of intercept



#Triangle plot
par(mar=c(3,3,1,1),oma=c(1,1,1,1))	
smoothScatter(NDVI[], LST[]
              ,colramp = colorRampPalette(c("white", 'darkgreen'))
              , transformation = function(x) x^.19		# the smaller the value the more colourful
              , ylim=c(cellStats(LST,min)-1,cellStats(LST,max)+1)
              , xlim=c((-0.05),1.2)
              , xlab="", ylab="", bty="n"
              , main=paste("NDVI-LST-plot for scene",image_date, sep=" "))

mtext("NDVI",1,padj=4, cex=0.8)
mtext("Land Surface Temperature [K]",2,padj=-4, cex=0.8)
abline(v=0,col="green3",lwd=3) # Vertikal, durch NDVI=0
abline(h=paste(y2) ,col="blue",lwd=3, lty=2) # horizontal, wet edge
abline(a=t, b=m,col="red",lwd=3,lty=2) # observed dry edge




#---------Calculation of phi max, phi min and phi i following Stisen et al (2008)--------

LSTi.max = m*class.ndvi + t  # Calculate LSTi.max for every NDVI class - important no negative values in NDVI_class

#Calculation of phi
phi_max = 1.26
phi_min_raster = phi_max * ((NDVI - x1) / (x2  -x1))^2


phi_raster = max(min((LSTi.max - LST)/(LSTi.max - y2) * (phi_max - phi_min_raster) + phi_min_raster, 1.26),0.1)


#---------Calculation of evaporative fraction (EF)----------

#first calculation of Delta= the slope of the saturated vapor pressure at given temperature, as defined by Murray (1967), suggested by Batra et al. (2006)

#input Air temp of EC data 2

T_EC_K =as.numeric(EC_overfly_hour_2$air_temperature) # in Kelvin

# important use Temperature in Kelvin
delta= (26297.76/(T_EC_K-29.65)^2)*exp((17.67*(T_EC_K-273.15)/(T_EC_K-29.65)))

#second: calculation of y with P(atmospheric pressure) and ?? (latent heat vaporization) and Cp (specific heat for dry air) following Wang, Li & Cribb (2006):

Cp=1012 #[Jkg-1K-1]
P0=1013.15 #[hPa] standard atmospheric pressure at sea surface level 
P=P0*10^((-z/18.400)*(T_EC_K/273))
lambda=4.2*(597-(0.6*(T_EC_K-273)))*1000

Gamma= (Cp*P)/(0.622*lambda)

EF = phi_raster*(delta/(delta+Gamma))


# ------ Calculation of ET -----

Rn = as.numeric(EC_overfly_hour_2$NET_rad)

#G = as.numeric(0.07*Rn)

G = 0.583*exp(-2.13*NDVI)*Rn

LE = phi_raster*((Rn-G)*(delta/(delta+Gamma)))



# ----- save raster ----
setwd(plot_path)
writeRaster(LE,filename=paste("rasterLE",image_date,'cargo_G', sep='_'),overwrite=F,format="GTiff")


## ----- calculate ET in [mm/h]------

t=3600 #timespan in seconds - 1h = 3600s

ET_mm_h = ((LE*t)/1000000)/(2.501-(0.002361*(T_EC_K-273.15)))



writeRaster(ET_mm_h,filename=paste("rasterET",image_date,'cargo_G', sep='_'),overwrite=T,format="GTiff")

writeRaster(G,filename=paste("G",image_date,'cargo_G', sep='_'),overwrite=T,format="GTiff")
