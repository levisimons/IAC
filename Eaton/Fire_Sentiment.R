rm(list=ls())
require(data.table)
require(raster)
require(sf)

#Set working directory
wd <- ""
setwd(wd)

#Read in Palisades fire dNBR
#Obtained via: wget https://popo.jpl.nasa.gov/pub/LA_Fires/dist/dnbr/COG_Palisades_AV3_provisional_dNBR_866nm_2198nm_20240905-20250116.tif
Palisades <- rast("COG_Palisades_AV3_provisional_dNBR_866nm_2198nm_20240905-20250116.tif")
#Reproject raster to EPSG:4326
Palisades <- terra::project(Palisades, "EPSG:4326")

#Read in Eaton fire dNBR
#Obtained via: wget https://popo.jpl.nasa.gov/pub/LA_Fires/dist/dnbr/COG_Eaton_AV3_provisional_dNBR_866nm_2198nm_20240905-20250116.tif
Eaton <- rast("COG_Eaton_AV3_provisional_dNBR_866nm_2198nm_20240905-20250116.tif")
#Reproject raster to EPSG:4326
Eaton <- terra::project(Eaton, "EPSG:4326")

#Read in CSUCI sentiment data from 2025
#Obtained from: https://github.com/levisimons/IAC/blob/main/Eaton/LAStronger-Rhetoric2025.csv
LAStronger_2025 <- fread(input="LAStronger-Rhetoric2025.csv",sep=",")
#Remove entries with missing spatial coordinates
LAStronger_2025 <- LAStronger_2025[!is.na(LAStronger_2025$x) & !is.na(LAStronger_2025$y),]
#Make into a spatial points object
LAStronger_2025 <- st_as_sf(LAStronger_2025,coords=c("x","y"),crs=4326)

#Read in CSUCI sentiment data from 2026
#Obtained from: https://github.com/levisimons/IAC/blob/main/Eaton/LAStronger-Rhetoric2026.csv
LAStronger_2026 <- fread(input="LAStronger-Rhetoric2026.csv",sep=",")
#Remove entries with missing spatial coordinates
LAStronger_2026 <- LAStronger_2026[!is.na(LAStronger_2026$x) & !is.na(LAStronger_2026$y),]
#Make into a spatial points object
LAStronger_2026 <- st_as_sf(LAStronger_2026,coords=c("x","y"),crs=4326)

#Extract Palisades dNBR data for 2025
LAStronger_2025$dNBR_Palisades <- raster::extract(Palisades,LAStronger_2025,df=T)[,2]
#Extract Eaton dNBR data for 2025
LAStronger_2025$dNBR_Eaton <- raster::extract(Eaton,LAStronger_2025,df=T)[,2]
#Create a unified dNBR column
LAStronger_2025$dNBR <-  rowSums(st_drop_geometry(LAStronger_2025[, c("dNBR_Palisades", "dNBR_Eaton")]), na.rm = TRUE)
#Convert 0 to NA in dNBR column
LAStronger_2025$dNBR <- ifelse(LAStronger_2025$dNBR==0,NA,LAStronger_2025$dNBR)
#Add year column
LAStronger_2025$year <- 2025

#Extract Palisades dNBR data for 2025
LAStronger_2026$dNBR_Palisades <- terra::extract(Palisades,LAStronger_2026,df=T)[,2]
#Extract Eaton dNBR data for 2025
LAStronger_2026$dNBR_Eaton <- terra::extract(Eaton,LAStronger_2026,df=T)[,2]
#Create a unified dNBR column
LAStronger_2026$dNBR <-  rowSums(st_drop_geometry(LAStronger_2026[, c("dNBR_Palisades", "dNBR_Eaton")]), na.rm = TRUE)
#Convert 0 to NA in dNBR column
LAStronger_2026$dNBR <- ifelse(LAStronger_2026$dNBR==0,NA,LAStronger_2026$dNBR)
#Add year column
LAStronger_2026$year <- 2026

#Create filtered data set for analysis
Sentiment <- rbind(st_drop_geometry(LAStronger_2025[,c("year","Sign Category Overall Sentiment-Bin 3","Sign Category-Bin1","dNBR")]),st_drop_geometry(LAStronger_2026[,c("year","Sign Category Overall Sentiment-Bin 3","Sign Category-Bin1","dNBR")]))
#Filter data set to only contain complete records]
Sentiment[Sentiment==""] <- NA
Sentiment <- Sentiment[complete.cases(Sentiment),]

#Plot and test sentiment versus burn severity
ggplot(Sentiment, aes(x = `Sign Category Overall Sentiment-Bin 3`, y = dNBR))+
  labs(x = "Sign sentiment category", y = "dNBR (Burn severity)")+
  geom_violin()
kruskal.test(`Sign Category Overall Sentiment-Bin 3`~dNBR,data=Sentiment)
ks.test(Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="negative"]$dNBR,Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="positive"]$dNBR)
ks.test(Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="negative"]$dNBR,Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="neutral"]$dNBR)
ks.test(Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="positive"]$dNBR,Sentiment[Sentiment$`Sign Category Overall Sentiment-Bin 3`=="neutral"]$dNBR)

#Focus on construction signs versus pro-community signs
Sentiment$SignType <- ifelse(grepl("pro-community",Sentiment$`Sign Category-Bin1`,ignore.case=T),"Pro-Community",
                             Sentiment$SignType <- ifelse(grepl("construction",Sentiment$`Sign Category-Bin1`,ignore.case=T),"Construction",NA))

#Plot and test sentiment versus burn severity
ggplot(Sentiment[!is.na(Sentiment$SignType),], aes(x = `SignType`, y = dNBR))+
  labs(x = "Sign sentiment category", y = "dNBR (Burn severity)")+
  geom_violin()
ks.test(Sentiment[!is.na(Sentiment$SignType) & Sentiment$SignType=="Construction"]$dNBR,Sentiment[!is.na(Sentiment$SignType) & Sentiment$SignType=="Pro-Community"]$dNBR,alternative="less")
