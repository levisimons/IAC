rm(list=ls())
require(data.table)
require(sf)
require(geodata)
require(ggplot2)
require(leaflet)
require(raster)
require(lubridate)
require(dismo)
require(randomForest)
require(DescTools)
require(plyr)

#Set random number string
set.seed(1)

#Set working directory
wd <- ""
setwd(wd)

#Set sf settings to reduce risk of errors
sf_use_s2(FALSE)

#Read in sample data
sample_input <- fread(input="EFA_XRF-demo.csv",sep=",")
#Set sample dates
sample_input$Date <- mdy(sample_input$Date)
#Set days since the Eaton fire started
sample_input$Days <- as.numeric(sample_input$Date-ymd("2025-01-07"))
#Set as spatial object
sample_input <- st_as_sf(sample_input,coords=c("lon","lat"),crs=4326)
#Get bounding box of sample points
sample_bounds <- st_as_sfc(st_bbox(sample_input))
#Set CRS to EPSG:3857 temporarily to use a system in meters
sample_bounds <- st_transform(sample_bounds,crs=3857)
#Create a 5 km buffer around the sampling area
sample_area <- st_bbox(c(xmin=as.numeric(st_bbox(sample_bounds)["xmin"]-5000),
                         ymin=as.numeric(st_bbox(sample_bounds)["ymin"]-5000),
                         xmax=as.numeric(st_bbox(sample_bounds)["xmax"]+5000),
                         ymax=as.numeric(st_bbox(sample_bounds)["ymax"]+5000)),crs=3857)
sample_area <- st_as_sfc(st_bbox(sample_area))
#Transform buffered area back to EPSG:4326
study_area <- st_as_sf(st_transform(sample_area,crs=4326))

#Read in Eaton fire boundary
Eaton <- st_read("Eaton_Perimeter_20250121.shp")

#Set CRS
Eaton <- st_transform(Eaton,crs=4326)

#Input LA county parcel data from https://data.lacounty.gov/datasets/lacounty::assessor-parcel-data-rolls-2021-present/about
Parcels <- fread(input="Parcel_Data_2021_Table_8468414499611436475.csv",sep=",")

#Perform geographic filtering on parcels data
Parcels <- Parcels[Parcels$`Location Latitude`>=st_bbox(study_area)[2] &
                     Parcels$`Location Latitude`<=st_bbox(study_area)[4] &
                     Parcels$`Location Longitude`>=st_bbox(study_area)[1] &
                     Parcels$`Location Longitude`<=st_bbox(study_area)[3],]

#Convert to a spatial object
Parcels <- st_as_sf(Parcels,coords=c("Location Longitude","Location Latitude"),crs=4326)

#Retain columns of interest from parcel data.
Parcel_columns <- c("Number of Buildings","Year Built","Square Footage","Number of Bedrooms","Number of Bathrooms","Number of Units","Total Value")
Parcel_subset <- Parcels[,Parcel_columns]

#Obtain bioclimatic data
region_bioclim <- worldclim_tile(var = "bio", res = 0.5, lon = mean(st_bbox(study_area)[c(1,3)]), lat = mean(st_bbox(study_area)[c(2,4)]),path=".")
region_bioclim <- crop(region_bioclim,study_area)
region_bioclim <- mask(region_bioclim,study_area)
#Extract bioclimatic data
bioclim_points <- raster::extract(region_bioclim,Parcel_subset)
bioclim_points$ID <- NULL

#Obtain precipitation data
region_precip <- worldclim_tile(var = "prec", res = 0.5, lon = mean(st_bbox(study_area)[c(1,3)]), lat = mean(st_bbox(study_area)[c(2,4)]),path=".")
region_precip <- crop(region_precip,study_area)
region_precip <- mask(region_precip,study_area)
#Extract precipitation data
precip_points <- raster::extract(region_precip,Parcel_subset)
precip_points$ID <- NULL

#Obtain solar radiation data
region_rad <-  worldclim_tile(var = "srad", res = 0.5, lon = mean(st_bbox(study_area)[c(1,3)]), lat = mean(st_bbox(study_area)[c(2,4)]),path=".")
region_rad <- crop(region_rad,study_area)
region_rad <- mask(region_rad,study_area)
#Extract solar radiation data
rad_points <- raster::extract(region_rad,Parcel_subset)
rad_points$ID <- NULL

#Obtain wind data
region_wind <-  worldclim_tile(var = "wind", res = 0.5, lon = mean(st_bbox(study_area)[c(1,3)]), lat = mean(st_bbox(study_area)[c(2,4)]),path=".")
region_wind <- crop(region_wind,study_area)
region_wind <- mask(region_wind,study_area)
#Extract wind data
wind_points <- raster::extract(region_wind,Parcel_subset)
wind_points$ID <- NULL

#Get Wildland Urban Interface raster of North America from https://geoserver.silvis.forest.wisc.edu/geodata/globalwui/NA.zip
#Obtain just the portion of the data from the Eaton fire.
if(!file.exists("WUI_Eaton.tif")){
  WUI_tiles <- list.dirs(paste(wd,"/NA",sep=""))[grepl(paste(wd,"/NA/X",sep=""),list.dirs(paste(wd,"/NA",sep="")))]
  i=1
  for(WUI_tile in WUI_tiles){
    # Obtain the bounding box in EPSG:4326 for the raster tile
    raster_bounds <- st_bbox(st_transform(st_as_sfc(st_bbox(rast(paste(WUI_tile,"/WUI.tif",sep=""))), crs = st_crs(rast(paste(WUI_tile,"/WUI.tif",sep="")))), 4326))
    print(paste(i,length(WUI_tiles),WUI_tile,as.numeric(raster_bounds)))
    
    if(raster_bounds[1]<=st_bbox(study_area)[3] & raster_bounds[3]>=st_bbox(study_area)[1] &
       raster_bounds[2]<=st_bbox(study_area)[4] & raster_bounds[4]>=st_bbox(study_area)[2]){
      input_raster <- raster(paste(WUI_tile,"/WUI.tif",sep=""))
      output_raster <- crop(input_raster,st_transform(study_area, crs=st_crs(rast(paste(WUI_tile,"/WUI.tif",sep="")))))
      output_raster <- mask(output_raster,st_transform(study_area, crs=st_crs(rast(paste(WUI_tile,"/WUI.tif",sep="")))))
      output_raster <- terra::project(rast(output_raster), "EPSG:4326")
      writeRaster(output_raster,"WUI_Eaton.tif",overwrite=T)
    }
    i=i+1
  }
}
#Extract WUI data
region_WUI <- raster("WUI_Eaton.tif")
WUI_points <- as.data.frame(terra::extract(region_WUI,Parcel_subset,method="simple"))
colnames(WUI_points) <- c("WUI")
#Round the returned numbers to make sure they correspond to the most common pixel value by the extraction point
WUI_points$WUI <- round(WUI_points$WUI)

#Get wildfire hazard potential for California from https://usfs-public.box.com/shared/static/et4mghz8sq0kxsk2ag6fpey5uec4f8fn.zip 
#Clip down to Eaton.
if(!file.exists("WHP_Eaton.tif")){
  input_raster <- raster("WHP_CA.tif")
  output_raster <- crop(input_raster,st_transform(study_area, crs=st_crs(rast("WHP_CA.tif"))))
  output_raster <- mask(output_raster,st_transform(study_area, crs=st_crs(rast("WHP_CA.tif"))))
  output_raster <- terra::project(rast(output_raster), "EPSG:4326")
  writeRaster(output_raster,"WHP_Eaton.tif",overwrite=T)
}
#Read in Eaton area WHP
region_WHP <- raster("WHP_Eaton.tif")
#Extract WHP data
WHP_points <- as.data.frame(terra::extract(region_WHP,Parcel_subset,method="simple"))
colnames(WHP_points) <- c("WHP")

#Get elevation data.
elevation_input <- elevation_3s(mean(st_bbox(study_area)[1],st_bbox(study_area)[1]),mean(st_bbox(study_area)[2],st_bbox(study_area)[4]),path=".")
#Clip elevation to study area.
region_elevation <- crop(elevation_input,study_area)
region_elevation <- mask(region_elevation,study_area)
#Calculate slope
region_slope <- terra::terrain(region_elevation,v="slope",unit="radians")
#Calculate aspect
region_aspect <- terra::terrain(region_elevation,v="aspect",unit="radians")
#Extract topographic data
elevation_points <- as.data.frame(terra::extract(region_elevation,Parcel_subset,method="simple"))
elevation_points$ID <- NULL
colnames(elevation_points) <- c("elevation")
slope_points <- as.data.frame(terra::extract(region_slope,Parcel_subset,method="simple"))
slope_points$ID <- NULL
colnames(slope_points) <- c("slope")
aspect_points <- as.data.frame(terra::extract(region_aspect,Parcel_subset,method="simple"))
aspect_points$ID <- NULL
colnames(aspect_points) <- c("aspect")

#Combine parcel and environmental data
env_data <- cbind(Parcel_subset,bioclim_points,precip_points,rad_points,wind_points,WUI_points,WHP_points,elevation_points,slope_points,aspect_points)
#Filter out parcels with missing data
env_data <- env_data[env_data$`Number.of.Buildings`>0 & env_data$`Year.Built`>0 & env_data$`Square.Footage`>0,]
#Store parcel coordinates
env_coordinates <- as.data.frame(st_coordinates(env_data))
colnames(env_coordinates) <- c("longitude","latitude")
#Remove points with missing extracted data, but retain remaining points as spatial objects
env_data <- st_drop_geometry(env_data)
env_data <- cbind(env_data,env_coordinates)
env_data <- env_data[complete.cases(env_data),]
env_data <- st_as_sf(env_data,coords=c("longitude","latitude"),crs=4326)

#Identify collinear variables
tmp <- st_drop_geometry(env_data[,!colnames(env_data) %in% c("WUI")])
tmp <- cor(tmp)
tmp[upper.tri(tmp)] <- 0
diag(tmp) <- 0
env_data_filtered <- env_data[, !apply(tmp, 2, function(x) any(abs(x) > 0.70, na.rm = TRUE))]   
env_vars_filtered <- colnames(st_drop_geometry(env_data_filtered))
#Retain only non-collinear variables
env_data <- env_data[,env_vars_filtered]

#Set CRS to EPSG:3857 temporarily to use a system in meters
sample_input <- st_transform(sample_input, crs = 3857)
#Create buffers around sample points of a set distance
sample_buffer <- st_buffer(sample_input,dist=1400)
#Set CRS back to EPSG:4326
sample_input <- st_transform(sample_input, crs = 4326)
sample_buffer <- st_transform(sample_buffer,crs=4326)
#Get the parcel points within the sample buffers
sampled_env <- st_intersection(sample_buffer,env_data)

#Store sample coordinates
sampled_coordinates <- as.data.frame(cbind(sample_input$EFA.ID,st_coordinates(sample_input)))
sampled_coordinates <- sampled_coordinates[!duplicated(sampled_coordinates),]
colnames(sampled_coordinates) <- c("EFA.ID","longitude","latitude")

#Get the mean of environmental and parcel metrics
sampled_env_summary <- st_drop_geometry(sampled_env) %>%
  dplyr::group_by(EFA.ID) %>%
  dplyr::summarise(dplyr::across(-c(No., Date, Days, WUI), median, na.rm = TRUE))

#Add back in days post-Eaton fire
sampled_days <- st_drop_geometry(sampled_env[,c("EFA.ID","Days")])
sampled_days <- sampled_days[!duplicated(sampled_days),]
sampled_env_summary <- dplyr::left_join(sampled_env_summary,sampled_days)

#Define a function to get the mode of a group of values.
getmode_na <- function(v, na.rm = TRUE) {
  if (na.rm) {
    v <- na.omit(v)
  }
  uniqv <- unique(v)
  if (length(uniqv) == 0) { # Handle case where all values are NA
    return(NA)
  }
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
#Get the mode of the WUI values within each sampling buffer
WUI_mode <- sampled_env %>%
  dplyr::group_by(EFA.ID) %>%
  dplyr::summarize(WUI_mode = getmode_na(WUI))
WUI_mode <- st_drop_geometry(WUI_mode)

#Add in the most common WUI value within the sampling buffers.
sampled_env_summary <- dplyr::left_join(sampled_env_summary,WUI_mode)

#Prepare modeling data for random forest
sampled_env_summary <- sampled_env_summary[,!colnames(sampled_env_summary) %in% c("O","Mg","Al","Si","P","S","K","Ca","Ti","Mn","Fe","Ni","Cu","Zn","Ga","As","Rb","Sr","Y","Zr","Nb")]
sampled_env_summary$WUI_mode <- as.factor(sampled_env_summary$WUI_mode)


#Run random forest models and evaluate them
i=1
i_max <- 1000
j=1
stat_tests <- c()
partial_plot_list <- c()
importance_list <- c()
if("WUI" %in% env_vars_filtered){env_vars_filtered <- gsub("WUI","WUI_mode",env_vars_filtered)}
for(i in 1:i_max){
  sampled_env_train <- sampled_env_summary[sample(nrow(sampled_env_summary),0.8*nrow(sampled_env_summary)),]
  sampled_env_test <- sampled_env_summary[!sampled_env_summary$EFA.ID %in% sampled_env_train$EFA.ID,]
  #Run a random forest model over this data subset.
  rf1 <- suppressWarnings(tuneRF(x=sampled_env_train[,!(colnames(sampled_env_train) %in% c("EFA.ID","Pb"))],y=sampled_env_train$Pb,stepFactor=1,plot=FALSE,doBest=TRUE))
  #Predict lead concentrations
  predicted_Pb <- predict(rf1,newdata=sampled_env_test[,!(colnames(sampled_env_test) %in% c("EFA.ID","Pb"))])
  #Store correlations and significance values for model evaluations
  tmp <- data.frame(matrix(nrow=1,ncol=2))
  colnames(tmp) <- c("rho","p")
  tmp$rho <- cor.test(sampled_env_test$Pb,predicted_Pb,method="spearman")$estimate
  tmp$p <- cor.test(sampled_env_test$Pb,predicted_Pb,method="spearman")$p.value
  stat_tests[[i]] <- tmp
  
  #Store relative importance of variable outputs as a temporary data frame.
  tmp <- as.data.frame(rf1$importance)
  #Set one column to store the variable names from the row names.
  tmp$VariableName <- rownames(tmp)
  #Store this importance data frame in the importance list.
  importance_list[[i]] <- tmp
  
  i=i+1
  
  #Store significant partial plot results
  if(cor.test(sampled_env_test$Pb,predicted_Pb,method="spearman")$p.value<=0.05){
    #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
    for(environmental_layer in env_vars_filtered){
      #Store partial plot chart data in a temporary data frame.
      tmp <- as.data.frame(partialPlot(rf1,as.data.frame(sampled_env_test[,!(colnames(sampled_env_test) %in% c("Pb","EFA.ID"))]),x.var=c(environmental_layer),plot=F))
      #Rename probability column
      colnames(tmp) <- c(environmental_layer,"Pb predicted")
      #Store partial plot data in a list of data frames.
      partial_plot_list[[j]] <- tmp
      j <- j+1
    }
  }
}
stat_tests <- rbindlist(stat_tests)
FisherZInv(mean(FisherZ(stat_tests[p<=0.05,rho])))
FisherZInv(sd(FisherZ(stat_tests[p<=0.05,rho])))
nrow(stat_tests[p<=0.05,])/nrow(stat_tests)


#Convert list of importance data frames to a single data frame.
importance_total <- rbind.fill(importance_list)
#Calculate the mean relative importance for each variable.
importance_total <- aggregate(x=importance_total$IncNodePurity,by = list(importance_total$VariableName),FUN = mean)
#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")
#Convert importance to rank importance.
importance_total$Importance <- rank(desc(importance_total$Importance))
#Save rank importance table.
write.table(importance_total,"Eaton_importance_uncalibrated_rank_importance.txt",quote=FALSE,sep="\t",row.names = FALSE)


#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,"Eaton_uncalibrated_partial_plots.txt",quote=FALSE,sep="\t",row.names = FALSE)
partial_plots <- read.table("Eaton_uncalibrated_partial_plots.txt", header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")
#Set WUI mode to categorical for plotting
partial_plots$WUI_mode <- as.factor(partial_plots$WUI_mode)

k <- 17
#Plot heat maps for continuous data.
if(!is.factor(partial_plots[,env_vars_filtered[k]])){
  ggplot(partial_plots, aes(x=!!as.name(env_vars_filtered[k]), y=`Pb predicted`) )+
    xlab(env_vars_filtered[k])+ylab("Predicted Pb")+
    geom_bin2d(bins = 50)+
    scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",nrow(stat_tests[p<=0.05,])," models)",sep=""))+
    stat_smooth(aes(y = `Pb predicted`, fill=`Pb predicted`),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,env_vars_filtered[k]])))+
    theme_bw(base_size=25)
}
#Plot violin plots for categorical data.
if(is.factor(partial_plots[,env_vars_filtered[k]])){
  ggplot(partial_plots, aes(x=!!as.name(env_vars_filtered[k]), y=`Pb predicted`) )+
    xlab(env_vars_filtered[k])+ylab("Predicted Pb")+
    geom_violin(aes(x=!!as.name(env_vars_filtered[k]), y=`Pb predicted`))+
    theme_bw(base_size=25)
}


#Map environmental clusters.
sample_map <- dplyr::left_join(sampled_env_summary,sampled_coordinates)
sample_map <- st_as_sf(sample_map,coords=c("longitude","latitude"),crs=4326)

pal <- colorNumeric(palette = "plasma", domain = sample_map$Pb)
leaflet(sample_map) %>%
  addTiles() %>%
  addCircleMarkers(
    radius = 5,
    fillColor = ~pal(Pb), # Color points by category
    fillOpacity = 0.8,
    stroke = TRUE,
    color = "black", # Border color
    weight = 1,
    popup = ~paste("Mean Pb:", Pb) # Optional: add popup info
  ) %>%
  addPolygons(data=Eaton,color = "blue", weight = 0.1, smoothFactor = 0.1, opacity = 1, fillOpacity = 0.1, fillColor = "blue", highlightOptions = highlightOptions(color = "white", weight = 2, sendToBack = TRUE)) %>%
  addLegend(pal = pal, values = ~Pb, opacity = 0.7, title = "Mean Pb",
            position = "bottomright")
