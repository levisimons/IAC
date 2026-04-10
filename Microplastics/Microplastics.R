rm(list=ls())
require(data.table)
require(dplyr)
require(FactoMineR)
require(vcd)
require(factoextra)
require(terrainr)
require(sf)
require(raster)
require(FedData)
require(whitebox)
require(exactextractr)
require(terra)
require(prism)
require(leaflet)
require(leaflegend)

#Set working directory
wd <- ""
setwd(wd)

#Read in microplastics data
microplastics_input <- fread(input="Microplastics.csv",sep=",")
microplastics_input[microplastics_input==""] <- NA

#Read in microplastics metadata
microplastics_metadata <- fread(input="Microplastics_metadata.csv",sep=",")
#Force longitude to being numeric
microplastics_metadata$Longitude <- as.numeric(gsub("[^0-9.-]", "", microplastics_metadata$Longitude))

#Create a unique location column
microplastics_metadata <- microplastics_metadata %>%
  dplyr::group_by(Latitude, Longitude) %>%
  dplyr::mutate(Location = cur_group_id())

#Get land usage fractions for each watershed around each sampling point.
NLCD_summary <- c()
for(i in 1:length(unique(microplastics_metadata$Location))){
  #Make a spatial points object of sample site
  tmp <- microplastics_metadata[microplastics_metadata$Location==unique(microplastics_metadata$Location)[i],c("Longitude","Latitude")]
  tmp <- tmp[!duplicated(tmp),]
  sample_site <- st_as_sf(tmp,coords = c("Longitude", "Latitude"), crs = 4326)
  st_write(sample_site,"sample_site.shp",append=F)
  
  #Get sampling centroid
  sampling_centroid <- sample_site
  #Define a box with 100km buffer around the sampling centroid
  sampling_region <- set_bbox_side_length(sampling_centroid, 5000)  # 5km buffer
  
  #Get 30m digital elevation model for sample locations
  DEM_tiles <- get_tiles(sampling_region, services = "elevation", resolution = 30)  
  
  #Merge DEM tiles into a single layer
  DEM <- rast(merge_rasters(DEM_tiles[["elevation"]], tempfile(fileext = ".tif")))
  writeRaster(DEM,"dem_30m.tif",overwrite=T)
  
  #Get 2019 NLCD data
  NLCD <- get_nlcd(
    template = sampling_region,
    label = "my_region",
    year = 2019,           # Available years: 2001, 2004, 2006, 2008, 2011, 2016, 2019
    dataset = "landcover", # or "impervious", "canopy"
    landmass = "L48",      # Lower 48 US states (default)
    extraction.dir = "./nlcd_data"
  )
  
  #Hydrologic preprocessing (DEM → flow)
  whitebox::wbt_fill_depressions(dem=paste(wd,"/dem_30m.tif",sep=""),output=paste(wd,"/dem_filled.tif",sep=""))
  
  #Flow direction raster generation.
  whitebox::wbt_d8_pointer(dem = "dem_filled.tif", output = "flow_dir.tif")
  
  #Generate watershed raster
  whitebox::wbt_watershed(d8_pntr = "flow_dir.tif",pour_pts="sample_site.shp",output="watershed.tif")
  
  #Generate watershed polygon
  watershed <- as.polygons(rast("watershed.tif"), dissolve = TRUE)
  #Transform watershed polygon CRS to match the NLCD one.
  watershed <- project(watershed,rast(NLCD))
  
  #Clipe NLCD raster to polygon
  NLCD_watershed <- mask(crop(NLCD, watershed), watershed)
  
  #Get fraction of NLCD classes per watershed
  NLCD_stats <- as.data.frame(table(as.data.frame(NLCD_watershed)$Class))
  NLCD_stats$percent <- 100*(NLCD_stats$Freq/sum(NLCD_stats$Freq))
  colnames(NLCD_stats) <- c("Class","Freq","percent")
  NLCD_stats <- NLCD_stats[,c("Class","percent")] %>%
    pivot_wider(
      names_from = Class,
      values_from = percent
    )
  #Add in Location
  NLCD_stats$Location <- unique(microplastics_metadata$Location)[i]
  
  #Store results
  NLCD_summary[[i]] <- NLCD_stats
  
  print(paste(i,length(unique(microplastics_metadata$Location))))
}
NLCD_summary <- rbindlist(NLCD_summary)
#Convert NaN values to NA
NLCD_summary[, (names(NLCD_summary)) := lapply(.SD, function(x) {
  if (is.numeric(x)) {
    x[is.nan(x)] <- NA
  }
  x
})]
#Keep only columns that contain at least one value NOT 0 and NOT NA
NLCD_summary <- NLCD_summary[, which(sapply(NLCD_summary, function(x) any(!is.na(x) & x != 0))), with = FALSE]

#Set points with missing NLCD data to just being 100% open water and 0% for everything else.
NLCD_summary[is.na(NLCD_summary[["Open Water"]]), `Open Water` := 100]
NLCD_summary[NLCD_summary[["Open Water"]] == 100, (names(NLCD_summary)) := lapply(.SD, function(x) fifelse(is.na(x), 0, x))]

#Get land class names
NLCD_classes <- colnames(NLCD_summary)[!(colnames(NLCD_summary) %in% c("Location"))]

#Merge in NLCD values into microplastics metadata
microplastics_metadata <- dplyr::left_join(microplastics_metadata,NLCD_summary,relationship="many-to-many")

#Set PRISM working directory for rainfall data
prism_set_dl_dir(paste(wd,"/PRISM",sep=""))

#Get local leading monthly rainfall totals around each sampling point.
microplastics_metadata$`Sample Collection Date` <- as.Date(microplastics_metadata$`Sample Collection Date`,"%m/%d/%y")
PRISM_summary <- c()
k <- 1
for(i in 1:length(unique(microplastics_metadata$`Sample Collection Date`))){
  sample_date <- unique(microplastics_metadata$`Sample Collection Date`)[i]
  #Get monthly rainfall data
  suppressWarnings(get_prism_monthlys(
    type = "ppt",
    years = as.numeric(format(sample_date,"%Y")),
    mon = as.numeric(format(sample_date,"%m")),
    keepZip = F,
    resolution = "4km"
  ))
  #Get rainfall file name
  PRISM_file <- prism_archive_ls() %>%
    grep("ppt", ., value = TRUE)
  #Get rainfall raster
  PRISM_raster <- rast(paste(wd,"/PRISM/",PRISM_file,"/",PRISM_file,".bil",sep=""))
  #Extract monthly rainfall data at sample locations
  tmp <- microplastics_metadata[microplastics_metadata$`Sample Collection Date`==sample_date,]
  for(j in unique(tmp$Location)){
    single_sample <- tmp[tmp$Location==j,c("Location","Latitude","Longitude")]
    single_sample <- single_sample[!duplicated(single_sample),]
    sample_site <- st_as_sf(single_sample,coords = c("Longitude", "Latitude"), crs = 4326)
    sample_site <- st_transform(sample_site, crs(PRISM_raster))
    PRISM_extracted <- data.frame(matrix(nrow=1,ncol=3))
    colnames(PRISM_extracted) <- c("Location","Sample Collection Date","Rainfall")
    PRISM_extracted$Location <- single_sample$Location
    PRISM_extracted$`Sample Collection Date` <- sample_date
    PRISM_extracted$Rainfall <- terra::extract(PRISM_raster, vect(sample_site))[, 2]
    PRISM_summary[[k]] <- PRISM_extracted
    k <- k+1
  }
}
PRISM_summary <- rbindlist(PRISM_summary)

#Merge in rainfall values into microplastics metadata
microplastics_metadata <- dplyr::left_join(microplastics_metadata,PRISM_summary,relationship="many-to-many")

#Merge metadata and sample data
microplastics <- dplyr::left_join(microplastics_input,microplastics_metadata,by=c("Environment","Media","System","Site","Sample ID"),relationship="many-to-many")

#Retain only plastics data
microplastics_filtered <- microplastics[microplastics$`Plastic or Not`=="plastic" & !is.na(microplastics$Latitude),]
#Create a numerical date variable
microplastics_filtered$`Sampling Date` <- as.Date(microplastics_filtered$`Sampling Date`, format = "%m/%d/%y")
microplastics_filtered$date_num <- as.numeric(microplastics_filtered$`Sampling Date`)
#Set certain variables as factors
microplastics_filtered$Color <- as.factor(microplastics_filtered$Color)
microplastics_filtered$Morphology <- as.factor(microplastics_filtered$Morphology)
microplastics_filtered$`1st Order Morph` <- as.factor(microplastics_filtered$`1st Order Morph`)
microplastics_filtered$Polymer <- as.factor(microplastics_filtered$Polymer)
microplastics_filtered$Media <- as.factor(microplastics_filtered$Media)
#Force certain variables to be numeric
microplastics_filtered$`Projected Surface Area (mm2)` <- as.numeric(microplastics_filtered$`Projected Surface Area (mm2)`)
microplastics_filtered$Longitude <- as.numeric(gsub("[^0-9.-]", "", microplastics_filtered$Longitude))

#Retain only variables of interest
microplastics_subset <- microplastics_filtered[,c("date_num","Location","Environment","Media","Discharge (cms)",..NLCD_classes,"Rainfall","Color","Morphology","1st Order Morph","Polymer","Projected Surface Area (mm2)","Feret Length (mm)","Min Feret (mm)","Fiber Length (mm)")]
microplastics_subset <- microplastics_subset[complete.cases(microplastics_subset)]

#Factorial Analysis of Mixed Data (FAMD)
microplastics_FAMD <- FAMD(microplastics_subset)

#Variable contributions to eigenvector dimensions
microplastics_FAMD$var$contrib

#How well each variable is represented?
#0 = poor representation, 1 = excellent representation.
microplastics_FAMD_cor <- data.frame(microplastics_FAMD$var$cos2)

#Map microplastics samples
leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircleMarkers(
    lng=microplastics$Longitude,
    lat=microplastics$Latitude,
    radius = 5,
    fillColor = "orange",
    fillOpacity = 0.33,
    stroke = TRUE,
    color = "black", # Border color
    weight = 1
  )

#Test how strongly correlated certain variables are to microplastics concentration
microplastics_stat_test <- microplastics_filtered[,c("MP Concentration","date_num","Discharge (cms)",..NLCD_classes,"Rainfall")]
microplastics_stat_test <- microplastics_stat_test[!duplicated(microplastics_stat_test),]
microplastics_stat_summary <- c()
i <- 1
for(var in c("date_num","Discharge (cms)",NLCD_classes,"Rainfall")){
  tmp <- data.frame(matrix(nrow=1,ncol=3))
  colnames(tmp) <- c("Variable","rho","p")
  tmp$Variable <- var
  tmp$rho <- cor.test(microplastics_stat_test$`MP Concentration`,microplastics_stat_test[[var]],method="spearman")$estimate
  tmp$p <- cor.test(microplastics_stat_test$`MP Concentration`,microplastics_stat_test[[var]],method="spearman")$p.value
  microplastics_stat_summary[[i]] <- tmp
  i <- i+1
}
microplastics_stat_summary <- rbindlist(microplastics_stat_summary)
