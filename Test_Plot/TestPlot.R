rm(list=ls())
require(data.table)
require(sf)
require(dplyr)
require(rgbif)
require(tidyr)
require(vegan)
require(ggplot2)
require(iNEXT)

#Set random number string
set.seed(1)

#Set working directory
wd <- ""
setwd(wd)

#Set sf settings to reduce risk of errors
sf_use_s2(FALSE)

#Read in geospatial data from Test Plot
plot_input <- do.call(
  rbind,
  lapply(st_layers("TEST_PLOT.kml")$name, function(lyr) st_read("TEST_PLOT.kml", layer = lyr, quiet = TRUE))
)
#Only retain polygon or multipolygon boundaries.
test_plots <- plot_input[st_geometry_type(plot_input) %in% c("POLYGON", "MULTIPOLYGON"),]
#Set description name to group test plot areas.
test_plots$Description <- ifelse(
  grepl("-", test_plots$Name),
  sub(".*-\\s*", "", test_plots$Name),
  test_plots$Name
)
test_plots$Description <- ifelse(
  grepl("RAINBOW ", test_plots$Description),
  sub("RAINBOW ", "RAINBOW", test_plots$Description),
  test_plots$Description
)
#Remove planned or northern California test plots
test_plots <- test_plots[!(test_plots$Description %in% c("OHLONE HILLSIDE","STAR KING","STARR KING","SAN BRUNO","(NEW PLOT) JEFFERSON BLVD")),]

#Get area of plots
test_plots$area <- as.numeric(st_area(test_plots))

#Get polygon buffers to use to query GBIF
for(plot_boundary in unique(test_plots$Description)){
  #Get test plot boundaries
  test_plot <- test_plots[test_plots$Description==plot_boundary,]
  #Merge plots for each restoration site
  test_plot <- test_plot |>
    group_by(Description) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
  test_plot <- st_as_sf(test_plot)
  #Use a 300m buffer for the test plot polygons
  test_plot_buffer <- test_plot |>
    st_transform(3857) |>      
    st_buffer(300) |>
    st_transform(4326)
  test_plot_buffer <- test_plot_buffer |>
    st_make_valid()
  print(paste(plot_boundary,st_as_text(st_geometry(st_as_sfc(st_bbox(test_plot_buffer))))))
}

#GBIF.org (1 January 2026) GBIF Occurrence Download https://doi.org/10.15468/dl.vm32m5
gbif_input <- fread(input="0070373-251120083545085.zip",quote="")
#Set unique sample IDs
gbif_input[, sampleid := .GRP, by = .(eventDate,decimalLongitude,decimalLatitude)]
#Create spatial points object
gbif_spatial <- st_as_sf(gbif_input,coords=c("decimalLongitude","decimalLatitude"),crs=4326)

#California native plant list from https://calscape.org/search?plant=&orderBy=&location_name=Los%20Angeles%20County%2C%20CA%2C%20USA&lat=34.3871821&lng=-118.1122679&page=1&perPage=60&height_from=&height_to=&width_from=&width_to=
ca_plants <- fread(input="LosAngelesNativePlants.csv",sep=",")

#California native animal list from: LASAN (2021) Los Angeles Native Fauna. Los Angeles
ca_animals <- fread(input="LosAngelesNativeAnimals.csv",sep=",")

#Determine if a species is a California native
gbif_spatial$Native_Plant <- ifelse(gbif_spatial$species %in% ca_plants$Botanical_Name,1,0)
gbif_spatial$Native_Animal <- ifelse(gbif_spatial$species %in% ca_animals$species,1,0)

i_max <- 11
j=1
biodiversity <- c()
metadata <- fread(input="Test_Plot_Metadata.csv",sep=',')
for(plot_boundary in unique(test_plots$Description)){
  #Get test plot boundaries
  test_plot <- test_plots[test_plots$Description==plot_boundary,]
  #Add in test plot metadata
  test_plot <- dplyr::left_join(test_plot,metadata,by=c("Description"="Name"))
  #Get start year for restoration
  Start_year <- as.numeric(format(as.Date(unique(test_plot$Start_date), format = "%m/%d/%y"),"%Y"))
  #Merge plots for each restoration site
  test_plot <- test_plot |>
    group_by(Description) |>
    summarise(geometry = st_union(geometry), .groups = "drop")
  test_plot <- st_as_sf(test_plot)
  #Get start year for restoration
  test_plot$Start_year <- Start_year
  #Define a buffer distance
  buffer_distance <- 30
  
  species_input <- c()
  species_spatial <- c()
  test_plot_buffer <- c()
  #Estimate species richness with a Chao1 estimator
  tmp <- data.frame(matrix(nrow=1,ncol=19))
  colnames(tmp) <- c("richness","richness_density","shannon","shannon_density","ca_richness","ca_richness_density","ca_shannon","ca_shannon_density","occurrences","ca_occurrences","inner_area","outer_area","area","kingdom","zone","distance","Year","Restoration_Year","Plot")
  for(i in 1:i_max){
    #Create a set of buffered zones around the test plot, including the test plot.
    test_plot_buffer[[i]] <- test_plot |>
      st_transform(3857) |>      
      st_buffer((i-1)*buffer_distance) |>
      st_transform(4326)
    test_plot_buffer[[i]] <- test_plot_buffer[[i]] |>
      st_make_valid()
    
    #Filter species occurrences by study area
    occ_tmp <- st_filter(gbif_spatial,test_plot_buffer[[i]])
    species_spatial[[i]] <- occ_tmp[!is.na(occ_tmp$species),]
    
    #Join study area data with species occurrences
    species_spatial[[i]] <- st_join(species_spatial[[i]],st_as_sf(test_plot_buffer[[i]]$geometry))
    #Only retain occurrences unique to the buffer ring around the test plot
    if(i>1){
      species_spatial[[i]] <- species_spatial[[i]][!species_spatial[[i]]$occurrenceID %in% species_spatial[[i-1]]$occurrenceID,]
    } else{
      species_spatial[[i]] <- species_spatial[[i]]
    }
    
    #Calculate species richness data over time, test plot, and surrounding areas.
    if(nrow(species_spatial[[i]])>0){
      for(year_selected in min(na.omit(species_spatial[[i]]$year)):max(na.omit(species_spatial[[i]]$year))){
        for(kingdom_selected in c("Animalia","Plantae")){
          tmp$Year <- year_selected
          tmp$Plot <- plot_boundary
          tmp$Restoration_Year <- tmp$Year-test_plot$Start_year
          tmp$zone <- i
          tmp$distance <- (i-1)*buffer_distance
          tmp$kingdom <- kingdom_selected
          
          #Get a count of occurrences per buffer ring
          if(kingdom_selected=="Plantae"){species_count <- as.data.frame(table(st_drop_geometry(species_spatial[[i]][species_spatial[[i]]$Native_Plant==0 & species_spatial[[i]]$year==year_selected & species_spatial[[i]]$kingdom==kingdom_selected,c("scientificName")])))}
          if(kingdom_selected=="Animalia"){species_count <- as.data.frame(table(st_drop_geometry(species_spatial[[i]][species_spatial[[i]]$Native_Animal==0 & species_spatial[[i]]$year==year_selected & species_spatial[[i]]$kingdom==kingdom_selected,c("scientificName")])))}
          
          if(nrow(species_count)>0){
            #Estimate alpha diversity metrics
            result <- iNEXT(species_count$Freq, q = 0, datatype = "abundance")$AsyEst
            
            #Estimate richness for the test plot
            tmp$richness <- result[rownames(result)=="Species Richness","Estimator"]
            #Estimate Shannon diversity for the test plot
            tmp$shannon <- result[rownames(result)=="Shannon diversity","Estimator"]
            #Number of occurrences per area per time
            tmp$occurrences <- sum(species_count$Freq)
            #Estimate chao1 richness density for the test plot and the surrounding rings.
            if(i>1){
              tmp$inner_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i-1]]$geometry,crs=3857)))
              tmp$outer_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i]]$geometry,crs=3857)))
              tmp$area <- tmp$outer_area - tmp$inner_area
              tmp$richness_density <- as.numeric(tmp$richness/(tmp$outer_area-tmp$inner_area))
              tmp$shannon_density <- as.numeric(tmp$shannon/(tmp$outer_area-tmp$inner_area))
            } else{
              tmp$inner_area <- 0
              tmp$outer_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i]]$geometry,crs=3857)))
              tmp$area <- tmp$outer_area - tmp$inner_area
              tmp$richness_density <- as.numeric(tmp$richness/(tmp$outer_area-tmp$inner_area))
              tmp$shannon_density <- as.numeric(tmp$shannon/(tmp$outer_area-tmp$inner_area))
            }
          } else{
            tmp$richness <- NA
            tmp$richness_density <- NA
            tmp$shannon <- NA
            tmp$shannon_density <- NA
          }
          
          #Get a count of California native occurrences per buffer ring
          if(kingdom_selected=="Plantae"){ca_species_count <- as.data.frame(table(st_drop_geometry(species_spatial[[i]][species_spatial[[i]]$Native_Plant==1 & species_spatial[[i]]$year==year_selected & species_spatial[[i]]$kingdom==kingdom_selected,c("scientificName")])))}
          if(kingdom_selected=="Animalia"){ca_species_count <- as.data.frame(table(st_drop_geometry(species_spatial[[i]][species_spatial[[i]]$Native_Animal==1 & species_spatial[[i]]$year==year_selected & species_spatial[[i]]$kingdom==kingdom_selected,c("scientificName")])))}
          
          if(nrow(ca_species_count)>0){
            #Estimate alpha diversity metrics
            result <- iNEXT(ca_species_count$Freq, q = 0, datatype = "abundance")$AsyEst
            
            #Estimate richness for the test plot
            tmp$ca_richness <- result[rownames(result)=="Species Richness","Estimator"]
            #Estimate Shannon diversity for the test plot
            tmp$ca_shannon <- result[rownames(result)=="Shannon diversity","Estimator"]
            #Number of occurrences per area per time
            tmp$ca_occurrences <- sum(ca_species_count$Freq)
            #Estimate richness density for the test plot and the surrounding rings.
            if(i>1){
              tmp$inner_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i-1]]$geometry,crs=3857)))
              tmp$outer_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i]]$geometry,crs=3857)))
              tmp$area <- tmp$outer_area - tmp$inner_area
              tmp$ca_richness_density <- as.numeric(tmp$ca_richness/(tmp$outer_area-tmp$inner_area))
              tmp$ca_shannon_density <- as.numeric(tmp$ca_shannon/(tmp$outer_area-tmp$inner_area))
            } else{
              tmp$inner_area <- 0
              tmp$outer_area <- as.numeric(st_area(st_transform(test_plot_buffer[[i]]$geometry,crs=3857)))
              tmp$area <- tmp$outer_area - tmp$inner_area
              tmp$ca_richness_density <- as.numeric(tmp$ca_richness/(tmp$outer_area-tmp$inner_area))
              tmp$ca_shannon_density <- as.numeric(tmp$ca_shannon/(tmp$outer_area-tmp$inner_area))
            }
          } else{
            tmp$ca_richness <- NA
            tmp$ca_richness_density <- NA
            tmp$ca_shannon <- NA
            tmp$ca_shannon_density <- NA
          }
          
          biodiversity[[j]] <- tmp
          j=j+1
        }
      }
    }
    print(paste(plot_boundary,i))
  }
}
biodiversity <- rbindlist(biodiversity)
biodiversity$Restored <- ifelse(biodiversity$Restoration_Year<0,"Pre-restoration","Post-restoration")
biodiversity$occurrence_density <- biodiversity$occurrences/biodiversity$area
biodiversity$ca_occurrence_density <- biodiversity$ca_occurrences/biodiversity$area
biodiversity_filtered <- biodiversity[!is.na(biodiversity$richness),]
biodiversity_filtered <- as.data.frame(biodiversity_filtered)
ca_biodiversity_filtered <- biodiversity[!is.na(biodiversity$ca_richness),]
ca_biodiversity_filtered <- as.data.frame(ca_biodiversity_filtered)
biodiversity_filtered$Restored <- factor(biodiversity_filtered$Restored, levels = c("Pre-restoration","Post-restoration"))   
ca_biodiversity_filtered$Restored <- factor(ca_biodiversity_filtered$Restored, levels = c("Pre-restoration","Post-restoration"))   

#Plot diversity density versus distance
ggplot(biodiversity_filtered, aes(x = distance, y = richness_density)) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10()+geom_point()+facet_grid(rows=vars(Restored),cols=vars(kingdom)) + labs(x = "distance (m)", y = "Non-native richness density\ninverse square meter")
ggplot(ca_biodiversity_filtered, aes(x = distance, y = ca_richness_density)) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10()+geom_point()+facet_grid(rows=vars(Restored),cols=vars(kingdom)) + labs(x = "distance (m)", y = "Native richness density\ninverse square meter")
ggplot(biodiversity_filtered, aes(x = distance, y = shannon_density)) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10()+geom_point()+facet_grid(rows=vars(Restored),cols=vars(kingdom)) + labs(x = "distance (m)", y = "Non-native Shannon diversity density\ninverse square meter")
ggplot(ca_biodiversity_filtered, aes(x = distance, y = ca_shannon_density)) + geom_smooth(method = "lm", se = FALSE) + scale_y_log10()+geom_point()+facet_grid(rows=vars(Restored),cols=vars(kingdom)) + labs(x = "distance (m)", y = "Native Shannon diversity density\ninverse square meter")

#Check for significant relationships between species density and distance from restoration site
stat_test <- c()
index_selected <- c("richness_density","shannon_density")[2]
i=1
for(kingdom_selected in unique(biodiversity_filtered$kingdom)){
  for(restored_selected in unique(biodiversity_filtered$Restored)){
    tmp <- data.frame(matrix(nrow=1,ncol=4))
    colnames(tmp) <- c("kingdom","Restored","rho","p")
    tmp$kingdom <- kingdom_selected
    tmp$Restored <- restored_selected
    tmp$rho <- cor.test(biodiversity_filtered[biodiversity_filtered$distance >= 0 & biodiversity_filtered$Restored==restored_selected & biodiversity_filtered$kingdom==kingdom_selected,"zone"],biodiversity_filtered[biodiversity_filtered$distance >= 0 & biodiversity_filtered$Restored==restored_selected & biodiversity_filtered$kingdom==kingdom_selected,index_selected],method="spearman")$estimate
    tmp$p <- cor.test(biodiversity_filtered[biodiversity_filtered$distance >= 0 & biodiversity_filtered$Restored==restored_selected & biodiversity_filtered$kingdom==kingdom_selected,"zone"],biodiversity_filtered[biodiversity_filtered$distance >= 0 & biodiversity_filtered$Restored==restored_selected & biodiversity_filtered$kingdom==kingdom_selected,index_selected],method="spearman")$p.value
    stat_test[[i]] <- tmp
    i=i+1
  }
}
stat_test <- rbindlist(stat_test)

#Check for significant relationships between California species density and distance from restoration site
ca_stat_test <- c()
index_selected <- c("ca_richness_density","ca_shannon_density")[2]
i=1
for(kingdom_selected in unique(ca_biodiversity_filtered$kingdom)){
  for(restored_selected in unique(ca_biodiversity_filtered$Restored)){
    tmp <- data.frame(matrix(nrow=1,ncol=4))
    colnames(tmp) <- c("kingdom","Restored","rho","p")
    tmp$kingdom <- kingdom_selected
    tmp$Restored <- restored_selected
    tmp$rho <- cor.test(ca_biodiversity_filtered[ca_biodiversity_filtered$distance >= 0 & ca_biodiversity_filtered$Restored==restored_selected & ca_biodiversity_filtered$kingdom==kingdom_selected,"zone"],ca_biodiversity_filtered[ca_biodiversity_filtered$distance >= 0 & ca_biodiversity_filtered$Restored==restored_selected & ca_biodiversity_filtered$kingdom==kingdom_selected,index_selected],method="spearman")$estimate
    tmp$p <- cor.test(ca_biodiversity_filtered[ca_biodiversity_filtered$distance >= 0 & ca_biodiversity_filtered$Restored==restored_selected & ca_biodiversity_filtered$kingdom==kingdom_selected,"zone"],ca_biodiversity_filtered[ca_biodiversity_filtered$distance >= 0 & ca_biodiversity_filtered$Restored==restored_selected & ca_biodiversity_filtered$kingdom==kingdom_selected,index_selected],method="spearman")$p.value
    ca_stat_test[[i]] <- tmp
    i=i+1
  }
}
ca_stat_test <- rbindlist(ca_stat_test)

#Add in plot areas
tmp <- st_drop_geometry(test_plots[,c("Description","area")])
tmp <- aggregate(area ~ Description, data = tmp, sum, na.rm = TRUE)
biodiversity_filtered <- dplyr::left_join(biodiversity_filtered,tmp)
ca_biodiversity_filtered <- dplyr::left_join(ca_biodiversity_filtered,tmp)

#How does chao1 species density vary with field variables?
anovas <- c()
tmp <- anova(lm(richness_density~Restored+Restoration_Year+distance+Plot+area,data=biodiversity_filtered[biodiversity_filtered$kingdom=="Plantae",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Non-native plants"
anovas[[1]] <- tmp
tmp <- anova(lm(richness_density~Restored+Restoration_Year+distance+Plot+area,data=biodiversity_filtered[biodiversity_filtered$kingdom=="Animalia",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Non-native animals"
anovas[[2]] <- tmp
tmp <- anova(lm(ca_richness_density~Restored+Restoration_Year+distance+Plot+area,data=ca_biodiversity_filtered[ca_biodiversity_filtered$kingdom=="Plantae",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Native plants"
anovas[[3]] <- tmp
tmp <- anova(lm(ca_richness_density~Restored+Restoration_Year+distance+Plot+area,data=ca_biodiversity_filtered[ca_biodiversity_filtered$kingdom=="Animalia",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Native animals"
anovas[[4]] <- tmp
anovas <- rbindlist(anovas)
anovas <- anovas[anovas$factor!="Residuals",]
anovas$`F value` <- ifelse(anovas$`Pr(>F)`>0.05,NA,anovas$`F value`)
#Plot F-values from the anovas with chao1 density
ggplot(anovas, aes(x = category, y = factor)) +
  geom_tile(aes(fill = `F value`), color = "white") +
  scale_fill_viridis_c(name = "F",limits = c(0, 150)) +
  labs(x = "category", y = "factor") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

#How does Shannon species density vary with field variables?
anovas <- c()
tmp <- anova(lm(shannon_density~Restored+Restoration_Year+distance+Plot+area,data=biodiversity_filtered[biodiversity_filtered$kingdom=="Plantae",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Non-native plants"
anovas[[1]] <- tmp
tmp <- anova(lm(shannon_density~Restored+Restoration_Year+distance+Plot+area,data=biodiversity_filtered[biodiversity_filtered$kingdom=="Animalia",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Non-native animals"
anovas[[2]] <- tmp
tmp <- anova(lm(ca_shannon_density~Restored+Restoration_Year+distance+Plot+area,data=ca_biodiversity_filtered[ca_biodiversity_filtered$kingdom=="Plantae",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Native plants"
anovas[[3]] <- tmp
tmp <- anova(lm(ca_shannon_density~Restored+Restoration_Year+distance+Plot+area,data=ca_biodiversity_filtered[ca_biodiversity_filtered$kingdom=="Animalia",]))
tmp$factor <- rownames(tmp)
tmp$category <- "Native animals"
anovas[[4]] <- tmp
anovas <- rbindlist(anovas)
anovas <- anovas[anovas$factor!="Residuals",]
anovas$`F value` <- ifelse(anovas$`Pr(>F)`>0.05,NA,anovas$`F value`)
#Plot F-values from the anovas with chao1 density
ggplot(anovas, aes(x = category, y = factor)) +
  geom_tile(aes(fill = `F value`), color = "white") +
  scale_fill_viridis_c(name = "F",limits = c(0, 150)) +
  labs(x = "category", y = "factor") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

#Compare species density distributions by kingdom and distance from restoration sites.
group_test <- c()
index_selected <- c("richness_density","shannon_density")[2]
i=1
for(distance_selected in min(biodiversity_filtered$distance):max(biodiversity_filtered$distance)){
  for(kingdom_selected in unique(biodiversity_filtered$kingdom)){
    if(length(unique(biodiversity_filtered[biodiversity_filtered$distance==distance_selected & biodiversity_filtered$kingdom==kingdom_selected,"Restored"]))>1){
      tmp <- data.frame(matrix(nrow=1,ncol=4))
      colnames(tmp) <- c("distance","kingdom","kolmogorov_smirnov","p")
      tmp$kingdom <- kingdom_selected
      tmp$distance <- distance_selected
      tmp$kolmogorov_smirnov <- ks.test(biodiversity_filtered[biodiversity_filtered$distance==distance_selected & biodiversity_filtered$kingdom==kingdom_selected & biodiversity_filtered$Restored=="Pre-restoration",index_selected],biodiversity_filtered[biodiversity_filtered$distance==distance_selected & biodiversity_filtered$kingdom==kingdom_selected & biodiversity_filtered$Restored=="Post-restoration",index_selected],alternative = "greater")$statistic
      tmp$p <- ks.test(biodiversity_filtered[biodiversity_filtered$distance==distance_selected & biodiversity_filtered$kingdom==kingdom_selected & biodiversity_filtered$Restored=="Pre-restoration",index_selected],biodiversity_filtered[biodiversity_filtered$distance==distance_selected & biodiversity_filtered$kingdom==kingdom_selected & biodiversity_filtered$Restored=="Post-restoration",index_selected],alternative = "greater")$p.value
      group_test[[i]] <- tmp
      i <- i+1
    }
  }
}
group_test <- rbindlist(group_test)

#Compare California species density distributions by kingdom and distance from restoration sites.
ca_group_test <- c()
index_selected <- c("ca_richness_density","ca_shannon_density")[2]
i=1
for(distance_selected in min(ca_biodiversity_filtered$distance):max(ca_biodiversity_filtered$distance)){
  for(kingdom_selected in unique(ca_biodiversity_filtered$kingdom)){
    if(length(unique(ca_biodiversity_filtered[ca_biodiversity_filtered$distance==distance_selected & ca_biodiversity_filtered$kingdom==kingdom_selected,"Restored"]))>1){
      tmp <- data.frame(matrix(nrow=1,ncol=4))
      colnames(tmp) <- c("distance","kingdom","kolmogorov_smirnov","p")
      tmp$kingdom <- kingdom_selected
      tmp$distance <- distance_selected
      tmp$kolmogorov_smirnov <- ks.test(ca_biodiversity_filtered[ca_biodiversity_filtered$distance==distance_selected & ca_biodiversity_filtered$kingdom==kingdom_selected & ca_biodiversity_filtered$Restored=="Pre-restoration",index_selected],ca_biodiversity_filtered[ca_biodiversity_filtered$distance==distance_selected & ca_biodiversity_filtered$kingdom==kingdom_selected & ca_biodiversity_filtered$Restored=="Post-restoration",index_selected],alternative = "greater")$statistic
      tmp$p <- ks.test(ca_biodiversity_filtered[ca_biodiversity_filtered$distance==distance_selected & ca_biodiversity_filtered$kingdom==kingdom_selected & ca_biodiversity_filtered$Restored=="Pre-restoration",index_selected],ca_biodiversity_filtered[ca_biodiversity_filtered$distance==distance_selected & ca_biodiversity_filtered$kingdom==kingdom_selected & ca_biodiversity_filtered$Restored=="Post-restoration",index_selected],alternative = "greater")$p.value
      ca_group_test[[i]] <- tmp
      i <- i+1
    }
  }
}
ca_group_test <- rbindlist(ca_group_test)
#Plot Kolmogorov-Smirnov test values
group_test$category <- ifelse(group_test$kingdom=="Plantae","Non-native plants","Non-native animals")
ca_group_test$category <- ifelse(ca_group_test$kingdom=="Plantae","Native plants","Native animals")
group_tests <- rbind(group_test,ca_group_test)
group_tests$kolmogorov_smirnov <- ifelse(group_tests$p > 0.05,NA,group_tests$kolmogorov_smirnov)
ggplot(group_tests, aes(x = distance, y = category)) +
  geom_tile(aes(fill = kolmogorov_smirnov), color = "white") +
  scale_fill_viridis_c(name = "Kolmogorov Smirnov\nscore",limits = c(0, 1)) +
  labs(x = "distance (m)", y = "category") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

#Map project data sets
require(terra)
require(leaflet)
require(leaflegend)
require(purrr)
require(mapview)

i <- 1
plot_maps <- c()
plot_buffer_maps <- c()
for(plot_boundary in list.files(path=paste(wd,"/Test_Plot",sep=""))){
  #Get test plot boundaries
  test_plot <- st_read(paste(wd,"/Test_Plot/",plot_boundary,sep=""))
  test_plot <- test_plot |>
    filter(!st_is_empty(geometry)) |>
    st_make_valid() |>
    filter(st_geometry_type(geometry) %in% c("POLYGON", "MULTIPOLYGON")) |>
    st_cast("MULTIPOLYGON") |>
    st_transform(4326)
  #Use a 300m buffer for the test plot polygons
  test_plot_buffer <- test_plot |>
    st_transform(3857) |>      
    st_buffer(300) |>
    st_transform(4326)
  test_plot_buffer <- test_plot_buffer |>
    st_make_valid()
  #Combine test plot boundaries
  plot_maps[[i]] <- test_plot
  plot_buffer_maps[[i]] <- test_plot_buffer
  i <- i+1
}
plot_maps <- dplyr::bind_rows(plot_maps)
plot_buffer_maps <- dplyr::bind_rows(plot_buffer_maps)
leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addPolygons(data=plot_maps,color="red") |>
  addPolygons(data=plot_buffer_maps,color="orange")
