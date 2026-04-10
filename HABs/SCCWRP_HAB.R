rm(list=ls())
require(data.table)
require(DescTools)
require(dismo)
require(dplyr)
require(plyr)
require(randomForest)
require(tidyr)
require(leaflet)
require(leaflegend)

#Set random number string
set.seed(1)

#Set working directory
wd <- ""
setwd(wd)

#Get the list of SCCOOS files
input_files <- list.files(path=paste(wd,"/SCCOOS",sep=""),pattern=".csv")
i <- 1
SCCOOS <- c()
for(input_file in input_files){
  #Read in southern California SCCOOS data from https://erddap.sccoos.org/erddap/tabledap/index.html?page=1&itemsPerPage=1000
  SCCOOS[[i]] <- fread(input=paste(wd,"/SCCOOS/",input_file,sep=""))
  i <- i+1
}
#Merge files into a single data frame.
SCCOOS <- as.data.frame(rbindlist(SCCOOS, fill=T))

#Define physicochemical variables
var_chem <- c("Temp (degree_C)","Phosphate (uM)","Silicate (uM)","Nitrite (uM)","Nitrate (uM)","Ammonium (uM)")

#Define biological variables
var_bio <- c("Akashiwo_sanguinea (cells/L)","Alexandrium_spp (cells/L)","Dinophysis_spp (cells/L)","Lingulodinium_polyedra (cells/L)","Prorocentrum_spp (cells/L)","Pseudo_nitzschia_delicatissima_group (cells/L)","Pseudo_nitzschia_seriata_group (cells/L)","Ceratium_spp (cells/L)","Cochlodinium_spp (cells/L)","Gymnodinium_spp (cells/L)","Other_Diatoms (cells/L)","Other_Dinoflagellates (cells/L)","Total_Phytoplankton (cells/L)")

#Store variable importances for all models
importance_total <- c()
#Store model accuracies for all models
accuracy_total <- c()
for(i in 1:length(var_bio)){
  #Select all physicochemical variables and one biological variable at a time for modelling.
  SCCOOS_filtered <- SCCOOS[,c("Location_Code",var_chem,var_bio[i])]
  #Remove entries with missing data
  SCCOOS_filtered <- SCCOOS_filtered[complete.cases(SCCOOS_filtered),]
  
  #Iterate over the number of model iterations
  model_num <- 10
  #Store variable importances
  importance_list <- c()
  #Store model accuracies
  accuracy_list <- c()
  #Run and evaluate models of harmful algal blooms
  for(j in 1:model_num){
    #Create a subset of the SCCOOS data composed of a randomly selected 80% of rows from SCCOOS_filtered.
    #Construct a training and testing set from this subsetted data.
    group <- kfold(SCCOOS_filtered,5)
    SCCOOS_training <- SCCOOS_filtered[group!=1,]
    SCCOOS_testing <- SCCOOS_filtered[group==1,]
    
    #Run a random forest model over this data subset.
    #Inputs are chemical data, output is the selected biological variable
    rf1 <- suppressWarnings(tuneRF(x=SCCOOS_training[,colnames(SCCOOS_training) %in% c("Location_Code",var_chem)],y=SCCOOS_training[,colnames(SCCOOS_training) %in% var_bio],stepFactor=1,plot=FALSE,doBest=TRUE))
    
    #Store relative importance of variable outputs as a temporary data frame.
    tmp <- as.data.frame(rf1$importance)
    #Set one column to store the variable names from the row names.
    tmp$VariableName <- rownames(tmp)
    #Store this importance data frame in the importance list.
    importance_list[[j]] <- tmp
    
    #Calculate the correlation between actual and predicted algal species variable to evaluate model accuracy.
    #Calculate the predicted algal concentration
    prediction <- predict(rf1,SCCOOS_testing)
    #Store significant correlation results
    print(paste(var_bio[i],j))
    if(cor.test(prediction,SCCOOS_testing[,var_bio[i]],method="spearman")$p.value < 0.05){
      accuracy_list[j] <- cor.test(prediction,SCCOOS_testing[,var_bio[i]],method="spearman")$estimate
    }
  }
  #Calculate the mean and standard deviation on correlations between predicted and actual algal concentrations
  accuracy_model <- data.frame(matrix(nrow=1,ncol=3))
  colnames(accuracy_model) <- c("algae","mean_rho","sd_rho")
  accuracy_model$algae <- var_bio[i]
  #Check to see if there's any data to summarize before summarizing it.
  if(length(na.omit(accuracy_list))>0){
    accuracy_model$mean_rho <- FisherZInv(mean(na.omit(FisherZ(accuracy_list))))
    accuracy_model$sd_rho <- FisherZInv(sd(na.omit(FisherZ(accuracy_list))))
  } else{
    accuracy_model$mean_rho <- NA
    accuracy_model$sd_rho <- NA
  }
  
  #Store model accuracies
  accuracy_total[[i]] <- accuracy_model
  
  #Convert list of importance data frames to a single data frame.
  importance_model <- rbind.fill(importance_list)
  #Calculate the mean relative importance for each variable.
  importance_model <- aggregate(x=importance_model$IncNodePurity,by = list(importance_model$VariableName),FUN = mean)
  #Rename columns.
  colnames(importance_model) <- c("VariableName","Importance")
  #Convert importance to rank importance.
  importance_model$Importance <- rank(desc(importance_model$Importance))
  #Add column per algae
  importance_model$algae <- var_bio[i]
  
  #Store model importance values
  importance_total[[i]] <- importance_model
}
#Summarize model accuracies and variable importance values
importance_total <- rbind.fill(importance_total,fill=T)
accuracy_total <- rbind.fill(accuracy_total,fill=T)


#Read in SCCWRP sample data from https://postfire.sccwrp.org/pages/analysis
SCCWRP_input <- as.data.frame(fread(input="SCCWRP_station_data.csv",sep=","))
#Set date variable
SCCWRP_input$sampledate <- as.Date(SCCWRP_input$sampledate,format="%Y-%m-%d")

#Set dates for the fires
fire_start <- as.Date("2025-01-07",format="%Y-%m-%d")
fire_end <- as.Date("2025-01-31",format="%Y-%m-%d")

#Set time from the firest
SCCWRP_input$from_fire_start <- as.numeric(SCCWRP_input$sampledate-fire_start)
SCCWRP_input$from_fire_end <- as.numeric(SCCWRP_input$sampledate-fire_end)

#Set No data values
SCCWRP_input[SCCWRP_input==-88] <- NA

#Select SCCWRP data with the same variables used to model HABs
#Only keep coastal data, that is samples taken from saltwater.
SCCWRP <- SCCWRP_input[SCCWRP_input$matrixname=="saltwater" & SCCWRP_input$analytename %in% c("Temperature","Ammonia","Ammonia as N","Nitrate","Nitrate as N","Nitrite","Nitrite as N","Phosphate","Orthophosphate","Orthophosphate as P","Silicate"),]

#Convert SCCWRP water temperatures
SCCWRP$`Temp (degree_C)` <- ifelse(SCCWRP$analytename %in% c("Temperature"),SCCWRP$result,NA)

#Convert SCCWRP ammonia concentrations
SCCWRP$`Ammonium (uM)` <- ifelse(SCCWRP$analytename %in% c("Ammonia","Ammonia as N"),SCCWRP$result/17.03,NA)

#Convert SCCWRP nitrate concentrations
SCCWRP$`Nitrate (uM)` <- ifelse(SCCWRP$analytename %in% c("Nitrate","Nitrate as N"),SCCWRP$result/62,NA)

#Convert SCCWRP nitrite concentrations
SCCWRP$`Nitrite (uM)` <- ifelse(SCCWRP$analytename %in% c("Nitrite","Nitrite as N"),SCCWRP$result/46.01,NA)

#Convert SCCWRP phosphate concentrations
SCCWRP$`Phosphate (uM)`<- ifelse(SCCWRP$analytename %in% c("Phosphate","Orthophosphate","Orthophosphate as P"),SCCWRP$result/94.97,NA)

#Convert SCCWRP silicate concentrations
SCCWRP$`Silicate (uM)` <- ifelse(SCCWRP$analytename %in% c("Silicate"),SCCWRP$result/60.08,NA)

#Collapse entries in the SCCWRP data which contain the model variables of interest
SCCWRP <- SCCWRP[,colnames(SCCWRP) %in% c("stationname","latitude","longitude","from_fire_start","from_fire_end",var_chem)]
SCCWRP <- SCCWRP %>%
  group_by(stationname,latitude,longitude,from_fire_start,from_fire_end) %>%
  dplyr::summarise(
    `Temp (degree_C)` = mean(`Temp (degree_C)`, na.rm = TRUE),
    `Phosphate (uM)` = mean(`Phosphate (uM)`, na.rm = TRUE),
    `Silicate (uM)` = mean(`Silicate (uM)`, na.rm = TRUE),
    `Nitrite (uM)` = mean(`Nitrite (uM)`, na.rm = TRUE),
    `Nitrate (uM)` = mean(`Nitrate (uM)`, na.rm = TRUE),
    `Ammonium (uM)` = mean(`Ammonium (uM)`, na.rm = TRUE),
    .groups = "drop"
  )

for(i in 1:length(var_bio)){
  #Select all physicochemical variables and one biological variable at a time for modelling.
  SCCOOS_filtered <- SCCOOS[,c(var_chem,var_bio[i])]
  #Remove entries with missing data
  SCCOOS_filtered <- SCCOOS_filtered[complete.cases(SCCOOS_filtered),]
  
  #Impute missing values in SCCWRP data using median values derived from SCCOOS data
  for (v in names(SCCWRP[,colnames(SCCWRP) %in% var_chem])) {
    SCCWRP[,colnames(SCCWRP) %in% var_chem][[v]][is.na(SCCWRP[,colnames(SCCWRP) %in% var_chem][[v]])] <-
      median(SCCOOS_filtered[,colnames(SCCOOS_filtered) %in% var_chem][[v]], na.rm = TRUE)
  }
  SCCWRP <- as.data.frame(SCCWRP)
  
  #Run a random forest model over this filtered data.
  #Inputs are chemical data, output is the selected biological variable
  rf_full <- suppressWarnings(tuneRF(x=SCCOOS_filtered[,colnames(SCCOOS_filtered) %in% var_chem],y=SCCOOS_filtered[,colnames(SCCOOS_filtered) %in% var_bio],stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Calculate the correlation between actual and predicted algal species variable to evaluate model accuracy.
  #Calculate the predicted algal concentration
  SCCWRP[,var_bio[i]] <- predict(rf_full,SCCWRP[,colnames(SCCWRP) %in% var_chem])
  
}

#Summarize how much predicted algal concentrations vary by location and time
model_summary <- c()
for(i in 1:length(var_bio)){
  tmp <- anova(lm(as.formula(paste("`",var_bio[i],"`", "~stationname*from_fire_start",sep="")), data=SCCWRP))
  tmp$factor <- rownames(tmp)
  tmp$algae <- var_bio[i]
  model_summary[[i]] <- tmp
}
model_summary <- rbindlist(model_summary,fill=T)

#Merge modeled algal concentrations back into SCCWRP data with metal concentrations
SCCWRP_model <- dplyr::left_join(SCCWRP,SCCWRP_input)

#Define metal variables
var_metal <- c("Zinc","Vanadium","Tin","Thallium","Silver","Selenium","Nickel","Molybendum","Mercury","Manganese","Magnesium","Lead","Iron","Hexavalent chromium","Copper","Chromium","Cadmium","Beryllium","Arsenic","Antimony","Aluminum")

#Test to see how strongly modeled algal concentrations are correlated with metal concentrations
k <- 1
metal_tests <- c()
for(i in 1:length(var_metal)){
  for(j in 1:length(var_bio)){
    #Get a subset of data focusing on a particular metal.
    SCCWRP_subset <- SCCWRP_model[SCCWRP_model$analytename==var_metal[i],]
    #Initialize an empty temporary matrix to store correlations
    tmp <- data.frame(matrix(nrow=1,ncol=4))
    colnames(tmp) <- c("metal","algae","rho","p")
    #Test correlations between metals and algal bloom concentrations
    tmp$algae <- var_bio[j]
    tmp$metal <- var_metal[i]
    #Check if there's enough data to test for a correlation
    if(nrow(SCCWRP_subset)>1){
      tmp$rho <- cor.test(SCCWRP_model[SCCWRP_model$analytename==var_metal[i],"result"],SCCWRP_model[SCCWRP_model$analytename==var_metal[i],var_bio[j]],method="spearman")$estimate
      tmp$p <- cor.test(SCCWRP_model[SCCWRP_model$analytename==var_metal[i],"result"],SCCWRP_model[SCCWRP_model$analytename==var_metal[i],var_bio[j]],method="spearman")$p.value
    } else{
      tmp$rho <- NA
      tmp$p <- NA
    }
    #Store results
    metal_tests[[k]] <- tmp
    k <- k+1
  }
}
metal_tests <- rbindlist(metal_tests,fill=T)

#Visualize correlations between metal and algal concentrations in a grid
#Clean up algal labels
metal_tests$algae <- gsub("_"," ",gsub(" \\(cells/L\\)","",metal_tests$algae))
ggplot(metal_tests, aes(x = algae, y = metal)) +
  ## Non-significant tiles (p > 0.05)
  geom_tile(
    data = subset(metal_tests, p > 0.05),
    fill = "gray80",
    color = "white"
  ) +
  ## Significant tiles (p ≤ 0.05), colored by rho
  geom_tile(
    data = subset(metal_tests, p <= 0.05),
    aes(fill = rho),
    color = "white"
  ) +
  scale_fill_viridis_c(name = "rho") +
  labs(x = "algal concentrations (cells/L)", y = "metal concentrations (uM)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    text = element_text(size = 12),
    panel.grid = element_blank()
  )

#Rearrange modeled HAB algal data for plotting over time
SCCWRP_HAB_time <- SCCWRP %>%
  pivot_longer(
    cols = all_of(var_bio),
    names_to = "Variable",
    values_to = "Value"
  )
#Clean up algal labels
SCCWRP_HAB_time$algae <- gsub("_"," ",gsub(" \\(cells/L\\)","",SCCWRP_HAB_time$Variable))
#Plot modeled algal concentrations over time
ggplot(SCCWRP_HAB_time, aes(x = from_fire_start, y = Value)) +
  geom_line() + 
  facet_wrap(~ algae, scales = "free_y") +
  theme_bw() +
  labs(x = "Days from fire start", y = "cells / L")

#Map project data sets
SCCOOS_tmp <- SCCOOS[,c("longitude (degrees_east)","latitude (degrees_north)")]
SCCOOS_tmp <- SCCOOS_tmp[!duplicated(SCCOOS_tmp),]
SCCWRP_tmp <- SCCWRP[,c("longitude","latitude")]
SCCWRP_tmp <- SCCWRP_tmp[!duplicated(SCCWRP_tmp),]
leaflet() |>
  addProviderTiles(providers$CartoDB.Positron) |>
  addCircleMarkers(
    lng=SCCOOS_tmp$`longitude (degrees_east)`,
    lat=SCCOOS_tmp$`latitude (degrees_north)`,
    radius = 6,
    fillColor = "orange",
    fillOpacity = 0.5,
    stroke = TRUE,
    color = "black", # Border color
    weight = 1
  ) |>
  addCircleMarkers(
    lng=SCCWRP_tmp$longitude,
    lat=SCCWRP_tmp$latitude,
    radius = 3,
    fillColor = "blue",
    fillOpacity = 0.5,
    stroke = TRUE,
    color = "black", # Border color
    weight = 1
  )

#
#Clean up algal labels
importance_total$algae <- gsub("_"," ",gsub(" \\(cells/L\\)","",importance_total$algae))
ggplot(importance_total, 
       aes(x = algae, 
           y = VariableName, 
           fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +   # nice continuous color scale
  theme_minimal() +
  labs(x = "cells / L",
       y = "Variable",
       fill = "Rank\nimportance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
