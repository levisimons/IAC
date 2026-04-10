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

#Create 10km sampling buffer around SCCOOS stations to capture nearby microplastics data
tmp <- SCCOOS[,c("Location_Code","latitude (degrees_north)","longitude (degrees_east)")]
tmp <- tmp[!duplicated(tmp),]
tmp <- st_as_sf(x = tmp, coords = c("longitude (degrees_east)", "latitude (degrees_north)"), crs = "+proj=longlat +datum=WGS84")
SCCOOS_buffer <- tmp |>
  st_transform(3857) |>      
  st_buffer(10000) |>
  st_transform(4326)
SCCOOS_buffer <- SCCOOS_buffer |>
  st_make_valid()

#Create a sampling day variable
SCCOOS$day <- as.numeric(as.Date(SCCOOS$`time (UTC)`)-min(as.Date(SCCOOS$`time (UTC)`)))

#Read in microplastics data
#Obtained from: https://github.com/levisimons/IAC/blob/main/Microplastics/Microplastics.csv
microplastics_input <- fread(input="Microplastics.csv",sep=",")
microplastics_input[microplastics_input==""] <- NA

#Read in microplastics metadata
#Obtained from: https://github.com/levisimons/IAC/blob/main/Microplastics/Microplastics_metadata.csv
microplastics_metadata <- fread(input="Microplastics_metadata.csv",sep=",")
#Force longitude to being numeric
microplastics_metadata$Longitude <- as.numeric(gsub("[^0-9.-]", "", microplastics_metadata$Longitude))

#Create a unique location column
microplastics_metadata <- microplastics_metadata %>%
  dplyr::group_by(Latitude, Longitude) %>%
  dplyr::mutate(Location = cur_group_id())

#Combine microplastics metadata with microplastics data
microplastics <- dplyr::left_join(microplastics_input,microplastics_metadata,relationship="many-to-many")

#Retain data with known microplastic concentrations
microplastics <- microplastics[microplastics$`Plastic or Not`=="plastic" & !is.na(microplastics$`MP Concentration`),c("Location","MP Concentration","Sample Collection Date","Longitude","Latitude")]
microplastics <- microplastics[!duplicated(microplastics),]

#Create numeric sampling day variable
microplastics$day <- as.numeric(as.Date(microplastics$`Sample Collection Date`,format="%m/%d/%y")-min(as.Date(SCCOOS$`time (UTC)`)))

#Join microplastics data with nearby SCCOOS sampling locations
microplastics <- st_as_sf(microplastics,coords=c("Longitude","Latitude"),crs = 4326)
microplastics <- st_join(microplastics,SCCOOS_buffer)

#Create a LOESS function between microplastic concentration and sampling day
microplastic_loess <- loess(microplastics$`MP Concentration`~microplastics$day)
#Predict microplastic concentrations over time
microplastics_predicted <- as.data.frame(predict(microplastic_loess,c(min(microplastics$day):max(microplastics$day))))
microplastics_predicted$day <- c(min(microplastics$day):max(microplastics$day))
colnames(microplastics_predicted) <- c("predicted_concentration","day")

#Filter SCCOOS data to only contain locations with interpolated microplastics concentrations
SCCOOS_subset <- SCCOOS[SCCOOS$Location_Code %in% na.omit(unique(microplastics$Location_Code)),]
#Add in interpolated microplastics concentrations
SCCOOS_subset <- dplyr::left_join(SCCOOS_subset,microplastics_predicted)
#Only retain present and positive interpolated microplastics concentrations
SCCOOS_subset <- SCCOOS_subset[!is.na(SCCOOS_subset$predicted_concentration) & SCCOOS_subset$predicted_concentration > 0,]

#Define physicochemical variables
var_chem <- c("Temp (degree_C)","Phosphate (uM)","Silicate (uM)","Nitrite (uM)","Nitrate (uM)","Ammonium (uM)","predicted_concentration")

#Define biological variables
var_bio <- c("Akashiwo_sanguinea (cells/L)","Alexandrium_spp (cells/L)","Dinophysis_spp (cells/L)","Lingulodinium_polyedra (cells/L)","Prorocentrum_spp (cells/L)","Pseudo_nitzschia_delicatissima_group (cells/L)","Pseudo_nitzschia_seriata_group (cells/L)","Ceratium_spp (cells/L)","Cochlodinium_spp (cells/L)","Gymnodinium_spp (cells/L)","Other_Diatoms (cells/L)","Other_Dinoflagellates (cells/L)","Total_Phytoplankton (cells/L)")

#Store variable importances for all models
importance_total <- c()
#Store model accuracies for all models
accuracy_total <- c()
for(i in 1:length(var_bio)){
  #Select all physicochemical variables and one biological variable at a time for modelling.
  SCCOOS_filtered <- SCCOOS_subset[,c(var_chem,var_bio[i])]
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
    if(!is.na(cor.test(prediction,SCCOOS_testing[,var_bio[i]],method="spearman")$p.value) & cor.test(prediction,SCCOOS_testing[,var_bio[i]],method="spearman")$p.value < 0.05){
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
importance_total$VariableName <- ifelse(importance_total$VariableName=="predicted_concentration","Estimated microplastics concentration\n(particles per cubic meter",importance_total$VariableName)
accuracy_total <- rbind.fill(accuracy_total,fill=T)

ggplot(importance_total, 
       aes(x = algae, 
           y = VariableName, 
           fill = Importance)) +
  geom_tile() +
  scale_fill_viridis_c(direction = -1) +   # good continuous color scale
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,   # right-align so labels don’t overlap ticks
      vjust = 1
    )
    )+
  labs(x = "Algae",
       y = "Variable Name",
       fill = "Importance")
