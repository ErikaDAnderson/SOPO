#=====================================================================================================
# Script Name: SOPO2020.R
# Script Author: Erika Anderson
# Script Start Date: 2021-02
# R Version: 3.6.1
#
# State of the Pacific Ocean meeting March 2021
# CPUE anomalies using swept volume
# Kriging of CPUE from Queen Charlotte Sound to Dixon Entrance
# Length to weight residuals by region (QCSD, HS, DE)
# Calorimetry results if available
# Genetic stock identification if available
# HSSALMON database limited to same area in fall only
# no IPES database needed since no fall surveys for comparison
# daynight only
#
#=====================================================================================================

# load libraries
library(tidyverse) # load core packages
library(here) # to use relative file names
library(lubridate) # dates
library(RODBC) # MS Access databases
library(patchwork) # combine plots
# library(viridis) # colors graphs
# library(sf) # spatial manipulation (newer than sp) so works with ggplot2
# library(sp) # spatial data manipulation
# library(rgdal) # to load shapefiles and rasters
# library(gstat) # model fit & Krige interpolation
# library(data.table) # bind data frames together from list
# library(raster) # load raster for grid (and predict function, alternative to gstat krige)
# library(modelr) # models length to weight, residuals
# library(rcompanion) # confident intervals
# library(readxl) # read excel files for GSI


#####################################
# database file path
# change these database file paths for your own version
#####################################

# connection to local high seas MS Access database
db_hs <- "C:/Users/andersoned/Documents/High Seas Salmon Database/HSSALMON.accdb"

outputFolder <- str_c("Output/2020/")

#####################################
# load high seas data
#####################################
# estalish connection to high sea Access database
myconn_hs <- odbcConnectAccess2007(db_hs)

# get bridge data from high seas
# regions QCSD, HS, DE (not QCI since HG West)
# fall only (SEP, OCT, NOV)
# headrope depth <=20 m to omit deeper tows
cpue_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.CRUISE, BRIDGE_COMPLETE.STATION_ID, BRIDGE_COMPLETE.REGION, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.Year, BRIDGE_COMPLETE.Month, BRIDGE_COMPLETE.START_LAT, BRIDGE_COMPLETE.START_LONG, BRIDGE_COMPLETE.DISTANCE, BRIDGE_COMPLETE.DUR, BRIDGE_COMPLETE.[SOG-KTS], BRIDGE_COMPLETE.START_BOT_DEPTH, BRIDGE_COMPLETE.END_BOT_DEPTH, BRIDGE_COMPLETE.HEAD_DEPTH, BRIDGE_COMPLETE.NET_OPENING_WIDTH, BRIDGE_COMPLETE.NET_OPENING_HEIGHT, CATCH_ALL.SPECIES_CODE, CATCH_ALL.JUVENILE_YN, CATCH_ALL.Count, CATCH_ALL.COUNT_ESTIMATED_YN
FROM BRIDGE_COMPLETE INNER JOIN CATCH_ALL ON BRIDGE_COMPLETE.STATION_ID = CATCH_ALL.STATION_ID
WHERE (((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.Month)='SEP' Or (BRIDGE_COMPLETE.Month)='OCT' Or (BRIDGE_COMPLETE.Month)='NOV') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20) AND ((CATCH_ALL.SPECIES_CODE)='108J' Or (CATCH_ALL.SPECIES_CODE)='118J' Or (CATCH_ALL.SPECIES_CODE)='112J' Or (CATCH_ALL.SPECIES_CODE)='115J' Or (CATCH_ALL.SPECIES_CODE)='124J'));
                          ")

tows_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.CRUISE, BRIDGE_COMPLETE.STATION_ID, BRIDGE_COMPLETE.REGION, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.Year, BRIDGE_COMPLETE.Month, BRIDGE_COMPLETE.START_LAT, BRIDGE_COMPLETE.START_LONG, BRIDGE_COMPLETE.DISTANCE, BRIDGE_COMPLETE.DUR, BRIDGE_COMPLETE.[SOG-KTS], BRIDGE_COMPLETE.START_BOT_DEPTH, BRIDGE_COMPLETE.END_BOT_DEPTH, BRIDGE_COMPLETE.HEAD_DEPTH, BRIDGE_COMPLETE.NET_OPENING_WIDTH, BRIDGE_COMPLETE.NET_OPENING_HEIGHT
FROM BRIDGE_COMPLETE
WHERE (((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.Month)='SEP' Or (BRIDGE_COMPLETE.Month)='OCT' Or (BRIDGE_COMPLETE.Month)='NOV') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20));
                         ")
# get length weight data from high seas
# regions QCSD, HS, DE
# fall only (SEP, OCT, NOV)
# headrope depth <=20 m to omit deeper tows
lw_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.YEAR, BRIDGE_COMPLETE.MONTH, BRIDGE_COMPLETE.REGION_CODE, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL_ALL.SPECIES_CODE, BIOLOGICAL_ALL.SHIP_LENGTH, BIOLOGICAL_ALL.SHIP_WT
FROM (BIOLOGICAL_JUNCTION INNER JOIN BIOLOGICAL_ALL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL_ALL.FISH_NUMBER) INNER JOIN BRIDGE_COMPLETE ON BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE_COMPLETE.STATION_ID
WHERE (((BRIDGE_COMPLETE.MONTH)='SEP' Or (BRIDGE_COMPLETE.MONTH)='OCT' Or (BRIDGE_COMPLETE.MONTH)='NOV') AND ((BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20) AND ((BIOLOGICAL_ALL.SPECIES_CODE)='108' Or (BIOLOGICAL_ALL.SPECIES_CODE)='112' Or (BIOLOGICAL_ALL.SPECIES_CODE)='115' Or (BIOLOGICAL_ALL.SPECIES_CODE)='118' Or (BIOLOGICAL_ALL.SPECIES_CODE)='124'));
                       ")

# get calorimetry data from high seas
# regions QCSD, HS, DE
# fall only (SEP, OCT, NOV)
# headrope depth <=20 m to omit deeper tows
cal_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.Year, BRIDGE_COMPLETE.MONTH, BIOLOGICAL_JUNCTION.FISH_NUMBER, CALORIMETRY.HEAT_RELEASED_CAL, CALORIMETRY.HEAT_RELEASED_KJ, CALORIMETRY.DUPLICATE, CALORIMETRY.DATA_ISSUE, BRIDGE_COMPLETE.HEAD_DEPTH, BRIDGE_COMPLETE.REGION_CODE
FROM (BIOLOGICAL_JUNCTION INNER JOIN CALORIMETRY ON BIOLOGICAL_JUNCTION.FISH_NUMBER = CALORIMETRY.FISH_NUMBER) INNER JOIN BRIDGE_COMPLETE ON BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE_COMPLETE.STATION_ID
WHERE (((BRIDGE_COMPLETE.MONTH)='SEP' Or (BRIDGE_COMPLETE.MONTH)='OCT' Or (BRIDGE_COMPLETE.MONTH)='NOV') AND ((CALORIMETRY.DATA_ISSUE)='N') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20) AND ((BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD'));
                        ")

# get calorimetry data for hisoric data
# regions QCSD, HS, DE
# fall only (SEP, OCT, NOV)
# headrope depth <=20 m to omit deeper tows
cal_historic_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.YEAR, BRIDGE_COMPLETE.MONTH, BIOLOGICAL_ALL.SPECIES_CODE, BIOLOGICAL_JUNCTION.FISH_NUMBER, PROXIMATE_FISH.ENERGY_BOMB, PROXIMATE_FISH.ENERGY_BOMB_BLIND_DUPL, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.HEAD_DEPTH
FROM ((BIOLOGICAL_JUNCTION INNER JOIN BRIDGE_COMPLETE ON BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE_COMPLETE.STATION_ID) INNER JOIN PROXIMATE_FISH ON BIOLOGICAL_JUNCTION.FISH_NUMBER = PROXIMATE_FISH.FISH_NUMBER) INNER JOIN BIOLOGICAL_ALL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL_ALL.FISH_NUMBER
WHERE (((BRIDGE_COMPLETE.MONTH)='SEP' Or (BRIDGE_COMPLETE.MONTH)='OCT' Or (BRIDGE_COMPLETE.MONTH)='NOV') AND ((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20));
                          ")

# close database
close(myconn_hs)

#####################################
# wrangle high seas data
#####################################
# need to add zero tows for individual salmon species
# zeros missing from CATCH_ALL table in 2020

# make function to add zero tows to each salmon speies
zero_fn <- function(cpue_hs_orig, tows_hs_orig, speciesName) {
  
  cpue_hs_orig %>%
    filter(SPECIES_CODE == speciesName) %>%
    full_join(., tows_hs_orig, by = c("CRUISE", "STATION_ID", "REGION", "REGION_CODE", "Year", "Month", "START_LAT", "START_LONG", "DISTANCE", "DUR", "SOG-KTS", "START_BOT_DEPTH", "END_BOT_DEPTH", "HEAD_DEPTH", "NET_OPENING_WIDTH", "NET_OPENING_HEIGHT")) %>%
    mutate(SPECIES_CODE = speciesName,
           COUNT0 = if_else(is.na(Count), 0, Count)) %>%
    select(CRUISE, STATION_ID, REGION_CODE, Year, Month, START_LAT, START_LONG,
           DISTANCE, DUR, `SOG-KTS`, NET_OPENING_WIDTH, NET_OPENING_HEIGHT,
           SPECIES_CODE, Count, COUNT0)
  
}

# cpue df with zero tows
cpue_hs_pk <- zero_fn(cpue_hs_orig, tows_hs_orig, "108J")
cpue_hs_cm <- zero_fn(cpue_hs_orig, tows_hs_orig, "112J")
cpue_hs_co <- zero_fn(cpue_hs_orig, tows_hs_orig, "115J")
cpue_hs_se <- zero_fn(cpue_hs_orig, tows_hs_orig, "118J")
cpue_hs_ck <- zero_fn(cpue_hs_orig, tows_hs_orig, "124J")

# check that tows between species are same 
# # ck, cm, co have 299 rows but pk = 302 and se = 300 why?
# one query had 15 m instead of 20 m head rope depths - fized
# anyDuplicated(cpue_hs_se$STATION_ID)
# anyDuplicated(cpue_hs_pk$ STATION_ID)
# 
# unique(cpue_hs_pk$STATION_ID[!cpue_hs_pk$STATION_ID %in% cpue_hs_ck$STATION_ID])
# unique(cpue_hs_se$STATION_ID[!cpue_hs_se$STATION_ID %in% cpue_hs_ck$STATION_ID])

# bind tpgether
cpue_hs <- rbind(cpue_hs_ck, cpue_hs_cm, cpue_hs_co, cpue_hs_pk, cpue_hs_se)


# check for empty net dimensions and distance values
missing_hs <- cpue_hs %>%
  filter(is.na(NET_OPENING_WIDTH)| is.na(NET_OPENING_HEIGHT) | is.na(DISTANCE)| 
           NET_OPENING_WIDTH == 0 | NET_OPENING_HEIGHT == 0 | DISTANCE == 0)


#### distance missing from cruise 9640 and 201466 
# zero NET_OPENING_HEIGHT in 201267 -> replace with average of 15 m
# go back and add duration and speed to original query to calculate distance
# calculate missing distance from duration and speed
cpue <- cpue_hs %>%
  mutate(DISTANCE_NM = if_else(is.na(DISTANCE), DUR*`SOG-KTS`, DISTANCE),
         NET_OPENING_HEIGHT = if_else(NET_OPENING_HEIGHT == 0, 15, NET_OPENING_HEIGHT),
         NET_AREA_KM = (NET_OPENING_WIDTH/1000) * (NET_OPENING_HEIGHT/1000),
         DISTANCE_KM = DISTANCE_NM * 1.852,
         SWEPT_VOLUME = NET_AREA_KM*DISTANCE_KM,
      # make CPUE by swept volume and log(cpue + 1)
        CPUE = COUNT0 / SWEPT_VOLUME,
         logCPUE1 = log(CPUE + 1), # natural log by default
        TRIP_YEAR = Year) %>%
  # remove station with no distance, nor speed to calculate it
  # only zero counts at this station anyway
  # can't check dataset as why missing and no end lat/long either
  filter(STATION_ID != "HS201466-T04")

#####################################
# cpue annual anomalies
#####################################

# number of tows per year for High seas
cpueNum_hs_tows <- cpue %>%
  group_by(TRIP_YEAR, STATION_ID) %>%
  count() %>%
  ungroup() %>%
  group_by(TRIP_YEAR) %>%
  count()

# create function to calculate anomalies for each salmon species
anomFn <- function(cpue, speciesCode, seperateReg) {
  
  if (seperateReg == FALSE) {
  cpue_select <- cpue %>%
    filter(SPECIES_CODE == speciesCode) %>%
    group_by(TRIP_YEAR) %>%
    summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE),
              .groups = "drop") 
  }
  
  if (seperateReg == TRUE) {
    cpue_select <- cpue %>%
      filter(SPECIES_CODE == speciesCode) %>%
      group_by(TRIP_YEAR, REGION_CODE) %>%
      summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE),
                .groups = "drop") 
  }
  
  # calculate mean and standard deviation for this time series
  meanCPUE_ts <- mean(cpue_select$meanCPUE, na.rm = TRUE)
  sdCPUE_ts <- sd(cpue_select$meanCPUE, na.rm = TRUE)
  
  # calculate anomalies
  cpue_select <- cpue_select %>%
    mutate(anom = (meanCPUE - meanCPUE_ts)/sdCPUE_ts,
           speciesCol = case_when(
             speciesCode == "108J" ~ "Pink",
             speciesCode == "112J" ~ "Chum",
             speciesCode == "115J" ~ "Coho",
             speciesCode == "118J" ~ "Sockeye",
             speciesCode == "124J" ~ "Chinook"))
  
  # make year factor for nice graph
  cpue_select$Year_fac <- as.factor(cpue_select$TRIP_YEAR)
  cpue_select$Year_fac <- factor(cpue_select$Year_fac, 
                                 levels = c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                            "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                            "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"))
  
  return(cpue_select)
  
}

# apply function to species overall regions
cpuePK_df <- anomFn(cpue, "108J", FALSE)
cpueCM_df <- anomFn(cpue, "112J", FALSE)
cpueCO_df <- anomFn(cpue, "115J", FALSE)
cpueCK_df <- anomFn(cpue, "124J", FALSE)
cpueSE_df <- anomFn(cpue, "118J", FALSE)

# by region
cpuePK_reg <- anomFn(cpue, "108J", TRUE)
cpueCM_reg <- anomFn(cpue, "112J", TRUE)
cpueCO_reg <- anomFn(cpue, "115J", TRUE)
cpueCK_reg <- anomFn(cpue, "124J", TRUE)
cpueSE_reg <- anomFn(cpue, "118J", TRUE)

# bind together
cpue_noReg <- rbind(cpueCK_df, cpueCO_df, cpueCM_df, cpuePK_df, cpueSE_df)
cpue_Reg <- rbind(cpueCK_reg, cpueCO_reg, cpueCM_reg, cpuePK_reg, cpueSE_reg)

# graph function
cpueGraphFN <- function(df, regionName) {
  
  if (is.na(regionName) == TRUE) {
    
    ggplot(data = df, aes(x = Year_fac, y = anom)) +
             geom_bar(stat = "identity", fill = "darkred") +
             theme_bw() +
             geom_hline(yintercept = 0) +
             #geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
             facet_wrap(~speciesCol) +
             labs(x = "Ocean Sample Year",
                  y = "ln(CPUE + 1) Anomalies") +
             theme(title = element_text(face = "bold", size = 14)) +
             #geom_vline(xintercept = c(2, 21, 23), linetype = "dotted") +
             #ylim(-2, 2) +
             scale_x_discrete(drop = FALSE,
                              breaks = c(2000, 2005, 2010, 2015, 2020)) +
             theme(
               axis.title = element_text(face = "bold", size = 14),
               axis.text = element_text(size = 12),
               strip.text = element_text(size = 14),
               panel.grid.minor.y = element_blank(), 
               panel.background = element_rect(fill = "white",colour = "black"),
               strip.background = element_rect(fill = "white")) 
           
  } else {
    
    regionNameLong <- if (regionName == "DE") {"Dixon Entrance"
    } else if (regionName == "HS") {"Hecate Strait" 
    } else if (regionName == "QCSD") {"Queen Charlotte Sound"
          }
    
    ggplot(data = filter(cpue_Reg, REGION_CODE == regionName),
           aes(x = Year_fac, y = anom)) +
      geom_bar(stat = "identity", fill = "darkred") +
      theme_bw() +
      geom_hline(yintercept = 0) +
      #geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
      facet_wrap(~speciesCol) +
      labs(x = "Ocean Sample Year",
           y = "ln(CPUE + 1) Anomalies") +
      theme(title = element_text(face = "bold", size = 14)) +
      #geom_vline(xintercept = c(2, 21, 23), linetype = "dotted") +
      #ylim(-2, 2) +
      scale_x_discrete(drop = FALSE,
                       breaks = c(2000, 2005, 2010, 2015, 2020)) +
      theme(
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white")) +
      labs(subtitle = regionNameLong)
  }
}

# apply function overall and to each region
cpueAll <- cpueGraphFN(cpue_noReg, NA) 
cpueAll
ggsave(str_c(outputFolder, "CPUE_NorthCoast2020", str_replace_all(Sys.Date(), "-", ""),".png"))

cpueDE <- cpueGraphFN(cpue_Reg, "DE") 
cpueDE
cpueHS <- cpueGraphFN(cpue_noReg, "HS") 
cpueHS
cpueQCSD <- cpueGraphFN(cpue_noReg, "QCSD") 
cpueQCSD

# combine into one plot
cpueRegions <- cpueDE / cpueHS / cpueQCSD

ggsave(str_c(outputFolder, "CPUE_Regions2020", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 9, units = "in")

#####################################
# use Kriging to see spatial distribution
#####################################
# folder name for Kriging outputs
OutputKriging <- str_c(outputFolder, "Kriging", Sys.Date())

# create directory for plots
dir.create(OutputKriging)

# save all tows as csv to display in ArcMap to make a raster
write_csv(tows_hs_orig, file = str_c(outputFolder, "tows_hs_orig.csv"))

# need grid for Kriging
saFilename <- here::here("Input/Spatial/ipes_wgs84.tif")

# load raster for study area to crop interpolation grid
sa <- raster(saFilename)
class(sa)

# convert to spatial point data frame
sa <- rasterToPoints(sa, spatial = TRUE)
class(sa)

# convert to pixel data frame
gridded(sa) = TRUE
class(sa)

# check grid
plot(sa)