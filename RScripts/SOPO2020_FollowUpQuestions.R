#=====================================================================================================
# Script Name: SOPO2020_FollowUpQuestions.R
# Script Author: Erika Anderson
# Script Start Date: 2021-03
# R Version: 3.6.1
#
# State of the Pacific Ocean meeting March 2021
# HSSALMON database limited to same areas (DE, HS, QCSD) in fall (SEP, OCT, NOV)
# CPUE anomalies by species and region using swept volume
# to address questions about target head rope depth of 0 or 15 m for chinook
# JK would like to see CPUE by region for 15 m only for chinook 
# graph chinook caught by region and depth for fun 
#
#=====================================================================================================

# load libraries
library(tidyverse) # load core packages
library(here) # to use relative file names
library(lubridate) # dates
library(RODBC) # MS Access databases
library(patchwork) # combine plots
library(modelr) # models length to weight, residuals
library(broom) # evaluate models
library(oce) # load ctd data
library(patchwork) # combine plots
library(readxl) # read excel files for GSI
library(sf) # spatial manipulation (newer than sp) so works with ggplot2
library(sp) # spatial data manipulation
library(rgdal) # to load shapefiles and rasters


#####################################
# database file path
# change these paths for your own version
#####################################

# connection to local high seas MS Access database
db_hs <- "C:/Users/andersoned/Documents/High Seas Salmon Database/HSSALMON.accdb"

outputFolder <- "Output/2020/Questions/"

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
WHERE (((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.Month)='SEP' Or (BRIDGE_COMPLETE.Month)='OCT' Or (BRIDGE_COMPLETE.Month)='NOV') AND ((CATCH_ALL.SPECIES_CODE)='108J' Or (CATCH_ALL.SPECIES_CODE)='118J' Or (CATCH_ALL.SPECIES_CODE)='112J' Or (CATCH_ALL.SPECIES_CODE)='115J' Or (CATCH_ALL.SPECIES_CODE)='124J'));
                        ")

tows_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.CRUISE, BRIDGE_COMPLETE.STATION_ID, BRIDGE_COMPLETE.REGION, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.Year, BRIDGE_COMPLETE.Month, BRIDGE_COMPLETE.START_LAT, BRIDGE_COMPLETE.START_LONG, BRIDGE_COMPLETE.DISTANCE, BRIDGE_COMPLETE.DUR, BRIDGE_COMPLETE.[SOG-KTS], BRIDGE_COMPLETE.START_BOT_DEPTH, BRIDGE_COMPLETE.END_BOT_DEPTH, BRIDGE_COMPLETE.HEAD_DEPTH, BRIDGE_COMPLETE.NET_OPENING_WIDTH, BRIDGE_COMPLETE.NET_OPENING_HEIGHT
FROM BRIDGE_COMPLETE
WHERE (((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.Month)='SEP' Or (BRIDGE_COMPLETE.Month)='OCT' Or (BRIDGE_COMPLETE.Month)='NOV'));
                        ")

# close database
close(myconn_hs)

#####################################
# wrangle high seas catch data
#####################################
# # save all tows as csv to display in ArcMap to make a raster
# write_csv(tows_hs_orig, file = str_c(outputFolder, "tows_hs_orig.csv"))

# after mapping csv of tows, recommend removing some tows as not really within regions
### might want to change these in database
removeStations <- c("HS9819-D06", "HS9819-D07", "HS9819-D08", "HS9819-D09",
                    "HS9819-D11", "HS9819-D13",
                    "BCSI-201778-QCSD17", "BCSI-201778-QCSD21", "HS201145-IBC19")

cpue_hs_orig <- cpue_hs_orig %>%
  filter(!(STATION_ID %in% removeStations))

tows_hs_orig <- tows_hs_orig %>%
  filter(!(STATION_ID %in% removeStations))

# need to add zero tows for individual salmon species
### zeros missing from CATCH_ALL table in 2020

# make function to add zero tows to each salmon speies
zero_fn <- function(cpue_hs_orig, tows_hs_orig, speciesName) {
  
  cpue_hs_orig %>%
    filter(SPECIES_CODE == speciesName) %>%
    full_join(., tows_hs_orig, by = c("CRUISE", "STATION_ID", "REGION", "REGION_CODE", "Year", "Month", "START_LAT", "START_LONG", "DISTANCE", "DUR", "SOG-KTS", "START_BOT_DEPTH", "END_BOT_DEPTH", "HEAD_DEPTH", "NET_OPENING_WIDTH", "NET_OPENING_HEIGHT")) %>%
    mutate(SPECIES_CODE = speciesName,
           COUNT0 = if_else(is.na(Count), 0, Count)) %>%
    select(CRUISE, STATION_ID, REGION_CODE, Year, Month, START_LAT, START_LONG,
           HEAD_DEPTH, DISTANCE, DUR, `SOG-KTS`, NET_OPENING_WIDTH, NET_OPENING_HEIGHT,
           SPECIES_CODE, Count, COUNT0)
  
}

# cpue df with zero tows
cpue_hs_pk <- zero_fn(cpue_hs_orig, tows_hs_orig, "108J")
cpue_hs_cm <- zero_fn(cpue_hs_orig, tows_hs_orig, "112J")
cpue_hs_co <- zero_fn(cpue_hs_orig, tows_hs_orig, "115J")
cpue_hs_se <- zero_fn(cpue_hs_orig, tows_hs_orig, "118J")
cpue_hs_ck <- zero_fn(cpue_hs_orig, tows_hs_orig, "124J")

# check that tows between species are same 
# # 312 is all with no limit on head rope depth

# bind tpgether
cpue_hs <- rbind(cpue_hs_ck, cpue_hs_cm, cpue_hs_co, cpue_hs_pk, cpue_hs_se)


# check for empty net dimensions and distance values or zeros
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
  group_by(TRIP_YEAR, HEAD_DEPTH, REGION_CODE, STATION_ID) %>%
  count() %>%
  ungroup() %>%
  group_by(TRIP_YEAR, HEAD_DEPTH, REGION_CODE) %>%
  count() %>%
  arrange(TRIP_YEAR)

# create function to calculate anomalies for each salmon species
# add ability to calculate regions seperately
anomFn <- function(cpue, speciesCode, seperateReg, headDepthMax, headDepthMin) {
  
  if (seperateReg == FALSE) {
    cpue_select <- cpue %>%
      filter(SPECIES_CODE == speciesCode) %>%
      filter(HEAD_DEPTH < headDepthMax & HEAD_DEPTH >= headDepthMin) %>%
      group_by(TRIP_YEAR) %>%
      summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE),
                .groups = "drop") 
    
    # calculate mean and standard deviation for time series
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
    
  }
  
  if (seperateReg == TRUE) {
    
    mylist <- list()
    loopVector <- 1:3
    
    for (i in seq_along(loopVector)) {
      
      regVector <- unique(cpue$REGION_CODE)
      regionName <- regVector[i]
      
      cpue_select <- cpue %>%
        filter(SPECIES_CODE == speciesCode & REGION_CODE == regionName) %>%
        filter(HEAD_DEPTH < headDepthMax & HEAD_DEPTH >= headDepthMin) %>%
        group_by(TRIP_YEAR, REGION_CODE) %>%
        summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE),
                  .groups = "drop") 
      
      # calculate mean and standard deviation for time series
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
      
      mylist[[i]] <- cpue_select
      
    }
    
    # make regions into one databframe
    cpue_select <- bind_rows(mylist)
  }
  
  # make year factor for nice graph
  cpue_select$Year_fac <- as.factor(cpue_select$TRIP_YEAR)
  cpue_select$Year_fac <- factor(cpue_select$Year_fac, 
                                 levels = c("1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                            "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                            "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"))
  
  return(cpue_select)
  
}

# create function for different head rope depths of all species
headDepthFn <- function(cpue, anomFn, headDepthMax, headDepthMin) {
# by region for < 1
cpuePK_reg <- anomFn(cpue, "108J", TRUE, headDepthMax, headDepthMin)
cpueCM_reg <- anomFn(cpue, "112J", TRUE, headDepthMax, headDepthMin)
cpueCO_reg <- anomFn(cpue, "115J", TRUE, headDepthMax, headDepthMin)
cpueCK_reg <- anomFn(cpue, "124J", TRUE, headDepthMax, headDepthMin)
cpueSE_reg <- anomFn(cpue, "118J", TRUE, headDepthMax, headDepthMin)

cpue_Reg <- rbind(cpueCK_reg, cpueCO_reg, cpueCM_reg, cpuePK_reg, cpueSE_reg)
}

# by region for 0-1 m
cpue_Reg1 <- headDepthFn(cpue, anomFn, 1, 0)
# by region for < 9 since one chinook caught at 8 m head depth
cpue_Reg9 <- headDepthFn(cpue, anomFn, 9, 0)
# by region for < 16 
cpue_Reg16 <- headDepthFn(cpue, anomFn, 16, 0)
# by region for < 20 
cpue_Reg20 <- headDepthFn(cpue, anomFn, 20, 0)
# try 10-20 m 
cpue_Reg19_20 <- headDepthFn(cpue, anomFn, 20, 10)

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
        #axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"),
        strip.background = element_rect(fill = "white")) 
    
  } else {
    
    regionNameLong <- if (regionName == "DE") {"Dixon Entrance"
    } else if (regionName == "HS") {"Hecate Strait" 
    } else if (regionName == "QCSD") {"Queen Charlotte Sound"
    }
    
    dataOld <- df %>%
      filter(REGION_CODE == regionName) %>%
      mutate(anom = if_else(Year_fac == 2020, NA_real_, anom))
    
    dataCurrent <- df %>%
      filter(REGION_CODE == regionName) %>%
      mutate(anom = if_else(Year_fac != 2020, NA_real_, anom))
    
    ggplot(data = dataOld,
           aes(x = Year_fac, y = anom)) +
      geom_bar(stat = "identity") +
      geom_bar(data = dataCurrent, aes(x = Year_fac, y = anom),
               stat = "identity", fill = "darkred") +
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
        #axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white"),
        plot.subtitle = element_text(hjust = 0.5)) +
      labs(subtitle = regionNameLong)
  }
}

# create function to graph all three regions with different head depths
cpueGraphFN2 <- function(df, cpueGraphFN, dfname) {
  
# create plots
cpueDE <- cpueGraphFN(df, "DE") 
cpueHS <- cpueGraphFN(df, "HS") 
cpueQCSD <- cpueGraphFN(df, "QCSD") 

# combine into one plot
cpueRegions <- cpueDE / cpueHS / cpueQCSD

ggsave(str_c(outputFolder, "CPUE_Regions2020_", dfname, "_", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 9, units = "in")
}

cpueGraphFN2(cpue_Reg1, cpueGraphFN, "1")
cpueGraphFN2(cpue_Reg9, cpueGraphFN, "9")
cpueGraphFN2(cpue_Reg16, cpueGraphFN, "16")
cpueGraphFN2(cpue_Reg20, cpueGraphFN, "20")

# try another way when facets are not all represented
#cpueGraphFN2(cpue_Reg19_20, cpueGraphFN, "10-20")

# create plots
#cpueDE <- cpueGraphFN(cpue_Reg19_20, "DE") # no tows enough to graph
cpueHS <- cpueGraphFN(cpue_Reg19_20, "HS") 
cpueQCSD <- cpueGraphFN(cpue_Reg19_20, "QCSD") 

# combine into one plot
cpueRegions <- cpueHS / cpueQCSD

ggsave(str_c(outputFolder, "CPUE_Regions2020_10-20_", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 6, units = "in")

###########################################
# graph chinook by depth and region
###########################################

depth_ck <- cpue_hs_ck %>%
  group_by(Year, REGION_CODE, HEAD_DEPTH) %>%
  summarise(totalNum = sum(COUNT0),
            .groups = "drop")

ggplot(depth_ck, aes(totalNum, HEAD_DEPTH, color = REGION_CODE)) +
  geom_point() +
  scale_y_reverse() +
  theme_bw() +
  labs(x = "Total Juvenile Chinook Caught",
       y = "Depth (m)",
       color = "Region")

ggsave(str_c(outputFolder, "ChinookCatchByDepth", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 4, units = "in")

###########################################

