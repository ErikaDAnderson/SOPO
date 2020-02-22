#=====================================================================================================
# Script Name: SOPO2019.R
# Script Author: Erika Anderson
# Script Date: 2020-02
# R Version: 3.6.1
#
# State of the Pacific Ocean meeting March 2020
# CPUE anomalies using swept volume
# Kriging of CPUE over IPES area in June/July
# Length to weight residuals
# Calorimetry results
# Genetic stock identification
# from IPES database
# from high seas salmon database limited spatially
# see email saved in Input for more details about previous years' info
#
#=====================================================================================================

# Load packages
library(tidyverse) # data wrangling
library(magrittr) # save as itself
library(RODBC) # MS Access databases
library(cowplot) # combine plots
library(lubridate) # dates
library(modelr) # models
library(viridis) # colors graphs

#####################################
# estalish connection to high sea Access database
db_hs <- "C:/Users/andersoned/Documents/BCSI/High Seas Salmon Database/HSSALMON.accdb"
myconn_hs <- odbcConnectAccess2007(db_hs)

# get bridge data from high seas
bridge_hs_orig <- sqlQuery(myconn_hs, "SELECT STATION_INFO.CRUISE, STATION_INFO.STATION_ID, STATION_INFO.REGION, STATION_INFO.REGION_CODE, STATION_INFO.SYNOPTIC_STATION, BRIDGE.Year, BRIDGE.Month, BRIDGE.Day, BRIDGE.START_LAT, BRIDGE.START_LONG, BRIDGE.DISTANCE, BRIDGE.START_BOT_DEPTH, BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, BRIDGE.HEAD_DEPTH, BRIDGE.PK_JUV, BRIDGE.CM_JUV, BRIDGE.SE_JUV, BRIDGE.CO_JUV, BRIDGE.CK_JUV
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290))
UNION
SELECT STATION_INFO.CRUISE, STATION_INFO.STATION_ID, STATION_INFO.REGION, STATION_INFO.REGION_CODE, STATION_INFO.SYNOPTIC_STATION, BRIDGE.Year, BRIDGE.Month, BRIDGE.Day, BRIDGE.START_LAT, BRIDGE.START_LONG, BRIDGE.DISTANCE, BRIDGE.START_BOT_DEPTH, BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, BRIDGE.HEAD_DEPTH, BRIDGE.PK_JUV, BRIDGE.CM_JUV, BRIDGE.SE_JUV, BRIDGE.CO_JUV, BRIDGE.CK_JUV
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID
WHERE (((STATION_INFO.REGION_CODE)='QCST') AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290));
")

# get bio data from high seas
bio_hs_orig <- sqlQuery(myconn_hs, "SELECT STATION_INFO.REGION_CODE, BRIDGE.STATION_ID, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL.SHIP_FL, BIOLOGICAL.SHIP_WT
FROM (STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID) INNER JOIN (BIOLOGICAL_JUNCTION INNER JOIN BIOLOGICAL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON (BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((STATION_INFO.REGION_CODE)='QCST') AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290))
UNION
SELECT STATION_INFO.REGION_CODE, BRIDGE.STATION_ID, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL.SHIP_FL, BIOLOGICAL.SHIP_WT
FROM (STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID) INNER JOIN (BIOLOGICAL_JUNCTION INNER JOIN BIOLOGICAL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON (BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290));
                        ")

# get calorimetry data from high seas
cal_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.STATION_ID, BIOLOGICAL_JUNCTION.FISH_NUMBER, CALORIMETRY.HEAT_RELEASED_CAL, CALORIMETRY.HEAT_RELEASED_KJ, CALORIMETRY.DUPLICATE, CALORIMETRY.DATA_ISSUE
FROM (BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN CALORIMETRY ON BIOLOGICAL_JUNCTION.FISH_NUMBER = CALORIMETRY.FISH_NUMBER
WHERE (((CALORIMETRY.DATA_ISSUE)='N'));
")

# get GSI data from high seas

               

# close database
close(myconn_hs)

#####################################
# wrangle high seas data

#####################################
# get IPES data

# CPUE data
# load as csv file since view built on views
cpue_ipes_orig <- read_csv("Input/2019/JB_VIEW_IPES_CPUE.csv")

# salmon only for CPUE
cpue_ipes <- cpue_ipes_orig %>%
  filter(SPECIES_CODE %in% c(108, 112, 115, 118, 124))

# estalish connection to IPES Access database
db_ipes <- "C:/Users/andersoned/Documents/GitHub/IPES_Report/Input/2019/IPES_TrawlDB_v19.07f_2017_18_19.mdb"
myconn_ipes <- odbcConnectAccess2007(db_ipes)

cal_ipes_orig <- sqlQuery(myconn_hs, "SELECT CALORIMETRY_IPES.FISH_NUMBER, CALORIMETRY_IPES.HEAT_RELEASED_CAL, CALORIMETRY_IPES.HEAT_RELEASED_KJ, CALORIMETRY_IPES.DUPLICATE, CALORIMETRY_IPES.DATA_ISSUE
FROM CALORIMETRY_IPES
WHERE (((CALORIMETRY_IPES.DATA_ISSUE)='N'));
")


# close database
close(myconn_ipes)

#####################################

# get data from IPES
# wrangle IPES data


#####################################
# CPUE anomalies 
# species specific changes in abndance

#####################################
# Kriging to see spatial distribution

#####################################
# Length to weight residuals as anomalies over time series
# Condition of salmon

#####################################
# Calorimetry results
# Change this presented as anomalies?
# or just limited years of data

#####################################
# GSI results for 2019

#####################################
# display other species 
# counts or biomass?


#####################################