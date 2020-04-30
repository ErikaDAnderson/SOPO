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
# from IPES database since 2017
# from high seas salmon database limited to same area in June/July only
# see email saved in Input for more details about previous year
# daynight definition needed to be changed slightly since two tows started at 9:55 pm
#
#=====================================================================================================

# load libraries
library(tidyverse) # load core packages
library(here) # to use relative file names
library(magrittr) # save as same name
library(lubridate) # dates
library(RODBC) # MS Access databases
library(egg) # combine plots
library(viridis) # colors graphs
library(sf) # spatial manipulation (newer than sp) so works with ggplot2
library(sp) # spatial data manipulation
library(rgdal) # to load shapefiles and rasters
library(gstat) # model fit & Krige interpolation
library(data.table) # bind data frames together from list
library(raster) # load raster for grid (and predict function, alternative to gstat krige)
library(modelr) # models length to weight, residuals
library(rcompanion) # confident intervals
library(readxl) # read excel files for GSI

#####################################
# load IPES data
#####################################

# CPUE data
# load as csv file since view built on views
# use for swept volume and join to catch 
# use BRIDGE_FIELD_ID because different database version

# March 2020 adjustments to view based on target depth averages instead of overall averages
# from IPES_TrawlDB_v20.02b database 
volume_ipes_orig <- read_csv("Input/2019/JB_VIEW_IPES_CPUE_BRIDGE_LOG_FIELD_ID.csv")

# estalish connection to IPES Access database
# used IPES Report Version since has extra GSI table that I included
db_ipes <- "C:/Users/andersoned/Documents/GitHub/IPES_Report/Input/2019/IPES_TrawlDB_v19.07f_2017_18_19.mdb"
myconn_ipes <- odbcConnectAccess2007(db_ipes)

tows_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, BRIDGE_LOG.EVENT_TYPE, 
BRIDGE_LOG.START_LATITUDE, BRIDGE_LOG.START_LONGITUDE, 
IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, 
BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, BRIDGE_LOG.USABILITY_CODE
FROM TRIP LEFT JOIN BRIDGE_LOG ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 
Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') 
AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((BRIDGE_LOG.USABILITY_CODE)<>5));
                          ")
# limited to daylight, usable tows within IPES area
cpue_ipes_orig <- sqlQuery(myconn_ipes, "SELECT BRIDGE_LOG.BRIDGE_LOG_ID, BRIDGE_LOG.BRIDGE_LOG_FIELD_ID, 
TRIP.TRIP_NAME, BRIDGE_LOG.EVENT_NUMBER, BRIDGE_LOG.EVENT_TYPE, BRIDGE_LOG.START_LATITUDE, 
BRIDGE_LOG.START_LONGITUDE, BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, CATCH.SPECIES_CODE, 
CATCH.JUVENILE_CATCH_COUNT, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, BRIDGE_LOG.USABILITY_CODE
FROM TRIP INNER JOIN (BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) 
ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND 
((BRIDGE_LOG.USABILITY_CODE)<>5));
       ")

# load bio data for length weight 
lw_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, SPECIMEN.UNIVERSAL_FISH_LABEL, 
CATCH.SPECIES_CODE, SPECIMEN.LENGTH, SPECIMEN.WEIGHT
FROM TRIP LEFT JOIN ((BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) 
LEFT JOIN SPECIMEN ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((CATCH.SPECIES_CODE)='108' Or (CATCH.SPECIES_CODE)='112' Or (CATCH.SPECIES_CODE)='115' Or 
(CATCH.SPECIES_CODE)='118' Or (CATCH.SPECIES_CODE)='124') AND ((BRIDGE_LOG.EVENT_TYPE)='midwater tow') 
AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND 
((BRIDGE_LOG.USABILITY_CODE)<>5));
                         ")

# GSI data from IPES for 2019
# IPES survey blocks only
# juvenile defined by length of coho, sockeye, chinook only
# daytime usable tows
gsi_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, SPECIMEN_COLLECTED.COLLECTED_ATTRIBUTE_CODE, 
SPECIMEN_COLLECTED.STORAGE_CONTAINER_SUB_ID, CATCH.SPECIES_CODE, BRIDGE_LOG.BRIDGE_LOG_ID, 
BRIDGE_LOG.BRIDGE_LOG_FIELD_ID, BRIDGE_LOG.TRIP_ID, BRIDGE_LOG.EVENT_DATE, BRIDGE_LOG.BLOCK_DESIGNATION, 
BRIDGE_LOG.USABILITY_CODE, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, SPECIMEN.LENGTH
FROM TRIP LEFT JOIN (BRIDGE_LOG LEFT JOIN (CATCH LEFT JOIN (SPECIMEN LEFT JOIN SPECIMEN_COLLECTED ON 
SPECIMEN.SPECIMEN_ID = SPECIMEN_COLLECTED.SPECIMEN_ID) ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON 
BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((TRIP.TRIP_YEAR)=2019) AND ((SPECIMEN_COLLECTED.COLLECTED_ATTRIBUTE_CODE)=4) AND 
((CATCH.SPECIES_CODE)='115' Or (CATCH.SPECIES_CODE)='118' Or (CATCH.SPECIES_CODE)='124') AND 
((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((BRIDGE_LOG.USABILITY_CODE)<>5) AND 
((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>20.9 Or 
DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND ((SPECIMEN.LENGTH)<350));
                          ")

# close database
close(myconn_ipes)

#####################################
# wrangle IPES data
#####################################

# find distinct tows for IPES data
tows_ipes <-  tows_ipes_orig %>%
  dplyr::select(TRIP_YEAR, BLOCK_DESIGNATION, START_LATITUDE, START_LONGITUDE)
  

# limit swept volume columns
volume_ipes <- volume_ipes_orig %>%
  filter(BLOCK_DESIGNATION > 0) %>%
  filter(DayNight == "Day") %>%
  filter(USABILITY_CODE == 1) %>%
  distinct(BRIDGE_LOG_ID, BRIDGE_LOG_FIELD_ID, TRIP_NAME, EVENT_NUMBER, EVENT_TYPE, BLOCK_DESIGNATION, STRATUM,
         OfficialVolumeSwept_km3)

# select salmon 
# calculate CPUE by swept volume
cpue_ipes_salmon <- cpue_ipes_orig %>%
  filter(SPECIES_CODE %in% c(108, 112, 115, 118, 124)) %>%
  left_join(., volume_ipes, by = c("BRIDGE_LOG_FIELD_ID", "TRIP_NAME", "EVENT_NUMBER", 
                                        "EVENT_TYPE", "BLOCK_DESIGNATION", "STRATUM")) %>%
  mutate(EVENT = str_c(TRIP_NAME, str_pad(EVENT_NUMBER, 3, side = "left", pad = 0), sep = "-"),
         CPUE = JUVENILE_CATCH_COUNT, 
         logCPUE1 = log(CPUE + 1),
         TRIP_YEAR = as.numeric(str_extract(TRIP_NAME, "[0-9]+"))) %>%
  dplyr::select(TRIP_YEAR, BLOCK_DESIGNATION, EVENT, START_LATITUDE, START_LONGITUDE, JUVENILE_CATCH_COUNT, 
         OfficialVolumeSwept_km3, SPECIES_CODE, CPUE, logCPUE1) %>%
  rename(SWEPT_VOLUME = OfficialVolumeSwept_km3)

#####################################
# load high seas data
#####################################
# estalish connection to high sea Access database
db_hs <- "C:/Users/andersoned/Documents/High Seas Salmon Database/HSSALMON.accdb"
myconn_hs <- odbcConnectAccess2007(db_hs)

# get bridge data from high seas
# synoptic stations only
# June and July
# headrope depth <20 m
cpue_hs_orig <- sqlQuery(myconn_hs, "SELECT STATION_INFO.CRUISE, STATION_INFO.STATION_ID, 
STATION_INFO.REGION, STATION_INFO.REGION_CODE, STATION_INFO.SYNOPTIC_STATION, BRIDGE.Year, 
BRIDGE.Month, BRIDGE.START_LAT, BRIDGE.START_LONG, BRIDGE.DISTANCE, BRIDGE.START_BOT_DEPTH, 
BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, BRIDGE.HEAD_DEPTH, 
BRIDGE.PK_JUV, BRIDGE.CM_JUV, BRIDGE.SE_JUV, BRIDGE.CO_JUV, BRIDGE.CK_JUV
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') 
AND ((BRIDGE.HEAD_DEPTH)<20));
                         ")

# get length weight data from high seas
# synoptic stations
# June and July
# headrope depth < 20 m
lw_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.Year, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL.SPECIES, 
BIOLOGICAL.SHIP_FL, BIOLOGICAL.SHIP_WT
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON 
BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN BIOLOGICAL ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON 
(STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) AND (STATION_INFO.STATION_ID = BRIDGE.STATION_ID)
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') AND 
(BRIDGE.HEAD_DEPTH)<20);
                       ")

# get calorimetry data from high seas
cal_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.YEAR, BIOLOGICAL_JUNCTION.FISH_NUMBER, 
CALORIMETRY.HEAT_RELEASED_CAL, CALORIMETRY.HEAT_RELEASED_KJ, CALORIMETRY.DUPLICATE, 
CALORIMETRY.DATA_ISSUE
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON 
BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN CALORIMETRY ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = CALORIMETRY.FISH_NUMBER) ON 
(STATION_INFO.STATION_ID = BRIDGE.STATION_ID) AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((CALORIMETRY.DATA_ISSUE)='N') AND ((STATION_INFO.SYNOPTIC_STATION)=True) AND 
((BRIDGE.MONTH)='JUN' Or (BRIDGE.MONTH)='JUL') AND (BRIDGE.HEAD_DEPTH)<20);
                        ")

# get calorimetry data for IPES from high seas table
# need to limit using bridge info to IPES blocks and usable tows
cal_ipes_orig <- sqlQuery(myconn_hs, "SELECT CALORIMETRY_IPES.FISH_NUMBER, CALORIMETRY_IPES.DUPLICATE, 
CALORIMETRY_IPES.HEAT_RELEASED_CAL, CALORIMETRY_IPES.HEAT_RELEASED_KJ, CALORIMETRY_IPES.DATA_ISSUE, 
CALORIMETRY_IPES.COMMENTS
FROM CALORIMETRY_IPES
WHERE (((CALORIMETRY_IPES.DATA_ISSUE)='N'));
                          ")

# get calorimetry data for hisoric data
# limited to synoptic stations
# June and July
# headrope depth <20 m
cal_historic_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.YEAR, BRIDGE.MONTH, BIOLOGICAL.SPECIES, 
BIOLOGICAL_JUNCTION.FISH_NUMBER, PROXIMATE_FISH.ENERGY_BOMB, PROXIMATE_FISH.ENERGY_BOMB_BLIND_DUPL
FROM STATION_INFO INNER JOIN (((BIOLOGICAL_JUNCTION INNER JOIN PROXIMATE_FISH ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = PROXIMATE_FISH.FISH_NUMBER) INNER JOIN BRIDGE ON 
BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE.STATION_ID) INNER JOIN BIOLOGICAL ON 
BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON (STATION_INFO.STATION_ID = BRIDGE.STATION_ID) 
AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((BRIDGE.MONTH)='JUN' Or (BRIDGE.MONTH)='JUL') AND ((STATION_INFO.SYNOPTIC_STATION)=True) AND 
(BRIDGE.HEAD_DEPTH)<20);
                          ")

# close database
close(myconn_hs)

#####################################
# wrangle high seas data
#####################################

# check for empty net dimensions and distance values
missing_hs <- cpue_hs_orig %>%
  filter(is.na(NET_OPENING_WIDTH)| is.na(NET_OPENING_HEIGHT) | is.na(DISTANCE))

# pull vector of cruises with missing info
missing_hs_vec <- missing_hs %>%
  pull(CRUISE)

#### all missing from cruise 201893 
# no missing distance values

# caluculate net averages for cruises at specific headrope depths
hs_net_avg <- cpue_hs_orig %>%
  group_by(CRUISE, HEAD_DEPTH) %>%
  summarize(AvgNET_OPENING_WIDTH = mean(NET_OPENING_WIDTH, na.rm = TRUE),
            AvgNET_OPENING_HEIGHT = mean(NET_OPENING_HEIGHT, na.rm = TRUE)) %>%
  filter(CRUISE %in% missing_hs_vec)

# use Sea Crest average net width and height from whole survey
# Jackie did this in past with net height = 14 and net width = 33 m
# gear comparison study was used initially 
# same vessel, Captain and net but the chain links and floats were different
# use average of 10.6 for headrope ~ 15 m


# replace missing net values
cpue_hs_net <- cpue_hs_orig %>%
  mutate(NET_OPENING_WIDTH = case_when(
    !(is.na(NET_OPENING_WIDTH)) ~ NET_OPENING_WIDTH,
      is.na(NET_OPENING_WIDTH) ~ 33), # from average for 201893 survey
    NET_OPENING_HEIGHT = case_when(
      !(is.na(NET_OPENING_HEIGHT)) ~ NET_OPENING_HEIGHT,
        is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH <= 7~ 14, # from average for 201893 survey
        is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH > 7 ~ 11), # average from same cruise rounded
    NET_AREA_KM = (NET_OPENING_WIDTH/1000) * (NET_OPENING_HEIGHT/1000),
    DISTANCE_KM = DISTANCE * 1.852,
    SWEPT_VOLUME = NET_AREA_KM*DISTANCE_KM)
  

# create function to rearrange format to match IPES for each salmon
wrangle_fn <- function(df, colName, speciesCode) {
  
df %>%
  mutate(BLOCK_DESIGNATION = NA) %>%
  dplyr::select(Year, BLOCK_DESIGNATION, STATION_ID, START_LAT, START_LONG, colName, SWEPT_VOLUME) %>%
  mutate(SPECIES_CODE = speciesCode) %>%
  rename(TRIP_YEAR = Year,
    EVENT = STATION_ID,
    START_LATITUDE = START_LAT,
    START_LONGITUDE = START_LONG,
    JUVENILE_CATCH_COUNT = colName)
}

# apply to all salmon species  
cpue_hs_pk <- wrangle_fn(cpue_hs_net, "PK_JUV", 108)
cpue_hs_ck <- wrangle_fn(cpue_hs_net, "CK_JUV", 124)
cpue_hs_co <- wrangle_fn(cpue_hs_net, "CO_JUV", 115)
cpue_hs_cm <- wrangle_fn(cpue_hs_net, "CM_JUV", 112)
cpue_hs_se <- wrangle_fn(cpue_hs_net, "SE_JUV", 118)

# bind species back together
cpue_hs <- rbind(cpue_hs_ck, cpue_hs_cm, cpue_hs_co, cpue_hs_pk, cpue_hs_se)

# make CPUE by swept volume and log(cpue + 1)
cpue_hs <- cpue_hs %>%
  mutate(CPUE = JUVENILE_CATCH_COUNT / SWEPT_VOLUME,
         logCPUE1 = log(CPUE + 1)) # natural log by default

#####################################
# cpue annual anomalies
#####################################
# need to add zero tows to IPES for individual salmon species

# make function to add zero tows to each salmon speies

zero_fn <- function(cpue_ipes_salmon, tows_ipes, speciesName) {
  
  cpue_ipes_salmon %>%
    filter(SPECIES_CODE == speciesName) %>%
    full_join(., tows_ipes, by = c("TRIP_YEAR", "BLOCK_DESIGNATION", "START_LATITUDE", "START_LONGITUDE")) %>%
    mutate(SPECIES_CODE = speciesName) %>%
    mutate_if(is.numeric, replace_na, replace = 0) %>%
    dplyr::select(TRIP_YEAR, BLOCK_DESIGNATION, EVENT, START_LATITUDE, START_LONGITUDE, JUVENILE_CATCH_COUNT,
            SWEPT_VOLUME, SPECIES_CODE, CPUE, logCPUE1)
    
}

# cpue df with zero tows
cpue_ipes_pk <- zero_fn(cpue_ipes_salmon, tows_ipes, 108)
cpue_ipes_cm <- zero_fn(cpue_ipes_salmon, tows_ipes, 112)
cpue_ipes_co <- zero_fn(cpue_ipes_salmon, tows_ipes, 115)
cpue_ipes_se <- zero_fn(cpue_ipes_salmon, tows_ipes, 118)
cpue_ipes_ck <- zero_fn(cpue_ipes_salmon, tows_ipes, 124)

# IPES tows check that IPES tows between species are same
# if not find duplicates
# originally ck and co had 177
# the others had 175 like the total tows
# then found 2018 had two night tows that started at 9:55 PM being labelled as day
# changed day night definition to omit them
tows_miss_ck <- cpue_ipes_ck %>%
  group_by(TRIP_YEAR, BLOCK_DESIGNATION, EVENT) %>%
  filter(n() > 1)

tows_miss_co <- cpue_ipes_co %>%
  group_by(TRIP_YEAR, BLOCK_DESIGNATION, EVENT) %>%
  filter(n() > 1)

tows_miss_se <- cpue_ipes_se %>%
  group_by(TRIP_YEAR, BLOCK_DESIGNATION, EVENT) %>%
  filter(n() > 1)

tows_miss_cm <- cpue_ipes_cm %>%
  group_by(TRIP_YEAR, BLOCK_DESIGNATION, EVENT) %>%
  filter(n() > 1)

# bind ipes tpgether
cpue_ipes <- rbind(cpue_ipes_ck, cpue_ipes_cm, cpue_ipes_co, cpue_ipes_pk, cpue_ipes_se)

# bind hs and ipes cpue together
cpue <- rbind(cpue_hs, cpue_ipes)

# number of tows per year for CPUE calculation for High seas
cpueNum_hs <- cpue_hs %>%
  group_by(TRIP_YEAR) %>%
  count() %>%
  ungroup() %>%
  mutate(source = "hs")

# number of tows per year for CPUE calculation for IPES from tows
cpueNum_ipes_tows <- tows_ipes %>%
  group_by(TRIP_YEAR) %>%
  count() %>%
  ungroup() %>%
  mutate(source = "ipes_tows")

# number of tows per year for CPUE calculation for IPES from catch
cpueNum_ipes_catch <- cpue_ipes %>%
  distinct(TRIP_YEAR, BLOCK_DESIGNATION) %>%
  group_by(TRIP_YEAR) %>%
  count() %>%
  ungroup() %>%
  mutate(source = "ipes_catch")

# check that catch and tow data are the same for IPES

# make year with number of tows vector for label
# add 2018 IPES and 2018 high seas together
# use hs for 1998 to 2017
# use IPES for 2017 and 2019
yearTowVec <- c("1998\n(65)", "1999\n(65)", "2000\n(80)", "2001\n(115)", "2002\n(65)", "2003\n(40)",
                "2004\n(75)", "2005\n(70)","2006\n(85)", "2007\n(160)", "2008\n(155)", "2009\n(140)",
                "2010\n(155)", "2011\n(155)", "2012\n(120)","2013\n(100)", "2014\n(35)", "2015\n(165)",
                "2016\n(0)", "2017\n(53)", "2018\n(91)", "2019\n(68)")

# create function to calculate anomalies for each salmon species
anom_fn <- function(df, speciesCode) {
  
  cpue_select <- cpue %>%
    filter(SPECIES_CODE == speciesCode) %>%
    group_by(TRIP_YEAR) %>%
    summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE)) %>%
    ungroup() 
  
  # calculate mean and standard deviation for this time series
  meanCPUE_ts <- mean(cpue_select$meanCPUE, na.rm = TRUE)
  sdCPUE_ts <- sd(cpue_select$meanCPUE, na.rm = TRUE)
  
  # calculate anomalies
  cpue_select <- cpue_select %>%
    mutate(anom = (meanCPUE - meanCPUE_ts)/sdCPUE_ts,
           speciesCol = case_when(
             speciesCode == 108 ~ "Pink",
             speciesCode == 112 ~ "Chum",
             speciesCode == 115 ~ "Coho",
             speciesCode == 118 ~ "Sockeye",
             speciesCode == 124 ~ "Chinook"))
  
  # make year factor for nice graph
  cpue_select$Year_fac <- as.factor(cpue_select$TRIP_YEAR)
  cpue_select$Year_fac <- factor(cpue_select$Year_fac, 
                                 levels = c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                            "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                            "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

  return(cpue_select)
  
}

# apply function to species
cpuePK_df <- anom_fn(cpue, 108)
cpueCM_df <- anom_fn(cpue, 112)
cpueCO_df <- anom_fn(cpue, 115)
cpueCK_df <- anom_fn(cpue, 124)
cpueSE_df <- anom_fn(cpue, 118)

# bind together
cpue_df <- rbind(cpueCK_df, cpueCO_df, cpueCM_df, cpuePK_df, cpueSE_df)

# graph
ggplot(cpue_df, aes(x = Year_fac, y = anom)) +
  geom_bar(stat = "identity", fill = "darkred") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  #geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
  facet_wrap(~speciesCol) +
  labs(x = "Ocean Sample Year",
       y = "ln(CPUE + 1) Anomalies") +
  theme(title = element_text(face = "bold", size = 14)) +
  geom_vline(xintercept = 19, linetype = "dotted") +
  ylim(-2, 2) +
  scale_x_discrete(drop = FALSE,
                   breaks = c(2000, 2005, 2010, 2015)) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(), 
        #panel.grid.major.x = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black")) 

ggsave(str_c("Output/2019/CPUE_AllSpecies.png"))

## rerun with no zero tows to see patterns between surveys without zero inflation
# create function to calculate anomalies for each salmon species
anom_no_zero_fn <- function(df, speciesCode) {
  
  cpue_select <- cpue %>%
    filter(CPUE != 0) %>%
    filter(SPECIES_CODE == speciesCode) %>%
    mutate(lncpue = log(CPUE)) %>% # natural log by default
                                   # no need for plus one without zero tows
    group_by(TRIP_YEAR) %>%
    summarize(meanCPUE = mean(lncpue, na.rm = TRUE)) %>%
    ungroup() 
  
  # calculate mean and standard deviation for this time series
  meanCPUE_ts <- mean(cpue_select$meanCPUE, na.rm = TRUE)
  sdCPUE_ts <- sd(cpue_select$meanCPUE, na.rm = TRUE)
  
  # calculate anomalies
  cpue_select <- cpue_select %>%
    mutate(anom = (meanCPUE - meanCPUE_ts)/sdCPUE_ts,
           speciesCol = case_when(
             speciesCode == 108 ~ "Pink",
             speciesCode == 112 ~ "Chum",
             speciesCode == 115 ~ "Coho",
             speciesCode == 118 ~ "Sockeye",
             speciesCode == 124 ~ "Chinook"))
  
  return(cpue_select)
  
}

# apply function to species
cpuePK_noZero <- anom_no_zero_fn(cpue, 108)
cpueCM_noZero <- anom_no_zero_fn(cpue, 112)
cpueCO_noZero <- anom_no_zero_fn(cpue, 115)
cpueCK_noZero <- anom_no_zero_fn(cpue, 124)
cpueSE_noZero <- anom_no_zero_fn(cpue, 118)

# bind together
cpue_noZero <- rbind(cpuePK_noZero, cpueCM_noZero, cpueCO_noZero, cpueCK_noZero, cpueSE_noZero)


# make year factor for nice graph
cpue_noZero$Year_fac <- as.factor(cpue_noZero$TRIP_YEAR)
cpue_noZero$Year_fac <- factor(cpue_noZero$Year_fac, 
                               levels = c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                          "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                          "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

# graph
ggplot(cpue_noZero, aes(x = Year_fac, y = anom)) +
  geom_bar(stat = "identity", fill = "darkred") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  facet_wrap(~speciesCol) +
  labs(x = "Ocean Sample Year",
       y = "ln(CPUE) Anomalies") +
  theme(title = element_text(face = "bold", size = 14)) +
  geom_vline(xintercept = 19, linetype = "dotted") +
  #ylim(-2, 2) +
  scale_x_discrete(drop = FALSE,
                  breaks = c(2000, 2005, 2010, 2015)) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 14),
        #panel.grid.minor.y = element_blank(), 
        #panel.grid.major.x = element_blank(), 
        #panel.background = element_rect(fill = "white",colour = "black")
        ) 

ggsave(str_c("Output/2019/CPUE_AllSpecies_NoZeros.png"))


### try CPUE with no volume correction
# compared min, max and mean of swept volumes
# they were comparable between hs and ipes
anom_noVolume_fn <- function(df, speciesCode) {
  
  cpue_select <- cpue %>%
    filter(SPECIES_CODE == speciesCode) %>%
    mutate(lncpue1 = log(JUVENILE_CATCH_COUNT + 1)) %>% # natural log by default
    group_by(TRIP_YEAR) %>%
    summarize(meanCPUE = mean(lncpue1, na.rm = TRUE)) %>%
    ungroup() 
  
  # calculate mean and standard deviation for this time series
  meanCPUE_ts <- mean(cpue_select$meanCPUE, na.rm = TRUE)
  sdCPUE_ts <- sd(cpue_select$meanCPUE, na.rm = TRUE)
  
  # calculate anomalies
  cpue_select <- cpue_select %>%
    mutate(anom = (meanCPUE - meanCPUE_ts)/sdCPUE_ts,
           speciesCol = case_when(
             speciesCode == 108 ~ "Pink",
             speciesCode == 112 ~ "Chum",
             speciesCode == 115 ~ "Coho",
             speciesCode == 118 ~ "Sockeye",
             speciesCode == 124 ~ "Chinook"))
  
  return(cpue_select)
  
}

# apply function to species
cpuePK_noVol <- anom_noVolume_fn(cpue, 108)
cpueCM_noVol <- anom_noVolume_fn(cpue, 112)
cpueCO_noVol <- anom_noVolume_fn(cpue, 115)
cpueCK_noVol <- anom_noVolume_fn(cpue, 124)
cpueSE_noVol <- anom_noVolume_fn(cpue, 118)

# bind together
cpue_noVol <- rbind(cpuePK_noVol, cpueCM_noVol, cpueCO_noVol, cpueCK_noVol, cpueSE_noVol)


# make year factor for nice graph
cpue_noVol$Year_fac <- as.factor(cpue_noVol$TRIP_YEAR)
cpue_noVol$Year_fac <- factor(cpue_noVol$Year_fac, 
                               levels = c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                          "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                          "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

# graph
ggplot(cpue_noVol, aes(x = Year_fac, y = anom)) +
  geom_bar(stat = "identity", fill = "darkred") +
  theme_bw() +
  geom_hline(yintercept = 0) +
  facet_wrap(~speciesCol) +
  labs(x = "Ocean Sample Year",
       y = "ln(CPUE + 1) Anomalies") +
  theme(title = element_text(face = "bold", size = 14)) +
  geom_vline(xintercept = 19, linetype = "dotted") +
  #ylim(-2, 2) +
  scale_x_discrete(drop = FALSE,
                   breaks = c(2000, 2005, 2010, 2015)) +
  theme(#axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(size = 14),
    #panel.grid.minor.y = element_blank(), 
    #panel.grid.major.x = element_blank(), 
    #panel.background = element_rect(fill = "white",colour = "black")
  ) 

ggsave(str_c("Output/2019/CPUE_AllSpecies_NoVolume.png"))

#####################################
# use Kriging to see spatial distribution in 2019
#####################################
# folder name for Kriging outputs
OutputFolder <- paste0("Output/2019/Kriging", Sys.Date())

# create directory for plots
dir.create(OutputFolder)

# load IPES study area for grid
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

# projection for grid
sa@proj4string

# read coast shapefile for coast line
coast <- readOGR("Input/Spatial/Land.shp")

# reproject into WGS 84
coast <- spTransform(coast, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "))

# convert to dats frame to plot in ggplot
coastdf <- fortify(coast, region = "OBJECTID")

# limit cpue data to 2019 only
df <- cpue %>%
  filter(TRIP_YEAR == 2019) 

# plot to see data
ggplot(df) +
  geom_point(aes(x = START_LONGITUDE, y = START_LATITUDE, color = logCPUE1)) +
  theme_bw()

# check for duplicates locations since they cause errors in interpolation
newdf_dups <- df %>%
  group_by(SPECIES_CODE, START_LATITUDE, START_LONGITUDE) %>%
  filter(n() > 1)

# create empty list
mylist <- list()

# create vector of salmon species excluding pink since there was no daylight catch in 2019
speciesVec <- c(112, 115, 118, 124)

    for (i in speciesVec) {
      
        newdf <- df %>%
          filter(SPECIES_CODE == i) %>%
          rename(lat = START_LATITUDE,
                 long = START_LONGITUDE,
                 logCPUE = logCPUE1) %>%
          mutate(SPECIES_NAME = case_when(
            SPECIES_CODE == 124 ~ "Chinook",
            SPECIES_CODE == 112 ~ "Chum",
            SPECIES_CODE == 115 ~ "Coho",
            SPECIES_CODE == 118 ~ "Sockeye",
            SPECIES_CODE == 108 ~ "Pink"))
        
        # create name of data frame
        nameDf <- unique(newdf$SPECIES_NAME)
        
        # print to console for troubleshooting
        print(nameDf)
        
        # duplicate locations cause issue with interpolations
        # calculate max CPUE and logCPUE for use in duplicate interpolations
        newdf %<>%
          group_by(lat, long) %>%
          summarise(
            maxCPUE = max(CPUE),
            maxlogCPUE =  max(logCPUE)) %>%
          ungroup() %>%
          rename(CPUE = maxCPUE,
                 logCPUE = maxlogCPUE) 
        
        # convert simple data frame into a spatial data frame object
        coordinates(newdf)= ~ long + lat
        
        # set CRS for data frame as geographic spatial coordiantes
        WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
        proj4string(newdf) <- WGS84

        # causes issues with plotting coast with interpolation grid
        # works in geographic coordinates though
        # # project to BC Albers to match study area grid and make pretty
        # newdf <- spTransform(newdf, "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")

        #######################
        # create variogram
        # could use logCPUE or CPUE
        # logCPUE minimizes impact of extra large tows
        dfVariogram = variogram(logCPUE~1, data = newdf)
        
        # gstat function that calculates nugget, sill, range from data
        dfVariogramModel <- fit.variogram(dfVariogram, vgm(c("Gau")))
        
        # plot variogram
        print(plot(dfVariogram, model = dfVariogramModel))
        
        # description of variogram
        print(summary(dfVariogramModel))
        
        #############################################
        
        # an interpolation function from gstat
        TheSurface <- gstat::krige((logCPUE) ~ 1, newdf, sa, model = dfVariogramModel)
        
        # transform surface to data frame for plotting
        TheSurface <- as.data.frame(TheSurface)
        
        # create title for plot
        plotTitle <- str_c(nameDf, " Distribution")
        
        # print to console
        print(plotTitle)
        
        # plot the individual interpolation surfaces
      plot <- ggplot() +
          geom_path(data = coast, aes(x = long, y = lat, group = group)) +
          geom_tile(data = TheSurface, 
                    aes(x = x, y = y, fill = var1.pred)) + 
            coord_equal() +
        xlim(-129.5, -123.5) +
        ylim(48, 51.5) +
          scale_fill_viridis_c() +
          theme_bw() +
          #theme(legend.position = "none") + # used for presentation
          labs(title = nameDf,
               x = "",
               y = "",
               fill = "ln(CPUE+1)") 
        
        plotName <- paste0(OutputFolder, "/Kriging", i, ".png")
        
        ggsave(plotName, plot)
        
        # add new data frame to list
        mylist[[nameDf]] <- plot
    }

# pull out of list
krig_cm <- mylist[["Chum"]]
krig_co <- mylist[["Coho"]]
krig_se <- mylist[["Sockeye"]]
krig_ck <- mylist[["Chinook"]]

# add Kriging plots together
krig_all <- egg::ggarrange(krig_ck, krig_cm, krig_co, krig_se)
krig_all

### added legend to Krigged map for publication
# changed title = nameDf from title = plotTitle
# copy plot from R window to avoid introducing white margins to plot

########################################
### try Kriging on high seas data to see distribution
#### ten years ago in 2009

# choose year to compare
yearID <- 2005

# limit cpue data to one year only
df_hs <- cpue %>%
  filter(TRIP_YEAR == yearID) 

# plot to see data
ggplot(df_hs) +
  geom_point(aes(x = START_LONGITUDE, y = START_LATITUDE, color = logCPUE1)) +
  theme_bw()

# check for duplicates locations since they cause errors in interpolation
newdf_dups <- df_hs %>%
  group_by(SPECIES_CODE, START_LATITUDE, START_LONGITUDE) %>%
  filter(n() > 1)

# create empty list
mylistHS <- list()

# create vector of salmon species excluding pink since there was no daylight catch in 2019
speciesVec <- c(112, 115, 118, 124)

for (i in speciesVec) {
  
  newdf <- df_hs %>%
    filter(SPECIES_CODE == i) %>%
    rename(lat = START_LATITUDE,
           long = START_LONGITUDE,
           logCPUE = logCPUE1) %>%
    mutate(SPECIES_NAME = case_when(
      SPECIES_CODE == 124 ~ "Chinook",
      SPECIES_CODE == 112 ~ "Chum",
      SPECIES_CODE == 115 ~ "Coho",
      SPECIES_CODE == 118 ~ "Sockeye",
      SPECIES_CODE == 108 ~ "Pink"))
  
  # create name of data frame
  nameDf <- unique(newdf$SPECIES_NAME)
  
  # print to console for troubleshooting
  print(nameDf)
  
  # duplicate locations cause issue with interpolations
  # calculate max CPUE and logCPUE for use in duplicate interpolations
  newdf %<>%
    group_by(lat, long) %>%
    summarise(
      maxCPUE = max(CPUE),
      maxlogCPUE =  max(logCPUE)) %>%
    ungroup() %>%
    rename(CPUE = maxCPUE,
           logCPUE = maxlogCPUE) 
  
  # convert simple data frame into a spatial data frame object
  coordinates(newdf)= ~ long + lat
  
  # set CRS for data frame as geographic spatial coordiantes
  WGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
  proj4string(newdf) <- WGS84
  
  # causes issues with plotting coast with interpolation grid
  # works in geographic coordinates though
  # # project to BC Albers to match study area grid and make pretty
  # newdf <- spTransform(newdf, "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
  
  #######################
  # create variogram
  # could use logCPUE or CPUE
  # logCPUE minimizes impact of extra large tows
  dfVariogram = variogram(logCPUE~1, data = newdf)
  
  # gstat function that calculates nugget, sill, range from data
  dfVariogramModel <- fit.variogram(dfVariogram, vgm(c("Gau")))
  
  # plot variogram
  print(plot(dfVariogram, model = dfVariogramModel))
  
  # description of variogram
  print(summary(dfVariogramModel))
  
  #############################################
  
  # an interpolation function from gstat
  TheSurface <- gstat::krige((logCPUE) ~ 1, newdf, sa, model = dfVariogramModel)
  
  # transform surface to data frame for plotting
  TheSurface <- as.data.frame(TheSurface)
  
  # create title for plot
  plotTitle <- str_c(nameDf, " Distribution ", yearID)
  
  # print to console
  print(plotTitle)
  
  # plot the individual interpolation surfaces
  plot <- ggplot() +
    geom_path(data = coast, aes(x = long, y = lat, group = group)) +
    geom_tile(data = TheSurface, 
              aes(x = x, y = y, fill = var1.pred)) + 
    coord_equal() +
    xlim(-129.5, -123.5) +
    ylim(48, 51.5) +
    scale_fill_viridis_c() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = plotTitle,
         x = "",
         y = "") 
  
  plotName <- paste0(OutputFolder, "/KrigingHS", i, ".png")
  
  ggsave(plotName, plot)
  
  # add new data frame to list
  mylistHS[[nameDf]] <- plot
}

# pull out of list
krig_cmHS <- mylistHS[["Chum"]]
krig_coHS <- mylistHS[["Coho"]]
krig_seHS <- mylistHS[["Sockeye"]]
krig_ckHS <- mylistHS[["Chinook"]]
#krig_pkHS <- mylistHS[["Pink"]] # no pinks

# add Kriging plots together
krig_HS <- egg::ggarrange(krig_ckHS, krig_cmHS, krig_coHS, krig_seHS)
krig_HS


#####################################
# Map IPES and high seas tows
#####################################
# use CPUE to map tows from high seas and IPES

# plot the high seas tows
map_hs <- ggplot() +
  geom_path(data = coast, aes(x = long, y = lat, group = group)) +
  geom_point(data = cpue_hs, 
            aes(x = START_LONGITUDE, y = START_LATITUDE),
            color = "darkred") + 
  coord_equal() +
  xlim(-129.5, -123.5) +
  ylim(48, 51.5) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"),
        title = element_text(face = "bold", size = 14)) +
  labs(title = "High Seas Salmon 1998-2015", # has data to 2018 but confusing in report so changed
       x = "",
       y = "") 

ggsave("Output/2019/HS_Tows.png", map_hs)

# plot the IPES tows facet by year
map_ipes <- ggplot() +
  geom_path(data = coast, aes(x = long, y = lat, group = group)) +
  geom_point(data = cpue_ipes, 
             aes(x = START_LONGITUDE, y = START_LATITUDE, color = factor(TRIP_YEAR))) + 
  coord_equal() +
  xlim(-129.5, -123.5) +
  ylim(48, 51.5) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"),
        title = element_text(face = "bold", size = 14)) +
  labs(title = "IPES 2017-2019",
       x = "",
       y = "",
       color = "Year") 

ggsave("Output/2019/IPES_Tows.png", map_ipes)

map_both <- egg::ggarrange(map_hs, map_ipes, 
                           nrow = 1)
map_both

ggsave("Output/2019/HSvsIPESSurveyMaps.png", map_both)

#####################################
# Length to weight residuals 
#####################################
# present as residuals over time series
# gives estimate of condition of salmon

# wrangle hs data
# remove adults
# use 350 mm as limit for all species
lw_hs <- lw_hs_orig %>%
  filter(SHIP_FL < 350) %>%
  mutate(SPECIES_CODE = case_when(
    SPECIES == "CHUM" ~ 112,
    SPECIES == "SOCKEYE" ~ 118,
    SPECIES == "CHINOOK" ~ 124,
    SPECIES == "COHO" ~ 115,
    SPECIES == "PINK" ~ 108
  )) %>%
  rename(TRIP_YEAR = Year,
         LENGTH = SHIP_FL,
         WEIGHT = SHIP_WT) %>%
  filter(!is.na(SPECIES_CODE)) %>%
  dplyr::select(TRIP_YEAR, FISH_NUMBER, SPECIES_CODE, LENGTH, WEIGHT)

# wrangle ipes data
# remove adults
lw_ipes <- lw_ipes_orig %>%
  filter(LENGTH < 350) %>%
  rename(FISH_NUMBER = UNIVERSAL_FISH_LABEL)

# combine IPES and high seas
lw <- rbind(lw_hs, lw_ipes)

# removed rows with blanks
lw <- lw %>%
  filter(!is.na(WEIGHT)) %>%
  filter(!is.na(LENGTH))

# plot prelim data
ggplot(lw, aes(log10(LENGTH), log10(WEIGHT))) +
  geom_point() +
  geom_smooth(method = "lm")

# fit model
modLW <- lm(log10(WEIGHT) ~ log10(LENGTH) * factor(SPECIES_CODE), data = lw)

# add residuals
residsLW <- lw %>%
  add_residuals(., modLW) %>%
  rename(Year = TRIP_YEAR,
         Residuals = resid)

# expand grid for zero years
yearGrid <- expand_grid(Year = 1998:2019) %>%
  mutate(ZeroCount = 0,
         Year = as.integer(Year))

# create function to count species by year for graph
lwCounts_fn <- function(df, speciesCode, yearGrid) {
  
# select species
 df_select <- df %>%
  filter(SPECIES_CODE == speciesCode)

# calculate number of fish per year
 df_select <- df_select %>%
  group_by(Year) %>%
  count() %>%
  rename(PrelimCount = n) %>%
  ungroup() %>%
  full_join(., yearGrid, by = "Year") %>%
  mutate(FinalCount = if_else(is.na(PrelimCount), as.integer(ZeroCount), PrelimCount)) %>%
  arrange(Year) %>%
  dplyr::select(Year, FinalCount)
 
 return(df_select)

}

# apply function to species
lwCount_pk <- lwCounts_fn(residsLW, 108, yearGrid)
lwCount_cm <- lwCounts_fn(residsLW, 112, yearGrid)
lwCount_co <- lwCounts_fn(residsLW, 115, yearGrid)
lwCount_se <- lwCounts_fn(residsLW, 118, yearGrid)
lwCount_ck <- lwCounts_fn(residsLW, 124, yearGrid)

# make years as factors for graphs
residsLW$Year <- factor(residsLW$Year, 
                        levels = c("1998", "1999", "2000", "2001", "2002","2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

# make function to graph LW as boxplots
LWboxplot_fn <- function(df, speciesCode, speciesName) {
  
  # select species
  df_select <- df %>%
    filter(SPECIES_CODE == speciesCode)
  
  # plot title
  plotTitle <- str_c(speciesName, " Length to Weight Relationship")

# graph
ggplot(data = df_select, aes(x = Year, y = Residuals, group = Year)) + 
  geom_boxplot(fill = "darkred") + 
  theme(panel.grid.minor.y = element_blank(), 
        #panel.grid.major.x = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black")) + 
  scale_y_continuous(expand = c(0, 0), limits = c(-.22, .22)) +
  theme(title = element_text(face = "bold", size = 14)) +
  labs(title = plotTitle,
       Y = "",
       x = "Year") +
  geom_hline(yintercept = 0, color = "black", size = 1) 
}

# pink LW plot
boxplot_pk <- LWboxplot_fn(residsLW, 108, "Pink")
boxplot_pk +
  scale_x_discrete(drop = FALSE,
                   labels = c("1998\n(39)", "1999\n(0)", "2000\n(128)", "2001\n(0)", "2002\n(0)","2003\n(0)", 
                 "2004\n(30)", "2005\n(15)", "2006\n(36)", "2007\n(7)", "2008\n(0)",
                 "2009\n(1)", "2010\n(20)", "2011\n(0)", "2012\n(68)", "2013\n(11)",
                  "2014\n(71)", "2015\n(2)", "2016\n(0)", "2017\n(1)", "2018\n(80)", "2019\n(0)"))

ggsave("Output/2019/LW_Pink.png")

# chum LW plot
boxplot_cm <- LWboxplot_fn(residsLW, 112, "Chum")
boxplot_cm +
  scale_x_discrete(drop = FALSE,
                   labels = c("1998\n(68)", "1999\n(206)", "2000\n(220)", "2001\n(36)", "2002\n(0)","2003\n(60)", 
                              "2004\n(44)", "2005\n(17)", "2006\n(137)", "2007\n(164)", "2008\n(0)",
                              "2009\n(110)", "2010\n(85)", "2011\n(66)", "2012\n(109)", "2013\n(64)",
                              "2014\n(41)", "2015\n(120)", "2016\n(0)", "2017\n(255)", "2018\n(164)", "2019\n(181)"))

ggsave("Output/2019/LW_Chum.png")

# coho LW plot
boxplot_co <- LWboxplot_fn(residsLW, 115, "Coho")
boxplot_co +
  scale_x_discrete(drop = FALSE,
                   labels = c("1998\n(42)", "1999\n(99)", "2000\n(150)", "2001\n(148)", "2002\n(0)","2003\n(17)", 
                              "2004\n(24)", "2005\n(51)", "2006\n(125)", "2007\n(150)", "2008\n(0)",
                              "2009\n(245)", "2010\n(81)", "2011\n(85)", "2012\n(178)", "2013\n(77)",
                              "2014\n(3)", "2015\n(235)", "2016\n(0)", "2017\n(84)", "2018\n(161)", "2019\n(89)"))

ggsave("Output/2019/LW_Coho.png")

# sockeye LW plot
boxplot_se <- LWboxplot_fn(residsLW, 118, "Sockeye")
boxplot_se +
  scale_x_discrete(drop = FALSE,
                   labels = c("1998\n(42)", "1999\n(119)", "2000\n(179)", "2001\n(46)", "2002\n(0)","2003\n(86)", 
                              "2004\n(128)", "2005\n(8)", "2006\n(86)", "2007\n(263)", "2008\n(0)",
                              "2009\n(138)", "2010\n(54)", "2011\n(120)", "2012\n(173)", "2013\n(8)",
                              "2014\n(2)", "2015\n(8)", "2016\n(0)", "2017\n(25)", "2018\n(90)", "2019\n(102)"))

ggsave("Output/2019/LW_Sockeye.png")

# chinook LW plot
boxplot_ck <- LWboxplot_fn(residsLW, 124, "Chinook")
boxplot_ck +
  scale_x_discrete(drop = FALSE,
                   labels = c("1998\n(43)", "1999\n(100)", "2000\n(61)", "2001\n(42)", "2002\n(0)","2003\n(53)", 
                              "2004\n(44)", "2005\n(7)", "2006\n(70)", "2007\n(65)", "2008\n(0)",
                              "2009\n(58)", "2010\n(110)", "2011\n(161)", "2012\n(248)", "2013\n(21)",
                              "2014\n(0)", "2015\n(71)", "2016\n(0)", "2017\n(37)", "2018\n(60)", "2019\n(85)"))

ggsave("Output/2019/LW_Chinook.png")

# graph LW as boxplots facted for presentation
# no tow numbers included
# use species names as facet titles

# add species Name
# remove pink since no pink in 2019
residsLW_all <- residsLW %>%
  filter(SPECIES_CODE != 108) %>%
  mutate(SPECIES_NAME = case_when(
    SPECIES_CODE == 124 ~ "Chinook",
    SPECIES_CODE == 112 ~ "Chum",
    SPECIES_CODE == 115 ~ "Coho",
    SPECIES_CODE == 118 ~ "Sockeye",
    SPECIES_CODE == 108 ~ "Pink"))


  # graph
  ggplot(data = residsLW_all, aes(x = Year, y = Residuals, group = Year)) + 
    geom_boxplot(fill = "darkred") + 
    facet_wrap(~SPECIES_NAME) +
    labs(x = "Ocean Sample Year") +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    theme_bw() +
    theme(#axis.text.x = element_text(angle = 90, hjust = 1),
          panel.grid.minor.y = element_blank(), 
          #panel.grid.major.x = element_blank(), 
          panel.background = element_rect(fill = "white",colour = "black"),
          strip.text = element_text(size = 14),
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 14)) +
    scale_x_discrete(drop = FALSE,
                     breaks = c(1995, 2000, 2005, 2010, 2015))
  
  ggsave("Output/2019/LW_all.png")

# # combine all LW plots
# lw_all <- egg::ggarrange(boxplot_ck, boxplot_cm, boxplot_co, boxplot_se)

#####################################
# Calorimetry results
#####################################
# wrangle IPES data

# look at comments for issues
unique(cal_ipes_orig$COMMENTS)

# join to bio data to limit to usable, daylight tows and juveniles
cal_ipes <- cal_ipes_orig %>%
  # average accross duplicates
  group_by(FISH_NUMBER) %>%
  summarize(HEAT_RELEASED_CAL = mean(HEAT_RELEASED_CAL, na.rm = TRUE),
            HEAT_RELEASED_KJ = mean(HEAT_RELEASED_KJ, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(., residsLW, by = "FISH_NUMBER") %>%
  filter(!is.na(SPECIES_CODE))

# wrangle high seas data
cal_hs <- cal_hs_orig %>%
  # average accross duplicates
  group_by(FISH_NUMBER) %>%
  summarize(HEAT_RELEASED_CAL = mean(HEAT_RELEASED_CAL, na.rm = TRUE),
            HEAT_RELEASED_KJ = mean(HEAT_RELEASED_KJ, na.rm = TRUE)) %>%
  ungroup() %>%
  # add species codes and length and weights
  left_join(., residsLW, by = "FISH_NUMBER") %>%
  filter(!is.na(SPECIES_CODE))

# bind data together
cal <- rbind(cal_hs, cal_ipes)

# confirm all juveniles using length column
max(cal$LENGTH)

# compare residuals to heat released values
ggplot(data = filter(cal, SPECIES_CODE != 108), 
                      aes(HEAT_RELEASED_CAL, Residuals, color = Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~SPECIES_CODE) +
   theme_bw()

# this does not look great for correlation
# don't have time now but look into this later
# perhaps need to untransform the Residuals to have clearer relationship
# perhaps relationship only in sockeye?

# calculate mean calorimetry results by year and species
# include confident intervals
cal_ci <- groupwiseMean(HEAT_RELEASED_CAL ~ SPECIES_CODE + factor(Year),
              data = cal,
              conf = 0.95,
              digits = 3)

# wrangle historic high seas calorimetry data
# reported in megajoules per kilogram
# equal to kilojoules per gram reported in newer calorimetry table
# table mean for duplicate samples
cal_historic <- cal_historic_orig %>%
  mutate(SPECIES_CODE = case_when(
    SPECIES == "CHINOOK" ~ 124,
    SPECIES == "COHO" ~ 115),
    HEAT_RELEASED_KJ = if_else(is.na(ENERGY_BOMB_BLIND_DUPL), 
                               ENERGY_BOMB, (ENERGY_BOMB + ENERGY_BOMB_BLIND_DUPL)/2),
    Year = YEAR) %>%
  dplyr::select(FISH_NUMBER, HEAT_RELEASED_KJ, Year, SPECIES_CODE)


# simplify current calorimetry data
cal_simp <- cal %>%
  dplyr::select(FISH_NUMBER, HEAT_RELEASED_KJ, Year, SPECIES_CODE)
  
# bind historic to current calorimetry data
cal_simp <- rbind(cal_simp, cal_historic)

# make year amd species as factors
cal_simp <- cal_simp %>%
  mutate(YearFac = as.factor(Year),
         SPECIES_NAME = case_when(
           SPECIES_CODE == 124 ~ "Chinook",
           SPECIES_CODE == 112 ~ "Chum",
           SPECIES_CODE == 115 ~ "Coho",
           SPECIES_CODE == 118 ~ "Sockeye",
           SPECIES_CODE == 108 ~ "Pink"))

# graph calorimetry results
ggplot(data = filter(cal_simp, SPECIES_CODE != 108),
       aes(x = YearFac, y = HEAT_RELEASED_KJ)) +
  #geom_boxplot(fill = "darkred") +
  geom_violin(color = "black", fill = "darkred") +
  geom_point(color = "black") +
  #geom_jitter(aes(color = YearFac), color = "darkred") +
   theme_bw() +
  facet_wrap(~SPECIES_NAME) +
  theme(panel.grid.minor.y = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"),
        strip.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 14),
        legend.position = "none") +
  #scale_color_viridis_d() +
  labs(x = "Year",
       y = "Heat Released (Kilojoules/gram)")

# save graph
ggsave("Output/2019/calorimetry.png")

#####################################
# GSI results for 2019
#####################################
# Note that there were no GSI data from high seas database
# the June survey was loaded into IPES with negative blocks (omited here)
# the October survey within high seas is not included

# load chinook data for 2019
gsi_ck_orig <- read_excel("Input/2019/chinookBCSI_2019_2020-03-02.xlsx",
                          sheet = "Individual IDs",
                          skip = 3)

reg_ck_orig <- read_excel("Input/2019/chinookBCSI_2019_2020-03-02.xlsx",
                     sheet = "Regional",
                     skip = 11)

# wrangle region data
reg_ck <- reg_ck_orig %>%
  mutate(Region_Code_1 = as.character(Code)) %>%
  dplyr::select(Region_Code_1, Region1)

# wrangle chinook data
gsi_ck <- gsi_ck_orig %>%
  # remove rows between samples and empty rows
  filter(Fish != "2019") %>%
  filter(Stock...3 != "2019") %>%
  rename(Stock = Stock...3,
         Region_Code_1 = Region...4,
         Prob_1 = `Prob 1`) %>%
  dplyr::select(Fish, Stock, Region_Code_1, Prob_1) %>%
  filter(Prob_1 > 0.5) %>%
  left_join(reg_ck, by = "Region_Code_1") %>%
  mutate(DNA_NUMBER = as.numeric(str_extract(Fish, "[0-9]+$")),
         BATCH_NUMBER = str_extract(Fish, "#[0-9]+"),
         BATCH_NUMBER = str_extract(BATCH_NUMBER, "[0-9]+"),
         BATCH_DNA_NUMBER = str_c(BATCH_NUMBER, DNA_NUMBER, 
                                  sep = "-"))

# limit to 2019 only for SOPO
gsi_ck_2019 <- gsi_ck %>%
  filter(BATCH_NUMBER != 68) %>%
  # join to sampling info to limit to usable, day, IPES area tows
  inner_join(gsi_ipes_orig, by = c("DNA_NUMBER" = "STORAGE_CONTAINER_SUB_ID")) 

# get number of samples for chinook
nCk <- nrow(gsi_ck_2019)

# relevel region for graph
gsi_ck_2019$Region1 <- factor(gsi_ck_2019$Region1,
                               levels = c("Snake-Sp/Su", "Up Col-Su/F", "Up Col-Sp",
                                 "North & Central Oregon", "WCVI", "Puget Sound"))

# graph it
gsiPlot_ck <- ggplot(gsi_ck_2019, 
       aes(x = Region1)) +
  geom_bar(fill = "darkred") +
  labs(x = " ",
       y = " ",
       title = "Chinook GSI in June and July 2019",
       caption = str_c("n = ", nCk)) +
  theme_bw() +
  theme(title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  scale_fill_viridis_d()
gsiPlot_ck

# load coho Snp data
gsi_co_orig1 <- read_excel("Input/2019/PID20190082_BCSI_B74(19)_2019-10-11.xlsx",
                           sheet = "collection_table_ids")

# load coho Snp data
gsi_co_orig2 <- read_excel("Input/2019/PID20190108_BCSI_B74(19)_2_sc239_2019-12-02.xlsx",
                           sheet = "collection_table_ids")

# load regional data
reg_co_orig <- read_excel("Input/2019/PID20190108_BCSI_B74(19)_2_sc239_2019-12-02.xlsx",
                          sheet = "repunits_estimates",
                          skip = 6)

# tidy regional data for Coho
reg_co <- reg_co_orig %>%
  dplyr::select(repunit, CU_NAME, Country)

# add vial column into first data set
gsi_co_orig1 <- gsi_co_orig1 %>%
  mutate(vial = as.numeric(str_extract(indiv, "[0-9]+$"))) %>%
  dplyr::select(indiv, vial, everything())

# bind coho data together for 2019  
gsi_co_orig <- rbind(gsi_co_orig1, gsi_co_orig2)

# wrangle coho data
gsi_co <- gsi_co_orig %>%
  dplyr::select(indiv, vial, repunit.1, collection.1, prob.1) %>%
  rename(Fish = indiv,
         DNA_NUMBER = vial,
         Region = repunit.1,
         Stock = collection.1) %>%
  mutate(BATCH_NUMBER = 74,
         BATCH_DNA_NUMBER = str_c(BATCH_NUMBER, DNA_NUMBER,
                                  sep = "-")) %>%
  left_join(., reg_co, by = c("Region" = "repunit")) %>%
  filter(prob.1 > 0.5) %>%
  # limit to daylight tows from IPES area with juvenile lengths
  inner_join(., gsi_ipes_orig, by = c("DNA_NUMBER" = "STORAGE_CONTAINER_SUB_ID"))

# get number of samples
nCo <- nrow(gsi_co)

# relevel region for graph
gsi_co$Region <- factor(gsi_co$Region,
                              levels = c("CR", "HOOD", "EVI+GStr", "JdF",
                                          "OR", "NOOK", "DOUG", "Nahwitti", 
                                         "CLAY", "SC+SFj"))

# graph it
gsiPlot_co <- ggplot(gsi_co, 
                     aes(x = Region)) +
  geom_bar(fill = "darkred") +
  labs(x = " ",
       y = "Count",
       title = "Coho GSI in June and July 2019",
       caption = str_c("n = ", nCo)) +
  theme_bw() +
  theme(title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  scale_fill_viridis_d()
gsiPlot_co  

# load sockeye gsi data
gsi_se_orig <- read_excel("Input/2019/sockeyeBCSIBatch75(19)_2019-11-28.xlsx",
                          sheet = "Individual IDs",
                          skip = 3)

# load SE region data
reg_se_orig <- read_excel("Input/2019/sockeyeBCSIBatch75(19)_2019-11-28.xlsx",
                          sheet = "Regional",
                          skip = 11)

# wrangle region data
reg_se <- reg_se_orig %>%
  mutate(Region_Code_1 = as.character(Code)) %>%
  dplyr::select(Region_Code_1, Region1)

# wrangle sockeye data
gsi_se <- gsi_se_orig %>%
  # remove rows between samples and empty rows
  filter(Fish != "2019") %>%
  filter(Stock...3 != "2019") %>%
  rename(Stock = Stock...3,
         Region_Code_1 = Region...4,
         Prob_1 = `Prob 1`) %>%
  dplyr::select(Fish, Stock, Region_Code_1, Prob_1) %>%
  filter(Prob_1 > 0.5) %>%
  left_join(reg_se, by = "Region_Code_1") %>%
  mutate(DNA_NUMBER = as.numeric(str_extract(Fish, "[0-9]+$")),
         # looked at data and it is all batch 75
         BATCH_NUMBER = 75,
         BATCH_NUMBER = str_extract(BATCH_NUMBER, "[0-9]+"),
         BATCH_DNA_NUMBER = str_c(BATCH_NUMBER, DNA_NUMBER, 
                                  sep = "-"))

# join to sampling info to limit to usable, day, IPES area tows
gsi_se_2019 <- gsi_se %>%
  inner_join(gsi_ipes_orig, by = c("DNA_NUMBER" = "STORAGE_CONTAINER_SUB_ID")) 

# get number of samples fof chinook
nSe <- nrow(gsi_se_2019)

# relevel region for graph
gsi_se_2019$Region1 <- factor(gsi_se_2019$Region1,
                              levels = c("VI", "Summer(Fr)", 
                                         "Washington", "Early Summer(Fr)"))
# graph it
gsiPlot_se <- ggplot(gsi_se_2019, 
                     aes(x = Region1)) +
  geom_bar(fill = "darkred") +
  labs(x = "Region of Origin",
       y = " ",
       title = "Sockeye GSI in June and July 2019",
       caption = str_c("n = ", nSe)) +
  theme_bw() +
  theme(title = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 14),
        strip.text = element_text(size = 14),
        axis.text = element_text(size = 14)) +
  scale_fill_viridis_d()
gsiPlot_se

# plot gsi together
gsi_all <- egg::ggarrange(gsiPlot_ck, gsiPlot_co, gsiPlot_se)
gsi_all
#####################################
# Other species? 
# counts or biomass?
# ran out of time
# get over view from IPES report

#####################################
