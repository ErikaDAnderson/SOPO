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
# from high seas salmon database limited to same area and June/July August?
# see email saved in Input for more details about previous year
#
#=====================================================================================================

# load libraries
library(tidyverse) # load core packages
library(here) # to use relative file names
library(magrittr) # save as same name
library(lubridate) # dates
library(RODBC) # MS Access databases
library(cowplot) # combine plots
library(viridis) # colors graphs
library(sf) # spatial manipulation (newer than sp) so works with ggplot2
library(sp) # spatial data manipulation
library(rgdal) # to load shapefiles and rasters
library(gstat) # model fit & Krige interpolation
library(data.table) # bind data frames together from list
library(raster) # load raster for grid (and predict function, alternative to gstat krige)
library(modelr) # models length to weight, residuals

#####################################
# load IPES data
#####################################

# CPUE data
# load as csv file since view built on views
# use for swept volume and join to catch
volume_ipes_orig <- read_csv("Input/2019/EA_JB_VIEW_IPES_CPUE.csv")

# estalish connection to IPES Access database
db_ipes <- "C:/Users/andersoned/Documents/GitHub/IPES_Report/Input/2019/IPES_TrawlDB_v19.07f_2017_18_19.mdb"
myconn_ipes <- odbcConnectAccess2007(db_ipes)

tows_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, BRIDGE_LOG.EVENT_TYPE, BRIDGE_LOG.START_LATITUDE, BRIDGE_LOG.START_LONGITUDE, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>21 Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, BRIDGE_LOG.USABILITY_CODE
FROM TRIP LEFT JOIN BRIDGE_LOG ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND ((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>21 Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((BRIDGE_LOG.USABILITY_CODE)<>5));
                          ")
# limited to daylight, usable tows within IPES area
cpue_ipes_orig <- sqlQuery(myconn_ipes, "SELECT BRIDGE_LOG.BRIDGE_LOG_ID, TRIP.TRIP_NAME, BRIDGE_LOG.EVENT_NUMBER, BRIDGE_LOG.EVENT_TYPE, BRIDGE_LOG.START_LATITUDE, BRIDGE_LOG.START_LONGITUDE, BRIDGE_LOG.BLOCK_DESIGNATION, BRIDGE_LOG.STRATUM, CATCH.SPECIES_CODE, CATCH.JUVENILE_CATCH_COUNT, IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>21 Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day') AS DayNight, BRIDGE_LOG.USABILITY_CODE
FROM TRIP INNER JOIN (BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>21 Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND ((BRIDGE_LOG.USABILITY_CODE)<>5));
       ")

# load bio data for length weight 
lw_ipes_orig <- sqlQuery(myconn_ipes, "SELECT TRIP.TRIP_YEAR, SPECIMEN.UNIVERSAL_FISH_LABEL, CATCH.SPECIES_CODE, SPECIMEN.LENGTH, SPECIMEN.WEIGHT
FROM TRIP LEFT JOIN ((BRIDGE_LOG LEFT JOIN CATCH ON BRIDGE_LOG.BRIDGE_LOG_ID = CATCH.BRIDGE_LOG_ID) LEFT JOIN SPECIMEN ON CATCH.CATCH_ID = SPECIMEN.CATCH_ID) ON TRIP.TRIP_ID = BRIDGE_LOG.TRIP_ID
WHERE (((CATCH.SPECIES_CODE)='108' Or (CATCH.SPECIES_CODE)='112' Or (CATCH.SPECIES_CODE)='115' Or (CATCH.SPECIES_CODE)='118' Or (CATCH.SPECIES_CODE)='124') AND ((BRIDGE_LOG.EVENT_TYPE)='midwater tow') AND ((BRIDGE_LOG.BLOCK_DESIGNATION)>0) AND ((IIf(DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])>21 Or DatePart('h',[BRIDGE_LOG].[END_DEPLOYMENT_TIME])<6,'Night','Day'))='Day') AND ((BRIDGE_LOG.USABILITY_CODE)<>5));
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
  distinct(BRIDGE_LOG_ID, TRIP_NAME, EVENT_NUMBER, EVENT_TYPE, BLOCK_DESIGNATION, STRATUM,
         OfficialVolumeSwept_km3)

# select salmon 
# calculate CPUE by swept volume
cpue_ipes_salmon <- cpue_ipes_orig %>%
  filter(SPECIES_CODE %in% c(108, 112, 115, 118, 124)) %>%
  left_join(., volume_ipes, by = c("BRIDGE_LOG_ID", "TRIP_NAME", "EVENT_NUMBER", 
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
db_hs <- "C:/Users/andersoned/Documents/BCSI/High Seas Salmon Database/HSSALMON.accdb"
myconn_hs <- odbcConnectAccess2007(db_hs)

# get bridge data from high seas
cpue_hs_orig <- sqlQuery(myconn_hs, "SELECT STATION_INFO.CRUISE, STATION_INFO.STATION_ID, STATION_INFO.REGION, STATION_INFO.REGION_CODE, STATION_INFO.SYNOPTIC_STATION, BRIDGE.Year, BRIDGE.Month, BRIDGE.START_LAT, BRIDGE.START_LONG, BRIDGE.DISTANCE, BRIDGE.START_BOT_DEPTH, BRIDGE.END_BOT_DEPTH, BRIDGE.NET_OPENING_WIDTH, BRIDGE.NET_OPENING_HEIGHT, BRIDGE.HEAD_DEPTH, BRIDGE.PK_JUV, BRIDGE.CM_JUV, BRIDGE.SE_JUV, BRIDGE.CO_JUV, BRIDGE.CK_JUV
FROM STATION_INFO INNER JOIN BRIDGE ON STATION_INFO.STATION_ID = BRIDGE.STATION_ID
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290) AND ((BRIDGE.HEAD_DEPTH)<22));
                         ")

# get length weight data from high seas
lw_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.Year, BIOLOGICAL_JUNCTION.FISH_NUMBER, BIOLOGICAL.SPECIES, BIOLOGICAL.SHIP_FL, BIOLOGICAL.SHIP_WT
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN BIOLOGICAL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL.FISH_NUMBER) ON (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) AND (STATION_INFO.STATION_ID = BRIDGE.STATION_ID)
WHERE (((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.Month)='JUN' Or (BRIDGE.Month)='JUL') AND ((BRIDGE.HEAD_DEPTH)<22) AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290));
                       ")

# get calorimetry data from high seas
cal_hs_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE.YEAR, BRIDGE.STATION_ID, BIOLOGICAL_JUNCTION.FISH_NUMBER, CALORIMETRY.HEAT_RELEASED_CAL, CALORIMETRY.HEAT_RELEASED_KJ, CALORIMETRY.DUPLICATE, CALORIMETRY.DATA_ISSUE
FROM STATION_INFO INNER JOIN ((BRIDGE INNER JOIN BIOLOGICAL_JUNCTION ON BRIDGE.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) INNER JOIN CALORIMETRY ON BIOLOGICAL_JUNCTION.FISH_NUMBER = CALORIMETRY.FISH_NUMBER) ON (STATION_INFO.STATION_ID = BRIDGE.STATION_ID) AND (STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID)
WHERE (((CALORIMETRY.DATA_ISSUE)='N') AND ((STATION_INFO.SYNOPTIC_STATION)=True) AND ((BRIDGE.MONTH)='JUN' Or (BRIDGE.MONTH)='JUL') AND ((BRIDGE.HEAD_DEPTH)<22) AND ((BRIDGE.START_BOT_DEPTH)<290) AND ((BRIDGE.END_BOT_DEPTH)<290));
                        ")

# get calorimetry data for IPES from high seas table
cal_ipes_orig <- sqlQuery(myconn_hs, "SELECT CALORIMETRY_IPES.FISH_NUMBER, CALORIMETRY_IPES.DUPLICATE, CALORIMETRY_IPES.HEAT_RELEASED_CAL, CALORIMETRY_IPES.HEAT_RELEASED_KJ, CALORIMETRY_IPES.DATA_ISSUE, CALORIMETRY_IPES.COMMENTS
FROM CALORIMETRY_IPES
WHERE (((CALORIMETRY_IPES.DATA_ISSUE)='N'));
                          ")


# # get GSI data from high seas
# 
# 
# # calorimetry data for IPES is stored in high seas
# cal_ipes_orig <- sqlQuery(myconn_hs, "SELECT CALORIMETRY_IPES.FISH_NUMBER, CALORIMETRY_IPES.HEAT_RELEASED_CAL, CALORIMETRY_IPES.HEAT_RELEASED_KJ, CALORIMETRY_IPES.DUPLICATE, CALORIMETRY_IPES.DATA_ISSUE
# FROM CALORIMETRY_IPES
# WHERE (((CALORIMETRY_IPES.DATA_ISSUE)="N"));
#")


# close database
close(myconn_hs)

#####################################
# wrangle high seas data
#####################################

# check for empty net dimensions and distance values
missing_hs <- cpue_hs_orig %>%
  filter(is.na(NET_OPENING_WIDTH)| is.na(NET_OPENING_HEIGHT)| is.na(DISTANCE))

#### all missing from cruise 201893 
# use Sea Crest values from gear comparison 
# according to head depth
# no missing distance values

# replace missing net values
cpue_hs_net <- cpue_hs_orig %>%
  mutate(NET_OPENING_WIDTH = case_when(
    !(is.na(NET_OPENING_WIDTH)) ~ NET_OPENING_WIDTH,
      is.na(NET_OPENING_WIDTH) & HEAD_DEPTH <= 7 ~ 41,
      is.na(NET_OPENING_WIDTH) & HEAD_DEPTH > 7 ~ 47),
    NET_OPENING_HEIGHT = case_when(
      !(is.na(NET_OPENING_HEIGHT)) ~ NET_OPENING_HEIGHT,
        is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH <= 7 ~ 19,
        is.na(NET_OPENING_HEIGHT) & HEAD_DEPTH > 7 ~ 13),
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
         logCPUE1 = log(CPUE + 1))

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
  
  # remove duplicate lines
  # they are zero tows anyway
    # likely due to blank juvenile counts with trace amounts in 2018
  distinct(TRIP_YEAR, BLOCK_DESIGNATION, EVENT, START_LATITUDE, START_LONGITUDE, JUVENILE_CATCH_COUNT,
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
tows_miss_ck <- cpue_ipes_ck %>%
  group_by(TRIP_YEAR, BLOCK_DESIGNATION, EVENT) %>%
  filter(n() > 1)

tows_miss_co <- cpue_ipes_cm %>%
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

# there are two tows misssing from 2018 and one from 2019 likely from deleted tows
# ignore those missing tows and use catch data

# make year with number of tows vector for label
# add 2018 IPES and 2018 high seas together
# use hs for 1998 to 2017
# use IPES for 2017 and 2019
yearTowVec <- c("1998\n(12)", "1999\n(12)", "2000\n(16)", "2001\n(23)", "2002\n(11)", "2003\n(7)",
                "2004\n(15)", "2005\n(14)","2006\n(17)", "2007\n(32)", "2008\n(0)", "2009\n(14)",
                "2010\n(19)", "2011\n(26)", "2012\n(23)","2013\n(16)", "2014\n(7)", "2015\n(33)",
                "2016\n(0)", "2017\n(54)", "2018\n(91)", "2019\n(68)")

# create function to make anomalies graph for each salmon species
anom_fn <- function(df, speciesCode, speciesName, yearTowVec, tows_ipes) {

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
    mutate(anom = (meanCPUE - meanCPUE_ts)/sdCPUE_ts)
  
  # make year factor for nice graph
  cpue_select$Year_fac <- as.factor(cpue_select$TRIP_YEAR)
  cpue_select$Year_fac <- factor(cpue_select$Year_fac, 
                                levels = c("1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005",
                                           "2006", "2007", "2008", "2009", "2010", "2011", "2012",
                                           "2013", "2014", "2015", "2016", "2017", "2018", "2019"))

  # graph
  plot <- ggplot(cpue_select, aes(x = Year_fac, y = anom)) +
    geom_bar(stat = "identity", fill = "darkred") +
    theme_bw() +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(-1, 1), linetype = "dashed") +
    labs(x = "Ocean Sample Year\n(Number of Tows)",
         y = "ln(CPUE + 1)",
         title = str_c("Juvenile ", speciesName, " Salmon Annual CPUE Anomalies")) +
    theme(title = element_text(face = "bold", size = 14)) +
    geom_vline(xintercept = c(11, 19), linetype = "dotted") +
    ylim(-2, 2) +
    scale_x_discrete(drop = FALSE,
                     labels = yearTowVec) + 
  theme(axis.title = element_text(face = "bold", size = 14))
  
  ggsave(str_c("Output/2019/CPUE_", speciesName, ".png"), plot)
  
  return(plot)

}

# apply function to species
cpuePK_plot <- anom_fn(cpue, 108, "Pink", yearTowVec, tows_ipes)
cpuePK_plot
cpueCM_plot <- anom_fn(cpue, 112, "Chum", yearTowVec, tows_ipes)
cpueCM_plot
cpueCO_plot <- anom_fn(cpue, 115, "Coho", yearTowVec, tows_ipes)
cpueCO_plot
cpueCK_plot <- anom_fn(cpue, 124, "Chinook", yearTowVec, tows_ipes)
cpueCK_plot
cpueSE_plot <- anom_fn(cpue, 118, "Sockeye", yearTowVec, tows_ipes)
cpueSE_plot

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
        plotTitle <- str_c(nameDf, " 2019 Distribution")
        
        # print to console
        print(plotTitle)
        
        # plot the individual interpolation surfaces
      ggplot() +
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
               x = "Longitude",
               y = "Latitude") 
        
        plotName <- paste0(OutputFolder, "/Kriging", i, ".png")
        
        ggsave(plotName)
        
        # add new data frame to list
        mylist[[nameDf]] <- TheSurface
      }

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

#####################################
# Calorimetry results
#####################################

# wrangle IPES data to correct fish_number
# join to bio data to limit to usable, daylight tows and juveniles only
cal_ipes <- cal_ipes_orig %>%
  left_join(., lw_ipes, by = "FISH_NUMBER")
# average accross duplicates




#####################################
# GSI results for 2019

#####################################
# display other species? 
# counts or biomass?


#####################################
