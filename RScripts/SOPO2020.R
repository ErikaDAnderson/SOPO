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
# Calorimetry results (if available)
# Genetic stock identification (if available)
# HSSALMON database limited to same area in fall only
# no IPES database needed since no fall surveys for comparison
# incorporate water profiles somehow as graphic?
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
library(modelr) # models length to weight, residuals
# library(rcompanion) # confident intervals
# library(readxl) # read excel files for GSI

#####################################
# database file path
# change these paths for your own version
#####################################

# connection to local high seas MS Access database
db_hs <- "C:/Users/andersoned/Documents/High Seas Salmon Database/HSSALMON.accdb"

outputFolder <- "Output/2020/"

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

# close database
close(myconn_hs)

#####################################
# wrangle high seas data
#####################################
# # save all tows as csv to display in ArcMap to make a raster
# write_csv(tows_hs_orig, file = str_c(outputFolder, "tows_hs_orig.csv"))

# after mapping csv of tows; I recommend removing some tows as not in regions
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
  
  # calculate mean and standard deviation for time series
  meanCPUE_ts <- mean(cpue_select$meanCPUE, na.rm = TRUE)
  sdCPUE_ts <- sd(cpue_select$meanCPUE, na.rm = TRUE)
  
  ###***** come back to this for anomalies by species???
  
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
    cpue_select <- cpue %>%
      filter(SPECIES_CODE == speciesCode) %>%
      group_by(TRIP_YEAR, REGION_CODE) %>%
      summarize(meanCPUE = mean(logCPUE1, na.rm = TRUE),
                .groups = "drop") 
    
    ###***** come back to this for anomalies by species and region
  }
  

  
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

# #####################################
# # use Kriging to see spatial distribution
# #####################################
# # folder name for Kriging outputs
# OutputKriging <- str_c(outputFolder, "Kriging", Sys.Date())
# 
# # create directory for plots
# dir.create(OutputKriging)
# 
# 
# # need grid for Kriging
# saFilename <- here::here("Input/Spatial/ipes_wgs84.tif")
# 
# # load raster for study area to crop interpolation grid
# sa <- raster(saFilename)
# class(sa)
# 
# # convert to spatial point data frame
# sa <- rasterToPoints(sa, spatial = TRUE)
# class(sa)
# 
# # convert to pixel data frame
# gridded(sa) = TRUE
# class(sa)
# 
# # check grid
# plot(sa)

#####################################
# load more high seas data
#####################################
# estalish connection to high sea Access database
myconn_hs <- odbcConnectAccess2007(db_hs)
                      
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


# #####################################
# # use Kriging to see spatial distribution
# #####################################
# # folder name for Kriging outputs
# OutputKriging <- str_c(outputFolder, "Kriging", Sys.Date())
# 
# # create directory for plots
# dir.create(OutputKriging)
# 
# 
# # need grid for Kriging

#####################################
# Length to weight residuals 
#####################################
# present as residuals over time series
# gives estimate of condition of salmon

# wrangle hs data
# remove adults
# use 350 mm as limit for all species
lw_hs <- lw_hs_orig %>%
  filter(SHIP_LENGTH < 350) %>%
  rename(TRIP_YEAR = YEAR,
         LENGTH = SHIP_LENGTH,
         WEIGHT = SHIP_WT) %>%
  filter(!is.na(SPECIES_CODE)) %>%
  dplyr::select(TRIP_YEAR, REGION_CODE, FISH_NUMBER, SPECIES_CODE, LENGTH, WEIGHT)

# removed rows with blanks
lw <- lw_hs %>%
  filter(!is.na(WEIGHT)) %>%
  filter(!is.na(LENGTH)) 

# plot prelim data
ggplot(lw, aes(log10(LENGTH), log10(WEIGHT), color = REGION_CODE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~SPECIES_CODE)

#### plot all regions together first (then repeat with region as factor in model)
# fit model
modLW <- lm(log10(WEIGHT) ~ log10(LENGTH) * factor(SPECIES_CODE), data = lw)

# add residuals
residsLW <- lw %>%
  add_residuals(., modLW) %>%
  rename(Year = TRIP_YEAR,
         Residuals = resid)

# expand grid for zero years
yearGrid <- expand_grid(Year = 1996:2020) %>%
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
                        levels = c("1996", "1997", "1998", "1999", "2000", "2001", 
                                   "2002","2003", "2004", "2005", "2006", "2007", 
                                   "2008", "2009", "2010", "2011", "2012", "2013", 
                                   "2014", "2015", "2016", "2017", "2018", "2019", 
                                   "2020"))

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
boxplot_pk 
# +
#   scale_x_discrete(drop = FALSE,
#                    labels = c("1996\n(0)", "1997\n(0)", "1998\n(2)", "1999\n(28)", "2000\n(339)", 
#                               "2001\n(310)", "2002\n(205)","2003\n(122)", "2004\n(110)", "2005\n(43)",
#                               "2006\n(20)", "2007\n(112)", "2008\n(207)","2009\n(104)", "2010\n(71)",
#                               "2011\n(229)", "2012\n(81)", "2013\n(1)","2014\n(0)", "2015\n(0)", 
#                               "2016\n(0)", "2017\n(49)", "2018\n(0)", "2019\n(13)", "2020\n(166)" ))

ggsave(str_c(outputFolder, "LW_Pink.png"))

# chum LW plot
boxplot_cm <- LWboxplot_fn(residsLW, 112, "Chum")
boxplot_cm 
# +
#   scale_x_discrete(drop = FALSE,
#                    labels = c("1996\n(0)", "1997\n(0)", "1998\n(6)", "1999\n(47)", "2000\n(342)", 
#                               "2001\n(86)", "2002\n(113)","2003\n(66)", "2004\n(135)", "2005\n(44)",
#                               "2006\n(11)", "2007\n(105)", "2008\n(198)","2009\n(96)", "2010\n(85)",
#                               "2011\n(121)", "2012\n(48)", "2013\n(5)","2014\n(0)", "2015\n(0)", 
#                               "2016\n(0)", "2017\n(158)", "2018\n(0)", "2019\n(18)", "2020\n(184)" ))

ggsave(str_c(outputFolder, "LW_Chum.png"))

# coho LW plot
boxplot_co <- LWboxplot_fn(residsLW, 115, "Coho")
boxplot_co 
ggsave(str_c(outputFolder, "LW_Coho.png"))

# sockeye LW plot
boxplot_se <- LWboxplot_fn(residsLW, 118, "Sockeye")
boxplot_se 
ggsave(str_c(outputFolder, "LW_Sockeye.png"))

# chinook LW plot
boxplot_ck <- LWboxplot_fn(residsLW, 124, "Chinook")
boxplot_ck 
ggsave(str_c(outputFolder, "LW_Chinook.png"))

##### repeat with region within model

# fit model
modLW_Reg <- lm(log10(WEIGHT) ~ log10(LENGTH) * factor(SPECIES_CODE) * factor(REGION_CODE), data = lw)

# look at model
broom::glance(modLW_Reg)

# add residuals
residsLW_Reg <- lw %>%
  add_residuals(., modLW_Reg) %>%
  rename(Year = TRIP_YEAR,
         Residuals = resid)


# make years as factors for graphs
residsLW_Reg$Year <- factor(residsLW_Reg$Year, 
                        levels = c("1996", "1997", "1998", "1999", "2000", "2001", 
                                   "2002","2003", "2004", "2005", "2006", "2007", 
                                   "2008", "2009", "2010", "2011", "2012", "2013", 
                                   "2014", "2015", "2016", "2017", "2018", "2019", 
                                   "2020"))

# make function to graph LW as boxplots
LWboxplotReg_fn <- function(df, speciesCode, speciesName) {
  
  # select species
  df_select <- df %>%
    filter(SPECIES_CODE == speciesCode)
  
  # plot title
  plotTitle <- str_c(speciesName, " Length to Weight Relationship")
  
  # graph
  ggplot(data = df_select, aes(x = Year, y = Residuals, group = Year)) + 
    geom_boxplot(fill = "darkred") + 
    scale_y_continuous(expand = c(0, 0), limits = c(-.22, .22)) +
    theme(panel.grid.minor.y = element_blank(), 
          panel.background = element_rect(fill = "white",colour = "black"), 
          title = element_text(face = "bold", size = 14),
          strip.background = element_rect(fill = "white")) +
    labs(title = plotTitle,
         Y = "",
         x = "Year") +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    facet_wrap(~REGION_CODE, ncol = 1) 
}

# pink LW plot
boxplot_pk_reg <- LWboxplotReg_fn(residsLW_Reg, 108, "Pink")
boxplot_pk_reg 
ggsave(str_c(outputFolder, "LW_Pink_Reg.png"))

# chum LW plot
boxplot_cm_reg <- LWboxplotReg_fn(residsLW_Reg, 112, "Chum")
boxplot_cm_reg 
ggsave(str_c(outputFolder, "LW_Chum_Reg.png"))

# coho LW plot
boxplot_co_reg <- LWboxplotReg_fn(residsLW_Reg, 115, "Coho")
boxplot_co_reg 
ggsave(str_c(outputFolder, "LW_Coho_Reg.png"))

# sockeye LW plot
boxplot_se_reg <- LWboxplotReg_fn(residsLW_Reg, 118, "Sockeye")
boxplot_se_reg 
ggsave(str_c(outputFolder, "LW_Sockeye_Reg.png"))

# chinook LW plot
boxplot_ck_reg <- LWboxplotReg_fn(residsLW_Reg, 124, "Chinook")
boxplot_ck_reg 
ggsave(str_c(outputFolder, "LW_Chinook_Reg.png"))

# How about we graph the residuals for 2020 faceted by region and grouped by species
residsLW_Reg %>%
  filter(Year == 2020) %>%
  ggplot(aes(x = as.factor(SPECIES_CODE), y = Residuals, group = SPECIES_CODE)) + 
  geom_boxplot(fill = "darkred") + 
  scale_y_continuous(expand = c(0, 0)) +
#, limits = c(-.22, .22)) +
  theme(panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"), 
        title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white")) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  facet_wrap(~REGION_CODE, ncol = 1) +
  labs(y = "Residuals",
       x = "Species") 

# make labels better for species and regions

       