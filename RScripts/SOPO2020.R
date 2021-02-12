#=====================================================================================================
# Script Name: SOPO2020.R
# Script Author: Erika Anderson
# Script Start Date: 2021-02
# R Version: 3.6.1
#
# State of the Pacific Ocean meeting March 2021
# HSSALMON database limited to same areas (DE, HS, QCSD) in fall (SEP, OCT, NOV)
# no IPES database needed since no fall surveys
# CPUE anomalies by species and region using swept volume
# Length to weight residuals by species and region (QCSD, HS, DE)

# Kriging of CPUE from Queen Charlotte Sound to Dixon Entrance
# Calorimetry results (if available)
# Genetic stock identification (if available)
# incorporate water profiles somehow as graphic?
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

# library(viridis) # colors graphs
# library(sf) # spatial manipulation (newer than sp) so works with ggplot2
# library(sp) # spatial data manipulation
# library(rgdal) # to load shapefiles and rasters
# library(gstat) # model fit & Krige interpolation
# library(data.table) # bind data frames together from list
# library(raster) # load raster for grid (and predict function, alternative to gstat krige)
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
# one query had 15 m instead of 20 m head rope depths - fixed
# anyDuplicated(cpue_hs_se$STATION_ID)
# anyDuplicated(cpue_hs_pk$ STATION_ID)
# 
# unique(cpue_hs_pk$STATION_ID[!cpue_hs_pk$STATION_ID %in% cpue_hs_ck$STATION_ID])
# unique(cpue_hs_se$STATION_ID[!cpue_hs_se$STATION_ID %in% cpue_hs_ck$STATION_ID])

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
  group_by(TRIP_YEAR, STATION_ID) %>%
  count() %>%
  ungroup() %>%
  group_by(TRIP_YEAR) %>%
  count()

# create function to calculate anomalies for each salmon species
# add ability to calculate regions seperately
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
    
    dataOld <- cpue_Reg %>%
      filter(REGION_CODE == regionName) %>%
      mutate(anom = if_else(Year_fac == 2020, NA_real_, anom))
    
    dataCurrent <- cpue_Reg %>%
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
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        strip.background = element_rect(fill = "white"),
        plot.subtitle = element_text(hjust = 0.5)) +
      labs(subtitle = regionNameLong)
  }
}

# apply function overall and to each region
cpueAll <- cpueGraphFN(cpue_noReg, NA) 
cpueAll
ggsave(str_c(outputFolder, "CPUE_NorthCoast2020", str_replace_all(Sys.Date(), "-", ""),".png"))

# combine plots
cpueDE <- cpueGraphFN(cpue_Reg, "DE") 
cpueDE
cpueHS <- cpueGraphFN(cpue_Reg, "HS") 
cpueHS
cpueQCSD <- cpueGraphFN(cpue_Reg, "QCSD") 
cpueQCSD

# combine into one plot
cpueRegions <- cpueDE / cpueHS / cpueQCSD

ggsave(str_c(outputFolder, "CPUE_Regions2020", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 9, units = "in")

#####################################
# load length weight high seas data
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


# close database
close(myconn_hs)

#####################################
# Length to weight residuals 
#####################################
# residuals over time series
# gives estimate of condition of salmon

# wrangle hs data
# remove adults
# definition of juveniles by species for fall from bridge table
# SE PK CM < 300 mm
# CK < 350 mm
# CO < 400 mm
lw_hs <- lw_hs_orig %>%
  mutate(AGE = case_when(
    SPECIES_CODE == 108 & SHIP_LENGTH < 300 ~ "J",
    SPECIES_CODE == 112 & SHIP_LENGTH < 300 ~ "J",
    SPECIES_CODE == 118 & SHIP_LENGTH < 300 ~ "J",
    SPECIES_CODE == 124 & SHIP_LENGTH < 350 ~ "J",
    SPECIES_CODE == 115 & SHIP_LENGTH < 400 ~ "J",
    TRUE ~ "A")) %>%
  filter(AGE == "J") %>%
  rename(TRIP_YEAR = YEAR,
         LENGTH = SHIP_LENGTH,
         WEIGHT = SHIP_WT) %>%
  filter(!is.na(SPECIES_CODE)) %>%
  dplyr::select(TRIP_YEAR, REGION_CODE, FISH_NUMBER, SPECIES_CODE, LENGTH, WEIGHT)

# removed rows with blanks
lw <- lw_hs %>%
  filter(!is.na(WEIGHT)) %>%
  filter(!is.na(LENGTH)) %>%
  
  # remove outliers from sockeye until they are fact checked
  ### TZ is checking in office for new updates
  ### already corrected BCSI-2019125-T02-007 from length 102 to 162 in local db
  filter(FISH_NUMBER != "BCSI-2020017-HN04-118J-001") %>%
  filter(FISH_NUMBER != "HS201145-T06-118-003")

# plot prelim data
ggplot(lw, aes(log10(LENGTH), log10(WEIGHT), 
               #color = as.factor(SPECIES_CODE), 
               group = SPECIES_CODE)) +
  geom_point() +
  geom_smooth(method = "lm", formula = 'y~x') +
  facet_wrap(~SPECIES_CODE) +
  theme_bw()

ggsave(str_c(outputFolder, "LW_linearModels", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 4, units = "in")

# fit models
modLW_simp <- lm(log10(WEIGHT) ~ log10(LENGTH), data = lw)
modLW <- lm(log10(WEIGHT) ~ log10(LENGTH) * factor(SPECIES_CODE), data = lw)
modLW_Reg <- lm(log10(WEIGHT) ~ log10(LENGTH) * factor(SPECIES_CODE) * factor(REGION_CODE), data = lw)
broom::glance(modLW,modLW_Reg)
anova(modLW_simp, modLW, modLW_Reg) 
# species and region improves to model
# see later that graphically has little difference between species and species-region models

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
    filter(SPECIES_CODE == speciesCode) %>%
  
  # calculate number of fish per year
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
          panel.background = element_rect(fill = "white",colour = "black")) + 
    scale_y_continuous(limits = c(-.22, .22)) +
    theme(title = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14)) +
    labs(title = plotTitle,
         Y = "",
         x = "Year") +
    geom_hline(yintercept = 0, color = "black", size = 1) 
}

# pink LW plot
boxplot_pk <- LWboxplot_fn(residsLW, 108, "Pink")
boxplot_pk 
## use if adding number of tows manually for paper
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

# add residuals and account for years without 
residsLW_Reg <- lw %>%
  add_residuals(., modLW_Reg) %>%
  rename(Year = TRIP_YEAR,
         Residuals = resid) 

# expand grid for year without values to be represented on graph
residsLW_Reg_exp <- residsLW_Reg %>%
  full_join(., yearGrid, by = "Year") 

# make function to graph LW as boxplots
LWboxplotReg_fn <- function(df, speciesCode, speciesName) {
  
  # select species
  df_select <- df %>%
    filter(SPECIES_CODE == speciesCode)
  
  df_historic <- df_select %>%
    mutate(Residuals = if_else(Year == 2020, NA_real_, Residuals))
  
  df_current <- df_select %>%
    mutate(Residuals = if_else(Year != 2020, NA_real_, Residuals))
    
  
  # plot title
  plotTitle <- str_c(speciesName, " Length to Weight Relationship")
  
  # graph
  ggplot(data = df_historic, aes(x = Year, y = Residuals, group = Year)) + 
    geom_boxplot(fill = "grey") + 
    geom_boxplot(data = df_current, aes(x = Year, y = Residuals, group = Year),
                 fill = "darkred") + 
    #scale_y_continuous(limits = c(-.22, .22)) +
    theme(title = element_text(face = "bold", size = 14),
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 14),
          panel.grid.minor.y = element_blank(), 
          panel.background = element_rect(fill = "white",colour = "black"), 
          strip.background = element_rect(fill = "white")) +
    labs(title = plotTitle,
         Y = "",
         x = "Year") +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    facet_wrap(~REGION_CODE, ncol = 1, scales = "free_y") 
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
# after graphing both models; they look similar so can use either
residsLW %>%
  filter(Year == 2020) %>%
  mutate(SpeciesName = case_when(
    SPECIES_CODE == 108 ~ "Pink",
    SPECIES_CODE == 112 ~ "Chum",
    SPECIES_CODE == 115 ~ "Coho",
    SPECIES_CODE == 118 ~ "Sockeye",
    SPECIES_CODE == 124 ~ "Chinook"),
    RegionName = case_when(
      REGION_CODE == "DE" ~ "Dixon Entrance",
      REGION_CODE == "HS" ~ "Hecate Strait",
      REGION_CODE == "QCSD" ~ "Queen Charlotte Sound")) %>%
  
  ggplot(aes(x = as.factor(SpeciesName), y = Residuals, group = SpeciesName)) + 
  geom_boxplot(fill = "darkred") + 
  theme(panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"), 
        title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  facet_wrap(~RegionName, ncol = 1, scales = "free_y") +
  labs(y = "Residuals",
       x = "Species") 

ggsave(str_c(outputFolder, "LW2020All", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 6, units = "in")


residsLW_Reg %>%
  filter(Year == 2020) %>%
  mutate(SpeciesName = case_when(
    SPECIES_CODE == 108 ~ "Pink",
    SPECIES_CODE == 112 ~ "Chum",
    SPECIES_CODE == 115 ~ "Coho",
    SPECIES_CODE == 118 ~ "Sockeye",
    SPECIES_CODE == 124 ~ "Chinook"),
    RegionName = case_when(
      REGION_CODE == "DE" ~ "Dixon Entrance",
      REGION_CODE == "HS" ~ "Hecate Strait",
      REGION_CODE == "QCSD" ~ "Queen Charlotte Sound")) %>%
  
  ggplot(aes(x = as.factor(SpeciesName), y = Residuals, group = SpeciesName)) + 
  geom_boxplot(fill = "darkred") + 
  theme(panel.grid.minor.y = element_blank(), 
        panel.background = element_rect(fill = "white",colour = "black"), 
        title = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "white"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  geom_hline(yintercept = 0, color = "black", size = 1) +
  facet_wrap(~RegionName, ncol = 1, scales = "free_y") +
  labs(y = "Residuals",
       x = "Species") 

ggsave(str_c(outputFolder, "LW2020_Reg", str_replace_all(Sys.Date(), "-", ""),".png"),
       height = 6, units = "in")

#####################################
# load calorimetry high seas data
#####################################
# estalish connection to high sea Access database
myconn_hs <- odbcConnectAccess2007(db_hs)

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
cal_historic_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.Year, BRIDGE_COMPLETE.Month, BIOLOGICAL_ALL.SPECIES_CODE, BIOLOGICAL_JUNCTION.FISH_NUMBER, PROXIMATE_FISH.ENERGY_BOMB, PROXIMATE_FISH.ENERGY_BOMB_BLIND_DUPL, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.HEAD_DEPTH, BIOLOGICAL_ALL.SHIP_LENGTH, BIOLOGICAL_ALL.SHIP_WT
FROM ((BIOLOGICAL_JUNCTION INNER JOIN BRIDGE_COMPLETE ON BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE_COMPLETE.STATION_ID) INNER JOIN PROXIMATE_FISH ON BIOLOGICAL_JUNCTION.FISH_NUMBER = PROXIMATE_FISH.FISH_NUMBER) INNER JOIN BIOLOGICAL_ALL ON BIOLOGICAL_JUNCTION.FISH_NUMBER = BIOLOGICAL_ALL.FISH_NUMBER
WHERE (((BRIDGE_COMPLETE.Month)='SEP' Or (BRIDGE_COMPLETE.Month)='OCT' Or (BRIDGE_COMPLETE.Month)='NOV') AND ((BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE' Or (BRIDGE_COMPLETE.REGION_CODE)='QCSD') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<=20));
                              ")

# close database
close(myconn_hs)

#####################################
# wrangle calorimetry data
######################################

# check that all juveniles from lengths

# wrangle high seas data
# only 6 sockeye completed from 2019 from QCSD
cal_bcsi <- cal_hs_orig %>%
  # average accross duplicates
  group_by(FISH_NUMBER) %>%
  summarize(HEAT_RELEASED_CAL = mean(HEAT_RELEASED_CAL, na.rm = TRUE),
            HEAT_RELEASED_KJ = mean(HEAT_RELEASED_KJ, na.rm = TRUE),
            .groups = "drop") %>%
  # add species codes and length and weights
  left_join(., residsLW, by = "FISH_NUMBER") %>%
  filter(!is.na(SPECIES_CODE)) %>%
  dplyr::select(FISH_NUMBER, HEAT_RELEASED_KJ, Year, SPECIES_CODE, REGION_CODE, Residuals)

# wrangle historic
cal_historic <- cal_historic_orig %>%
  mutate(HEAT_RELEASED_KJ = if_else(is.na(ENERGY_BOMB_BLIND_DUPL), 
                               ENERGY_BOMB, (ENERGY_BOMB + ENERGY_BOMB_BLIND_DUPL)/2),
    LENGTH = SHIP_LENGTH,
    WEIGHT = SHIP_WT) %>%
  add_residuals(., modLW) %>%
  rename(Residuals = resid) %>%
  dplyr::select(FISH_NUMBER, HEAT_RELEASED_KJ, Year, SPECIES_CODE, REGION_CODE, Residuals)

residsLW <- lw %>%
  add_residuals(., modLW) 

# bind historic to current calorimetry data
cal_all <- rbind(cal_bcsi, cal_historic)

# compare residuals to heat released values
ggplot(data = cal_all, 
       aes(HEAT_RELEASED_KJ, Residuals, color = REGION_CODE)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~SPECIES_CODE) +
  theme_bw()

# historic samples show length to weight condition factor not correlated to ED

#####################################
# get GSI available data
#####################################

# estalish connection to high sea Access database
myconn_hs <- odbcConnectAccess2007(db_hs)

# get gsi data
# regions QCSD, HS, DE
# fall only (SEP, OCT, NOV)
# headrope depth <=20 m to omit deeper tows
gsi_orig <- sqlQuery(myconn_hs, "SELECT BRIDGE_COMPLETE.CRUISE, BRIDGE_COMPLETE.YEAR, BRIDGE_COMPLETE.REGION_CODE, BRIDGE_COMPLETE.HEAD_DEPTH, DNA_ALL_RESULTS_STOCKS.Species, BIOLOGICAL_JUNCTION.FISH_NUMBER, DNA_ALL_RESULTS_STOCKS.STOCK_FINAL, DNA_ALL_RESULTS_STOCKS.STOCK_1, DNA_ALL_RESULTS_STOCKS.PROB_1, DNA_ALL_RESULTS_STOCKS.STOCK_OVER50, DNA_STOCK_STOCK_FINAL_REGION.HS_REGION
FROM (DNA_ALL_RESULTS_STOCKS INNER JOIN (BIOLOGICAL_JUNCTION INNER JOIN BRIDGE_COMPLETE ON BIOLOGICAL_JUNCTION.STATION_ID = BRIDGE_COMPLETE.STATION_ID) ON DNA_ALL_RESULTS_STOCKS.FISH_NUMBER = BIOLOGICAL_JUNCTION.FISH_NUMBER) LEFT JOIN DNA_STOCK_STOCK_FINAL_REGION ON DNA_ALL_RESULTS_STOCKS.STOCK_FINAL = DNA_STOCK_STOCK_FINAL_REGION.STOCK_FINAL
WHERE (((BRIDGE_COMPLETE.YEAR)=2019 Or (BRIDGE_COMPLETE.YEAR)=2020) AND ((BRIDGE_COMPLETE.REGION_CODE)='QCSD' Or (BRIDGE_COMPLETE.REGION_CODE)='HS' Or (BRIDGE_COMPLETE.REGION_CODE)='DE') AND ((BRIDGE_COMPLETE.HEAD_DEPTH)<20));
                          ")
# close database
close(myconn_hs)

# # could use excel files since not complete in database
# ## or add regions to link North Coast stocks
# PID20200096_BCSI_B90(20)_sc242_2020-12-09.xlsx
# PID20200096_BCSI_Batch_92(20)_sc40_2021-01-20.xlsx

# find unique stocks and regions
gsi_for_db <- gsi_orig %>%
  filter(is.na(HS_REGION)) %>%
  select(Species, STOCK_1, HS_REGION) %>%
  distinct()

# # write into csv to enter HS_REGIONS
# write_csv(gsi_for_db, str_c(outputFolder, "gsi_for_db.csv"),
#           na = "")

# pull in high sea regions for new stocks from csv
gsi_regions <- read_csv(str_c(outputFolder, "gsi_for_db.csv"),
                        col_types = 
                          cols(
                            Species = col_character(),
                            STOCK_1 = col_character(),
                            STOCK_1_UPPER = col_character(),
                            HS_REGION = col_character()
                          ))
                    
#####################################
# wrangle GSI data
#####################################


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
