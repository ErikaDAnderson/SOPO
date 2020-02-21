#=====================================================================================================
# Script Name: SOPO2019.R
# Script Author: Erika Anderson
# Script Date: 2020-02
# R Version: 3.6.1
#
# State of the Pacific Ocean March 2020
# swept volume CPUE and genetic stock identification
# from high seas salmon database for June along migratory corridor  
# Johnstone Strait, Queen Charlotte Strait and southern Queen Charlotte Sound,
# see email saved in Input for more details about previous years' info
#
#=====================================================================================================

# Load packages
library(tidyverse) # core data manipulation & visualization
library(magrittr) # save as same dataframe
library(RODBC) # database connection to query
#library(stringr) # string manipulation
#library(plyr) # used in older code for data manipulation
library(bootstrap) # error bars
library(viridis) # colors for color blind people

#####################################

# estalish connection to high sea Access database
db <- "C:/Users/andersoned/Documents/BCSI/High Seas Salmon Database/HSSALMON.accdb"
myconn <- odbcConnectAccess2007(db)

#####################################
# triangle transect
# create GSI profile along longitudes to see if there is natural break

# parameters for species
SE <- "'sockeye'"
CO <- "'coho'"
CK <- "'chinook'"
PK <- "'pink'"
CM <- "'chum'"

# pull triangle transect GSI data in June from Access using query
triangle_original <- sqlQuery(myconn, "SELECT DNA_CHINOOK_STOCK_ID.FISH_NUMBER, DNA_CHINOOK_STOCK_ID.STOCK_1, DNA_CHINOOK_STOCK_ID.REGION_CODE_1, DNA_CHINOOK_STOCK_ID.REGION_1, DNA_CHINOOK_STOCK_ID.PROB_1, TriangleTransectStations_June.STATION_ID, TriangleTransectStations_June.START_LAT, [START_LONG]*-1 AS START_LONG_NEG, ${CK} AS SPECIES
FROM (STATION_INFO INNER JOIN TriangleTransectStations_June ON STATION_INFO.STATION_ID = TriangleTransectStations_June.STATION_ID) INNER JOIN (BIOLOGICAL_JUNCTION INNER JOIN DNA_CHINOOK_STOCK_ID ON (BIOLOGICAL_JUNCTION.FISH_NUMBER = DNA_CHINOOK_STOCK_ID.FISH_NUMBER) AND (BIOLOGICAL_JUNCTION.FISH_NUMBER = DNA_CHINOOK_STOCK_ID.FISH_NUMBER)) ON STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID
UNION
SELECT DNA_CHUM_STOCK_ID.FISH_NUMBER, DNA_CHUM_STOCK_ID.STOCK_1, DNA_CHUM_STOCK_ID.REGION_CODE_1, DNA_CHUM_STOCK_ID.REGION_1, DNA_CHUM_STOCK_ID.PROB_1, TriangleTransectStations_June.STATION_ID, TriangleTransectStations_June.START_LAT, [START_LONG]*-1 AS START_LONG_NEG, ${CM} AS SPECIES
FROM DNA_CHUM_STOCK_ID INNER JOIN ((STATION_INFO INNER JOIN TriangleTransectStations_June ON STATION_INFO.STATION_ID = TriangleTransectStations_June.STATION_ID) INNER JOIN BIOLOGICAL_JUNCTION ON STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) ON DNA_CHUM_STOCK_ID.FISH_NUMBER = BIOLOGICAL_JUNCTION.FISH_NUMBER
UNION
SELECT DNA_COHO_STOCK_ID.FISH_NUMBER, DNA_COHO_STOCK_ID.STOCK_1, DNA_COHO_STOCK_ID.REGION_CODE_1, DNA_COHO_STOCK_ID.REGION_1, DNA_COHO_STOCK_ID.PROB_1, TriangleTransectStations_June.STATION_ID, TriangleTransectStations_June.START_LAT, [START_LONG]*-1 AS START_LONG_NEG, ${CO} AS SPECIES
FROM DNA_COHO_STOCK_ID INNER JOIN ((STATION_INFO INNER JOIN TriangleTransectStations_June ON STATION_INFO.STATION_ID = TriangleTransectStations_June.STATION_ID) INNER JOIN BIOLOGICAL_JUNCTION ON STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) ON DNA_COHO_STOCK_ID.FISH_NUMBER = BIOLOGICAL_JUNCTION.FISH_NUMBER
UNION
SELECT DNA_PINK_STOCK_ID.FISH_NUMBER, DNA_PINK_STOCK_ID.STOCK_1, DNA_PINK_STOCK_ID.REGION_CODE_1, DNA_PINK_STOCK_ID.REGION_1, DNA_PINK_STOCK_ID.PROB_1, TriangleTransectStations_June.STATION_ID, TriangleTransectStations_June.START_LAT, [START_LONG]*-1 AS START_LONG_NEG, ${PK} AS SPECIES
FROM DNA_PINK_STOCK_ID INNER JOIN ((STATION_INFO INNER JOIN TriangleTransectStations_June ON STATION_INFO.STATION_ID = TriangleTransectStations_June.STATION_ID) INNER JOIN BIOLOGICAL_JUNCTION ON STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) ON DNA_PINK_STOCK_ID.FISH_NUMBER = BIOLOGICAL_JUNCTION.FISH_NUMBER
UNION
SELECT DNA_SOCKEYE_STOCK_ID.FISH_NUMBER, DNA_SOCKEYE_STOCK_ID.STOCK_1, DNA_SOCKEYE_STOCK_ID.REGION_CODE_1, DNA_SOCKEYE_STOCK_ID.REGION_1, DNA_SOCKEYE_STOCK_ID.PROB_1, TriangleTransectStations_June.STATION_ID, TriangleTransectStations_June.START_LAT, [START_LONG]*-1 AS START_LONG_NEG, ${SE} AS SPECIES
FROM DNA_SOCKEYE_STOCK_ID INNER JOIN ((STATION_INFO INNER JOIN TriangleTransectStations_June ON STATION_INFO.STATION_ID = TriangleTransectStations_June.STATION_ID) INNER JOIN BIOLOGICAL_JUNCTION ON STATION_INFO.STATION_ID = BIOLOGICAL_JUNCTION.STATION_ID) ON DNA_SOCKEYE_STOCK_ID.FISH_NUMBER = BIOLOGICAL_JUNCTION.FISH_NUMBER
;
")

# remove NAs
triangle <- triangle_original %>%
  filter(., !is.na(PROB_1)) %>%
  
  # remove rows less than 50% probability
  filter(PROB_1 > 0.5) 

# load allocation table
allocation <- read.csv("Input/stockAllocation3.csv")

# join to assign stocks to prob dataset
triangle %<>%
  left_join(., allocation, by = c("STOCK_1" = "Stock"))

# find stocks within region of origin
allocateMissing <- triangle %>%
  filter(., is.na(Origin)) %>%
  group_by(STOCK_1, REGION_1, SPECIES) %>%
  dplyr::count(., SPECIES) %>%
  arrange(desc(n))

#### repeat until allocateMissing has zero data
### use stock names to search google maps if not in Beacham's supplementary tables

# drop levels not in data
triangle$Origin <- droplevels(triangle$Origin)

# check to see unique Origins
unique(triangle$Origin)

# relevel factor to have N-S order in plot
triangle$Origin <- factor(triangle$Origin,
                          levels = c("North Coast", "Central Coast",
                                     "ECVI", "Fraser", "South Coast", "Puget Sound", "WCVI",
                                     "Washington", "Columbia", "Snake", "Oregon"))

# plot 
ggplot(triangle, aes(x = START_LONG_NEG, fill = Origin)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~ SPECIES, nrow = 3, scales = "free_y") +
  labs(y = "Number of juvenile salmon",
       x = "Longitude",
       fill = "Origin",
       title = "Genetic Stock ID for Triangle Transect") +
  scale_x_continuous(breaks = c(-129.5, -129.25, -129.0, -128.75, -128.5, -128.25)) +
  theme_bw()

# save plot
ggsave("Output/June_GSI_traingle.png")

# plot for color blind people
ggplot(triangle, aes(x = START_LONG_NEG, fill = Origin)) +
  geom_histogram(binwidth = 0.05) +
  facet_wrap(~ SPECIES, nrow = 3, scales = "free_y") +
  labs(y = "Number of juvenile salmon",
       x = "Longitude",
       fill = "Origin",
       title = "Genetic Stock ID for Triangle Transect") +
  scale_x_continuous(breaks = c(-129.5, -129.25, -129.0, -128.75, -128.5, -128.25)) +
  scale_fill_viridis_d() +
  theme_bw()

# save plot
ggsave("Output/June_GSI_traingle_viridis.png")

#####################################