setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load necessary libraries , Set Census API key
library(sf)
library(nominatimlite)
library(tidycensus)
library(tidyverse)
library(tmap)
library(tigris)
library(fastDummies)
library(dplyr)
census_api_key('7deb96b03ae11f23a2fb544839ee195ccd646ac3')


# 1. Socio-economic Data
## Get Census data for block groups
blockgroup <- get_acs(
  geography = "block_group", 
  variables = c(
    'tot_popE' = 'B01001_001E', 
    'male_popE' = 'B01001_002E', 
    'female_popE' = 'B01001_026E', 
    'race_totalE' = 'B02001_001E', 
    'whiteE' = 'B02001_002E', 
    'blackE' = 'B02001_003E', 
    'nativeE' = 'B02001_004E', 
    'asianE' = 'B02001_005E', 
    'pacific_islanderE' = 'B02001_006E', 
    'otherE' = 'B02001_007E'
  ),
  year = 2022,
  state = "GA", 
  county = c("Fulton", "DeKalb", "Clayton"), 
  survey = "acs5", 
  geometry = TRUE, 
  output = "wide"
)

## Calculate race ratios and population density
blockgroup <- blockgroup %>%
  select(GEOID, NAME, geometry, ends_with("E")) %>%
  mutate(
    pop_den = tot_popE / st_area(.),
    white_ratio = whiteE / race_totalE,
    black_ratio = blackE / race_totalE,
    other_ratio = (nativeE + asianE + pacific_islanderE + otherE) / race_totalE,
  )

## Get tract-level data
## There are not data on education level and median income for blockgroup, thus Get tract level data and put it to corresponding block groups.
tract_data <- get_acs(
  geography = "tract", 
  variables = c(
    'median_incomeE' = 'B06011_001E', 
    'edu_totalE' = 'B06009_001E', 
    'less_than_hsE' = 'B06009_002E'
  ),
  year = 2022,
  state = "GA", 
  county = c("Fulton", "DeKalb", "Clayton", "Cobb"), 
  survey = "acs5", 
  output = "wide"
) %>%
  mutate(less_than_hs_ratio = less_than_hsE / edu_totalE) %>%
  select(GEOID, median_incomeE, less_than_hs_ratio)

## Match and join to block group data
blockgroup <- blockgroup %>%
  mutate(tract_id = substr(GEOID, 1, 11)) %>%
  left_join(tract_data, by = c("tract_id" = "GEOID"))

## Clean up by dropping the helper tract_id column
blockgroup <- blockgroup %>%
  select(-tract_id)



# 2. Filter block groups within Atlanta City Region

## Get City of Atlanta boundary
atlanta <- nominatimlite::geo_lite_sf('Atlanta, GA', points_only = FALSE)
atlanta <- st_transform(atlanta, st_crs(blockgroup))

## Filter block group data within Atlanta city region
blockgroup_atlanta <- blockgroup[st_intersects(blockgroup, atlanta, sparse = FALSE), ]





# 3. Crime Data

## Load crime data and categorize violent and non-violent crimes
crime.data <- read.csv("Data/Crime_Data.csv")
crime.data <- crime.data %>%
  mutate(violent = ifelse(NIBRS.Code.Name %in% c(
    "Aggravated Assault", "Murder & Nonnegligent Manslaughter", "Sodomy", 
    "Animal Cruelty", "Statutory Rape", "Rape", "Fondling", "Arson", 
    "Kidnapping/Abduction", "Intimidation", "Simple Assault", 
    "Weapon Law Violations"), 1, 0))

## Convert crime data to sf object
crime.sf <- st_as_sf(crime.data, coords = c("Longitude", "Latitude"), crs = 4326)

## Ensure both objects have the same CRS
blockgroup_atlanta <- st_transform(blockgroup_atlanta, crs = st_crs(crime.sf))

## Perform spatial join
bg_crime <- st_join(blockgroup_atlanta, crime.sf %>% mutate(count = ifelse(violent == 0, 1, 10000)))

## Summarize crime data by block group
bg_crime_count <- bg_crime %>%
  group_by(GEOID) %>%
  summarise(count = sum(count, na.rm = TRUE)) %>% 
  st_drop_geometry()

blockgroup_atlanta <- blockgroup_atlanta %>%
  left_join(bg_crime_count, by = "GEOID") %>% 
  mutate(violent_count = count %/% 10000, nonviolent_count = count %% 10000)




# 4. Land-use Data

## Load zoning data and classify zones
zoning <- read.csv("Data/Official_Zoning_Districts_-_Open_Data.csv")
zoning <- zoning %>%
  mutate(class = case_when(
    str_detect(ZONECLASS, "^(R-(1|2|3|4|5)|PD-H|FCR-)") ~ "Low-density Residential",
    str_detect(ZONECLASS, "^(R-LC|MRC-|PD-MU)") ~ "Residential-Commercial",
    str_detect(ZONECLASS, "^(RG-|MR-)") ~ "High-density Residential",
    str_detect(ZONECLASS, "^I-") ~ "Industrial",
    str_detect(ZONECLASS, "^(C-|LD|LW|NC|PD-OC)") ~ "Commercial",
    str_detect(ZONECLASS, "^O-") ~ "Institutional",
    TRUE ~ ""
  ))

## Load land use data and calculate land use percentages by block group
### NOTE : It took too much time. We made geojson file already by this command.
'''
landuse <- st_read("Data/Current_Land_Use.geojson") %>%
  select(FID, ACRES, LandUse, geometry) %>%
  st_make_valid() %>%
  st_transform(st_crs(blockgroup))

landuse <- landuse %>%
  st_intersection(blockgroup) %>%
  group_by(GEOID) %>%
  summarise(area = sum(st_area(geometry)))
'''

landuse <- st_read("Data/bg_landuse.geojson")
columns_to_keep <- c("GEOID", "NAME", "area", "Commercial", "HighdensityResidential", 
                     "Industrial", "Institutional", "LowdensityResidential", 
                     "ResidentialCommercial", "Percentage.sum")
landuse <- landuse[, columns_to_keep]

landuse_dff <- landuse[landuse$Percentage.sum>80,]
landuse_isna <- !is.na(landuse_dff$Percentage.sum)
landuse <- landuse_dff[landuse_isna,]

## Remove spatial class for joining
blockgroup_atlanta_with_geometry <- blockgroup_atlanta
blockgroup_atlanta <- st_set_geometry(blockgroup_atlanta, NULL)
landuse_df <- st_set_geometry(landuse, NULL)

## Join data
blockgroup_joined <- blockgroup_atlanta %>%
  left_join(landuse_df %>% select(-NAME), by = "GEOID")

## Filter out unmatched GEOIDs between blockgroup_atlanta and landuse to remove NaN values
missing_geoids <- setdiff(blockgroup_atlanta$GEOID, landuse_df$GEOID)

blockgroup_joined <- blockgroup_joined %>%
  filter(!GEOID %in% missing_geoids) %>%
  left_join(blockgroup_atlanta_with_geometry %>% select(GEOID, geometry), by = "GEOID") %>%
  st_as_sf()




# 5. Police Station Data

## Calculate the distance from each block group centroid to nearest police station
stations <- read_csv("Data/police_stations.csv")
stations_sf <- st_as_sf(stations, coords = c("x", "y"), crs = 3857) %>%
  st_transform(crs = st_crs(blockgroup_joined))

block_centroids <- blockgroup_joined %>% st_centroid()
block_centroids <- st_transform(block_centroids, st_crs(stations_sf))

data_final <- blockgroup_joined %>%
  mutate(min_station_dist = apply(st_distance(block_centroids, stations_sf), 1, min))


# 6. Final Data Cleaning and Save

## Keep Columns we will use for regression (it can be changed later)
## Keep GEOID and Name to match for index
columns_to_keep <- c("GEOID", "NAME", "tot_popE", "pop_den", "white_ratio", "black_ratio", "other_ratio", "median_incomeE", "less_than_hs_ratio", 
                     "Commercial", "HighdensityResidential", "Industrial", "Institutional", "LowdensityResidential", "ResidentialCommercial", "min_station_dist",
                     "geometry", "violent_count", "nonviolent_count")
data_final <- data_final[, columns_to_keep]

## Remove Null Rows
data_final <- data_final %>%
  drop_na()

## Export final dataset to GeoJSON
st_write(data_final, "data_final.geojson", driver = "GeoJSON")

## Remove geometry
data_final_csv <- data_final %>%
  st_drop_geometry()

## Export to CSV
write.csv(data_final_csv, "data_final.csv", row.names = FALSE)

