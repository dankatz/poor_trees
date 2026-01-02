# Associations between urban tree damage and neighborhood demography in 20 American cities 
# Maya Mangala Munuma, Alexander Young, Eli Robinson, Daniel S.W. Katz

# Data assembly for UFIA and census 
# This version of the analysis relies on buffering the public location of the FIA plot
# previous versions of this script are available at
# https://github.com/bearsofthemoss/tree_census/tree/main/code

### set up work environment
#load all required packages
library(tidycensus)
library(ggplot2)
library(dplyr)
library(sf)
library(units)
#library(plyr)
library(stringr)
library(scales)
library(here)
library(stringdist)
library(tidyr)
library(purrr)
library(readr)
library(lme4)
library(lmerTest)
library(ggsignif)
#rm(list=ls())

wd <- here::here()
setwd(file.path(wd)) #getwd()

options(scipen=999)
options(tigris_use_cache = TRUE)


# categorizing cities by location in the country
NE_cities <- c("BaltimoreMD2022Curr", "BurlingtonVT2022Curr",  "ChicagoIL2022Curr", "ClevelandOH2022Curr", 
               "DesMoinesIA2022Curr",
               "MadisonWI2022Curr",  "MilwaukeeWI2022Curr",  "MinneapolMN2022Curr",
               "PittsburghPA2022Curr", "PortlandME2022Curr", 
               "ProvidenceRI2022Curr", "RochesterNY2022Curr", "TrentonNJ2022Curr", "WashingtonDC2022Curr" )
SW_cities <- c("AustinTX", "HoustonTX", "SanAntonioTX" )

NE_cities_not_eval <- c("BaltimoreMD", "BurlingtonVT",  "ChicagoIL", "ClevelandOH",   "DesMoinesIA",  
                        "MadisonWI",  "MilwaukeeWI",  "MinneapolMN","PittsburghPA", "PortlandME", 
                        "ProvidenceRI", "RochesterNY", "TrentonNJ", "WashingtonDC" )

#including additional cities that were available on 11/23/2024
evals <- c("AustinTX2022Curr", "BaltimoreMD2022Curr", "BurlingtonVT2022Curr",  "ChicagoIL2022Curr",
           "ClevelandOH2022Curr",  "DesMoinesIA2022Curr", "HoustonTX2022Curr", "KansasCityMO2022Curr",
           "MadisonWI2022Curr",  "MilwaukeeWI2022Curr",  "MinneapolMN2022Curr",  "PittsburghPA2022Curr",
           "PortlandME2022Curr", "PortlandOR2022Curr", "ProvidenceRI2022Curr", "RochesterNY2022Curr", 
           "SanAntonioTX2022Curr", "SanDiegoCA2022Curr", "SpringfielMO2020Curr", "StLouisMO2022Curr",    
           "TrentonNJ2022Curr", "WashingtonDC2022Curr" )

### adding the Urban FIA data from datamart and creating derived variables #####################################
    # this version of UFIA was downloaded Nov 2024
    psca <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_PLOT_STRAT_CALC_ASSGN.csv"))
    psc <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/POP_STRATUM_CALC.csv"))
    plt <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_PLOT.csv")) #summary(plt)
    ref_plot_status <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/REF_PLOT_STATUS.csv"))
    mtre <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_MOTHER_TREE.csv"))
    indiv_tree <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_TREE.csv"))
    subp <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_SUBPLOT.csv"))
    cnd <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_COND.csv"))
    spcnd <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/ID_SUBP_COND.csv"))
    ref_species <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/REF_SPECIES.csv"))
    ref_species_group <- read.csv(file.path("data","FIADB_URBAN_ENTIRE_CSV/REF_SPECIES_GROUP.csv"))
    plt$PLOT_STATUS_CD_LAB <- ref_plot_status$ABBR[match(plt$PLOT_STATUS_CD, ref_plot_status$VALUE)]
    
  ## summarizing tree data at the plot level 
   #creating a damage category from the ID_TREE file for adding to mtre before summarizing
      indiv_tree_damage <- indiv_tree %>% 
        mutate(damaged = case_when(DMG_ROOT_STEM_GIRDLING == 1 ~ 1,
                                   DMG_TRUNK_BARK_INCLUSION == 1 ~ 1,
                                   DMG_EXCESS_MULCH == 1 ~ 1,
                                   DMG_TOPPING_PRUNING == 1 ~ 1,
                                   DMG_SIDEWALK_ROOT_CONFLICT == 1 ~ 1,
                                   DMG_OVERHEAD_WIRES == 1 ~ 1,
                                   DMG_IMPROPER_PLANTING== 1 ~ 1, 
                                   .default = 0)) %>% 
        group_by(PLOTID, TREE) %>% 
        dplyr:: summarise(damaged = max(damaged, na.rm = TRUE))

  ## summarizing mtre tree data at the plot level 
    plot_tree_summary <- 
      left_join(mtre, indiv_tree_damage) %>% #adding in the tree damage from the indiv_tree (ID_TREE) dataset
      group_by(PLT_CN) %>% 
      filter(SUBP == 1) %>%  #restricting to non-sapling trees (DBH > 5 in)
      filter(STATUSCD == 1 | STATUSCD == 2) %>% #removing trees that weren't measured due to no longer being in the sample
      #STATUSCD 0 == tree is not in the remeasured plot, STATUSCD 3 == cut and utilized, STATUSCD 4 == removed
      mutate(trees_alive = case_when(STATUSCD == 2 ~ 0, #STATUSCD 1 == live tree, STATUSCD 2 == dead tree
                                     STATUSCD == 1 ~ 1),
             trees_planted = case_when(IS_PLANTED == 1 ~ 1, #1 == planted
                                       IS_PLANTED == 2 ~ 0, #2 == natural origin
                                       IS_PLANTED == 3 ~ 0)) %>% #3 == "not sure"; we are lumping unsure with natural origin
      dplyr::summarize(plot_BA = sum(BASAL_AREA, na.rm = TRUE),
                       # The 'TPA_UNADJ' column value (trees per acre) 
                       # might be convenient as the units are sq ft per acre 
                       plot_n_trees = n(),
                       plot_prop_alive = mean(trees_alive, na.rm = TRUE),
                       plot_mean_dieback = mean(CROWN_DIEBACK_CD, na.rm = TRUE),
                       plot_prop_damaged = mean(damaged, na.rm = TRUE),
                       plot_mean_no_foliage = mean(FOLIAGE_ABSENT, na.rm = TRUE),
                       plot_crown_diam = mean(CROWN_DIA_90, na.rm = TRUE),
                       plot_street_tree = mean(IS_STREET_TREE, na.rm = TRUE),
                       plot_planted = mean(trees_planted, na.rm = TRUE),
                       plot_leaf_area = sum(LEAF_AREA_ITREE, na.rm = TRUE),
                       plot_mean_LAI_raw = mean(LEAF_AREA_INDEX_ITREE, na.rm = TRUE), #this is an average of LAI for all trees at a plot,
                          #it is not the average LAI of canopy within the plot though
                       plot_mean_LAI = weighted.mean(x = LEAF_AREA_INDEX_ITREE, w = CROWN_GROUND_AREA_ITREE, na.rm = TRUE), #this is 
                          #a mean of LAI at the plot weighted by tree canopy area
                       plot_compensatory_value = mean(COMPENSATORY_VALUE_ITREE, na.rm = TRUE),
                       plot_spp_richness = n_distinct(SPCD))
    
    

    
### start the city loop for extracting US Census (ACS) data for each UFIA plot ##############
pcv_out <- list() 
for(i in c(1:length(evals))){   # to run all cities. 
    #for(i in c(1:2)){   #
  
  ### For each city prepare UFIA data ======================================================
  
  # Specify the EVALID (evaluation ID)
  city_choose <- evals[i]  #city_choose <- evals[1]
  
  # subset the Pop stratum Calc (psc) table for a single evaluation
  ## The rows provide the NLCD-based area expansion factors, each plot is assigned to one.
  Stratum <- psc[psc$EVALID == city_choose ,]
  Stratum$STRATUM_CN <- as.factor(Stratum$CN)
  
  ##  Select only the city estimation unit for 3 cities that have multiple estimation units
  if(city_choose == "StLouis2021Curr"){
    Stratum <- Stratum[Stratum$ESTN_UNIT_NAME=="City of St. Louis, MO",] }
  
  if(city_choose == "KansasCity2021Curr"){
    Stratum <- Stratum[Stratum$ESTN_UNIT_NAME=="City of Kansas City, MO",] }
  
  if(city_choose == "SanAntonio2021Curr"){
    Stratum <- Stratum[Stratum$ESTN_UNIT_NAME=="City of San Antonio, TX",] }
  
  # parse Stratum dataframe to identify the city
  parsed_city_name <- sub(" [0-9]{4}.*", "", unique(Stratum$EVAL_NAME))
  
  ## now find which plots are associated with each strata using the psca table.
  ##  Plot stratum calc assign- has plot (PLT_CN) and Stratum (PSC_CN) pairings.
  Plot <- psca[psca$PSC_CN %in% Stratum$CN ,]
  
  # bring the lat and lon over from the plot  (plt) table.
  Plot$LON <- plt$LON[match(Plot$PLT_CN, plt$CN)]
  Plot$LAT <- plt$LAT[match(Plot$PLT_CN, plt$CN)]
  
  # include label for plot status
  Plot$PLOT_STATUS_CD_LAB <- plt$PLOT_STATUS_CD_LAB[match(Plot$PLT_CN, plt$CN)]
  
  # Add the eval ID to the Plot dataframe
  Plot$EVALID <- city_choose
  
  # Set CRS of Plot df to nad 83
  sf_plots <-  st_as_sf(Plot, coords = c("LON", "LAT"), crs = 4269) 
  
  # bring in state code for census API
  sf_plots$state <- plt$STATECD[match(sf_plots$PLT_CN, plt$CN)]
  
  
  ### get  relevant US census data and calculate summary for each buffered plot =======================
  ## download census data for the relevant state
  bg_data <- get_acs(geography = "block group", 
                     variables = c(
                       #c_total_pop = "B03002_001",
                       #c_pop = "B01003_001", #total population
                       
                       #Race/ethnicity # see descriptions here: https://censusreporter.org/topics/race-hispanic/
                       c_total_pop = "B03002_001", c_white = "B03002_003", c_black = "B03002_004", c_latinx = "B03002_012",
      
                       #poverty: for more detail on how poverty is assessed, look here: https://censusreporter.org/topics/poverty/
                       c_poverty_1 = "C17002_001", #total population of poverty table #https://censusreporter.org/tables/C17002/
                       c_poverty_2 = "C17002_002", #population less than 50% of poverty line
                       c_poverty_3 = "C17002_003", #population 50-99% of poverty line
                       
                       #c_building_age = "B25035_001",
                       medincome = "B19013_001"), 
                    
                     state = unique(sf_plots$state), 
                     #summary_var = "B03002_001",
                     year = 2020, 
                     survey = "acs5", #use the 5-yr survey for greater reliability
                     output = "wide",
                     geometry = TRUE) %>% 
            #calculate derived variables for poverty and race/ethnicity
                    mutate(estimate_c_poverty = (c_poverty_2E + c_poverty_3E)/c_poverty_1E,
                           estimate_c_perc_poverty = estimate_c_poverty * 100) %>% 
                    mutate(estimate_c_white = c_whiteE/c_total_popE,
                           estimate_c_perc_white = estimate_c_white * 100)
  

  ## filter block groups to the relevant place
    #get the list of places in the focal state
    places <- get_acs(geography = "place", 
                      variables = c(medincome = "B19013_001"), 
                      state = unique(sf_plots$state), 
                      year = 2020,
                      survey = "acs5",
                      geometry = TRUE) 
    
    # Find the index of the best match
    best_match_index <- stringdist::amatch(parsed_city_name, places$NAME , maxDist = Inf)
    
    # Get the best matching row from the dataframe
    just_the_place <- places[best_match_index, ]
    
    # spatial subset of the state block groups by the city place boundary
    city_block_groups <- bg_data[just_the_place, ]
    
    # this catches an error with Kansas City and Houston
    if(just_the_place$NAME =="Warsaw city, Missouri"){
      just_the_place <- places[places$NAME=="Kansas City city, Missouri" ,]
      city_block_groups <- bg_data[just_the_place, ]
    }
    
    # this catches an error with Houston 2021
    if(just_the_place$NAME =="Howe town, Texas"){
      just_the_place <- places[places$NAME=="Houston city, Texas" ,]
      city_block_groups <- bg[just_the_place, ]
    }
  
  
  ## buffer UFIA plots in focal city
    sf_plots_1km <- st_buffer(sf_plots, 1000) #plot(sf_plots_1km)
    sf_plots_1km_bbox <- st_bbox(sf_plots_1km) #for clipping block groups to plot area# plot(sf_plots_1km_bbox)
  
  ## clip block groups to buffered UFIA plots
    bg_data_city <- st_crop(bg_data, sf_plots_1km_bbox) %>% 
      mutate(bg_area_original = st_area(geometry),  #for use in weighted average of demographic variables 
             people_per_area = as.numeric(c_total_popE /bg_area_original))
    bg_data_clipped <- st_intersection(bg_data_city, sf_plots_1km) #takes ~10 sec to run
    
    #removing any block groups that no one lives in (poverty and race/income are irrelevant there) 
    bg_data_clipped <- bg_data_clipped %>%  
                        filter(c_total_popE > 0) %>%  #plot(bg_data_clipped[1])
                        filter(c_poverty_1E > 0) 
      
    bg_data_clipped_c <- bg_data_clipped %>% 
         filter(st_is_valid(.) == TRUE) #remove any potential errors in geometry
    
     print(paste("removed this many block groups with bad geometry:", nrow(bg_data_clipped) - nrow(bg_data_clipped_c), 
                 "out of a total of", nrow(bg_data_clipped)))
    
     #is census data missing for any areas?
       print(paste("number of NA values in poverty estimate:", sum(is.na(bg_data_clipped_c$estimate_c_perc_poverty))))
       print(paste("number of NA values in whiteness estimate:", sum(is.na(bg_data_clipped_c$estimate_c_perc_white))))
             
     
  ##calculate demographic response variable weighted means within 1 km of each plot
   bg_focal_city <- bg_data_clipped_c %>% 
        mutate(area_subset_m2 = st_area(geometry), #calculating the area of the polygon subset
               people_in_subset = people_per_area * as.numeric(area_subset_m2)) %>% #calculate the people within that polygon subset
        st_drop_geometry()   %>% 
        dplyr::group_by(PLT_CN, CN, EVALID) %>% 
        dplyr::summarize( estimate_c_perc_poverty = round(weighted.mean(x = estimate_c_perc_poverty, w = people_in_subset), 3),
                          estimate_c_perc_white = round(weighted.mean(x = estimate_c_perc_white, w = people_in_subset), 3),
                          people_per_ha = weighted.mean(x = people_per_area, w = people_in_subset)*10000) %>% 
        mutate(city = city_choose) %>% 
        ungroup()
  
   
   ## export results from city loop
   if(i == 1){bg_out <- bg_focal_city}else{bg_out <- bind_rows(bg_out, bg_focal_city)}
   
   print(paste("finished with city:", city_choose))
} #end city loop for census data

    sum(is.na(bg_out))
    
### combine UFIA and US Census data and export for analysis #######################################################  
    # bg_out 
    # plot_tree_summary
    ufia_acs <- left_join(plot_tree_summary, bg_out) 

  
    ### save the file used in the analysis 
    csv_out_path <- file.path(here::here(),"out")
    formatted_date <- format(Sys.Date(), "%y%m%d") # Gives "02-01-2025"

    write_csv(x = ufia_acs, file = file.path(csv_out_path, paste0("ufia_acs_for_analysis_", formatted_date, ".csv")))
