# Associations between urban tree damage and neighborhood demography in 20 American cities 
# Maya Mangala Munuma, Alexander Young, Eli Robinson, Daniel S.W. Katz

# This script contains data analysis and visualization
# It uses data that were assembled in the "UFIA_data_assembly.R" script in this same repo
# This version of the analysis relies on buffering the public location of the FIA plot
# previous versions of this project are available at
# https://github.com/bearsofthemoss/tree_census/tree/main/code


### set up work environment
#load all required packages
library(ggplot2)
library(dplyr)
library(units)
library(here)
library(tidyr)
library(purrr)
library(readr)
library(lme4)
library(lmerTest)
library(ggsignif)
#rm(list=ls())

wd <- here::here()
setwd(file.path(wd)) #getwd()



### load the file for the analysis 
csv_out_path <- file.path(here::here(),"out")
ufia_acs <- read_csv( file.path(csv_out_path, "ufia_acs_for_analysis_260102.csv")) %>% 
  mutate(plot_perc_damaged = plot_prop_damaged * 100) %>% 
  filter(!is.na(city)) #remove the plots that did not connect to census data


## data visualization and analysis #################################################################


### Fig 1: comparison of tree damage with poverty and whiteness ############################
  m1 <- lmer(plot_mean_LAI ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
  m1_summary <- summary(m1)
  print(m1_summary)
  
  
  m1_slope <- round(m1_summary$coefficients[2,1], 2)
  m1_intercept <- round(m1_summary$coefficients[1,1], 2)
  m1_slope_p <- round(m1_summary$coefficients[2,5], 2)

  ### PICK BACK UP HERE: FIGURE OUT HOW TO CREATE CI'S FOR MODEL AND PLOT THOSE ACROSS THE RAW DATA

### Panel a: damage and poverty


#panel_a <- 
  ggplot(ufia_acs, aes(x = estimate_c_perc_poverty, y = plot_perc_damaged, group = city)) + 
  geom_point(color = "gray", alpha = 0.3) +
  geom_smooth(method = "lm", se = FALSE, color = "gray70") +
  theme_bw() + 
  theme(legend.position = "none") + 
  xlab("poverty (%)") + ylab("damaged trees (%)")
  


### SI X: correlation between poverty and race/ethnicity ###########################################
  ufia_acs %>% 
    ggplot(aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white)) + geom_point()+ theme_bw() +
    xlab("poverty (%)") + ylab("white (percent)")


#creating a figure for a single comparison

m1 <- lmer(plot_perc_damaged ~ estimate_c_perc_poverty + estimate_c_perc_white + (1|city), data = ufia_acs)
m1_summary <- summary(m1)
m1_slope <- round(m1_summary$coefficients[2,1], 2)
m1_intercept <- round(m1_summary$coefficients[1,1], 2)
m1_slope_p <- round(m1_summary$coefficients[2,5], 2)

#sink("lm.txt")
print(m1_summary)
#sink()
str(m1_summary)
confint(m1) #confidence intervals for each coefficient
ranef(m1) #random effects for each city

plot(m1)
# ee <- effects::Effect(c("estimate_c_perc_poverty","estimate_c_perc_white"),m1) 
# plot(ee)

sjPlot::plot_model(m1, type = "pred") +
    geom_hex(aes(x = estimate_c_perc_poverty, y = plot_perc_damaged),  data = ufia_acs,
        bins = 10, alpha = 0.2, show.legend = TRUE) +  
    scale_fill_viridis_c(lim = c(1, 1000), na.value = NA, trans = "log10") + 
    #annotate("text", x = 20, y = 95, label = paste0("y = ", m1_slope, " * x + ", m1_intercept, ", p = ", m1_slope_p)) +
    guides(fill = guide_legend(title="observations (n)")) +
   # theme(legend.position = "inside", legend.position.inside = c(0.8, 0.2)) + 
    ggthemes::theme_few() + ggtitle("") + xlab("poverty (%)") + ylab("damaged trees (%)")



qqnorm(resid(m1))
hist(resid(m1))

poverty_list <- data.frame(estimate_c_perc_poverty = 1:100)
white_list <- data.frame(estimate_c_perc_white = 1:100)
city_list = data.frame(city = sort(unique(ufia_acs$city)))
            
            
pred_df <- expand_grid(city_list, poverty_list, white_list)
pred_df2 <- pred_df %>% 
          mutate( #pred_val = merTools::predictInterval(merMod = m1, newdata = ., level = 0.95) #this went from -50 - 120
                  pred_val = predict(m1, .)
                 )






# # categorizing cities by location in the country
# NE_cities <- c("BaltimoreMD2022Curr", "BurlingtonVT2022Curr",  "ChicagoIL2022Curr", "ClevelandOH2022Curr", 
#                "DesMoinesIA2022Curr",
#                "MadisonWI2022Curr",  "MilwaukeeWI2022Curr",  "MinneapolMN2022Curr",
#                "PittsburghPA2022Curr", "PortlandME2022Curr", 
#                "ProvidenceRI2022Curr", "RochesterNY2022Curr", "TrentonNJ2022Curr", "WashingtonDC2022Curr" )
# SW_cities <- c("AustinTX", "HoustonTX", "SanAntonioTX" )
# 
# NE_cities_not_eval <- c("BaltimoreMD", "BurlingtonVT",  "ChicagoIL", "ClevelandOH",   "DesMoinesIA",  
#                         "MadisonWI",  "MilwaukeeWI",  "MinneapolMN","PittsburghPA", "PortlandME", 
#                         "ProvidenceRI", "RochesterNY", "TrentonNJ", "WashingtonDC" )