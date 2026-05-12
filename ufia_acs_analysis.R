# Associations between urban tree damage and neighborhood demography in 28 American cities 
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
#library(units)
library(here)
library(tidyr)
library(purrr)
library(readr)
library(lme4)
library(lmerTest)
#library(ggsignif)
library(broom)
library(broom.mixed)
library(performance)
#rm(list=ls())

wd <- here::here()
setwd(file.path(wd)) #getwd()



### load the file for the analysis 
csv_out_path <- file.path(here::here(),"out")
ufia_acs <- read_csv( file.path(csv_out_path, "ufia_acs_for_analysis_260511.csv")) %>% 
  mutate(plot_perc_damaged = plot_prop_damaged * 100,
         plot_perc_alive = plot_prop_alive * 100,
         plot_perc_dead = 100 - plot_perc_alive,
         plot_planted_perc = plot_planted * 100,
         plot_street_tree_perc = plot_street_tree * 100) %>% 
  filter(!is.na(city)) #remove the plots that did not connect to census data


## summary statistics for paper ####################################################################
#number of trees censused per city
ufia_acs %>% 
  group_by(city) %>% 
  summarize(n = sum(plot_n_trees),
            n_plots = n()) %>% 
  arrange(n_plots) %>% 
  print(n = 50)

#number of trees 
ufia_acs %>% 
  summarize(n_trees = sum(plot_n_trees))
  
#number of plots
nrow(ufia_acs)
  
#number of "cities"
ufia_acs %>% 
  select(city) %>% 
  distinct() %>% 
  summarize(n = n())


## data visualization and analysis #################################################################
### coarse map of study sites #################################
ufia_acs %>% 
  mutate(lat_1 = round(LAT, 1),
         lon_1 = round(LON, 1)) %>% 
  ggplot(aes(x = lon_1, y = lat_1)) + geom_point()




### Fig 1: comparison of tree damage and LAI with poverty and whiteness ############################

  ## Panel A: the combined effect of both poverty and whiteness on tree damage
    m1 <- lmer(plot_perc_damaged ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
    m1_summary <- summary(m1)
    print(m1_summary)
    check_collinearity(m1)
    ranef(m1)

    pov_df <- data.frame(estimate_c_perc_poverty = 0:round(max(ufia_acs$estimate_c_perc_poverty), 0))
    whi_df <- data.frame(estimate_c_perc_white = 0:100)
    
    pred_df2 <- expand_grid(pov_df, whi_df) %>% 
      mutate(city = "Washington, DC")
    
    pred_df3 <-  pred_df2 %>% 
      mutate(pred = predict(m1, newdata = pred_df2) - 0.5555428 ) #account for the effect of Washington DC 
    
    fig1_a <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
          geom_tile() + scale_fill_viridis_c(name = "damage (%)", option = "inferno") + xlab("poverty (%)") + ylab("White (%)") +
            ggthemes::theme_few() + theme(legend.position = "bottom") +
          geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                            inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")

    
    #extract some statistics about what differences this would result in at different percentiles of race and poverty
      study_wide_mean_poverty <- mean(ufia_acs$estimate_c_perc_poverty)
      study_wide_mean_white <- mean(ufia_acs$estimate_c_perc_white)
      study_wide_poverty1 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.01)
      study_wide_poverty99 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.99)
      study_wide_white1 <- quantile(ufia_acs$estimate_c_perc_white, 0.01)
      study_wide_white99 <- quantile(ufia_acs$estimate_c_perc_white, 0.99)
      
      min(ufia_acs$estimate_c_perc_white)
      min(ufia_acs$estimate_c_perc_poverty)
      max(ufia_acs$estimate_c_perc_poverty)
      max(ufia_acs$estimate_c_perc_white)
      
      
      pred_df <- data.frame(estimate_c_perc_poverty = c(study_wide_poverty1, study_wide_poverty99,
                                                        study_wide_mean_poverty, study_wide_mean_poverty,
                                                        study_wide_poverty99, study_wide_poverty1),
                            estimate_c_perc_white = c(study_wide_mean_white, study_wide_mean_white,
                                                      study_wide_white1, study_wide_white99,
                                                      study_wide_white1, study_wide_white99)) %>% 
        mutate(city = "Washington, DC")
     
      ranef(m1) # What is the random effect for the city-level intercept
      
      preds <- predict(m1, newdata = pred_df) - 0.5555428 #account for the effect of the city 
      pred_df %>% mutate(preds = preds)
      
      
    
    ## Panel B: the combined effect of both poverty and whiteness on LAI
    m1_b <- lmer(plot_mean_LAI ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
    m1_bsummary <- summary(m1_b)
    print(m1_bsummary)
    check_collinearity(m1_b)
    ranef(m1_b)
    
    pov_df <- data.frame(estimate_c_perc_poverty = 0:round(max(ufia_acs$estimate_c_perc_poverty), 0))
    whi_df <- data.frame(estimate_c_perc_white = 0:100)
    
    pred_df2 <- expand_grid(pov_df, whi_df) %>% 
      mutate(city = "Washington, DC")
    
    pred_df3 <-  pred_df2 %>% 
      mutate(pred = predict(m1_b, newdata = pred_df2) - 0.5555428 ) #account for the effect of Washington DC 
    
    fig1_b <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
      geom_tile() + scale_fill_viridis_c(name = expression("mean LAI (m"^2*"/m"^2*")"), 
                                         option = "viridis") + xlab("poverty (%)") + ylab("White (%)") +
      ggthemes::theme_few() + theme(legend.position = "bottom") +
      geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                 inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")
      
  
        fig1 <- cowplot::plot_grid(fig1_a, fig1_b, ncol = 2, labels = c("A", "B"))
        fig1
      
        #ggsave(fig1, filename = "fig1_260511.jpeg", dpi = 300, width = 8, height = 5, units = "in")
      

        #LAI: extract some statistics about what differences this would result in at  percentiles of race and poverty
                study_wide_mean_poverty <- mean(ufia_acs$estimate_c_perc_poverty)
                study_wide_mean_white <- mean(ufia_acs$estimate_c_perc_white)
                study_wide_poverty1 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.01)
                study_wide_poverty99 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.99)
                study_wide_white1 <- quantile(ufia_acs$estimate_c_perc_white, 0.01)
                study_wide_white99 <- quantile(ufia_acs$estimate_c_perc_white, 0.99)
                
                min(ufia_acs$estimate_c_perc_white)
                min(ufia_acs$estimate_c_perc_poverty)
                max(ufia_acs$estimate_c_perc_poverty)
                max(ufia_acs$estimate_c_perc_white)
        
        pred_df <- data.frame(estimate_c_perc_poverty = c(study_wide_poverty1, study_wide_poverty99,
                                                          study_wide_mean_poverty, study_wide_mean_poverty,
                                                          study_wide_poverty99, study_wide_poverty1),
                              estimate_c_perc_white = c(study_wide_mean_white, study_wide_mean_white,
                                                        study_wide_white1, study_wide_white99,
                                                        study_wide_white1, study_wide_white99)) %>% 
          mutate(city = "Washington, DC")
        
        ranef(m1_b) # What is the random effect for the city-level intercept
        
        preds <- predict(m1_b, newdata = pred_df) - 0.59060457 #account for the effect of the city 
        pred_df %>% mutate(preds = preds)
    
### Table 1: model results for each variable #######################################################
ufia_response_vars <- c(#"plot_spp_richness", #"plot_n_trees",
                        "plot_planted", #"plot_street_tree", #"people_per_ha",
                        #"plot_leaf_area", #"plot_BA", "plot_crown_diam", #"plot_compensatory_value",
                        "plot_mean_LAI", 
                        "plot_perc_damaged",  "plot_perc_alive" ) #"plot_mean_no_foliage", "plot_mean_dieback"
  
  
  # what to include in the model
  predictors <- "estimate_c_perc_poverty +  estimate_c_perc_white + (1|city)"
  
  # Fit models for each response variable
  coefficients_df <- ufia_response_vars %>%
    set_names() %>%  # keeps names for .id
    map(~ lmer(as.formula(paste(.x, "~", predictors)), data = ufia_acs)) %>%
    map_dfr(broom.mixed::tidy, .id = "response_variable") %>% 
    filter(effect == "fixed", term != "(Intercept)") %>%
    select(response_variable, term, estimate, std.error, p.value)
  print(coefficients_df)
  
  coefficients_df <- ufia_response_vars %>%
    set_names() %>%
    map(~ lmer(as.formula(paste(.x, "~", predictors)), data = ufia_acs)) %>%
    map_dfr(broom.mixed::tidy, conf.int = TRUE, .id = "response_variable") %>%
    filter(effect == "fixed", term != "(Intercept)") %>%
    select(response_variable, term, estimate, std.error, conf.low, conf.high, p.value)
  
  # ggplot(coefficients_df, aes(x = response_variable, y = estimate, ymin = conf.low, ymax = conf.high)) +
  #   geom_pointrange() + facet_wrap(~term)+
  #   theme_bw()
  
  # To get a wide format (coefficients as columns)
  coefficients_wide <- coefficients_df %>%
    mutate(estimate_conf = paste0(round(estimate, 4), " (", round(conf.low, 4), " - ", round(conf.high, 4), " )", round(p.value, 4))) %>% 
    select(response_variable, term, estimate_conf) %>%
    tidyr::pivot_wider(names_from = term, values_from = estimate_conf)
  
  print(coefficients_wide)
  
  coefficients_wide <- coefficients_df %>%
    tidyr::pivot_wider(
      names_from = term,
      values_from = c(estimate, std.error, p.value),
      names_glue = "{term}_{.value}"
    )
  
  print(coefficients_wide)
  
  # library(knitr)
  # library(kableExtra)
  # 
  # Format the table with rounded values
  coefficients_df %>%
    mutate(
      estimate = round(estimate, 4),
      std.error = round(std.error, 4),
      p.value = format.pval(p.value, digits = 4)
    ) %>%
    kable() %>%
    kable_styling()
  
  
  mi <- lmer(plot_perc_alive ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
  mi_summary <- summary(mi)   #
  print(mi_summary)
  poverty_slope_i <- round(mi_summary$coefficients[2], 3)
  
  
  
  
### SI 1: correlation between poverty and race/ethnicity ###########################################
  ufia_acs %>% 
    ggplot(aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white)) + geom_point()+ theme_bw() +
    xlab("neighborhood poverty (%)") + ylab("neighborhood Whiteness (%)")
  
  cor(ufia_acs$estimate_c_perc_poverty, ufia_acs$estimate_c_perc_white, method = "pearson")



  
### SI 2: comparison of plot trees alive, tree species richness, planted trees, and street trees with poverty and whiteness ################
study_wide_max_poverty <- max(ufia_acs$estimate_c_perc_poverty)
    
### panel A, B: city level relationships between damage and poverty and whiteness
  ufia_acs_focal <- ufia_acs %>% 
    mutate(focal_var = plot_perc_damaged )
  focal_y_lab <- expression("damage (%)")
  m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
  m_si2_focal_summary <- summary(m_si2_focal)
  check_collinearity(m_si2_focal)
  print(m_si2_focal_summary)
  
  poverty_p <- paste("p =", round( m_si2_focal_summary$coefficients[2,5], 3))
  white_p <- paste("p =", round( m_si2_focal_summary$coefficients[3,5], 3))
  
  ### poverty 
      ## individual city level
      pred_grid_city_pov <- ufia_acs %>% 
        group_by(city) %>% 
        summarise(estimate_c_perc_poverty = 0:study_wide_max_poverty, .groups = "drop") %>% 
        mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
      pred_focal_pov <- pred_grid_city_pov %>% 
        mutate(fit = predict(m_si2_focal, pred_grid_city_pov))
      
      #overall effect  
      effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_a_df <- as.data.frame(effect_var_a)
      
      #create plot    
      fig_si2_pov_a <- 
        ggplot() +
        geom_line(data = pred_focal_pov, aes(x = estimate_c_perc_poverty, y = fit, group = city), alpha = 0.2) +
        
        geom_line(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit), color = "blue") +
        geom_ribbon(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit, 
                                                ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
        annotate("text", label = poverty_p, x = 15, y = 80)+
        ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))

  
  ### whiteness
      ## individual city level
      pred_grid_city_white <- ufia_acs %>% 
        group_by(city) %>% 
        summarise(estimate_c_perc_white = 0:100, .groups = "drop") %>% 
        mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
      pred_focal_white <- pred_grid_city_white %>% 
        mutate(fit = predict(m_si2_focal, pred_grid_city_white))
      
      #overall effect  
      effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_b_df <- as.data.frame(effect_var_b)
      
      #create plot    
      fig_si2_white_b <- 
        ggplot() +
        geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit, group = city), alpha = 0.2) +
        
        geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit), color = "blue") +
        geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit, 
                                                ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
        annotate("text", label = white_p, x = 85, y = 63)+
        
        ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
      
      
### panel C, D: city level relationships between LAI and poverty and whiteness
      ufia_acs_focal <- ufia_acs %>% 
        mutate(focal_var = plot_mean_LAI )
      focal_y_lab <- expression("mean LAI (m"^2*"/m"^2*")")
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      poverty_p <- paste("p =", round( m_si2_focal_summary$coefficients[2,5], 3))
      white_p <- paste("p =", round( m_si2_focal_summary$coefficients[3,5], 3))
      
      
      ### poverty 
        ## individual city level
        pred_grid_city_pov <- ufia_acs %>% 
          group_by(city) %>% 
          summarise(estimate_c_perc_poverty = 0:study_wide_max_poverty, .groups = "drop") %>% 
          mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
        pred_focal_pov <- pred_grid_city_pov %>% 
          mutate(fit = predict(m_si2_focal, pred_grid_city_pov))
        
        #overall effect  
        effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
        effect_var_a_df <- as.data.frame(effect_var_a)
        
        #create plot    
        fig_si2_pov_c <- 
          ggplot() +
          geom_line(data = pred_focal_pov, aes(x = estimate_c_perc_poverty, y = fit, group = city), alpha = 0.2) +
          
          geom_line(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit), color = "blue") +
          geom_ribbon(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit, 
                                                  ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
          annotate("text", label = poverty_p, x = 15, y = 7.4)+
          ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
        
      
      ### whiteness
        ## individual city level
        pred_grid_city_white <- ufia_acs %>% 
          group_by(city) %>% 
          summarise(estimate_c_perc_white = 0:100, .groups = "drop") %>% 
          mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
        pred_focal_white <- pred_grid_city_white %>% 
          mutate(fit = predict(m_si2_focal, pred_grid_city_white))
        
        #overall effect  
        effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
        effect_var_b_df <- as.data.frame(effect_var_b)
        
        #create plot    
        fig_si2_white_d <- 
          ggplot() +
          geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit, group = city), alpha = 0.2) +
          
          geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit), color = "blue") +
          geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit, 
                                                  ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
          annotate("text", label = white_p, x = 83, y = 7.4)+
          ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
        
        
### panel E, F: city level relationships between dead trees and poverty and whiteness
      ufia_acs_focal <- ufia_acs %>% 
        mutate(focal_var = plot_perc_dead  )
      focal_y_lab <- expression("dead (%)")
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      poverty_p <- paste("p =", round( m_si2_focal_summary$coefficients[2,5], 3))
      white_p <- paste("p =", round( m_si2_focal_summary$coefficients[3,5], 3))
      
      
      ### poverty 
        ## individual city level
        pred_grid_city_pov <- ufia_acs %>% 
          group_by(city) %>% 
          summarise(estimate_c_perc_poverty = 0:study_wide_max_poverty, .groups = "drop") %>% 
          mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
        pred_focal_pov <- pred_grid_city_pov %>% 
          mutate(fit = predict(m_si2_focal, pred_grid_city_pov))
        
        #overall effect  
        effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
        effect_var_a_df <- as.data.frame(effect_var_a)
        
        #create plot    
        fig_si2_pov_e <- 
          ggplot() +
          geom_line(data = pred_focal_pov, aes(x = estimate_c_perc_poverty, y = fit, group = city), alpha = 0.2) +
          
          geom_line(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit), color = "blue") +
          geom_ribbon(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit, 
                                                  ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
          annotate("text", label = poverty_p, x = 15, y = 18) +
          ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
        
      
      ### whiteness
        ## individual city level
        pred_grid_city_white <- ufia_acs %>% 
          group_by(city) %>% 
          summarise(estimate_c_perc_white = 0:100, .groups = "drop") %>% 
          mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
        pred_focal_white <- pred_grid_city_white %>% 
          mutate(fit = predict(m_si2_focal, pred_grid_city_white))
        
        #overall effect  
        effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
        effect_var_b_df <- as.data.frame(effect_var_b)
        
        #create plot    
        fig_si2_white_f <- 
          ggplot() +
          geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit, group = city), alpha = 0.2) +
          
          geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit), color = "blue") +
          geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit, 
                                                  ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
          annotate("text", label = white_p, x = 85, y = 19)+
          ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
        #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
        
### save fig si2 
    fig_si2 <- cowplot::plot_grid(fig_si2_pov_a, fig_si2_white_b,
                                  fig_si2_pov_c, fig_si2_white_d,
                                  fig_si2_pov_e, fig_si2_white_f,
                                  ncol = 2, labels = c("A", "B", "C", "D", "E", "F"))
    
    ggsave(fig_si2, filename = "fig_SI_2_260512.jpeg",  width = 7, height = 10, units = "in", dpi = 300) #may need some resizing
    
        