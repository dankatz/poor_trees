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
ufia_acs <- read_csv( file.path(csv_out_path, "ufia_acs_for_analysis_260102.csv")) %>% 
  mutate(plot_perc_damaged = plot_prop_damaged * 100,
         plot_perc_alive = plot_prop_alive * 100,
         plot_perc_dead = 100 - plot_perc_alive,
         plot_planted_perc = plot_planted * 100,
         plot_street_tree_perc = plot_street_tree * 100) %>% 
  filter(!is.na(city)) #remove the plots that did not connect to census data


## data visualization and analysis #################################################################


### Fig 1: comparison of tree damage with poverty and whiteness ############################
  m1 <- lmer(plot_perc_damaged ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
  m1_summary <- summary(m1)
  print(m1_summary)
  check_collinearity(m1)
  ranef(m1)
  
  ##  Panel A: the effect of poverty on tree damage
  # Calculate the effects 
  effect_var1 <- effects::effect("estimate_c_perc_poverty", m1, xlevels = 100)
  effect_var1_df <- as.data.frame(effect_var1)
  #effect_var1_df
  fig1_a <- 
    ggplot(effect_var1_df, aes(x = estimate_c_perc_poverty, y = fit)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      theme_bw() + ylab("damaged trees (%)") + xlab("neighborhood poverty (%)") +
    geom_rug(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = plot_perc_damaged), alpha = 0.05, position = position_jitter(height = 2))
  
  ##  Panel B: the effect of whiteness on tree damage
  # Calculate the effects 
  effect_var2 <- effects::effect("estimate_c_perc_white", m1, xlevels = 100)
  effect_var2_df <- as.data.frame(effect_var2)
  fig1_b <-
    ggplot(effect_var2_df, aes(x = estimate_c_perc_white, y = fit)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      theme_bw() + ylab("damaged trees (%)") + xlab("neighborhood whiteness (%)") +
    geom_rug(data = ufia_acs, aes(x = estimate_c_perc_white, y = plot_perc_damaged), alpha = 0.05, position = position_jitter(height = 2))
  
  ## Panel C: the combined effect of both poverty and whiteness
        pov_df <- data.frame(estimate_c_perc_poverty = 0:100)
        whi_df <- data.frame(estimate_c_perc_white = 0:100)
        
        pred_df2 <- expand_grid(pov_df, whi_df) %>% 
          mutate(city = "Baltimore, MD")
        
        pred_df3 <-  pred_df2 %>% 
          mutate(pred = predict(m1, newdata = pred_df2) + 4.1197) #account for the effect of Baltimore 
        
      
        fig1 <- cowplot::plot_grid(fig1_a, fig1_b, ncol = 1, labels = c("A", "B"))
        fig1
      
        ggsave(fig1, filename = "fig1_260108.jpeg", dpi = 300, width = 8, height = 4, units = "in")
      
        fig1_c <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
          geom_tile() + scale_fill_viridis_c(name = "damage (%)", option = "inferno") + xlab("poverty (%)") + ylab("white (%)") +
            ggthemes::theme_few() + theme(legend.position = "bottom") +
          geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                            inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")
        
        fig1_next <- cowplot::plot_grid(fig1, fig1_c, ncol = 2, labels = c(" ", "C"), rel_widths = c(1, 1.6))
  
  ## save file    
  ggsave(fig1_next, filename = "fig1_260108.jpeg", dpi = 300, width = 7, height = 5, units = "in")
  
  
  #extract some statistics about what differences this would result in
    study_wide_mean_poverty <- mean(ufia_acs$estimate_c_perc_poverty)
    study_wide_mean_white <- mean(ufia_acs$estimate_c_perc_white)
    min(ufia_acs$estimate_c_perc_white)
    min(ufia_acs$estimate_c_perc_poverty)
    max(ufia_acs$estimate_c_perc_poverty)
    
    pred_df <- data.frame(estimate_c_perc_poverty = c(study_wide_mean_poverty, study_wide_mean_poverty, 0, 84,0,84),
                          estimate_c_perc_white = c( 0, 100, study_wide_mean_white, study_wide_mean_white, 100, 0)) %>% 
              mutate(city = "Baltimore, MD")
    ranef(m1) # What is the random effect for the city-level intercept
    
    predict(m1, newdata = pred_df) + 4.1197 #account for the effect of Baltimore 

  
    
    
# ### Fig 2: comparison of plot tree mortality with poverty and whiteness ############################
#   m2 <- lmer(plot_perc_dead  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
#   m2_summary <- summary(m2)
#   print(m2_summary)
#   check_collinearity(m2)
#   
#   ##  Panel A: the effect of poverty on plot mean LAI
#   # Calculate the effects 
#   effect_var3 <- effects::effect("estimate_c_perc_poverty", m2, xlevels = 100)
#   effect_var3_df <- as.data.frame(effect_var3)
#   #effect_var1_df
#   fig2_a <- 
#     ggplot(effect_var3_df, aes(x = estimate_c_perc_poverty, y = fit)) +
#     geom_line() +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#     theme_bw() + ylab("dead trees (%)") + xlab("neighborhood poverty (%)") +
#     geom_rug(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = plot_mean_LAI), alpha = 0.05, position = position_jitter(height = 2))
#   
#   ##  Panel B: the effect of whiteness on tree LAI
#   # Calculate the effects 
#   effect_var4 <- effects::effect("estimate_c_perc_white", m2, xlevels = 100)
#   effect_var4_df <- as.data.frame(effect_var4)
#   fig2_b <-
#     ggplot(effect_var4_df, aes(x = estimate_c_perc_white, y = fit)) +
#     geom_line() +
#     geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
#     theme_bw() + ylab(expression("dead trees (%)")) + xlab("neighborhood whiteness (%)") +
#     geom_rug(data = ufia_acs, aes(x = estimate_c_perc_white, y = plot_mean_LAI), alpha = 0.05, position = position_jitter(height = 2)) 
#   
#   
#   fig2 <- cowplot::plot_grid(fig2_a, fig2_b, ncol = 2, labels = c("A", "B"))
#   fig2
#   #ggsave(fig2, filename = "fig2_260107.jpeg", dpi = 300, width = 7, height = 3.5, units = "in")
#   
#   #extract some statistics about what differences this would result in
#   study_wide_mean_poverty <- mean(ufia_acs$estimate_c_perc_poverty)
#   study_wide_mean_white <- mean(ufia_acs$estimate_c_perc_white)
#   min(ufia_acs$estimate_c_perc_white)
#   min(ufia_acs$estimate_c_perc_poverty)
#   max(ufia_acs$estimate_c_perc_poverty)
#   
#   pred_df <- data.frame(estimate_c_perc_poverty = c(study_wide_mean_poverty, study_wide_mean_poverty, 0, 84),
#                         estimate_c_perc_white = c( 0, 100, study_wide_mean_white, study_wide_mean_white)) %>% 
#     mutate(city = "Baltimore, MD")
#   ranef(m2) # What is the random effect for the city-level intercept
#   
#   predict(m2, newdata = pred_df) + 0.1073 #account for the effect of Baltimore 
#   
#   #a new figure including both poverty and whiteness
#       pov_df <- data.frame(estimate_c_perc_poverty = 0:100)
#       whi_df <- data.frame(estimate_c_perc_white = 0:100)
#       
#       pred_df2 <- expand_grid(pov_df, whi_df) %>% 
#         mutate(city = "Baltimore, MD")
#       
#       pred_df3 <-  pred_df2 %>% 
#         mutate(pred = predict(m2, newdata = pred_df2) - 1.749) #account for the effect of Baltimore #ranef(m2) 
#       
#       ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
#         geom_tile() + scale_fill_viridis_c(name = "dead trees (%)", option = "inferno") + xlab("poverty (%)") + ylab("white (%)") +
#         ggthemes::theme_few()
#   
#       #exploring a new figure including both poverty and whiteness
#         pov_df <- data.frame(estimate_c_perc_poverty = 0:100)
#         whi_df <- data.frame(estimate_c_perc_white = 0:100)
#         
#         pred_df2 <- expand_grid(pov_df, whi_df) %>% 
#           mutate(city = "Baltimore, MD")
#         
#         pred_df3 <-  pred_df2 %>% 
#           mutate(pred = predict(m2, newdata = pred_df2) + 0.1073) #account for the effect of Baltimore 
#         
#         ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
#           geom_tile() + scale_fill_viridis_c(name = "damage (%)", option = "inferno") + xlab("poverty (%)") + ylab("white (%)") +
#           ggthemes::theme_few()
      
### Table 1: model results for each variable #######################################################
ufia_response_vars <- c("plot_spp_richness", #"plot_n_trees",
                        "plot_planted", "plot_street_tree", #"people_per_ha",
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
  
  
  mi <- lmer(plot_perc_damaged ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
  mi_summary <- summary(mi)   #print(m1_summary)
  poverty_slope_i <- round(mi_summary$coefficients[2], 3)
  
  
  
  
### SI X: correlation between poverty and race/ethnicity ###########################################
  ufia_acs %>% 
    ggplot(aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white)) + geom_point()+ theme_bw() +
    xlab("neighborhood poverty (%)") + ylab("neighborhood whiteness (%)")
  
  cor(ufia_acs$estimate_c_perc_poverty, ufia_acs$estimate_c_perc_white, method = "pearson")



  
### SI 2: comparison of plot trees alive, tree species richness, planted trees, and street trees with poverty and whiteness ################
  ### panel A, B: LAI 
      ufia_acs_focal <- ufia_acs %>% 
          mutate(focal_var = plot_mean_LAI )
      focal_y_lab <- expression("mean LAI (m"^2*"/m"^2*")")
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      
      ## Poverty panel for focal variable
      effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_a_df <- as.data.frame(effect_var_a)
      fig_si2_pov_a <- 
        ggplot(effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))+
        coord_cartesian(ylim = c(5.5, 8))
      
      ## whiteness panel for focal variable
      effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_b_df <- as.data.frame(effect_var_b)
      fig_si2_whi_b <- 
        ggplot(effect_var_b_df, aes(x = estimate_c_perc_white, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2)) +
        coord_cartesian(ylim = c(5.5, 8))
      
  ### panel C, D: tree species richness 
      ufia_acs_focal <- ufia_acs %>% 
        mutate(focal_var = plot_spp_richness )
      focal_y_lab <- "species richness (n)"
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      
      ## Poverty panel for focal variable
      effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_a_df <- as.data.frame(effect_var_a)
      fig_si2_pov_c <- 
        ggplot(effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 0.5))
      
      ## whiteness panel for focal variable
      effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_b_df <- as.data.frame(effect_var_b)
      fig_si2_whi_d <- 
        ggplot(effect_var_b_df, aes(x = estimate_c_perc_white, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 0.5))
      
      
  ### panel E, F: percent of trees that are planted
      ufia_acs_focal <- ufia_acs %>% 
        mutate(focal_var = plot_planted_perc)
      focal_y_lab <- "planted trees (%)"
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      
      ## Poverty panel for focal variable
      effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_a_df <- as.data.frame(effect_var_a)
      fig_si2_pov_e <- 
        ggplot(effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
      
      ## whiteness panel for focal variable
      effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_b_df <- as.data.frame(effect_var_b)
      fig_si2_whi_f <- 
        ggplot(effect_var_b_df, aes(x = estimate_c_perc_white, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
      
      
  ### panel G, H: percent of trees that are street trees
      ufia_acs_focal <- ufia_acs %>% 
        mutate(focal_var = plot_street_tree_perc  )
      focal_y_lab <- "street trees (%)"
      m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
      m_si2_focal_summary <- summary(m_si2_focal)
      check_collinearity(m_si2_focal)
      print(m_si2_focal_summary)
      
      ## Poverty panel for focal variable
      effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_a_df <- as.data.frame(effect_var_a)
      fig_si2_pov_g <- 
        ggplot(effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
      
      ## whiteness panel for focal variable
      effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
      effect_var_b_df <- as.data.frame(effect_var_b)
      fig_si2_whi_h <- 
        ggplot(effect_var_b_df, aes(x = estimate_c_perc_white, y = fit)) +
        geom_line() +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
        theme_bw() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") +
        geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
      
  ### panel I, J: dead trees
    ufia_acs_focal <- ufia_acs %>% 
      mutate(focal_var =  plot_perc_dead  )
    focal_y_lab <- "dead trees (%)"
    m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
    m_si2_focal_summary <- summary(m_si2_focal)
    check_collinearity(m_si2_focal)
    print(m_si2_focal_summary)
    
    ## Poverty panel for focal variable
    effect_var_a <- effects::effect("estimate_c_perc_poverty", m_si2_focal, xlevels = 100)       # Calculate the effects 
    effect_var_a_df <- as.data.frame(effect_var_a)
    fig_si2_pov_i <- 
      ggplot(effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      theme_bw() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") +
      geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2)) +
      coord_cartesian(ylim = c(0,10))
    
    ## whiteness panel for focal variable
    effect_var_b <- effects::effect("estimate_c_perc_white", m_si2_focal, xlevels = 100)       # Calculate the effects 
    effect_var_b_df <- as.data.frame(effect_var_b)
    fig_si2_whi_j <- 
      ggplot(effect_var_b_df, aes(x = estimate_c_perc_white, y = fit)) +
      geom_line() +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      theme_bw() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") +
      geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))+
      coord_cartesian(ylim = c(0,10))
    
  fig_si2 <- cowplot::plot_grid(fig_si2_pov_a, fig_si2_whi_b,
                                fig_si2_pov_c, fig_si2_whi_d,
                                fig_si2_pov_e, fig_si2_whi_f,
                                fig_si2_pov_g, fig_si2_whi_h,
                                fig_si2_pov_i, fig_si2_whi_j,
                                ncol = 2, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"))

  ggsave(fig_si2, filename = "fig_SI_2_260108.jpeg", dpi = 300, width = 6, height = 9, units = "in")

  
  
  
  
  
##############################################

  
  m3 <- lmer(plot_perc_damaged   ~ estimate_c_perc_poverty +  estimate_c_perc_white +  
               plot_planted_perc  + plot_street_tree_perc + (1|city), data = ufia_acs)
  m3_summary <- summary(m3)
  print(m3_summary)
  check_collinearity(m3)
  
  #exploring a new figure including both poverty and whiteness
  pov_df <- data.frame(estimate_c_perc_poverty = 0:100)
  whi_df <- data.frame(estimate_c_perc_white = 0:100)
  
  pred_df2 <- expand_grid(pov_df, whi_df) %>% 
    mutate(city = "Baltimore, MD")
  
  pred_df3 <-  pred_df2 %>% 
    mutate(pred = predict(m3, newdata = pred_df2) + 0.1073) #account for the effect of Baltimore 
  
  ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
    geom_tile() + scale_fill_viridis_c(name = "LAI", option = "inferno") + xlab("poverty (%)") + ylab("white (%)") +
    ggthemes::theme_few()
  
  
  
  
  
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




# Calculate the effects for var2
effect_var2 <- effects::effect("estimate_c_perc_white", m1)
plot(effect_var2)

# Calculate the interaction effect
# Replace "var1*var2" with your actual interaction term
interaction_effect <- effect("estimate_c_perc_poverty*estimate_c_perc_white", m1, confidencelevel = 0.95)


# Plot the interaction
plot(interaction_effect, multiline = TRUE, confint = TRUE, ci.style = "bars",
     main = "Effect of Var1 and Var2 Interaction",
     xlab = "Variable 1", ylab = "Response Variable")


# Convert the effect object to a data frame
effect_df <- as.data.frame(interaction_effect)

# Use ggplot2 to create a custom plot
ggplot(effect_df, aes(x = var1, y = fit, color = var2, group = var2)) +
  geom_point() +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = var2), alpha = 0.2) +
  theme_classic() +
  labs(title = "Custom Interaction Plot",
       x = "Variable 1",
       y = "Predicted Response")





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