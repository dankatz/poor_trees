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
library(DHARMa)
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

ufia_acs_indiv <- read_csv( file.path(csv_out_path, "ufia_acs_indiv_for_analysis_260720.csv")) %>% 
  filter(!is.na(city)) #remove the plots that did not connect to census data

#length(unique(ufia_acs_indiv$PLT_CN))

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
  
#histogram of number of trees per plot
hist(ufia_acs$plot_n_trees, breaks = 200)
ufia_acs |> group_by(plot_n_trees) |> summarize(n = n())
median(ufia_acs$plot_n_trees)

#number of plots
nrow(ufia_acs)
  
#number of "cities"
ufia_acs %>% 
  select(city) %>% 
  distinct() %>% 
  summarize(n = n())

#proportion of plots that are in forested areas
ufia_acs |> 
  group_by(PLOT_STATUS_CD_LAB) |> 
  summarize(n = n())
  #589/(589 + 2574)

#proportion of trees that are street trees
ufia_acs_indiv |> 
  group_by(tree_street_tree ) |> 
  summarize(n = n())
  1387 /(1387 + 20925)

#proportion of trees that are planted
ufia_acs_indiv |> 
  group_by(tree_planted ) |> 
  summarize(n = n())
  4896/(4896 + 17422)

  #proportion of trees that are damaged
  ufia_acs_indiv |> 
    group_by(damaged   ) |> 
    summarize(n = n())
  5131/(5131 + 18528)
  
#proportion of trees that are dead
  ufia_acs_indiv |> 
    group_by(tree_alive  ) |> 
    summarize(n = n())
  2987/(2987 + 20672)
  

#mean LAI
  mean(ufia_acs$plot_mean_LAI, na.rm = T )
  sd(ufia_acs$plot_mean_LAI, na.rm = T )
  
## data visualization and analysis #################################################################
### coarse map of study sites #################################
ufia_acs %>% 
  mutate(lat_1 = round(LAT, 1),
         lon_1 = round(LON, 1)) %>% 
  ggplot(aes(x = lon_1, y = lat_1)) + geom_point()




### Fig 1: comparison of tree damage and LAI with poverty and whiteness ############################
# formerly SI 2
  
  #helper function to extract p-values from models
    get_term_stats <- function(model, term_name) {
      tidy(model, effects = "fixed", conf.int = TRUE) %>%
        filter(term == term_name) %>%
        mutate(p_label = paste0("p = ", format.pval(p.value, digits = 2, eps = 0.001))) |> 
        mutate(p_label = case_when(p_label == "p = <0.001" ~ "p ≤ 0.001",
                                   .default = p_label))
    }
     
study_wide_max_poverty <- max(ufia_acs$estimate_c_perc_poverty)

### panel A, B: city level relationships between damage and poverty and whiteness
    ufia_acs_indiv_focal <- ufia_acs_indiv %>% 
      mutate(focal_var = damaged )
    focal_y_lab <- expression("damage (%)")
    m1_a <- glmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city/PLT_CN), 
                  data = ufia_acs_indiv_focal, family = binomial(link = "logit"))
    m1_a_summary <- summary(m1_a)
    
    print(m1_a_summary)
    check_collinearity(m1_a)
    VarCorr(m1_a)
    r2(m1_a)
    
    #check DHARMa residuals
    sim_res <- simulateResiduals(m1_a, n = 1000)
    plot(sim_res)
    
        # Test whether residual pattern differs systematically by city — 
        # directly relevant to your random-slope question
        testCategorical(sim_res, catPred = ufia_acs_indiv_focal$city)
    
        # Formal test for non-linearity/misspecification
        testResiduals(sim_res)
        testOutliers(sim_res)
        plot(sim_res)  # outliers are flagged in the residual-vs-predicted panel
        outliers(sim_res) 
        ufia_acs_indiv_focal[outliers(sim_res), ]
        ranef(m1_a)
        performance::r2_tjur(m1_a)
        
    poverty_p <-  get_term_stats(m1_a, "estimate_c_perc_poverty") |> pull(p_label)
    white_p <- get_term_stats(m1_a, "estimate_c_perc_white") |> pull(p_label)
    
    ### poverty 
    ## individual city level
    pred_grid_city_pov <- 
      expand_grid(
        city = unique(ufia_acs$city),
        estimate_c_perc_poverty = 0:study_wide_max_poverty
      ) |> 
      mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
    
    
    pred_focal_pov <- pred_grid_city_pov %>% 
      mutate(fit = predict(m1_a, pred_grid_city_pov, re.form = ~(1|city), type = "response"))
    
    #overall effect  
    effect_var_a <- effects::effect("estimate_c_perc_poverty", m1_a, xlevels = 100)       # Calculate the effects 
    effect_var_a_df <- as.data.frame(effect_var_a)
    
    #create plot    
    fig_1_pov_a <- 
      ggplot() +
      geom_line(data = pred_focal_pov, aes(x = estimate_c_perc_poverty, y = fit * 100, group = city), alpha = 0.2) +
      
      geom_line(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit* 100), color = "blue") +
      geom_ribbon(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit * 100, 
                                              ymin = lower * 100, ymax = upper * 100), fill = "blue", alpha = 0.2) +
      annotate("text", label = poverty_p, x = 15, y = 80)+
      ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") 
    #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_poverty, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))


    ### whiteness
    ## individual city level
    pred_grid_city_white <- 
      expand_grid(
        city = unique(ufia_acs$city),
        estimate_c_perc_white = 0:100
      ) |> 
      mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
    
    
    pred_focal_white <- pred_grid_city_white %>% 
      mutate(fit = predict(m1_a, pred_grid_city_white, re.form = ~(1|city), type = "response"))
    
    #overall effect  
    effect_var_b <- effects::effect("estimate_c_perc_white", m1_a, xlevels = 100)       # Calculate the effects 
    effect_var_b_df <- as.data.frame(effect_var_b)
    
    #create plot    
    fig_1_white_b <- 
      ggplot() +
      geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit * 100, group = city), alpha = 0.2) +
      
      geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit* 100), color = "blue") +
      geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit* 100, 
                                              ymin = lower* 100, ymax = upper* 100), fill = "blue", alpha = 0.2) +
      annotate("text", label = white_p, x = 85, y = 63)+
      
      ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
    #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))

    
    

  ### panel C, D: city level relationships between trees alive and poverty and whiteness
    ufia_acs_indiv_focal <- ufia_acs_indiv_focal %>% 
      mutate(focal_var = tree_alive ,
             obs_id = 1:nrow(.))
    focal_y_lab <- expression("alive (%)")
    m1_c_city <- glmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city/PLT_CN), data = ufia_acs_indiv_focal,
                  family = binomial(link = "logit"))
    isSingular(m1_c_city)
    ranef(m1_c_city)$city
    
    m1_c <- glmer(focal_var ~ estimate_c_perc_poverty + estimate_c_perc_white + (1 | PLT_CN),
                        data = ufia_acs_indiv_focal, family = binomial(link = "logit"))

    #should a city-level effect be included?
    anova(m1_c_city, m1_c) #no, do not include a city-level random effect for this model
    m1_c_summary <- summary(m1_c)
    check_collinearity(m1_c)
    print(m1_c_summary)
    r2(m1_c)
    
    #check DHARMa residuals
        sim_res <- simulateResiduals(m1_c, n = 1000)
        plot(sim_res)
        
        # Test whether residual pattern differs systematically by city — 
        # directly relevant to your random-slope question
        testCategorical(sim_res, catPred = ufia_acs_indiv_focal$city)
        
        # Formal test for non-linearity/misspecification
        testResiduals(sim_res)
        testOutliers(sim_res)
        testDispersion(sim_res) #the result is significant statistically, but the difference is very small (1.00 vs 1.06)
        plot(sim_res)  # outliers are flagged in the residual-vs-predicted panel
        outliers(sim_res) 
        ufia_acs_indiv_focal[outliers(sim_res), ]
        ranef(m1_c)
    
    poverty_p <-  get_term_stats(m1_c, "estimate_c_perc_poverty") |> pull(p_label)
    white_p <- get_term_stats(m1_c, "estimate_c_perc_white") |> pull(p_label)
    
    
    ### poverty 
    ## individual plot level
    # pred_grid_city_pov <- 
    #   expand_grid(
    #     PLT_CN = unique(ufia_acs$PLT_CN),
    #     estimate_c_perc_poverty = 0:study_wide_max_poverty
    #   ) |> 
    #   mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
    # 
    # pred_focal_pov <- pred_grid_city_pov %>% 
    #   mutate(fit = predict(m1_c, pred_grid_city_pov, re.form = ~(1|PLT_CN), type = "response"))
    # 
    
    #overall effect  
    effect_var_a <- effects::effect("estimate_c_perc_poverty", m1_c, xlevels = 100)       # Calculate the effects 
    effect_var_a_df <- as.data.frame(effect_var_a)
    
    #create plot    
    fig_1_pov_e <- 
      ggplot() +
      #geom_line(data = pred_focal_pov, aes(x = estimate_c_perc_poverty, y = fit * 100, group = PLT_CN), alpha = 0.2) +
      
      geom_line(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit* 100), color = "blue") +
      geom_ribbon(data = effect_var_a_df, aes(x = estimate_c_perc_poverty, y = fit* 100, 
                                              ymin = lower * 100, ymax = upper* 100), fill = "blue", alpha = 0.2) +
      annotate("text", label = poverty_p, x = 15, y = 98) +
      ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood poverty (%)") 

    
    ### whiteness
    # ## individual city level
    # pred_grid_city_white <- 
    #   expand_grid(
    #     city = unique(ufia_acs$city),
    #     estimate_c_perc_white = 0:100
    #   ) |> 
    #   mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
    # 
    # pred_focal_white <- pred_grid_city_white %>% 
    #   mutate(fit = predict(m1_c, pred_grid_city_white, re.form = ~(1|city), type = "response"))
    
    #overall effect  
    effect_var_b <- effects::effect("estimate_c_perc_white", m1_c, xlevels = 100)       # Calculate the effects 
    effect_var_b_df <- as.data.frame(effect_var_b)
    
    #create plot    
    fig_1_white_f <- 
      ggplot() +
      # geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit * 100, group = city), alpha = 0.2) +
      
      geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit* 100), color = "blue") +
      geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit* 100, 
                                              ymin = lower * 100, ymax = upper * 100), fill = "blue", alpha = 0.2) +
      annotate("text", label = white_p, x = 15, y = 97)+
      ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
    #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
    
    
  ### panel C, D: city level relationships between LAI and poverty and whiteness
    ufia_acs_focal <- ufia_acs %>% 
      mutate(focal_var = plot_mean_LAI )
    focal_y_lab <- expression("mean LAI (m"^2*"/m"^2*")")
    m1_b <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs_focal)
    m1_b_summary <- summary(m1_b)
    check_collinearity(m1_b)
    print(m1_b_summary)
    poverty_p <- paste("p =", round( m1_b_summary$coefficients[2,5], 3))
    white_p <- paste("p =", round( m1_b_summary$coefficients[3,5], 3))
    r2(m1_b)
    performance::icc(m1_b, by_group = TRUE)
    performance::check_model(m1_b)
    
    
    ### poverty 
    ## individual city level
    pred_grid_city_pov <- 
      expand_grid(
        city = unique(ufia_acs$city),
        estimate_c_perc_poverty = 0:study_wide_max_poverty
      ) |> 
      mutate(estimate_c_perc_white = mean(ufia_acs$estimate_c_perc_white))
    
    pred_focal_pov <- pred_grid_city_pov %>% 
      mutate(fit = predict(m1_b, pred_grid_city_pov))
    
    #overall effect  
    effect_var_a <- effects::effect("estimate_c_perc_poverty", m1_b, xlevels = 100)       # Calculate the effects 
    effect_var_a_df <- as.data.frame(effect_var_a)
    
    #create plot    
    fig_1_pov_c <- 
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
    pred_grid_city_white <- 
      expand_grid(
        city = unique(ufia_acs$city),
        estimate_c_perc_white = 0:100
      ) |> 
      mutate(estimate_c_perc_poverty = mean(ufia_acs$estimate_c_perc_poverty))
    
    pred_focal_white <- pred_grid_city_white %>% 
      mutate(fit = predict(m1_b, pred_grid_city_white))
    
    #overall effect  
    effect_var_b <- effects::effect("estimate_c_perc_white", m1_b, xlevels = 100)       # Calculate the effects 
    effect_var_b_df <- as.data.frame(effect_var_b)
    
    #create plot    
    fig_1_white_d <- 
      ggplot() +
      geom_line(data = pred_focal_white, aes(x = estimate_c_perc_white, y = fit, group = city), alpha = 0.2) +
      
      geom_line(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit), color = "blue") +
      geom_ribbon(data = effect_var_b_df, aes(x = estimate_c_perc_white, y = fit, 
                                              ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
      annotate("text", label = white_p, x = 15, y = 7.4)+
      ggthemes::theme_few() + ylab(focal_y_lab) + xlab("neighborhood whiteness (%)") 
    #geom_rug(data = ufia_acs_focal, aes(x = estimate_c_perc_white, y = focal_var), alpha = 0.05, position = position_jitter(height = 2))
    
    
    
    
    ### save fig 1
    fig_1 <- cowplot::plot_grid(fig_1_pov_a, fig_1_white_b,
                                fig_1_pov_e, fig_1_white_f,
                                fig_1_pov_c, fig_1_white_d,
                                
                                ncol = 2, labels = c("A", "B", "C", "D", "E", "F"),
                                rel_heights = c(1, 1, 1.3))
    fig_1
    
    ggsave(fig_1, filename = "fig_1_260810.jpeg",  width = 7, height = 10, units = "in", dpi = 300) #may need some resizing
    
    
    # 
    # ufia_acs_focal <- ufia_acs %>% 
    #   mutate(focal_var = plot_perc_damaged )
    # focal_y_lab <- expression("damage (%)")
    # m_si2_focal <- lmer(focal_var  ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city/PLT_CN), data = ufia_acs_indiv_focal)
    # m_si2_focal_summary <- summary(m_si2_focal)
    # check_collinearity(m_si2_focal)
    # print(m_si2_focal_summary)
    # VarCorr(m_si2_focal)
    # 

    
### Fig 2: predicted combined effect of both poverty and whiteness on tree damage, % alive, and LAI ########################
   #need to run the models from Fig. 1 before running this
    
    ## Panel A: the combined effect of both poverty and whiteness on tree damage
     m2_a <- m1_a#lmer(damaged ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city/PLT_CN), data = ufia_acs_indiv_focal)
    m2_a_summary <- summary(m2_a)
    print(m2_a_summary)
    
    pov_df <- data.frame(estimate_c_perc_poverty = 0:round(max(ufia_acs$estimate_c_perc_poverty), 0))
    whi_df <- data.frame(estimate_c_perc_white = 0:100)
    
    pred_df2 <- expand_grid(pov_df, whi_df)
    
    pred_df3 <-  pred_df2 %>% 
      mutate(pred = predict(m2_a, newdata = pred_df2, re.form = NA, type = "response") ) # no random effects — population-level
    
    fig2_a <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred * 100)) +
          geom_tile() + scale_fill_viridis_c(name = "damage (%)", option = "plasma") + xlab("poverty (%)") + ylab("White (%)") +
            ggthemes::theme_few() + theme(legend.position = "bottom") +
          geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                            inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")

    ## Panel B: the combined effect of both poverty and whiteness on tree death
    m2_c <- m1_c #lmer(tree_alive ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city/PLT_CN), data = ufia_acs_indiv_focal)
    
    pov_df <- data.frame(estimate_c_perc_poverty = 0:round(max(ufia_acs$estimate_c_perc_poverty), 0))
    whi_df <- data.frame(estimate_c_perc_white = 0:100)
    
    pred_df2 <- expand_grid(pov_df, whi_df)
    
    pred_df3 <-  pred_df2 %>% 
      mutate(pred = predict(m2_c, newdata = pred_df2, re.form = NA, type = "response") ) # no random effects — population-level
    
    fig2_b <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred * 100)) +
      geom_tile() + scale_fill_viridis_c(name = "alive (%)", option = "rocket", direction = -1) + xlab("poverty (%)") + ylab("White (%)") +
      ggthemes::theme_few() + theme(legend.position = "bottom") +
      geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                 inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")
    
    
    
    ## Panel C: the combined effect of both poverty and whiteness on LAI
    m2_b <- m1_b#lmer(plot_mean_LAI ~ estimate_c_perc_poverty +  estimate_c_perc_white + (1|city), data = ufia_acs)
    
    pov_df <- data.frame(estimate_c_perc_poverty = 0:round(max(ufia_acs$estimate_c_perc_poverty), 0))
    whi_df <- data.frame(estimate_c_perc_white = 0:100)
    
    pred_df2 <- expand_grid(pov_df, whi_df) 
    
    pred_df3 <-  pred_df2 %>% 
      mutate(pred = predict(m2_b, newdata = pred_df2, re.form = NA) ) #account for the effect of Washington DC 
    
    fig2_c <- ggplot(pred_df3, aes(x= estimate_c_perc_poverty, y = estimate_c_perc_white, fill = pred)) +
      geom_tile() + scale_fill_viridis_c(name = expression("mean LAI (m"^2*"/m"^2*")"), 
                                         option = "viridis", direction = -1) + xlab("poverty (%)") + ylab("White (%)") +
      ggthemes::theme_few() + theme(legend.position = "bottom") +
      geom_point(data = ufia_acs, aes(x = estimate_c_perc_poverty, y = estimate_c_perc_white), 
                 inherit.aes = FALSE, alpha = 0.4, size = 0.2, color = "gray")
    
    
    ### save fig 2 
    fig_2 <- cowplot::plot_grid( fig2_a, fig2_b, fig2_c,
                                 ncol = 3, labels = c("A", "B", "C"))
    fig_2
    ggsave(fig_2, filename = "fig_2_260723.jpeg",  width = 10, height = 4, units = "in", dpi = 300) #may need some resizing
    
    
    
    
### Table 1: model results for each variable #######################################################
    options(pillar.sigfig = 4)
    
    #damage #note: converted to odds ratios to percent to match what's in the table and figs
    tidy(m1_a, effects = "fixed", conf.int = TRUE, conf.method = "Wald", exponentiate = TRUE) %>%
      filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
      select(term, estimate, conf.low, conf.high)
    
        #showing for a 10% increase
        tidy(m1_a, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
          filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
          mutate(
            or_10unit = round(exp(estimate * 10), 3),
            ci_low_10unit = round(exp(conf.low * 10), 3),
            ci_high_10unit = round(exp(conf.high * 10), 3)
          ) %>%
          select(term, or_10unit, ci_low_10unit, ci_high_10unit)
            
    
    #alive #note: convert to percent to match what's in the table and figs
        test <- tidy(m1_c, effects = "fixed", conf.int = TRUE, conf.method = "Wald", exponentiate = TRUE) %>%
          filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
          select(term, estimate, conf.low, conf.high)
        
        #showing for a 10% increase
        tidy(m1_c, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
          filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
          mutate(
            or_10unit = exp(estimate * 10),
            ci_low_10unit = exp(conf.low * 10),
            ci_high_10unit = exp(conf.high * 10)
          ) %>%
          select(term, or_10unit, ci_low_10unit, ci_high_10unit)
        
    
    #lai
    tidy(m1_b, effects = "fixed", conf.int = TRUE, conf.method = "profile") %>%
      filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
      select(term, estimate, conf.low, conf.high)
    
    tidy(m1_b, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
      filter(term %in% c("estimate_c_perc_poverty", "estimate_c_perc_white")) %>%
      mutate(
        model = "LAI",
        metric = "β (raw coefficient)",
        estimate_final = round(estimate * 10, 3),
        ci_low_final = round(conf.low * 10, 3),
        ci_high_final = round(conf.high * 10, 3)
      )
    
    
### some other stats for the paper ##################################################    
    #extract some statistics about what differences this would result in at different percentiles of race and poverty
      study_wide_median_poverty <- median(ufia_acs$estimate_c_perc_poverty)
      study_wide_median_white <- median(ufia_acs$estimate_c_perc_white)
      study_wide_poverty1 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.01)
      study_wide_poverty99 <- quantile(ufia_acs$estimate_c_perc_poverty, 0.99)
      study_wide_white1 <- quantile(ufia_acs$estimate_c_perc_white, 0.01)
      study_wide_white99 <- quantile(ufia_acs$estimate_c_perc_white, 0.99)
      
      min(ufia_acs$estimate_c_perc_white)
      min(ufia_acs$estimate_c_perc_poverty)
      max(ufia_acs$estimate_c_perc_poverty)
      max(ufia_acs$estimate_c_perc_white)
      
      ### some predictions for the 1st and 99th percentile
      pred_df <- data.frame(estimate_c_perc_poverty = c(study_wide_poverty1, study_wide_poverty99,
                                                        study_wide_median_poverty, study_wide_median_poverty,
                                                        study_wide_poverty99, study_wide_poverty1),
                            estimate_c_perc_white = c(study_wide_median_white, study_wide_median_white,
                                                      study_wide_white1, study_wide_white99,
                                                      study_wide_white1, study_wide_white99)) 
     
      ## predictions for damage  
      preds <- predict(m1_a, newdata = pred_df, re.form = NA, type = "response") # no random effects — population-level
      pred_df %>% mutate(preds = preds)
      
      
      #predictions for trees alive
      preds <- predict(m1_c, newdata = pred_df, re.form = NA, type = "response") # no random effects — population-level
      pred_df %>% mutate(preds = preds)
            
      #predictions for LAI
      preds <- predict(m1_b, newdata = pred_df, re.form = NA) # no random effects — population-level
      pred_df %>% mutate(preds = preds)
      
      
  
  

 
  
### SI 1: does this pattern hold true only within individual cities or in aggregate? ##############################
  tree_level_data <- ufia_acs_indiv_focal %>%
    group_by(city) %>%
    mutate(poverty_city_mean = mean(estimate_c_perc_poverty),
           poverty_within = estimate_c_perc_poverty - poverty_city_mean,
           white_city_mean = mean(estimate_c_perc_white),
           white_within = estimate_c_perc_white - white_city_mean
           ) %>%
    ungroup()
  
  ### damage  #note: convert to percent to match what's in the table and figs
      model_decomp_damage <- glmer(damaged  ~ 
                                     poverty_within + poverty_city_mean + 
                                    white_within + white_city_mean + 
                                    (1|city/PLT_CN), 
                    data = tree_level_data, family = binomial(link = "logit"))
      model_decomp_damage
      summary(model_decomp_damage)
      
      #export coefficients 
      results_to_clipboard <- 
        tidy(model_decomp_damage, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
        filter(term %in% c("poverty_within", "poverty_city_mean", "white_city_mean", "white_within")) %>%
        mutate(
          estimate = round(estimate, 3),
          conf.low = round(conf.low, 3),
          conf.high = round(conf.high, 3)
        ) %>%
          select(term, estimate, conf.low, conf.high)%>%
        mutate(estimate_ci = paste0(estimate, " [", conf.low, ", ", conf.high, "]")) %>%
        select(term, estimate_ci)
      write.table(results_to_clipboard, "clipboard", sep = "\t", row.names = FALSE)   # Windows
        
      
  
  ### trees alive  #note: convert to percent to match what's in the table and figs
  model_decomp_alive <- glmer(tree_alive ~ poverty_within + poverty_city_mean + 
                                white_within + white_city_mean +
                                (1 | PLT_CN ),
                              data = tree_level_data, family = binomial(link = "logit"))
  model_decomp_alive
  summary(model_decomp_alive)
  
  #export coefficients 
  results_to_clipboard <- 
    tidy(model_decomp_alive, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    filter(term %in% c("poverty_within", "poverty_city_mean", "white_city_mean", "white_within")) %>%
    mutate(
      estimate = round(estimate, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3)
    ) %>%
    select(term, estimate, conf.low, conf.high)%>%
    mutate(estimate_ci = paste0(estimate, " [", conf.low, ", ", conf.high, "]")) %>%
    select(term, estimate_ci)
  write.table(results_to_clipboard, "clipboard", sep = "\t", row.names = FALSE)   # Windows
  
  
  
  
  ### tree LAI
  model_decomp_lai <- lmer(plot_mean_LAI  ~ poverty_within + poverty_city_mean + 
                               white_within + white_city_mean +
                               (1 | city ),
                             data = tree_level_data)
  model_decomp_lai
  summary(model_decomp_lai)
  
  #export coefficients 
  results_to_clipboard <- 
    tidy(model_decomp_lai, effects = "fixed", conf.int = TRUE, conf.method = "Wald") %>%
    filter(term %in% c("poverty_within", "poverty_city_mean", "white_city_mean", "white_within")) %>%
    mutate(
      estimate = round(estimate, 3),
      conf.low = round(conf.low, 3),
      conf.high = round(conf.high, 3)
    ) %>%
    select(term, estimate, conf.low, conf.high)%>%
    mutate(estimate_ci = paste0(estimate, " [", conf.low, ", ", conf.high, "]")) %>%
    select(term, estimate_ci)
  write.table(results_to_clipboard, "clipboard", sep = "\t", row.names = FALSE)   # Windows
  
  
  
  
### SI 3: do models that fit a random slope per city work better? ###########################################################################
  
### tree damage
  
  # Random slope for poverty across cities
  model_poverty_slope <- lmer(
    damaged ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 + estimate_c_perc_poverty | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  # Random slope for race across cities
  model_white_slope <- lmer(
    damaged ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 + estimate_c_perc_white | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  #base model
  model_base <- lmer(
    damaged ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  anova(model_base, model_poverty_slope)
  anova(model_base, model_white_slope)
  
  
### tree alive
  
  # Random slope for poverty across cities
  model_poverty_slope <- lmer(
    tree_alive ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 + estimate_c_perc_poverty | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  # Random slope for race across cities
  model_white_slope <- lmer(
    tree_alive ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 + estimate_c_perc_white | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  #base model
  model_base <- lmer(
    tree_alive ~ estimate_c_perc_poverty + estimate_c_perc_white + 
      (1 | city) + (1 | city:PLT_CN),
    data = ufia_acs_indiv_focal,
    REML = FALSE
  )
  
  anova(model_base, model_poverty_slope)
  anova(model_base, model_white_slope)
  