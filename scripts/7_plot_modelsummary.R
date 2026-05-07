# ---
# title: Plot visual summary of final models
# author: Bryan Alpecho
# date: 2025-
# output: Figure 3, 4, and 5
# ---

#setup
library(tidyverse)
library(patchwork)
library(ggpointdensity)
library(marginaleffects)
library(glmmTMB)

date <- "13042026"

plot_model_summary_omnivores <- function(date){
  
  #Supplementary Figure 1
  #setup: read in df and model
  df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))
  
  print("Visual summary of glm for omnivorous zooplankton in preparation")
  #01 to plot predictions by (A) Chl-a 
  #Omnivores proportion vs. Chl-a
  pop_preds_omni_chla <- predictions(Omni_mdl_zib,
                                     newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.001)),
                                     re.form = NA) # Zeros out random effects
  
  omni_plot_chla <- ggplot() + pub_theme +
    geom_point(data = df_noNA,
               aes(x = chla_sqrt, y = ROC_SVT_zib), alpha = 0.1) +
    geom_ribbon(data = pop_preds_omni_chla, aes(x = chla_sqrt, y = estimate,
                                                ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) +
    geom_line(data = pop_preds_omni_chla, aes(x = chla_sqrt, y = estimate),colour="blue") +
    labs(y = "Proportion of omnivores", x = expression(bold(sqrt("Chl-a"))))
  
  #02 plot intercepts and slope of (A) Longhurst Provinces and (B) Tow No. and Days
  print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
  REs <- ranef(Omni_mdl_zib, condVar = TRUE) 
  TG <- "Omni"
  
  re_lh <- REs$cond$longhurst
  
  ### Longhurst Province ###
  qq <- attr(re_lh, "condVar")
  
  # Extract intercepts
  rand.interc <- REs$cond$longhurst
  
  # Make a dataframe for plotting
  df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                        sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                        lev.names = factor(rownames(rand.interc))) %>% 
    arrange(Intercepts) %>% 
    within({  # Reorder levels
      lev.names <- factor(as.character(lev.names),
                          as.character(lev.names))
    })
  
  re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
    geom_hline(yintercept=0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                      ymax=Intercepts+sd.interc),
                  width = 0,color="black") +
    geom_point(color = "black", size = 2) +
    guides(size="none",shape="none") + 
    theme(axis.text.x=element_text(size=10), 
          axis.title.x=element_text(size=13),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    coord_flip() + 
    labs(y = "Intercept", x = "Longhurst Provinces")
  
  #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
  
  ### Tow slope and intercept ###
  
  # Create a dataframe for plotting
  re_tow <- REs$cond$`survey:tow_no`
  
  towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
  
  # Specify frequency breaks
  my_breaks <- c(1,5,30,175,1000)
  
  #
  re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
    pub_theme +
    geom_pointdensity() +
    scale_colour_viridis_c() +
    labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
    theme(
      axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
    ) 
  
  #03 Combine plots per trophic group. 
  design <- "
                      AAAAA#
                      AAAAA#
                      BBBCCC
                      BBBCCC
                      BBBCCC
                      BBBCCC
                    "
  
  #Figure 3. Omnivores
  print("Patching plots together")
  (final_patch <- omni_plot_chla + re_lh_plot + re_tow_plot_point_density +
      plot_layout(design = design) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )    
  ggsave(paste("output/plots/Omni_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  print(paste("Saved plot: output/plots/Omni_",date,".png",sep=""))
}

#FIGURE 1
plot_model_summary_carnivores <- function(date){
  # Create a publication-ready theme (Adapted from 2025 UQ MME Lab Winter R Workshop)
  # #01 read in model
  # load("output/previousModels/revision/fMdl_carni.RData") 
  # df <- read_csv(paste0("data_input/global_df_complete_",date,".rds"))
  df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))
  
  print("Visual summary of glm for carnivorous zooplankton in preparation")
  #B. Carnivores
  #01 to generate predictions for Carnivores proportion vs. Chl-a
  pop_preds_carni_chla <- predictions(Carni_mdl_zib, 
                                      newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.001)), 
                                      re.form = NA) # Zeros out random effects
  
  carni_plot_chla <- ggplot(data = pop_preds_carni_chla) + pub_theme + 
    geom_point(data = df_noNA, 
               aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
    geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                    ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype = 2) + 
    geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
    labs(y = "Proportion of carnivores", x = expression(bold(sqrt("Chl-a")))) 
  
  #02 plot intercepts and slope of (A) Longhurst Provinces and (B) Tow No. and Days
  REs <- ranef(Carni_mdl_zib, condVar = TRUE) 
  TG <- "Carni"
  re_lh <- REs$cond$longhurst
  
  ### Longhurst Province ###
  qq <- attr(re_lh, "condVar")
  
  # Extract intercepts
  rand.interc <- REs$cond$longhurst
  
  # Make a dataframe for plotting
  df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                        sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                        lev.names = factor(rownames(rand.interc))) %>% 
    arrange(Intercepts) %>% 
    within({  # Reorder levels
      lev.names <- factor(as.character(lev.names),
                          as.character(lev.names))
    })
  
  re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
    geom_hline(yintercept=0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                      ymax=Intercepts+sd.interc),
                  width = 0,color="black") +
    geom_point(color = "black", size = 2) +
    guides(size="none",shape="none") + 
    theme(axis.text.x=element_text(size=10), 
          axis.title.x=element_text(size=13),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    coord_flip() + 
    labs(y = "Intercept", x = "Longhurst Provinces")
  
  #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
  
  ### Tow slope and intercept ###
  print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
  # Create a dataframe for plotting
  re_tow <- REs$cond$`survey:tow_no`
  
  towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
  
  # Specify frequency breaks
  my_breaks <- c(1,5,30,175,1000)
  
  #
  re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
    pub_theme +
    geom_pointdensity() +
    scale_colour_viridis_c() +
    labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
    theme(
      axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
    ) 
  
  #03 Combine plots per trophic group. 
  design <- "
                      AAAAA#
                      AAAAA#
                      BBBCCC
                      BBBCCC
                      BBBCCC
                      BBBCCC
                    "
  
  #Figure 4. Carnivores
  print("Patching plots together")
  (final_patch <- carni_plot_chla + re_lh_plot + re_tow_plot_point_density +
      plot_layout(design = design) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )    
  
  ggsave(paste("output/plots/Carni_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  print(paste("Plot saved: output/plots/Carni_",date,".png",sep=""))
}


plot_model_summary_filterfeeders <- function(date){
  # #01 read in model
  # load("output/previousModels/revision/fMdl_filterfeeder.RData") 
  # df_filter <- read_csv(paste0("data_input/global_df_complete_Filter_",date,".rds"))
  df_noNA <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))
  
  #FIGURE 2
  print("Visual summary of glm for gelatinous filter-feeders in preparation")
  #C. Filter-feeders
  #01 to generate predictions for Filter-feeders proportion vs. Chl-a
  pop_preds_FF_chla <- predictions(Filter_mdl_zib, 
                                   newdata = datagrid(chla_sqrt = seq(floor(min(df_noNA$chla_sqrt)),ceiling(max(df_noNA$chla_sqrt)),0.001)), 
                                   re.form = NA) # Zeros out random effects
  
  ff_plot_chla <- ggplot(data = pop_preds_FF_chla) + pub_theme + 
    geom_point(data = df_noNA, 
               aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
    geom_ribbon(aes(x = chla_sqrt, y = estimate, 
                    ymin = conf.low, ymax = conf.high), alpha = 0.3, colour = "blue", linetype=2) + 
    geom_line(aes(x = chla_sqrt, y = estimate),colour="blue") + 
    labs(y = "Proportion of filter-feeders", x = expression(bold(sqrt("Chl-a")))) 
  
  #02 plot intercepts and slope of (A) Longhurst Provinces and (B) Tow No. and Days
  print("Plotting intercepts and slope of Longhurst Provinces and Tow within Survey")
  ###Longhurst Province
  # Extract random effects
  REs <- ranef(Filter_mdl_zib, condVar = TRUE) #Update the model input
  TG <- "Filter"
  
  
  
  re_lh <- REs$cond$longhurst
  
  ### Longhurst Province ###
  qq <- attr(re_lh, "condVar")
  
  # Extract intercepts
  rand.interc <- REs$cond$longhurst
  
  # Make a dataframe for plotting
  df_plot <- data.frame(Intercepts = REs$cond$longhurst[,1],
                        sd.interc = 2*sqrt(qq[,,1:length(qq)]),
                        lev.names = factor(rownames(rand.interc))) %>% 
    arrange(Intercepts) %>% 
    within({  # Reorder levels
      lev.names <- factor(as.character(lev.names),
                          as.character(lev.names))
    })
  
  re_lh_plot <- ggplot(df_plot, aes(lev.names, Intercepts)) + pub_theme + 
    geom_hline(yintercept=0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin=Intercepts-sd.interc, 
                      ymax=Intercepts+sd.interc),
                  width = 0,color="black") +
    geom_point(color = "black", size = 2) +
    guides(size="none",shape="none") + 
    theme(axis.text.x=element_text(size=10), 
          axis.title.x=element_text(size=13),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank()) +
    coord_flip() + 
    labs(y = "Intercept", x = "Longhurst Provinces")
  
  #annotate("text", y = -2.3, x = 28, label = "(D)", size = 5)
  
  ### Tow slope and intercept ###
  
  # Create a dataframe for plotting
  re_tow <- REs$cond$`survey:tow_no`
  
  towDf <- as.data.frame(re_tow) %>% setNames(c("Intercept", "Slope"))
  
  # Specify frequency breaks
  my_breaks <- c(1,5,30,175,1000)
  
  #
  re_tow_plot_point_density <- ggplot(towDf, aes(Intercept, Slope)) +
    pub_theme +
    geom_pointdensity() +
    scale_colour_viridis_c() +
    labs(colour = "Density", x = "Tow within survey intercept", y = "Tow within survey slope") +
    theme(
      axis.text = element_text(size = 14), legend.position.inside = c(0.85, 0.85)    # Adjust the size for axis numbers/text
    ) 
  
  #03 Combine plots per trophic group. 
  design <- "
                      AAAAA#
                      AAAAA#
                      BBBCCC
                      BBBCCC
                      BBBCCC
                      BBBCCC
                    "
  print("Patching plots together")
  
  (final_patch <- ff_plot_chla + re_lh_plot + re_tow_plot_point_density +
      plot_layout(design = design) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.title = element_text(size = 12, face = "bold"))
      )
  )    
  
  ggsave(paste("output/plots/Filter_",date,".png",sep=""), plot = final_patch,
         width = 8, height = 10, dpi = 300)
  print(paste("Plot saved: output/plots/Filter_",date,".png",sep=""))
}  

#to plot the model estimates per survey in response scale 
#Supplementary Figure 2
plot_model_perSurvey_responseScale <- function(date){
  
  #Carnivores proportion per CPR Survey
  df_noNA_carni <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RCO_SVT_zib)) %>% filter(!is.na(tow_days))
  
  pop_preds_carni_surveys <- predictions(Carni_mdl_zib,
                                         newdata = datagrid(survey = unique(df_noNA_carni$survey),
                                                            chla_sqrt = seq(floor(min(df_noNA_carni$chla_sqrt)),ceiling(max(df_noNA_carni$chla_sqrt)),0.01)),
                                         re.form = NA) # Zeros out random effects
  print("Plotting for carnivores")
  
  carni_plot_surveys <- ggplot(data = pop_preds_carni_surveys) + pub_theme +
    geom_point(data = df_noNA_carni,
               aes(x = chla_sqrt, y = RCO_SVT_zib), alpha = 0.1) +
    geom_ribbon(aes(x = chla_sqrt, y = estimate,
                    ymin = conf.low, ymax = conf.high,
                    colour = survey, fill = survey), alpha = 0.3, colour = NA, linetype = 2) +
    geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) +
    labs(y = "Proportion of carnivores", x = expression(bold(sqrt("Chl-a")))) +
    theme(legend.position = "none", axis.title = element_text(size = 12))
  
  #Omnivores proportion per CPR Survey
  df_noNA_omni <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(ROC_SVT_zib)) %>% filter(!is.na(tow_days))
  
  pop_preds_omni_surveys <- predictions(Omni_mdl_zib,
                                        newdata = datagrid(survey = unique(df_noNA_omni$survey),
                                                           chla_sqrt = seq(floor(min(df_noNA_omni$chla_sqrt)),ceiling(max(df_noNA_omni$chla_sqrt)),0.001)),
                                        re.form = NA) # Zeros out random effects
  print("Plotting for omnivores")
  omni_plot_surveys <- ggplot(data = pop_preds_omni_surveys) + pub_theme +
    geom_point(data = df_noNA_omni,
               aes(x = chla_sqrt, y = ROC_SVT_zib), alpha = 0.1) +
    geom_ribbon(aes(x = chla_sqrt, y = estimate,
                    ymin = conf.low, ymax = conf.high,
                    colour = survey, fill = survey), alpha = 0.3, colour = NA) +
    geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) +
    labs(fill = "Survey", y = "Proportion of omnivores", x = expression(bold(sqrt("Chl-a")))) +
    theme(legend.position = 'none', axis.title = element_text(size = 12)) 
  
  #Filter-feeders proportion per CPR Survey
  df_noNA_FF <- df %>% filter(!is.na(chla_sqrt)) %>% filter(!is.na(RFF_SVT_zib)) %>% filter(!is.na(tow_days))
  
  pop_preds_FF_surveys <- predictions(Filter_mdl_zib,
                                      newdata = datagrid(survey = unique(df_noNA_FF$survey),
                                                         chla_sqrt = seq(floor(min(df_noNA_FF$chla_sqrt)),ceiling(max(df_noNA_FF$chla_sqrt)),0.001)),
                                      re.form = NA) # Zeros out random effects
  print("Plotting for filter-feeders")
  ff_plot_surveys <- ggplot(data = pop_preds_FF_surveys) + pub_theme +
    geom_point(data = df_noNA_FF,
               aes(x = chla_sqrt, y = RFF_SVT_zib), alpha = 0.1) +
    geom_ribbon(aes(x = chla_sqrt, y = estimate,
                    ymin = conf.low, ymax = conf.high,
                    colour = survey, fill = survey), alpha = 0.3, colour = NA) +
    geom_line(aes(x = chla_sqrt, y = estimate, colour = survey), show.legend = F) +
    labs(fill = "Survey", y = "Proportion of filter-feeders", x = expression(bold(sqrt("Chl-a")))) +
    theme(legend.position = "right", axis.title = element_text(size = 12)) +
    scale_fill_discrete(labels=c('Australian CPR', 'Atlantic CPR','North Pacific CPR','SCAR Southern Ocean CPR'))
  
  #to combine plots per trophic group. 
  design <- "AAAABBB
               AAAABBB
               #CCCC##
               #CCCC##"
  
  #Figure
  print("Patching plots together")
  
  (final_patch <- carni_plot_surveys + omni_plot_surveys + ff_plot_surveys +
      plot_layout(design = design) +
      plot_annotation(
        tag_levels = "A",
        theme = theme(plot.title = element_text(size = 18, face = "bold"),
                      plot.margin = unit(c(0.2,0.2,0.2,0.2),"cm"))
      )
  )    
  
  ggsave(paste("output/plots/ResponsePerSurvey_",date,".png",sep=""), plot = final_patch,
         width = 10, height = 6, dpi = 300)
  print(paste("Plot saved: output/plots/ResponsePerSurvey_",date,".png",sep=""))
  
}
  
plot_model_summary_omnivores(date)

plot_model_summary_carnivores(date)

plot_model_summary_filterfeeders(date)

plot_model_perSurvey_responseScale(date)

## End