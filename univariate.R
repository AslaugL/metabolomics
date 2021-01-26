#Set working directory
#setwd('C:/Users/Aslaug L/Desktop/Master')
setwd('C:/Users/aslau/Documents/Uni/Master/gitHub')

#Load libraries
library(tidyverse)
library(broom)
library(stringi)
library(readxl)
library(janitor)
library(funModeling)
library(rstatix)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(ggfortify)
library(gghalves)
library(ggbeeswarm)
library(ggiraph)
library(ggrepel)
library(DT)
library(knitr)
library(kableExtra)

#ggplot graphics
theme_set(theme_minimal())

#Read data----
data <- readRDS('tidy_data.Rds')
#List for various df's to keep environment clean
various <- list()

#Functions
filterP <- function(df, filterp, num = 0.05){
  #Filter p.adj value
  if(filterp == 'yes'){
    df %>%
      dplyr::filter(p.adj < num)
  }else if(filterp == 'no'){
    df <- df
  }else {
    stop("filterp must be yes or no")
  }
}
filterCompounds <- function(df, compounds){
  #Filter the compound column for any compounds provided
  if(any(df$compound %in% compounds)){
    df %>% filter(compound %in% compounds) #If any compounds are found in the Compound column, filter the column by these
  }else if(!any(df$compound %in% compounds) & length(compounds) > 0){ #If none of the provided compounds are found, give notice
    print("None of the provided compounds were found in the data")
    df <- df
  }else if(is.null(compounds)){ #If no compounds are provided, do nothing
    df <- df
  }
}
extract_stats <- function(df, filterp = 'yes') {
#Get a string of stats for plots, from a rstatix-results df (t-test/wilcox + respective effect sizes)
  
stats <- df %>%
  filterP(., filterp = filterp) %>% #Filter out significantly different features
  mutate_at(., 'p.adj', ~formatC(., format = 'e', digits = 2)) %>% #Format scientific notation
  mutate_at(., c('estimate', 'conf.low', 'conf.high', 'effsize'), ~round(., digits = 2)) %>% #Format digits
    
  #Write the string based on test type, t-test or wilcox 
  dplyr::mutate(., string = case_when(
  method == "Welch's t-test" ~ sprintf('%s, Welch t test, FDR: %s\n Estimate: %s, ci: %s; %s, Effect size: %s', .$compound, .$p.adj, .$estimate, .$conf.low, .$conf.high, .$effsize),
  method == 'Wilcoxon' ~ sprintf('%s, Wilcox, FDR: %s\n Estimate: %s, ci: %s;  %s,  Effect size: %s', .$compound, .$p.adj, .$estimate, .$conf.low, .$conf.high, .$effsize))
  ) %>% select(compound, string)
}
changeGGplotTxtSize <- function(plot, txt_size = 12){
  #Change the text size of the axis title and text of a ggplot
  plot + theme(axis.title = element_text(size = txt_size),
               axis.text = element_text(size = txt_size),
               legend.title = element_text(size = txt_size),
               legend.text = element_text(size = txt_size))
}
plotViolinBox <- function(df) {

  plot <- ggplot(df, aes(x = health_status, y = value, color = health_status)) +
    
    #Add mock linetypes for the legend
    geom_line(aes(linetype = 'dashed')) +
    geom_line(aes(linetype = 'solid')) +
    
    #Build the jitterplot or dotplot, half violin and half boxplot
    geom_half_violin(side = 'l') +
    geom_half_boxplot(side = 'r', outlier.size = 0) +
    geom_quasirandom(width = 0.2, alpha = 0.5) +
    
    #Add the stippled line for the mean
    stat_summary(fun = mean, geom = 'errorbar',
                 aes(ymax = ..y.., ymin = ..y..),
                 linetype = 'dashed',
                 width = 0.38, position = position_nudge(x = 0.185)) + #Add dashed line to indicate mean
    
    #Make it nicer
    #scale_color_manual(values = )
    #scale_shape_manual(values = ) +
    scale_linetype_manual(name = 'Boxplot summary',
                          values = c('solid', 'dashed'),
                          labels = c('Median', 'Mean')) +
    
    #Set legend to the right
    theme(legend.position = 'right')
}
createViolinBoxPlot <- function(df, stats, filterp = 'yes', compounds = NULL, print = NULL) {
#Create violinboxplots with stats, df = tidy dataframe with , stats = rstatix results list or df
  
  #metabolite data
  data <- df %>% filterCompounds(., compounds = compounds)
  
  #Extract string of stat info, and join with metabolite value data
  plot_df <- stats %>%
    extract_stats(., filterp = filterp) %>%
    inner_join(., data)
  
  #plot
  temp_plot <- plotViolinBox(plot_df)
  
  #Facet_wrap paginate or facet_wrap based on print argument
  if(is.null(print)){
    plot <- temp_plot + facet_wrap(~string, scales = 'free',
                                   ncol = 2)
  } else if (print == 'yes') {
    plot <- temp_plot + facet_wrap_paginate(~string, 
                                            ncol = 1,
                                            nrow = 1,
                                            scales = 'free',
                                            page = 1)
  }else {
    stop("Print must be null or 'yes'")
  }
  
}
createDTable <- function(datatable){
  #Create an interactive table with ability to filter columns and download results
  table <- datatable(datatable, 
                     rownames = FALSE, 
                     filter = 'top',
                     extensions = c('Buttons'),
                     options =
                       list(
                         pageLength = 50, #Set number of features shown at default
                         dom = 'Bfrtip', #Where to put the different table elements
                         buttons = c('copy', 'csv', 'excel', 'print') #Add buttons for saving
                       )
  )
}

### Univariate tests----
#Shap wilk
various$shap_wilk <- data %>%
  group_by(health_status, compound) %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>%
  
  #Does result imply t or wilcox test?
  group_by(compound) %>%
  mutate(test_to_run = case_when(
    any(p.value <= 0.05) ~ 'Wilcox',
    TRUE ~ 't-test'
  )) %>%
  ungroup()

#Run univariate tests using rstatix, t-test or wilcox depending on shapiro wilk result
temp <- list(
  't_test' = various$shap_wilk %>%
    filter(test_to_run == 't-test') %>%
    select(compound) %>%
    inner_join(data) %>% #Get the metabolite values
    group_by(compound) %>%
    
    #run t-test
    t_test(value ~ health_status, detailed = TRUE) %>%
    adjust_pvalue(method = "BH") %>%
    mutate(method = "Welch's t-test"),
  
  #Effect size
  't_test_effsize' = various$shap_wilk %>%
    filter(test_to_run == 't-test') %>%
    select(compound) %>%
    inner_join(data) %>%
    group_by(compound) %>%
    
    #Get cohen's d
    cohens_d(value ~ health_status),
  
  'wilcox' = various$shap_wilk %>%
    filter(test_to_run == 'Wilcox') %>%
    select(compound) %>%
    inner_join(data) %>%
    group_by(compound) %>%
    
    #Run wilcox
    rstatix::wilcox_test(value ~ health_status, detailed = TRUE, paired = FALSE) %>%
    adjust_pvalue(method = "BH") %>%
    ungroup(),
  
  'wilcox_effsize' = various$shap_wilk %>%
    filter(test_to_run == 'Wilcox') %>%
    select(compound) %>%
    inner_join(data) %>%
    group_by(compound) %>%
    
    #Get r
    rstatix::wilcox_effsize(value ~ health_status, ci = FALSE) #ci takes forever
)
univariate_tests <- bind_rows(
  full_join(temp$t_test, temp$t_test_effsize %>% select(-magnitude)),
  full_join(temp$wilcox, temp$wilcox_effsize %>% select(-magnitude))
)

#Half violin/half boxplot
plot <- createViolinBoxPlot(data, univariate_tests, filterp = 'no', compounds = 'SM(18:1)')
plot +
  
  #Fix labs
  labs(
    x = 'Health status',
    y = 'Value'
  ) +
  #Fix tickmarks
  scale_x_discrete(breaks = c('case', 'control'), labels = c('ME/CFS', 'Healthy')) +
  scale_color_manual(values = c('#450A5C', '#7CD04F'), labels = c('ME/CFS', 'Healthy'), name = 'Health status', #Viridis colors
                    guide = guide_legend(override.aes = list(shape = c(21,21), fill = c('#450A5C', '#7CD04F')))) #Filled circular shapes
  
ggsave('example_violinboxplot.png')

#Fold change
fold_change <- tidy_data %>%
  pivot_wider(names_from = 'health_status',
              values_from = 'value') %>%
  mutate(fold_change = case/control)

