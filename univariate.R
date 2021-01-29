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
univariate_tests <- readRDS('univariate_test_results.Rds')
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
FCthreshold <- function(df, p_threshold = 0.05, FC_threshold = 1.5, color_by){
  #Add a 'threshold' column to change the shape of datapoints in the volcano-plot
  #and a 'color' column to color by based on which type of metabolite is over/underexpressed. 
  #df must have a log2(FC) column
  #p.adj = adjusted p.value threshold, FC = fold_change threshold, color_by = column in which information to color by is found
  
  #Get the log2 of the threshold for fold_change to filter by
  fc <- log2(FC_threshold)
  
  #Add threshold column
  added_columns <- df %>%
    mutate(threshold = case_when(
      p.adj < p_threshold & `log2(FC)` > fc ~ 'Overrepresented',
      p.adj < p_threshold & `log2(FC)` < -fc ~ 'Underrepresented',
      TRUE ~ 'Below threshold'
    )) %>%
    
    #Add color_by column
    mutate(color_by = case_when(
      threshold %in% c('Overrepresented', 'Underrepresented') ~ !!as.name(color_by),
      threshold == 'Below threshold' ~ 'Below threshold'
    ))
  
  
}
plotVolcano <- function(df, title, fold_change_threshold = 1.5) {
  #Use ggplot2 to plot a volcano plot
  #Creates vertical lines to separate the different fold change categories, lines are drawn with a default fold change of 2
  
  #Tooltip for ggiraph
  #df$tooltip <- c(paste0("Name: ", df$Item, "\n FDR: ", formatC(df$p.adj, format = 'e', digits = 2), "\n Fold change ", df$fold_change))
  
  #Plot
  plot <- ggplot(df, aes(`log2(FC)`, -log10(p.adj), color = color_by, shape = threshold)) +
    
    #Add the scatters
    geom_point(alpha = 0.5, size = 2, aes(fill = color_by, shape = threshold)) +
    #geom_point_interactive(tooltip=df$tooltip, alpha = 0.5, size = 2) + #Add the  interactive tooltips
    
    #Add lines to show where the significant threshold is
    geom_vline(xintercept = log2(fold_change_threshold), linetype = 'dashed', alpha = 0.5) +
    geom_vline(xintercept = log2(1/fold_change_threshold), linetype = 'dashed', alpha = 0.5) +
    geom_hline(yintercept = 1.3, linetype = 'dashed', alpha = 0.5) +
    
    #Make it nicer, change colors and labels
    #scale_color_manual(values = c('Below threshold' = '#DDDDDD')) + 
    scale_shape_manual(values = color_shape_manuals$shape$Volcanoplot) +
    theme(legend.title = element_blank()) + #Remove legend title
    
    #Change labs
    labs(
      title = title,
      x = "log2 fold change",
      y = "-log10 adjusted p-value")
} #Need to add dynamic color scale
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

temp 

#Plot color/shape manuals----
color_shape_manuals = list(
  'color' = list(
    
    #Color PCA loadings or volcano plots by type of metabolite
    #'volcano' = setNames(c('#DDDDDD', '#88CCEE', '#44AA99'), c('Below threshold', unique(data$SUPER_PATHWAY))),
    
    #Color R2Y and Q2 in permutation plot
    'permutations' = setNames(c('#21A884', '#482374'), c('R2Y', 'Q2Y'))
  ),
  'shape' = list(
    'Volcanoplot' = setNames(c('\u25B2', '\u25CF', '\u25BC'), c('Overrepresented', 'Below threshold', 'Underrepresented')) #Upward/downward triangle for over/underrepresented features, circle for 'below threshold'
  )
)

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
#save univariate tests
#saveRDS(univariate_tests, 'univariate_test_restuls.Rds')

#Fold change
various$fold_change <- data %>%
  
  #Get the median for each compound for both case and control
  group_by(compound, health_status) %>%
  summarise(median = median(value)) %>%
  ungroup() %>%
  
  #Get the super pathway for filtering later
  inner_join(., data %>% select(compound, SUPER_PATHWAY, SUB_PATHWAY)) %>% unique() %>%
  
  #Find fold change
  pivot_wider(names_from = 'health_status',
              values_from = 'median') %>%
  
  mutate(fold_change = case/control) %>% #1 Inf (Cotinine), some NaN (0 values in both case and control groups)
  mutate(`log2(FC)` = log2(fold_change))

#Add fold change to tidy_data
#data <- inner_join(data, various$fold_change %>%
#                     pivot_longer(cols = c(case, control),
#                                  names_to = 'health_status',
#                                  values_to = 'median'))
  
#saveRDS(data, 'tidy_data.Rds')

### Plots----

###
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

#ggsave('example_violinboxplot.png')

###
# Volcano plots


temp <- inner_join(data, univariate_tests) %>%
  FCthreshold(., color_by = 'SUPER_PATHWAY') %>%
  #Add 'Below threshold' at the end of the factor for prettier plots
  mutate_at('color_by', ~as.factor(.)) %>%
  mutate_at(., 'color_by', ~fct_relevel(., 'Below threshold', after = Inf)) %>%
  #Only columns necessary or plot
  select(compound, `log2(FC)`, p.adj, color_by, threshold) %>% unique()

plot <- plotVolcano(temp, 'ME/CFS vs Control')
plot +
  #Change color scale
  scale_color_manual(values = c('Below threshold' = '#DDDDDD', 'Amino Acid' = '#88CCEE', 'Peptide' = '#44AA99')) +
  
  #Add labels when so few mets are different
  geom_label_repel(data = temp %>% filter(threshold != 'Below threshold'), aes(label = compound), show.legend = FALSE)

