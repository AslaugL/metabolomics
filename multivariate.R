#Set working directory
#setwd('C:/Users/Aslaug L/Desktop/Master')
setwd('C:/Users/aslau/Documents/Uni/Master/gitHub/MECFS')

#Load libraries
library(tidyverse)
library(readxl)
library(stringi)
library(janitor)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(ggfortify)
library(gghalves)
library(ggbeeswarm)
library(ggiraph)
library(ggrepel)
library(psych)
library(ropls)
library(DT)
library(knitr)
library(reactable)
library(kableExtra)
library(cowplot)

#ggplot graphics
theme_set(theme_minimal() + theme(legend.position = 'bottom'))

#Read data
raw_data <- read_xlsx('germain2018.xlsx')

tidy_data <- raw_data %>%
  
  #Remove unnecessary columns
  select(-c(`PATHWAY SORTORDER`, `COMP ID`, `PLATFORM`, `CHEMICAL ID`, RI, MASS, CAS, PUBCHEM, KEGG, HMDB)) %>%
  
  #Tidy
  pivot_longer(.,
               cols = -c(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`),
               names_to = 'sample',
               values_to = 'value') %>%
  #Rename BIOCHEMICAL to feature for plotting functions
  rename(feature = BIOCHEMICAL, feature_anno = `SUPER PATHWAY`) %>%
  #Create a column for control/sample status
  mutate(group = case_when(
    str_detect(sample, 'C') ~ 'Control',
    str_detect(sample, 'P') ~ 'ME/CFS'
  ))


#List for various df's to keep environment clean
various <- list()

#Create metadata for annotations----
meta <- list(
  
  #Sample metadata#group is used to annotate samples in plots by default
  'samples' = tidy_data %>% select(sample, group) %>% unique(),
  
  #Feature metadata, for coloring loading plots etc
  'features' = tidy_data %>% select(feature, feature_anno) %>% unique())

#Color manuals
color_manuals <- list(
  
  #color for annotations
  'sample_group' = setNames(c('#44AA99', '#AA4499'),
                    unique(meta$samples$group)),
  'feature_anno' = setNames(c('#DDDDDD', '#88CCEE', '#44AA99', '#117733', '#DDCC77', '#999933', '#882255', '#CC6677', '#AA4499'),
                            c('Below threshold', unique(meta$features$feature_anno)))
  
)


### Multivariate----
changeGGplotTxtSize <- function(plot, txt_size = 10){
  #Change the text size of the axis title and text of a ggplot
  plot + theme(axis.title = element_text(size = txt_size),
               #axis.text = element_text(size = txt_size),
               legend.title = element_text(size = txt_size),
               legend.text = element_text(size = txt_size))
}
addLinesToPlot <- function(plot, yintercept = 0, xintercept = 0){
  #Adds a gray dashed line to the plot at the intercepts
  
  plot <- plot +
    geom_vline(xintercept = xintercept, color = 'gray', linetype = 'dashed') +
    geom_hline(yintercept = yintercept, color = 'gray', linetype = 'dashed')
  
}

#PCA
PCAprep <- function(df){
  #For (r)opls make a wide dataframe with only numericals and subject ID as rownames, and a separate df with metadata and subject as rownames
  #ropls preprocesses the dataframe before analysis, no need to do it in advance
  
  #Data to analyse
  prepped <- df %>%
    select(sample, group, feature, value) %>% #Get necessary columns
    unite('sample_group', c(sample, group), sep = '_') %>% #Create new subject ID with health status to use to separate samples later
  
  #Turn wide
    pivot_wider(., #Turn wide
                names_from = 'feature',
                values_from = 'value') %>%
    column_to_rownames('sample_group') #Turn subject ID to rownames for opls
  
}
plotPCA <- function(df, title, prc){
  #The principal components to be plotted must be the first two columns of the df
  #prc = vector with the % the principal components contribute
  
  #Get the names of the principal component columns to be used as x and y axis
  pc_x <- names(df[1])
  pc_y <- names(df[2])
  
  #Plot PCA colored by annotation column
  plot <- ggplot(df, aes(x = !!ensym(pc_x), y = !!ensym(pc_y), color = group)) +
    geom_point() +
  
  #Set color by annotation
  scale_color_manual(values = color_manuals$sample_group,
                    guide = guide_legend(order = 1)) #+ #Set as first legend
  
  #Add title, correct labels
  plot <- plot +
    labs(title = title,
         color = 'Group',
         x = prc[1],
         y = prc[2])
  
  #Add a dashed line at x = 0, y = 0 to show center of plot
  plot %>% addLinesToPlot()
  
}
getPCpercentages <- function(prcomp_object){
  #Takes a prcomp object and calculates the percentages each PC contributes, pulls out the pc's of interest, default 1 and 2
  percentage <- round(((prcomp_object$sdev^2) / (sum(prcomp_object$sdev^2)))*100,2)
  percentage <- paste0(colnames(prcomp_object$x), " (", as.character(percentage), "%)")
}
getPCAloadings <- function(prcomp_object, principal_components) {
  #Get the loading scores for each feature
  
  #Get the feature names
  features <- as_tibble(rownames(prcomp_object$rotation)) %>%
    rename(feature = value) %>% #Call the column with names 'Item'
    rownames_to_column(., 'number')
  
  #Get the loadings and plot
  loadings <- as_tibble(prcomp_object$rotation, rownames = NA) %>% #The loadings, keep rownames
    dplyr::select(principal_components) %>% #Filter out the principal components of interest
    rownames_to_column('number') %>% #rownames are not kept when doing this, workaround by adding the features df
    inner_join(., features) %>%
    select(-number) #Remove unnecessary column
  
}
plotPCAloadings <- function(df, cutoff = 0.2, addLabels = FALSE){
  #Plot PCA loadings, df must be a dataframe with the principal components in the first two columns,
  #after getPCAloadings
  
  #Get the names of the principal component columns to be used as x and y axis
  pc_x <- names(df[1])
  pc_y <- names(df[2])
  
  #Plot the loadings colored by feature_annotation
  plot <- ggplot(df, aes(x = !!ensym(pc_x), y = !!ensym(pc_y), color = feature_anno)) + 
    scale_color_manual(values = color_manuals$feature_anno) + #Change colors
    geom_point(alpha = 0.8, size = 2)
  
  #Set names of axis ad fix color legend so it doesn't show letters
  plot <- plot + labs(
    color = '',
    x = paste0('Loadings ', pc_x),
    y = paste0('Loadings ', pc_y)
  ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  
  #Add labels/names to the features above/below the threshold
  if(isTRUE(addLabels)){
    plot <- plot +
      geom_label_repel(data = df %>% filter(df[1] > cutoff | df[2] > cutoff | df[1] < -cutoff | df[2] < -cutoff),
                       aes(label = feature), #Only show labels for the features above/below the threshold
                       size = 2.5, box.padding = unit(0.2, "lines"), show.legend = FALSE) #Change size of labels and not show legend for this layer
    
  }else if(isFALSE(addLabels)){
    plot
  } else {
    stop('addLabels must be TRUE or FALSE.')
  }
  
  plot %>% addLinesToPlot(.)  
  
}
createPCA <- function(df, title = '', principal_components = c(1,2), plots = 'scores', addLabels = FALSE, features = NULL) {
  #df = A wide df, title = title of the plot
  #principal_components = which principal components to draw, default is 1 and 2
  #plots = type of plots to make, scores, loadings, or both ('scores', 'loadings')
  
  #Do the PCA analysis
  prcomp_data <- df %>% prcomp(., scale. = FALSE)
  
  #Get the metadata
  subject_meta <- meta$samples
  
  ## Get the data needed to plot in one dataframe ##
  
  #Principal components
  plot_subjects <- as_tibble(prcomp_data$x) %>% #The principal components
    dplyr::select(principal_components) %>% #Filter out the principal components of interest
    
    #Add subject metadata and week
    mutate(sample = subject_meta$sample) %>%
    inner_join(., subject_meta)
  
  #Get percentage the different PC's contribute
  percentages <- getPCpercentages(prcomp_data)
  
  #Variables/Items contributing the most to the PC
  loadings <- getPCAloadings(prcomp_data, principal_components = principal_components) %>%
    #Add annotatons to loadings
    inner_join(., meta$features)

  #Plot principal components
  plot_subject_PCA <- plotPCA(plot_subjects, title = title, prc = percentages[principal_components])
  
  #Plot loadings
  plot_loadings <- plotPCAloadings(loadings, addLabels = addLabels)
  
  #Which plot to plot?
  if(plots == 'scores'){
    plot <- plot_subject_PCA
  } else if (plots == 'loadings') {
    plot <- plot_loadings
  } else if (plots == 'both') { #Stopped working for some reason
    plot <- plot_grid(plot_subject_PCA, plot_loadings,
                      labels = 'AUTO',
                      ncol = 2) 
  } else {
    stop("Please define which plot you want. 'scores', 'loadings' or 'both'")
  }
  
  plot %>% changeGGplotTxtSize(.)
  
}

#Do PCA and create scores and loading plots
#temp <- tidy_data %>% PCAprep(.) %>% createPCA(., plots = 'both', principal_components = c(1,2), addLabels = TRUE)
#save_plot('PCA_scores_and_loadings.png', temp, ncol = 2)

#OPLS-DA----
oplsMeta <- function(df) {
  #Get the metadata necessary for (r)opls, a df with subjectID as rownames and factors in the other columns

  #Meta
  meta <- df %>%
    select(sample, group) %>% #Get necessary columns
    unique() %>% #Remove duplicates
    
    #unite sample and group identifiers, keep group to separate
    unite('sample_group', c(sample, group), remove = FALSE) %>%
    column_to_rownames('sample_group') %>% #Turn subject ID into rownames for opls
    mutate_at('group', ~as.factor(.))
  
}
doOPLS <- function(df, separator){
  #separator = column with factor that separates the groups, such as case/control
  
  #Values to analyse
  values <- PCAprep(df)
  #Metadata
  meta <- oplsMeta(df)
  
  #Run (o)pls
  values.opls <- opls(values, meta[[separator]],
                      predI = NA, orthoI = NA,
                      permI = 1000,
                      scaleC = 'standard') #Scale = how should ropls scale the data? 'Standard' is autoscaling.
  
  values.opls
} 
getOPLSscores <- function(opls_object) {
  #Get the scores from the ropls results
  
  df <- data.frame(opls_object@scoreMN, opls_object@orthoScoreMN) %>%
    
    #Get the normal subject IDs back and get metadata
    rownames_to_column('sample_group') %>% #Get the united sample_group column back
    separate(sample_group, c('sample', 'group'), sep = '_') #split into sample and group
  
}
getOPLSpercentages <- function(opls_object){
  #Get the % of variance explained from opls components, not finished
  temp <- getOPLSscores(opls_object = opls_object) %>%
    select(matches('\\d$')) %>%
    summarise_all(sd) %>%
    pivot_longer(.,
                 cols = matches('\\d$'),
                 names_to = 'component',
                 values_to = 'sd')
}
getOPLSloadings <- function(opls_object) {
  #Get the OPLS loadings results and add type metadata
  
  loadings <- data.frame(opls_object@loadingMN, opls_object@orthoLoadingMN) %>%
    rownames_to_column('feature') %>% #Get the feature names back
    
    #Add type metadata
    inner_join(., tidy_data %>% select(feature, feature_anno) %>% unique()) %>%
    
    #Create abbreviations for type to use less space in plot
    mutate_at(., 'feature_anno', ~as.factor(.))
}
getOPLSpermutations <- function(opls_object){
  #Get the data needed to create a plot of the permutations done when building the opls model
  permutation_data <- data.frame(opls_object@suppLs$permMN) %>%
    
    #rename relevant columns
    rename(R2Y = R2Y.cum., Q2Y = Q2.cum.) %>%
    
    #Select columns relevant to plot and turn into a tidy/long df
    select(R2Y, Q2Y, sim) %>%
    pivot_longer(.,
                 cols = c(R2Y, Q2Y),
                 names_to = 'Measurement',
                 values_to = 'Value')
}
getOPLSsummary <- function(opls_object){
  #Summary statistics from ropls
  
  summary <- data.frame(opls_object@summaryDF) %>%
    pivot_longer(everything(),
                 names_to = 'statistics',
                 values_to = 'value') %>%
    
    #Combine the statistics and values to later add to plots
    mutate(combined = paste0(statistics, ':\n', value)) %>%
    
    #Clean up
    mutate_at('combined', ~str_replace_all(., '.cum.', ' cum'))
}
getOPLSvip <- function(opls_object) {
  #Get the oplsVIP scores
  
  #Get VIP and the correct colnames
  temp <- as.data.frame(opls_object@vipVn)
  colnames(temp) <- 'VIP'
  
  pVIP <- temp %>%
    rownames_to_column('feature') %>%
    inner_join(., meta$features)
}
getOPLSorthovip <- function(opls_object){
  #Get the opls ortho VIP scores by 
  
  #Get VIP and the correct colnames
  temp <- as.data.frame(opls_object@orthoVipVn)
  colnames(temp) <- 'VIP'
  
  oVIP <- temp %>%
    rownames_to_column('feature') %>%
    inner_join(., tidy_data %>% select(feature, feature_anno))
}
getOPLScorrcov <- function(opls_object, org_df) {
  #Get the covariance and correlations between the scores and the features used to calculate them
  
  #Dataframe with samples in rows and features in columns, data should be scaled
  opls_df <- PCAprep(org_df) %>% scale(.)

  #Scores
  scores <- getOPLSscores(opls_object = opls_object) %>%
    select(p1)
  
  #Correlation, use psych's corr.test and tibble to keep the names
  temp <- corr.test(scores, opls_df, ci = FALSE, method = 'kendall')
  
  corr <- as_tibble(temp$r) %>%
    pivot_longer(everything(),
                 names_to = 'feature',
                 values_to = 'Correlation')
  
  #Covariance 
  cov <- data.frame(cov(scores, opls_df)) %>%
    pivot_longer(everything(),
                 names_to = 'feature2',
                 values_to = 'Covariance')
  
  #Get original loadings with feature names
  temp <- getOPLSloadings(opls_object = opls_object) %>%
    select(p1, feature) 
  
  #Turn into one df
  corrcov <- bind_cols(corr, cov) %>%
    select(-feature2) %>%
    inner_join(., meta$features) %>% #Add feature metadata
    
    #Add VIPvn to color things by later
    inner_join(getOPLSvip(opls_object = opls_object)) %>%
    
    #Add original score back
    inner_join(., temp)
}
plotOPLS <- function(df, title = ''){
  
  df <- df %>% mutate_at('group', ~as.factor(.))
  
  #Plot the opls data with a 95% ci confidence ellipse
  plot <- ggplot(df, aes(x = p1, y = o1, colour = group)) +
    
    #Set color by group
    scale_color_manual(values = color_manuals$sample_group,
                       guide = guide_legend(order = 1)) + #Set as first legend
  
  #Set shape by different sample metadata, remove guidelines as it is messy
  #scale_shape_manual(values = color_shape_manuals$shape$PCAs, guide = 'none') + #Add shapes for gender and interaction Gender/Week
  #
  #Add fake guide with proper Gender shapes
  #scale_alpha_manual(name = 'Gender',
  #                  values = c(1,1),
  #                 guide = guide_legend(override.aes = list(shape = c(17,15)))) #color_shape_manuals$shape$PCAs)))
  #For some reason it won't read the manual without error sometimes..
  
  #Plot OPLS score, color datapoints by group
      geom_point(size = 2) + #, aes(shape = Gender)) + #Add gender as extra information
      #geom_point_interactive(tooltip = df$tooltip, size = 2) +#, aes(shape = Gender)) + #Add interactivity
      stat_conf_ellipse(inherit.aes = FALSE, aes(x = p1, y = o1, colour = group, fill = group), alpha = 0.1, geom = 'polygon') + #Add eliptical around values
      
      #Fix the legends
      guides(color = guide_legend(order = 1, override.aes = list(fill = NA, linetype = 'blank')), #Remove gray background and border of color key
             fill = FALSE,
             linetype = FALSE) #No fill/linetype legend keys
   
  #Plotting an OPLS model where measurements were taken at different weeks, have different shapes from different weeks but keep the color for separating the groups   
  #}else if(separator == 'Week'){ #Plot two weeks, differentiating between weeks with filles/empty shapes
  #  plot <- plot +
  #    geom_point(aes(shape = Week), size = 2) +#aes(shape = interaction(Week, Gender), alpha = Gender, size = Week)) + #Filled shapes for week 20, borders for week 20 colored by intervention, alpha is used to create a Gender legend
      #geom_point_interactive(tooltip = df$tooltip, aes(shape = interaction(Week, Gender)), size = 1) + #Add interactivity
      
      #Add  95% ci elipticals around values
  #    stat_conf_ellipse(inherit.aes = FALSE, aes(x = p1, y = o1,
  #                                               colour = Intervention, linetype = Week, fill = Intervention),
  #                      alpha = 0.1, size = 0.2, geom = 'polygon') +
      
      #Add a filled shape for week 20, an empty for week 0 using an unneccessary/empty aestetics
  #    scale_size_manual(values = c(2, 2), #Set size
  #                      guide = guide_legend(order = 2, #Set as second legend
   #                                          override.aes = list(shape = c(1,16)))) + #Change into empty/filled circles
      #Fix the legends
    #  guides(color = guide_legend(order = 1, override.aes = list(fill = NA, linetype = 'blank')),
     #        linetype = FALSE,#guide_legend(order = 3, 
             #            override.aes = list(name ='', colour = 'black', fill = NA, shape = c(1,16))), #Color must be added here as linetypes by default are added to color legend and are empty in their own
      #       fill = FALSE) + #No fill or lines in the color legend
      #scale_linetype_manual(name = '',
      #                     values = c('solid', 'longdash')) + #Remove second 'Week' from legend
      #scale_shape_manual(values = c(1, 16))
  #}
  
  #Add title, correct labels
  plot <- plot +
    labs(title = title,
         colour = 'Group',
         x = 'OPLS-DA predictive component',
         y = 'Orthogonal component') 
  
  #Add lines to show origo
  plot %>% addLinesToPlot(.)
  
}
plotOPLSloadings <- function(df, cutoff = 0.08, addLabels = FALSE, interactive = FALSE){
  #Plot OPLS loadings, df must be a dataframe with p1 and o1 in the second and third column respectively
  
  #Get the names of the principal component columns to be used as x and y axis
  pc_x <- names(df[2])
  pc_y <- names(df[3])
  
  #Plot loadings colored by type
  plot <- ggplot(df, aes(x = !!ensym(pc_x), y = !!ensym(pc_y), color = feature_anno)) +
    scale_color_manual(values = color_manuals$feature_anno) #Change colors
  
  #Use geom point or geom point interactive based in if interactive is true or false
  if(isFALSE(interactive)){
    plot <- plot + geom_point(alpha = 0.8, size = 3)
  } else if (isTRUE(interactive)) {
    plot <- plot + geom_point_interactive(tooltip = df$tooltip, alpha = 0.8, size = 3)
  } else {
    stop('Interactive must be TRUE or FALSE')
  }
  
  #Set names of axis ad fix color legend so it doesn't show letters
  plot <- plot + labs(
    color = '',
    x = paste0('Loadings ', pc_x),
    y = paste0('Loadings ', pc_y)
  ) +
    guides(color = guide_legend(override.aes = list(alpha = 1, size = 3)))
  
  #Add labels/names to the features above/below the threshold
  if(isTRUE(addLabels)){
    plot <- plot +
      geom_label_repel(data = df %>% filter(df[2] > cutoff | df[3] > cutoff | df[2] < -cutoff | df[3] < -cutoff),
                       aes(label = feature), #Only show labels for the features above/below the threshold
                       size = 2.5, box.padding = unit(0.2, "lines"), show.legend = FALSE) #Change size of labels and not show legend for this layer
    
  }else if(isFALSE(addLabels)){
    plot
  } else {
    stop('addLabels must be TRUE or FALSE.')
  }
  
  plot %>% addLinesToPlot(.)
  
}
plotOPLSpermutations <- function(df, title = '') {
  #Create a permutation plot like from ropls
  #df = dataframe with the permutation data from opls model
  
  #Build plot
  plot <- ggplot(df, aes(x = sim, y = Value, colour = Measurement)) +
    
    #Reduce alpha of the overlapping points
    geom_point(alpha = 0.7) +
    
    #Add a horizontal line to show the original R2 and Q2 values
    geom_hline(data = df %>% filter(sim == 1 & Measurement == 'R2Y'), aes(yintercept = Value)) +#, colour = color_manual$permutations[1]) +
    geom_hline(data = df %>% filter(sim == 1 & Measurement == 'Q2Y'), aes(yintercept = Value)) +#, colour = color_manual$permutations[2]) +
    
    #Add legend in lower left corner
    theme(legend.position = c(0.85, 0.2),
          legend.title = element_blank(),
          legend.spacing.y = unit(0.05, 'cm'),
          
          #Add a white background with black lines
          legend.background = element_rect(fill="white",
                                           size=0.3, linetype="solid", 
                                           colour ="black")) +
    
    #Add title and correct axis text
    labs(
      title = title,
      colour = 'Stat', #change legend key title
      x = 'Correlation coefficient between original and permuted data',
      y = 'R2Y and Q2Y'
    ) #+
    
    #Change colors and shapes of points
    #scale_color_manual(values = color_manual$permutations)
  
}
plotOPLSsummary <- function(df){
  #Create a white plot element with the summary statistics of ropls opls model
  
  plot <- ggplot(df, aes(x = 1:length(combined), y = 0, label = combined)) +
    geom_text(size = 3.5) +
    theme_void() #Set theme with no axes or anything
  
}
plotOPLSsplot <- function(df, cov_cutoff = NULL, corr_cutoff = NULL, VIP_cutoff = NULL, addLabels = FALSE, interactive = FALSE) {
  #Creates an s-plot for the opls-da model
  #df = corrcov df
  
  #Add tooltip for ggiraph and mutate feature_anno so features can be marked 'Below threshold'
  df <- df %>% 
    #mutate(tooltip = addFeatureMetadataTooltip(df)) %>%
    mutate_at('feature_anno', ~as.character(.)) 
  
  #Turn every point below threshold to gray color, either by cov/corr or by VIP
  if(!is.null(corr_cutoff) & is.null(cov_cutoff)){
    df <- df %>%
      mutate(feature_anno = case_when(Correlation < corr_cutoff & Correlation > -corr_cutoff ~ 'Below threshold',
                              TRUE ~ feature_anno)) %>%
      mutate_at('feature_anno', ~as.factor(.))}
  else if(!is.null(cov_cutoff) & is.null(corr_cutoff)){
    df <- df %>%
      mutate(feature_anno = case_when(Covariance < cov_cutoff & Covariance > -cov_cutoff ~ 'Below threshold',
                              TRUE ~ feature_anno)) %>%
      mutate_at('feature_anno', ~as.factor(.))
  } else if(!is.null(cov_cutoff) & !is.null(corr_cutoff)){
    df <- df %>%
      mutate(feature_anno = case_when(Covariance < cov_cutoff & Covariance > -cov_cutoff ~ 'Below threshold',
                              Correlation < corr_cutoff & Correlation > -corr_cutoff ~ 'Below threshold',
                              TRUE ~ feature_anno)) %>%
      mutate_at('feature_anno', ~as.factor(.))
  }
  else if (!is.null(VIP_cutoff)) {
    df <- df %>%
      mutate(feature_anno = case_when(VIP < VIP_cutoff & VIP > -VIP_cutoff ~ 'Below threshold',
                              TRUE ~ feature_anno)) %>%
      mutate_at('feature_anno', ~as.factor(.))
  }
  
  #Turn feature_anno into factor and relevel so 'Below threshold' is the last
  df <- df %>%
    mutate_at('feature_anno', ~as.factor(.)) %>%
    mutate_at(., 'feature_anno', ~fct_relevel(., 'Below threshold', after = Inf))
  
  #Plot
  plot <- ggplot(df, aes(x = Covariance, y = Correlation, color = feature_anno)) +
    scale_color_manual(values = color_manuals$feature_anno) + #Change colors
    geom_point() +
    
    #Set axis labs
    labs(
      x = 'Covariance with group label',
      y = 'Correlation with group label'
    ) +
    
    #Increase size of color keys
    guides(colour = guide_legend(override.aes = list(size = 4)))
  
  #Add labels to important points or not
  if(isTRUE(addLabels)){
    plot <- plot + geom_label_repel(data = df %>% filter(feature_anno != 'Below threshold'), aes(label = feature), #Only show labels for the features above/below the threshold
                                    size = 2, box.padding = unit(0.2, "lines"), show.legend = FALSE) + #Change size of labels and not show legend for this layer) 
      labs(
        color = 'feature_anno')}
  else if (isFALSE(addLabels)){
    plot <- plot
  }else{
    stop("Do you want to add labels to the plot? TRUE/FALSE")
  }
  
  #Add interactivity or not
  #if(isTRUE(interactive)){
  #  plot <- plot + geom_point_interactive(tooltip = df$tooltip)
  #}else if(isFALSE(interactive)){
  #  plot <- plot
  #}else{
  #  stop("Should the plot have interactive points? TRUE/FALSE")
  #}
  
  plot %>% addLinesToPlot()
  
}
plotALLopls <- function(opls_object, separator, org_df, addLabels = FALSE, interactive = FALSE, cov_cutoff = NULL, corr_cutoff = NULL, VIP_cutoff = NULL) {
  #Plot the OPLS score plot, loadings and s-plot, opls_object = ropls results, separator = the groups used running ropls
  
  #Get the scores, loadings, covariance/correlations, permutations and model summary
  scores <- getOPLSscores(opls_object = opls_object)
  loadings <- getOPLSloadings(opls_object = opls_object)
  permutations <- getOPLSpermutations(opls_object = opls_object)
  corrcov <- getOPLScorrcov(opls_object = opls_object, org_df = org_df)
  summary <- getOPLSsummary(opls_object = opls_object)
  
  #Plot
  scores_plot <- plotOPLS(scores) %>%
    changeGGplotTxtSize()
  loadings_plot <- plotOPLSloadings(loadings, addLabels = addLabels, interactive = interactive) %>%
    changeGGplotTxtSize()
  permutations_plot <- plotOPLSpermutations(permutations) %>%
    changeGGplotTxtSize()
  splot <- plotOPLSsplot(corrcov, cov_cutoff = cov_cutoff, corr_cutoff = corr_cutoff, VIP_cutoff = VIP_cutoff) %>%
    changeGGplotTxtSize()
  summary_plot <- plotOPLSsummary(summary)
  
  #Plot scores and loadings together first
    scores_loadings <- plot_grid(
      
      #Score plot
      scores_plot + theme( 
        
        #Set the legend of the OPLS-score plot in the lower left inner corner
        legend.position = c(.27, .07),
        legend.box = 'vertical',
        legend.direction = 'horizontal',
        legend.margin = margin(-7, -6, 0, 0),
        legend.spacing.y = unit(0.07, 'cm'),
        
        #Add a white background with black lines
        legend.background = element_rect(color="white")), 
      
      
      #Loading plot
      loadings_plot + theme(legend.position = 'none'),
      rel_widths= c(1,1),
      labels = 'AUTO',
      axis = 'l',
      align = 'hv',
      ncol = 2)
  
  #Add s-plot and summary on the bottom
  sl_splot <- plot_grid(scores_loadings,
                        
                        #Fix the legend position
                        splot + theme(
                          legend.box.spacing = unit(0.2, 'cm')
                        ),
                        summary_plot,
                        nrow = 3,
                        labels = c('', 'C'),
                        axis = 'l',
                        rel_heights = c(1,2,0.16))
  
}  
organizeOPLSvips <- function(opls_object, org_df = tidy_data){
  #Gets the VIP, covariance and correlation scores and create a table showing how many VIPs >1 are either low or high in the different types of features
  
  #Get VIPs
  VIPs <- getOPLSvip(opls_object = opls_object)
  
  #Get the scores
  df <- getOPLScorrcov(opls_object = opls_object, org_df = org_df) %>%
    filter(VIP > 1) %>%
    
    #Coun the number of increased or decreased features in case vs control by feature annotation
    #group_by(feature_anno) %>%
    #mutate(`Increased` = sum(Covariance < 0)) %>%
    #mutate(`Decreased` = sum(Covariance > 0)) %>%
    #ungroup() %>%
    #select(feature, feature_anno, contains('creased'), contains('total'), Correlation) %>% unique() %>%
    
    inner_join(., VIPs) %>%
    
    #Organize
    select(feature_anno, feature, VIP, Correlation)
  
}

#Do OPLS
opls <- tidy_data %>% doOPLS(., 'group')

#Plots and tables
#Scores and loadings
plot1 <- getOPLSscores(opls) %>% plotOPLS(.)
plot2 <- getOPLSloadings(opls) %>% plotOPLSloadings(., addLabels = FALSE)
save_plot('OPLS_scores_and_loadings.png', plot_grid(plot1, plot2, ncol = 2), nrow = 2, ncol = 4)

#S-plot
temp <- getOPLScorrcov(opls, tidy_data) %>% plotOPLSsplot(.)

#Large plot with everything
temp <- plotALLopls(opls, separator = 'group', org_df = tidy_data, VIP_cutoff = 1)

#VIP and correlations table, correlation is with cases
temp <- organizeOPLSvips(opls)
