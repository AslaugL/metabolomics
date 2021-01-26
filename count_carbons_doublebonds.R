#Count number of carbons and double bonds in lipids from Biocrates or zora biosciences for plotting

#Set working directory
setwd('C:/Users/aslau/Documents/Uni/Master/gitHub')

#Load libraries
library(tidyverse)
library(stringi)

#Read data
data <- readRDS('tidy_data.Rds')

#Function
countCarbonsDoublebonds <- function(df){
  #Count the number of carbons and double bonds in lipid species, the lipid species must be in the first column of df
  
  counts <- df %>% 
    
    #Count carbonds and double bonds
    mutate(temp_carbons = str_extract_all(pull(., var = 1), '\\d{1,2}(?=:)')) %>%
    mutate(temp_double_bonds = str_extract_all(pull(., var = 1), '(?<=:)\\d{1}')) %>%
    
    #If there are multiple numbers, it will create a list column that must be unnested before turned into numberics for summation
    unnest(cols = c(temp_carbons, temp_double_bonds)) %>%
    mutate_at(c('temp_carbons', 'temp_double_bonds'), ~as.numeric(as.character(.))) %>%
    group_by(pull(., var = 1)) %>%
    mutate(carbons = sum(temp_carbons), double_bonds = sum(temp_double_bonds)) %>%
    ungroup(.) %>%
    
    #remove temp columns and duplicates created by unnest
    select(-c(temp_carbons, temp_double_bonds, `pull(., var = 1)`))  %>%
    unique()
}
getType <- function(df){
  #Get the type of lipid (DG, TG, LPC, AC etc) in one column for plotting
  
  with_type <- df %>%
    separate(col = 1, into = c('Type', 'Numbers'), sep = "-(?!\\w)| |\\(", remove = FALSE) %>% #Separate type of lipids from the numbers
    mutate_at('Numbers', ~str_replace(., '\\)', '')) %>% #Remove unnecessary characters
    
    #Get the tag such as -OH, -DH for some lipds
    separate(col = Numbers, into = c('Numbers', 'Tag'), sep = '-(?=\\w{1,2})') %>%#This creates an 'NA' in lipids with no tag
    unite(col = Type, c(Type, Tag), sep = '-') %>%
    #Add tag to type
    mutate_at('Type', ~str_replace(., '-NA', '')) %>% #Remove the NA tag
    
    #Remove unnecessary columns
    select(-Numbers)
}
changeGGplotTxtSize <- function(plot, txt_size = 12){
  #Change the text size of the axis title and text of a ggplot
  plot + theme(axis.title = element_text(size = txt_size),
               axis.text = element_text(size = txt_size),
               legend.title = element_text(size = txt_size),
               legend.text = element_text(size = txt_size))
}
plotCarbonsDoublebondsLipids <- function(df) {
#Creates a heatmap of logfold change of the lipids with carbon length on y-axis and double_bonds on x-axis
  #df = output from countCarbonsDoublebonds and getType
  
  #Temporary plot
  plot <- ggplot(df, aes(x = double_bonds, y = as.factor(carbons))) +
    geom_tile(aes(fill = `log2(FC)`)) +
    geom_text(aes(label = `log2(FC)`), color = 'white', size = 3, show.legend = FALSE) +
  
  #Facet by type of lipid
   facet_grid(rows = vars(Type),
              scales = 'free',
              space = 'free') +
  
  #Set color and axis scales
    scale_fill_continuous(type = 'viridis') +
    scale_x_continuous(breaks = unique(df$double_bonds)) +
    
    #Add borders around the plot and have the names of the grid elements horizontal
    theme_light() +
    theme(strip.text.y = element_text(angle = 0),
          strip.text = element_text(color = 'black'),
          strip.background = element_rect(color = 'white', fill= 'white'),
          legend.position = 'right',
          legend.key.width = unit(0.5, 'cm')) + #Reduce width of legend key
    
    #Fix axis, legend title and keys
    labs(
      y = 'Total carbons',
      x = 'Double bonds'
    ) +
    guides(fill = guide_colorbar(title = 'Log2(FC)',
                                 label.position = 'left', 
                                 label.hjust = 1))
  
  #Change txt size
  plot <- plot %>% changeGGplotTxtSize(txt_size = 10)
}


#Sample some of the data for a plot
temp <- data %>%
  select(compound, `log2(FC)`) %>%
  filter(str_detect(compound, 'CE\\(|DAG|LPC')) %>%
  unique() %>%
  
  #Count carbons and double bonds and get type of lipid
  countCarbonsDoublebonds() %>% getType() %>%
  
  #Some lipid molecules have multiple configurations of fatty acids with similar carbon/double bond number, take average fold change
  select(-compound) %>%
  group_by(Type, carbons, double_bonds) %>%
  mutate(meanLog2 = mean(`log2(FC)`)) %>%
  select(-`log2(FC)`) %>%
  rename(`log2(FC)` = meanLog2) %>% unique() %>%
  #Round for easier reading
  mutate(`log2(FC)` = round(`log2(FC)`, 2))

plot <- plotCarbonsDoublebondsLipids(temp)
plot

ggsave('LipidCarbonDoublebonds.png')
