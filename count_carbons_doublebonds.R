#Count number of carbons and double bonds in lipids from Biocrates or zora biosciences for plotting

#Set working directory
#setwd('C:/Users/Aslaug L/Desktop/Master')
setwd('C:/Users/aslau/Documents/Uni/Master/gitHub')

#Load libraries
library(tidyverse)
library(stringi)

#Read data
data <- readRDS('lipids_BiocratesZorabioscience.Rds')

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
