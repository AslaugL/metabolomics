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

#Read data----
data <- excel_sheets('germain2020.xlsx') %>% #Get names of sheets
  set_names() %>%
  map(read_excel,
      path = 'germain2020.xlsx')
#List for various df's to keep environment clean
various <- list()

#Get to know the data----
#No info on the samples

#Missing values, low variance metabolites in control and cases
#wide df that can be filterest based on case/control status
temp <- t(as_tibble(data$data_matrix)) %>% #As tibble to keep row/columnnames
  row_to_names(., 1) %>% #Turn compound_id to column names
  as_tibble(rownames = NA) %>% #Turn into tibble, keep sample_id rownames
  rownames_to_column(., 'sample_id') %>%
  inner_join(., data$sample_metadata %>% select(-sample_id_internal))

#Get info on missing values and low variance metabolites in controls
various$missing_and_zeros_control <- temp %>%
  filter(health_status == 'control') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/52*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5) %>%
  #Rename variable to compound_id
  rename(compound_id = variable)
#Cases
various$missing_and_zeros_cases <- temp %>%
  filter(health_status == 'case') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/52*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5) %>%
  #Rename variable to compound_id
  rename(compound_id = variable)

#Compounds removed from either groups in last step
various$removed_lowVar_controls <- temp %>% 
  pivot_longer(cols = -sample_id,
               names_to = 'compound_id',
               values_to = 'value') %>%
  select(compound_id) %>%
  unique() %>%
  anti_join(., various$missing_and_zeros_control) %>%
  inner_join(data$data_dictionary) #Add compound metadata

various$removed_lowVar_cases <- temp %>% 
  pivot_longer(cols = -sample_id,
               names_to = 'compound_id',
               values_to = 'value') %>%
  select(compound_id) %>%
  unique() %>%
  anti_join(., various$missing_and_zeros_cases) %>%
  inner_join(data$data_dictionary) #Add compound metadata

#Remove variables with <5% variance, turn into Tidy format----
temp <- full_join(various$removed_lowVar_controls, various$removed_lowVar_cases) #Metabolites with <5% variance in either group

tidy_data <- data$data_matrix %>%
  pivot_longer(.,
               cols = -compound_id,
               names_to = 'sample_id') %>%
  anti_join(., temp %>% select(compound_id)) %>% #Remove metabolites with <5% variance in either group 
  
  #Add sample metadata (case/control status)
  inner_join(., data$sample_metadata %>% select(-sample_id_internal)) %>%
  
  #get compound names, not just id, and some metadata
  inner_join(., data$data_dictionary %>% select(compound_id, BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY)) %>%
  rename(compound = BIOCHEMICAL)

