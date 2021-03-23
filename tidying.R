#Set working directory
setwd('C:/Users/aslau/Documents/Uni/Master/gitHub/MECFS')

#Load libraries
library(tidyverse)
library(broom)
library(stringi)
library(readxl)
library(janitor)
library(funModeling)

#Goal
#Tidy dataframe with a 'feature' column containing features, a 'feature_anno'
#explaining what the features are, a 'sample_id' column containing the samples,
#a 'value' column containing the values of the features for the samples.
#Low variance or high missing number features should be discarded.

#Read data----
#2020 dataset
data2020 <- excel_sheets('germain2020.xlsx') %>% #Get names of sheets
  set_names() %>%
  map(read_excel,
      path = 'germain2020.xlsx')
#List for various df's to keep environment clean
various <- list()

#Get to know the data----
#No info on the samples

#Missing values, low variance metabolites in control and cases
#wide df that can be filtered based on case/control status
temp <- t(as_tibble(data2020$data_matrix)) %>% #As tibble to keep row/columnnames
  row_to_names(., 1) %>% #Turn compound_id to column names
  as_tibble(rownames = NA) %>% #Turn into tibble, keep sample_id rownames
  rownames_to_column(., 'sample_id') %>% #'sample_id' is used in later functions to analyze the data
  inner_join(., data2020$sample_metadata %>% select(-sample_id_internal)) %>%
  rename(group = health_status) #group is used in later functions

#Get info on missing values and low variance metabolites in controls
various$missing_and_zeros_control <- temp %>%
  filter(group == 'control') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/sum(temp$group == 'control')*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5 & unique > 1) %>%
  #Rename variable to compound_id
  rename(compound_id = variable)
#Cases
various$missing_and_zeros_cases <- temp %>%
  filter(group == 'case') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/sum(temp$sample_id == 'case')*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5 & unique > 1) %>%
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
  inner_join(data2020$data_dictionary) #Add compound metadata

various$removed_lowVar_cases <- temp %>%
  pivot_longer(cols = -sample_id,
               names_to = 'compound_id',
               values_to = 'value') %>%
  select(compound_id) %>%
  unique() %>%
  anti_join(., various$missing_and_zeros_cases) %>%
  inner_join(data2020$data_dictionary) #Add compound metadata

#Remove variables with <5% variance, turn into Tidy format----
temp <- full_join(various$removed_lowVar_controls, various$removed_lowVar_cases) #Metabolites with <5% variance in either group

tidy_data <- data2020$data_matrix %>%
  pivot_longer(.,
               cols = -compound_id,
               names_to = 'sample_id') %>%
  anti_join(., temp %>% select(compound_id)) %>% #Remove metabolites with <5% variance in either group

  #Add sample metadata (case/control status)
  inner_join(., data2020$sample_metadata %>% select(-sample_id_internal)) %>%

  #get compound names, not just id, and some metadata
  inner_join(., data2020$data_dictionary %>% select(compound_id, BIOCHEMICAL, SUPER_PATHWAY, SUB_PATHWAY)) %>%

  #Change column names to work with other functions later,
  #feature names column is named 'feature', if there is any column that can be used to annotate features call it feature_anno
  #and the column that show which group participants are in called 'group'
  rename(feature = BIOCHEMICAL,
         feature_anno = SUPER_PATHWAY,
         group = health_status) %>%

  #Change names of control/cases
  mutate(group = case_when(
    group == 'case' ~ 'ME/CFS',
    group == 'control' ~ 'Control'
  ))

saveRDS(tidy_data, 'tidy_data2020.Rds')

#2018 dataset
#Read data
data2018 <- read_xlsx('germain2018.xlsx')

tidy_data <- data2018 %>%

  #Remove unnecessary columns
  select(-c(`PATHWAY SORTORDER`, `COMP ID`, `PLATFORM`, `CHEMICAL ID`, RI, MASS, CAS, PUBCHEM, KEGG, HMDB)) %>%

  #Tidy
  pivot_longer(.,
               cols = -c(BIOCHEMICAL, `SUPER PATHWAY`, `SUB PATHWAY`),
               names_to = 'sample_id',
               values_to = 'value') %>%
  #Rename BIOCHEMICAL to feature for plotting functions
  rename(feature = BIOCHEMICAL, feature_anno = `SUPER PATHWAY`) %>%
  #Create a column for control/sample status
  mutate(group = case_when(
    str_detect(sample_id, 'C') ~ 'Control',
    str_detect(sample_id, 'P') ~ 'ME/CFS'
  ))

#Get to know the data
#wide df that can be filtered based on case/control status
temp <- tidy_data %>%
  select(-c(feature_anno, `SUB PATHWAY`)) %>%
  pivot_wider(names_from = feature,
              values_from = value)

#Get info on missing values and low variance metabolites in controls
various$missing_and_zeros_control <- temp %>%
  filter(group == 'Control') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/sum(temp == 'Control')*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5 & unique > 1) %>%
  #Rename variable to compound_id
  rename(feature = variable)
#Cases
various$missing_and_zeros_cases <- temp %>%
  filter(group == 'ME/CFS') %>%
  df_status(.) %>%
  #Percentage of unique values for each metabolite
  mutate(p_unique = unique/sum(temp == 'ME/CFS')*100) %>%
  #Filter out metabolites with <5% unique values
  filter(p_unique > 5 & unique > 1) %>%
  #Rename variable to compound_id
  rename(feature = variable)

#Compounds removed from either groups in last step
various$removed_lowVar_controls <- tidy_data %>%
  anti_join(., various$missing_and_zeros_control) %>%
  select(feature) %>%
  unique()

various$removed_lowVar_cases <- tidy_data %>%
  anti_join(., various$missing_and_zeros_cases) %>%
  select(feature) %>%
  unique()
temp <- full_join(various$removed_lowVar_controls, various$removed_lowVar_cases) #Metabolites with <5% variance in either group

#Remove from final df
tidy_data <- tidy_data %>%
  anti_join(temp)

saveRDS(tidy_data, 'tidy_data2018.Rds')
