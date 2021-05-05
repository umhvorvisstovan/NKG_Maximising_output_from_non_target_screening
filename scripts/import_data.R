# importing data

library(tidyverse)
library(readxl)

effluent_features <- read_csv("data/second_report_effluent_features.csv", col_types = cols())
effluent_results <- read_csv("data/second_report_effluent_results.csv", col_types = cols())
effluent_substances <- read_csv("data/second_report_effluent_unique_substances.csv", col_types = cols())

sample_info <- read_csv("data/first_report_sample_info.csv", col_types = cols())

initial_results <- read_csv("data/first_report_results.csv", col_types = cols())
initial_substances <- read_csv("data/first_report_unique_substances.csv", col_types = cols())

gcgc <- read_csv("data/second_report_requantified_biota_gcgc.csv", col_types = cols())


# table with identified features in effluent samples
identified <- effluent_features %>% 
  filter(IDLevel < 3) %>% 
  filter(is.na(Comment)) %>% 
  distinct(ID, IDLevel, Name) %>% 
  left_join(effluent_substances, by = "Name")

# effluent decon id table
effluent_table <- identified %>% 
  select(ID, CAS, Name, Norman_SusDat_ID, IDLevel) %>% 
  distinct_all() %>% 
  left_join(effluent_results %>% select(ID, sample_no), by = "ID") %>% 
  left_join(sample_info %>% select(sample_no, country) %>% distinct_all(), by = "sample_no") %>% 
  select(-ID, -sample_no)

gcgc_table <- gcgc %>% 
  left_join(sample_info %>% select(sample_no, country) %>% distinct_all(), by = "sample_no") %>% 
  select(-sample_no) %>% 
  left_join(initial_substances %>% select(Name, CAS, Molecular_Formula, Norman_SusDat_ID), by = "Name") %>% 
  mutate(Molecular_Formula = if_else(is.na(Molecular_Formula), MolecularFormula, Molecular_Formula),
         Molecular_Formula = str_remove_all(Molecular_Formula, "-"))


initial_table <- initial_results %>% 
  select(-CAS) %>% 
  left_join(sample_info %>% select(sample_no, country, matrix) %>% distinct_all(), by = "sample_no") %>% 
  select(-sample_no) %>% 
  left_join(initial_substances %>% select(Name, CAS, Molecular_Formula, Norman_SusDat_ID,
                                          CompoundInfo, DrugType, PotentialSource), by = "Name") %>% 
  mutate(Molecular_Formula = if_else(is.na(Molecular_Formula), MolecularFormula, Molecular_Formula),
         Molecular_Formula = str_remove_all(Molecular_Formula, "-"))

# different between the initial effluent data and the new effluent data

initial_effluent <- initial_results %>%
  left_join(sample_info, by = "sample_no") %>% 
  filter(matrix == "Effluent") %>% 
  distinct(Name) %>%  
  mutate(study = "initial")

second_effluent <- effluent_features %>%  
  distinct(Name) %>% 
  mutate(study = "second")

effluent_substances_only_initial <- initial_effluent %>% 
  filter(!(Name %in% effluent_features$Name))

effluent_substances_level2_only_second <- effluent_features %>%
  filter(IDLevel < 3) %>% 
  distinct(Name) %>%
  filter(!(Name %in% initial_effluent$Name))
