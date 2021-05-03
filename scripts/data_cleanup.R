# getting CAS numbers for IDs

library(tidyverse)
library(readxl)

sample_info <- read_excel("data/raw_data/originals/samples.xlsx")
write_csv(sample_info %>% arrange(sample_no), "data/first_report_sample_info.csv")

# ID score equal ID levels
df_idlevel <- data.frame(IDScore = c(-1, 0, 1, 7),
                         IDLevel = c(3.5, 5, 2, NA))

featureIDs <- read_excel("data/raw_data/originals/Aligned_Peak_List_WithID_Extended_2709018_BA.xlsx", sheet = "IDs") %>% 
  select(- ...1) %>% 
  #clean Name column
  mutate(Name_original = Name,
         Name = str_remove(Name, "\\.$"), #deletes trailing .
         Name = str_remove(Name, "\\;$")) %>% #deletes trailing ;
  rename(IDScore = `ID Score`) %>% 
  left_join(df_idlevel, by = "IDScore") %>% 
  select(-IDScore)

effluent_substances <- featureIDs %>% 
  distinct(Name_original, Name) %>% 
  filter(!is.na(Name_original)) %>% 
  filter(Name_original != "unknown") %>% 
  filter(Name_original != "Unknown") %>% 
  arrange(Name_original)

normansusdat <- read_csv("data/raw_data/originals/susdat_2020-06-18-123341_20210223.csv", 
                         col_types = cols(`Lowest_P-PNEC_(QSAR)_ug/L` = col_skip(),
                                          ChemSpiderID = col_skip(),
                                          ChemSpider_ID_based_on_InChIKey_19032018 = col_skip(),
                                          alogp_ChemSpider = col_skip(),
                                          xlogp_ChemSpider = col_skip())) %>% 
  distinct_all()
normansusdat_downloaded <- read_delim("data/raw_data/originals/norman_website.txt", delim = "\t", col_types = cols(PubChem_CID = col_character())) %>% 
  distinct_all() %>% 
  filter(!Norman_SusDat_ID %in% normansusdat$Norman_SusDat_ID)
normansusdat <- normansusdat %>% 
  bind_rows(normansusdat_downloaded)
rm(normansusdat_downloaded)

# clean the norman susdat list, take out duplicate with less information

# the IDs to keep
normansusdat_replicate_to_keep <- normansusdat %>%
  group_by(Name) %>%
  filter(n()>1) %>%
  ungroup() %>%
  mutate(info = rowSums (!is.na(.[,1:52]))) %>% 
  select(Norman_SusDat_ID, Name, info) %>% 
  group_by(Name) %>%
  arrange(desc(info)) %>% 
  mutate(replicates = paste0(Norman_SusDat_ID, collapse = ","),
         replicates = word(replicates, start = 2, end = -1, sep = ",")) %>% 
  top_n(1, info) %>% 
  select(-info)

# the IDs to delete
normansusdat_replicate_to_delete <- normansusdat %>%
  group_by(Name) %>%
  filter(n()>1) %>%
  ungroup() %>%
  filter(!(Norman_SusDat_ID %in% normansusdat_replicate_to_keep$Norman_SusDat_ID)) %>% 
  select(Norman_SusDat_ID)

# filter out the Norman_SusDat_ID that are replicates
normansusdat <- normansusdat %>% 
  filter(!(Norman_SusDat_ID %in% normansusdat_replicate_to_delete$Norman_SusDat_ID)) 

# all possible name synonyms in susdat
normannames <- bind_rows(
  normansusdat %>% 
    pivot_longer(cols = c(Name, Name_Dashboard, Name_ChemSpider, Name_IUPAC),
                 values_to = "Name", names_to = "Name_Source") %>% 
    filter(!is.na(Name)) %>% 
    distinct(Norman_SusDat_ID, Name_Source, Name, CAS_RN_Dashboard), 
  normansusdat %>%
    mutate(count_syn = str_count(Synonyms_ChemSpider, ";"),
           count_rel = str_count(Reliability_of_Synonyms_ChemSpider, ";"),
           diff = if_else(count_syn == count_rel, TRUE, FALSE)) %>% 
    filter(diff == TRUE) %>% 
    separate_rows(Synonyms_ChemSpider, Reliability_of_Synonyms_ChemSpider, sep = ";") %>% 
    filter(!is.na(Synonyms_ChemSpider)) %>% 
    distinct(Norman_SusDat_ID, Synonyms_ChemSpider, Reliability_of_Synonyms_ChemSpider, CAS_RN_Dashboard) %>% 
    rename(Name = Synonyms_ChemSpider) %>% 
    mutate(Name_Source = "Synonyms_ChemSpider") 
) %>% 
  arrange(Norman_SusDat_ID) %>% 
  filter(!is.na(Name)) %>% 
  mutate(Reliability_of_Synonyms_ChemSpider = as.numeric(Reliability_of_Synonyms_ChemSpider))

# match validation ranking
namerank <- data.frame(Name_Source = c("Name", "Name_Dashboard", "Name_ChemSpider", "Name_IUPAC", "Synonyms_ChemSpider"),
                       Namerank = c(1:5))

normannames <- normannames %>% 
  left_join(namerank, by = "Name_Source")

# all possible CAS in susdat
normanCASs <- bind_rows(
  normansusdat %>% 
    pivot_longer(cols = c(CAS_RN, CAS_RN_Cactus, CAS_RN_Dashboard, CAS_RN_PubChem),
                 values_to = "CAS", names_to = "CAS_Source") %>% 
    filter(!is.na(CAS)) %>% 
    mutate(CAS = str_remove_all(CAS, "CAS_RN: ")) %>% 
    mutate(CAS = str_remove_all(CAS, "NOCAS_")) %>%
    distinct(Norman_SusDat_ID, CAS_Source, CAS) %>% 
    separate_rows(CAS, sep =";"), 
  normansusdat %>%
    mutate(count_syn = str_count(CAS_RN_ChemSpider, ";"),
           count_rel = str_count(Reliability_of_CAS_ChemSpider, ";"),
           diff = if_else(count_syn == count_rel, TRUE, FALSE)) %>% 
    filter(diff == TRUE) %>% 
    separate_rows(CAS_RN_ChemSpider, Reliability_of_CAS_ChemSpider, sep = ";") %>% 
    filter(!is.na(CAS_RN_ChemSpider)) %>% 
    distinct(Norman_SusDat_ID, CAS_RN_ChemSpider, Reliability_of_CAS_ChemSpider) %>% 
    rename(CAS = CAS_RN_ChemSpider) %>% 
    mutate(CAS_Source = "CAS_RN_ChemSpider") 
) %>% 
  arrange(Norman_SusDat_ID) %>% 
  filter(!is.na(CAS)) %>% 
  mutate(Reliability_of_CAS_ChemSpider = as.numeric(Reliability_of_CAS_ChemSpider))

# match validation ranking
CASrank <- data.frame(CAS_Source = c("CAS_RN", "CAS_RN_Dashboard", "CAS_RN_PubChem", "CAS_RN_Cactus", "CAS_RN_ChemSpider"),
                      CASrank = c(1:5))

normanCASs <- normanCASs %>% 
  left_join(CASrank, by = "CAS_Source")

# match to all possible names in SusDat
IDs_matchnames <- effluent_substances %>% 
  left_join(normannames %>% filter, by = "Name") %>%
  filter(!is.na(Norman_SusDat_ID)) %>% 
  mutate(syn_rel = if_else(is.na(Reliability_of_Synonyms_ChemSpider), 10, Reliability_of_Synonyms_ChemSpider)) %>% 
  arrange(Namerank, syn_rel) %>% 
  group_by(Name) %>% 
  slice_head(n = 1) %>% 
  mutate(Norman_matched = Name_Source,
         Norman_matched_Synonym_Reliability = Reliability_of_Synonyms_ChemSpider) %>% 
  distinct(Name_original, Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability)

# missing Norman_SusDat_ID
missing <- effluent_substances %>%
  filter(!(Name %in% IDs_matchnames$Name)) %>% 
  .$Name

# batch search on CompTox -  https://comptox.epa.gov/dashboard
# split the missing data set in two because of a max 5000 line search on the comptox dashboard
comtox <- read_csv(file = "data/raw_data/CompToxChemicalsDashboard-Batch-Search_2021-04-25_10_43_44.csv", col_types = cols(), na = "-") %>% 
  bind_rows(read_csv(file = "data/raw_data/CompToxChemicalsDashboard-Batch-Search_2021-04-25_11_10_17.csv", col_types = cols(), na = "-")) %>%
  filter(FOUND_BY != "NO_MATCH") %>% 
  filter(INPUT %in% missing) %>% 
  mutate(Norman_matched = paste0("comptox batch - ", FOUND_BY)) %>% 
  rename(Name = INPUT,
         Name_preferred = PREFERRED_NAME,
         Name_IUPAC = IUPAC_NAME,
         CAS_RN = CASRN,
         MS_Ready_SMILES = MS_READY_SMILES,
         Monoiso_Mass = MONOISOTOPIC_MASS,
         StdInChIKey = INCHIKEY,
         StdInChI = INCHI_STRING,
         Molecular_Formula = MOLECULAR_FORMULA) %>% 
  mutate(Monoiso_Mass = as.numeric(Monoiso_Mass))


# check for SusDat_ID
comtox_exists <- comtox %>%
  select(Name, Name_preferred, Norman_matched) %>%
  left_join(normannames %>% select(-CAS_RN_Dashboard), by = c("Name_preferred" = "Name")) %>% 
  filter(!is.na(Norman_SusDat_ID)) %>% 
  arrange(Namerank) %>% 
  group_by(Norman_SusDat_ID) %>% 
  slice_min(Namerank, n = 1) %>%
  ungroup() %>% 
  group_by(Name) %>% 
  slice_min(Namerank, n = 1) %>%
  ungroup() %>% 
  mutate(Norman_matched = Name_Source,
         Norman_matched_Synonym_Reliability = Reliability_of_Synonyms_ChemSpider) %>%
  distinct(Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability) %>% 
  left_join(effluent_substances %>% select(Name, Name_original), by = "Name")

# add rows with comptox matches that are already in the susdat list
IDs_matchnames <- IDs_matchnames %>% 
  bind_rows(comtox_exists)

effluent_substances <- effluent_substances %>% 
  left_join(IDs_matchnames, by = c("Name", "Name_original")) %>% 
  left_join(normansusdat, by = "Norman_SusDat_ID") %>% 
  mutate(Name = if_else(is.na(Norman_SusDat_ID), Name.x, Name.y))

# filter comptox search that is not in SusDat
comtox_add <- comtox %>% 
  filter(!(Name %in% comtox_exists$Name)) %>% 
  select(Name, Norman_matched, CAS_RN, Molecular_Formula, Monoiso_Mass, DTXSID, Name_preferred) %>% 
  left_join(effluent_substances %>% select(Name_original, Name), by = "Name") %>% 
  select(-Name) %>% 
  rename(Name = Name_preferred)
# remove them and then add rows in feature ID
effluent_substances <- effluent_substances %>% 
  filter(!(Name_original %in% comtox_add$Name_original)) %>% 
  bind_rows(comtox_add)  %>% 
  distinct_all()

effluent_substances <- effluent_substances %>% 
  mutate(CAS = str_remove_all(CAS_RN, "[a-zA-Z _:]"),
         CAS = if_else(is.na(CAS), CAS_RN_Dashboard, CAS)) %>% 
  select(Name_original, Name, CAS, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability, 
         CAS_RN, DTXSID, Molecular_Formula) %>% 
  distinct_all() %>% 
  bind_rows(data.frame(Name_original = "unknown",
                       Name = "unknown"))


# save the unique effluent substances detected
write_csv(effluent_substances %>% arrange(Name) %>% distinct_all(), "data/second_report_effluent_unique_substances.csv")


do_match <- FALSE

if(do_match == TRUE){
  
  ### partial match to the list, but results have to be sifted manually....
  Names_left <- effluent_substances %>% 
    filter(is.na(Norman_SusDat_ID)) %>% 
    distinct(Name)
  
  #match function
  matchnorman <- function(x, df) {
    tryCatch({
      df %>% 
        filter(str_detect(Name, x))
    }, error = function(e){
      data.frame(Name = "error")
    })
  }
  
  df <- normannames
  x <- Names_left$Name
  
  possiblematches_to_SusDat <- slice(df, 0) %>% mutate(lookup = NA_character_)
  
  for(i in x){
    
    df1 <-matchnorman(i, df)
    
    if(nrow(df1) > 0) {
      df1$lookup <- i
      possiblematches_to_SusDat <- possiblematches_to_SusDat %>% 
        bind_rows(df1)
    }
  }
  
  write_csv(possiblematches_to_SusDat, "data/possible_matches_to_Norman_SusDat.csv")
  
}


# save the feature IDs detected in effluent
featureIDs <- featureIDs %>% 
  select(-Name) %>% 
  left_join(effluent_substances %>% select(Name_original, Name), by = "Name_original") %>% 
  select(-Name_original)
  
write_csv(featureIDs %>% 
            arrange(Name) %>% 
            distinct_all(), "data/second_report_effluent_features.csv")


# features detected in water samples, transformed to longform.
features <- read_excel("data/raw_data/originals/Aligned_Peak_List_WithID_Extended_2709018_BA.xlsx") %>% 
  select(- ...1) %>% 
  rename(Denmark = `VANN-DANMARK01_Peak_List`,
         `Faroe Islands` = `VANN-FAROYENE01_Peak_List`,
         Norway = `VANN-NORGE01_Peak_List`,
         Sweden = `VANN-SVERIGE01_Peak_List`,
         Greenland = `VANN-GREENLAND01_Peak_List`,
         Iceland = `VANN-ISLAND01_Peak_List`,
         Finland = `VANN-FINLAND01_Peak_List`) %>% 
  pivot_longer(cols = c(Denmark, `Faroe Islands`, Norway, Sweden, Greenland, Iceland, Finland),
               names_to = "country_name", values_to = "QuantArea") %>% 
  filter(QuantArea > 0) %>% 
  mutate(matrix = "Effluent",
         Method = "decon+ULSA") %>% 
  left_join(sample_info %>% select(country_name, matrix, sample_no), by = c("country_name", "matrix")) %>% 
  select(-country_name, matrix)

write_csv(features %>% arrange(ID, sample_no), "data/second_report_effluent_results.csv")


# clean alldata (initial_report), the data from the first report
allData <- read_excel("data/raw_data/originals/AllData.xlsx", na = "NA") %>% 
  rename(Name_original = CommonName) %>% 
  rename(matrix = Sample,
         species = Species) %>% 
  mutate(matrix = if_else(matrix == "Water", "Effluent", matrix),
         species = if_else(species == "Water", "Effluent", species),
         species = str_to_title(species)) %>% 
  rename(country_name = Country) %>% 
  # get sample_no
  left_join(sample_info, by = c("country_name", "matrix", "species")) %>%
  rename(Monoiso_Mass = Weight,
         DrugType = `Drug type`,
         PotentialSource = `Potential Sources`) %>% 
  select(sample_no, SampleID, Method, Rt1, Rt2, QuantSN, QuantMass, QuantArea, Concentration,
         MolecularFormula, Monoiso_Mass, Smiles, CAS, Name_original, IDLevel,
         CompoundInfo, DrugType, PotentialSource)

# initial report unique substances
allData_substances <- allData %>% 
  distinct(Name_original, CAS)


# match to all possible names in SusDat
allData_namematch <- allData_substances %>% 
  mutate(Name = Name_original) %>% 
  filter(!is.na(Name_original)) %>% 
  left_join(normannames %>% filter, by = "Name") %>%
  filter(!is.na(Norman_SusDat_ID)) %>% 
  mutate(syn_rel = if_else(is.na(Reliability_of_Synonyms_ChemSpider), 10, Reliability_of_Synonyms_ChemSpider)) %>% 
  arrange(Namerank, syn_rel) %>% 
  group_by(Name) %>% 
  slice_head(n = 1) %>% 
  mutate(Norman_matched = Name_Source,
         Norman_matched_Synonym_Reliability = Reliability_of_Synonyms_ChemSpider) %>% 
  distinct(Name_original, Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability)

# match to all possible names in SusDat
allData_CASmatch <- allData_substances %>% 
  mutate(Name = Name_original) %>% 
  filter(!is.na(CAS)) %>% 
  filter(!(Name %in% allData_namematch$Name)) %>% 
  left_join(normanCASs %>% filter, by = "CAS") %>%
  filter(!is.na(Norman_SusDat_ID)) %>% 
  mutate(syn_rel = if_else(is.na(Reliability_of_CAS_ChemSpider ), 10, Reliability_of_CAS_ChemSpider )) %>% 
  arrange(CASrank, syn_rel) %>% 
  group_by(Name) %>% 
  slice_head(n = 1) %>% 
  rename(CAS_original = CAS) %>% 
  mutate(Norman_matched = CAS_Source,
         Norman_matched_Synonym_Reliability = Reliability_of_CAS_ChemSpider ) %>% 
  distinct(Name_original, Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability)

allData_match <- allData_namematch %>% 
  bind_rows(allData_CASmatch)

rm(allData_namematch, allData_CASmatch)

# missing info
initial_missing <- allData_substances %>% 
  filter(!(Name_original %in% allData_match$Name_original))
#write_csv(initial_missing %>% select(Name_original), "data/raw_data/initial_report_missing_info.csv")

initial_missing_cas <- allData_substances %>% 
  filter(!(Name_original %in% allData_match$Name_original)) %>% 
  filter(!is.na(CAS))
#write_csv(initial_missing_cas %>% select(CAS), "data/raw_data/initial_report_missing_info_cas.csv")


# batch search on CompTox -  https://comptox.epa.gov/dashboard
comtox <- read_csv(file = "data/raw_data/CompToxChemicalsDashboard-Batch-Search_2021-04-30_20_57_32.csv", col_types = cols(), na = "-") %>% 
  distinct_all() %>% 
  filter(FOUND_BY != "NO_MATCH") %>% 
  mutate(Norman_matched = paste0("comptox batch - ", FOUND_BY)) %>% 
  rename(Name = INPUT,
         Name_preferred = PREFERRED_NAME,
         Name_IUPAC = IUPAC_NAME,
         CAS_RN_Dashboard = CASRN,
         MS_Ready_SMILES = MS_READY_SMILES,
         Monoiso_Mass = MONOISOTOPIC_MASS,
         StdInChIKey = INCHIKEY,
         StdInChI = INCHI_STRING,
         Molecular_Formula = MOLECULAR_FORMULA) %>% 
  mutate(Monoiso_Mass = as.numeric(Monoiso_Mass))

comtox1 <- read_csv(file = "data/raw_data/CompToxChemicalsDashboard-Batch-Search_2021-04-30_21_12_49.csv", col_types = cols(), na = "-") %>% 
  distinct_all() %>% 
  filter(FOUND_BY != "NO_MATCH") %>% 
  filter(FOUND_BY != "Checksum Failed") %>% 
  mutate(Norman_matched = paste0("comptox batch - ", FOUND_BY)) %>% 
  rename(CAS_original = INPUT,
         Name = PREFERRED_NAME,
         Name_IUPAC = IUPAC_NAME,
         CAS_RN_Dashboard = CASRN,
         MS_Ready_SMILES = MS_READY_SMILES,
         Monoiso_Mass = MONOISOTOPIC_MASS,
         StdInChIKey = INCHIKEY,
         StdInChI = INCHI_STRING,
         Molecular_Formula = MOLECULAR_FORMULA) %>% 
  mutate(Monoiso_Mass = as.numeric(Monoiso_Mass)) %>% 
  filter(!(Name %in% comtox$Name_preferred))

# check for SusDat_ID
comtox_exists <- comtox %>%
  select(Name, Name_preferred, Norman_matched) %>%
  left_join(normannames %>% select(-CAS_RN_Dashboard), by = c("Name_preferred" = "Name")) %>% 
  filter(!is.na(Norman_SusDat_ID)) %>% 
  arrange(Namerank) %>% 
  group_by(Norman_SusDat_ID) %>% 
  slice_min(Namerank, n = 1) %>%
  ungroup() %>% 
  group_by(Name) %>% 
  slice_min(Namerank, n = 1) %>%
  ungroup() %>% 
  mutate(Norman_matched = Name_Source,
         Norman_matched_Synonym_Reliability = Reliability_of_Synonyms_ChemSpider) %>%
  distinct(Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability) %>% 
  left_join(effluent_substances %>% select(Name, Name_original), by = "Name")

# add rows with comptox matches that are already in the susdat list
allData_match <- allData_match %>% 
  bind_rows(comtox_exists)

allData_substances <- allData_substances %>% 
  mutate(Name = Name_original) %>% 
  mutate(CAS_original = CAS) %>% 
  left_join(allData_match, by = c("Name", "Name_original")) %>% 
  left_join(normansusdat, by = "Norman_SusDat_ID") %>%
  mutate(Name = if_else(is.na(Norman_SusDat_ID), Name.x, Name.y)) %>% 
  distinct(CAS_original, Name_original, Name, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability,
           CAS_RN, CAS_RN_Dashboard, Molecular_Formula, Monoiso_Mass, DTXSID)

# filter comptox search that is not in SusDat
comtox_add <- comtox %>% 
  filter(!(Name %in% comtox_exists$Name)) %>% 
  select(Name, Norman_matched, CAS_RN_Dashboard, Molecular_Formula, Monoiso_Mass, DTXSID, Name_preferred) %>% 
  left_join(allData_substances %>% select(CAS_original, Name_original, Name), by = "Name") %>% 
  select(-Name) %>% 
  rename(Name = Name_preferred)
# remove them and then add rows in feature ID
allData_substances <- allData_substances %>% 
  filter(!(Name_original %in% comtox_add$Name_original)) %>% 
  bind_rows(comtox_add)  %>% 
  distinct_all()

# filter comptox search that is not in SusDat
comtox_add <- comtox1 %>% 
  filter(!(Name %in% comtox_exists$Name)) %>% 
  select(CAS_original, Name, Norman_matched, CAS_RN_Dashboard, Molecular_Formula, Monoiso_Mass, DTXSID) %>% 
  left_join(allData_substances %>% select(CAS_original, Name_original), by = "CAS_original")
# remove them and then add rows in feature ID
allData_substances <- allData_substances %>% 
  filter(!(Name_original %in% comtox_add$Name_original)) %>% 
  bind_rows(comtox_add)  %>% 
  distinct_all()

initial_report_substances <- allData_substances %>% 
  left_join(allData %>% rename(CAS_original = CAS) %>% select(-Monoiso_Mass), by = c("Name_original", "CAS_original")) %>% 
  mutate(CAS = str_remove_all(CAS_RN, "[a-zA-Z _:]"),
         CAS = if_else(is.na(CAS), CAS_original, CAS)) %>% 
  select(CAS_original, Name_original, Name, CAS, Norman_SusDat_ID, Norman_matched, Norman_matched_Synonym_Reliability, 
         CAS_RN, DTXSID, Molecular_Formula, CompoundInfo, DrugType, PotentialSource) %>% 
  distinct_all()

write_csv(initial_report_substances %>% arrange(Name), "data/first_report_unique_substances.csv")

initial_report_results <- allData %>% 
  select(-CAS, -CompoundInfo, -DrugType, -PotentialSource) %>% 
  left_join(initial_report_substances %>% select(Name_original, Name, CAS), by = "Name_original") %>% 
  select(-Name_original) %>% 
  distinct_all() %>% 
  arrange(Name)

write_csv(initial_report_results %>% arrange(Name), "data/first_report_results.csv")


# clean the re-quantified GCGC data
gcgc <- read_excel("data/raw_data/originals/ReQuantifiedBiotaGCGC.xlsx") %>% 
  rename(Name_original = CommonName) %>% 
  rename(country_name = Country,
         matrix = Sample) %>% 
  select(-CAS) %>% 
  left_join(initial_report_substances %>% select(Name_original, Name, -CAS), by = c("Name_original")) %>% 
  select(-Name_original) %>% 
  left_join(sample_info %>% select(country_name, matrix, sample_no), by = c("country_name", "matrix")) %>% 
  select(-country_name, -matrix)

write_csv(gcgc %>% arrange(Name), "data/second_report_requantified_biota_gcgc.csv")


