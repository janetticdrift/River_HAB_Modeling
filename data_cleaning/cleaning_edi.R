###########################
#Environmental Data Initative (EDI): Cleaned Files
###########################
#Read in cleaned and organized datasets from cleaning.R, and reorganize and rename 
#dataframes for export as CSV files to upload to EDI

#Files include water chemistry data collected by the Blaszczak Lab (University of 
#Nevada, Reno), benthic percent cover estimates, microscopy analyses, 
#anatoxin concentrations normalized by organic matter, and qPCR gene copy estimates.

#Read in the cleaned data files
source(here::here("data_cleaning/cleaning_HAB.R"))

###########################
#Water Chemistry
###########################
#In cleaning_HAB.R:
  #1) Values that were below the minimum detection value were replaced with 
      #half the minimum detection value
  #2) For 2024, ammonium was calculated using ammonia, pH, and temperature
  #3) Removed an outlier nitrate value in 2023 and replaced it with NA
  #4) Replaced an instrument error in conductivity reading with the lowest 
      #equivalent HOBO estimate
  
water_chemistry_data <- nutrients_raw %>% 
  dplyr::select(!c(nitrate_replace, ammonia_replace, year, ammonia)) %>% 
  dplyr::relocate(field_date, site_reach, site, reach, 
                  nitrate_mg_N_L, oPhos_ug_P_L, ammonium_mg_N_L, DIN_mg_N_L,
                  temp_C, cond_uS_cm, pH, DO_mg_L)


###########################
#Benthic Algae Percent Cover
###########################

percent_cover_data <- percoverdata %>% 
  dplyr::relocate(other_nfixers, .before = bare_biofilm) %>% 
  dplyr::select(!year) %>% 
  dplyr::arrange(field_date)

###########################
#Mat Microscopy Analysis
###########################

microscopy_data <- microscopy1 %>% 
  ungroup() %>% 
  dplyr::select(!c(year, slide_rep)) %>% 
  unite("site", site:location, sep = "-", remove = F) %>% 
  unite("site_reach", c(site, reach), sep = "-", remove = F) %>% 
  dplyr::select(!location) %>% 
  dplyr::arrange(desc(sample_type), field_date) %>% 
  pivot_wider(names_from = Species, values_from = Abundance) %>% 
  dplyr::relocate(date_analyzed, .after = field_date) %>% 
  dplyr::relocate(leptolyngbya_and_geitlerinema, .after = calothrix) %>% 
  dplyr::rename(analyzed_date = date_analyzed, mat_sampled = sample_type) %>% 
  dplyr::mutate(mat_sampled = case_when(mat_sampled == "TM" ~ "Microcoleus",
                                        mat_sampled == "TAC" ~ "Anabaena"))



###########################
#Anatoxin Concentrations
###########################

anatoxin_data <- toxindf %>% 
  dplyr::mutate(site = "SFE-M") %>% 
  dplyr::select(!year) %>% 
  unite("site_reach", c(site, reach), sep = "-", remove = F) %>% 
  dplyr::relocate(site, .after = site_reach) %>% 
  dplyr::rename(mat_sampled = sample_type) %>% 
  dplyr::mutate(per_org_matter = ifelse(per_org_matter>1, 1, per_org_matter))

###########################
#qPCR Gene Copies
###########################

genecopy <- read.csv(here::here("data/qPCR_2024.csv"))
