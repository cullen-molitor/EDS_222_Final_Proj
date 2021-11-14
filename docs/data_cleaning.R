#
#
#
#
#     Use for Manipulating data frames 
#
#
#     Use ALT + O to see outline
#
#

{ # Library   ----
 
  library(tidyverse)
  library(here)
  library(lubridate)
  library(vegan)
  library(rnoaa)
  
}

{ # Metadata    ----
  
  Species_Info <- readr::read_csv(
    here::here(
    "data", "meta", "Species_Complete.csv"
    ))
  
  Site_Info <- readr::read_csv(
    here::here(
      "data", "meta", "Site_Info.csv"
    ))
 
  VFT_Species <- c(
    "blacksmith", 
    "black_surfperch", 
    "striped_surfperch", 
    "opaleye",                     
    "garibaldi", 
    "senorita",
    "kelp_bass", 
    "pile_perch",                 
    "kelp_rockfish", 
    "blue_rockfish", 
    "olive_rockfish", 
    "California_sheephead_female",
    "California_sheephead_male", 
    "rock_wrasse_female", 
    "rock_wrasse_male")

}

{ # SST Anomaly Index (ONI)   ----
  
    oni <- read.table( # Read in  ONI to be added to all data
      "https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/detrend.nino34.ascii.txt",
      header = T) %>%
      dplyr::mutate(Date = as.Date(ISOdate(YR, MON, 1)),
                    DateStart = as.Date(ISOdate(YR, MON, 1)),
                    DateEnd = ceiling_date(DateStart, "month")) %>%
      dplyr::rename(ONI_ANOM = ANOM,
                    Month = MON,
                    SurveyYear = YR) %>% 
      dplyr::select(SurveyYear, Month, Date, DateStart, DateEnd, ONI_ANOM) 
 
}

{ # Secchi  ----
  
  Secchi <- 
    readr::read_csv(
      file = here::here(
        "data", "raw", "KFM_Secchi_1985-2019.txt")) %>%
    dplyr::rename(SiteCode = Site_Code) %>% 
    dplyr::left_join(Site_Info) %>%
    tidyr::separate(Date, into = c('Date','Time'), sep = ' ') %>%
    dplyr::mutate(Date = lubridate::mdy(Date),
                  SurveyYear = lubridate::year(Date)) %>% 
    dplyr::select(SiteNumber, IslandCode, IslandName, SiteCode, SiteName,
                  Date, SurveyYear, HorizontalSecci, Surge, MeanDepth) 
  
}

{ # RDFC    ----
  RDFC_Density <- 
    readr::read_csv(
      file = here::here("data", "raw", "KFM_RovingDiverFishCount_RawData_1982-2019.txt"),
      col_names = TRUE, locale = locale(encoding = "ISO-8859-1"),
      col_types = cols(Count = col_double())) %>%
    tidyr::separate(SurveyDate, c('Date','Time'),' ') %>%
    dplyr::mutate(Date = lubridate::mdy(Date),
                  Count = as.double(Count)) %>% 
    dplyr::left_join(Site_Info) %>%
    dplyr::group_by(SiteCode, SurveyYear) %>% 
    dplyr::filter(IslandCode != "CL",
                  Date == base::max(Date), SurveyYear > 2003, 
                  ExperienceLevel == "E") %>%
    dplyr::ungroup() %>% 
    dplyr::mutate(CommonName = gsub('ñ', 'n', CommonName),
                  Abundance = gsub("c", "C", Abundance),
                  Abundance = gsub("f", "F", Abundance),
                  Abundance = gsub("s", "S", Abundance),
                  Abundance = gsub("m", "M", Abundance),
                  Abundance = gsub("^$", NA, Abundance),
                  Abundance = gsub("-", NA, Abundance),
                  Score = ifelse(is.na(Score), 0, Score),
                  Abundance = ifelse(is.na(Abundance), "N", Abundance),
                  Count = ifelse(Score == 0, 0, Count)) %>% 
    dplyr::filter(!CommonName %in% c(
      "black surfperch, all", "blacksmith, all", "blue rockfish, all", "kelp bass, all", 
      "kelp rockfish, all", "olive rockfish, all",  "opaleye, all", "pile perch, all", 
      "senorita, all", "striped surfperch, all")) %>% 
    dplyr::group_by(SiteNumber, IslandCode, IslandName, SiteCode, SiteName, SurveyYear, Date, 
                    ScientificName, CommonName, ReserveStatus, Reference) %>%
    dplyr::summarise(Count = base::mean(Count, na.rm = TRUE)) %>% 
    dplyr::mutate(Count = base::ifelse(Count > 0 & Count < 1, 1, base::round(Count, 0))) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(Secchi) %>% 
    readr::write_csv(here::here("data", "tidy", "rdfc.csv"))
  
}

{ # VFT   ----
  VFT_Density <- 
    readr::read_csv(
      here::here("data", "raw", "KFM_VisualFishTransect_RawData_1985-2019.txt"),
      locale = locale(encoding = "ISO-8859-1"), col_names = TRUE, col_types = NULL) %>% 
    tidyr::separate(SurveyDate, c('Date','Time'),' ') %>%
    dplyr::mutate(Date = lubridate::mdy(Date),
                  CommonName = gsub('ñ', 'n', CommonName)) %>%
    dplyr::left_join(Site_Info) %>%
    dplyr::rename(Count = CountA) %>%
    dplyr::group_by(SiteNumber, SurveyYear) %>%
    dplyr::filter(Date == max(Date), 
                  IslandCode != "CL") %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(SiteNumber, SurveyYear, Species) %>% 
    tidyr::pivot_wider(names_from = "Transect_Number", names_prefix = "T", values_from = "Count") %>%
    dplyr::mutate(T1 = ifelse(SurveyYear < 1997, T1, T1 + T2),
                  T2 = ifelse(SurveyYear < 1997, T2, T3 + T4)) %>% 
    dplyr::select(-T3, -T4) %>%
    tidyr::pivot_longer(cols = c(T1, T2), values_to = "Count", names_to = "Transect_Number") %>%
    dplyr::ungroup() %>% 
    dplyr::group_by(SiteNumber, IslandCode, IslandName, SiteCode, SiteName, SurveyYear, Date,
                    CommonName, ScientificName, ReserveStatus, Reference, ReserveYear) %>%
    dplyr::summarise(Mean_Density = round(sum(Count, na.rm = TRUE) / 600, 4)) %>% 
    dplyr::ungroup() %>%  
    dplyr::distinct(SiteNumber, IslandCode, IslandName, SiteCode, SiteName, SurveyYear, Date,
                    CommonName, ScientificName, Mean_Density, ReserveStatus, Reference)  %>% 
    dplyr::left_join(Secchi) %>% 
    readr::write_csv(here::here("data", "tidy", "vft.csv"))
  
}

{ # Diversity   ----

  VFT_Diversity <- VFT_Density %>% 
    dplyr::group_by(SiteNumber, SurveyYear) %>% 
    dplyr::summarise(shannon = vegan::diversity(Mean_Density)) %>% 
    dplyr::ungroup() %>% 
    dplyr::left_join(Secchi) %>% 
    dplyr::left_join(Site_Info) %>%
    dplyr::select(SiteNumber, IslandCode, IslandName, SiteCode, SiteName, SurveyYear, Date,
                  shannon, ReserveStatus, Reference, HorizontalSecci, Surge, MeanDepth) %>% 
    dplyr::mutate(survey = "VFT") 

  RDFC_Diversity <- RDFC_Density %>% 
    dplyr::group_by(SiteNumber, SurveyYear) %>% 
    dplyr::summarise(shannon = vegan::diversity(Count)) %>%  
    dplyr::ungroup() %>% 
    dplyr::left_join(Secchi) %>%  
    dplyr::left_join(Site_Info) %>%
    dplyr::select(SiteNumber, IslandCode, IslandName, SiteCode, SiteName, SurveyYear, Date,
                  shannon, ReserveStatus, Reference, HorizontalSecci, Surge, MeanDepth) %>% 
    dplyr::mutate(survey = "RDFC")
  
  fish_diversity <- rbind(VFT_Diversity, RDFC_Diversity) %>% 
    readr::write_csv(here::here("data", "tidy", "fish_diversity.csv"))

}

{ # Temperature   ----
  
  temp <- readr::read_csv( 
    here::here("data", "raw","Temperature_RawData_1994-2019.txt")) %>%
    dplyr::filter(IslandCode != "CL", Site_Number < 38, !base::is.na(Temp_C)) %>%
    dplyr::select(-Date, -Time) %>%
    tidyr::separate(DateTime, into = c('Date','Time'), sep = ' ') %>%
    dplyr::mutate(Date = lubridate::mdy(Date)) %>% 
    dplyr::rename(SiteNumber = Site_Number) %>% 
    dplyr::group_by(Date, SiteNumber) %>%
    dplyr::summarise(
      Daily_Min = base::min(Temp_C),
      Daily_Max = base::max(Temp_C),
      Daily_Mean = base::mean(Temp_C)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(Site_Info) %>%
    dplyr::select(SiteNumber, IslandCode, IslandName, SiteCode, SiteName,
                  Date, Daily_Mean, Daily_Min, Daily_Max, MeanDepth) %>%
    readr::write_csv(here::here("data", "tidy","kfm_temp.csv"))
  
}

{ # Buoy  ----
  
  buoy_sm <- rnoaa::buoy(
    dataset = 'stdmet',
    buoyid = '46011', 
    year = 1980)[["data"]] %>% 
    tidyr::separate(time, into = c('date','time'), sep = 'T') %>%
    dplyr::mutate(date = lubridate::ymd(date)) %>% 
    dplyr::group_by(date) %>% 
    dplyr:: summarise(sst_mean = base::mean(sea_surface_temperature, na.rm = T))
  
  for (yr in c(1981:1992, 1994:2020)) {
    temp <- rnoaa::buoy(
      dataset = 'stdmet',
      buoyid = '46011', 
      year = yr)[["data"]] %>% 
      tidyr::separate(time, into = c('date','time'), sep = 'T') %>%
      dplyr::mutate(date = lubridate::ymd(date)) %>% 
      dplyr::group_by(date) %>% 
      dplyr::summarise(sst_mean = base::mean(sea_surface_temperature, na.rm = T))
    buoy_sm <- rbind(buoy_sm, temp)
  } 
  
  readr::write_csv(
    buoy_sm, 
    file = here::here("data", "tidy", "buoy_sst.csv"))
  
}









