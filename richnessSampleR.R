# This function was started by James B Dorey on the 2nd of November 2024 in order to sample 
# iNEXT or iChao with input data varying as drawn from a curve derived from sampling the 
# taxonomic literature
# For help, contact James at jbodrey@me.com


richnessSampleR <- function(
    # From X:Y iterations
  iterations = NULL,
  richnessInputs = richnessInputs,
  # What scale?
  scale = NULL, # Country, Continent, or Global
  # Functions
  ChaoWrapper = ChaoWrapper,
  iNEXTwrapper = iNEXTwrapper
){
  
  #### 0.0 Prep ####
  ##### 0.1 richnessInputs ####
  # Extract the datasets from richnessInputs
  literatureCurve = richnessInputs$litTibble
  taxonomy_valid = richnessInputs$taxonomy_valid
  noPointSpecies_diff = richnessInputs$noPointSpecies_diff
  taxonomyFile = richnessInputs$taxonomyFile
  beeData_counts = richnessInputs$beeData_counts
  continentChecklist = richnessInputs$continentChecklist
  worldMap = richnessInputs$worldMap
  countryNOTcontinent = richnessInputs$countryNOTcontinent
  country_speciesCounts = richnessInputs$country_speciesCounts
  
    # Fix some names in the country_speciesCounts
  country_speciesCounts <- country_speciesCounts %>% 
    dplyr::mutate(country_suggested = dplyr::if_else(country_suggested == "Russian Federation",
                                                      "Russia", country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Bosnia and Herzegovina",
                                                      "Bosnia", country_suggested),
                  country_suggested = dplyr::if_else(stringr::str_detect(country_suggested, "eSwatini"),
                                                     "Kingdom of eSwatini", country_suggested),
                  country_suggested = dplyr::if_else(stringr::str_detect(country_suggested, "Eswatini"),
                                                     "Kingdom of eSwatini", country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Russia",
                                                     "Russian Federation", country_suggested)
                  ) %>%
    # Shorten some country names
    dplyr::mutate(country_suggested = dplyr::if_else(country_suggested == "United States of America","USA",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Russian Federation","Russia",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Democratic Republic of the Congo","DRC",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Lao People's Democratic Republic","Lao",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "United Arab Emirates","UAE",country_suggested),
                  country_suggested = dplyr::if_else(stringr::str_detect(country_suggested,"Central African Republic"), "CAR",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Saint Vincent and the Grenadines","SV & Grenadines",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Republic of Cabo Verde","Cabo Verde",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Northern Mariana Islands","N. Mariana",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "São Tomé and Principe","São Tomé",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Bosnia and Herzegovina","Bosnia",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Kingdom of eSwatini","eSwatini",country_suggested),
                  country_suggested = dplyr::if_else(country_suggested == "Republic of the Congo","Rep. Congo",country_suggested))
    # combine some country names
  checklistFile = richnessInputs$checklistFile %>% 
    # Change country names
    dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "United States",
                                                      "United States of America", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "USA",
                                                      "United States of America", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Aland Islands",
                                                      "Åland Islands", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels",
                                                      "Belgium", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "DRC",
                                                      "Democratic Republic of the Congo", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Swaziland",
                                                      "Kingdom of eSwatini", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "The Gambia",
                                                      "Gambia", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Republic of Congo",
                                                      "Republic of the Congo", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Rep. Congo",
                                                      "Republic of the Congo", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cape Verde",
                                                      "Republic of Cabo Verde", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "West Bank",
                                                      "Palestine", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gaza",
                                                      "Palestine", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Kosovo",
                                                      "Serbia", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Federation of Bosnia and Herzegovina",
                                                      "Bosnia and Herzegovina", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Barbuda",
                                                      "Antigua and Barbuda", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name,
                                                                          "Brunei"),
                                                      "Brunei Darussalam", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Lao PDR",
                                                      "Lao People's Democratic Republic", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Lao",
                                                      "Lao People's Democratic Republic", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Saint-Martin",
                                                      "Saint Martin", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Somaliland",
                                                      "Somalia", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Wallis and Futuna Islands",
                                                      "Wallis and Futuna", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "CAR",
                                                      "Central African Republic", rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Faeroe Islands",
                                                      "Faroe Islands", rNaturalEarth_name)) %>%
    
    # Shorten some country names
    dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "United States of America","USA",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Russian Federation","Russia",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Democratic Republic of the Congo","DRC",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Lao People's Democratic Republic","Lao",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "United Arab Emirates","UAE",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Central African Republic","CAR",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Saint Vincent and the Grenadines","SV & Grenadines",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Republic of Cabo Verde","Cabo Verde",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Northern Mariana Islands","N. Mariana",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "São Tomé and Principe","São Tomé",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Bosnia and Herzegovina","Bosnia",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Kingdom of eSwatini","eSwatini",rNaturalEarth_name),
                  rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Republic of the Congo","Rep. Congo",rNaturalEarth_name))
  
  ##### 0.2 Functions ####
  
  # Make a function to sample the country sample sizes
  countrySampler <- function(df = NULL){ 
    # Set a random seed
    set.seed(iterations*7)
    # Get the maximum value for literatureCurve$x_coords
    max_x <- max(literatureCurve$x_coords)
    # Test that there is sufficient length (i.e., at least one records with n > NA)
    lengthTest <- df$n[complete.cases(df$n)] %>% 
      length()
    # Get the maximum value for the country's sample size or set the minium to n = 1
    if(lengthTest > 0){
      max_n = max(df$n[complete.cases(df$n)])
    }else{
      max_n = 1
    }
    # Find the length of the indput data
    inputLength <- length(df$n) + 100
    # Take the smaller of these two values
    maxSample <- dplyr::if_else(max_x > max_n, max_n, max_x)
    # Now apply the maxxed distribution to the country's checklist data
    df <- df %>%
      dplyr::mutate(n = dplyr::if_else(is.na(n), 
                                       # Draw from the literature distribution
                                       sample(x = literatureCurve$x_coords[1:maxSample],
                                              prob = literatureCurve$y_coords[1:maxSample], 
                                              replace = TRUE,
                                              size = inputLength)[dplyr::row_number()] %>%
                                         round(0),
                                       # Or keep n
                                       n),
                    dataFrom = dplyr::if_else(is.na(dataFrom), "checklist", dataFrom))
    return(df)
  }# END countrySampler
  
  
  #### 1.0 Sample the curve ####
  # Sample the curve at the start of each iteration to generate new values
  # Get a random sample based on the empirical dataset
  litSample <- sample(x = literatureCurve$x_coords, prob = literatureCurve$y_coords, 
                      replace = TRUE,
                      size = 10000) %>%
    round(0) 
  
  
  ###### 1.1 No occ. points #####
  # Find the species that are not represented in the occurrence dataset
  noPointSpecies <- noPointSpecies_diff  %>%
    # Add a sample size that is drawn from the same distriubtion of the literature and no-occurrences
    # dataset and generated in litSample
    dplyr::tibble(scientificName = .,
                  n = litSample[1:length(.)],
                  dataFrom = "randomisedSample")
  
  # 
  beeData_totalCounts <- beeData_counts %>% 
    dplyr::bind_rows(noPointSpecies)
  
  #### 2.0 Data prep ####
  ##### 2.1 Country ####
  if(scale == "Country"){
    # Now, for those species-country combinations without ANY points, but that should have points
    # from the checklist, make n = a number drawn from the literature records distribution
    # and state that they came from the checklist
    country_speciesChecklistCounts <- checklistFile %>%
      dplyr::select(validName, rNaturalEarth_name, year, author, authorCount) %>%
      dplyr::full_join(country_speciesCounts, by = c("validName" = "scientificName",
                                                     "rNaturalEarth_name" = "country_suggested")) %>%
      # Group by country and then split into a list per group
      dplyr::group_by(rNaturalEarth_name) %>%
      dplyr::group_split() %>%
      # For each country, apply the literature distribution but with an n max of the empirical maximum
      # for that country
      lapply(., countrySampler) %>%
      dplyr::bind_rows() %>%
      # Fix up the year, author, and authorCount columns by removing and re-adding them
      dplyr::select(!c(year, author, authorCount)) %>% 
      dplyr::left_join(taxonomyFile %>%
                         dplyr::select(validName, year, author, authorCount),
                       by = "validName") %>%
      dplyr::distinct(validName, rNaturalEarth_name, .keep_all = TRUE)
  
    
    # Create the Chao input data including checklist data
    countryChaoData_checklist <- country_speciesChecklistCounts %>%
      dplyr::select(validName, rNaturalEarth_name, n) %>%
      tidyr::pivot_wider(names_from = rNaturalEarth_name,
                         values_from = n,
                         values_fill = 0) %>%
      # Create the rownames
      tibble::column_to_rownames("validName") %>%
      dplyr::tibble()
    
  } # END if country
  
  
  ##### 2.2 Continent ####
  if(scale == "Continent"){
    # Turn the country occ data into a continent one
    continentOccs <- country_speciesCounts %>%
      dplyr::ungroup() %>%
      # Change some country names to better match the continent data
      dplyr::mutate(country_suggested = dplyr::if_else(country_suggested == "eSwatini",
                                                       "Kingdom of eSwatini", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "Faroe Islands",
                                                       "Faeroe Islands", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "French Guiana",
                                                       # Renamed Brazil to match to South America
                                                       "Brazil", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "Lao People's Democratic Republic",
                                                       "Laos", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "Macedonia",
                                                       "North Macedonia", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "Saint Martin",
                                                       "Saint-Martin", country_suggested),
                    country_suggested = dplyr::if_else(country_suggested == "Wallis and Futuna",
                                                       "Wallis and Futuna Islands", country_suggested),
                    # Names to a nearby island
                    country_suggested = dplyr::if_else(country_suggested == "Martinique",
                                                       "Barbados", country_suggested)) %>%
      # Join first by name
      dplyr::left_join(worldMap %>% dplyr::select(name, continent) %>%
                         sf::st_drop_geometry(), by = c("country_suggested" = "name")) %>%
      # Then join by name_long
      dplyr::left_join(worldMap %>% dplyr::select(name_long, continent) %>%
                         sf::st_drop_geometry(), by = c("country_suggested" = "name_long")) %>%
      # Now merge these continent columns
      dplyr::mutate(continent = dplyr::if_else(is.na(continent.x),
                                               continent.y, continent.x)) %>%
      # drop interim and old count columns
      dplyr::select(!c("continent.x", "continent.y")) %>%
      dplyr::group_by(scientificName, continent) %>%
      dplyr::mutate(sum = sum(n)) %>%
      dplyr::mutate(name_continent = stringr::str_c(scientificName, continent, sep = "__"))
    
    
    # Combine the occurrence and checklist data 
    continentCounts <- continentOccs %>%
      dplyr::select(!c("country_suggested", "n", "name_continent")) %>%
      dplyr::rename(n = sum) %>%
      dplyr::bind_rows(countryNOTcontinent) %>%
      # Group by country and then split into a list per group
      dplyr::group_by(continent) %>%
      dplyr::group_split() %>%
      # For each country, apply the literature distribution but with an n max of the empirical maximum
      # for that country
      lapply(., countrySampler) %>%
      dplyr::bind_rows() %>%
      dplyr::distinct(scientificName, continent, .keep_all = TRUE) 
    
    # # Pivot wider 
     continentWider <- continentCounts %>%
       dplyr::select(scientificName, continent, n) %>%
       tidyr::drop_na() %>% 
       tidyr::pivot_wider(names_from = continent,
                          values_from = n,
                          values_fill = 0) %>%
       # Create the rownames
       tibble::column_to_rownames("scientificName") %>%
       dplyr::tibble()
    
  } # END if continent
  
  
  
  #### 3.0 Analysis ####
  ##### 3.1 Country ####
  if(tolower(scale) == "country"){
    ###### a. Chao ####
    ChaoOut <- ChaoWrapper(data = countryChaoData_checklist,
                           k = 5,
                           datatype = "abundance",
                           conf = 0.95,
                           mc.cores = 1)
    
  } # END if country
  
  if(tolower(scale) == "country"){
    ###### b. iNext ####
    # Run iNEXT for each country using the iNEXT wrapper
    iNEXTout <- iNEXTwrapper(data = country_speciesChecklistCounts,
                             variableColumn = "rNaturalEarth_name",
                             valueColumn = "n",
                             datatype = "abundance",
                             mc.cores = 1)
  } # END if continent
  
  
  ##### 3.2 Continent ####
  if(tolower(scale) == "continent"){
    ###### a. Chao ####
    ChaoOut <- ChaoWrapper(data = continentWider,
                           k = 5,
                           datatype = "abundance",
                           conf = 0.95,
                           mc.cores = 1)
    
  } # END if country
  
  if(tolower(scale) == "continent"){
    ###### b. iNext ####
    # Run iNEXT for each country using the iNEXT wrapper
    iNEXTout <- iNEXTwrapper(data = continentCounts,
                             k = 5,
                             variableColumn = "continent",
                             datatype = "abundance",
                             conf = 0.95,
                             mc.cores = 1)
  } # END if continent
  
  
  ##### 3.3 Global ####
  if(tolower(scale) == "global"){
    
    ###### a. iChao ####
    # Get the Chao values for global bee species diversity
    ChaoOut <- SpadeR::ChaoSpecies(beeData_totalCounts$n,
                                   datatype = "abundance", 
                                   k = 10, conf = 0.95)
    
    ##### b. iNEXT ####
    global_iNEXT <- beeData_totalCounts %>%
      dplyr::select(scientificName, n) %>%  # Create the rownames
      tibble::column_to_rownames("scientificName") %>%
      dplyr::tibble()
    
    iNEXTout <- iNEXT::iNEXT(x = global_iNEXT$n, datatype = "abundance",
                             # 0 for species richness Shannon diversity (q = 1, the exponential of Shannon entropy) and Simpson diversity (q = 2, the inverse of Simpson concentration).
                             q = 0)
  } # END if global
  
  # Return a list of the objects produced
  return(dplyr::lst(ChaoOut, iNEXTout) %>%
           stats::setNames(c("ChaoOut", "iNEXTout")) %>%
           dplyr::lst() %>%
           # Set the name as the iteration number 
           stats::setNames(iterations))
} # END function




