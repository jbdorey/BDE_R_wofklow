# This function was started on the 20th of January 2025 by James Dorey in order to make some 
# country countryColumns more consistent for his paper titled "How many bee species are there? A 
# quantitative global estimate".
  # For questions, email him on jbdorey@me.com

countryHarmoniseR <- function(
    data = NULL,
    countryColumn = NULL,
      # set to TRUE in order to match small entities to continents
    continentAnalysis = FALSE){
  
  
  
  data <- data %>%
    dplyr::rename(countryColumn = tidyselect::any_of(countryColumn)) %>%
    # Shorten some country countryColumns
    dplyr::mutate(countryColumn = dplyr::if_else(countryColumn %in% c("United States of America","United States"),
                                        "USA",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Russia", "Russian Federation", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Democratic Republic of the Congo","DRC",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Lao People's Democratic Republic", "Lao PDR"),
                                        "Lao",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "United Arab Emirates","UAE",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Central African Republic","CAR",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Saint Vincent and the Grenadines","SV & Grenadines",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Republic of Cabo Verde","Cape Verde"),
                                        "Cabo Verde",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Northern Mariana Islands","N. Mariana",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("São Tomé and Principe","Sao Tome and Principe",
                                                                      "São Tomé and Príncipe"),
                                        "São Tomé",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Republic of the Congo","Republic of Congo"),
                                        "Rep. Congo",countryColumn),
                  countryColumn = dplyr::if_else(stringr::str_detect(countryColumn, "Northern") & stringr::str_detect(countryColumn, "Cyprus"),
                                        "Northern Cyprus", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("N. Cyprus"),
                                        "Northern Cyprus",countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "The Gambia", "Gambia", countryColumn),
                  countryColumn = dplyr::if_else(stringr::str_detect(countryColumn, "North") & stringr::str_detect(countryColumn, "Macedonia"),
                                        "Macedonia", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Kingdom of eSwatini","Swaziland", "eSwatini"), 
                                        "Eswatini", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Barbuda", "Antigua"),
                                        "Antigua and Barbuda", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Brussels","Flemish Region","Walloon Region"),
                                        "Belgium", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Caribbean Netherlands","Martinique"),
                                                      # RecountryColumnd Curaçao to match to continent
                                                      "Curaçao", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Cocos Islands","Indian Ocean Territories", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Federation of Bosnia and Herzegovina",
                                                    "Republic Srpska","Vojvodina", "Bosnia"),
                                        "Bosnia and Herzegovina", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Gibraltar","Spain", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Guadeloupe","Barbados", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Mayotte","Madagascar", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Réunion", "Reunion"),"Mauritius", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Svalbard Islands","Bouvet"),
                                        "Norway", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("West Bank","Gaza"),"Palestine", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "French Southern Territories",
                                        "French Southern and Antarctic Lands", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "South Georgia and South Sandwich Islands","South Georgia and the Islands", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Aland Islands","Åland Islands", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Falkland Islands","Falkland Islands / Malvinas", countryColumn),
                  countryColumn = dplyr::if_else(stringr::str_detect(countryColumn, "Darussalam"),"Brunei Darussalam", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn %in% c("Scotland","Wales","England"),
                                        "United Kingdom", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Tokelau","New Zealand", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Czech Republic", "Czechia", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "South Korea", "Republic of Korea", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Wallis and Futuna Islands", "Wallis and Futuna", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Faeroe Islands","Faroe Islands", countryColumn),
                  countryColumn = dplyr::if_else(countryColumn == "Saint Martin","Saint-Martin", countryColumn),
                  )
  
  if(continentAnalysis == TRUE){
    data = data %>% 
      dplyr::mutate(
        countryColumn = dplyr::if_else(countryColumn == "French Guiana","Brazil", countryColumn),
      )
  }
  
    # Re-apply the original column name
  names(data)[names(data) == "countryColumn"] <- countryColumn
  
  # Return the data
  return(data)
  
}



