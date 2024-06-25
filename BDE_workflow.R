# This R script was written by James Dorey, starting on the 7th of March 2024. 
# Here I will use DiscoverLife + BeeBDC data to begin to estimate species diversity of bees.
# For queries, please feel free to contact James Dorey at jbdorey@me.com


#### 0.0 Script preparation ####
##### 0.1 Working directory ####
# Choose the path to the root folder in which all other folders can be found (or made by dirMaker)
RootPath <- "/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow"

# Set the working directory
setwd(RootPath)

# Initialise renv the project if needed
# renv::init(project = RootPath) 
renv::activate(project = RootPath)


# Install BeeBDC from CRAN
install.packages("BeeBDC")
install.packages("bdc")
install.packages("iNEXT")
install.packages("readxl")
install.packages("SpadeR")
install.packages("mosaic")
install.packages("mosaicCalc")
  # install phyloseq before installing breakaway
BiocManager::install("phyloseq")
  # Generate PAT if needed
usethis::create_github_token()
devtools::install_github("adw96/breakaway", 
                         auth_token = "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")
# You could also install BeeBDC's development version using the below: 
# WARNING the development version may not pass all CRAN or GitHub tests.
# remotes::install_github("https://github.com/jbdorey/BeeBDC.git", user="jbdorey", 
#                         # To use the development version, do below, otherwise choose "main"
#                         ref = "devel", 
#                         force = TRUE)


##### 0.2 Load packages ####
# Save a snapshot of the environment
renv::activate(project = RootPath)
renv::snapshot(project = RootPath)
# Load all packages from the list specified above,
lapply(c("BeeBDC", "magrittr", "bdc", "SpadeR", "ggplot2", "dplyr", "phyloseq", "breakaway"), 
       library, character.only = TRUE)

# Load in an extra wrapper around SpadeR
source("ChaoWrapper.R")
source("ChaoWrapper.R")
source("iNEXTwrapper.R")



#### 1.0 Data prep ####
##### 1.1 Download data ####
# Download the taxonomy
taxonomyFile <- BeeBDC::beesTaxonomy(URL = "https://open.flinders.edu.au/ndownloader/files/47089969")
# Download the checklist
checklistFile <- beesChecklist(URL = "https://figshare.com/ndownloader/files/47092720")


##### 1.2 Taxonomy manipulation ####
# Extract author year, name, and number of authors on a species
taxonomyFile <- taxonomyFile %>% 
  # Extract the year of description
  dplyr::mutate(year = stringr::str_extract(authorship, "[0-9]{4}") %>%
                  as.numeric) %>%
  # Extract the authors
  dplyr::mutate(author = stringr::str_remove(authorship, "[0-9]{4}") %>% 
                  stringr::str_remove_all(., "\\(|\\)") %>%
                  stringr::str_squish() %>%
                  stringr::str_remove_all(",$") %>%
                  # Remove the "ands"
                  stringr::str_replace(" and", ", and") %>%
                  stringr::str_replace(",,", ",") %>%
                  stringr::str_remove_all(" and") %>%
                  # Remove the random numbers that have hung around
                  stringr::str_remove_all(", 1919|, 1900|, 1|-|935|, 7| sensu Robertson| et al") %>%
                  # Make a few things consistent
                  stringr::str_replace_all("E. A. B. ", "EAB ") %>%
                  stringr::str_replace_all(" &", ", ") %>%
                  stringr::str_replace_all("H. S. |H S ", "HS ") %>%
                  stringr::str_replace_all("J. R. ", "JR ") %>%
                  stringr::str_replace_all("M. C. |M C ", "MC ") %>%
                  stringr::str_replace_all("C. W. ", "CW ") %>%
                  stringr::str_replace_all("L. D. ", "LD ") %>%
                  stringr::str_replace_all("F. F. ", "FF ") %>%
                  stringr::str_replace_all("M. L. ", "ML ") %>%
                  stringr::str_replace_all("M. P. ", "MP ") %>%
                  stringr::str_replace_all("M. A. ", "MA ") %>%
                  stringr::str_replace_all("W. F. |W F ", "WF ") %>%
                  stringr::str_replace_all("R. P. |R P ", "R ") 
  ) %>% 
  # Get a count of author names per species
  dplyr::mutate(authorCount = stringr::str_count(author, ",")+1)

# Get this just for valid species 
taxonomy_valid <- taxonomyFile %>%
  dplyr::filter(taxonomic_status == "accepted")

# Get this just for synonyms species 
taxonomy_synonyms <- taxonomyFile %>%
  dplyr::filter(taxonomic_status == "synonym")


##### 1.3 Checklist manipulation ####
# Now, for each country pass on the newly-extracted info to the country checklist
checklistFile <- checklistFile %>% 
  dplyr::left_join(taxonomyFile %>%
                     dplyr::select(validName, year, author, authorCount),
                   by = "validName")


##### 1.4 Read BeeBDC occurrences ####
###### a. read data ####
# Read in the occurrence dataset
beeData <- readr::read_csv("/Users/jamesdorey/Desktop/Uni/My_papers/Bee_SDM_paper/Data_acquisition_workflow/Output/Intermediate/05_unCleaned_database_2024-02-15.csv",
                           col_types = BeeBDC::ColTypeR())
# Filter the data 
beeData <- beeData %>%
  BeeBDC::summaryFun(data = .,
                     # Don't filter for these
                     dontFilterThese = c(".gridSummary", ".lonFlag", ".latFlag", ".uncer_terms",
                                         ".sequential", ".year_outOfRange",
                                         # 
                                         ".GBIFflags",".coordinates_outOfRange",".sequential",
                                         ".val", ".cen", ".inst", ".equ",".rou", ".uncertaintyThreshold",
                                         ".sea", ".zer", ".otl", ".cap",  ".gbf", ".eventDate_empty",
                                         ".year_outOfRange", ".coordinates_empty"),
                     removeFilterColumns = TRUE,
                     filterClean = TRUE) %>%
  dplyr::mutate(country_suggested = dplyr::if_else(country == "MX",
                                                   "Mexico",
                                                   country_suggested))


##### 1.5 Literature correction ####
###### a. read lit occs ####
# Because those species without any occurrence records are known entities, but their sample 
# sizes are not, bring in data from taxonomic revisions to see if we can correct for this
# and apply the empirical distributions from a sample of papers across the unknown species
litRecords <- readxl::read_xlsx("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/DescriptionData/DiversityLiteratureRecords.xlsx",
                                col_types = "text") %>%
  # Fix some taxonomy 
  # Fill the scientificName if it is empty using species genus
  dplyr::mutate(scientificName = stringr::str_squish(scientificName),
                scientificName = dplyr::if_else(is.na(scientificName),
                                                stringr::str_c(genus, specificEpithet, sep = " "),
                                                scientificName)) %>%
  # Drop no names
  dplyr::filter(complete.cases(scientificName)) %>%
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     taxonomy = taxonomyFile,
                     path = getwd()) %>%
  # Find out if they are represented in the point data
  dplyr::mutate(hasPointRecords = dplyr::if_else(scientificName %in% unique(beeData$scientificName),
                                                 TRUE, FALSE))
# Turn these data into counts
litCounts <- litRecords %>%
  dplyr::group_by(scientificName) %>%
  dplyr::count() %>%
  dplyr::mutate(source = "litRecords",
                n = as.character(n))

###### b. read lit counts ####
# Read in the counts from literature data
# Counts from John Ascher
ascherCounts <- readxl::read_xlsx("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/DescriptionData/Ascher_comments_on_pointless.xlsx",
                                  col_types = "text") %>%
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     taxonomy = taxonomyFile,
                     path = getwd()) %>%
  # Drop an extra "n" column
  dplyr::select(!n) %>%
  # Rename columns
  dplyr::rename(n = `Minimum estimate from all sources`) %>%
  dplyr::select(n, scientificName) %>%
  dplyr::mutate(source = "Ascher2024")


# Read in counts from https://doi.org/10.1111/gcb.15879
fireCounts <- readxl::read_xlsx("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/DescriptionData/Combined_LifeHistory_Sep20.xlsx",
                                col_types = "text") %>%
  # Remove duplicate columns
  dplyr::select(!tidyselect::starts_with("Notes")) %>% 
  dplyr::select(!tidyselect::starts_with("REFERENCE")) %>%
  # Drop records that lack a count
  dplyr::filter(complete.cases(`Number_of_specimens_in_most-recent_description`)) %>%
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     speciesColumn = "Species_full",
                     taxonomy = taxonomyFile,
                     path = getwd()) %>%
  # Rename the columns of interest
  dplyr::rename(n = `Number_of_specimens_in_most-recent_description`,
                scientificName = Species_full) %>%
  # Combine the counts of males and females
  dplyr::group_by(scientificName) %>%
  dplyr::mutate(n = sum(as.numeric(n)) %>% as.character()) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(scientificName, .keep_all = TRUE) %>% 
  # Get the columns of interest
  dplyr::select(scientificName, n) %>%
  dplyr::mutate(source = "Dorey_etal_2021")

###### c. read from random species extractions ####

randomLit <- dplyr::bind_rows(
  readr::read_csv("CompletedTables/100RandomPoitnlessSpeecies_1_MCO.csv") %>% dplyr::mutate(Auth = "MCO"),
  readr::read_csv("CompletedTables/100RandomPoitnlessSpeecies_2_MCO.csv") %>% dplyr::mutate(Auth = "MCO"),
  readr::read_csv("CompletedTables/100RandomPoitnlessSpeecies_3_NJ_v3.csv") %>% dplyr::mutate(Auth = "NJ"),
  readr::read_csv("CompletedTables/100RandomPoitnlessSpeecies_4_DEG.csv") %>% dplyr::mutate(Auth = "DEG"),
  readr::read_csv("CompletedTables/100 Random pointsless species_5_AMG.csv") %>% dplyr::mutate(Auth = "AMG")
) %>%
  #dplyr::mutate(scientificName = scientificName %>% 
  #                  # Take the genus species name to be run thru HarmoniseR
  #                stringr::str_extract("[A-Za-z]+\\s[a-z]+")) %>% 
  # Only keep records where a sample size was found
  dplyr::filter(complete.cases(n)) %>% 
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     taxonomy = taxonomyFile,
                     path = getwd()) %>%
  dplyr::mutate(n = as.character(n))
  

  

###### c. combine lit #####
# Combine the lit record counts
combinedLit <- dplyr::bind_rows(litCounts, ascherCounts, fireCounts, randomLit) %>%
  # Drop NA counts
  dplyr::filter(!n == "NA") %>% 
  dplyr::mutate(n = n %>% as.numeric()) %>%
  dplyr::arrange(scientificName, desc(n)) %>%
  # Group by species name
  dplyr::group_by(scientificName) %>%
  # Take the highest n value from all sources
  dplyr::filter(dplyr::row_number() == 1) %>%
  dplyr::mutate(dataFrom = "literature") %>%
  # Remove species name not matching the taxonomy
  dplyr::filter(!scientificName =="Euryglossa rubricata")


###### d. quick plot ####
(TEST <- ggplot2::ggplot(data = combinedLit %>%
                           dplyr::group_by(scientificName) %>%
                           # Find out if they are represented in the point data
                           dplyr::mutate(hasPointRecords = dplyr::if_else(scientificName %in% unique(beeData$scientificName),
                                                                          TRUE, FALSE)),
                         ggplot2::aes(n, fill = hasPointRecords)) +
   ggplot2::geom_bar(position = "fill"))


##### 1.6 Prepare all counts ####
###### a. total counts ####
# Get counts of species in the dataset
beeData_counts <- beeData %>% 
  dplyr::group_by(scientificName) %>%
  dplyr::count() %>%
  dplyr::mutate(dataFrom = "points")


###### b. combine counts #####
# Combine combinedLit and beeData_counts
allActualCounts <- dplyr::full_join(combinedLit, beeData_counts, by = "scientificName",
                                    suffix = c("_lit", "_occ"))

# Get a complete overlap dataset (no NA)
overlap_actualCounts <- allActualCounts %>%
  dplyr::filter(complete.cases(n_lit)) %>%
  dplyr::filter(complete.cases(n_occ)) 

# Test the relationship
# Test for normality 
# NOT NORMAL
shapiro.test(overlap_actualCounts$n_lit)
shapiro.test(overlap_actualCounts$n_occ)
# Test for association
cor.test(overlap_actualCounts$n_lit, overlap_actualCounts$n_occ,
         alternative = "two.sided", method = "spearman")


# Plot the relationship between these datasets, where they overlap
(litOcc_plot <- ggplot2::ggplot(overlap_actualCounts,
                                ggplot2::aes(x = n_occ, y = n_lit)) +
    ggplot2::geom_point(ggplot2::aes(), colour = "darkgrey") +                                # scatter plot, coloured by sex
    ggplot2::labs(x = "Number of occurrences", y = "Number of literature records") +
    ggplot2::stat_smooth(
      method = "glm", ggplot2::aes(fill = "#FFC125", 
                                   colour = "black"
      )) +    # adding regression lines for each sex
    ggplot2::scale_colour_manual(values = c("#FFC125")) +
    ggplot2::scale_fill_manual(values = c("#FFC125")) +
    # ggplot2::ylim(c(0,max(overlap_actualCounts$n_occ))) +
    ggplot2::theme(#legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "transparent",
                                               colour = "black",
                                               linetype = NULL)) )
# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","1.6b_litOcc_plot.pdf"), 
                plot = litOcc_plot, 
                width = 6, height = 6, units = "in", dpi = 300)

###### c. lit only ####
# Get the literature ONLY counts
literature_only <- allActualCounts %>%
  dplyr::filter(is.na(n_occ))

(litHist <- ggplot2::ggplot(literature_only) +
    ggplot2::geom_histogram(stat = "bin", binwidth = 1,
                            ggplot2::aes(x = n_lit)) +
    ggplot2::xlab("Literature records") + ggplot2::ylab("Count") +
    ggplot2::theme(#legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "transparent",
                                               colour = "black",
                                               linetype = NULL)) )
# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","1.6c_litHistogram.pdf"), 
                plot = litHist, 
                width = 6, height = 6, units = "in", dpi = 300)

  ###### d. BUILD curve ####
  # Select only the relevant data
litCurveData <- literature_only %>%
  dplyr::group_by(n_lit) %>% 
  dplyr::mutate(count = dplyr::n())

  # Use mosaic to find the best linear model
f2 <- mosaic::fitModel( # 5.958
  count ~ (A*n_lit * n_lit^-log(B)), 
  #A*n_lit * n_lit^-sqrt(8),
  #A*n_lit + X + C *sqrt(n_lit),
  #start = list(A = 200, B = 10, C = 1),
  control = nls.control(maxiter = 5000,
                        minFactor = 0.000000001),
  data = litCurveData)
  # View the summary
summary(f2)
  # View the model parameters
f2


  # Build a plot of the points and the model
(litCountCurve <- ggformula::gf_point(
  count ~ n_lit, data = litCurveData) %>%
  mosaicCalc::slice_plot(f2(n_lit) ~ n_lit,
                         npts = 500, color = "#FFC125") +
  ggplot2::geom_point(ggplot2::aes(), colour = "darkgrey") +                                
  ggplot2::labs(x = "Number of occurrences", y = "Count")+    
  ggplot2::scale_colour_manual(values = c("#FFC125")) +
  ggplot2::scale_fill_manual(values = c("#FFC125")) +
  ggplot2::theme(#legend.position = "none",
    panel.background = ggplot2::element_rect(fill = NA,
                                             colour = "black",
                                             linetype = "solid"),
    panel.border = ggplot2::element_rect(fill = NA,
                                         colour = "black",
                                         linetype = "solid")) 
)

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","1.6c_litCountCurve.pdf"), 
                plot = litCountCurve, 
                width = 6, height = 6, units = "in", dpi = 300)



# Get the coordinates from the gam model in ggplot2 and put them in a tibble
litTibble <- dplyr::tibble(
  x_coords = ggplot2::ggplot_build(litCountCurve)$data[[2]]$x,
  y_coords = ggplot2::ggplot_build(litCountCurve)$data[[2]]$y)

# Set the seed for a consistent result
set.seed(1234)
# Get a random sample based on the empirical dataset
litSample <- sample(x = litTibble$x_coords, prob = litTibble$y_coords, 
                    replace = TRUE,
                    size = 10000) %>%
  round(0)


###### e. no occ. points #####
# Find the species that are not represented in the occurrence dataset
noPointSpecies <- taxonomy_valid %>% 
  dplyr::pull(validName) %>%
  setdiff(beeData$scientificName)  %>%
  # Add a sample size that is drawn from the same distriubtion of the literature and no-occurrences
  # dataset and generated in litSample
  dplyr::tibble(scientificName = .,
                n = litSample[1:length(.)],
                dataFrom = "randomisedSample")

# 
beeData_totalCounts <- beeData_counts %>% 
  dplyr::bind_rows(noPointSpecies)

  # This was completed and integrated above.
    #   #     # Extract 100 random species with no points for literature work
    #   ranomPointless <- noPointSpecies %>%
    #     # Remove the species that already have been counted
    #     dplyr::filter(!scientificName %in% literature_only$scientificName) %>%
    #     dplyr::slice_sample(n = 500) %>%
    #     dplyr::select(scientificName) %>%
    #     dplyr::mutate(n = NA_integer_,
    #                   doiORurl = NA_character_,
    #                   publicationYear = NA_integer_)

    #   # Split the ranomd 500 sepcies into a list of 100 rows each
    #   grps <- (split(ranomPointless, (seq(500)-1) %/% 100))
    #   # Loop to save these files for co-authors to work on
    #   for (i in seq_along(grps)) {
    #     write.csv(grps[[i]], paste0("100RandomPoitnlessSpeecies_", i, ".csv"))
    #   }


##### 1.7 Country counts ####
# Get the counts of species per country based on point occurrences
country_speciesCounts <- beeData %>%
  dplyr::mutate(country_suggested = dplyr::if_else(country_suggested == "United States",
                                                   "United States of America", 
                                                   country_suggested),
                country_suggested = dplyr::if_else(country_suggested == "Brunei Darussalam",
                                                   "Brunei Darussalam", country_suggested)) %>%
  dplyr::group_by(scientificName, country_suggested) %>%
  dplyr::count() %>%
  dplyr::mutate(dataFrom = "points") %>%
  dplyr::filter(!is.na(country_suggested))

# Harmonise country names with the checklist
checklistFile <- checklistFile %>%
  # Change country names
  dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "United States",
                                                    "United States of America", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Aland Islands",
                                                    "Åland Islands", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels",
                                                    "Belgium", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels",
                                                    "Belgium", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Swaziland",
                                                    "Eswatini", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "The Gambia",
                                                    "Gambia", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Republic of Congo",
                                                    "Republic of the Congo", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cape Verde",
                                                    "Republic of Cabo Verde", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "West Bank|Gaza",
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
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Saint-Martin",
                                                    "Saint Martin", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Somaliland",
                                                    "Somalia", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Wallis and Futuna Islands",
                                                    "Wallis and Futuna", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Faeroe Islands",
                                                    "Faroe Islands", rNaturalEarth_name)) 
#  # Check to see if checklist country names match
#TEST <- checklistFile %>% 
#  dplyr::select(rNaturalEarth_name) %>% 
#  dplyr::distinct() %>%
#  dplyr::mutate(checklist = "yes") %>%
#  dplyr::full_join(country_speciesCounts %>% 
#                     dplyr::ungroup() %>%
#                     dplyr::select(country_suggested) %>%
#                     dplyr::distinct() %>%
#                     dplyr::mutate(occurrences = "yes"),
#                   by = c("rNaturalEarth_name" = "country_suggested"))


# Make a function to sample the country sample sizes
countrySampler <- function(df = NULL){ 
  # Get the maximum value for litTibble$x_coords
  max_x <- max(litTibble$x_coords)
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
                                     sample(x = litTibble$x_coords[1:maxSample],
                                            prob = litTibble$y_coords[1:maxSample], 
                                            replace = TRUE,
                                            size = inputLength)[dplyr::row_number()] %>%
                                       round(0),
                                     # Or keep n
                                     n),
                  dataFrom = dplyr::if_else(is.na(dataFrom), "checklist", dataFrom))
  return(df)
}# END countrySampler

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

# Create the Chao input data including ONLY occurrence data
countryChaoData_occs <- country_speciesCounts %>%
  dplyr::select(scientificName, country_suggested, n) %>%
  tidyr::pivot_wider(names_from = country_suggested,
                     values_from = n,
                     values_fill = 0) %>%
  # Create the rownames
  tibble::column_to_rownames("scientificName") %>%
  dplyr::tibble()

  ##### 1.8 Save ####
  # Save these files 
readr::write_excel_csv(countryChaoData_checklist, file = "Table_outputs/1.8_countryChaoData_checklist.csv")
readr::write_excel_csv(country_speciesChecklistCounts, file = "Table_outputs/1.8_country_speciesChecklistCounts.csv")
readr::write_excel_csv(beeData_totalCounts, file = "Table_outputs/1.8_beeData_totalCounts.csv")


  # Read them in if needed
if(!exists("countryChaoData_checklist")){
  countryChaoData_checklist <- readr::read_csv("Table_outputs/1.8_countryChaoData_checklist.csv")
}
if(!exists("country_speciesChecklistCounts")){
  country_speciesChecklistCounts <- readr::read_csv("Table_outputs/1.8_country_speciesChecklistCounts.csv")
}
if(!exists("beeData_totalCounts")){
  beeData_totalCounts <- readr::read_csv("Table_outputs/1.8_beeData_totalCounts.csv")
}


#### 2.0 Statistics checklist ####
##### 2.1 Taxonomy accumulations ####
###### a. accepted names ####
# Make a dataset showing count of species per year according to the taxonomy
validCounts <- taxonomy_valid %>%
  dplyr::group_by(year) %>% 
  dplyr::count() %>%
  dplyr::arrange(year) %>%
  dplyr::ungroup() %>%
  # Add the cunmulative sum
  dplyr::mutate(cumulative = cumsum(n))

# Make a dataset showing count of species per year according to when the specimen was first found
occYear_counts <- beeData %>%
  dplyr::select(scientificName, year) %>%
  dplyr::arrange(year) %>%
  dplyr::group_by(scientificName) %>%
  dplyr::distinct(scientificName, .keep_all = TRUE) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(n = 1) %>%
  # Add the cunmulative sum
  dplyr::mutate(cumulative = cumsum(n))

# Quick and dirty histogram of species accumulation
(TaxoOccPlot <- ggplot2::ggplot(validCounts,
                ggplot2::aes(x = year)) + 
  ggplot2::geom_line(ggplot2::aes(y = cumulative, colour = "black")) +
  # Occurrence accumulation from first collection
  ggplot2::geom_line(data = occYear_counts, ggplot2::aes(y = cumulative, colour = "red")) +
  ggplot2::scale_x_continuous("Year", breaks = seq(from = 1760, to = 2024, by = 20), 
                              limits = c(1758, 2024)) +
  ggplot2::ylab("Cumulative species") + 
  ggplot2::theme(#legend.position = "none",
    panel.background = ggplot2::element_rect(fill = "transparent",
                                             colour = "black",
                                             linetype = NULL),
    panel.border = ggplot2::element_rect(fill = "transparent",
                                         colour = "black",
                                         linetype = NULL),
    legend.position = "inside",
    legend.position.inside = c(0.15,0.8),
    legend.key=element_blank()) +
  scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
  scale_colour_manual(name = 'Data source', 
                      values =c('black','red'), labels = c('Taxonomy','First occurrence')))

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","2.1ac_TaxoOccCurve.pdf"), 
                plot = TaxoOccPlot, 
                width = 8, height = 4, units = "in", dpi = 300)

###### b. synonyms ####
# Make a dataset showing count of species per year
INvalidCounts <- taxonomy_synonyms %>%
  dplyr::group_by(year) %>% 
  dplyr::count() %>%
  dplyr::arrange(year) %>%
  dplyr::ungroup() %>%
  # Add the cunmulative sum
  dplyr::mutate(cumulative = cumsum(n))

# Quick and dirty histogram of species accumulation
(synonymAccumPlot <- ggplot2::ggplot(INvalidCounts,
                ggplot2::aes(x = year)) + 
  ggplot2::geom_line(data = validCounts, ggplot2::aes(y = cumulative, colour = "black")) +
  ggplot2::geom_line(data = INvalidCounts, ggplot2::aes(y = cumulative, colour = "red")) +
    ggplot2::scale_x_continuous("Year", breaks = seq(from = 1760, to = 2024, by = 20), 
                                limits = c(1758, 2024)) +
    ggplot2::ylab("Cumulative species") + 
    ggplot2::theme(#legend.position = "none",
      panel.background = ggplot2::element_rect(fill = "transparent",
                                               colour = "black",
                                               linetype = NULL),
      panel.border = ggplot2::element_rect(fill = "transparent",
                                           colour = "black",
                                           linetype = NULL),
      legend.position = "inside",
      legend.position.inside = c(0.15,0.8),
      legend.key=element_blank()) +
    scale_fill_identity(name = 'the fill', guide = 'legend',labels = c('m1')) +
    scale_colour_manual(name = 'Species name', 
                        values =c('black','red'), labels = c('New valid name','New synonym')))

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","2.1b_SynonymCurve.pdf"), 
                plot = synonymAccumPlot, 
                width = 8, height = 4, units = "in", dpi = 300)

  # Combine the plots 
(accumCurves <- cowplot::plot_grid(TaxoOccPlot,
                               synonymAccumPlot, 
                               labels = c("(a)","(b)"),
                               ncol = 1, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","2.1b_Curves.pdf"),
                   plot = accumCurves,
                   base_width = 8,
                   base_height = 8)


###### c. USA names ####
US_counts <- checklistFile %>%
  dplyr::filter(rNaturalEarth_name == "United States of America") %>%
  dplyr::group_by(year) %>% 
  dplyr::count() %>%
  dplyr::arrange(year) %>%
  dplyr::ungroup() %>%
  # Add the cunmulative sum
  dplyr::mutate(cumulative = cumsum(n))

# Quick and dirty histogram of species accumulation
ggplot2::ggplot(US_counts,
                ggplot2::aes(x = year)) + 
  ggplot2::geom_line(ggplot2::aes(y = cumulative)) +
  ggplot2::scale_x_continuous("Year", breaks = seq(from = 1760, to = 2022, by = 20), 
                              limits = c(1758, 2022)) +
  ggplot2::ylab("Cumulative US described species") + 
  ggplot2::theme(#legend.position = "none",
    panel.background = ggplot2::element_rect(fill = "transparent",
                                             colour = "black",
                                             linetype = NULL)) 


##### 2.2 Chao ####
###### a. global records ####
# Get the Chao values for global bee species diversity
globalChao <- SpadeR::ChaoSpecies(beeData_totalCounts$n,
                                  datatype = "abundance", 
                                  k = 10, conf = 0.95)

  # Save this file
base::saveRDS(globalChao, 
              file = "2.2a_globalChao.Rda")

###### b. country records ####
countryChao_n1 <- ChaoWrapper(data = countryChaoData_checklist,
                                   k = 5,
                                   datatype = "abundance",
                                   conf = 0.95,
                                   mc.cores = 10)

# Save this file
base::saveRDS(countryChao_n1, 
              file = "2.2a_countryChao.Rda")


##### 2.3 iNEXT ####
  ###### a. global records ####
global_iNEXT <- beeData_totalCounts %>%
  dplyr::select(scientificName, n) %>%  # Create the rownames
  tibble::column_to_rownames("scientificName") %>%
  dplyr::tibble()

global_iNEXT_out <- iNEXT::iNEXT(x = global_iNEXT$n, datatype = "abundance",
                                 # 0 for species richness Shannon diversity (q = 1, the exponential of Shannon entropy) and Simpson diversity (q = 2, the inverse of Simpson concentration).
                                 q = 0)

# Save this file
base::saveRDS(global_iNEXT_out, 
              file = "2.3a_globaliNEXT.Rda")

iNEXT::ggiNEXT(global_iNEXT_out, type=1, facet.var="Assemblage")


  ###### b. country records ####
source("iNEXTwrapper.R")
  # Run iNEXT for each country using the iNEXT wrapper
country_iNEXT <- iNEXTwrapper(data = country_speciesChecklistCounts,
                              variableColumn = "rNaturalEarth_name",
                              valueColumn = "n",
                              datatype = "abundance",
                              mc.cores = 10)

# Save this file
base::saveRDS(country_iNEXT, 
              file = "2.3b_country_iNEXT.Rda")
  # Read back in if needed
if(!exists(country_iNEXT)){
  country_iNEXT <- base::readRDS("2.3b_country_iNEXT.Rda")
}

country_iNEXT$iNextEst$iNextEst

  # Get the sample sizes for each country ( get those < 30)
smallSampleSizes <- country_iNEXT$DataInfo %>%
  dplyr::arrange(n) %>%
  dplyr::filter(n < 30)

# Plot all country plots using ggiNEXTwrapper
source("ggiNEXTwrapper.R")
ggiNEXTwrapper(data = country_iNEXT,
                # Remove the countries under a certain sample size
               filterOut = smallSampleSizes %>%
                 dplyr::pull(Assemblage),
               legendPerPlot = FALSE,
               nrow = 4,
               ncol = 3,
               labels = NULL,
               fileName = "iNEXTplots",
               outPath = "/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
               base_width = 8.3,
               base_height = 11.7, 
               dpi = 300)


##### 2.4 Breakaway ####
  # Make sure you have the dataset
if(!exists("beeData_totalCounts")){
  beeData_totalCounts <- readr::read_csv("Table_outputs/1.8_beeData_totalCounts.csv")
}

  # Estimate species richness using the Bayesian package, breakaway.
library(breakaway)
data(toy_otu_table)
otu_data <- toy_otu_table

frequencytablelist <- build_frequency_count_tables(otu_data)

  # Transform the data into the required format
totalCounts_breakaway <- beeData_totalCounts %>%
  dplyr::arrange(n) %>%
  dplyr::group_by(n) %>%
  dplyr::count(n, name = "count") %>%
  dplyr::rename(index = "n",
                frequency = "count")

breakaway(totalCounts_breakaway,
          cutoff = 6)

plot(breakaway(totalCounts_breakaway))





