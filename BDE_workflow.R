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
install.packages("lmerTest")
install.packages("mosaicCalc")
install.packages("geodata")
install.packages("tidystringdist")
install.packages("igraph")
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
lapply(c("BeeBDC", "magrittr", "bdc", "SpadeR", "ggplot2", "dplyr", "phyloseq", "breakaway",
         "geodata", "lmerTest"), 
       library, character.only = TRUE)

# Load in an extra wrapper around SpadeR
source("ChaoWrapper.R")
source("ChaoWrapper.R")
source("iNEXTwrapper.R")
source("richnessSampleR.R")
source("countryHarmoniseR.R")


  ##### 0.3 Re-read data ####
  # Read in the saved datasets that may be needed in this script

if(!exists("beeData_counts")){
  beeData_counts <- readr::read_csv(paste0("Table_outputs/","1.6a_beeData_counts.csv"))}
if(!exists("litCurveData")){
  litCurveData <- readr::read_csv(paste0("Table_outputs/","1.6d_litPoints.csv"))}
if(!exists("country_speciesCounts")){
  country_speciesCounts <- readr::read_csv("Table_outputs/1.7_country_speciesCounts.csv")}
if(!exists("countryChaoData_checklist")){
  countryChaoData_checklist <- readr::read_csv("Table_outputs/1.9_countryChaoData_checklist.csv")}
if(!exists("country_speciesChecklistCounts")){
  country_speciesChecklistCounts <- readr::read_csv("Table_outputs/1.9_country_speciesChecklistCounts.csv")}
if(!exists("beeData_totalCounts")){
  beeData_totalCounts <- readr::read_csv("Table_outputs/1.9_beeData_totalCounts.csv")}
if(!exists("continentWider")){
  continentWider <- readr::read_csv("Table_outputs/1.9_continentWider.csv")}
if(!exists("taxonomy_valid")){
  taxonomy_valid <- readr::read_csv("Table_outputs/2.1a_taxonomy_valid.csv",
                                    guess_max = 10000)}
if(!exists("validCounts")){
  validCounts <- readr::read_csv("Table_outputs/2.1a_validCounts.csv")}
if(!exists("INvalidCounts")){
  INvalidCounts <- readr::read_csv("Table_outputs/2.1a_INvalidCounts.csv")}
if(!exists("occYear_counts")){
  occYear_counts <- readr::read_csv("Table_outputs/2.1a_occYear_counts.csv")}
if(!exists("globalChao")){
  globalChao <- readRDS("2.2a_globalChao.Rda")}
if(!exists("countryChao_n1")){
  countryChao_n1 <- readRDS("2.2b_countryChao.Rda")}
if(!exists("continentChao_n1")){
  continentChao_n1 <- readRDS("2.2c_continentChao.Rda")}
if(!exists("country_iNEXT")){
  country_iNEXT <- base::readRDS("2.3b_country_iNEXT.Rda")}
if(!exists("global_iNEXT_out")){
  global_iNEXT_out <- readRDS("2.3a_globaliNEXT.Rda")}
if(!exists("country_iNEXT")){
  country_iNEXT <- base::readRDS("2.3b_country_iNEXT.Rda")}
if(!exists("continent_iNEXT")){
  continent_iNEXT <- base::readRDS("2.3c_continent_iNEXT.Rda")}
if(!exists("continent_summary")){
  continent_summary <- readr::read_csv("Figure_outputs/country_iNEXT/2.3c_continentTable.csv")}
if(!exists("richnessInputs")){
  richnessInputs <- base::readRDS("2.4a_richnessInputs.Rda")}
# if(!exists("globalSampled")){
#   globalSampled <- base::readRDS("2.4b_globalSampled.Rda")}
if(!exists("countrySampled")){
  countrySampled <- base::readRDS("2.4c_countrySampled.Rda")}
if(!exists("continentSampled")){
  continentSampled <- base::readRDS("2.4d_continentSampled.Rda")}
if(!exists("combinedStatistics")){
  combinedStatistics <- readr::read_csv("Table_outputs/3.1a_combinedStatistics.csv")}
if(!exists("country_summary")){
  country_summary <- readr::read_csv(paste0("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
                                            "/2.3b_countryTable.csv"))}
if(!exists("countryElevations_var")){
  countryElevations_var <- readr::read_csv("CorrelatesData/countryElevations_var.csv")}
if(!exists("countryDistances")){
  countryDistances <- readr::read_csv("Table_outputs/5.1e_countryDistances.csv")}
if(!exists("combined_explanatory")){
  combined_explanatory <- readr::read_csv("CorrelatesData/combined_explanatory.csv")}
if(!exists("beeData")){
  BeeData <- readRDS("beeData.Rda")
}



#### 1.0 Data prep ####
##### 1.1 Download data ####
# Download the taxonomy
taxonomyFile <- BeeBDC::beesTaxonomy()
# Download the checklist
checklistFile <- BeeBDC::beesChecklist()


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
                           col_types = BeeBDC::ColTypeR()) %>%
  dplyr::mutate(country_suggested = dplyr::if_else(country == "MX",
                                                   "Mexico",
                                                   country_suggested))
# Filter the data 
beeData <- beeData %>%
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     taxonomy = taxonomyFile,
                     path = getwd(),
                     checkVerbatim = TRUE) %>%
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
# Save the dataset
base::saveRDS(beeData, 
              file = "beeData.Rda")
  # Read back in if it's not already in the environment 
if(!exists("beeData")){
  readRDS("beeData.Rda")
}

readr::write_excel_csv(beeData %>% 
                         dplyr::filter(country %in% c("Fiji", "Uganda", "Vietnam", "Zambia")) %>%
                         dplyr::select(c("scientificName", "country_suggested")),
                       "/Users/jamesdorey/Desktop/Uni/Packages/BeeBDC_development/beesCountry.csv")


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

randomOnly <- dplyr::bind_rows(ascherCounts, randomLit) %>%
  readr::write_csv(paste0("Table_outputs/","1.6c_randomLitTogether.csv"))
  

###### d. combine lit #####
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

combinedLit%>% 
  # Run through harmoniseR
  BeeBDC::harmoniseR(data = .,
                     taxonomy = taxonomyFile,
                     path = getwd())


###### e. quick plot ####
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

  # Get counts by species and country (for the occCurve2)
beeData_counts_spCountry <- beeData %>% 
  dplyr::group_by(scientificName, country_suggested) %>%
  dplyr::count() %>%
  dplyr::mutate(dataFrom = "points") %>%
  dplyr::filter(is.na(country_suggested))

readr::write_csv(beeData_counts, paste0("Table_outputs/","1.6a_beeData_counts.csv"))


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

  # Save this file 
readr::write_excel_csv(litCurveData, paste0("Table_outputs/","1.6d_litPoints.csv"))
if(!exists("litCurveData")){
  litCurveData <- readr::read_csv(paste0("Table_outputs/","1.6d_litPoints.csv"))
}

  # Use mosaic to find the best linear model
f2 <- mosaic::fitModel( # 5.958
  count ~ (A*n_lit * n_lit^-log(B)), 
  #A*n_lit * n_lit^-sqrt(8),
  #A*n_lit + X + C *sqrt(n_lit),
  start = list(A = 200, B = 10),
  control = nls.control(maxiter = 5000,
                        minFactor = 0.000000001),
  data = litCurveData)
  # View the summary# View tlog2()he summary
summary(f2)
  # View the model parameters
f2

  # Look at the points on the below figure
SummaryOfSampleSize <- litCurveData %>% 
  dplyr::group_by(n_lit) %>%
  dplyr::count()

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

  # Extract the occurrence only records
allOccCounts <- allActualCounts %>% 
  dplyr::filter(complete.cases(n_occ)) %>%
  dplyr::group_by(n_occ) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::distinct(n_occ, count) %>% 
  dplyr::mutate(count = (count*(max(litCurveData$count)/max(.$count)))
                )

readr::write_excel_csv(allOccCounts, paste0(RootPath, "/Table_outputs/1.6d_allOccCounts.csv"))



# Use mosaic to find the best linear model of the OCCURRENCE data
f_occ <- mosaic::fitModel( # 5.958
  count ~ (A+(B/(n_occ))), 
  #A*n_lit * n_lit^-sqrt(8),
  #A*n_lit + X + C *sqrt(n_lit),
  #start = list(A = 200, B = 10),
  control = nls.control(maxiter = 5000,
                        minFactor = 0.000000001),
  data = allOccCounts )
# View the summary# View tlog2()he summary
summary(f_occ)


# Build a plot of the points and the model
(litOccCurve <-ggplot2::ggplot(data = allOccCounts,
                        ggplot2::aes(x = n_occ, y = count)) +
    ggplot2::geom_function(fun = function(x) (-0.1891 + (226.4555/(x))), 
                           ggplot2::aes(color = "#3FB8AF"), linetype = 1, size = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = "#3FB8AF"), 
                        alpha = 0.65) +
    ggplot2::labs(x = "Number of occurrences", y = "Count")+ 
    
    ggplot2::geom_function(fun = function(x) (228.7531 * x * x^-log(12.1593)), 
                           ggplot2::aes(colour = "#FFC125"), linetype = 1, size = 1) +
    ggplot2::geom_point(data = litCurveData,
                        ggplot2::aes(x = n_lit, y = count, colour = "#FFC125"),
                                     alpha = 0.65) +
    ggplot2::scale_colour_manual(values = c("#3FB8AF", "#FFC125"),
                                 label = c("Occurrences", "Literature"),
                                 name = "Data source") +
      # Limit the x-axis to see the curve
    ggplot2::xlim(1,100) + 
    ggplot2::theme(#legend.position = "none",
      panel.background = ggplot2::element_rect(fill = NA,
                                               colour = "black",
                                               linetype = "solid"),
      panel.border = ggplot2::element_rect(fill = NA,
                                           colour = "black",
                                           linetype = "solid")) 
)
# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","1.6c_litOccCurve.pdf"), 
                plot = litOccCurve, 
                width = 6, height = 6, units = "in", dpi = 300)


  # Now build this curve for the country-species combinations
beeData_counts_spCountry
# Extract the occurrence only records
allOccCounts_spCountry <- beeData_counts_spCountry %>% 
  dplyr::rename(n_occ_spCou = n) %>% 
  dplyr::group_by(n_occ_spCou) %>% 
  dplyr::mutate(count = dplyr::n()) %>% 
  dplyr::distinct(n_occ_spCou, count) %>% 
  dplyr::mutate(count = (count*(max(litCurveData$count)/max(.$count)))
  )



# Use mosaic to find the best linear model of the OCCURRENCE data
f_occ_spCou <- mosaic::fitModel( # 5.958
  count ~ (A+(B/(n_occ_spCou))), 
  #A*n_lit * n_lit^-sqrt(8),
  #A*n_lit + X + C *sqrt(n_lit),
  #start = list(A = 200, B = 10),
  control = nls.control(maxiter = 5000,
                        minFactor = 0.000000001),
  data = allOccCounts_spCountry )
# View the summary# View tlog2()he summary
summary(f_occ_spCou)

# Build a plot of the points and the model
(litOccCurve_spCou <- ggplot2::ggplot(data = allOccCounts,
                               ggplot2::aes(x = n_occ, y = count)) +
    # Literature curve data 
    ggplot2::geom_function(fun = function(x) (228.7531 * x * x^-log(12.1593)), 
                           ggplot2::aes(colour = "#FFC125"), linetype = 1, size = 1) +
    ggplot2::geom_point(data = litCurveData,
                        ggplot2::aes(x = n_lit, y = count, colour = "#FFC125"),
                        alpha = 0.65) +
      # Occ global data
    ggplot2::geom_function(fun = function(x) (-0.1891 + (226.4555/(x))), 
                           ggplot2::aes(color = "#1BB6AFFF"), linetype = 1, size = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = "#1BB6AFFF"), 
                        alpha = 0.65) +
      # Occ species and country data 
    ggplot2::geom_function(fun = function(x) (-0.68833 + (217.51375/(x))), 
                           ggplot2::aes(color = "#172869FF"), linetype = 1, size = 1) +
    ggplot2::geom_point(ggplot2::aes(colour = "#172869FF"), 
                        alpha = 0.65) +
    ggplot2::scale_colour_manual(values = c("#FFC125", "#1BB6AFFF", "#172869FF"),
                                 label = c("Literature",
                                           "Occurrences global",
                                           "Occurrences country"),
                                 name = "Data source") +
    ggplot2::labs(x = "Number of occurrences", y = "Count")+ 
    
    # Limit the x-axis to see the curve
    ggplot2::xlim(1,100) + 
    ggplot2::theme(legend.position.inside = c(0.76, 0.86),
                   legend.position = "inside",
      panel.background = ggplot2::element_rect(fill = NA,
                                               colour = "black",
                                               linetype = "solid"),
      panel.border = ggplot2::element_rect(fill = NA,
                                           colour = "black",
                                           linetype = "solid")) 
)
# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","1.6c_litOccCurve_speciesCountry.pdf"), 
                plot = litOccCurve_spCou, 
                width = 6, height = 6, units = "in", dpi = 300)

# Get the coordinates from the gam model in ggplot2 and put them in a tibble
litTibble <- dplyr::tibble(
  x_coords = ggplot2::ggplot_build(litCountCurve)$data[[2]]$x,
  y_coords = ggplot2::ggplot_build(litCountCurve)$data[[2]]$y)

# Set the seed for a consistent result
set.seed(1234)
# Get a random sample based on the empirical dataset
if(TRUE) {
  stop("Be certain that you want to run the below line. It has the power to change the results slightly as it's a random sample.")
}
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

symdiff(beeData_totalCounts$scientificName,
        taxonomy_valid$validName)

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

readr::write_excel_csv(country_speciesCounts, "Table_outputs/1.7_country_speciesCounts.csv")
if(!exists("country_speciesCounts")){
  country_speciesCounts <- readr::read_csv("Table_outputs/1.7_country_speciesCounts.csv")
}


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


  ##### 1.8 Continent ####
# Download a world map to convert countries to continents
worldMap <- rnaturalearth::ne_countries(returnclass = "sf",
                                        scale = 50, type = "countries") 

# Turn the country occ data into a continent one
continentOccs <- country_speciesCounts %>%
  dplyr::ungroup() %>%
    # Change some country names to better match the continent data
  dplyr::mutate(country_suggested = dplyr::if_else(country_suggested == "Eswatini",
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


  # Turn the country checklist data into a continent one
continentChecklist <- beesChecklist(URL = "https://figshare.com/ndownloader/files/47092720") %>% 
  dplyr::select(validName, rNaturalEarth_name) %>%
dplyr::ungroup() %>%
  # Change some country names to better match the continent data
  dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Barbuda",
                                                   "Antigua and Barbuda", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels",
                                                   "Belgium", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Caribbean Netherlands",
                                                   # Renamed Curaçao to match to continent
                                                   "Curaçao", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cocos Islands",
                                                   "Indian Ocean Territories", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Federation of Bosnia and Herzegovina",
                                                   "Bosnia and Herzegovina", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "French Guiana",
                                                   "Brazil", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gaza",
                                                   "Palestine", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gibraltar",
                                                   "Spain", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Guadeloupe",
                                                   "Barbados", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Martinique",
                                                   "Curaçao", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Mayotte",
                                                   "Madagascar", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name,"Cyprus"),
                                                   "N. Cyprus", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Réunion",
                                                   "Mauritius", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Svalbard Islands",
                                                   "Norway", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "West Bank",
                                                   "Palestine", rNaturalEarth_name)) %>%
  # Join first by name
  dplyr::left_join(worldMap %>% dplyr::select(name, continent) %>%
                     sf::st_drop_geometry(), by = c("rNaturalEarth_name" = "name")) %>%
  # Then join by name_long
  dplyr::left_join(worldMap %>% dplyr::select(name_long, continent) %>%
                     sf::st_drop_geometry(), by = c("rNaturalEarth_name" = "name_long")) %>%
  # Now merge these continent columns
  dplyr::mutate(continent = dplyr::if_else(is.na(continent.x),
                                           continent.y, continent.x)) %>%
  # drop interim and old count columns
  dplyr::select(!c("continent.x", "continent.y")) %>%
  dplyr::mutate(name_continent = stringr::str_c(validName, continent, sep = "__"))

  # Find the names in the checklist that are missing from the occurrence records
  countryNOTcontinent <- symdiff(continentOccs$name_continent, continentChecklist$name_continent) %>%
      # Extract continent and species names into tibble columns and remove the . column
    dplyr::tibble(scientificName = stringr::str_extract(., ".*__") %>% stringr::str_remove("__"),
                  continent = stringr::str_extract(., "__.*") %>% stringr::str_remove("__"),
                  n = NA,
                  dataFrom = "litEstimate") %>%
    dplyr::select(!.) 
  
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
  
  # Pivot wider 
  continentWider <- continentCounts %>%
    dplyr::select(scientificName, continent, n) %>%
    tidyr::drop_na() %>% 
    tidyr::pivot_wider(names_from = continent,
                       values_from = n,
                       values_fill = 0) %>%
    # Create the rownames
    tibble::column_to_rownames("scientificName") %>%
    dplyr::tibble()
    
  ##### 1.9 Save ####
  # Save these files 
readr::write_excel_csv(countryChaoData_checklist, file = "Table_outputs/1.9_countryChaoData_checklist.csv")
readr::write_excel_csv(country_speciesChecklistCounts, file = "Table_outputs/1.9_country_speciesChecklistCounts.csv")
readr::write_excel_csv(beeData_totalCounts, file = "Table_outputs/1.9_beeData_totalCounts.csv")
readr::write_excel_csv(continentCounts, file = "Table_outputs/1.9_continentCounts.csv")
readr::write_excel_csv(continentWider, file = "Table_outputs/1.9_continentWider.csv")


  # Read them in if needed
if(!exists("countryChaoData_checklist")){
  countryChaoData_checklist <- readr::read_csv("Table_outputs/1.9_countryChaoData_checklist.csv")
}
if(!exists("country_speciesChecklistCounts")){
  country_speciesChecklistCounts <- readr::read_csv("Table_outputs/1.9_country_speciesChecklistCounts.csv")
}
if(!exists("beeData_totalCounts")){
  beeData_totalCounts <- readr::read_csv("Table_outputs/1.9_beeData_totalCounts.csv")
}
if(!exists("continentWider")){
  continentWider <- readr::read_csv("Table_outputs/1.9_continentWider.csv")
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
# Make a dataset showing count of species per year
INvalidCounts_pre <- taxonomy_synonyms %>%
    # Simplify the dataframe
  dplyr::select(accid, id, species, authorship, year) %>%
  # Remove numbers from the species name column
  dplyr::mutate(species = species %>% stringr::str_remove_all(., "[0-9]")) %>% 
    # Start by making columns to identify authographic variants (rather than synonyms)
  # Remove brackets
  dplyr::mutate(authorYear_brackRM = stringr::str_remove_all(authorship, "\\(|\\)"),
                .after = authorship) %>%
  dplyr::group_by(accid, authorYear_brackRM) %>%
    # Count rows > 1 = a match
  dplyr::mutate(grpNum = dplyr::row_number(), .after = authorship,
                grpSize = dplyr::n()) 

  # Get the soundex matches for invalid names 
iGraphInvalid <- INvalidCounts_pre %>%
  dplyr::filter(grpSize > 1) %>% 
    # Re-combine these data with themselves to get all combinations
    # within accid and authorYear_brackRM
  dplyr::left_join(., INvalidCounts_pre %>%
                     dplyr::select(accid, id, authorYear_brackRM, species) %>%
                     dplyr::rename(species2 = species,
                                   id2 = id),
                   by = c("accid" = "accid", "authorYear_brackRM" = "authorYear_brackRM"),
                   multiple = "all", relationship = "many-to-many") %>%
    # Remove duplicates
  dplyr::distinct() %>%
  dplyr::relocate(species, .before = species2) %>% 
  # Compare the distances within the species name column 
  tidystringdist::tidy_stringdist(v1 = species, v2 = species2) %>%
    # Remove self-matches
  dplyr::filter(!id == id2) %>%
  dplyr::filter(soundex == 0) %>%
  dplyr::filter(!species %in% c("sp","sp.")) %>% 
  dplyr::select(id, id2, species, species2, authorYear_brackRM) %>%
    # Use igraph now 
  igraph::graph_from_data_frame(., directed = FALSE) %>% 
  igraph::components()
    # Make this into a tibble with id and membership
grouped_ids <- dplyr::tibble(
  id = names(iGraphInvalid$membership) %>% as.numeric(),
  group = iGraphInvalid$membership
)

  # Add in the group numbers and then label everything accordingly for future grouping and selection
INvalidCounts <- INvalidCounts_pre %>%
  dplyr::left_join(., grouped_ids,
                   by = "id") %>%
  dplyr::ungroup() %>% 
  dplyr::mutate(groupInclusive = dplyr::if_else(is.na(group),
                                                      stringr::str_c("syn_", dplyr::row_number()),
                                                stringr::str_c("Orth_", group %>% as.character()))) %>% 
  dplyr::group_by(groupInclusive) %>% 
    # Now take the first of every group
  dplyr::filter(dplyr::row_number() == 1) %>% 
  dplyr::filter(stringr::str_count(species) > 1) %>%
    # Make counts
  dplyr::group_by(year) %>% 
  dplyr::count() %>%
  dplyr::arrange(year) %>%
  dplyr::ungroup() %>%
  # Add the cunmulative sum
  dplyr::mutate(cumulative = cumsum(n))


readr::write_csv(taxonomy_valid, "Table_outputs/2.1a_taxonomy_valid.csv")
readr::write_csv(validCounts, "Table_outputs/2.1a_validCounts.csv")
readr::write_csv(INvalidCounts, "Table_outputs/2.1a_INvalidCounts.csv")


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

readr::write_csv(occYear_counts, "Table_outputs/2.1a_occYear_counts.csv")

# Quick and dirty histogram of species accumulation
(TaxoOccPlot <- ggplot2::ggplot(validCounts,
                ggplot2::aes(x = year)) + 
    # Occurrence accumulation from synonym
    ggplot2::geom_line(data = INvalidCounts, ggplot2::aes(y = cumulative, colour = "Syn"),
                       linewidth = 1.2) +
        # First occurrence 
    ggplot2::geom_line(data = occYear_counts, ggplot2::aes(y = cumulative, colour = "Occ"),
                       linewidth = 1.2) +
      # New valid name
  ggplot2::geom_line(data = validCounts, ggplot2::aes(y = cumulative, colour = "Valid"),
                     linewidth = 1.2) +
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
                      breaks =  c('Syn',"Valid",'Occ'),
                      values =c('#d62828','#f77f00', "#fcbf49"), 
                      labels = c('New synonym', 'New valid name',"First occurrence record")))

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","2.1ac_TaxoOccCurve.pdf"), 
                plot = TaxoOccPlot, 
                width = 8, height = 4, units = "in", dpi = 300)

###### b. synonyms ####

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
                               labels = c("(A)","(B)"),
                               ncol = 1, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","2.1b_Curves.pdf"),
                   plot = accumCurves,
                   base_width = 8,
                   base_height = 8)

  ###### c. rate of increase ####
  # Extract the rate of increase for the full year range 
CompleteEnds <- validCounts %>% 
  dplyr::mutate(meanIncrease = mean(n),
                medianIncrease = median(n))
  # 90 per year average, 78 median

# Extract the rate of increase from 1960
endRange <- validCounts %>% 
  dplyr::filter(year > 1959)  %>% 
  dplyr::mutate(meanIncrease = mean(n),
                medianIncrease = median(n))

# 117 per year
# Year range to describe this many species
    #   > 3884/ 117
    #   [1] 33.19658
    #   > 4839/ 117
    #   [1] 41.35897

Aus2000 <- readxl::read_excel("Papers_since_2000.xlsx") %>% 
  dplyr::group_by(Phylogeny_used) %>%
  dplyr::mutate(NumSp = sum(Number_of_new_species),
                numSyn = sum(Num_synonyms, na.rm = TRUE),
                numTotal = sum(Total_taxa, na.rm = TRUE),
                numArticles = dplyr::n()) %>% 
  dplyr::distinct(Phylogeny_used, .keep_all = TRUE)
  # What percentage of articles use:
# Phylogenetic techniques
(Aus2000[Aus2000$Phylogeny_used == TRUE,]$numArticles / 
  sum(Aus2000[Aus2000$Phylogeny_used == TRUE,]$numArticles,
      Aus2000[Aus2000$Phylogeny_used == FALSE,]$numArticles)) * 100 

  # What percentage of species are described using:
# Phylogenetic techniques
(Aus2000[Aus2000$Phylogeny_used == TRUE,]$NumSp / 
    sum(Aus2000[Aus2000$Phylogeny_used == TRUE,]$NumSp,
        Aus2000[Aus2000$Phylogeny_used == FALSE,]$NumSp)) * 100 


###### d. USA names ####
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

  # Extract the relevant data
globalChao_iChao1 <- globalChao$Species_table %>%
  as.data.frame() %>% 
  tibble::rownames_to_column() %>%
  dplyr::filter(stringr::str_detect(rowname, "iChao1"))

###### b. country records ####
countryChao_n1 <- ChaoWrapper(data = countryChaoData_checklist,
                                   k = 5,
                                   datatype = "abundance",
                                   conf = 0.95,
                                   mc.cores = 10)

# Save this file
base::saveRDS(countryChao_n1, 
              file = "2.2b_countryChao.Rda")

###### c. continental ####
  # Run species estimates per continent, one core per continent
continentChao_n1 <- ChaoWrapper(data = continentWider,
                              k = 5,
                              datatype = "abundance",
                              conf = 0.95,
                              mc.cores = 8)

# Save this file
base::saveRDS(continentChao_n1, 
              file = "2.2c_continentChao.Rda")

    ###### d. read in if needed ####
if(!exists("globalChao")){
  globalChao <- readRDS("2.2a_globalChao.Rda")
}
if(!exists("countryChao_n1")){
  countryChao_n1 <- readRDS("2.2b_countryChao.Rda")
}
if(!exists("continentChao_n1")){
  continentChao_n1 <- readRDS("2.2c_continentChao.Rda")
}


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



chaoGlobal_est <- globalChao$Species_table %>% as.data.frame() %>% dplyr::tibble()
chaoGlobal_est <- chaoGlobal_est %>%
  dplyr::bind_cols(
    globalChao$Species_table %>% as.data.frame() %>% rownames() %>% stringr::str_squish()
    ) %>%
  setNames(c("Estimate", "se", "95%Lower", "95%Upper", "Statistic")) %>%
  dplyr::filter(Statistic == "iChao1 (Chiu et al. 2014)")



(global_iNEXTplot <- iNEXT::ggiNEXT(global_iNEXT_out, type=1, facet.var="None", color.var = "Order.q") +
  ggplot2::theme_classic() +
  ggplot2::ggtitle(paste0("Global bee species richness","\n",
                          "n = ", format(global_iNEXT_out$DataInfo$n, big.mark = ","),
                          "; obs. = ", format(global_iNEXT_out$AsyEst$Observed[[1]] %>% round(0),big.mark = ","),
                          # iNEXT
                          "\niNEXT = ", format(global_iNEXT_out$AsyEst$Estimator[[1]] %>% round(0),
                                               big.mark = ","),
                          " (", format(global_iNEXT_out$AsyEst$`95% Lower`[[1]] %>% round(0),
                                       big.mark = ","),
                          "-", format(global_iNEXT_out$AsyEst$`95% Upper`[[1]] %>% round(0),
                                      big.mark = ","),
                          "; +",
                          ((1 - (global_iNEXT_out$AsyEst$Observed[[1]] / global_iNEXT_out$AsyEst$Estimator[[1]]))*100) %>% 
                            round(0),"%)",
                          # i CHAO
                   "\niChao = ", format(chaoGlobal_est$Estimate %>% round(0),
                                        big.mark = ","),
                   " (", format(chaoGlobal_est$`95%Lower` %>% round(0),
                                big.mark = ","),
                   "-", format(chaoGlobal_est$`95%Upper` %>% round(0),
                               big.mark = ","),
                     "; +",
                   ((1 - (global_iNEXT_out$AsyEst$Observed[[1]] / chaoGlobal_est$Estimate))*100) %>% 
                              round(0), "%)"
                   ))+ 
  ggplot2::theme(legend.position="none") +
    # iNEXT
  ggplot2::geom_hline(yintercept = global_iNEXT_out$AsyEst$Estimator[[1]], linetype="solid", color = "black") +
  ggplot2::geom_hline(yintercept = global_iNEXT_out$AsyEst$`95% Lower`[[1]], linetype="dashed", color = "#FD9B63") +
  ggplot2::geom_hline(yintercept = global_iNEXT_out$AsyEst$`95% Upper`[[1]], linetype="dashed", color = "#FD9B63") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = global_iNEXT_out$AsyEst$`95% Upper`[[1]],  
                                    ymax = global_iNEXT_out$AsyEst$`95% Lower`[[1]]),
                       fill = "#FD9B63", alpha = 0.1, colour = NA) +
  # iChao
  ggplot2::geom_hline(yintercept = chaoGlobal_est$Estimate, linetype="solid", color = "black") +
  ggplot2::geom_hline(yintercept = chaoGlobal_est$`95%Lower`, linetype="dashed", color = "#55AD9B") +
  ggplot2::geom_hline(yintercept = chaoGlobal_est$`95%Upper`, linetype="dashed", color = "#55AD9B")+
  ggplot2::geom_ribbon(ggplot2::aes(ymin = chaoGlobal_est$`95%Upper`,  
                                    ymax = chaoGlobal_est$`95%Lower`),
                       fill = "#55AD9B", alpha = 0.1, colour = NA) +
  # Have the plot start at zero on X and Y 
  scale_x_continuous(expand = c(0.001, 0.001)) + scale_y_continuous(expand = c(0.001, 0.001))
)
  
ggplot2::ggsave(global_iNEXTplot, device = "pdf", path = paste0(getwd(), "/Figure_outputs/country_iNEXT"),
                file = "global_iNEXT.pdf", width = 6, height = 5)



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
if(!exists("country_iNEXT")){
  country_iNEXT <- base::readRDS("2.3b_country_iNEXT.Rda")
}

# Get the sample sizes for each country ( get those < 30)
smallSampleSizes_conti <- country_iNEXT$DataInfo %>%
  dplyr::arrange(n) %>%
  dplyr::filter(n < 30)
  # Northern Cyprus is always a problem to remove... probably some weird space character or similar
country_iNEXT$iNextEst$iNextEst <- 
  country_iNEXT$iNextEst$iNextEst[!stringr::str_detect(names(country_iNEXT$iNextEst$iNextEst), "\\sCyprus")]

# Plot all country plots using ggiNEXTwrapper
source("ggiNEXTwrapper.R")
country_summary <- ggiNEXTwrapper(data = country_iNEXT,
               # Remove the countries under a certain sample size
               filterOut = smallSampleSizes_conti %>%
                 dplyr::pull(Assemblage) %>%
                 c(., "Côte d'Ivoire", "Dem. Rep. Korea", "Equatorial Guinea", "Eritrea",
                   "North Macedonia", "Iceland", "Northern Cyprus",
                    "Western Sahara", "West Bank", "Guinea", "Mauritania"),
               legendPerPlot = FALSE,
               nrow = 4,
               ncol = 3,
               iChao_in = countryChao_n1,
               labels = NULL,
               fileName = "iNEXTplots",
               outPath = "/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
               base_width = 8.3,
               base_height = 11.7, 
               dpi = 300)

# Save this data table
readr::write_excel_csv(country_summary, file = paste0("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
                                                      "/2.3b_countryTable.csv"))


###### c. continent records ####
source("iNEXTwrapper.R")
# Run iNEXT for each country using the iNEXT wrapper
continent_iNEXT <- iNEXTwrapper(
  # Should this be continentCounts instead with variableColumn = "continent"?
  data = continentWider, 
                              k = 5,
                              datatype = "abundance",
                              conf = 0.95,
                              mc.cores = 8)

# Save this file
base::saveRDS(continent_iNEXT, 
              file = "2.3c_continent_iNEXT.Rda")


  # Get the sample sizes for each country ( get those < 30)
smallSampleSizes_conti <- continent_iNEXT$DataInfo %>%
  dplyr::arrange(n) %>%
  dplyr::filter(n < 500)

# Plot all country plots using ggiNEXTwrapper
source("ggiNEXTwrapper.R")
continent_summary <- ggiNEXTwrapper(data = continent_iNEXT,
                # Remove the countries under a certain sample size
               filterOut = smallSampleSizes_conti %>%
                 dplyr::pull(Assemblage),
               legendPerPlot = FALSE,
               iChao_in = continentChao_n1,
               nrow = 3,
               ncol = 2,
               labels = NULL,
               fileName = "iNEXTplots_Continent",
               outPath = "/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
               base_width = 8.3,
               base_height = 11.7, 
               dpi = 300)

# Save this data table
readr::write_excel_csv(continent_summary, file = paste0("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
                                                      "/2.3c_continentTable.csv"))

  ###### d. read in if needed ####
if(!exists("global_iNEXT_out")){
  global_iNEXT_out <- readRDS("2.3a_globaliNEXT.Rda")
}
if(!exists("country_iNEXT")){
  country_iNEXT <- base::readRDS("2.3b_country_iNEXT.Rda")
}
if(!exists("continent_iNEXT")){
  continent_iNEXT <- base::readRDS("2.3c_continent_iNEXT.Rda")
}
if(!exists("continent_summary")){
  continent_summary <- base::readRDS("2.3c_continentTable.csv")
}


  ##### 2.4 Iterative sampling ####
source("richnessSampleR.R")
    ###### a. data prep ####
  # Get the difference between the taxonomy and occurrence records
noPointSpecies_diff <- taxonomy_valid %>% 
  dplyr::pull(validName) %>%
  setdiff(beeData$scientificName)

  # Make a list of the data inputs so that this can be run easily on another computer
richnessInputs <- dplyr::lst(
  # The curve created by implementing mosaic::fitModel over the literature samples  
  litTibble,
  checklistFile,
  taxonomyFile,
  taxonomy_valid,
  noPointSpecies_diff,
  beeData_counts,
  continentChecklist,
  countryNOTcontinent,
  worldMap,
  # Get the counts of species per country based on point occurrences
  country_speciesCounts) 
    # Save the file for use on another computer
richnessInputs %>%
  base::saveRDS(., 
                file = "2.4a_richnessInputs.Rda")
if(!exists("richnessInputs")){
  richnessInputs <- base::readRDS("2.4a_richnessInputs.Rda")
}



    # Use the richnessSampleR function to iteratively generate richness estimates
    ###### b. global ####
globalSampled <- parallel::mclapply(
  X = 1:20,
  FUN = richnessSampleR,
  # FUNCTION INPUTS:
  # Input datasets
  richnessInputs = richnessInputs,
  # What scale?
  scale = "Global", # Country, Continent, or Global
  # Functions
  ChaoWrapper = ChaoWrapper,
  iNEXTwrapper = iNEXTwrapper,
  mc.cores = 10
)
  # Save
globalSampled %>%
  base::saveRDS(., file = "2.4b_globalSampled.Rda")

  # Interim read in the separately-calculated parts
globalSampled <- c(base::readRDS("2.4b_globalSampled_p1.Rda"), 
                        base::readRDS("2.4b_globalSampled_p2.Rda"),
                        base::readRDS("2.4b_globalSampled_p3.Rda"))

# Read in if needed
if(!exists("globalSampled")){
  globalSampled <- base::readRDS("2.4b_globalSampled.Rda")
}


counter <- 1
# Extract all of the Chao data
ChaoGlobalSampled <- globalSampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$ChaoOut$Species_table %>%
             as.data.frame() %>%
             tibble::rownames_to_column(var = "Name") %>% 
             dplyr::mutate(Name = stringr::str_squish(Name))
           # Extract the sample size
           sampleSize <- x[[1]]$ChaoOut$Basic_data_information %>% 
             tibble::rownames_to_column(var = "Name") %>% 
             dplyr::mutate(Name = stringr::str_squish(Name)) %>%
             t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% dplyr::tibble() %>%
             dplyr::select(rowname, V1) %>% dplyr::filter(!V1 %in% c("Sample size", "n")) %>%
             dplyr::rename(n = V1)
           # Add the sample size 
           extraction <- extraction %>% 
             dplyr::mutate(n = sampleSize$n)
         }) 
# Combine these together
global_iChaoSamples <- ChaoGlobalSampled %>%
  dplyr::bind_rows() %>%
  # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(Name, "iChao1 \\(Chiu et al\\. 2014\\)")) %>%
  dplyr::mutate(Name = "iChao") %>% 
  dplyr::mutate(variable = "Global", .before = 1)

# Extract all of the iNEXT data
counter <- 1
iNEXTGlobalSampled <- globalSampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$iNEXTout$AsyEst %>%
             dplyr::mutate(n = x[[1]]$iNEXTout$DataInfo$n %>%
                             as.character())
         }) 

# Combine these together
global_iNextSamples <- iNEXTGlobalSampled %>%
  dplyr::bind_rows() %>%
  tibble::rownames_to_column("statistic") %>% 
  # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(statistic, "Species Richness")) %>%
  dplyr::mutate(statistic = "iNEXT") %>% 
  # Rename to match with iChao
  dplyr::rename(Name = statistic,
                Estimate = Estimator,
                `95%Lower` = `95% Lower`,
                `95%Upper` = `95% Upper`) %>%
  dplyr::mutate(variable = "Global", .before = 1)

combined_global_ChaoiNext <- global_iNextSamples %>%
  dplyr::bind_rows(global_iChaoSamples) %>%
  dplyr::arrange( `95%Upper`) %>%
  # Add a column with group numbers, and then make them into groups of ten
  dplyr::mutate(groupNo = ceiling(cur_group_id()/10)) %>%
  dplyr::group_by(groupNo)


# Combine the data to calculate further statistics
# Combine the iNEXT data
globalcombined_iNEXT <- combined_global_ChaoiNext %>%
  dplyr::filter(Name == "iNEXT") %>%
  # Group by country
  dplyr::group_by(variable) %>%
  # Get the medians of the sampled values
  dplyr::mutate(observedRichness = Observed,
                niNEXT = median(n %>% as.numeric(.)),
                iNEXT_est = median(Estimate),
                iNEXT_lower = median(`95%Lower`),
                iNEXT_upper = median(`95%Upper`),
                iNEXT_increasePercent = ((iNEXT_est/observedRichness)-1)*100,
                iNEXT_increase = iNEXT_est - observedRichness,
                level = "Global") %>%
  dplyr::select(variable, observedRichness, level, niNEXT,
                iNEXT_est, iNEXT_lower, iNEXT_upper, iNEXT_increasePercent, iNEXT_increase) %>%
  dplyr::distinct()

# Combine the iChao data
globalcombined_iChao <- combined_global_ChaoiNext %>%
  dplyr::filter(Name == "iChao") %>%
  # Group by country
  dplyr::group_by(variable) %>%
  # Get the medians of the sampled values
  dplyr::mutate(nChao = median(as.numeric(n)),
                iChao_est = median(Estimate),
                iChao_lower = median(`95%Lower`),
                iChao_upper = median(`95%Upper`)) %>%
  dplyr::select(variable, nChao, iChao_est, iChao_lower, iChao_upper) %>%
  dplyr::distinct() 

# Combine them all
combinedStatistics_sampled_global <- globalcombined_iChao %>%
  dplyr::left_join(globalcombined_iNEXT, by = "variable") %>% 
  # Calculate the iChao increases
  dplyr::mutate(iChao_increasePercent = ((iChao_est/observedRichness)-1)*100,
                iChao_increase = iChao_est - observedRichness)





    ###### c. country ####
countrySampled <- parallel::mclapply(
  X = 1:100,
  FUN = richnessSampleR,
  # FUNCTION INPUTS:
  # Input datasets
  richnessInputs = richnessInputs,
  # What scale?
  scale = "Country", # Country, Continent, or Global
  # Functions
  ChaoWrapper = ChaoWrapper,
  iNEXTwrapper = iNEXTwrapper,
  mc.cores = 3
)

# Save
countrySampled %>%
  base::saveRDS(., file = "2.4c_countrySampled.Rda")
# Read in if needed
if(!exists("countrySampled")){
  countrySampled <- base::readRDS("2.4c_countrySampled.Rda")
}


counter <- 1
# Extract all of the Chao data
ChaoCountrySampled <- countrySampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$ChaoOut$diversityTable %>%
             dplyr::mutate(Name = stringr::str_squish(Name))
            # Extract the sample size
           sampleSize <- x[[1]]$ChaoOut$basicTable %>% 
             t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% dplyr::tibble() %>%
             dplyr::select(rowname, V1) %>% dplyr::filter(!V1 %in% c("Sample size", "n")) %>%
             dplyr::rename(n = V1)
            # Add the sample size 
           extraction <- extraction %>% 
             dplyr::left_join(sampleSize, by = c("variable" = "rowname"))
         }) 
# Combine these together
country_iChaoSamples <- ChaoCountrySampled %>%
  dplyr::bind_rows() %>%
  # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(Name, "iChao1 \\(Chiu et al\\. 2014\\)")) %>%
  dplyr::mutate(Name = "iChao") %>% 
  # Remove undesired continents
  dplyr::filter(!stringr::str_detect(variable, "Antarctica|Seven")) %>%
  # Group by continent
  dplyr::group_by(variable)

# Extract all of the iNEXT data
counter <- 1
iNEXTCountrySampled <- countrySampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$iNEXTout$AsyEst 
           extraction <- extraction %>% 
             dplyr::left_join(x[[1]]$iNEXTout$DataInfo %>% 
                                dplyr::select(Assemblage, n) %>%
                                dplyr::mutate(n = as.character(n)),
                              by = c("rNaturalEarth_name" = "Assemblage"))
           
         }) 

# Combine these together
country_iNextSamples <- iNEXTCountrySampled %>%
  dplyr::bind_rows() %>%
  # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(statistic, "Species Richness")) %>%
  dplyr::mutate(statistic = "iNEXT") %>% 
  # Remove undesired continents
  dplyr::filter(!stringr::str_detect(rNaturalEarth_name, "Antarctica|Seven")) %>%
  # Group by rNaturalEarth_name
  dplyr::group_by(rNaturalEarth_name) %>%
  # Rename to match with iChao
  dplyr::rename(variable = rNaturalEarth_name,
                Name = statistic,
                Estimate = Estimator,
                `95%Lower` = `95% Lower`,
                `95%Upper` = `95% Upper`)

combined_count_ChaoiNext <- country_iNextSamples %>%
  dplyr::bind_rows(country_iChaoSamples) %>%
    # Remove low sample countries
  dplyr::filter(!variable %in% c("Côte d'Ivoire", "Dem. Rep. Korea", "Equatorial Guinea", "Eritrea",
                                  "North Macedonia", "Iceland",
                                  "Western Sahara", "West Bank", "Guinea", "Mauritania"
                                 )) %>% 
  dplyr::filter(!stringr::str_detect(variable, "Northern" # for Northern Cyprus
                                     )) %>%
  dplyr::arrange( `95%Upper`) %>%
    # Add a column with group numbers, and then make them into groups of ten
  dplyr::mutate(groupNo = ceiling(cur_group_id()/10)) %>%
  dplyr::group_by(groupNo)

# Make a plot of the iChao values and the confidence intervals
(count_sampledPlot <- ggplot2::ggplot(data = combined_count_ChaoiNext) + 
    ggplot2::geom_violin(position="dodge", alpha=0.5,
                         ggplot2::aes(fill=Name, 
                                      y=`95%Lower`, x=variable), colour =  NA) +
    ggplot2::geom_violin(position="dodge", alpha=0.5,
                         ggplot2::aes(fill=Name,
                                      y=`95%Upper`, x=variable), colour =  NA) +
    ggplot2::geom_violin(position="dodge", alpha=1,
                         ggplot2::aes(fill=Name,
                                      y=Estimate, x=variable), colour =  "black") +
    ggplot2::scale_fill_manual(values=c("#55AD9B", "#FD9B63")) +
    ggplot2::theme_classic() + ggplot2::xlab("Country") + ggplot2::ylab("iChao") +
    ggplot2::facet_wrap(vars(variable), scales = "free", ncol = 5)+ 
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
)
# Save the plot
ggplot2::ggsave(file = "2.4c_countrySampled.pdf", 
                path = paste0(RootPath, "/Figure_outputs"),
                plot = count_sampledPlot, device = "pdf",
                width = 10, height = 50,
                limitsize = FALSE)

  # Combine the data to calculate further statistics
  # Combine the iNEXT data
countrycombined_iNEXT <- combined_count_ChaoiNext %>%
  dplyr::filter(Name == "iNEXT") %>%
    # Group by country
  dplyr::group_by(variable) %>%
    # Get the medians of the sampled values
  dplyr::mutate(observedRichness = Observed,
                niNEXT = median(n %>% as.numeric(.)),
                iNEXT_est = median(Estimate),
                iNEXT_lower = median(`95%Lower`),
                iNEXT_upper = median(`95%Upper`),
                iNEXT_increasePercent = ((iNEXT_est/observedRichness)-1)*100,
                iNEXT_increase = iNEXT_est - observedRichness,
                level = "Country") %>%
  dplyr::select(variable, observedRichness, level, niNEXT,
                iNEXT_est, iNEXT_lower, iNEXT_upper, iNEXT_increasePercent, iNEXT_increase) %>%
  dplyr::distinct()

# Combine the iChao data
countrycombined_iChao <- combined_count_ChaoiNext %>%
  dplyr::filter(Name == "iChao") %>%
  # Group by country
  dplyr::group_by(variable) %>%
  # Get the medians of the sampled values
  dplyr::mutate(nChao = median(as.numeric(n)),
                iChao_est = median(Estimate),
                iChao_lower = median(`95%Lower`),
                iChao_upper = median(`95%Upper`)) %>%
  dplyr::select(variable, nChao, iChao_est, iChao_lower, iChao_upper) %>%
  dplyr::distinct()
  
  # Combine them all
combinedStatistics_sampled_country <- countrycombined_iChao %>%
  dplyr::left_join(countrycombined_iNEXT, by = "variable") %>% 
    # Calculate the iChao increases
  dplyr::mutate(iChao_increasePercent = ((iChao_est/observedRichness)-1)*100,
                iChao_increase = iChao_est - observedRichness)




    ###### d. continent ####
  # Iteratively sample the continent-level richness
continentSampled <- parallel::mclapply(
  X = 1:100,
  FUN = richnessSampleR,
                          # FUNCTION INPUTS:
    # Input datasets
  richnessInputs = richnessInputs,
  # What scale?
  scale = "Continent", # Country, Continent, or Global
  # Functions
  ChaoWrapper = ChaoWrapper,
  iNEXTwrapper = iNEXTwrapper,
  mc.cores = 10
)

# Save
continentSampled %>%
  base::saveRDS(., file = "2.4d_continentSampled.Rda")
# Read in if needed
if(!exists("continentSampled")){
  continentSampled <- base::readRDS("2.4d_continentSampled.Rda")
}

counter <- 1
# Extract all of the Chao data
ChaoContinentSampled <- continentSampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$ChaoOut$diversityTable
           # Extract the sample size
           sampleSize <- x[[1]]$ChaoOut$basicTable %>% 
             t() %>% as.data.frame() %>% tibble::rownames_to_column() %>% dplyr::tibble() %>%
             dplyr::select(rowname, V1) %>% dplyr::filter(!V1 %in% c("Sample size", "n")) %>%
             dplyr::rename(n = V1)
           # Add the sample size 
           extraction <- extraction %>% 
             dplyr::left_join(sampleSize, by = c("variable" = "rowname"))
         }) 
  # Combine these together
continent_iChaoSamples <- ChaoContinentSampled %>%
  dplyr::bind_rows() %>%
    # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(Name, "iChao1 \\(Chiu et al\\. 2014\\)")) %>%
  dplyr::mutate(Name = "iChao") %>% 
    # Remove undesired continents
  dplyr::filter(!stringr::str_detect(variable, "Antarctica|Seven")) %>%
    # Group by continent
  dplyr::group_by(variable)

# Extract all of the iNEXT data
counter <- 1
iNEXTContinentSampled <- continentSampled %>%
  lapply(X = .,
         function(x){
           counter <<- counter + 1
           extraction <- x[[1]]$iNEXTout$AsyEst 
           extraction <- extraction %>% 
             dplyr::left_join(x[[1]]$iNEXTout$DataInfo %>% 
                                dplyr::select(Assemblage, n) %>%
                                dplyr::mutate(n = as.character(n)),
                              by = c("continent" = "Assemblage"))
             
         }) 

# Combine these together
continent_iNextSamples <- iNEXTContinentSampled %>%
  dplyr::bind_rows() %>%
  # Select the desired statistic estimate
  dplyr::filter(stringr::str_detect(statistic, "Species Richness")) %>%
  dplyr::mutate(statistic = "iNEXT") %>% 
  # Remove undesired continents
  dplyr::filter(!stringr::str_detect(continent, "Antarctica|Seven")) %>%
  # Group by continent
  dplyr::group_by(continent) %>%
    # Rename to match with iChao
  dplyr::rename(variable = continent,
                Name = statistic,
                Estimate = Estimator,
                `95%Lower` = `95% Lower`,
                `95%Upper` = `95% Upper`)

combined_cont_ChaoiNext <- continent_iChaoSamples %>%
  dplyr::bind_rows(continent_iNextSamples)


# Combine the data to calculate further statistics
# Combine the iNEXT data
contCombined_iNEXT <- combined_cont_ChaoiNext %>%
  dplyr::filter(Name == "iNEXT") %>%
  # Group by country
  dplyr::group_by(variable) %>%
  # Get the medians of the sampled values
  dplyr::mutate(observedRichness = Observed,
                niNEXT = median(n %>% as.numeric(.)),
                iNEXT_est = median(Estimate),
                iNEXT_lower = median(`95%Lower`),
                iNEXT_upper = median(`95%Upper`),
                iNEXT_increasePercent = ((iNEXT_est/observedRichness)-1)*100,
                iNEXT_increase = iNEXT_est - observedRichness,
                level = "Continent") %>%
  dplyr::select(variable, observedRichness, level, niNEXT,
                iNEXT_est, iNEXT_lower, iNEXT_upper, iNEXT_increasePercent, iNEXT_increase) %>%
  dplyr::distinct()

# Combine the iChao data
contCombined_iChao <- combined_cont_ChaoiNext %>%
  dplyr::filter(Name == "iChao") %>%
  # Group by country
  dplyr::group_by(variable) %>%
  # Get the medians of the sampled values
  dplyr::mutate(nChao = median(as.numeric(n)),
                iChao_est = median(Estimate),
                iChao_lower = median(`95%Lower`),
                iChao_upper = median(`95%Upper`)) %>%
  dplyr::select(variable, nChao, iChao_est, iChao_lower, iChao_upper) %>%
  dplyr::distinct()

# Combine them all
combinedStatistics_sampled_continent <- contCombined_iChao %>%
  dplyr::left_join(contCombined_iNEXT, by = "variable") %>% 
  # Calculate the iChao increases
  dplyr::mutate(iChao_increasePercent = ((iChao_est/observedRichness)-1)*100,
                iChao_increase = iChao_est - observedRichness)


# Make a plot of the iChao values and the confidence intervals
(cont_sampledPlot <- ggplot2::ggplot(data = combined_cont_ChaoiNext) + 
  ggplot2::geom_violin(position="dodge", alpha=0.5,
                       ggplot2::aes(fill=Name, 
                         y=`95%Lower`, x=variable)#,
                       #fill = "blue"
                       , colour =  NA) +
  ggplot2::geom_violin(position="dodge", alpha=0.5,
                       ggplot2::aes(fill=Name,
                         y=`95%Upper`, x=variable)#,
                       #fill = "red"
                       , colour =  NA) +
  ggplot2::geom_violin(position="dodge", alpha=1,
                       ggplot2::aes(fill=Name,
                         y=Estimate, x=variable),
                       #fill = "purple"
                       colour =  "black") +
    ggplot2::scale_fill_manual(values=c("#55AD9B", "#FD9B63")) +
  ggplot2::theme_classic() + ggplot2::xlab("Continent") + ggplot2::ylab("Species estimate") +
    guides(fill=guide_legend(title="Statistic"))) 
# Save the plot
ggplot2::ggsave(file = "2.4d_continentSampled.pdf", 
                path = paste0(RootPath, "/Figure_outputs"),
                plot = cont_sampledPlot, device = "pdf",
                width = 7, height = 5)

  ##### 2.5 Violin plots ####
# Build a legend
(violinLegend <- ggplot2::ggplot(dplyr::tibble(name = c("yes","yes"),
                                               statistic = c("iChao", "iNEXT") %>%
                                                 factor(levels = c("iChao", "iNEXT")),
                                               est = c(1,1)) ,
                                 aes(x = name, y = est)) + 
   ggplot2::geom_point(data = dplyr::tibble(name = c("yes","yes"),
                                            statistic = c("iChao", "iNEXT") %>%
                                              factor(levels = c("iChao", "iNEXT")),
                                            est = c(1,1)) ,
                       ggplot2::aes(y=est, x = name, colour = "red")) + 
   scale_color_manual(labels = c('Observed'), 
                      values = c('grey30'))  +
   ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT"),
                              values = c("iChao" = "#55AD9B", "iNEXT" = "#FD9B63")) +
   ggplot2::theme_classic() +
   theme(legend.title = element_blank(),
         legend.position = 'right',
         legend.margin = margin(0, 0, 0, 0),
         legend.spacing.y = unit(0, "pt")) +
   ggplot2::guides(fill = ggplot2::guide_legend(ncol = 1, byrow = TRUE,
                                                reverse = TRUE,
                                                 order = 1))
)


    ###### a. country ####
  # Find the top and bottom 10 countries for median richness
Top10_countries <- combined_count_ChaoiNext %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(median = median(Estimate)) %>%
  dplyr::arrange(median)  %>%
  dplyr::ungroup() %>% 
  # Take the top half of the data
  dplyr::slice((nrow(.)-9):nrow(.))

Top10_data <- combined_count_ChaoiNext %>%
  dplyr::arrange(.by_group = TRUE, Estimate, variable) %>%
  dplyr::mutate(variable = factor(variable, 
                                  levels = rev(unique(Top10_countries$variable)))) %>% 
  # Filter for only the top and bottom ten countries
  dplyr::filter(variable %in% Top10_countries$variable)  %>%
  dplyr::mutate(level = "Country")

# plot the top half of countries 
(country_top10_violin <- ggplot2::ggplot(Top10_data,
                                              aes(x = name, y = est)) + 
   ggplot2::geom_violin(position="dodge", alpha=0.5,
                        ggplot2::aes(fill=Name, 
                                     y=`95%Lower`, x=variable)#,
                        #fill = "blue"
                        , colour =  NA) +
   ggplot2::geom_violin(position="dodge", alpha=0.5,
                        ggplot2::aes(fill=Name,
                                     y=`95%Upper`, x=variable)#,
                        #fill = "red"
                        , colour =  NA) +
   ggplot2::geom_violin(position="dodge", alpha=1,
                        ggplot2::aes(fill=Name,
                                     y=Estimate, x=variable),
                        #fill = "purple"
                        colour =  "black") +
   ggplot2::scale_fill_manual(values=c("#55AD9B", "#FD9B63")) +
    ggplot2::geom_point(data = Top10_data %>%
                          dplyr::distinct(variable, Observed), 
                        ggplot2::aes(x = variable, y = Observed), col = "grey40") + 
   ggplot2::theme_classic() + ggplot2::xlab("Continent") + ggplot2::ylab("Species estimate") +
   guides(fill=guide_legend(title="Statistic")) +
   ggplot2::theme_classic() +
   ggplot2::theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
   #ggplot2::ylim(c(0, 4500)) +  
  ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species")) +
    annotation_custom(cowplot::ggdraw(cowplot::get_legend(violinLegend)) %>%
                        ggplot2::ggplotGrob(), xmin = 10, xmax = 8, 
                      ymin = 3500, ymax = 4900)
)


    ###### b. continent #####
  # Get the order of largest continents
continent_order <- combined_cont_ChaoiNext %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(median = median(Estimate)) %>%
  dplyr::arrange(median)  %>%
  dplyr::ungroup()

  # Build a continent dataset to plot
continentViolinData <- combined_cont_ChaoiNext %>%
  dplyr::arrange(.by_group = FALSE, Observed) %>%
  dplyr::mutate(variable = factor(variable, 
                                  levels = rev(unique(continent_order$variable))) ) %>%
  dplyr::mutate(level = "Continent")

(continent_violin <- ggplot2::ggplot(continentViolinData,
                                         aes(x = name, y = Estimate)) + 
   ggplot2::geom_violin(position="dodge", alpha=0.5,
                        ggplot2::aes(fill=Name, 
                                     y=`95%Lower`, x=variable)#,
                        #fill = "blue"
                        , colour =  NA) +
   ggplot2::geom_violin(position="dodge", alpha=0.5,
                        ggplot2::aes(fill=Name,
                                     y=`95%Upper`, x=variable)#,
                        #fill = "red"
                        , colour =  NA) +
   ggplot2::geom_violin(position="dodge", alpha=1,
                        ggplot2::aes(fill=Name,
                                     y=Estimate, x=variable),
                        #fill = "purple"
                        colour =  "black") +
   ggplot2::scale_fill_manual(values=c("#55AD9B", "#FD9B63")) +
   ggplot2::geom_point(data = continentViolinData %>%
                         dplyr::distinct(variable, Observed), 
                       ggplot2::aes(x = variable, y = Observed), col = "grey40") + 
   ggplot2::theme_classic() + ggplot2::xlab("Continent") + ggplot2::ylab("Species estimate") +
   guides(fill=guide_legend(title="Statistic")) +
   ggplot2::theme_classic() +
   ggplot2::theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
   #ggplot2::ylim(c(0, 9500)) +  
   ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species")) #+
 #annotation_custom(cowplot::ggdraw(cowplot::get_legend(barLegend)) %>%
 #                    ggplot2::ggplotGrob(), xmin = 17, xmax = 20, 
 #                  ymin = 2000, ymax = 4500)
)

    ###### c. global ####
# Build a continent dataset to plot
GlobalViolinData <- combined_global_ChaoiNext %>%
  dplyr::arrange(Observed) %>%
  dplyr::mutate(variable = factor(variable,
                                  levels = "Global")) %>%
  dplyr::mutate(level = "Global")

(global_violin <- ggplot2::ggplot(GlobalViolinData,
                                     aes(x = name, y = Estimate)) + 
    ggplot2::geom_violin(position="dodge", alpha=0.5,
                         ggplot2::aes(fill=Name, 
                                      y=`95%Lower`, x=variable)#,
                         #fill = "blue"
                         , colour =  NA) +
    ggplot2::geom_violin(position="dodge", alpha=0.5,
                         ggplot2::aes(fill=Name,
                                      y=`95%Upper`, x=variable)#,
                         #fill = "red"
                         , colour =  NA) +
    ggplot2::geom_violin(position="dodge", alpha=1,
                         ggplot2::aes(fill=Name,
                                      y=Estimate, x=variable),
                         #fill = "purple"
                         colour =  "black") +
    ggplot2::scale_fill_manual(values=c("#55AD9B", "#FD9B63")) +
    ggplot2::geom_point(data = GlobalViolinData %>%
                          dplyr::distinct(variable, Observed), 
                        ggplot2::aes(x = variable, y = Observed), col = "grey40") + 
    ggplot2::theme_classic() + ggplot2::xlab("Continent") + ggplot2::ylab("Species estimate") +
    guides(fill=guide_legend(title="Statistic")) +
    ggplot2::theme_classic() +
    #ggplot2::facet_grid(cols = vars(level), scales = "free", space = "free") +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    #ggplot2::ylim(c(0, 9500)) +  
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species")) #+
  #annotation_custom(cowplot::ggdraw(cowplot::get_legend(barLegend)) %>%
  #                    ggplot2::ggplotGrob(), xmin = 17, xmax = 20, 
  #                  ymin = 2000, ymax = 4500)
)

    ###### d. combine ####
# Combine the plots 
(violinPlots_iter <- 
   cowplot::plot_grid(
   cowplot::plot_grid(global_violin,
                                     continent_violin + ggplot2::ylab(""),
                                     country_top10_violin + ggplot2::ylab(""),
                                     cols = 3,
                                     labels = c("(A)","(B)", "(C)"),
                                     rel_widths = c(2, 7, 10), align = "h"),
  
                     TaxoOccPlot +
                       ggplot2::theme(legend.position.inside = c(0.15,0.7))+
                       ggplot2::ylim(0,25000),
                     labels = c("", "(D)"),
                     nrow = 2, ncol = 1, align = 'none', rel_heights = c(1.5, 1)
                     )
)



# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","2.5dII_violinIterativePlots_Taxo.pdf"),
                   plot = violinPlots_iter,
                   base_width = 10,
                   base_height = 7)

 #### 3.0 Downstream statistics ####
  ##### 3.1 Combine metrics ####
    ###### a. Combine data ####
  # Combine the iterative sampled data above
combinedStatistics <- dplyr::bind_rows(combinedStatistics_sampled_country %>%
                                         dplyr::mutate(scale = "Country", .after = 1),
                                       combinedStatistics_sampled_continent %>%
                                         dplyr::mutate(scale = "Continent", .after = 1),
                                       combinedStatistics_sampled_global %>%
                                         dplyr::mutate(scale = "Global", .after = 1)
                                       ) %>%
  dplyr::rename(name = variable)
# Save this file
readr::write_excel_csv(combinedStatistics, file = "Table_outputs/3.1a_combinedStatistics.csv")
if(!exists("combinedStatistics")){
  combinedStatistics <- readr::read_csv("Table_outputs/3.1a_combinedStatistics.csv")
}

  # Thge below code is for a SINGLE RUN; not the iterative sampled version above
# Make a combined global datasheet to match the rest
      #     global_combined <- dplyr::tibble(
      #       # Add iNEXT
      #       name = "Global",
      #       level = "Global",
      #       iNEXT_est = dplyr::pull(global_iNEXT_out$AsyEst[1,], Estimator),
      #       iNEXT_lower = dplyr::pull(global_iNEXT_out$AsyEst[1,] , `95% Lower`),
      #       iNEXT_upper = dplyr::pull(global_iNEXT_out$AsyEst[1,] , `95% Upper`)) %>%
      #       # Add iChao
      #       dplyr::mutate(
      #         n = globalChao$Basic_data_information %>% dplyr::filter(Variable == "n") %>% dplyr::pull(Value) %>%
      #           as.numeric(),
      #         observedRichness = globalChao$Basic_data_information %>% 
      #           dplyr::filter(Variable == "D") %>% dplyr::pull(Value) %>% as.numeric(),
      #         iChao_est = globalChao_iChao1 %>% dplyr::pull(Estimate),
      #         iChao_lower = globalChao_iChao1 %>% dplyr::pull(`95%Lower`),
      #         iChao_upper = globalChao_iChao1 %>% dplyr::pull(`95%Upper`),
      #         iNEXT_increasePercent = ((iNEXT_est/observedRichness)-1)*100,
      #         iChao_increasePercent = ((iChao_est/observedRichness)-1)*100 )
      #     
      #     
      #     # Combine all levels of metrics into one table
      #     combinedStatistics <- continent_summary %>%
      #       dplyr::rename(name = level) %>%
      #       dplyr::mutate(level = "Continental") %>%
      #       dplyr::bind_rows(country_summary %>% dplyr::rename(name = level) %>% 
      #                          dplyr::mutate(level = "Country")) %>%
      #       dplyr::bind_rows(global_combined) %>%
      #         # Turn level into a factor
      #       dplyr::mutate(level = level %>% factor(x = ., levels = c("Global", "Continental", "Country"),
      #                                              ordered = TRUE)) %>%
      #       dplyr::group_by(level) %>%
      #       dplyr::mutate(iNEXT_increase = iNEXT_est - observedRichness,
      #                     iChao_increase = iChao_est - observedRichness) 
      #     



    ###### b. pull apart estimates ####
  # Make an iNEXT data frame
combined_iNEXT <- combinedStatistics %>%
  dplyr::select(name, level,
                niNEXT,observedRichness, iNEXT_est, iNEXT_lower, iNEXT_upper, iNEXT_increasePercent) %>%
  dplyr::rename(est = iNEXT_est, lower = iNEXT_lower, upper = iNEXT_upper, 
                increasePercentage = iNEXT_increasePercent) %>%
  dplyr::mutate(statistic = "iNEXT")

# Make an iChao data frame
combined_iChao <- combinedStatistics %>%
  dplyr::select(name, level, nChao,observedRichness, iChao_est, iChao_lower, iChao_upper, iChao_increasePercent) %>%
  dplyr::rename(est = iChao_est, lower = iChao_lower, upper = iChao_upper, 
                increasePercentage = iChao_increasePercent)%>%
  dplyr::mutate(statistic = "iChao")

  ###### c. combine longer ####
# Combine these longer
longerCombined <- dplyr::bind_rows(combined_iNEXT, combined_iChao) %>%
    # Shorten some country names
  dplyr::mutate(name = dplyr::if_else(name == "United States of America","USA",name),
                name = dplyr::if_else(name == "Russian Federation","Russia",name),
                name = dplyr::if_else(name == "Democratic Republic of the Congo","DRC",name),
                name = dplyr::if_else(name == "Lao People's Democratic Republic","Lao",name),
                name = dplyr::if_else(name == "United Arab Emirates","UAE",name),
                name = dplyr::if_else(name == "Central African Republic","CAR",name),
                name = dplyr::if_else(name == "Saint Vincent and the Grenadines","SV & Grenadines",name),
                name = dplyr::if_else(name == "Republic of Cabo Verde","Cabo Verde",name),
                name = dplyr::if_else(name == "Northern Mariana Islands","N. Mariana",name),
                name = dplyr::if_else(name == "São Tomé and Principe","São Tomé",name),
                name = dplyr::if_else(name == "Bosnia and Herzegovina","Bosnia",name),
                name = dplyr::if_else(name == "Republic of the Congo","Rep. Congo",name))
  # save
readr::write_excel_csv(combinedStatistics, file = "Table_outputs/3.1c_combinedStatistics_longer.csv")

    ###### d. test CI overlap ####
CIoverlap <- combinedStatistics %>%
  dplyr::select(iChao_lower, iChao_upper, iNEXT_lower, iNEXT_upper, level) %>% 
  dplyr::mutate(TEST = dplyr::if_else(iChao_lower <= iNEXT_upper && iNEXT_lower <= iChao_upper,
                                      "yes", "no"))

    ###### e. read in if needed ####
if(!exists("country_summary")){
  country_summary <- readr::read_csv(paste0("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
                                            "/2.3b_countryTable.csv"))
}
if(!exists("continent_summary")){
  continent_summary <- readr::read_csv(paste0("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
                                              "/2.3c_continentTable.csv"))
}
if(!exists("combinedStatistics")){
  combinedStatistics <- readr::read_csv(paste0("Table_outputs/3.1c_combinedStatistics_longer.csv"))
}

  ##### 3.2 Bar plots ####
    ###### a. Global ####
(globalBoxplot <- ggplot2::ggplot(longerCombined %>% dplyr::filter(level == "Global"),
                           aes(x = name, y = est)) + 
  ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                      position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(values = c("#FD9B63", "#55AD9B") %>% rev()) +
    # Add in the observed richness
   ggplot2::geom_bar(aes(fill = statistic, y = observedRichness),
                     position = position_dodge(0.90),  stat = "identity", fill = "grey", 
                     alpha = 0.5) +
    # Add error bars 
  ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                          position = position_dodge(0.90), width = 0.2, linewidth = 1, 
                         col = c("black", "black")) +

   ggplot2::theme_classic() +
   ggplot2::theme(legend.position = "none") +
   #ggplot2::ylim(c(20000, 27000)) +
  ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species"))
  )

###### b. Continent ####
(continentBoxplot <- ggplot2::ggplot(longerCombined %>% dplyr::filter(level == "Continent") %>%
                                       dplyr::arrange(.by_group = FALSE, observedRichness) %>%
                                       dplyr::mutate(name = factor(name, 
                                                                   levels = rev(unique(.$name)))),
                                     aes(x = name, y = est)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(values = c("#FD9B63", "#55AD9B") %>% rev()) +
   # Add in the observed richness
   ggplot2::geom_bar(aes(fill = statistic, y = observedRichness),
                     position = position_dodge(0.90),  stat = "identity", fill = "grey", 
                     alpha = 0.5) +
   # Add error bars 
   ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                          position = position_dodge(0.90), width = 0.2, linewidth = 1#, 
                          #col = c("black", "black")
                          ) +
   
   ggplot2::theme_classic() +
   ggplot2::theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
   #ggplot2::ylim(c(20000, 27000)) +
   ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species"))
)



###### c. Country ####

# Make a custom legend
(barLegend <- ggplot2::ggplot(dplyr::tibble(name = c("yes","yes","yes"),
                                            statistic = c("iChao", "iNEXT", "Observed") %>%
                                              factor(levels = c("iChao", "iNEXT", "Observed")),
                                            est = c(1,1,1)) ,
                              aes(x = name, y = est)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT", "Observed"),
                              values = c("iChao" = "#55AD9B", "iNEXT" = "#FD9B63",
                                         "Observed" = "grey")) 
)



  # plot the top half of countries 
(countryBoxplot_TOP <- ggplot2::ggplot(longerCombined %>% dplyr::filter(level == "Country") %>%
                                         dplyr::ungroup() %>% 
                                       dplyr::arrange(.by_group = FALSE, observedRichness, name) %>%
                                       dplyr::mutate(name = factor(name, 
                                                                   levels = rev(unique(.$name)))) %>%
                                      # Take the top half of the data
                                     dplyr::slice_tail(n= (nrow(.)/2)),
                                     aes(x = name, y = est)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT"),
                              values = c("iNEXT" = "#FD9B63", "iChao" = "#55AD9B")) +
   # Add in the observed richness
   ggplot2::geom_bar(aes(fill = statistic, y = observedRichness),
                     position = position_dodge(0.90),  stat = "identity", fill = "grey", 
                     alpha = 0.5) +
   # Add error bars 
   ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                          position = position_dodge(0.90), width = 0.2, linewidth = 1#, 
                          #col = c("black", "black")
                          ) +
   ggplot2::theme_classic() +
   ggplot2::theme(legend.position = "none",
                  axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
   ggplot2::ylim(c(0, 4500)) +   ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species"))+
    annotation_custom(cowplot::ggdraw(cowplot::get_legend(barLegend)) %>%
                        ggplot2::ggplotGrob(), xmin = 170, xmax = 3, 
                      ymin = 2000, ymax = 4500)
)

  # plot the bottom half
(countryBoxplot_Bottom <- ggplot2::ggplot(longerCombined %>% dplyr::filter(level == "Country") %>%
                                            dplyr::ungroup() %>% 
                                            dplyr::arrange(.by_group = FALSE, observedRichness, name) %>%
                                            dplyr::mutate(name = factor(name, 
                                                                        levels = rev(unique(.$name)))) %>%
                                            # Take the bottom half of the data
                                            dplyr::slice_head(n= (nrow(.)/2)),
                                       aes(x = name, y = est)) + 
    ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                      position = position_dodge(0.90),  stat = "identity") +
    ggplot2::scale_fill_manual(values = c("#FD9B63", "#55AD9B") %>% rev()) +
    # Add in the observed richness
    ggplot2::geom_bar(aes(fill = statistic, y = observedRichness),
                      position = position_dodge(0.90),  stat = "identity", fill = "grey", 
                      alpha = 0.5) +
    # Add error bars 
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                           position = position_dodge(0.90), width = 0.2, linewidth = 1#, 
                           #col = c("black", "black")
                           ) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    ggplot2::ylim(c(0, 4500)) +
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species")) 
)

# plot the top half of countries 
(countryBoxplot_TOP_BOTTOM <- cowplot::plot_grid(countryBoxplot_TOP,
                                                 countryBoxplot_Bottom +
                                                   #ggplot2::theme(
                                                   #  axis.text.x = element_text(angle = 60, vjust = 1, hjust=0.5)) +
                                                   ggplot2::ylab(c("")) ,
                                                 nrow = 2,
                                                 labels = c("",""))
)


  ###### d. output ####
  # MAKE plots with a subset 
# Combine the plots 
(barPlots <- cowplot::plot_grid(
  cowplot::plot_grid(globalBoxplot,
                     continentBoxplot +
                       ggplot2::theme( axis.text.x = element_text(angle = 0, vjust = 1, hjust=0.5)) +
                       ggplot2::ylab(c("")) ,
                     labels = c("(A)","(B)"),
                     rel_widths = c(0.3, 1), align = "h"),
   countryBoxplot_TOP_BOTTOM, 
  TaxoOccPlot +
    ggplot2::theme(legend.position.inside = c(0.15,0.7))+
    ggplot2::ylim(0,32000),
  labels = c("", "(C)", "(D)"),
   nrow = 3, ncol = 1, align = 'none', rel_widths = c(0.5, 0.5))
)

# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","3.2d_BarPlots.pdf"),
                   plot = barPlots,
                   base_width = 8,
                   base_height = 9)


  # MAKE the supp plot
# Combine the plots 
(barPlots_supp <- cowplot::plot_grid(globalBoxplot,
                                     countryBoxplot_TOP, continentBoxplot, countryBoxplot_Bottom,
                                     labels = c("(A)","(C)", "(B)", "(D)"),
                                     ncol = 2, align = 'v', axis = 'l', 
                                     rel_widths = c(0.5, 2,0.5,2))
)
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","3.2d_suppBarPlots.pdf"),
                   plot = barPlots_supp,
                   base_width = 15,
                   base_height = 10)



#### 4.0 Country maps ####
  ##### 4.1 Map functions ####
countryMapFunction <- function(
    map_in = NULL,
    mapColumn = "observedRichness",
    class_n = 10,
    class_Style = "fisher",
    legendTitle = "Class count",
    title = "your title here",
    naColour = "grey"){
  
  # Download a world map to convert countries to continents
  worldMap <- rnaturalearth::ne_countries(returnclass = "sf",
                                          scale = 50, type = "countries") 

    # Make class intervals.
# Class intervals from ?classIntervals: fixed", "sd", "equal", "pretty", "quantile", 
# "kmeans", "hclust", "bclust", "fisher", "jenks", "dpih" or "headtails"
classes <- classInt::classIntervals(map_in[[mapColumn]], n = class_n, 
                                    style = class_Style, dig.lab=20,
                                    dataPrecision=0)
# Next we'll create a new column in our sf object using the base R cut() function to cut up our 
# percent variable into distinct groups:
map_in <- map_in %>%
  dplyr::mutate(class_count = cut(map_in[[mapColumn]], 
                                  classes$brks,
                                  include.lowest = T, dig.lab = 10)) %>%
  # format the class_count column to remove spaces, add comma break, and join min and max
  dplyr::mutate(class_count2 = class_count %>%
                  stringr::str_remove("\\[|\\]|\\(|\\)") %>%
                  stringr::str_remove_all("\\.[0-9]+") %>%
                  stringr::str_remove("\\]") %>%
                  stringr::str_replace(",", "-")) %>%
  tidyr::separate(col = class_count2, into = c("min", "max"), sep = "-") %>%
  dplyr::mutate(min = min %>% as.numeric() %>% format(big.mark = ",") %>% 
                  stringr::str_remove("\\s+"),
                max = max %>% as.numeric() %>% format(big.mark = ",") %>% 
                  stringr::str_remove("\\s+"),
                class_count2 = stringr::str_c(min, max, sep = "-") ) 

# Join the map and occurrence data
fullMap <-  dplyr::full_join(worldMap,  map_in %>% sf::st_drop_geometry(),
                             by = c("name_long" = "name_long")) %>%
  # Remove na rows
  tidyr::drop_na(tidyselect::any_of(mapColumn))

# Make the map
(spCountryMap <- ggplot2::ggplot(data = fullMap ) +
   # Add in a blank base-map to highlight countries with no data
   ggplot2::geom_sf(data = worldMap, size = 0.15, fill = naColour)+ 
   # Plot and colour the terrestrial base map
   ggplot2::geom_sf(ggplot2::aes(fill = class_count), size = 0.15)+ 
   # Set map limits, if wanted
   ggplot2::coord_sf(expand = FALSE, ylim = c(-60,90), lims_method = "box",
                     clip = "on") + 
   # Map formatting
   # Add in the map's north arrow
    # ggspatial::annotation_north_arrow(location = "tl", which_north = "true", 
    #                                   pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm"), 
    #                                   style = ggspatial::north_arrow_fancy_orienteering()
    #                                   ) + # Add in NORTH ARROW
   ggplot2::theme(panel.grid.major = ggplot2::element_line(color = grDevices::gray(.1, alpha = 0.1), 
                                                           linetype = "dashed", linewidth = 0.5), # Add grid lines
                  panel.border = ggplot2::element_rect(color = grDevices::gray(.1, alpha = 1), 
                                                       linetype = "solid", linewidth = 0.5,
                                                       fill = NA), # add panel border
                  # Add background - colour in the ocean
                  panel.background = ggplot2::element_rect(fill = "aliceblue") )+ 
   # For Dorey colour scheme use the below
   ggplot2::scale_fill_viridis_d(option = "inferno",
                                 na.value = naColour,
                                 name = legendTitle,
                                 labels = fullMap %>%
                                   dplyr::arrange( dplyr::across(dplyr::matches(mapColumn) )) %>% 
                                   dplyr::distinct(class_count2) %>%
                                   # options = "magma", "inferno", "plasma", "cividis"
                                   dplyr::pull(class_count2)) + 
   # Add in X and Y labels
   ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
   # Add in the title
   ggplot2::ggtitle( title)  )

return(spCountryMap)
} # END countryMapFunction



continentMapFunction <- function(
    map_in = NULL,
    mapColumn = "observedRichness",
    legendTitle = "Species",
    title = "your title here",
    naColour = "grey"){
  
  map_in$plotVar <- map_in[[mapColumn]]
  
  # Make the map
  (spCountryMap <- ggplot2::ggplot(data = map_in ) +
      # Add in a blank base-map to highlight countries with no data
      ggplot2::geom_sf(data = map_in, size = 0.15, fill = naColour) + 
      # Plot and colour the terrestrial base map
      ggplot2::geom_sf(ggplot2::aes(fill = plotVar), size = 0.15) + 
      # Set map limits, if wanted
      ggplot2::coord_sf(expand = FALSE, ylim = c(-60,90), lims_method = "box") + 
      # Map formatting
      # Add in the map's north arrow
      # ggspatial::annotation_north_arrow(location = "tl", which_north = "true", 
      #                                   pad_x = unit(0.1, "cm"), pad_y = unit(0.1, "cm"), 
      #                                   style = ggspatial::north_arrow_fancy_orienteering()
      #                                   ) + # Add in NORTH ARROW
      ggplot2::theme(panel.grid.major = ggplot2::element_line(color = grDevices::gray(.1, alpha = 0.1), 
                                                              linetype = "dashed", linewidth = 0.5), # Add grid lines
                     panel.border = ggplot2::element_rect(color = grDevices::gray(.1, alpha = 1), 
                                                          linetype = "solid", linewidth = 0.5,
                                                          fill = NA), # add panel border
                     # Add background - colour in the ocean
                     panel.background = ggplot2::element_rect(fill = "aliceblue") )+ 
      # For Dorey colour scheme use the below
      ggplot2::scale_fill_viridis_c(option = "inferno",
                                    na.value = naColour,
                                    name = legendTitle) + 
      # Add in X and Y labels
      ggplot2::xlab("Longitude") + ggplot2::ylab("Latitude") + 
      # Add in the title
      ggplot2::ggtitle( title)  )
  
  return(spCountryMap)
} # END continentMapFunction


  ##### 4.2 Add statistics to map ####
    # For the sake of matching the statistics to the countries, make the names longer again
combinedStatistics_map <- combinedStatistics %>%
  # LENGTHEN some country names
  dplyr::mutate(name = dplyr::if_else(name == "USA",             "United States of America",         name),
                name = dplyr::if_else(name == "Russia",          "Russian Federation",               name),
                name = dplyr::if_else(name == "DRC",             "Democratic Republic of the Congo", name),
                name = dplyr::if_else(name == "Lao",             "Lao People's Democratic Republic", name),
                name = dplyr::if_else(name ==  "UAE",            "United Arab Emirates",             name),
                name = dplyr::if_else(name == "CAR",             "Central African Republic",         name),
                name = dplyr::if_else(name == "SV & Grenadines", "Saint Vincent and the Grenadines", name),
                name = dplyr::if_else(name == "Cabo Verde",      "Republic of Cabo Verde",           name),
                name = dplyr::if_else(name == "N. Mariana",      "Northern Mariana Islands",         name),
                name = dplyr::if_else(name == "São Tomé",        "São Tomé and Principe",            name),
                name = dplyr::if_else(name == "Bosnia",          "Bosnia and Herzegovina",           name),
                name = dplyr::if_else(name == "Rep. Congo",      "Republic of the Congo",            name))

    ###### a. country ####
  # Add the statistic data to the map 
    # Match by name_long andthen by name
match1 <- worldMap %>%
  dplyr::select(name, name_long) %>%
  dplyr::full_join(combinedStatistics_map,
                   by = c("name_long" = "name")) %>%
  tidyr::drop_na(observedRichness)

match2 <- worldMap %>%
  dplyr::select(name, name_long) %>%
  dplyr::full_join(combinedStatistics_map,
                   by = c("name" = "name")) %>%
  tidyr::drop_na(observedRichness)

  # Combine these two match methods to get all countries included
dataMap <- match1 %>%
  dplyr::filter(!name %in% match2$name) %>%
  dplyr::bind_rows(match2)

    ###### b. continent ####
sf::sf_use_s2(FALSE)

continentDataMap <- rnaturalearth::ne_countries(returnclass = "sf",
                                                scale = 50, type = "map_units")  %>%
  sf::st_make_valid() %>% 
  dplyr::group_by(continent) %>% 
  dplyr::summarise(across(geometry, ~ sf::st_union(.)), .groups = "keep") %>%
  dplyr::summarise(across(geometry, ~ sf::st_combine(.))) %>%
  dplyr::left_join(combinedStatistics,
                   by = c("continent" = "name")) 

##### 4.3 Country maps ####
    ###### a. observed richness ####
(obsRich_map <- countryMapFunction(
    map_in = dataMap %>% 
      dplyr::filter(level == "Country"),
    mapColumn = "observedRichness",
    class_n = 10, naColour = "white",
    class_Style = "fisher",
    legendTitle = "Species class",
    title = "Observed species richness"))

    ###### b. iChao richness ####
(iChaoRich_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iChao_est",
  class_n = 10, naColour = "white",
  legendTitle = "Species class",
  class_Style = "fisher",
  title = "iChao species richness"))

    ###### c. iNEXT richness ####
(iNEXTRich_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iNEXT_est",
  class_n = 10, naColour = "white",
  legendTitle = "Species class",
  class_Style = "fisher",
  title = "iNEXT species richness"))

###### d. iChao +sp. ####
(iChaoIncrease_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iChao_increase",
  class_n = 7, naColour = "white",
  legendTitle = "Increase species class",
  class_Style = "fisher",
  title = "iChao increase (sp.)"))

###### e. iNEXT +sp. ####
(iNextIncrease_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iNEXT_increase",
  class_n = 7, naColour = "white",
  legendTitle = "Increase species class",
  class_Style = "fisher",
  title = "iNEXT increase (sp.)"))

    ###### f. iChao % ####
(iChaoIncreasePer_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iChao_increasePercent",
  class_n = 10, naColour = "white",
  legendTitle = "Increase class (%)",
  class_Style = "fisher",
  title = "iChao increase (%)"))

    ###### g. iNEXT % ####
(iNextIncreasePer_map <- countryMapFunction(
  map_in = dataMap %>% 
    dplyr::filter(level == "Country"),
  mapColumn = "iNEXT_increasePercent",
  class_n = 10, naColour = "white",
  legendTitle = "Increase class (%)",
  class_Style = "fisher",
  title = "iNEXT increase (%)"))

    ###### h. combine ####
  # Combine the plots 
  allCountryMaps <- cowplot::plot_grid(obsRich_map, NULL,
                                 iChaoRich_map, iNEXTRich_map,
                                 iChaoIncrease_map, iNextIncrease_map,
                                 iChaoIncreasePer_map, iNextIncreasePer_map,
                                     labels = c("(A)","","(B)", "(C)", "(D)","(E)",
                                                "(F)", "(G)"),
                                     ncol = 2, align = 'v', axis = 'l')
  # Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","4.4h_countryMaps.pdf"),
                   plot = allCountryMaps,
                   base_width = 12,
                   base_height = 10)




  ##### 4.5 Continent maps ####
    ###### a. observed richness ####
(obsRich_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "observedRichness",
  naColour = "gray90",
  legendTitle = "Species richness",
  title = "Observed species richness"))

    ###### b. iChao richness ####
(iChaoRich_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iChao_est",
  naColour = "gray90",
  legendTitle = "Species estimate",
  title = "iChao species richness"))

    ###### c. iNEXT richness ####
(iNEXTRich_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iNEXT_est",
  naColour = "gray90",
  legendTitle = "Species estimate",
  title = "iNEXT species richness"))

    ###### d. iChao +sp. ####
(iChaoIncrease_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iChao_increase",
  naColour = "gray90",
  legendTitle = "Increase species",
  title = "iChao increase (sp.)"))

    ###### e. iNEXT +sp. ####
(iNextIncrease_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iNEXT_increase",
  naColour = "gray90",
  legendTitle = "Increase species",
  title = "iNEXT increase (sp.)"))

    ###### d. iChao % ####
(iChaoIncreasePer_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iChao_increasePercent",
  naColour = "gray90",
  legendTitle = "Increase (%)",
  title = "iChao increase (%)"))

    ###### e. iNEXT % ####
(iNextIncreasePer_mapContinent <- continentMapFunction(
  map_in = continentDataMap %>% 
    dplyr::filter(level == "Continent"),
  mapColumn = "iNEXT_increasePercent",
  naColour = "gray90",
  legendTitle = "Increase (%)",
  title = "iNEXT increase (%)"))


    ###### f. combine ####
# Combine the plots 
allcontinentMaps <- cowplot::plot_grid(obsRich_mapContinent, NULL,
                                      iChaoRich_mapContinent, iNEXTRich_mapContinent,
                                      iChaoIncrease_mapContinent, iNextIncrease_mapContinent,
                                      iChaoIncreasePer_mapContinent, iNextIncreasePer_mapContinent,
                                      labels = c("(A)","","(B)", "(C)", "(D)","(E)",
                                                 "(F)","(G)"),
                                      ncol = 2, align = 'v', axis = 'l')
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","4.5h_continentMaps.pdf"),
                   plot = allcontinentMaps,
                   base_width = 12,
                   base_height = 9)

  ##### 4.6 Publication map ####
# Combine the plots 
(contCountry_map <- cowplot::plot_grid( iChaoIncrease_map + ggplot2::ggtitle("Undescribed species (iChao)"),    
                                        iChaoIncrease_mapContinent + ggplot2::ggtitle("Undescribed species (iChao)"),   
                                        iChaoIncreasePer_map + ggplot2::ggtitle("Undescribed species (%)"),  
                                        iChaoIncreasePer_mapContinent + ggplot2::ggtitle("Undescribed species (%)"),
                                        labels = c("(A)","(B)", "(C)", "(D)"),
                                        ncol = 2, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","4.6_contCountry_Maps.pdf"),
                   plot = contCountry_map,
                   base_width = 15,
                   base_height = 7)


#### 5.0 Correlates ####
  ##### 5.1 Data inputs ####
    ###### a. GDP_c ####
  # GPD-c data from https://ourworldindata.org/grapher/gdp-per-capita-maddison on the 6th of August 2024
GDP_c <- readr::read_csv("CorrelatesData/gdp-per-capita-maddison.csv", guess_max = 22000) %>% 
  dplyr::select(!`900793-annotations`) %>%
  dplyr::rename(year_GDPc = Year, GPD_per_cap = `GDP per capita`) %>%
    # Select for years prior to 2024 then select the most recent year 
  dplyr::filter(year_GDPc < 2025) %>%
  dplyr::group_by(Entity) %>% dplyr::arrange(year_GDPc %>% dplyr::desc()) %>%
  dplyr::filter(dplyr::row_number() == 1) %>% 
    # Add in missing GDP-C for Oceania countries
  dplyr::bind_rows(
      # Accessed on 22/Aug/2024 -- 
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=FJ 
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=PF
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=GU
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=NC
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=MP
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=PW
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=PG
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=WS
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=SB
    # https://data.worldbank.org/indicator/NY.GDP.PCAP.CD?locations=VU
    dplyr::tibble(Entity = c("Fiji", "French Polynesia", "Guam","New Caledonia",
                             "Northern Mariana Islands","Palau","Papua New Guinea",
                             "Samoa", "Solomon Islands", "Vanuatu"),
                  year_GDPc = c(2023, 2022, 2022,2022,
                                2020,2023,2023, 2023,
                                2023, 2023),
                  GPD_per_cap = c(5868, 18984, 40227,35745,
                                  17303,14563,2995,
                                  4139,2203,3367))
  )

    ###### b. TertiaryEd ####
  # Tertiary education from https://ourworldindata.org/grapher/share-of-the-population-with-completed-tertiary-education on the 6th of August 2024
TertiaryEd <- readr::read_csv("CorrelatesData/share-of-the-population-with-completed-tertiary-education.csv",
                              guess_max = 5000) %>%
  dplyr::rename(year_TerEd = Year, 
                percent_ed = `Combined - percentage of 25-64 years adults with incomplete tertiary education`) %>%
  # Select for years prior to 2024 then select the most recent year 
  dplyr::filter(year_TerEd < 2025) %>%
  dplyr::group_by(Entity) %>% dplyr::arrange(year_TerEd %>% dplyr::desc()) %>%
  dplyr::filter(dplyr::row_number() == 1)

    ###### c. cleanRecords ####
  # Calculate the number of clean records
    # BY COUNTRY
cleanRecords_country <- country_speciesCounts %>%
  dplyr::group_by(country_suggested) %>% 
  dplyr::summarise(country_n = sum(n)) 
    # BY CONTINENT
# Calculate the number of clean records
cleanRecords_continent <- continentOccs %>%
  dplyr::group_by(continent) %>% 
  dplyr::summarise(continent_n = sum(n))
  
    ###### d. globalDEM ####
  # Global DEM at 30s resolution
    # Initial download:
#globalDEM <- geodata::elevation_global(res = 0.5, path = "CorrelatesData")
    # Re-read in again
globalDEM <- terra::rast("CorrelatesData/elevation/wc2.1_30s/wc2.1_30s_elev.tif")
  # Extract the raster data to the world map data
countryElevations <- terra::vect(worldMap) %>%
  terra::extract(globalDEM, .)
# Save these data
readr::write_excel_csv(countryElevations, "CorrelatesData/countryElevations.csv")
  # Get the variance per group
countryElevations_var <- countryElevations %>%
  tidyr::drop_na() %>%
  dplyr::group_by(ID) %>%
  dplyr::mutate(min = min(wc2.1_30s_elev, na.rm = TRUE),
                max = max(wc2.1_30s_elev, na.rm = TRUE),
                elevationalRange = max - min) %>%
  dplyr::distinct(ID, .keep_all = TRUE) %>%
    # Add in the country names
  dplyr::left_join( worldMap %>%
                      dplyr::select(name_long, name, iso_a3_eh) %>%
                      sf::st_drop_geometry() %>%
                      dplyr::mutate(ID = dplyr::row_number()),
                    by = "ID") %>%
  # Rename some countries
  dplyr::mutate(
    name_long = dplyr::if_else(name_long == "Cape Verde", "Republic of Cabo Verde", name_long),
    name_long = dplyr::if_else(name_long == "The Gambia", "Gambia", name_long),
    name_long = dplyr::if_else(name_long == "Lao PDR", "Lao People's Democratic Republic", name_long),
    name_long = dplyr::if_else(name_long == "North Macedonia", "Macedonia", name_long),
    name_long = dplyr::if_else(name_long == "United States", "United States of America", name_long),
    name_long = dplyr::if_else(name_long == "Kingdom of eSwatini", "Eswatini", name_long)
  ) %>%
  dplyr::bind_rows(
    dplyr::tibble(name_long = c("French Guiana", "Martinique"),
                    elevationalRange = c(851,1397 )))
  # Save these data
readr::write_excel_csv(countryElevations_var, "CorrelatesData/countryElevations_var.csv")
if(!exists("countryElevations_var")){
  countryElevations_var <- readr::read_csv("CorrelatesData/countryElevations_var.csv")
}

  ###### e. roads ####
  # Read in the roads dataset
roads <- terra::vect("CorrelatesData/groads-v1-global-gdb/gROADS_v1.gdb") %>%
  sf::st_as_sf() %>%
    # Simplify the roads dataset
  dplyr::select(SOURCEID) 

  # Build a grid of points from which to find the nearest neighbour
  # point resolution
densgrid = 0.08 # c.a. 1 km at equator
  # Get a chunk size 
chunkSize = 100000
  # Make a 1 km grid cell of the world
grid <- sf::st_make_grid(worldMap, cellsize = densgrid, what = "centers") 
  # Find the points that are over land and extract their index
landPoints <- sf::st_intersects(grid, worldMap, sparse = TRUE) %>% 
  # return a tibble with the index of each match or NA where there was no match
  dplyr::tibble(indexMatch = .)  %>%
  dplyr::mutate(indexMatch = indexMatch %>% as.character() %>%
                  # deal with problems - Take the first number where two are provided
                  stringr::str_extract("[0-9]+") %>% 
                  # Remove zero to NA
                  stringr::str_replace("^[0]$", NA_character_) %>%
                  # Make numeric
                  as.numeric()
  ) %>%
  dplyr::mutate(point_num = dplyr::row_number()) %>%
  tidyr::drop_na() %>%
  dplyr::pull(point_num)
  # Filter to the land-only points 
grid <- grid[landPoints] %>% 
  sf::st_as_sf()
  # Get the nearest road to the points
nearest <- grid  %>%
  # Group by the row number and step size
  dplyr::group_by(group = ceiling(dplyr::row_number()/chunkSize)) %>%
  # Split the dataset up into a list by group
  dplyr::group_split(.keep = TRUE) %>%
  # Run the actual function
  parallel::mclapply(X = ., 
                     FUN = sf::st_nearest_feature,
                     roads,
                     mc.cores = 2) %>%
  # Combine the lists 
  unlist()

    #     # Find the index of the nearest road to each grid point
    #   nearest <- sf::st_nearest_feature(grid, roads) 

  # Make a list of pairs of points and their nearest road
pointPolyList <- grid %>%
  dplyr::bind_cols(roads[nearest,]) %>% 
  dplyr::mutate(number = dplyr::row_number()) %>% 
  dplyr::group_by(number) %>%
  dplyr::group_split()

  # Find the countries that don't have any roads accoridng to the roads dataset
noRoadCountries_index <- sf::st_intersects(roads, worldMap, sparse = TRUE) %>% 
  # return a tibble with the index of each match or NA where there was no match
  dplyr::tibble(indexMatch = .)  %>%
  dplyr::mutate(indexMatch = indexMatch %>% as.character() %>%
                  # deal with problems - Take the first number where two are provided
                  stringr::str_extract("[0-9]+") %>% 
                  # Remove zero to NA
                  stringr::str_replace("^[0]$", NA_character_) %>%
                  # Make numeric
                  as.numeric()
  ) 
  # Find the countries where roads don't overlap with them
noRoadCountries <- worldMap[dplyr::symdiff(TEST$indexMatch, 1:241),] %>%
  dplyr::pull(name_long)


  # Make a function to go over each pair in the list and find the distance between the two
nearest_distance <- function(inList){
  point = inList %>% dplyr::select(x)
  poly = inList %>% sf::st_drop_geometry() %>% sf::st_as_sf()
  dist <- sf::st_distance(point, poly) %>% 
    dplyr::tibble()
  return(dist)
}

  # Run the function in parallel to get the nearest neighbors 
distances <- parallel::mclapply(pointPolyList, 
                                FUN = nearest_distance,
                                mc.cores = 5) %>%
  dplyr::bind_rows()

  # Add the distance to nearest road onto the grid
roadDistanceGrid <- grid %>%
  dplyr::bind_cols(distances %>% setNames("distance")) 

# Simplify the world map ONCE to be used later
simplePoly <- worldMap %>% 
  dplyr::select(name_long) %>% 
  sf::st_drop_geometry() %>%
  dplyr::mutate(indexMatch = dplyr::row_number())

#Extract polygon information to points
extracted <- sf::st_intersects(roadDistanceGrid, worldMap, sparse = TRUE) %>% 
  # return a tibble with the index of each match or NA where there was no match
  dplyr::tibble(indexMatch = .) 
# If first element is full, unlist each one
extracted <- extracted %>%
  dplyr::mutate(indexMatch = indexMatch %>% as.character() %>%
                  # deal with problems - Take the first number where two are provided
                  stringr::str_extract("[0-9]+") %>% 
                  # Remove zero to NA
                  stringr::str_replace("^[0]$", NA_character_) %>%
                  # Make numeric
                  as.numeric()
  ) %>%
  # drop geometry
  sf::st_drop_geometry() 

# rejoin
countryDistances <- extracted %>%
  dplyr::left_join(simplePoly,
                   by = "indexMatch") %>%
  # Add in the database_id
  dplyr::bind_cols(roadDistanceGrid) %>%
  dplyr::group_by(name_long) %>%
  dplyr::mutate(mean_roadDist = mean(distance),
                median_roadDist = median(distance)) %>% 
  dplyr::distinct(name_long, .keep_all = TRUE) %>% 
  dplyr::select(!c(indexMatch, distance, x)) %>%
  dplyr::mutate(
    name_long = dplyr::if_else(name_long == "Kingdom of eSwatini","Eswatini", name_long),
    name_long = dplyr::if_else(name_long == "Lao PDR", "Lao People's Democratic Republic", name_long),
    name_long = dplyr::if_else(name_long == "North Macedonia", "Macedonia", name_long),
    name_long = dplyr::if_else(name_long == "The Gambia", "Gambia", name_long),
    name_long = dplyr::if_else(name_long == "United States", "United States of America", name_long)) %>%
    # Remove the noRoadCountries
  dplyr::filter(!name_long %in% noRoadCountries)
  

# Save these data
readr::write_csv(countryDistances, "Table_outputs/5.1e_countryDistances.csv")
if(!exists("countryDistances")){
  countryDistances <- readr::read_csv("Table_outputs/5.1e_countryDistances.csv")
}


  ###### f. lit records ####

country_speciesChecklistCounts <- country_speciesChecklistCounts %>% 
# Change some country names to better match the continent data
dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Barbuda",
                                                  "Antigua and Barbuda", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels",
                                                  "Belgium", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Caribbean Netherlands",
                                                  # Renamed Curaçao to match to continent
                                                  "Curaçao", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cocos Islands",
                                                  "Indian Ocean Territories", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Federation of Bosnia and Herzegovina",
                                                  "Bosnia and Herzegovina", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "French Guiana",
                                                  "Brazil", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gibraltar",
                                                  "Spain", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Guadeloupe",
                                                  "Barbados", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Martinique",
                                                  "Curaçao", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Mayotte",
                                                  "Madagascar", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Réunion",
                                                  "Mauritius", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Svalbard Islands",
                                                  "Norway", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "West Bank",
                                                  "Palestine", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "North Macedonia",
                                                  "Macedonia", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Kingdom of eSwatini",
                                                  "Eswatini", rNaturalEarth_name),
              rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "N. Cyprus",
                                                  "Cyprus", rNaturalEarth_name)) 

  # Get the count of no-occurrence-record (i.e., literature) species per country
litCountryCount <- country_speciesChecklistCounts %>% 
  dplyr::filter(dataFrom == "checklist") %>%
  countryHarmoniseR(countryColumn = "rNaturalEarth_name") %>% 
  dplyr::count(rNaturalEarth_name, name = "n_lit")
occSpecies_CountryCountry <- country_speciesChecklistCounts %>% 
  dplyr::filter(dataFrom == "points") %>%
  dplyr::count(rNaturalEarth_name, name = "n_occ") 
  # Combine these counts and then get a proportion of occurrence-based species
litOcc_counts <- litCountryCount %>% 
  dplyr::left_join(occSpecies_CountryCountry %>%
                     countryHarmoniseR(countryColumn = "rNaturalEarth_name"),
                     by = "rNaturalEarth_name")  %>%
    # Turn NA into zero
  dplyr::mutate(n_lit = dplyr::if_else(is.na(n_lit),
                                       0, n_lit),
                n_occ = dplyr::if_else(is.na(n_occ),
                                       0, n_occ)) %>% 
  dplyr::mutate(propOccurrences = n_occ/(n_occ+n_lit))
  

    ###### g. type specimens ####
  # Read in the John Ascher Type specimen database
typeSpecimens <- readxl::read_excel("Ascher bee files Dec 2024.xlsx", sheet = "North America",
                                    guess_max = 12000) %>%
  dplyr::bind_rows(readxl::read_excel("Ascher bee files Dec 2024.xlsx", sheet = "South America",
                                      guess_max = 12000),
                   readxl::read_excel("Ascher bee files Dec 2024.xlsx", sheet = "Old World",
                                      guess_max = 12000))
  # Read in Michael Orr's translation document
stateCipher <- readxl::read_excel("Ascher primary division political area codes.xlsx", sheet = "Primary divisions",
                                    guess_max = 12000) %>%
    # Create a country:state column to make the codes unique
  dplyr::mutate(countryStateCode = stringr::str_c(country, abv, sep = "_"))
  # Read in the BeeBDC Ashcer country cipher
countryCipher <- readxl::read_excel("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/DL country codes according to Ascher Nov 2020.xlsx",
                                    sheet = "Country codes used by DL accord",
                                    guess_max = 12000)
# Read in a manually-checked comparison sheet
Manual_ISO <- read.csv("/Users/jamesdorey/Desktop/Uni/Packages/BeeBDC_development/restrictedScripts/ISO_matching.csv") %>%
  # Modify some country names to match rnaturalearth
  dplyr::mutate(rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cape Verde","Republic of Cabo Verde", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Republic of Congo","Republic of the Congo", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "French Southern Territories","French Southern and Antarctic Lands", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Macedonia"),"North Macedonia", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Reunion"),"Réunion", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Swaziland"),"Kingdom of eSwatini", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "South Georgia and South Sandwich Islands"),"South Georgia and the Islands", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Aland Islands"),"Åland Islands", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Bouvet"),"Norway", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Falkland Islands"),"Falkland Islands / Malvinas", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Darussalam"),"Brunei Darussalam", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Scotland|Wales|England"),"United Kingdom", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Antigua"),"Antigua and Barbuda", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Barbuda"),"Antigua and Barbuda", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Brussels|Flemish Region|Walloon Region"),"Belgium", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Republic Srpska|Federation of Bosnia and Herzegovina"),"Bosnia and Herzegovina", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Northern") & 
                                                      stringr::str_detect(rNaturalEarth_name, "Cyprus"),
                                                    "Northern Cyprus", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(stringr::str_detect(rNaturalEarth_name, "Gibraltar"),"Spain", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Barbuda","Antigua and Barbuda", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Brussels","Belgium", rNaturalEarth_name),
                  # Renamed Curaçao to match to continent
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Caribbean Netherlands","Curaçao", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Cocos Islands","Indian Ocean Territories", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "French Guiana","Brazil", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gaza","Palestine", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Gibraltar","Spain", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Guadeloupe","Barbados", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Martinique","Curaçao", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Mayotte","Madagascar", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Réunion","Mauritius", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Svalbard Islands","Norway", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "West Bank","Palestine", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Tokelau","New Zealand", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "USA","United States", rNaturalEarth_name),
                rNaturalEarth_name = dplyr::if_else(rNaturalEarth_name == "Vojvodina","Bosnia and Herzegovina", rNaturalEarth_name)
                )

CheckLmap <-  full_join(worldMap,  Manual_ISO,
                        by = c("name_long" = "rNaturalEarth_name")) %>%
  sf::st_drop_geometry()

  # Translate the type country codes
typeSpecimens_edit <- typeSpecimens %>% 
  # Make an index column 
  dplyr::mutate(index = dplyr::row_number(), .before = 1) %>%
    # Remove special characters from type country code 
  dplyr::mutate(`Type country` = stringr::str_extract(`Type country`, "[A-Z]+")) %>% 
  dplyr::left_join(CheckLmap %>% dplyr::select(Code, name_long),
                    by = c("Type country" = "Code"),
                   relationship = "many-to-many") %>% 
  dplyr::rename(type_country = name_long) %>%
    # Remove a data quirk; Siachen Glacier
  dplyr::mutate(type_country = dplyr::if_else(type_country == "Siachen Glacier", NA_character_, type_country)) %>%
    # make a countryStateCode column
  dplyr::mutate(countryStateCode = stringr::str_c(`Type country`, `Type state`, sep = "_")) %>%
    # Match the state names
  dplyr::left_join(stateCipher %>% 
                     dplyr::select(countryStateCode, `PRIMARY DIVISION FIELD USE THIS`) %>%
                     dplyr::filter(complete.cases(countryStateCode)),
                    by = c("countryStateCode"),
                   relationship = "many-to-many") %>%
    # Keep unique only by index
  dplyr::distinct(.keep_all = TRUE)
  
  # Extract the failed typeSpecimens_
typeSpecimens_failed <- typeSpecimens_edit %>% 
  # Select the types without a Type country in text
  dplyr::filter(is.na(type_country)|type_country == "") %>%
  dplyr::filter(complete.cases(Lat)) %>% 
  sf::st_as_sf(coords = c("Lat", "Lon"), crs = sf::st_crs(worldMap)) %>% 
  sf::st_intersection(worldMap %>% dplyr::select(name_long)) %>%
  dplyr::select(!type_country) %>% 
  dplyr::rename(type_country = name_long) %>%
  sf::st_drop_geometry()


TEST <- typeSpecimens_edit %>% 
  # Drop the types without a Type country in text
  dplyr::filter(!is.na(type_country)|type_country == "") 
  
  
###### h. area ####
# Get country areas
countryArea <- dplyr::tibble(area_m = sf::st_area(worldMap)) %>% 
  dplyr::mutate(area_m = as.numeric(area_m)) %>% 
  dplyr::bind_cols(worldMap %>%
                     sf::st_drop_geometry()) %>% 
  countryHarmoniseR(countryColumn = "name_long") %>% 
  dplyr::select(area_m, name_long)


  ###### i. combine ####
  # Combine the our world in data
ourWorldinData <- GDP_c %>%
  dplyr::ungroup() %>%
    # Harmonise the country names
  countryHarmoniseR(
    data = .,
    countryColumn = "Entity"
  ) %>%
  dplyr::left_join(., TertiaryEd %>%
                     countryHarmoniseR(
                       data = .,
                       countryColumn = "Entity"
                       ) %>% 
                     dplyr::select(!Code),
                   by = "Entity") 

  # combine these statistics
combined_explanatory <- combinedStatistics %>%
    # Add in the number of clean records
  dplyr::left_join(cleanRecords_country %>% 
                     countryHarmoniseR(countryColumn = "country_suggested") %>%
                     dplyr::rename(cleanRecords = country_n), 
                   by = c("name" = "country_suggested")) %>%
    # Add in the elevational variation
  dplyr::left_join(countryElevations_var %>% dplyr::ungroup() %>% 
                     countryHarmoniseR(countryColumn = "name_long") %>%
                     dplyr::select(elevationalRange, name_long),
                   by = c("name" = "name_long")) %>%
    # Add in our world in data
  dplyr::left_join(ourWorldinData,
                   by = c("name" = "Entity")) %>%
    # Remove continent and global markers
  dplyr::filter(!level %in% c("Continental", "Global")) %>% 
    # Add in continent names
  dplyr::left_join(checklistFile %>% 
                     dplyr::select(rNaturalEarth_name, continent) %>%
                     countryHarmoniseR(countryColumn = "rNaturalEarth_name") %>%
                     dplyr::distinct() %>%
                     dplyr::filter(complete.cases(continent)),
                   by = c("name" = "rNaturalEarth_name")) %>%
  dplyr::left_join(countryDistances %>% 
                     countryHarmoniseR(countryColumn = "name_long"),
                   by = c("name" = "name_long")) %>%
  dplyr::left_join(litOcc_counts %>% 
                     countryHarmoniseR(countryColumn = "rNaturalEarth_name"),
                   by = c("name" = "rNaturalEarth_name")) %>%
  countryHarmoniseR(countryColumn = "name") %>% 
  dplyr::left_join(., countryArea, by = c("name" = "name_long")) 
  # Add in continental data that didn't match
combined_explanatory <- combined_explanatory %>%
  dplyr::mutate(continent = dplyr::if_else(is.na(continent) & name == "Antigua and Barbuda",
                                          "North America", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Saint-Martin",
                                          "North America", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Curaçao",
                                          "South America", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Eswatini",
                                          "Africa", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Cyprus",
                                          "Europe", continent),
                continent = dplyr::if_else(is.na(continent) & name == "French Guiana",
                                          "South America", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Bosnia and Herzegovina",
                                          "Europe", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Belgium",
                                          "Europe", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Palestine",
                                          "Asia", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Czechia",
                                          "Europe", continent),
                continent = dplyr::if_else(is.na(continent) & name == "Russian Federation",
                                            "Europe", continent)
                )


    # Save this file
  readr::write_excel_csv(combined_explanatory, "CorrelatesData/combined_explanatory.csv")

if(!exists("combined_explanatory")){
  combined_explanatory <- readr::read_csv("CorrelatesData/combined_explanatory.csv")
}

  # Remove extra continents
combined_explanatory <- combined_explanatory %>%
  dplyr::filter(!stringr::str_detect(continent, "Seven"))

      ###### j. explore in 3D ####
# Examine the relationships in 3D
plotly::plot_ly(data = combined_explanatory,
                x = ~(propOccurrences),  
                y = ~log(cleanRecords), 
                z = ~log(iChao_est), 
                type = "scatter3d", mode = "markers",
                name = ~name,
                color = ~continent, colors = "Dark2")

    ###### k. area and richness ####
  # Make and harmonise a list of island states
islandStates <- dplyr::tibble(islands = c("Antigua and Barbuda", "Bahamas", "Bahrain", "Barbados", "Brunei", "Cape Verde", "Comoros", "Cuba", "Cyprus", "Dominica", "Dominican Republic", "East Timor", "Fiji", "Grenada", "Haiti", "Iceland", "Indonesia", "Ireland", "Jamaica", "Japan", "Kiribati", "Madagascar", "Maldives", "Malta", "Marshall Islands", "Mauritius", "Micronesia", "Nauru", "New Zealand(Aotearoa)", "Palau", "Papua New Guinea", "Philippines", "Saint Kitts and Nevis", "Saint Lucia", "Saint Vincent and the Grenadines", "Samoa", "São Tomé and Príncipe", "Seychelles", "Singapore", "Solomon Islands", "Sri Lanka", "Tonga", "Trinidad and Tobago", "Tuvalu", "United Kingdom", "Vanuatu", "Northern Cyprus", "Taiwan", "Cook Islands", "Niue", "Åland", "American Samoa", "Andaman and Nicobar Islands", "Anguilla", "Aruba", "Bermuda", "Bouvet Island", "British Indian Ocean Territory", "British Virgin Islands", "Caribbean Netherlands", "Cayman Islands", "Christmas Island", "Cocos (Keeling) Islands", "Curaçao", "Easter Island", "Falkland Islands", "Faroe Islands", "French Polynesia", "French Southern and Antarctic Lands", "Galápagos Islands", "Greenland", "Guadeloupe", "Guam", "Guernsey", "Heard and McDonald Islands", "Hong Kong", "Isle of Man", "Jan Mayen", "Jersey", "Juan Fernández Islands", "Lakshadweep", "Macau", "Martinique", "Mayotte", "Montserrat", "New Caledonia", "Norfolk Island", "Northern Mariana Islands", "Pitcairn", "Henderson", "Ducie", "and Oeno Islands", "Puerto Rico", "Réunion", "San Andrés and Providencia", "Saint Barthélemy", "Saint Helena", "Ascension", "and Tristan da Cunha", "Saint Martin", "Saint Pierre and Miquelon", "Sint Maarten", "South Georgia and the South Sandwich Islands", "Sovereign Base Areas of Akrotiri and Dhekelia", "Svalbard", "Tokelau", "Turks and Caicos Islands", "United States Minor Outlying Islands", "U.S. Virgin Islands", "Wallis and Futuna", "Timor-Leste")) %>%
  countryHarmoniseR(countryColumn = "islands") %>%
  dplyr::distinct()
  # Add in the species per area
combined_explanatory <- combined_explanatory %>%
  dplyr::mutate(sp_per_km = (iChao_est/(area_m/1000000))) %>%
    # Add in island status
  dplyr::mutate(island = dplyr::if_else(name %in% islandStates$islands,
                                        "Island", "Mainland")) 
  # Plot the data
IslMain_plot <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = island, y = sp_per_km %>% log())) +
  ggplot2::geom_boxplot() +
  ggplot2::theme_classic() + ggplot2::ylab("Species per square kilometer") + 
  ggplot2::xlab("Country status") 
  # Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","5.1k_Sp_per_km_IslMain.pdf"), 
                  plot = IslMain_plot, 
                  width = 4, height = 4, units = "in", dpi = 300)
  # Non-parametric test of differences
wilcox.test(combined_explanatory %>% dplyr::filter(island == "Island") %>% dplyr::pull(sp_per_km),
       combined_explanatory %>% dplyr::filter(island == "Mainland") %>% dplyr::pull(sp_per_km),
       alternative = "two.sided")


  ##### 5.2 GLMM iChao Increase ####
shapiro.test( log(combined_explanatory$iChao_increasePercent ))
hist(         log(combined_explanatory$iChao_increasePercent ))
plot(         log(combined_explanatory$iChao_increasePercent ))

    ###### a. build model ####
  # Make a lmer for the increase in the estimated number of species
(lmerIncrease <- lmerTest::lmer(data = combined_explanatory ,
                   formula = log(iChao_increase) ~ 
                      # Fixed effects
                     log(GPD_per_cap) + log(percent_ed) + 
                     log(median_roadDist) +
                     log(elevationalRange)  + 
                     log(area_m) +
                     log(observedRichness)  + 
                     # Removing cleanRecords below  will change the sign of propOccurrences to negative
                     log(cleanRecords) *
                     log(propOccurrences) +
                      # Random effect
                     (1|continent),
                   REML = TRUE)
)

      ###### b. examine outputs ####
# Get the output summary
summary(lmerIncrease)

  # PERFORM SOME MODEL TESTS AND EXAMINE THE OUTPUTS
# Extract the residuals 
tdat_in <- dplyr::tibble(predicted= predict(lmerIncrease), 
                          residual = residuals(lmerIncrease), 
                          continent= combined_explanatory  %>% 
                           dplyr::ungroup() %>% 
                           dplyr::filter(dplyr::row_number() %in% c(names(predict(lmerIncrease)) %>%
                                           as.numeric())) %>% 
                           dplyr::pull(continent))
  # Check the model's residuals
ggplot2::ggplot(tdat_in,
                ggplot2::aes(x=predicted,y=residual, colour=continent)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_hline(yintercept=0, lty=3)
  # Histogram
ggplot2::ggplot(tdat_in,aes(x=residual)) + ggplot2::geom_histogram(bins=20, color="black")
  # qqplot
ggplot2::ggplot(tdat_in,aes(sample=residual)) + ggplot2::stat_qq() + ggplot2::stat_qq_line()


      ###### c. plot predictions ####
  # Extract the data that was predicted (no missing values) and add the model predictions
combinedPred <- combined_explanatory %>%
  dplyr::ungroup() %>% 
  dplyr::filter(dplyr::row_number() %in% c(names(predict(lmerIncrease)) %>%
                                             as.numeric())) %>%
  dplyr::mutate(prediction = predict(lmerIncrease))


# Plot GPD_per_cap
(GPD_per_cap_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(GPD_per_cap),
                                                               y=prediction,
                                                              colour=continent, fill = continent,
                                                              group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot percent_ed
(percent_ed_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(percent_ed),
                                                              y=prediction,
                                                             colour=continent, fill = continent,
                                                             group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot median_roadDist
(median_roadDist_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(median_roadDist),
                                                               y=prediction,
                                                               colour=continent, fill = continent,
                                                               group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot elevationalRange
(elevationalRange_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(elevationalRange),
                                                               y=prediction,
                                                               colour=continent, fill = continent,
                                                               group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot area_m
(area_m_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(area_m),
                                                               y=prediction,
                                                         colour=continent, fill = continent,
                                                         group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot observedRichness
(observedRichness_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(observedRichness),
                                                               y=prediction,
                                                               colour=continent, fill = continent,
                                                               group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot clean records
(CleanRecords_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(cleanRecords),
                                          y=prediction,
                                          colour=continent, fill = continent,
                                          group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot propOccurrences
(propOccurrences_5_2 <- ggplot2::ggplot(combinedPred,ggplot2::aes(x=log(propOccurrences),
                                                               y=prediction,
                                                               colour=continent, fill = continent,
                                                               group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao richness increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))


# Combine the plots 
(effects_in_predicted <- cowplot::plot_grid(GPD_per_cap_5_2, percent_ed_5_2,
                                  median_roadDist_5_2, elevationalRange_5_2, 
                                  area_m_5_2, observedRichness_5_2,  
                                  CleanRecords_5_2, propOccurrences_5_2,
                                  label_y = 1.015,
                                  label_x = 0.08,
                                  labels = c("(A)*","(B)",
                                             "(C)","(D)",
                                             "(E)", "(F)***",
                                             "(G)***", "(H)***"),
                                  ncol = 2, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.2c_effectsIncrease_predicted.pdf"),
                   plot = effects_in_predicted,
                   base_width = 10,
                   base_height = 13)


      ###### d. plot variables ####
showLabels = FALSE
  # Plot the variables
# GDP-c
(GPDc_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(GPD_per_cap), 
                                                                      y = log(iChao_increase), 
                                        group = continent, color = continent, fill = continent,
                                        label = name)) + 
  ggplot2::geom_point() + 
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log gross domestic capital (GDP)") +
  ggplot2::geom_smooth(method = "lm") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# percent_ed
(perEd_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(percent_ed), 
                                                                       y = log(iChao_increase), 
                                                                  group = continent, 
                                                                  color = continent,
                                                                  fill = continent)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log percentage higher education") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
  ggplot2::theme_classic())
# Roads
(roads_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(median_roadDist %>%
                                                                                 as.numeric()), 
                                                                       y = log(iChao_increase), 
                                                                       group = continent, 
                                                                       color = continent,
                                                                       fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log median distance to nearest road") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# elevationalRange
(elev_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(elevationalRange), 
                                                                      y = log(iChao_increase), 
                                                                      group = continent, 
                                                                      color = continent,
                                                                      fill = continent)) + 
    ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log elevational range") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
  ggplot2::theme_classic())
# area_m
(area_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(area_m), 
                                                                      y = log(iChao_increase), 
                                                                      group = continent, 
                                                                      color = continent,
                                                                      fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log area (m)") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# cleanRecords
(clean_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(cleanRecords), 
                                                                       y = log(iChao_increase), 
                                                                  group = continent, 
                                                                  color = continent, fill = continent,
                                                                  label = name)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log clean bee occurrence records") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
  ggplot2::theme_classic())
# observedRichness
(obsRich_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(observedRichness), 
                                                                         y = log(iChao_increase), 
                                                                         group = continent, 
                                                                         color = continent,
                                                                         fill = continent)) + 
    ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log observed species richness") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
  ggplot2::theme_classic())
# propOccurrences
(propOcc_in <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(propOccurrences), 
                                                                         y = log(iChao_increase), 
                                                                         group = continent, 
                                                                         color = continent,
                                                                         fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao richness increase") + ggplot2::xlab("Log proportion of occurrence-based species") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())

# Combine the plots 
(effects_in <- cowplot::plot_grid(GPDc_in, perEd_in,
                                  roads_in, elev_in,
                                  area_in, obsRich_in,  
                                  clean_in, propOcc_in, #obsRichprop_in,
                                  label_y = 1.015,
                                  label_x = 0.08,
                                   labels = c("(A)*","(B)",
                                              "(C)","(D)",
                                              "(E)", "(F)***",
                                              "(G)***", "(H)***"),
                                   ncol = 2, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.2d_effectsIncrease.pdf"),
                   plot = effects_in,
                   base_width = 10,
                   base_height = 13)

                                                                                                  

    ##### 5.3 GLMM iChao percent ####
      ###### a. build model ####
  # Make a lmer for the increase in the estimated number of species
(lmerPercentiChaoPer <- lmerTest::lmer(data = combined_explanatory,
                                formula = log(iChao_increasePercent) ~ 
                                  # Fixed effects
                                 log(GPD_per_cap) + log(percent_ed) + log(median_roadDist) +
                                  log(elevationalRange) + 
                                  log(area_m) +
                                  log(observedRichness) +
                                  # Removing cleanRecords below  will change the sign of propOccurrences to negative
                                 log(cleanRecords) * 
                                  log(propOccurrences) +
                                  # Random effect
                                 (1|continent))
 
)

      ###### b. examine outputs ####
  # Get the output summary
summary(lmerPercentiChaoPer)




# PERFORM SOME MODEL TESTS AND EXAMINE THE OUTPUTS
# Extract the residuals 
tdat_per <- dplyr::tibble(predicted= predict(lmerPercentiChaoPer), 
                         residual = residuals(lmerPercentiChaoPer), 
                         continent= combined_explanatory  %>% 
                           dplyr::ungroup() %>% 
                           dplyr::filter(dplyr::row_number() %in% c(names(predict(lmerPercentiChaoPer)) %>%
                                                                      as.numeric())) %>% 
                           dplyr::pull(continent))
# Check the model's residuals
ggplot2::ggplot(tdat_per,
                ggplot2::aes(x=predicted,y=residual, colour=continent)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_hline(yintercept=0, lty=3)
# Histogram
ggplot2::ggplot(tdat_per,aes(x=residual)) + ggplot2::geom_histogram(bins=20, color="black")
# qqplot
ggplot2::ggplot(tdat_per,aes(sample=residual)) + ggplot2::stat_qq() + ggplot2::stat_qq_line()


    ###### c. plot predictions ####
# Extract the data that was predicted (no missing values) and add the model predictions
combinedPred_perc <- combined_explanatory %>%
  dplyr::ungroup() %>% 
  dplyr::filter(dplyr::row_number() %in% c(names(predict(lmerPercentiChaoPer)) %>%
                                             as.numeric())) %>%
  dplyr::mutate(prediction = predict(lmerPercentiChaoPer))


# Plot GPD_per_cap
(GPD_per_cap_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(GPD_per_cap),
                                                              y=prediction,
                                                              colour=continent, fill = continent,
                                                              group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot percent_ed
(percent_ed_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(percent_ed),
                                                             y=prediction,
                                                             colour=continent, fill = continent,
                                                             group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot median_roadDist
(median_roadDist_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(median_roadDist),
                                                                  y=prediction,
                                                                  colour=continent, fill = continent,
                                                                  group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot elevationalRange
(elevationalRange_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(elevationalRange),
                                                                   y=prediction,
                                                                   colour=continent, fill = continent,
                                                                   group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot area_m
(area_m_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(area_m),
                                                         y=prediction,
                                                         colour=continent, fill = continent,
                                                         group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot observedRichness
(observedRichness_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(observedRichness),
                                                                   y=prediction,
                                                                   colour=continent, fill = continent,
                                                                   group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot clean records
(CleanRecords_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(cleanRecords),
                                                               y=prediction,
                                                               colour=continent, fill = continent,
                                                               group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))
# Plot propOccurrences
(propOccurrences_5_3 <- ggplot2::ggplot(combinedPred_perc,ggplot2::aes(x=log(propOccurrences),
                                                                  y=prediction,
                                                                  colour=continent, fill = continent,
                                                                  group=continent)) + 
    ggplot2::geom_point() + ggplot2::theme_classic() + 
    ggplot2::ylab("Log predicted iChao percent increase") +
    ggplot2::geom_smooth(method = "lm") + 
    ggplot2::theme(legend.position="right"))


# Combine the plots 
(effects_perc_predicted <- cowplot::plot_grid(GPD_per_cap_5_3, percent_ed_5_3,
                                            median_roadDist_5_3, elevationalRange_5_3, 
                                            area_m_5_3, observedRichness_5_3,  
                                            CleanRecords_5_3, propOccurrences_5_3,
                                            label_y = 1.015,
                                            label_x = 0.08,
                                            labels = c("(A)*","(B)",
                                                       "(C)","(D)",
                                                       "(E)", "(F)***",
                                                       "(G)***", "(H)***"),
                                            ncol = 2, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.3c_effectsIncrease_predicted.pdf"),
                   plot = effects_perc_predicted,
                   base_width = 10,
                   base_height = 13)


      ###### d. plot variables ####
# Plot the variables
# GDP-c
(GPDc_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(GPD_per_cap), 
                                                                       y = log(iChao_increasePercent), 
                                                                       group = continent, 
                                                                       color = continent,
                                                                       fill = continent)) + 
   ggplot2::geom_point() + 
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log gross domestic capital (GDP)") +
    ggplot2::geom_smooth(method = "lm") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# percent_ed
(perEd_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(percent_ed), 
                                                                        y = log(iChao_increasePercent), 
                                                                        group = continent, 
                                                                        color = continent,
                                                                        fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log percentage higher education") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# road
(road_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(mean_roadDist %>% as.numeric()), 
                                                                       y = log(iChao_increasePercent), 
                                                                       group = continent, 
                                                                       color = continent,
                                                                       fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log median distance to nearest road") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# elevationalRange
(elev_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(elevationalRange), 
                                                                       y = log(iChao_increasePercent), 
                                                                       group = continent, 
                                                                       color = continent,
                                                                       fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log elevational range") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# area_m
(area_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(area_m), 
                                                                       y = log(iChao_increasePercent), 
                                                                       group = continent, 
                                                                       color = continent,
                                                                       fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log area (m)") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# cleanRecords
(clean_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(cleanRecords), 
                                                                        y = log(iChao_increasePercent), 
                                                                        group = continent, 
                                                                        color = continent,
                                                                        fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log clean bee occurrence records") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())
# observedRichness-c
(obsRich_per <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(observedRichness), 
                                                                          y = log(iChao_increasePercent), 
                                                                          group = continent, 
                                                                          color = continent,
                                                                          fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log observed bee species richness") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
        } +
    ggplot2::theme_classic() )
# propOccurrences
(propOcc_inPer <- ggplot2::ggplot(data = combined_explanatory, ggplot2::aes(x = log(propOccurrences), 
                                                                         y = log(iChao_increasePercent), 
                                                                         group = continent, 
                                                                         color = continent,
                                                                         fill = continent)) + 
    ggplot2::geom_point() + 
    ggplot2::geom_smooth(method = "lm") +
    ggplot2::ylab("Log iChao percent increase") + ggplot2::xlab("Log proportion of occurrence-based species") +
    {if(showLabels == TRUE)
      ggplot2::geom_text(hjust=0, vjust=0) 
    } +
    ggplot2::theme_classic())



# Combine the plots 
(effects_per <- cowplot::plot_grid(GPDc_per, perEd_per,
                                   road_per, elev_per,
                                   area_per, obsRich_per, 
                                   clean_per, propOcc_inPer, #obsRichprop_per,
                                   label_y = 1.015,
                                   label_x = 0.08,
                                   labels = c("(A)*","(B)",
                                              "(C)","(D)",
                                              "(E)", "(F)***",
                                              "(G)***", "(H)***"),
                                   ncol = 2, align = 'v', axis = 'l'))
# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.3d_effectsPercent.pdf"),
                   plot = effects_per,
                   base_width = 10,
                   base_height = 13)


    ##### 5.3 Increase + percent ####
# Make a lmer for the increase in the estimated number of species
(lmer_inPer <- lmerTest::lmer(data = combined_explanatory %>%
                                  # Remove zero percent increase country (Aruba) for this analysis
                                dplyr::filter(iChao_increasePercent > 0),
                                formula = log(iChao_increase) ~ 
                                  # Fixed effects
                          log(iChao_increasePercent) +
                                  # Random effect
                                  (1|continent))
)
# Get the output summary
(summary_lmin <- summary(lmer_inPer))
# Plot the interaction
(lm_inPer_plot <- ggplot2::ggplot(data = combined_explanatory, 
                                  ggplot2::aes(x = log(iChao_increase), y = log(iChao_increasePercent), 
                                               group = continent, color = continent, 
                                               label = name)) + 
  ggplot2::geom_point() + 
  ggplot2::geom_smooth(method = "lm") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(
      "p = ", summary_lmin$coefficients %>% tibble::as_tibble() %>% dplyr::pull(`Pr(>|t|)`) %>% .[2] %>% 
        round(10)
    )) + ggplot2::xlab("Log of species increase (iChao)") +
    ggplot2::ylab("Log of percentage species increase (iChao)")
)


# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","5.3_lm_inPer_plot.pdf"), 
                plot = lm_inPer_plot, 
                width = 6, height = 6, units = "in", dpi = 300)

  ##### 5.4 Clean + propOcc ####
# cleanRecords and propOccurrences
(obsRichprop <- ggplot2::ggplot(data = combined_explanatory, 
                                    ggplot2::aes(x = log(propOccurrences), 
                                                 y = log(cleanRecords),
                                                 group = continent, color = continent, fill = continent,
                                                 label = name)) + 
   ggplot2::geom_point() + 
   ggplot2::geom_smooth(method = "lm") +
   ggplot2::ylab("Log proportion of species w. occurrences") + ggplot2::xlab("Log clean bee occurrence records") +
   {if(showLabels == TRUE)
     ggplot2::geom_text(hjust=0, vjust=0) 
   } +
   ggplot2::theme_classic())

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","5.4_obsRichprop.pdf"), 
                plot = obsRichprop, 
                width = 6, height = 6, units = "in", dpi = 300)

  ##### 5.5 Taxonomic gap ####
    ###### a. plot percent ####
# plot the taxonomic gap
# Make a custom legend
(barLegend_gap <- ggplot2::ggplot(dplyr::tibble(name = c("yes","yes"),
                                            statistic = c("iChao", "iNEXT") %>%
                                              factor(levels = c("iChao", "iNEXT")),
                                            est = c(1,1)) ,
                              aes(x = name, y = est)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = est), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT"),
                              values = c("iChao" = "#55AD9B", "iNEXT" = "#FD9B63")) 
)

combined_explanatory_longer <- dplyr::bind_rows(
      # ichao
    combined_explanatory %>% dplyr::select(c("name", "observedRichness", 
                                             "continent",
                                             tidyselect::contains("iChao"))) %>%
      dplyr::mutate(lowerPercent = ((iChao_lower/observedRichness)-1)*100 %>%
                      dplyr::if_else(. <= 0, 0.0001, .),
                    upperPercent = ((iChao_upper/observedRichness)-1)*100) %>%
      dplyr::mutate(statistic = "iChao") %>%
      dplyr::mutate(Chao_increasePercent = iChao_increasePercent) %>% 
      setNames(colnames(.) %>% stringr::str_remove(., "iChao_|iNEXT_")),
      # iNEXT
    combined_explanatory %>% dplyr::select(c("name", "observedRichness", "iChao_increasePercent",
                                             "continent",
                                             tidyselect::contains("iNEXT"))) %>%
      dplyr::mutate(lowerPercent = ((iNEXT_lower/observedRichness)-1)*100,
                    lowerPercent = dplyr::if_else(lowerPercent <= 0, 0.0001, lowerPercent),
                    upperPercent = ((iNEXT_upper/observedRichness)-1)*100) %>%
      dplyr::select(!"niNEXT") %>% 
      dplyr::mutate(statistic = "iNEXT") %>% 
      dplyr::mutate(Chao_increasePercent = iChao_increasePercent)  %>% 
      dplyr::select(!iChao_increasePercent) %>% 
      setNames(colnames(.) %>% stringr::str_remove(., "iChao_|iNEXT_"))
  ) %>%
    # Remove countries where only one statistic was estimated
  dplyr::filter(!name %in% c("Kuwait", "Bahrain", "São Tomé", "El Salvador"))%>% 
  dplyr::filter(!stringr::str_detect(continent, "Seven seas")) %>% 
  dplyr::mutate(colourBar = -10)
  
# Set up colour scale for continents
myColors <- c("#F9766D", "#B79E01", "#14BA39", "#02C0C4", "#609CFF", "#F564E3")
names(myColors) <- levels(combined_explanatory_longer_increase$continent %>% as.factor)
colScale <- ggplot2::scale_colour_manual(name = "grp",values = myColors)

# plot the top half of countries 
(gapBoxplot_TOP_perc <- ggplot2::ggplot(combined_explanatory_longer %>% #dplyr::filter(level == "Country") %>%
                                         dplyr::ungroup() %>% 
                                         dplyr::arrange(.by_group = FALSE, Chao_increasePercent, name) %>%
                                         dplyr::mutate(name = factor(name, 
                                                                     levels = rev(unique(.$name)))) %>%
                                         # Take the top half of the data
                                         dplyr::slice_tail(n= ((nrow(.))/2)),
                                       aes(x = name, y = increasePercent)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = increasePercent), #width = 0.5,
                     position = position_dodge(0.90),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT"),
                              values = c("iNEXT" = "#FD9B63", "iChao" = "#55AD9B")) +
   # Add error bars 
   ggplot2::geom_errorbar(aes(ymin = lowerPercent, ymax = upperPercent, group = statistic), 
                          position = position_dodge(0.90), width = 0.2, linewidth = 0.5,
                          colour = "grey20"
                          #col = c("black", "black")
   ) +
    # Add in a new colour scale for the continent colour bars along the x-axis
    ggnewscale::new_scale_fill() +
    ggplot2::geom_bar(aes(fill = continent, y = colourBar),  stat = "identity") +
    ggplot2::scale_fill_manual(values = myColors) +
    ggplot2::geom_hline(yintercept=c(100,50), linetype="dashed", 
                        color = "black", linewidth=0.5) +
   ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust=1),
                   panel.border = element_blank()) +
   #ggplot2::ylim(c(0, 500)) +   
    scale_y_continuous(limits = c(-10, 500), oob = scales::squish) +
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Increase percentage"))
)

# plot the bottom half
# plot the top half of countries 
(gapBoxplot_BOTTOM_perc <- ggplot2::ggplot(combined_explanatory_longer %>% #dplyr::filter(level == "Country") %>%
                                     dplyr::ungroup() %>% 
                                     dplyr::arrange(.by_group = FALSE, Chao_increasePercent, name) %>%
                                     dplyr::mutate(name = factor(name, 
                                                                 levels = rev(unique(.$name)))) %>%
                                     # Take the top half of the data
                                     dplyr::slice_head(n= ((nrow(.))/2)),
                                   aes(x = name, y = increasePercent)) + 
    ggplot2::geom_bar(aes(fill = statistic, y = increasePercent), #width = 0.5,
                      position = position_dodge(0.90),  stat = "identity") +
    ggplot2::scale_fill_manual(name = "Statistic",
                               labels = c("iChao", "iNEXT"),
                               values = c("iNEXT" = "#FD9B63", "iChao" = "#55AD9B")) +
    # Add error bars 
    ggplot2::geom_errorbar(aes(ymin = lowerPercent, ymax = upperPercent, group = statistic), 
                           position = position_dodge(0.90), width = 0.2, linewidth = 0.5,
                           colour = "grey20"
                           #col = c("black", "black")
    ) +
    ggplot2::geom_hline(yintercept=c(100,50), linetype="dashed", 
                        color = "black", linewidth=0.5) +
    # Add in a new colour scale for the continent colour bars along the x-axis
    ggnewscale::new_scale_fill() +
    ggplot2::geom_bar(aes(fill = continent, y = colourBar),  stat = "identity") +
    ggplot2::scale_fill_manual(values = myColors) +
    ggplot2::geom_hline(yintercept=c(100,50), linetype="dashed", 
                        color = "black", linewidth=0.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position=c(.9,.75),
                   axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust=1),
                   panel.border = element_blank()) +
    #ggplot2::ylim(c(0, 500)) +   
    scale_y_continuous(limits = c(-10, 500), oob = scales::squish) +
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Increase percentage"))
)


###### b. plot increase ####

combined_explanatory_longer_increase <- dplyr::bind_rows(
  # ichao
  combined_explanatory %>% dplyr::select(c("name", "observedRichness", "continent",
                                           tidyselect::contains("iChao"))) %>%
    dplyr::mutate(statistic = "iChao") %>%
    dplyr::mutate(Chao_increase = iChao_increase) %>% 
    setNames(colnames(.) %>% stringr::str_remove(., "iChao_|iNEXT_")) %>%
    dplyr::mutate(lower = lower - observedRichness,
                  lower = dplyr::if_else(lower < 0, 000.1, lower),
                  upper = upper - observedRichness),
  # iNEXT
  combined_explanatory %>% dplyr::select(c("name", "observedRichness", "iChao_increase",
                                           "continent",
                                           tidyselect::contains("iNEXT"))) %>%
    dplyr::select(!"niNEXT") %>% 
    dplyr::mutate(statistic = "iNEXT") %>% 
    dplyr::mutate(Chao_increase = iChao_increase)  %>% 
    dplyr::select(!iChao_increase) %>% 
    setNames(colnames(.) %>% stringr::str_remove(., "iChao_|iNEXT_"))%>%
    dplyr::mutate(lower = lower - observedRichness,
                  lower = dplyr::if_else(lower < 0, 000.1, lower),
                  upper = upper - observedRichness)
) %>%
  # Remove countries where only one statistic was estimated
  dplyr::filter(!name %in% c("Kuwait", "Bahrain", "São Tomé", "El Salvador")) %>% 
  dplyr::filter(!stringr::str_detect(continent, "Seven seas")) %>% 
  dplyr::mutate(colourBar = -15)

# plot the top half of countries 
(gapBoxplot_TOP_in <- ggplot2::ggplot(combined_explanatory_longer_increase %>% #dplyr::filter(level == "Country") %>%
                                          dplyr::ungroup() %>% 
                                          dplyr::arrange(.by_group = FALSE, Chao_increase, name) %>%
                                          dplyr::mutate(name = factor(name, 
                                                                      levels = rev(unique(.$name)))) %>%
                                          # Take the top half of the data
                                          dplyr::slice_tail(n= ((nrow(.))/2)),
                                        aes(x = name, y = increase)) + 
   ggplot2::geom_bar(aes(fill = statistic, y = increase), width = 0.8,
                     position = position_dodge(0.8),  stat = "identity") +
   ggplot2::scale_fill_manual(name = "Statistic",
                              labels = c("iChao", "iNEXT"),
                              values = c("iNEXT" = "#FD9B63", "iChao" = "#55AD9B")) +
   # Add error bars 
   ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                          position = position_dodge(0.90), width = 0.2, linewidth = 0.5,
                          colour = "grey20"
                          #col = c("black", "black")
   ) +
    
    # Add in a new colour scale for the continent colour bars along the x-axis
    ggnewscale::new_scale_fill() +
    ggplot2::geom_bar(aes(fill = continent, y = colourBar),  stat = "identity") +
    ggplot2::scale_fill_manual(values = myColors) + 
  theme_classic() +
    ggplot2::theme(legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust=1),
                   panel.border = element_blank()) +
    ggplot2::scale_y_continuous(limits = c(-15,1000), expand = c(0, 0)) +  
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Increased richness")) #+
    #annotation_custom(cowplot::ggdraw(cowplot::get_legend(barLegend_gap)) %>%
    #                    ggplot2::ggplotGrob(), xmin = 170, xmax = 3, 
    #                  ymin = 800, ymax = 1000)
)

# plot the bottom half
# plot the top half of countries 
(gapBoxplot_BOTTOM_in <- ggplot2::ggplot(combined_explanatory_longer_increase %>% #dplyr::filter(level == "Country") %>%
                                             dplyr::ungroup() %>% 
                                             dplyr::arrange(.by_group = FALSE, Chao_increase, name) %>%
                                             dplyr::mutate(name = factor(name, 
                                                                         levels = rev(unique(.$name)))) %>%
                                             # Take the top half of the data
                                             dplyr::slice_head(n= ((nrow(.))/2)),
                                           aes(x = name, y = increase)) + 
    ggplot2::geom_bar(aes(fill = statistic, y = increase), width = 0.8,
                      position = position_dodge(0.80),  stat = "identity") +
    ggplot2::scale_fill_manual(name = "Statistic",
                               labels = c("iChao", "iNEXT"),
                               values = c("iNEXT" = "#FD9B63", "iChao" = "#55AD9B")) +
    # Add error bars 
    ggplot2::geom_errorbar(aes(ymin = lower, ymax = upper, group = statistic), 
                           position = position_dodge(0.90), width = 0.2, linewidth = 0.5,
                           colour = "grey20"
                           #col = c("black", "black")
    ) +
      # Add in a new colour scale for the continent colour bars along the x-axis
    ggnewscale::new_scale_fill() +
    ggplot2::geom_bar(aes(fill = continent, y = colourBar),  stat = "identity") +
    ggplot2::scale_fill_manual(name = "Continent",
                               values = myColors) + 
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position=c(.9,.75),
                   axis.text.x = ggplot2::element_text(angle = 60, vjust = 1, hjust=1),
                   panel.border = element_blank()) +
    ggplot2::scale_y_continuous(limits = c(-15,1000), expand = c(0, 0)) +   
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Increased richness"))
)


      ###### c. combine ####
# Combine INCREASE
(gapPlots_in <- cowplot::plot_grid(gapBoxplot_TOP_in, gapBoxplot_BOTTOM_in,
                                     labels = c("(A)","(B)"),
                                     ncol = 1, align = 'v', axis = 'l', 
                            label_y = 1.005,
                            label_x = 0.04)
)

# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.5c_gapPlots_in.pdf"),
                   plot = gapPlots_in,
                   base_width = 12,
                   base_height = 10)

# Combine PERCENTAGE
(gapPlots_per <- cowplot::plot_grid(
                                gapBoxplot_TOP_perc, gapBoxplot_BOTTOM_perc, 
                                labels = c("(A)","(B)"),
                                ncol = 1, align = 'v', axis = 'l', 
                                label_y = 1.005,
                                label_x = 0.04)
)

# Save the plot
cowplot::save_plot(filename = paste0("Figure_outputs/","5.5c_gapPlots_per.pdf"),
                   plot = gapPlots_per,
                   base_width = 12,
                   base_height = 10)


#### 6.0 Compare curves ####
    ##### 6.1 read data ####
  # read in the three datasets produced from the three different curves
litCurve_combine <- readr::read_csv("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Table_outputs/3.1c_combinedStatistics_longer.csv") %>%
  dplyr::filter(scale == "Country") %>% 
  dplyr::select(name, iChao_est, iChao_lower, iChao_upper, observedRichness) %>% 
  dplyr::mutate(source = "Literature") %>%
  dplyr::arrange(dplyr::desc(iChao_est))
occCurveGlobal_combine <- readr::read_csv("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/OccCurve/Table_outputs/3.1c_combinedStatistics_longer.csv") %>%
  dplyr::select(name, iChao_est, iChao_lower, iChao_upper, iChao_increasePercent, iChao_increase) %>% 
  dplyr::mutate(source = "Global") 
occCurveCountry_combine <- readr::read_csv("/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/OccCurve_spCou/Table_outputs/3.1c_combinedStatistics_longer.csv") %>%
  dplyr::select(name, iChao_est, iChao_lower, iChao_upper, iChao_increasePercent, iChao_increase) %>% 
  dplyr::mutate(source = "Country") 

  # Combine the curve data together
combinedCurveData <- litCurve_combine %>%
  dplyr::bind_rows(occCurveGlobal_combine, occCurveCountry_combine) %>% 
  dplyr::group_by(source) %>% dplyr::arrange(dplyr::desc(iChao_est)) %>%
  # Take the top 20 countries from the literature curve
  dplyr::filter(name %in% (litCurve_combine %>%
                             dplyr::slice_head(n = 20) %>%
                             dplyr::pull(name)))

  ##### 6.2 Plot ####
# Plot the outputs for the top 20 countries from each curve
(curveBarPlot <- ggplot2::ggplot(combinedCurveData %>% 
                                   dplyr::arrange(.by_group = FALSE, iChao_est) %>%
                                   dplyr::mutate(name = factor(name, 
                                                               levels = rev(unique(.$name)))),
                                 aes(x = name, y = iChao_est, fill = source)) + 
    ggplot2::geom_bar(aes(fill = source, y = iChao_est), #width = 0.5,
                      position = position_dodge(0.90),  stat = "identity") +
    ggplot2::scale_fill_manual(values = c("#FFC125", "#1BB6AFFF", "#172869FF"),
                               label = c("Literature",
                                         "Occurrences global",
                                         "Occurrences country"),
                               name = "Data source")  +
    # Add in the observed richness
    ggplot2::geom_bar(aes(fill = source, y = observedRichness),
                      position = position_dodge(0.90),  stat = "identity", fill = "grey", 
                      alpha = 0.5) +
    # Add error bars 
    ggplot2::geom_errorbar(aes(ymin = iChao_lower, ymax = iChao_upper, group = source), 
                           position = position_dodge(0.90), width = 0.2, linewidth = 1,
                           colour = "grey50"
    ) +
    
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position.inside = c(0.76, 0.86),
                   legend.position = "inside",
                   axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
    #ggplot2::ylim(c(20000, 27000)) +
    ggplot2::xlab(c( "")) + ggplot2::ylab(c("Species"))
)

# Save the plot
ggplot2::ggsave(paste0("Figure_outputs/","6.0_curveBarPlot.pdf"), 
                plot = curveBarPlot, 
                width = 6, height = 6, units = "in", dpi = 300)


  ##### 6.3 lmer ####
    ###### a. global #####
occCurveGlobal_combine <- occCurveGlobal_combine %>%
  countryHarmoniseR(countryColumn = "name") %>% 
  dplyr::left_join(., combined_explanatory %>% dplyr::select(!c("iChao_est", "iChao_lower",
                                                             "iChao_upper", "iChao_increasePercent",
                                                             "iChao_increase")),
                                                          by = "name")

# Make a lmer for the increase in the estimated number of species
lmerTest::lmer(data = occCurveGlobal_combine %>%
                 dplyr::filter(iChao_increase > 0),
                                formula = log(iChao_increase) ~ 
                                  # Fixed effects
                                  log(GPD_per_cap) + log(percent_ed) + 
                                  log(median_roadDist) +
                                  log(elevationalRange)  + 
                                  log(area_m) +
                                  log(observedRichness)  + 
                                  # Removing cleanRecords below  will change the sign of propOccurrences to negative
                                  log(cleanRecords) *
                                  log(propOccurrences) +
                                  # Random effect
                                  (1|continent),
                                REML = TRUE) %>% 
  # Get the output summary
  summary(.)

# Make a lmer for the increase in the estimated PERCENTAGE of species
lmerTest::lmer(data = occCurveGlobal_combine %>%
                 dplyr::filter(iChao_increasePercent > 0),
               formula = log(iChao_increasePercent) ~ 
                 # Fixed effects
                 log(GPD_per_cap) + log(percent_ed) + 
                 log(median_roadDist) +
                 log(elevationalRange)  + 
                 log(area_m) +
                 log(observedRichness)  + 
                 # Removing cleanRecords below  will change the sign of propOccurrences to negative
                 log(cleanRecords) *
                 log(propOccurrences) +
                 # Random effect
                 (1|continent),
               REML = TRUE) %>% 
  # Get the output summary
  summary(.)

###### b. country #####
occCurveCountry_combine <- occCurveCountry_combine %>%
  countryHarmoniseR(countryColumn = "name") %>% 
  dplyr::left_join(., combined_explanatory %>% dplyr::select(!c("iChao_est", "iChao_lower",
                                                                "iChao_upper", "iChao_increasePercent",
                                                                "iChao_increase")),
                   by = "name")

# Make a lmer for the increase in the estimated number of species
lmerTest::lmer(data = occCurveCountry_combine %>%
                 dplyr::filter(iChao_increase > 0),
               formula = log(iChao_increase) ~ 
                 # Fixed effects
                 log(GPD_per_cap) + log(percent_ed) + 
                 log(median_roadDist) +
                 log(elevationalRange)  + 
                 log(area_m) +
                 log(observedRichness)  + 
                 # Removing cleanRecords below  will change the sign of propOccurrences to negative
                 log(cleanRecords) *
                 log(propOccurrences) +
                 # Random effect
                 (1|continent),
               REML = TRUE) %>% 
  # Get the output summary
  summary(.)

# Make a lmer for the increase in the estimated PERCENTAGE of species
lmerTest::lmer(data = occCurveCountry_combine %>%
                 dplyr::filter(iChao_increasePercent > 0),
               formula = log(iChao_increasePercent) ~ 
                 # Fixed effects
                 log(GPD_per_cap) + log(percent_ed) + 
                 log(median_roadDist) +
                 log(elevationalRange)  + 
                 log(area_m) +
                 log(observedRichness)  + 
                 # Removing cleanRecords below  will change the sign of propOccurrences to negative
                 log(cleanRecords) *
                 log(propOccurrences) +
                 # Random effect
                 (1|continent),
               REML = TRUE) %>% 
  # Get the output summary
  summary(.)





#### TEST ####
  # Get the Fiji map from natural earth

# shift coordinates to recenter worldmap
worldmap <- ggplot2::map_data("world", wrap = c(0, 360))
# Download the Fijian DEM rasters — one for each side of the dateline and convert to EPSG:3460
# — Fiji 1986
FijiRastWest <- geodata::elevation_3s(country='FJI', 
                                      path = RootPath,
                                      mask = TRUE,
                                      lat = c( -16),
                                      lon = c(180),
                                      res = "0.5") %>%
  terra::project(., terra::crs("EPSG:3460"),
                 gdal = TRUE,   method = "near",
                 threads = TRUE, res = 90)
FijiRastEast <- geodata::elevation_3s(country='FJI', 
                                      path = RootPath,
                                      mask = TRUE,
                                      lat = c( -16),
                                      lon = c(-180),
                                      res = "0.5") %>%
  terra::project(., terra::crs("EPSG:3460"),
                 gdal = TRUE,   method = "near",
                 threads = TRUE, res = 90)
# Merge the fiji map halves
FijiMap <- terra::merge(FijiRastWest, FijiRastEast) %>%
  terra::classify(c(0,200,400,600,800,1000,1200,1400),
                  right = FALSE)

# Crop roads to Fiji dataset
FijiRoads <- roads %>%
  sf::st_transform( crs = sf::st_crs("EPSG:3460")) %>% 
  sf::st_crop(., terra::ext(FijiMap)) 


# Set up limits
yLimInput = c(3591385+200000, 4062964+30000)
xLimInput = c(1871432- 70000, 2298672 - 100000.0)
# Reproject the Fiji 1986 extent to WGS84 for the inset
WGS84extent <- sf::st_bbox(c(xmin = xLimInput[1], xmax = xLimInput[2], 
                             ymin = yLimInput[1], ymax = yLimInput[2]),
                           crs = terra::crs("EPSG:3460")) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(., crs = terra::crs("EPSG:4326")) %>%
  st_bbox()
# Extract the limits and format
WGS84extent2 <- tibble::tibble(
  point = c("xmin", "xmax", "ymin", "ymax"),
  coords = c(WGS84extent[1], WGS84extent[3],WGS84extent[2], WGS84extent[4]) %>%
    as.numeric()) %>%
  dplyr::mutate(coords = dplyr::if_else(coords < 0 & 
                                          stringr::str_detect(point, "^x"),
                                        coords + 360,
                                        coords)  )





