# This function was written by James B. Dorey starting on the 19th of April 2023 to wrap iNEXT
# and output table for multiple species/countries/variables at once

iNEXTwrapper <- function(data = NULL,
                         variableColumn = "groupVariable",
                         valueColumn = "n",
                         q = 0,
                         datatype = "abundance",
                         conf = 0.95,
                         se = TRUE,
                         nboot = 50,
                         size = NULL,
                         endpoint = NULL,
                         knots = 40,
                         mc.cores = 1){
  
  
  #### 0.0 Prep ####
  ##### 0.1 Functions ####
  # Duplicate a function within ChaoSpecies to not run data-poor species
  #f <- function(i, data) {
  #  length(data[which(data == i)])
  #}
  
  # wrap the iNEXT function
  wrapper <- function(inputData_i = df_list[[1]]){
      ###### a. data prep ####
      # Format the data coming in
    inputData_i <- inputData_i[[1]] %>%
      dplyr::select(groupVariable, groupCount) %>%  # Create the rownames
      dplyr::tibble()
    
     ###### b. iNEXT ####
    # Don't run species where there are no counts that are less than k
    # i.e. ONLY run species that won't throw an error due to weird/poor cases of data
    # if(any(inputData_i[[1]] > 0 & inputData_i[[1]] < 5) & 
    #    !(f(1, inputData_i[[1]]) == sum(inputData_i[[1]], na.rm = TRUE))
    # ){
      # Run the ChaoSpecies function
      spOutput <- iNEXT::iNEXT(
        x = inputData_i$groupCount,
        datatype = datatype,
        q = q,
        conf = conf,
        se = se,
        nboot = nboot,
        size = size,
        endpoint = endpoint,
        knots = knots)

        ###### c. DataInfo ####
        # Extract the dataInfo and then replace the Assemblage with the groupVariable value
      dataInfo_out <- spOutput$DataInfo %>%
        dplyr::mutate(Assemblage = inputData_i$groupVariable %>% unique())
      
      ###### d. AsyEst ####
        # Save AsyEst and replace the rownames with a column 
      AsyEst_out <- spOutput$AsyEst %>%
        dplyr::mutate(statistic = rownames(.), .before = 1) %>% 
        dplyr::tibble() %>%
          # Add the grouping variable to the start
        dplyr::mutate(groupVariable = inputData_i$groupVariable %>% unique(),
                      .before = 1)
        
      
      ###### e. iNextEst ####
    # Save the iNextEst as a list that's named per groupVariable
      iNextEst <- spOutput #%>%
        #dplyr::lst() %>%
        #stats::setNames(inputData_i$groupVariable %>% unique())
      
      
      # Return the data as a list
      return(list( dataInfo_out, AsyEst_out, iNextEst) %>%
               setNames(c( "DataInfo", "AsyEst", "iNextEst")))
    # }
  }
  
  #### 1.0 Data prep ####
    ##### 1.1 Prep tibble ####
    # Prep the correct column
  # Temporarily rename the variableColumn to "groupVariable" within the function
  data <- data %>%
    dplyr::rename("groupVariable" = tidyselect::any_of(variableColumn))
  # Temporarily rename the variableColumn to "groupVariable" within the function
  data <- data %>%
    dplyr::rename("groupCount" = tidyselect::any_of(valueColumn))
  
    # Turn the data into a list by the country
  data_list <- data %>%
    tidyr::drop_na(groupVariable) %>% 
    dplyr::group_by(groupVariable) %>%
    dplyr::group_split()
    # Create an empty list to add to
  df_list <- dplyr::lst()
  variableNames <- c()
  
    ##### 1.2 Loop prep ####
  # Loop the data to make it into a list (it's a two-level list)
  for(i in 1:length(data_list)){
    # Extract the list and give it the country name
    loopList <- data_list[[i]] %>%
      dplyr::lst() %>%
      setNames(data_list[[i]]$groupVariable[[1]])
    # Add to the variable list to the greater list
    df_list <- df_list %>%
      append( dplyr::lst(loopList))
    # Add the variable name to a running list to ensure name order is maintained
    variableNames <- c(variableNames, data_list[[i]]$groupVariable[[1]])
  }
  # Set the names for the list
  df_list <- df_list %>%
    stats::setNames(variableNames)
  
 
  
  
  #### 2.0 Run functions ####    
  # Run the function per species or level
  iNEXTOutput <- parallel::mclapply(
    X = df_list,
    FUN = wrapper,
    mc.cores = mc.cores
  ) 
  
  
  #### 3.0 Process outputs ####
  ##### 3.1 seperate F + P ####
  # Find the non-null variables and extract those
  non_empty_list_test <- !sapply(iNEXTOutput <- iNEXTOutput, is.null)
  # Save a list of variables that could not be analysed
  failures <- iNEXTOutput[!non_empty_list_test] %>%
    # Get the failure names 
    names()
  # Remove the failed list
  iNEXTOutput <- iNEXTOutput[non_empty_list_test]
  
  
  ##### 3.2 dataInfo outputs ####
  # Now, combine each level of the list across variables
  # Extract the diversity table stuff
  dataInfoOut <- lapply(iNEXTOutput, function(x) x[["DataInfo"]]) %>%
    # Bind together with the original two columns
    dplyr::bind_rows()
  
  
  ##### 3.3 AsyEst outputs ####
  # Row bind the diversity statistics into a single table
  iNEXTOutputOut <- lapply(iNEXTOutput, function(x) x[["AsyEst"]]) %>%
    dplyr::bind_rows()
  # Return the variableColumn name to its original state
  names(iNEXTOutputOut)[names(iNEXTOutputOut) == "groupVariable"] <- variableColumn
  
  ##### 3.4 iNEXT outputs ####
    # Combine the iNextEsts ino a single list
  iNextEstOut <- iNEXTOutput %>%
    lapply(., function(x) x) %>%
    dplyr::lst() %>%
    stats::setNames("iNextEst")
  
  
  ##### 3.5 Combine ####
  output <- dplyr::lst(dataInfoOut, iNEXTOutputOut, iNextEstOut) %>% 
    setNames(c("DataInfo", "AsyEst", "iNextEst"))
  
  #### 4.0 User output ####
  # provide some user output
  if(length(failures > 0)){
    writeLines(paste0(
      " - We could not examine the following variables (because of insufficent data or sample size): ",
      paste(failures, collapse = ", ")
    ))}
  
  writeLines(paste0(" - Outputs can be found in a list with two tibbles called 'DataInfo' and",
                    " 'AsyEst' and a list of iNext outputs per groupVariable in iNextEst'."))

  
  # Return the output
  return(output)
  
} # END function

  

  
  