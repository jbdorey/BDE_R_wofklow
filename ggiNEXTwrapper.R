# This function was written starting on the 24th of June 2024 by James B Dorey in order
# to wrap iNEXT::ggiNEXT to work with iNEXTwrapper and plot multiple accumulation curves
# at once across pages.

#' ggplot2 extension for an iNEXT object
#' 
#' `BeeBDC::ggiNEXTwrapper()` wraps `iNEXT::ggiNEXT()` and extends its functionality to work with 
#' `BeeBDC::iNEXTwrapper()` and produce plots for multiple groups in one go. The plots can be 
#' multiple per page and across multiple pages.
#' 
#' @param data a list of \code{iNEXT} objects computed by `BeeBDC::iNEXTwrapper()`.
#' @param type three types of plots: sample-size-based rarefaction/extrapolation curve (\code{type = 1}); 
#' sample completeness curve (\code{type = 2}); coverage-based rarefaction/extrapolation curve (\code{type = 3}).            
#' @param se a logical variable to display confidence interval around the estimated sampling curve.
#' @param facet.var create a separate plot for each value of a specified variable: 
#'  no separation \cr (\code{facet.var="None"}); 
#'  a separate plot for each diversity order (\code{facet.var="Order.q"}); 
#'  a separate plot for each assemblage (\code{facet.var="Assemblage"}); 
#'  a separate plot for each combination of order x assemblage (\code{facet.var="Both"}).              
#' @param color.var create curves in different colors for values of a specified variable:
#'  all curves are in the same color (\code{color.var="None"}); 
#'  use different colors for diversity orders (\code{color.var="Order.q"}); 
#'  use different colors for sites (\code{color.var="Assemblage"}); 
#'  use different colors for combinations of order x assemblage (\code{color.var="Both"}).  
#' @param grey a logical variable to display grey and white ggplot2 theme. 
#' @param ... other arguments passed on to methods. Not currently used.
#' @return a ggplot2 object
#' 

ggiNEXTwrapper <- function(
    data = country_iNEXT,
    filterOut = NULL,
    type = 1,
    se = TRUE,
    facet.var = "None",
    color.var = "Order.q",
    grey = FALSE,
    legendPerPlot = FALSE,
    nrow = 3,
    ncol = 4,
    labels = NULL,
    fileName = "iNEXTplots",
    outPath = "/Users/jamesdorey/Desktop/Uni/My_papers/BeeDiversityEstimates/BDE_R_wofklow/Figure_outputs/country_iNEXT",
    base_width = 8.3,
    base_height = 11.7, 
    dpi = 300,
    ...){
  
  decimalLatitude <- decimalLongitude <- database_id <- scientificName <- NULL
  
  requireNamespace("magrittr")
  requireNamespace("iNEXT")
  requireNamespace("cowplot")
  requireNamespace("ggplot2")
  
  
  #### 0.0 Prepare function ####
    ##### 0.1 Fatal warnings ####
  if (is.null(data)) {
    stop("No data provided. Please provide some data produced from BeeBDC::ggiNEXTwrapper().")
  }
  
    ##### 0.2 Prepare values ####
    # If no labels were provided, take them as the length of nrow and ncol a:z
  if(is.null(labels)){
    labels = letters[1:(nrow*ncol)]
  }
  
  #### 1.0 Prepare data ####
    ##### 1.0 Extract listed countries ####
  dataExtracted <- data$iNextEst$iNextEst
  
    # Filter out levels in filterOut
  if(!is.null(filterOut)){
    dataExtracted <- dataExtracted[!names(dataExtracted) %in% filterOut]
  }
  
  #### 1.0 Make plots ####
    ##### 1.1 Build plots ####
  ggiNEXT_fun <- function(dataIn = dataExtracted,
                      type = type,
                      se = se,
                      facet.var = facet.var,
                      color.var = color.var,
                      grey = grey){
      # extract the iNEXT data itself
    dataWithin <- dataIn$iNextEst
    iNEXT_plot <- iNEXT::ggiNEXT(x = dataWithin,
                                 type = type,
                                 se = se,
                                 facet.var = facet.var,
                                 color.var = color.var,
                                 grey = grey) +
      ggplot2::theme_classic() +
      ggplot2::ggtitle(paste0(dataIn$AsyEst$groupVariable %>% unique(),"\n",
                              "n = ", format(dataIn$DataInfo$n, big.mark = ","),"\n",
                              "Obs. = ", format(dataIn$AsyEst$Observed %>% round(0),
                                                    big.mark = ","),
                              "; Est. = ", format(dataIn$AsyEst$Estimator %>% round(0),
                                                     big.mark = ","),"\n",
                              "Lower = ", format(dataIn$AsyEst$`95% Lower` %>% round(0),
                                                 big.mark = ","),
                              "; Upper = ", format(dataIn$AsyEst$`95% Upper` %>% round(0),
                                                 big.mark = ",")
                              ))
    
      # Remove the legend from each plot
    if(legendPerPlot == FALSE){
      iNEXT_plot <- iNEXT_plot + 
        ggplot2::theme(legend.position="none")
    }
    
    # Return the plot 
    return(iNEXT_plot)
  } # END ggiNEXT_fun


  # Make the plots per country using lapply
  countryPlots <- dataExtracted %>%
    lapply(X = .,
           FUN = ggiNEXT_fun,
           type = type,
           se = se,
           facet.var = facet.var,
           color.var = color.var,
           grey = grey)

  
  
  
  #### 2.0 Combine plots ####
    ##### 2.1 Chunk plot list ####
    # Get the number of chunks to be output
  numberOfChunks <- ceiling(length(countryPlots)/(nrow*ncol) )
  plotsPerPage <- (nrow*ncol)
  
  countryPlot_chunks <- dplyr::lst()
  j = 1
    # Put each chunk of plots into its own list
  for(i in 1:numberOfChunks){
    countryPlot_chunk_i <- countryPlots[j:(plotsPerPage + j)]
    countryPlot_chunks <- append(countryPlot_chunks, list(countryPlot_chunk_i))
    j = j + plotsPerPage 
  }
  
  
  ##### 2.2 Combine plots ####
    # A fucntion to do the cowplots
  multiPlotFun <- function(plotData = countryPlot_chunks,
                           nrow = nrow,
                           ncol = ncol,
                           labels = labels,
                           base_width = base_width,
                           base_height = base_height, 
                           dpi = dpi){
      # Make the plot grid
    plotOut <- do.call("plot_grid",
                       args = list(plotlist = plotData,
                                   nrow = nrow,
                                   ncol = ncol,
                                   labels = labels[1:length(plotData[lengths(plotData) > 0])]))
    
    # Save the plots
    cowplot::save_plot(
      plot = plotOut,
      filename = paste(outPath, "/", fileName, "_", names(plotData)[[1]], "_to_",
                       names(plotData)[[length(plotData)]], ".pdf", sep = ""),
      base_width = base_height,
      base_height = base_height, 
      dpi = dpi)
    
    # Return the plot grid
    return(plotOut)
  }
  
  
    # Combine the plots with cowplot
 combinedPlots <- lapply(X = countryPlot_chunks,
                         FUN = multiPlotFun,
                         nrow = nrow,
                         ncol = ncol,
                         labels = labels,
                         base_width = base_height,
                         base_height = base_height, 
                         dpi = dpi )

} # END ggiNEXTwrapper





