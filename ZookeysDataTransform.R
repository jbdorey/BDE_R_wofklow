


Indata <- readxl::read_excel("/Users/jamesdorey/Downloads/zookeys-755-001-s001.xlsx") %>%
  dplyr::mutate(Total = as.numeric(Total)) %>%
  dplyr::mutate(Total = dplyr::if_else(is.na(Total), 1, Total))

outData <- Indata %>%
  tidyr::uncount(Total) %>%
  dplyr::mutate(sex = dplyr::if_else(Females>0 & Males == 0,"F",
                                     dplyr::if_else(Males > 0 & Females == 0,
                                                    "M", NA_character_))) %>%
    # Save the output
  readr::write_excel_csv("/Users/jamesdorey/Downloads/Outfile.csv")
  

