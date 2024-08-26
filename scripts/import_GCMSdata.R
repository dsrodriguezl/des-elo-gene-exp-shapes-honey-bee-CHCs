tic("Script to import GCMS integration data")

# Packages ----
# Install devtools package if it is not yet installed
install_new_packs("devtools")

# Install analyzeGC v0.2.1 if the package is not yet installed
if ("analyzeGC" %in% 
    rownames(installed.packages()) == FALSE) {
  devtools::install_github("dsrodriguezl/analyzeGC"
                           , dependencies = TRUE
                           , ref = "v0.2.1")
}

# Vector holding the list of packages that will be used by the script
script_packs <- c("analyzeGC", 
                  "GCalignR", 
                  "readxl")

# install and/or load packages
load_my_packs(script_packs)

# Import data files ----
## samples ####
# Vector defining what string patterns to remove form the CSV file names to
# obtain the individual sample label
patterns_2_delete <- paste("DR-PL_",
                           "DR_PL_",
                           "_c",
                           paste(paste0("_",
                                        1:20) |> 
                                   rev()
                                 , collapse = "|"),
                           ".CSV",
                           ".csv"
                           , sep = "|")

## Import the GCMS integration data for all the samples of the batch
samples_data_list <- 
  # List the CSV files holding the inegration data of the batch
  list.files(path = here("data", "raw"
                         , "gcms-integration"
                         , gcms_batch
                         )
             # Get all CSV files in the folder
             , pattern = ".CSV" 
             , full.names = T
             , recursive = T) |>
  # Do not include standards
  str_subset('STD', negate = T) |> 
  # Import the integration data
  import_mh_data(patterns_2_delete = patterns_2_delete)

samples_data_list |> summary()

## Standards (STD) ####
patterns_2_delete <- paste("STD",
                           "_SSL_",
                           ".CSV",
                           ".csv"
                           , sep = "|")

## Import the GCMS integration data for all the standards of the batch
standards_list <- list.files(path = here("data", "raw"
                                         , "gcms-integration"
                                         , gcms_batch
                                         )
                             # Get all CSV files in the folder
                             , pattern = ".CSV" 
                             , full.names = T
                             , recursive = T) |>
  # Do not include samples
  str_subset('STD') |> 
  import_mh_data(patterns_2_delete = patterns_2_delete)

standards_list |> summary()

# group data frames by sampling groups within a list ----
## Load group membership information
load(here("data", "processed", "grouping_info.Rdata"))

# Add group label column combining desired factors to analyse
samples_info <- grouping_info |> 
  # Group label is compound by all the factors
  unite(group_label, where(is.factor), sep = "_", remove = F) |> 
  relocate(group_label, .after = last_col())

## Group GCMS data files based on unique group labels
samples_data_list <- mg_list(sample.info = samples_info
                             , group.label = "group_label"
                             , samples.data.list = samples_data_list)

# Input check-up for GCalignR ----
pdf(here("output"
         , paste0("GCalignR_input-check_"
                  , gcms_batch
                  , '.pdf'))
    , width = 10, height = 5)

lapply(samples_data_list
       , check_input
       , plot = T)

lapply(samples_data_list
       , peak_interspace
       , rt_col_name = "RT"
       , quantile_range = c(0, 0.8)
       , quantiles = 0.05)


# Close graphic device to export plots into the PDF file
dev.off()
print("Input check-up plots were exported")

# Export data ----
folder <- here("data", "tmp"
               , "data2align")
if (file.exists(folder)) {
  cat("The folder already exists")
} else {
  dir.create(folder)
}
save(list = c("samples_data_list", "standards_list")
     , file = here(folder
                   , paste0(gcms_batch
                            , ".Rdata")))

print("The data list(s) to align have been exported to the tmp data folder")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "SInf_import-gcms-script.txt"))

## Detach/unload packages
pacman::p_unload(char = script_packs)

## Clear environment
objects_list <- ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|"))
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

gc()

toc()


