tic("Group tables script")

# Packages ----
# Vector holding the list of packages that will be used
script_packs <- c("openxlsx"
                  , "purrr"
                  , "analyzeGC")

# install and/or load packages
load_my_packs(script_packs)

# Import data files ----
## Load aligned GC-MS data frames
load(here("data", "tmp"
          , paste0(gcms_batch
                   , "_aligned_gcms-data.Rdata")))

## Load the compound information, result from the identification process
## This was manually added to the "identification" page of the XLSX files

# List of file paths for the XLSX with the corrected alignments
path_comps_id <- list.files(path = here("data", "tmp")
                            , pattern = gcms_batch
                            , full.names = T) |> 
  str_subset("table.xlsx") |> 
  str_subset("master_", negate = T) |> 
  str_subset("batch_", negate = T)

# Import the compound identification form the "identification" page of the XLSX
# files
comps_id <- lapply(path_comps_id
                   , read.xlsx 
                   , sheet = "identification"
                   ) |> 
  lapply(as_tibble)
summary(comps_id)

# Name the entries of the list as the group the belong to
names(comps_id) <- str_split(path_comps_id
                                 , "/"
                                 , simplify = T) |> 
  str_subset(".xlsx") |> 
  str_remove_all(paste("_table.xlsx"
                       , "_RT"
                       , paste0(gcms_batch, "_")
                       , sep = "|"))
summary(comps_id)
comps_id

## Load group membership information
load(here("data", "processed", "grouping_info.Rdata"))
grouping_info

# Shape the data frames ----
## Standards ####
### Export PFD with plots of the abundance of standards vs mean RT
### They visualize the  differences in performance of the total ion counts
### of hydrocarbons along the run for a given batch
pdf(here("output"
         , paste0("standards-area-v-RT_"
                  , gcms_batch
                  , '.pdf'))
    , width = 18, height = 8)
# shape the data frame with the compound information of the standards
std_info <- shape_hcstd_info(comps_id.STD = comps_id$STD
                             , aligned_std = aligned_STD_data_list
                             , short_std_pattern = "L"
                             , long_std_pattern = "H")
dev.off()

## Compounds information ####
### Extract the compounds information for the groups
comps_info <- comps_id |>
  # Iterate through groups
  lapply(get_hc_info)

comps_info

## Samples ####

### Adjust the abundance (area) of the peaks
### given the observed  differences in performance of the total ion counts
### in the standards
#### Get plots to compare before and after correction
pdf(here("output"
         , paste0("samples_correction_plots_"
                  , gcms_batch
                  , '.pdf'))
    , width = 18, height = 8)

adjusted_samples_list <- corrected_samples_list2 |>
  lapply(adjust_abundance, std.info = std_info)

dev.off()

### Add CHC identification to the aligned sample data frames
unfiltered_samples_list <- add_comps_info(samples.list = adjusted_samples_list
                                          ,comps.info.list = comps_info)

### grouping_info ####
### Define grouping variables as factors
grouping_info[, colnames(grouping_info)] <-
  lapply(grouping_info[, colnames(grouping_info)], as.factor)
grouping_info

# Clean the data ----
## Delete trace compounds ####
filtered_samples_list <- unfiltered_samples_list |>
  # iterate through groups
  lapply(trace_comps
         # Remove from each sample any compound that represents less than 0.01%
         # of the total ion count of the corresponding sample
         , threshold = 0.01)

## Delete non-CHC compounds ####
### Samples 
filtered_samples_list2 <- filtered_samples_list |>
  lapply(drop_na_compounds)

# Kováts Retention index RI) ----
samples_plus_ri_list <- filtered_samples_list2 |>
  #  Iterate through groups, calculating the RIs of the compounds
  lapply(kovats_retention_index, std.info = std_info)

# Export final data frames ----
## Prepare the final group tables
group_tables <- samples_plus_ri_list |>
  lapply(shape_group_table)

## Export the data frames
save(list = c("group_tables"
              , "grouping_info")
     , file = here("data"
                   , "tmp"
                   , paste0(gcms_batch, "_group_tables.Rdata")))
print("The data frames for analysis were exported")

## Export the RT data frames to the XSLX files
## This eases the evaluation of the group tables
for (df in names(group_tables)) {
  group_table_file <- here("data", "tmp"
                           , paste0(gcms_batch, "_", df, "_table.xlsx"))
  
  wb <- loadWorkbook(group_table_file)
  if (sheets(wb) |> str_detect("group_table") |> any()) {
    wb |> 
      removeWorksheet(sheet = "group_table")
  }
  wb |> 
    addWorksheet(sheetName = "group_table")
  wb |> 
    conditionalFormatting(sheet = "group_table"
                          , cols = 1
                          , rows = 1:nrow(group_tables |> 
                                            pluck(df) |> 
                                            pluck("RT"))
                          , type = "duplicates")
  wb |> 
    writeData(x = group_tables |> 
                pluck(df) |> 
                pluck("Area")
              , sheet = "group_table")
  wb |> 
    saveWorkbook(file = group_table_file
                 , overwrite = T)
  print(paste0(gcms_batch, "_", df, "_table.xlsx", " was modified"))
  rm(df)
}

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_group-samples-script.txt"))


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

