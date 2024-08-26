tic("Master table")

# Packages ----
# Vector holding the list of packages that will be used
script_packs <- c("purrr"
                  , "magrittr"
                  , "openxlsx"
                  , "ggdist"
                  , "ggplot2"
                  , "analyzeGC")

# install and/or load packages
load_my_packs(script_packs)

## Vector listing objects in the environment that should not be erased during 
## or immediately after running the for loop
important_script_objects <- ls() |> 
  str_remove_all(analysis_objects |> 
                   rev() |>
                   paste(collapse = "|")) |> 
  # select the non-empty entries
  str_subset("^\\S")

load(here("data", "processed",  "grouping_info.Rdata"))

# Add samples_info to important_script_objects
important_script_objects <- c(important_script_objects, "samples_info")

# Batch tables ----
## FOR loop to assemble the "master tables" of each gcms batch (batch table)
for (gcms_batch in gcms_batch_list) {
  # Load group tables 
  load(here("data", "tmp"
            , paste0(gcms_batch
                     , "_group_tables.Rdata")))
  # Firts batch table
  ## Define columns that hold the information of the compounds, not samples'
  ## composition data
  comps_vars <- group_tables |> 
    pluck(1) |> 
    pluck(1) |> 
    select(RI:Mod.position) |> 
    colnames()
  comps_vars
  
  ## Merge all the group tables of the batch into a single batch table
  master_table <- group_tables |> 
    build_master_table()
  str(master_table)
  
  ## Export the batch table into an XLSX file, for check-ups
  wb <- createWorkbook()
  wb |> 
    addWorksheet(sheetName = "master_table")
  wb |> 
    conditionalFormatting(sheet = "master_table"
                          , cols = 2
                          , rows = 1:nrow(master_table |> 
                                            pluck("RT"))
                          , type = "duplicates")
  wb |> 
    writeData(x = master_table |> 
                pluck("Area")
              , sheet = "master_table")
  wb |> 
    saveWorkbook(file = here("data", "tmp"
                             , paste0(gcms_batch
                                      , "_batch_table.xlsx"))
                 , overwrite = T)
  
  ## Extract the abundance composition data from the batch table, in a data frame
  ## These data frames are labeled with "daten"
  master_daten <- master_table |> 
    pluck("Area") |> 
    column_to_rownames("Peak") |> 
    select(!all_of(c(comps_vars))) |> 
    t() |> 
    as.data.frame()
  
  ## NAs are replaced with 0s as they represent the absence of a compound in a 
  ## sample
  master_daten[is.na(master_daten)] <- 0
  str(master_daten)
  
  ## Rename master_daten. The new name is the batch name with the "daten" label
  assign(paste(gcms_batch, quote(master_daten), sep = "_")
         , master_daten)
  
  ## Extract the compound information from the batch table, in a data frame
  ## These data frames are labeled with "Comps"
  master_Comps <- master_table |> 
    pluck("Area") |> 
    select(all_of(comps_vars))
  
  # Ensure that RI reamins as an integer
  master_Comps$RI <- master_Comps$RI |> as.integer()
  str(master_Comps)
  
  ## Rename master_Comps. The new name is the batch name with the "Comps" label
  assign(paste(gcms_batch, quote(master_Comps), sep = "_")
         , master_Comps)
  
  ## Rename grouping_info. The new name is the batch name with the 
  ## "grouping_info" label
  assign(paste(gcms_batch, quote(grouping_info), sep = "_")
         , grouping_info)
  
  ## Export the batch data frames into an Rdata file
  save(list = c(paste(gcms_batch, quote(master_daten), sep = "_")
                , paste(gcms_batch, quote(master_Comps), sep = "_")
                , paste(gcms_batch, quote(grouping_info), sep = "_"))
       , file = here("data", "tmp"
                     , paste0(gcms_batch
                              , "_data-frames.Rdata")))
}

# Decluter memory by removing objects created by the for loop that are not
# listed in important_script_objects nor analysis_objects
objects_list <- ls() |> 
  str_remove_all(c(analysis_objects 
                   , important_script_objects) |>
                   rev() |> 
                   paste(collapse = "|")) |> 
  # select the non-empty entries
  str_subset("^\\S")
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

# First master table ----
## List of Rdata files paths
mastertables_list <- list.files(here("data", "tmp")
                                , pattern = "data-frames.Rdata"
                                , full.names = T) 

## Load the batch data frames
for (mastertable in mastertables_list) {
  load(mastertable)
}

## Data frame indicating the name of the data frames of each batch by type
mastertables_list <- sapply(mastertables_list, load) |> 
  as.data.frame()
colnames(mastertables_list) <- gcms_batch_list
rownames(mastertables_list) <- mastertables_list |> 
  pull(1) |> 
  lapply(str_remove_all
         , gcms_batch_list |> 
           paste0(collapse = "|")) |> 
  str_split("_") |>  
  lapply(str_subset
         , "[:alpha:]") |> 
  lapply(paste, collapse = "_") |> 
  unlist()

## Bind comps and daten objects together ----
# empty vector for listing the batch tables
batch_tables <- c()
# Iterate through batches
for (batch in colnames(mastertables_list)) {
  # Get batch "daten" data frame
  daten_name <- mastertables_list |> 
    select(all_of(batch)) |> 
    filter(row.names(mastertables_list) == "master_daten") |> 
    as.character()
  
  # Get batch "Comps" data frame
  daten_Comps <- mastertables_list |> 
    select(all_of(batch)) |> 
    filter(row.names(mastertables_list) == "master_Comps") |> 
    as.character()
  
  # transpose "daten" data frame
  tmp_daten <- get(daten_name) |> 
    t() |> 
    as.data.frame() 
  
  # Bind the "Comps" and "daten" data frames together to get the batch table
  tmp_table <- get(daten_Comps) |> 
    cbind(tmp_daten) |>
    as_tibble()
  
  # Define batch table name
  name_tmp_table <- paste0(batch, "_table")
  
  # Rename name_tmp_table with the batch table name
  assign(name_tmp_table, tmp_table) 
  
  # Add batch table name to batch_tables
  batch_tables <- c(batch_tables, name_tmp_table)
  
  # declutter memory
  rm(tmp_daten, tmp_table, daten_Comps, daten_name)
}

# Create a list that contaisn all batch tables
batch_tables_list <- lapply(batch_tables, get)

# Name entries in the list as the corresponding batch tables
names(batch_tables_list) <- batch_tables

## Build the master table by merging all the batch tables
master_table <- batch_tables_list |> 
  build_master_table()
str(master_table)

## Assemble grouping_info ----
## Merge the grouping_info data frames of all batches into one
grouping_info <- mastertables_list |> 
  filter(row.names(mastertables_list) == "grouping_info") |> 
  as.character() |> 
  lapply(get) |> 
  reduce(merge
         , all = T
         , sort = F)

# Set sample labels as the row names
rownames(grouping_info) <- grouping_info$Individual

# Shape the grouping_info
grouping_info <- grouping_info |> 
  # reorder columns
  select(Individual, Location:Task) |> 
  # create a group_label column by appending all factors
  unite(group_label, !Individual, sep = "_", remove = F) 

## Export the master table ----
## Export master table to a XLSX file for check-ups
wb <- createWorkbook()
wb |> 
  addWorksheet(sheetName = "master_table_1")
wb |> 
  conditionalFormatting(sheet = "master_table_1"
                        , cols = 2
                        , rows = 1:nrow(master_table)
                        , type = "duplicates")
wb |> 
  writeData(x = master_table
            , sheet = "master_table_1")
wb |> 
  saveWorkbook(file = here("data", "tmp", "master_table.xlsx")
               , overwrite = T)

# Second master table ----

## Define columns that hold the information of the compounds, not samples'
## composition data
comps_vars <- master_table |> 
  select(Peak:Mod.position) |> 
  colnames()

## Get the "daten" data frame of the master table
## It holds the samples' composition data
master_daten <- master_table |> 
  select(!all_of(comps_vars)) |> 
  t() |> 
  as.data.frame()

## Retrieve the group tables from the master table ----
### A list holding the group tables, made by dissaesembling the master table
group_tables_list <- retrieve_group_tables(group.label = "group_label"
                                          , master.table = master_table
                                          , grouping.info = grouping_info)

## Fuse peaks ----
### Duplicated compounds information ----
#### Check for duplicated compounds
pdf(here("output", 'density-distribution_duplicated-compounds.pdf')
    , width = 18, height = 8)
duplicated_compounds_presence <-
  group_tables_list |>
  assess_duplicated_compounds()
dev.off()

duplicated_comps <- master_table |> 
  filter(duplicated(Compound)) |> 
  pull(Compound) |> 
  unique()

#### Export tables reporting the presence/absence of peaks with duplicated 
#### compounds in the different groups, as XLSX files.
#### This helps in defining what peaks should be fused together.
wb <- createWorkbook()
for (dup_comp in duplicated_comps) {
  dup_comp_table <- group_tables_list |> 
    lapply(filter
           , Compound == dup_comp) |> 
    lapply(select, Peak:present) |> 
    lapply(mutate
           , count = ifelse(present == T, 1, 0)
           , .keep = "unused") |> 
    reduce(rbind) |>
    arrange(RI)
  
  for (peak in unique(dup_comp_table$Peak)) {
    dup_comp_table["count"][dup_comp_table["Peak"] == peak] <- 
      dup_comp_table |> 
      filter(Peak == peak) |> 
      pull(count) |> 
      sum() |> 
      as.integer()
  }
  
  dup_comp_table <- dup_comp_table |> 
    filter(!duplicated(Peak)) 
 
  dup_comp_table <- group_tables_list |> 
    lapply(filter, Compound == dup_comp) |> 
    lapply(select, present) |> 
    reduce(cbind) |> 
    set_colnames(names(group_tables_list)) |> 
    mutate(Peak = group_tables_list |> 
             lapply(filter, Compound == dup_comp) |> 
             pluck(1) |> 
             pull(Peak)
           , RI = group_tables_list |> 
             lapply(filter, Compound == dup_comp) |> 
             pluck(1) |> 
             pull(RI)) |> 
    select(Peak, RI, everything()) |> 
    t() |> 
    as.data.frame()
  
  wb |> 
    addWorksheet(sheetName = dup_comp |> str_replace(":","="))
  wb |> 
    conditionalFormatting(sheet = dup_comp |> str_replace(":","=")
                          , cols = 2:(ncol(dup_comp_table) + 1)
                          , rows = 3:nrow(dup_comp_table)
                          , rule = "FALSE"
                          , type = "contains")
  wb |> 
    setColWidths(sheet = dup_comp |> str_replace(":","=")
                 , cols = 1:(ncol(dup_comp_table) +  1)
                 , widths = "auto")
  wb |> 
    freezePane(sheet = dup_comp |> str_replace(":","=")
               , firstActiveRow = 2
               , firstActiveCol = 2)
  wb |> 
    writeData(x = dup_comp_table
              , sheet = dup_comp |> str_replace(":","=")
              , rowNames = T
              , colNames = F)
}

wb |> 
  saveWorkbook(file = here("data", "tmp"
                           , "group-presence_duplicated-compounds.xlsx")
               , overwrite = T)

### Fuse the peaks ----
### Each entry of the list corresponds to the names of two or more peaks that
### have to be fused together across the data set (master table).
fusion_list <- {list(
  # C19_ene
  c(paste0("P"
           , c(5, 6, 7)))
  # C21_ene
  , c(paste0("P"
             , c(12, 13)))
  # C22_ene
  , c(paste0("P"
             , c(15, 16)))
  , c(paste0("P"
             , c(17, 18)))
  # C23_diene
  ,c(paste0("P"
            , c(20, 21, 22)))
  # C23_ene
  ,c(paste0("P"
            , c(23, 24, 25)))
  ,c(paste0("P"
            , c(26, 27, 28, 29, 30)))
  # 9-;11-MeC23
  ,c(paste0("P"
            , c(32, 33, 34)))
  # C24_ene
  ,c(paste0("P"
            , c(36, 37)))
  # C25_diene
  ,c(paste0("P"
            , c(39, 40)))
  ,c(paste0("P"
            , c(41, 42)))
  # C25_ene
  ,c(paste0("P"
            , c(43, 44, 45)))
  ,c(paste0("P"
            , c(46, 47, 48)))
  # 11-;13-;15-Me25
  ,c(paste0("P"
            , c(50, 51, 52)))
  # C27_ene
  ,c(paste0("P"
            , c(57, 58)))
  ,c(paste0("P"
            , c(59, 60, 61)))
  # 11-;13-MeC27
  ,c(paste0("P"
            , c(63, 64, 65)))
  # 11-,15--DimeC27
  ,c(paste0("P"
            , c(66, 67)))
  # 13-;14-MeC28
  ,c(paste0("P"
            , c(69, 70)))
  # C29_ene
  ,c(paste0("P"
            , c(71, 72, 73, 74, 75)))
  # 11-;13-15-MeC29
  ,c(paste0("P"
            , c(77, 78, 79)))
  # 11-,17-DimeC29
  ,c(paste0("P"
            , c(80, 81, 82)))
  # C31-diene
  ,c(paste0("P"
            , c(85, 86, 87)))
  # C31_ene
  ,c(paste0("P"
            , c(89, 90, 91, 92)))
  ,c(paste0("P"
            , c(93, 94)))
  # 13-;15-MeC31
  ,c(paste0("P"
            , c(96, 97, 98)))
  #C32_ene
  ,c(paste0("P"
            , c(100, 101)))
  # C33_diene
  ,c(paste0("P"
            , c(103, 104, 105)))
  ,c(paste0("P"
            , c(106, 107)))
  # C33_ene
  ,c(paste0("P"
            , c(108, 109, 110, 111)))
  #C35_diene
  ,c(paste0("P"
            , c(114, 115, 116)))
  #C35_ene
  ,c(paste0("P"
            , c(117, 118, 119, 120)))
)}

fusion_list

### Fuse the peaks in the master table
master_table <- fuse_all_peaks(master_table, fusion_list)

## Export the master table ----
## Export master table to a XLSX file for check-ups
wb <- loadWorkbook(here("data", "tmp", "master_table.xlsx"))
if (sheets(wb) |> str_detect("master_table_fusedpeaks") |> any()) {
  wb |> 
    removeWorksheet(sheet = "master_table_fusedpeaks")
}
wb |> 
  addWorksheet(sheetName = "master_table_fusedpeaks")
wb |> 
  conditionalFormatting(sheet = "master_table_fusedpeaks"
                        , cols = 2
                        , rows = 1:nrow(master_table)
                        , type = "duplicates")
wb |> 
  writeData(x = master_table
            , sheet = "master_table_fusedpeaks")
wb |> 
  saveWorkbook(file = here("data", "tmp", "master_table.xlsx")
               , overwrite = T)

# Third master table ----
## Retrieve the group tables from the master table ----
group_tables_list <- retrieve_group_tables("group_label"
                                           , master.table = master_table
                                           , grouping.info = grouping_info)

### Filter rare compounds of each group ----
# Use the group_frequency_filter function to remove from each group table
# the compounds that are not present in at least 50% (default threshold) of the
# samples within the corresponding group
group_tables_list <- group_tables_list |> 
  lapply(group_frequency_filter)

## rebuild the master table
master_table <- group_tables_list |> 
  build_master_table()

## Export the master table ----
## Export master table to a XLSX file for check-ups
wb <- loadWorkbook(here("data", "tmp", "master_table.xlsx"))
if (sheets(wb) |> str_detect("master_table_50%filter") |> any()) {
  wb |> 
    removeWorksheet(sheet = "master_table_50%filter")
}
wb |> 
  addWorksheet(sheetName = "master_table_50%filter")
wb |> 
  conditionalFormatting(sheet = "master_table_50%filter"
                        , cols = 2
                        , rows = 1:nrow(master_table)
                        , type = "duplicates")
wb |> 
  writeData(x = master_table
            , sheet = "master_table_50%filter")
wb |> 
  saveWorkbook(file = here("data", "tmp", "master_table.xlsx")
               , overwrite = T)

# Final master table ----
## Use the abundance_transformation to transform the abundance of the compounds
## to percentage (default transformation)
master_table <- master_table |> 
  abundance_transformation()

## Repair compound info ----
### As we fused some peaks together, the name of the compound(s) within the 
### peaks needs to be repaired
### At this point the double bond positions for the alkenes was found via DMDS 
### derivatization. So, we could add this to the data set.
master_table <- master_table |> 
  rows_update(tribble(~Peak, ~RI, ~Compound, ~Mod.position
                      , "P3", 1874, "9-C19:ene", "9-"
                      , "P7", 2073, "9-C21:ene", "9-"
                      , "P14", 2273, "9-C23:ene", "9-"
                      , "P15", 2280, "7-C23:ene", "7-"
                      , "P17", 2336, "9-;11-MeC23", "9-;11-"
                      , "P18", 2373, "9-C24:ene", "9-"
                      , "P19", 2382, "7-C24:ene", "7-"
                      , "P23", 2473, "9-C25:ene", "9-"
                      , "P24", 2480, "7-C25:ene", "7-"
                      , "P26", 2534, "11-;13-MeC25", "11-;13-"
                      , "P30", 2674, "9-C27:ene", "9-"
                      , "P31", 2682, "7-C27:ene", "7-"
                      , "P33", 2732, "11-;13-MeC27", "11-;13-;15-"
                      , "P34", 2760, "11-,15-DimeC27", "11-,15-"
                      , "P36", 2832, "13-;14-MeC28", "13-;14-"
                      , "P37", 2879, "9-C29:ene", "9-"
                      , "P39", 2931, "11-;13-;15-MeC29", "11-;13-;15-"
                      , "P40", 2960, "11-,17-DimeC29", "11-,17-"
                      , "P41", 2982, "9-C30:ene", "9-"
                      , "P45", 3078, "10-C31:ene", "10-"
                      , "P46", 3082, "8-C31:ene", "9-"
                      , "P48", 3130, "13-;15-MeC31", "13-;15-"
                      , "P49", 3157, "13-,17-DimeC31", "13-,17-"
                      , "P50", 3178, "8-;10-C32:ene", "8-;10-"
                      , "P54", 3278, "10-;12-C33:ene", "10-;12-"
                      , "P56", 3327, "13-;15-;17-MeC33", "13-;15-;17-"
                      , "P58", 3470, "10-;12-C35:ene", "10-;12-"))

## Export the final master table ----
## Export master table to a XLSX file for check-ups
wb <- loadWorkbook(here("data", "tmp", "master_table.xlsx"))
if (sheets(wb) |> str_detect("master_table_final") |> any()) {
  wb |> 
    removeWorksheet(sheet = "master_table_final")
}
wb |> 
  addWorksheet(sheetName = "master_table_final")
wb |> 
  conditionalFormatting(sheet = "master_table_final"
                        , cols = 2
                        , rows = 1:nrow(master_table)
                        , type = "duplicates")
wb |> 
  writeData(x = master_table
            , sheet = "master_table_final")
wb |> 
  saveWorkbook(file = here("data", "tmp", "master_table.xlsx")
               , overwrite = T)

## Export master data frames to be used by subsequent analysis scripts
### Define columns that hold the information of the compounds, not samples'
### composition data
comps_vars <- master_table |> 
  select(Peak:Mod.position) |> 
  colnames()

### Extract the "daten" data frame from the master table
master_daten <- master_table |> 
  select(!all_of(comps_vars)) |> 
  t() |> 
  as.data.frame()

### NAs in master_daten represent the absence of a compound in a given sample
### So, they are replaced by 0s
master_daten[is.na(master_daten)] <- 0
str(master_daten)

### Extract teh "Comps" data frame from the master table
master_Comps <- master_table |> 
  select(all_of(comps_vars))

### Ensure the RIs are integers
master_Comps$RI <- master_Comps$RI |> as.integer()
str(master_Comps)

### Set peak labels as column names for master_daten
colnames(master_daten) <- master_Comps$Peak
str(master_daten)

### Grouping_info now holds complete (non-abbreviated) names
grouping_info <- samples_info |> 
  filter(Individual %in% grouping_info$Individual)

# export the master data frames into an Rdata file
save(list = c("master_daten"
              , "master_Comps"
              , "grouping_info")
     , file = here("data", "processed", "mastertable_data-frames.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_mastertable-script.txt"))


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
