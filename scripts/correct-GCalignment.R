tic("Script to correct alignments")

# Packages ----
# Vector holding the list of packages that will be used
script_packs <- c("ggplot2", "openxlsx", "analyzeGC")

# install and/or load packages
load_my_packs(script_packs)

# Load data ----
load(here("data", "tmp", "data2align"
          , paste0("uncorrected-alignment_"
                   , gcms_batch
                   , ".Rdata")))

# Correct miss-alignments ----
## Create empty rows within the data frames ####
### List with instructions for adding empty "peaks" to the data frames ####
empty_peaks <- {list("0061" = tribble(~position.reference, ~direction,
                                      "P140", "before")
                     , "0156" = tribble(~position.reference, ~direction,
                                    "P35", "after"))
}

### Add empty peaks to the data frames ####
aligned_samples_data_list <- aligned_samples_data_list |>
lapply(add_empty_peaks
      , empty.peaks = empty_peaks)

## List with the instructions for moving the peaks within samples ####
#batch 20022023
peaks_movements <- {list(#batch 20022023  
                        "0013" = data.frame(peaks_list = c("P89")
                                            , movement_dirs = c('up'))
                        , "0051" = data.frame(peaks_list = c("P119")
                                              , movement_dirs = c('down'))
                        , "0053" = data.frame(peaks_list = c("P119")
                                              , movement_dirs = c('down'))
                        , "0060" = data.frame(peaks_list = c("P119")
                                              , movement_dirs = c('down'))
                        , "0062" = data.frame(peaks_list = c("P202")
                                              , movement_dirs = c('down'))
                        , "0065" = data.frame(peaks_list = c("P141")
                                              , movement_dirs = c('up'))
                        , "0066" = data.frame(peaks_list = c("P141"
                                                             , "P202")
                                              , movement_dirs = c('up'
                                                                  , 'down'))
                        , "0067" = data.frame(peaks_list = c("P141")
                                              , movement_dirs = c('up'))
                        , "0068" = data.frame(peaks_list = c("P141")
                                              , movement_dirs = c('up'))
                        , "0069" = data.frame(peaks_list = c("P141")
                                              , movement_dirs = c('up'))
                        , "0070" = data.frame(peaks_list = c("P140"
                                                             , "P141"
                                                             , "P202")
                                              , movement_dirs = c('up'
                                                                  , 'up'
                                                                  , 'down'))
                        , "0093" = data.frame(peaks_list = c("P246")
                                              , movement_dirs = c('up'))
                        , "0095" = data.frame(peaks_list = c("P172")
                                              , movement_dirs = c('down'))
                        , "0096" = data.frame(peaks_list = c("P172")
                                              , movement_dirs = c('down'))
                        , "0097" = data.frame(peaks_list = c("P109")
                                              , movement_dirs = c('down'))
                        , "0099" = data.frame(peaks_list = c("P172","P246")
                                              , movement_dirs = c('down','up'))
                        , "0100" = data.frame(peaks_list = c("P246")
                                              , movement_dirs = c('up'))
                        , "0102" = data.frame(peaks_list = c("P246"
                                                             , "P259"
                                                             , "P258")
                                              , movement_dirs = c('up'
                                                                  , 'down'
                                                                  , 'down'))
                        , "0140" = data.frame(peaks_list = c("P158")
                                              , movement_dirs = c('down'))
                        , "0144" = data.frame(peaks_list = c("P211")
                                              , movement_dirs = c('up'))
                        , "0145" = data.frame(peaks_list = c("P211")
                                              , movement_dirs = c('up'))
                        , "0150" = data.frame(peaks_list = c("P211")
                                              , movement_dirs = c('up'))
                        #batch 06032023
                        , "0154" = data.frame(peaks_list = c("P169")
                                              , movement_dirs = c('up'))
                        , "0156" = data.frame(peaks_list = c("P35"
                                                             , "P34")
                                              , movement_dirs = c('down'
                                                                  , 'down'))
                       )}

## correct alignments ####
corrected_samples_data_list <- aligned_samples_data_list |>
  # Iterate through alignments in the list, to correct them all
  lapply(correct_alignment
         , movements_list = peaks_movements)

print("Alignment was corrected. Check heat map to verify that it is correct.")

# Check the corrected alignment ----
## PDF to export diagnostic heat maps of the corrected alignments
{pdf(here("output"
         , paste0(gcms_batch
                  , "_corrected-alignment-heatmap"
                  , '.pdf'))
    , width = 20, height = 10)}

## samples ####
for (df in names(corrected_samples_data_list)) {
  diagnostic_heatmap(corrected_samples_data_list[[df]]
                     , title = paste0("corrected alignment of "
                                      , df)
                     , alignment.type = "corrected")
}

dev.off()
print("Heat map to verify corrected alignment was exported")

# Export the aligned data for CHC identification ----
# correct the mean RT in all the alignments within the list
corrected_samples_list2 <- lapply(corrected_samples_data_list
                                  , recalculate_meanRT)

## Save the data frames 
save(list = c("corrected_samples_list2"
              , "aligned_STD_data_list"
              )
     , file = here("data", "tmp"
                   , paste0(gcms_batch
                            , "_aligned_gcms-data.Rdata")))
print("The aligned data frames were exported")

# Extract RT
corrected_samples_list_RT <- corrected_samples_list2 |> 
  lapply(extract_RT)

STD_RT <- extract_RT(aligned_STD_data_list) |> 
  rownames_to_column("Peak")

# Export RT data frames to XLSX files
## This is useful to check up the alignments, while looking at the chromatograms
## in the GC-MS analysis software.
## In addition, the compounds identification can be manually added in the 
## "identification" page of the XLSX file
for (df in names(corrected_samples_list_RT)) {
  group_table_file <- here("data", "tmp"
                           , paste0(gcms_batch, "_", df, "_table.xlsx"))
  if (file.exists(group_table_file)) {
    wb <- loadWorkbook(group_table_file)
    wb |> 
      removeWorksheet(sheet = "corrected_alignment")
    wb |> 
      addWorksheet(sheetName = "corrected_alignment")
    wb |> 
      writeData(x = corrected_samples_list_RT[[df]]
                , sheet = "corrected_alignment")
    wb |> 
      saveWorkbook(file = group_table_file
                   , overwrite = T)
    print(paste0(gcms_batch, "_", df, "_table.xlsx", " was modified"))
  } else {
    wb <- createWorkbook()
    wb |> 
      addWorksheet(sheetName = "corrected_alignment")
    wb |> 
      writeData(x = corrected_samples_list_RT[[df]]
                , sheet = "corrected_alignment")
    if (!str_detect(paste(names(wb), collapse = "_"), "identification")) {
      wb |>
        addWorksheet(sheetName = "identification")
    }
    wb |> 
      saveWorkbook(file = group_table_file)
    print(paste0(gcms_batch, "_", df, "_table.xlsx", " was exported"))
  }
  rm(df)
}

standards_table_file <- here("data", "tmp"
                             , paste0(gcms_batch
                                      , "_STD_RT_table.xlsx"))
if (file.exists(standards_table_file)) {
  wb <- loadWorkbook(standards_table_file)
  wb |> 
    removeWorksheet(sheet = "corrected_alignment")
  wb |> 
    addWorksheet(sheetName = "corrected_alignment")
  wb |> 
    writeData(x = STD_RT
              , sheet = "corrected_alignment")
  wb |> 
    saveWorkbook(file = standards_table_file
                 , overwrite = T)
  print(paste0(gcms_batch, "_STD_RT_table.xlsx", " was modified"))
} else {
  wb <- createWorkbook()
  wb |> 
    addWorksheet(sheetName = "corrected_alignment")
  wb |> 
    writeData(x = STD_RT
              , sheet = "corrected_alignment")
  if (!str_detect(paste(names(wb), collapse = "_"), "identification")) {
    wb |>
      addWorksheet(sheetName = "identification")
  }
  wb |> 
    saveWorkbook(file = standards_table_file
                 , overwrite = T)
  print(paste0(gcms_batch, "_STD_RT_table.xlsx", " was modified"))
}

print("The RT tables for CHC identification was exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_alignment-correction-script.txt"))

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
