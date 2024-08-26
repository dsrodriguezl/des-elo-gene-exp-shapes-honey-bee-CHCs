tic("MS alignment script")
# Packages ----

install_new_packs("viridis")
# Vector holding the list of packages that will be used
script_packs <- c("analyzeGC"
                  , "ggplot2"
                  , "GCalignR"
                  , "doParallel")

# install / load packages
load_my_packs(script_packs)

# Load the data ----
load(here("data", "tmp", "data2align"
          , paste0(gcms_batch
                   , ".Rdata")))
  
# Data alignment ----
## Samples 
## Let's parallelize this
## register the parallel backend
no_cores <- detectCores() / 2 
cl <- makeCluster(no_cores)  
clusterExport(cl, varlist = c("linear_shift_criteria"
              , "partial_alignment_threshold"
              , "row_merging_threshold"))
registerDoParallel(cl)   

# Perform the alignment(s) of samples
aligned_samples_data_list <- 
  parLapply(cl
            , samples_data_list
            , align_chromatograms2
            , linear_shift_criteria = linear_shift_criteria
            , partial_alignment_threshold = partial_alignment_threshold
            , row_merging_threshold = row_merging_threshold)
print(aligned_samples_data_list)

## stop the cluster of the parallel back-end to free-up resources
stopCluster(cl)
registerDoSEQ() 

## STD 
# perform alignment of standards
aligned_STD_data_list <- standards_list |> 
  align_chromatograms2(linear_shift_criteria = linear_shift_criteria
                       , partial_alignment_threshold = partial_alignment_threshold
                       , row_merging_threshold = row_merging_threshold)
print(aligned_STD_data_list)

# Check precision of the alignment(s) ----
## Get aligned RT data frames
### Samples 
samples_list_RT <- aligned_samples_data_list |> lapply(extract_RT)

### STD 
STD_RT <- aligned_STD_data_list |> extract_RT()

# Export aligned data frames  ####
## CSV of RT data frames ####
## This is useful to check up the alignments, while looking at the chromatograms
## in the GC-MS analysis software.

print("Exporting CSV files of aligned RT data frames of samples")
for (df in names(samples_list_RT)) {
  write.csv(samples_list_RT[[df]]
            , here("data", "tmp", "data2align"
                   , paste0("aligned-RT_", gcms_batch, "_", df, ".csv")))
  print(paste0("aligned-RT_", gcms_batch, "_", df, ".csv", " was exported"))
  rm(df)
}

print("Exporting CSV file of aligned RT data frame of standards")
write.csv(STD_RT
          , here("data", "tmp", "data2align"
                 , paste0("aligned_RT_"
                          , gcms_batch
                          ,"_STD.csv"))
          , row.names = T)
print(paste0("aligned_RT_", gcms_batch,"_STD.csv", " was exported"))

## data frames in a Rdata file ####
save(list = c("aligned_samples_data_list"
              , "aligned_STD_data_list")
     , file = here("data", "tmp", "data2align"
                   , paste0("uncorrected-alignment_"
                            , gcms_batch
                            , ".Rdata")))
print("The aligned data frames were exported")

# Export diagnostic plots ####

#### Generate PDF file to contain plots
print("Exporting plots to diagnose the alignments")
{pdf(here("output"
          , paste0("uncorrected-alignment-plots_"
                   , gcms_batch
                   , '.pdf'))
    , width = 30
    , height = 15)}

## samples ####
for (df in names(aligned_samples_data_list)) {
  diagnostic_heatmap(aligned_samples_data_list[[df]]
                     , title = paste0("Alignment of "
                                      , df)
                     , alignment.type = "automatic")
}

## STD ####
diagnostic_heatmap(aligned_STD_data_list
                   , title = "Alignment of standards"
                   , alignment.type = "automatic")

#### Close graphic device to export plots into the PDF file
dev.off()
print("Diagnostic plots for the alignments were exported")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_align-gcms-script.txt"))

## Detach/unload packages
pacman::p_unload(char = script_packs)

toc()
## Clear environment
objects_list <- ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|"))
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

gc()
