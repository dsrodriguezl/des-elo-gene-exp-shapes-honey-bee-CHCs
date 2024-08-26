tic("Script to calculate the CNRQs values from the qPCR data")

# Packages ----
# Install packages if they are not yet
install_new_packs(c("devtools", "matrixStats", "broom"))

# Install easyqpcr2 v0.1.0 if the package is not installed
if ("easyqpcr2" %in% 
    rownames(installed.packages()) == FALSE) {
  devtools::install_github("dsrodriguezl/easyqpcr2"
                           , dependencies = TRUE
                           , build_vignettes = TRUE
                           , ref = "V0.1.0")
}

# Vector holding the list of packages that will be used by the script
script_packs <- c("purrr"
                  , "easyqpcr2")

# install and/or load packages
load_my_packs(script_packs)

# Load the data ----
load(here("data", "tmp", "qPCR", "CT-data.Rdata"))

# Calculate the amplification efficiencies ----
efficiencies <- amp_efficiency(data = standard_curves_df
                               , q = c(50000000
                                       , 5000000
                                       , 500000
                                       , 50000
                                       , 5000
                                       , 500)
                               , r = 3
                               , na.rm = TRUE)

# calibrator factors ----
## Get the names of the reference genes in a vector
reference_genes_names <- run_data_df |> 
  # Both of my reference genes' names begin with "Rpl"
  select(starts_with("Rpl")) |> 
  colnames()

## Calculate the calibration factors
cal_factors_list <- run_data_df |> 
  calibration_factors(n_replicates = 3
                      , reference_genes = reference_genes_names
                      , amp_efficiencies = efficiencies)

# calculate Normalized RQs  ----
nrqs_list <- run_data_df |> 
  nrqs(n_replicates = 3
       , reference_genes = reference_genes_names
       , amp_efficiencies = efficiencies
       , cal.factors.list = cal_factors_list)

# Calibrated NRQs (for each run) ----
qpcr_rel_data <- nrqs_list |> 
  cnrqs(run.data.df = run_data_df)

# # Normalize with respect to the control group (geometrical mean) ----
load(here("data", "processed", "grouping_info.Rdata"))

groups_qpcr <- samples_info |> 
  filter(Individual %in% row.names(qpcr_rel_data))

# Export the data frames ----
# remove the calibrator
qpcr_rel_data <- qpcr_rel_data |> 
  filter(!str_detect(row.names(qpcr_rel_data) 
                              , "Cal")) 

# Transform into a tibble
qpcr_rel_data <- qpcr_rel_data |> 
  mutate(Individual = row.names(qpcr_rel_data)) |> 
  as_tibble()

destination_folder <- here("data", "processed")
# Create destination folder if it doesn't exist
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

save(list = c("qpcr_rel_data", "groups_qpcr")
     , file = here(destination_folder, "CNQRs.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_CNRQs-scaling-script.txt"))


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