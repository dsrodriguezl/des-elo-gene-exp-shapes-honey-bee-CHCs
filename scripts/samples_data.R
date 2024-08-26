tic("Script to process samples data")

# Packages ----
# Vector holding the list of packages that will be used by the script
script_packs <- c("readr", "lubridate", "forcats")

# install and/or load packages
load_my_packs(script_packs)

# Load the data set ----
samples_info <- read_csv(here("data", "raw", "samples_data.csv")) |> 
  # Do not include empty rows
  filter(!is.na(Collection.date)) |> 
  # Rename column
  mutate(Date = Collection.date
         , .keep = "unused") |> 
  # Extract year and month from date
  mutate(Year = year(Date)
         , Month = month(Date)
         , .keep = "unused") |> 
  # Reorder columns
  select(Individual, Year, Month, everything())

# Create a data frame, where the categorical variables are defined as factors
# This will be used for data wrangling operations
grouping_info <- samples_info |> 
  mutate_if(is.character, as.factor) |> 
  mutate(Individual = Individual |> as.character())

# Transform abbreviations into the real group names, then define the categorical
# variables as factors
# This data frame will be used for plotting operations and statistical analyses 
samples_info <- samples_info |> 
  mutate(Subspecies = ifelse(Subspecies == "Ca"
                               , "A. m. carnica"
                               , "A. m. iberiensis")
         , Task = ifelse(Task == "Nu"
                         , "Nurse bees"
                         , "Forager bees")) |> 
  mutate_at(vars(!contains("Individual"))
            , as.factor) |> 
  mutate(Task = fct_relevel(Task,"Nurse bees"))

# Export data ----
save(list = c("grouping_info", "samples_info")
     , file = here("data", "processed", "grouping_info.Rdata"))

print("The samples data has been exported to the processed data folder")

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "sInf_data-script.txt"))

## Detach/unload packages
pacman::p_unload(char = script_packs)

## Clear environment
objects_list <- ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|"))
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

toc()
