tic("Script to import qPCR data")

# Packages ----
# Vector holding the list of packages that will be used by the script
script_packs <- c("purrr", "readr", "ggplot2", "magrittr")

# install and/or load packages
load_my_packs(script_packs)

# FUNCTIONS ----

## A function to import and begin to structure the raw Ct
## from the CSV files obtained from the Rotor-Gene software
import_qpcr_data <- function(data.path.list, type) {
  # Procedure if importing samples data
  if (type == "samples") {
    # Get the gene_run names
    names_list <- str_split(data.path.list
                            , "/"
                            , simplify = T) |> 
      str_subset(".csv") |> 
      str_remove(".csv") |> 
      str_split("-", simplify = T) |> 
      str_subset("[:alpha:]")
    
    # Import the data
    run.data.list <- data.path.list |> 
      lapply(read_csv
             , col_select = c(Name
                              , Ct))
  }
  
  # Procedure if importing data from standard curves
  if (type == "standards") {
    names_list <- str_split(data.path.list
                            , "/"
                            , simplify = T) |> 
      str_subset(".csv") |> 
      str_remove(".csv") |> 
      str_split("_", simplify = T) |> 
      str_subset("[:alpha:]") |> 
      str_subset("SK", negate = T)
    
    run.data.list <- data.path.list |> 
      lapply(read_csv
             , col_select = c(Name
                              , Ct
                              , 7))
  }
  
  # Define names of list entries
  names(run.data.list) <- names_list
  run.data.list
  
  # Iterate through runs structuring the data frames with the Ct data 
  # of each gene
  for (run in names(run.data.list)) {
    # print(run)
    run_df <- run.data.list[[run]]
    
    if (type == "standards") {
      gene <- run 
      
      colnames(run_df)[3] <- "amount"
      
      run_df <- run_df |> 
        filter(!is.na(amount) & Name != "ntc") |> 
        select(Ct, amount)
      
      colnames(run_df) <- c(gene, "amount")
      
      run_df <- run_df |> 
        mutate(Replicate = ifelse(duplicated(amount)
                               , paste0(amount, "_2")
                               , paste0(amount, "_1"))) |>
        mutate(Replicate = ifelse(duplicated(Replicate)
                               , paste0(Replicate |> 
                                          str_remove("_2")
                                        , "_3")
                               , Replicate))
    }
    
    if (type == "samples") {
      # Extract gene name
      gene <- run |>
        str_split("_", simplify = T) |> 
        str_subset("[:alpha:]")
      
      # Extract run number
      run_number <- run |>
        str_split("_", simplify = T) |>
        str_subset("[:digit:]") |>
        str_subset("[:alpha:]", negate = T)
      
      run_df <- run_df |>
        filter(get("Name") != "ntc") |>
        filter(get("Name") != "nct") |>
        mutate("Name" = ifelse(Name == "Cal"
                               , paste0(Name, run_number)
                               , Name))
      colnames(run_df) <- c("sample", gene)
      run_df <- run_df |>
        mutate(sample = ifelse(duplicated(sample)
                               , paste0(sample, "_2")
                               , paste0(sample, "_1"))) |>
        mutate(sample = ifelse(duplicated(sample)
                               , paste0(sample |>
                                          str_remove("_2")
                                        , "_3")
                               , sample))
    }
    
    run.data.list[[run]] <- run_df
    
  }
  
  if (type == "samples") {
    runs <- run.data.list |> names() |> 
      str_split("_", simplify = T) |> 
      str_subset("[:alpha:]", negate = T) |> 
      unique()
    
    mg_list <- list()
    for (run in runs) {
      # print(run)
      run_genes <- run.data.list |> 
        names() |> 
        str_subset(paste0("_", run))
      run_list <- run.data.list |> 
        purrr::keep(run.data.list |> 
                      names() %in% run_genes) |> 
        as.list()
      
      mg_list[[run]] <- run_list
    }
    run.data.list <- mg_list
  }
  
  run.data.list
}

## A function to clean the imported qpcr data
### It removes technical outliers from the replicates of a sample in a qPCR run
### for a given gene.
### The CT value from replicates of a sample are removed from the data set,
### if their CT value is farther than 2 to the median CT value of the sample.
clean_qpcr_data <- function(run.data.list) {
  run.data.list |> 
    imap(function(run_list, run_id) {
      # print(paste("Processing run:", run_id))
      
      # Iterate through genes within run
      run_list <- run_list |> 
        imap(function(gene_df, run_gene) {
          # Define gene name
          gene_id <- run_gene |> str_remove_all("_[:digit:]")
          print(paste("Gene:", gene_id))
          
          # Shape gene data frame
          gene_df <- gene_df |> 
            mutate(Individual = sample |> str_remove_all("_[:digit:]")) |> 
            pivot_longer(cols = all_of(gene_id)
                         , names_to = "gene"
                         , values_to = "ct") 
          
          # Calculate median CT value per sample, using the CT value of the replicates
          summary_df <- gene_df|>
            group_by(Individual, gene) |>
            summarise(median_ct = median(ct, na.rm = T))
          
          # Determine which replicates are outliers and correct their CT values 
          # by turning them into NAs
          gene_df <- gene_df |> 
            full_join(summary_df) |> 
            select(-gene) |> 
            mutate(outlier = abs(median_ct - ct) > 2
                   , corrected_ct = ifelse(abs(median_ct - ct) > 2
                                           , NA
                                           , ct))
          
          # Plot illustrating median, range for considering a CT value an outlier
          # and the replicates per sample
          plot <- gene_df |>
            ggplot(aes(x = Individual
                       , y = ct)) +
            geom_pointrange(aes(y = median_ct
                                , ymin = median_ct - 2
                                , ymax = median_ct + 2)
                            , shape = 3
                            , linewidth = 1.5
                            , size = 2
                            , color = "grey40") +
            geom_point(fill = "red"
                       , shape = 21
                       , size = 2
            ) +
            labs(title = paste0("Run ", run_id, ": ", gene_id)
                 , subtitle = paste("CT values per replicate are visualized"
                                    , "as red dots."
                                    , "\nMedian CT and criteria to remove"
                                    , "technical outliers are shown as"
                                    , "bar-ranges.")) +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90)
                  , plot.title = element_text(face = "bold"))
          print(plot)
          
          gene_df
        })
      
      # Iterate through genes within run
      run_list |>
        imap(function(gene_df, run_gene) {
          # Define gene name
          gene_id <- run_gene |> str_remove_all("_[:digit:]")
          
          # Count the number of replicates left per sample
          summary_df <- gene_df |> 
            group_by(Individual) |> 
            summarise(count = sum(!is.na(corrected_ct)))
          
          # Remove the CT values of any sample that has only one replicate left
          gene_df |> 
            full_join(summary_df) |> 
            mutate(corrected_ct = ifelse(count > 1
                                         , corrected_ct
                                         , NA)) |> 
            # Shape the data frame
            select(sample, corrected_ct) |>
            set_colnames(c("sample", gene_id))
        })
    })
}

## A function to add the well and Run information to each run within 
## run_data_list
add_run_info <- function(run.data.list) {
   run_names <- names(run.data.list)
   for (run in run_names) {
     run.data.list[[run]] <- run.data.list[[run]] |>
       mutate(Well = 1:nrow(run.data.list[[run]])
              , Run = run
              , Sample = sample
              , .keep = "unused")
   }
   run.data.list
}

# Import data files ----
## Samples
run_data_list <- 
  # Make a list of the individual sample files paths
  list.files(path = here("data", "raw", "qPCR", "samples_runs")
             # Get all CSV files in the folder
             , pattern = ".csv" 
             , full.names = T
             , recursive = T) |>
  # Import the data in the csv within a structured list
  import_qpcr_data(type = "samples")

pdf(here("output"
         , paste0("ct-data-filtering-plot.pdf"))
    , width = 18, height = 8)

# Remove technical outliers
run_data_list <- run_data_list |> 
  clean_qpcr_data()
dev.off()


## standard curves
standard_curves_list <- 
  # Make a list of the individual sample files paths
  list.files(path = here("data", "raw", "qPCR")
             # Get all CSV files in the folder
             , pattern = ".csv" 
             , full.names = T
             , recursive = F) |>
  # Import the data in the csv within a structured list
  import_qpcr_data(type = "standards")

# Shape the data frames ----
## Samples
run_data_df <- run_data_list |> 
  # Merge gene data frames within each run 
  lapply(reduce
         , merge
         , by = "sample"
         , all = T) |> 
  # Add well and run information to each run data frame
  add_run_info() |>
  # Bind all run data frames into one
  list_rbind()

run_data_df <- run_data_df |> 
  # Place the calibrators at the end of the data frame
  filter(Sample |> 
           str_detect("Cal"
                      , negate = T)) |> 
  rbind(run_data_df |> 
          filter(Sample |> 
                   str_detect("Cal"))) |> 
  # Remove replicate number from the sample label
  mutate(Sample = Sample |> 
           str_remove_all(paste("_1"
                                , "_2"
                                , "_3"
                                , sep = "|"))) |> 
  # Reorganize the columns
  select(Well:Sample, starts_with("Rpl"), everything()) |> 
  # Transform in a tibble
  as_tibble()

run_data_df

## standard curves
standard_curves_df <-  standard_curves_list |> 
  # Merge all gene data frames into one
  reduce(merge
         , all = T) |> 
  # Sort rows by amount in a descending order
  arrange(desc(amount)) |> 
  # Set dilution replicate label
  mutate(Replicate = rep(3:8
                         , each = 3)) |> 
  # Remove column amount
  select(-amount) |> 
  # transform into a tibble
  as_tibble()

standard_curves_df  

# Export the data frames ----
destination_folder <- here("data", "tmp", "qPCR")
if (!dir.exists(destination_folder)) {
  dir.create(destination_folder)
}

save(list = c("run_data_df", "standard_curves_df")
     , file = here(destination_folder, "CT-data.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_import-qpcr-script.txt"))

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


