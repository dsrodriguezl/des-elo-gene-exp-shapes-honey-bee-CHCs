tic("Script producing NMDS plots")

# set-up ----
# Vector holding the list of packages that will be used
script_packs <- c('vegan'
                  , 'goeveg'
                  , "purrr"
                  , "viridis"
                  , "ggplot2"
                  , "patchwork"
                  , "ragg"
                  , "ggtext")

# install / load packages
load_my_packs(script_packs)

# get some ggplot presets
source(here("scripts", "custom-ggplot-settings.R"))

# Load the data
load(here("data", "processed", "mastertable_data-frames.Rdata"))
load(here("data", "processed", "chc-v-genes_correlations.Rdata"))

# Perform the NMDS ----
# Set a seed to ensure reproducibility
set.seed(123456) 
sol <- metaMDS(master_daten
                # Use Bray-Curtis dissimilarity index
               , distance = "bray"
               # Two dimensional NMDS
               , k = 2 
               # Minimum number of random starts
               , try = 200
               # Maximum number of random starts
               , trymax = 400
               # Avoid automatic transformation of the data
               , autotransform = F)
sol

## Evaluating the NMDS
stressplot(sol)

## Export NMDS object in an Rdata file
save("sol", file = here("data", "processed", "nmds.Rdata"))

# Make a new data frame with the NMDS coordinates of the samples
NMDS <- sol |> 
  # Get NMDS coordinates
  scores(tidy = T) |> 
  # Preserve only coordinates of samples, not variables
  filter(score == "sites") |> 
  # Rename column label as Individual
  mutate(Individual = label
         , .keep = "unused") |> 
  # Remove column score
  select(-score) |> 
  # Merge with data frame containing samples metadata
  merge(grouping_info, all = T) |> 
  # Transform into a tibble
  as_tibble()
NMDS

# Make a new data frame with the NMDS coordinates of the compounds (varaibles)
NMDS_comps <- sol |> 
  # Get NMDS coordinates
  scores(tidy = T) |> 
  # Preserve only coordinates of variables, not samples
  filter(score == "species") |> 
  # Rename column label as Peak
  mutate(Peak = label
         , .keep = "unused") |> 
  # Merge with data frame containing compounds information
  merge(master_Comps |> 
          select(Peak, Compound, RI)
        , all = T
        , sort = F) |> 
  # Transform Peak into a factor, set the order of the levels to the inverse of 
  # the ordinal number of the peaks.
  # This is important to ensure the order of the peaks in the plots
  mutate(Peak = factor(Peak
                       , levels = unique(Peak) |> 
                         rev())) |>
  # Merge with data on the correlation of the expression of the genes and the
  # abundance of compounds.
  ## Desaturase genes
  merge(des_comps_cor |> 
          ## Remove columns with the confidence interval limits
          select(-(lower:upper)) |> 
          ## Pivot data frame to have a column per gene, with the measured
          ## correlation 
          pivot_wider(names_from = gene
                      , values_from = estimate)
        , all = T
        , sort = F) |> 
  ## Elongase genes
  merge(elo_comps_cor |> 
          ## Remove columns with the confidence interval limits
          select(-(lower:upper)) |> 
          ## Pivot data frame to have a column per gene, with the measured
          ## correlation 
          pivot_wider(names_from = gene
                      , values_from = estimate)
        , all = T
        , sort = F) |> 
  # Transform into a tibble
  as_tibble()

# NMDS plots ----
## Plot depicting the NMDS coordinate of the samples, as well as their task and 
## subspecies.
nmdsPlot <- {ggplot(data = NMDS, aes(x = NMDS1, y = NMDS2
                                     , fill = Subspecies)) +
    geom_point(aes(shape = Task), alpha = 0.75, size = 2) +
    task_shape_scale +
    subspecies_fill_scale +
    scale_linetype_discrete() +
    coord_equal() +
    nmds_theme}
nmdsPlot

## NMDS plots including the coordinates of the compounds with the strongest
## correlation (positive and negative) to the expression of each gene, as well
## as the magnitude of the measured correlation.
nmds_plus_list <- c("Des", "Elo") |> 
  #  iterating through gene types (desaturases and elongases)
  lapply(function(gene_fam) {
    temp_NMDS_comps <- NMDS_comps |> 
      # Transform into NA any correlation that is not > 0.4
      mutate_at(NMDS_comps |> 
                  select(-(Peak:RI))|> 
                  colnames()
                , function(x)(ifelse(is.na(x)
                                     , x
                                     , ifelse(abs(x) > 0.4
                                              , x
                                              , NA)))) |> 
      # Remove columns corresponding to genes that do not belong to the type of 
      # the current iteration
      select(Peak:RI, contains(gene_fam)) |> 
      # Pivot to have a column stating the gene name and another with the 
      # measured correlation to the compound in the peak of the corresponding 
      # row.
      pivot_longer(-(Peak:RI)
                   , names_to = "gene"
                   , values_to = "rel_exp") |> 
      # Remove NMDS coordinates of peaks containing compounds for which the
      # correlation with the expression of a given gene is NA.
      mutate(NMDS1 = ifelse(is.na(rel_exp)
                            , NA
                            , NMDS1)
             , NMDS2 = ifelse(is.na(rel_exp)
                              , NA
                              , NMDS2))
    
    # Add to the NMDS plots the coordinate of the compounds and the magnitude of
    # their correlation to the expression of the genes.
    p <- nmdsPlot +
      geom_segment(data = temp_NMDS_comps
                   , aes(x = NMDS1
                         , y = NMDS2
                         , xend = 0
                         , yend = 0
                         , color = rel_exp)
                   , linewidth = 1
                   , inherit.aes = F) +
      scale_color_gradient2(low = "#3aa2fcff"
                            , high = "#c82803ff"
                            , mid = "white"
                            , limits = c(-1,1)
                            , name = paste("Spearman's"
                                           , "\ncorrelation"
                            )
                            , guide = 
                              guide_colorbar(title.position = "top"
                                             , barheight = 5)) +
      coord_equal() +
      # Divide the plot into facets, one per gene
      facet_wrap(~gene
                 , ncol = 1) +
      theme(legend.position = "left")
    })

# Name the entries in the list as the corresponding gene types
names(nmds_plus_list) <- c("Des", "Elo")

# Assemble compound plots by binding together the correlograms for the 
# correlation of the individual compounds with expression of the genes and the 
# NMDS plots. Then, export the plots as a PNG.
nmds_plus_list |> 
  # Iterate through gene types (x)
  imap(function(x, y) {
    # Bind together correlogram and NMDS plots
    p <- (wrap_elements(full = 
                          get(paste0(str_to_lower(y)
                                          , "_comps_plot")) + 
                          # labs(title = "A") +
                          theme(legend.position = "none")) |
      x) +
      plot_annotation(tag_levels = "A") +
      plot_layout(guides = "collect"
                  , widths = c(5, 3))
    
    # Generate PNG image with the compound plot
    ggsave(here("figs", paste0(y, "_chc-gexp_plot", ".png"))
           , plot = p
           , dpi = "print"
           , width = 6
           , height = 6.2
           , units = "in"
           , device = agg_png
           , scaling = 0.65)
    })

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "SInf_nmds-script.txt"))

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
