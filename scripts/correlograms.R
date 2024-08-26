# Packages ----
script_packs <- c("readr"
                  , "purrr"
                  , "magrittr"
                  , "patchwork"
                  , "ggtext"
                  , "gghalves"
                  , "ggh4x"
                  , "ggside"
                  , "GGally"
                  , "ragg"
                  , "ggplot2"
                  , "parsnip"
                  , "broom"
                  , "forcats"
                  , "confintr"
                  )

# Install and/or Load packages
load_my_packs(script_packs)

# ggplot settings ----
source(here("scripts", "custom-ggplot-settings.R"))

# Load the data ----
load(here("data", "processed", "mastertable_data-frames.Rdata"))
load(here("data", "processed", "univariate_chc-stats.Rdata"))
load(here("data", "processed", "CNQRs.Rdata"))

# Shape data frames ----
grouping_info

# merge together data frame with data of compound classes abundance and 
# data frame with data of gene expression
cclass_rel_df <- merge(Prop_CompsClass
                       , qpcr_rel_data
                       , all = T
                       , by = "Individual") |>
  as_tibble()

# merge together data frame with data of mean hydrocarbon chain length and 
# data frame with data of gene expression
cl_rel_df <- merge(Prop_chain.length 
                   , qpcr_rel_data
                   , all = T
                   , by = "Individual") |>
  as_tibble()

# Data frame with information of unsaturated compounds
comps.olefins <- master_Comps |> 
  filter(Class == "Alkene" | Class == "Alkadiene")

# Data frame with the abundance of unsaturated hydrocarbons in all samples
daten.olefins <- master_daten |> 
  select(all_of(comps.olefins$Peak))

# Data frame with abundance of unsaturated hydrocarbons and expression of
# desaturase-encoding genes
eda_olefinsVdes <- merge(cbind(grouping_info
                               , daten.olefins)
                         , qpcr_rel_data |> 
                           select(-starts_with(c("Elo"
                                            , "Rpl")))
                         , by = "Individual") |> 
  as_tibble()

# Data frame with abundance of hydrocarbons and expression of
# elongase-encoding genes
eda_compsVelos <- merge(cbind(grouping_info
                              , master_daten)
                        , qpcr_rel_data |> 
                          select(-starts_with(c("Des"
                                           , "Rpl")))
                        , by = "Individual") |> 
  as_tibble()

# Desaturases-encoding genes vs CHC composition ----
# Correlogram for expression of desaturases-encoding genes and abundance of 
# unsaturated hydrocarbon classes (alkenes and alkadienes)
des_class_plot <- {
  cclass_rel_df |> 
    # Only numeric variables remain
    select(where(is.numeric)) |> 
    # Remove expression of elongases and reference genes
    select(-starts_with(c("El", "Rpl"))) |> 
    # Re-order columns
    select(Alkane
           , everything()) |> 
    # remove rows with NAs
    drop_na() |> 
    # Function to calculate the Spearman's correlation index between
    # the expression of the genes and the abundance of the compound classes
    (function(x) {
      x |> 
        select(starts_with("Des")) |> 
        imap(function(.x, id) {
          print(paste("gene:", id))
          x |> 
            select(-starts_with("Des")) |> 
            imap(function(y, id_y) {
              print(paste("class:", id_y))
              tmp <- ci_cor(x = .x
                            , y = y
                            , method = "spearman"
                            , type = "bootstrap"
                            , seed = 12345)
              
              tmp |> 
                pluck("estimate") |> 
                cbind.data.frame(tmp |> 
                                   pluck("interval") |> 
                                   as.data.frame() |> 
                                   t() |> 
                                   as.data.frame()) |> 
                set_colnames(c("estimate", "lower", "upper")) |> 
                set_rownames(NULL)
            })
        })
    })() |> 
    list_transpose() |> 
    # keep only correlation between expression of genes and abundance of 
    # alkenes or alkadienes
    keep_at(c("Alkene", "Alkadiene")) |> 
    # Set gene names as a variable (column) of the data frames in the list
    imap(function(x = .x, y = .y){
      x |> 
        mutate(Class = y) |> 
        rownames_to_column("gene")
    }) |> 
    # Merge all data frames in the list into one
    reduce(merge, all = T) |> 
    as_tibble() |> 
    # Plot correlogram
    ggplot(aes(y = Class
               , x = gene)) +
    geom_raster(aes(fill = estimate)) +
    geom_text(aes(label = estimate |> 
                    round(3) |> 
                    format(digits = 2
                           , nsmall = 3) |> 
                  paste(paste0("("
                               , lower |> 
                                 round(3) |> 
                                 format(digits = 2
                                        , nsmall = 3)
                               , ", "
                               , upper |> 
                                 round(3) |> 
                                 format(digits = 2
                                        , nsmall = 3)
                               , ")"))
                  )
              ) +
    scale_fill_gradient2(low = "#3aa2fcff"
                         , high = "#c82803ff"
                         , mid = "white"
                         , midpoint = 0
                         , limits = c(-1,1)
                         , name = "Spearman's \ncorrelation") +
    labs(y = "Unsaturated CHCs' classes"
         , x = "Desaturase genes"
    ) +
    theme_classic() +
    theme(axis.text = element_text(face = "bold"))
}
# ggsave(here("figs", "des_cor_class.png")
#        , plot = des_class_plot 
#        , dpi = "print"
#        , width = 5
#        , height = 3
#        , units = "in"
#        , device = agg_png
#        , scaling = 0.7)

# Speraman's correlation of desaturases-encoding genes and abundance of
# unsaturated hydrocarbons
des_comps_cor <- {
  eda_olefinsVdes |> 
    # Only numeric variables remain
    select(where(is.numeric)) |> 
    # remove rows with NAs
    drop_na() |>
    # Function to calculate the Spearman's correlation index between
    # the expression of the genes and the abundance of the compounds
    (function(x) {
      x |> 
        select(starts_with("Des")) |> 
        imap(function(.x, id) {
          print(paste("gene:", id))
          x |> 
            select(-starts_with("Des")) |> 
            imap(function(y, id_y) {
              print(paste("peak:", id_y))
              tmp <- ci_cor(x = .x
                            , y = y
                            , method = "spearman"
                            , type = "bootstrap"
                            , seed = 12345)
              
              tmp |> 
                pluck("estimate") |> 
                cbind.data.frame(tmp |> 
                                   pluck("interval") |> 
                                   as.data.frame() |> 
                                   t() |> 
                                   as.data.frame()) |> 
                set_colnames(c("estimate", "lower", "upper")) |> 
                set_rownames(NULL)
            })
        })
    })() |> 
    list_transpose() |> 
    # Set gene names as a variable (column) of the data frames in the list
    imap(function(x = .x, y = .y){
      x |> 
        mutate(Peak = y) |> 
        rownames_to_column("gene")
    }) |> 
    # Merge all data frames in the list into one
    reduce(merge, all = T, sort = F) |> 
    # Peak as factor, ensure correct order of levels by peak number
    mutate(Peak = factor(Peak
                         , levels = unique(Peak) |> 
                           rev())
           , gene = as.factor(gene)) |> 
    as_tibble()
}

# Correlogram for expression of desaturases-encoding genes and abundance of 
# unsaturated hydrocarbons
des_comps_plot <- {
  des_comps_cor |> 
    # plot correlogram
    (function(df) {
      df |> 
        ggplot(aes(y = Peak
                   , x = gene)) +
        geom_raster(aes(fill = estimate)) +
        scale_fill_gradient2(low = "#3aa2fcff"
                               , high = "#c82803ff"
                               , mid = "white"
                               , limits = c(-1,1)
                             , name = "Spearman's \ncorrelation") +
        geom_text(aes(label = estimate |> 
                        round(3) |> 
                        format(digits = 2
                               , nsmall = 3) |> 
                        paste(paste0("("
                                     , lower |> 
                                       round(3) |> 
                                       format(digits = 2
                                              , nsmall = 3)
                                     , ", "
                                     , upper |> 
                                       round(3) |> 
                                       format(digits = 2
                                              , nsmall = 3)
                                     , ")"))
                      )
                  , size = 2.5
        ) +
        scale_y_discrete(labels = master_Comps |>
                           filter(Peak %in% df$Peak) |>
                           arrange(desc(RI)) |>
                           pull(Compound)) +
        labs(y = "Unsaturated CHCs"
             , x = "Desaturase genes"
             ) +
        theme_classic() +
        theme(axis.text = element_text(face = "bold"))
    })()
}
# ggsave(here("figs", "des_cor_comps.png")
#        , plot = des_comps_plot 
#        , dpi = "print"
#        , width = 4
#        , height = 6
#        , units = "in"
#        , device = agg_png
#        , scaling = 0.65)

# Elongases-encoding genes vs CHC composition ----
# Speraman's correlation of elongases-encoding genes and mean hydrocarbon
# chain length
elo_cl_cor <- {
  cl_rel_df |> 
    # Only numeric variables remain
    select(where(is.numeric)) |> 
    # Remove expression of elongases and reference genes
    select(-starts_with(c("Des", "Rpl"))) |> 
    # remove rows with NAs
    drop_na() |>  
    # Function to calculate the Spearman's correlation index between
    # the expression of the genes and the mean hydrocarbon chain length
    (function(x) {
      x |> 
        select(starts_with("El")) |> 
        imap(function(.x, id) {
          print(paste("gene:", id))
          x |> 
            select(-starts_with("El")) |> 
            imap(function(y, id_y) {
              print(paste("var:", id_y))
              tmp <- ci_cor(x = .x
                            , y = y
                            , method = "spearman"
                            , type = "bootstrap"
                            , seed = 12345)
              
              tmp |> 
                pluck("estimate") |> 
                cbind.data.frame(tmp |> 
                                   pluck("interval") |> 
                                   as.data.frame() |> 
                                   t() |> 
                                   as.data.frame()) |> 
                set_colnames(c("estimate", "lower", "upper")) |> 
                set_rownames(NULL)
            })
        })
    })() |> 
    list_transpose() |> 
    # Set gene names as a variable (column) of the data frames in the list
    imap(function(x = .x, y = .y){
      x |> 
        mutate(var = y) |> 
        rownames_to_column("gene")
    }) |> 
    # Merge all data frames in the list into one
    reduce(merge, all = T) |> 
    as_tibble()
}

# Correlogram for expression of elongases-encoding genes and mean hydrocarbon
# chain length
elo_cl_plot <- {elo_cl_cor |> 
    # Plot correlogram
  ggplot(aes(x = gene, y = var)) +
  geom_raster(aes(fill = estimate)) +
  scale_fill_gradient2(low = "#3aa2fcff"
                       , high = "#c82803ff"
                       , mid = "white"
                       , limits = c(-1,1)
                       , name = "Spearman's \ncorrelation") +
  geom_text(aes(label = estimate |> 
                  round(3) |> 
                  format(digits = 2
                         , nsmall = 3) |> 
                  paste(paste0("("
                               , lower |> 
                                 round(3) |> 
                                 format(digits = 2
                                        , nsmall = 3)
                               , ", "
                               , upper |> 
                                 round(3) |> 
                                 format(digits = 2
                                        , nsmall = 3)
                               , ")"))
                  )
  ) +
  labs(x = "Elongase genes"
       , y = "CHCs' mean chain length") +
  theme_classic() +
  theme(axis.text.y = element_blank()
        , axis.ticks.y = element_blank()
        , axis.line.y = element_blank()
        , axis.text.x = element_text(face = "bold"))}
# ggsave(here("figs", "elo_cor_cl.png")
#        , plot = elo_cl_plot 
#        , dpi = "print"
#        , width = 4
#        , height = 2
#        , units = "in"
#        , device = agg_png
#        , scaling = 0.7)


# Speraman's correlation of elongases-encoding genes and abundance of
# hydrocarbons
elo_comps_cor <- {eda_compsVelos |> 
    # Only numeric variables remain
    select(where(is.numeric)) |> 
    # remove rows with NAs
    drop_na() |>
    # Function to calculate the Spearman's correlation index between
    # the expression of the genes and the abundance of the compounds
    (function(x) {
      x |> 
        select(starts_with("El")) |> 
        imap(function(.x, id) {
          print(paste("gene:", id))
          x |> 
            select(-starts_with("El")) |> 
            imap(function(y, id_y) {
              print(paste("peak:", id_y))
              tmp <- ci_cor(x = .x
                            , y = y
                            , method = "spearman"
                            , type = "bootstrap"
                            , seed = 12345)
              
              tmp |> 
                pluck("estimate") |> 
                cbind.data.frame(tmp |> 
                                   pluck("interval") |> 
                                   as.data.frame() |> 
                                   t() |> 
                                   as.data.frame()) |> 
                set_colnames(c("estimate", "lower", "upper")) |> 
                set_rownames(NULL)
            })
        })
    })() |> 
    list_transpose() |> 
    # Set gene names as a variable (column) of the data frames in the list
    imap(function(x = .x, y = .y){
      x |> 
        mutate(Peak = y) |> 
        rownames_to_column("gene")
    }) |> 
    # Merge all data frames in the list into one
    reduce(merge, all = T, sort = F) |> 
    # Peak as factor, ensure correct order of levels by peak number
    mutate(Peak = factor(Peak
                         , levels = unique(Peak) |> 
                           rev())
           , gene = as.factor(gene)) |> 
    as_tibble()
}

# Correlogram for expression of elongases-encoding genes and abundance of
# hydrocarbons
elo_comps_plot <- {
  elo_comps_cor |>
    # plot correlogram
    (function(df){
      df |> 
        ggplot(aes(y = Peak
                   , x = gene)) +
        geom_raster(aes(fill = estimate)) +
        scale_fill_gradient2(low = "#3aa2fcff"
                             , high = "#c82803ff"
                             , mid = "white"
                             , limits = c(-1,1)
                             , name = "Spearman's \ncorrelation") +
        geom_text(aes(label = estimate |> 
                        round(3) |> 
                        format(digits = 2
                               , nsmall = 3) |> 
                        paste(paste0("("
                                     , lower |> 
                                       round(3) |> 
                                       format(digits = 2
                                              , nsmall = 3)
                                     , ", "
                                     , upper |> 
                                       round(3) |> 
                                       format(digits = 2
                                              , nsmall = 3)
                                     , ")"))
                        )
                  , size = 2.5) +
        scale_y_discrete(labels = master_Comps |> 
                           filter(Peak %in% df$Peak) |> 
                           arrange(desc(RI)) |> 
                           pull(Compound)) +
        labs(y = "Cuticular hydrocarbons"
             , x = "Elongase genes"
             ) +
        theme_classic() +
        theme(axis.text = element_text(face = "bold"))
    })()
}
# ggsave(here("figs", "elo_cor_comps.png")
#        , plot = elo_comps_plot 
#        , dpi = "print"
#        , width = 4
#        , height = 6.5
#        , units = "in"
#        , device = agg_png
#        , scaling = 0.65)

# Compound plot ----
# Design matrix for compound plot
design <- "AA
          AA
          AA
          AA
          BB
          BB
          BB"

# Compound plot with correlograms of gene expression vs abundance of alkenes 
# and alkadienes, and gene expression vs mean hydrocarbon chain length
des_elo_plot <- (des_class_plot / free(elo_cl_plot)) +
  plot_layout(guides = "collect"
              , design = design) +
  plot_annotation(tag_levels = "A")

ggsave(here("figs", "des&elo_cor-plot.png")
       , plot = des_elo_plot
       , dpi = "print"
       , width = 6
       , height = 4.8
       , units = "in"
       , device = agg_png
       , scaling = 0.8)

# Export data ----
save(list = c("des_comps_cor"
              ,  "des_comps_plot"
              , "elo_cl_cor"
              , "elo_comps_cor"
              , "elo_comps_plot")
     , file = here("data"
                   , "processed"
                   ,  "chc-v-genes_correlations.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_correlograms-script.txt"))

## Detach/unload packages
pacman::p_unload(char = script_packs)

## Clear environment
objects_list <- ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|"))
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

gc()
