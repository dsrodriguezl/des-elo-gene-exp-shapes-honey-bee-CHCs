# Script to generate custom ggplot settings

# Themes ----
## Theme for cluster analysis dendrogram plots ----
cluster_theme <- theme(legend.position = "left"
                       , legend.key = element_rect(fill = "white"
                                                   , color = "white")
                       , legend.key.size = unit(3.5, "mm")
                       , legend.spacing.y = unit(0.2, "lines")
                       , legend.title = element_text(size = 10)
                       , legend.text = element_text(size = 8)
                       , panel.background = element_rect(fill = "white")
                       , plot.margin = unit(c(1, 10, 1, 1), "lines")
                       , axis.text.x = element_text(size = 7)
                       , axis.text.y = element_blank()
                       , axis.ticks = element_blank()
                       , axis.title = element_blank())

dispersion_theme <- theme_classic() +
  theme(strip.text.x = element_text(face = "bold"
                                    , size = 10)
        , strip.text.y = element_markdown(face = "bold"
                                          , size = 10)
        , strip.background = element_blank()
        #, legend.position = "left"
        , legend.key.size = unit(6, "mm")
        , legend.spacing.y = unit(0.3
                                  , "lines")
        , legend.title = element_text(size = 10)
        , panel.border = element_rect(fill = NA
                                      , color = "black"
                                        , linetype = "dotted"
                                      , linewidth = 0.1)
        , ggside.panel.scale = 0.25
        , ggside.axis.text = element_blank()
        , ggside.axis.ticks = element_blank()
  )

## Theme for NMDS plots ----
nmds_theme <- theme_classic() +
  theme(strip.text.x = element_text(face = "bold"
                                  , size = 10)
        , strip.text.y = element_markdown(face = "bold"
                                             , size = 10)
        , strip.background = element_blank()
        , axis.text = element_text(face = "bold")
        , legend.key.size = unit(3.5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , panel.border = element_rect(fill = NA
                                      , color = "black"
                                      , linetype = "dotted"
                                      , linewidth = 0.5))

## Theme for box plots ----
boxplot_theme <- theme_classic() +
  theme(title = element_text(face = "bold"
                             , size = 11)
        , strip.text = element_text(face = "bold"
                                    , size = 10)
        , strip.background = element_blank()
        , legend.position = "top"
        # , legend.justification = "center"
        , legend.key.size = unit(7, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(face = "bold"
                                      , size = 10)
        , legend.text = element_text(size = 10)
        , axis.text = element_text(face = "bold"
                                   , size = 9)
        , axis.title = element_text(face = "bold"
                                    , size = 10)
        # , axis.title.y = element_blank()
        # , axis.text.y = element_blank()
        # , axis.ticks.y = element_blank()
        )

## ----
chc_richness_theme <- theme_classic() +
  theme(legend.key.size = unit(3.5, "mm")
        , legend.spacing.y = unit(0.2, "lines")
        , legend.title = element_text(size = 10)
        , legend.text = element_text(size = 8)
        , axis.text = element_text(size = 8)
        , axis.title = element_text(size = 10))

# Scales ----
## colors ----
### Tasks
task_fill_scale <- 
  scale_fill_manual(values = c("Nurses" = "#5D177FFF"
                               , "Pollen foragers" = "#D1426FFF"
                               , "Nurse bees" = "#5D177FFF"
                               , "Forager bees" = "#D1426FFF"
                               , "Foragers" = "#D1426FFF")
                    , guide = guide_legend(title = "Task"
                                           , override.aes = list(shape = 21)))

task_color_scale <- 
  scale_color_manual(values = c("Nurses" = "#5D177FFF"
                                , "Pollen foragers" = "#D1426FFF"
                                , "Nurse bees" = "#5D177FFF"
                                , "Forager bees" = "#D1426FFF"
                                , "Foragers" = "#D1426FFF")
                     , guide = guide_legend(title = "Task"))

### Subspecies
subspecies_fill_scale <-
  scale_fill_manual(values = c("A. m. carnica" = "#F1E51DFF"
                               , "A. m. iberiensis" = "#46085CFF")
                    , guide = guide_legend(title = "Subspecies"
                                           , label.theme = 
                                             element_text(face = "italic"
                                                          , size = 8
                                                          )
                                           , byrow = T
                                           # Force the shape in the legend to
                                           # include filling, so it can show
                                           # the color of the fill
                                           , override.aes = list(shape = 21)))

subspecies_color_scale <- 
  scale_color_manual(values = c("A. m. carnica" = "#F1E51DFF"
                                , "A. m. iberiensis" = "#46085CFF")
                     , guide = guide_legend(title = "Subspecies"
                                            , label.theme = 
                                              element_text(face = "italic"
                                                           , size = 8
                                                           )
                                            , byrow = T))

## shapes ----
task_shape_scale <- 
  scale_shape_manual(values = c("Pollen foragers" = 21
                                , "Nurses" = 24
                                , "Forager bees" = 21
                                , "Foragers" = 21
                                , "Nurse bees" = 24)
                     , name = "Task"
                     , guide = guide_legend(override.aes =
                                              list(color = "black")
                                            , label.theme = 
                                              element_text(size = 8)))


# plotting functions ----
my_lp_plot <- function(ggplot_object, jitter_height = 0) {
  p <- ggplot_object +
    geom_half_violin(color = "grey60"
                     , linewidth = 0.2
                     , alpha = 0.55
                     , side = "r"
                     , scale = "width"
                     , nudge = 0.06) +
    geom_half_boxplot(color = "grey35"
                      , size = 1
                      , alpha = 0.9
                      , side = "r"
                      #, center = T
                      , outlier.shape = NA
                      , errorbar.draw = F
                      , nudge = 0.02) +
    geom_half_point(color = "grey40"
                    , side = "l"
                    , shape = 21
                    , size = 1.5
                    , stroke = 0.8
                    # , alpha = 0.8
                    , range_scale = 0.45
                    , transformation = position_jitter(seed = 12345
                                                       , width = 0.04
                                                       , height = jitter_height)
                    ) +
    boxplot_theme +
    coord_flip()
  
  p
}

custom_ggs <- ls() |> 
  str_remove_all(analysis_objects |> 
                   rev() |> 
                   paste(collapse = "|")) |> 
  str_remove_all("script_packs") |> 
  str_subset("[:alpha:]")


