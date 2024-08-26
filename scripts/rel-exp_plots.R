tic("Script for plotting mresults from model inference on gene expression")

# Packages ----
script_packs <- c("ggplot2"
                  , "ggtext"
                  , "purrr"
                  , "forcats"
                  , "gghalves"
                  , "patchwork"
                  , "ragg"
)

## Install and/or load packages to be used by the script
load_my_packs(script_packs)

# ggplot pre-sets ----
source(here("scripts", "custom-ggplot-settings.R"))

# Load the data ----
load(here("data", "processed", "rel-exp_boot-models.Rdata"))

# Plots ----

## Generate plots per gene
plots_list <- models_list |> 
  # Iterate through genes of interest
  imap(function(x, gene) {
    # Extract the data for the plots
    ## raw data
    model_data <- x$data_modeling$model$data
    ## Bootstrapped GLM predictions
    boot_predictions <- x$bootstrap_modeling$boot_predictions
    ## Bootstrapped coefficients
    boot_coefficients <- x$bootstrap_modeling$coefficients |>
      mutate(grp_term = term |> 
               str_remove_all(paste("`"
                                    , "Task"
                                    , "Subspecies"
                                    , sep = "|")) |> 
               str_replace_all("A. m. iberiensis"
                               , "*A. m. iberiensis*")
             , term = term |>
               str_remove_all(paste("`"
                                    , "A. m. iberiensis"
                                    , "Forager bees"
                                    , sep = "|")) |> 
               as.factor())
    
    # Data frame to define reference groups
    labels_df <- data.frame(term = c("Task", "Subspecies")
                            , .estimate = 1
                            , grp_term = c("Nurse bees", "*A. m. carnica*")
                            , reference = T
                            , p_value = NA) |> 
      rbind(data.frame(term = unique(boot_coefficients$term)
                       , .estimate = unique(boot_coefficients$.estimate)
                       , grp_term = unique(boot_coefficients$grp_term)
                       , reference = F
                       , p_value = unique(boot_coefficients$boot_p_value)))
    
    # Plot bootstrapped predictions
    boot_pred_plot <- {model_data |>
        merge(boot_predictions
              , all = T
              , sort = F)  |>
        mutate(p_value = labels_df |> 
                 filter(!is.na(p_value) & term == "Task") |> 
                 pull(p_value)) |>
        as_tibble() |>
        ggplot(aes(x = Subspecies
                   , y = predicted
                   , color =  Task
                   , fill = Task)) +
        geom_point(aes(y = get(gene))
                   , alpha = 0.3
                   , size = 2
                   , position = position_jitterdodge(jitter.width = 0.5
                                                     , jitter.height = 0
                                                     , dodge.width = 0.75
                                                     , seed = 12345)) +
        geom_pointrange(aes(ymin = conf.low
                            , ymax = conf.high)
                        , size = 0.6
                        , linewidth = 1.2
                        , position = position_dodge(width = 0.75)) +
        labs(y = "Relative expresion") +
        coord_flip() +
        task_fill_scale +
        task_color_scale +
        boxplot_theme +
        theme(axis.text.y = element_text(face = "bold.italic")
        )}
    
    # Plot bootstrapped coefficients
    boot_coef_plot <- {boot_coefficients |>
        ggplot(aes(x = term
                   , y = estimate
                   )) +
        geom_half_violin(color = "grey50"
                         , fill = "grey50"
                         , linewidth = 0.2
                         , side = "r"
                         , scale = "width"
                         , nudge = 0.2
                         ) +
        geom_pointrange(aes(ymin = .lower
                            , ymax = .upper
                            , y = .estimate)
                        , position = position_nudge(0.09)
                        , fill = "black"
                        , size = 0.6
                        , linewidth = 1.2) +
        geom_richtext(data = labels_df |> 
                    filter(!reference)
                  , aes(x = term
                        , y = .estimate 
                        , label = grp_term
                        )
                  , fontface = "bold"
                  , size = (boxplot_theme$axis.text$size / .pt) * 0.95
                  , label.color = NA
                  , fill = NA
                  , nudge_x = 0.015) +
        geom_pointrange(data = labels_df |> 
                          filter(reference)
                        , aes(ymin = .estimate
                            , ymax = .estimate
                            , y = .estimate)
                        , position = position_nudge(-0.15)
                        , fill = "black"
                        , size = 0.6
                        , linewidth = 1.2) +
        geom_richtext(data = labels_df |> 
                        filter(reference)
                      , aes(x = term
                          , y = .estimate
                          , label = grp_term)
                      , fontface = "bold"
                      , size = (boxplot_theme$axis.text$size / .pt) * 0.95
                      , label.color = NA
                      , fill = NA
                      , nudge_x = -0.21) +
        geom_hline(yintercept = 1
                   , linewidth = 1
                   , linetype = "dashed"
                   , alpha = 0.75
                   , color = "grey20") +
        geom_text(aes(x = 0.5
                      , y = 1.
                      , label = "Reference")
                  , fontface = "bold"
                  , size = boxplot_theme$axis.text$size / .pt
                  , color = "gray40") +
        boxplot_theme +
        ylim(0, NA) +
        coord_flip() +
        labs(y = "Coefficients"
             , x = NULL) +
        guides(fill = "none") +
        facet_wrap(~fct_rev(term)
                   , scales = "free_y"
                   , nrow = 2) +
        theme(axis.text.y = element_blank()
              , axis.ticks.y = element_blank()
              , legend.position = "none"
              , panel.border = element_rect(fill = NA
                                            , color = "black"
                                            , linetype = "dotted"
                                            , linewidth = 0.4))}
    
    list(pred_plot = boot_pred_plot
         , coef_plot = boot_coef_plot)
  })

## Design matrix for compound plots
design <- "#EE#
            AACC
            AACC
            AACC
            AACC
            BBDD
            BBDD
            BBDD
            BBDD"

## Compound plot for the bootstrapped predictions for the four genes
(plots_list |> 
  lapply(pluck, "pred_plot") |> 
  wrap_plots() + 
  plot_layout(guides = "collect"
              , ncol = 2
              , design = design) +
   guide_area() +
   plot_annotation(tag_levels = list(names(plots_list))) &
   theme(legend.position = "top"
         , legend.background = element_blank()
         , plot.tag.position = c(0.03, 1.01))
 ) |> 
  ggsave(here("figs", "relexp_pred_plot.png")
         , plot = _
         , dpi = "print"
         , width = 6
         , height = 3.8
         , units = "in"
         , device = agg_png
         , scaling = 0.7
         )

## Compound plot for the bootstrapped coefficients for the four genes
(plots_list |> 
  lapply(pluck, "coef_plot") |> 
  wrap_plots() + 
  plot_layout(guides = "collect"
              , ncol = 2
              , byrow = F) +
  plot_annotation(tag_levels = list(names(plots_list))) &
  theme(legend.position = "top"
        , legend.background = element_blank()
        , plot.tag.position = c(0.03, 1.01))) |> 
  ggsave(here("figs", "relexp_coef_plot.png")
         , plot = _
         , dpi = "print"
         , width = 6
         , height = 6
         , units = "in"
         , device = agg_png
         , scaling = 0.65)

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "SInf_rel-exp_tsk&sbsp_plots.txt"))

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
