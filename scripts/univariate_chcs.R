# Vector holding the list of packages that will be used
script_packs <- c("forcats"
                  , "patchwork"
                  , "ggtext"
                  , "gghalves"
                  , "ggh4x"
                  , "ragg"
                  , "ggplot2"
                  , "tidymodels"
                  , "emmeans"
                  , "ggeffects"
                  , "quantreg")

# Install and/or Load packages
load_my_packs(script_packs)

# FUNCTIONS ----
## In this section the functions that will be used latter in the script are defined
## Without running the code on this section in advance, some parts of the script
## will not work as they rely on these functions

## Data frame by hydrocarbon classes ----
cclasses_df <- function(daten, grouping.info, Comps.data){
  # Extract unique compound classes from Comps.data
  cclasses <- unique(Comps.data$Class)
  
  # Create a list of data frames for each class
  # Each data frame contains the rows from daten where Comps.data$Class == c.class
  class_dfs <- lapply(cclasses, function(c.class) {
    as.data.frame(t(daten)) |>
      filter(Comps.data$Class == c.class)
  })
  
  # Assign names to list elements
  names(class_dfs) <- cclasses
  
  # Create a data frame with sample names and ranks
  Sampls_index <- data.frame('Sample name' = rownames(daten),
                             rank = rank(grouping.info$Individual),
                             stringsAsFactors = FALSE)
  
  # Create a data frame with the sum of columns for each class
  # The column names are taken from the names of the list elements (i.e. the compound classes)
  Prop.CompsClass <- data.frame(grouping.info,
                                sapply(class_dfs, function(df) colSums(df)),
                                row.names = Sampls_index$Sample.name)
  
  # Convert to tibble and return result
  Prop.CompsClass |> as_tibble()
}

## Chain length data set ----
## Function to comps_classrate a data set with the weighted average chain length per
## individual
cl_df <- function(daten, Comps.data, grouping.info) {
  # Transpose daten and convert to data frame
  daten_t <- as.data.frame(t(daten))
  
  # Calculate sum of abundance for each unique chain length
  abundance_sums <- sapply(unique(Comps.data$Chain.length)
                           , function(chain_length) {
                             daten_t |> 
                               filter(Comps.data$Chain.length == 
                                        chain_length) |>
                               colSums() |> 
                               as.numeric()
                             })
  
  # Convert abundance_sums to data frame and set row names
  cl.df <- as.data.frame(abundance_sums)
  rownames(cl.df) <- rownames(daten)
  
  # Set column names of cl.df to chain lengths
  colnames(cl.df) <- unique(Comps.data$Chain.length)
  
  # Merge grouping.info with cl.df
  cl.df <- cbind(grouping.info, cl.df)
  
  # Calculate weighted mean chain length per individual and return result as tibble
  result <- merge(grouping.info,
                  cl.df |> 
                    pivot_longer(cols = !where(is.factor) & !Individual,
                                 names_to = "Chain.length",
                                 values_to = "Abundance") |>  
                    mutate(Chain.length = as.numeric(Chain.length)) |>
                    group_by(Individual) |>
                    summarise(cl_w.mean = weighted.mean(Chain.length
                                                        , Abundance))
                  , by = "Individual") |> 
    as_tibble()
  
  result
}

# ggplot themes ----
source(here("scripts", "custom-ggplot-settings.R"))

# Load the data frames ----
load(here("data", "processed", "mastertable_data-frames.Rdata"))

# Hydrocarbon classes' abundance per group ----
## Creates a data frame with the abundance of each compounds class per sample 
## and the factors to group the samples
Prop_CompsClass <- cclasses_df(master_daten
                               , grouping_info |> 
                                 select_if(function(col)(!is.logical(col)))
                               , master_Comps)  |> 
  mutate_if(is.factor, fct_drop)

# Vector listing the name of mono- and di-unsaturated hydrocarbon classes
comps_classes <- Prop_CompsClass |> 
  select(Alkene, Alkadiene) |> 
  colnames()

## Fit models ----
# empty list for holding model results
cclasses_models_list <- list()

# PDF for the diagnostic plots of the models
pdf(here("output"
         , paste0("model-diagnotsics_chc-classes"
                  , '.pdf'))
    , width = 10, height = 8)

# Fit quantile regression models for the abundance of alkenes and alkadienes
cclasses_models_list <- lapply(comps_classes, function(comps_class) {
  
  # Define model formula
  model_formula <- paste(comps_class
                         , "~"
                         , paste("Subspecies"
                                 , "Task"
                                 # Uses an interaction model for the alkenes 
                                 # and anadditive model for the alkadienes
                                 , sep = ifelse(comps_class == "Alkene"
                                                , " * "
                                                ,  " + ")))
  
  # Subset data frame to not include the other comps_classs
  # This just to make the final object lighter
  model_data <- Prop_CompsClass |>
    select(contains(comps_class), Individual, Subspecies:Task)
  
  # Fit the quantile regression model for the 50% quantile
  model <- rq(as.formula(model_formula), data = model_data, tau = 0.5)
  
  # Print diagnostic plots
  plot(c(0, 1)
       , c(0, 1)
       , ann = F
       , bty = 'n'
       , type = 'n'
       , xaxt = 'n'
       , yaxt = 'n')
  text(x = 0.5, y = 0.5
       , paste("rq - tau = 0.5:", model_formula)
       , cex = 1.6, col = "black")
  
  # Inspect interaction
  model |> 
    emmip(Task ~ Subspecies) |> 
    print()
  
  # Get estimated marginal means from the GLM
  pairwise_emmeans <- model |> 
    emmeans(pairwise ~ Task * Subspecies, type = "response") |> 
    map(tidy, conf.int = T)
  
  # Plot p values per pairwise contrast
  p_plot <- pairwise_emmeans$contrasts |> 
    ggplot(aes(x = contrast
               , y = adj.p.value)) +
    geom_col(color = "lightblue2"
             , fill = "lightblue2") +
    geom_text(aes(label = adj.p.value |> 
                    format(digits = 3
                           , nsmall = 2)
                  , fontface = "bold")
              , angle = 270
              , color = "grey20") +
    theme_classic() +
    coord_flip() +
    geom_hline(aes(yintercept = 0.05)
               , linetype = "dashed", linewidth = 1) +
    labs(title = "P values of model's pairwise contrasts"
         , subtitle = comps_class)
    print(p_plot)
    
    # Get a preliminary prediction plot
    model |> 
      ggpredict(c("Subspecies", "Task")) |> 
      plot() |> 
      print()
    
    prediction_plot <- pairwise_emmeans |> 
      pluck(1) |> 
      select(Task:estimate, starts_with("conf.")) |> 
      merge(model_data
            , all = T
            , sort = F) |> 
      mutate(Task = Task |> as.factor() |> fct_relevel("Nurse bees")
             , Subspecies = Subspecies |> as.factor()) |> 
      as_tibble() |> 
      ggplot(aes(x = Subspecies
                 , y = get(comps_class)
                 , color = Task
                 , fill = Task)) +
      geom_point(position = position_jitterdodge(seed = 12345
                                                 , jitter.width = 0.3
                                                 , jitter.height = 0
                                                 , dodge.width = 0.75)
                 , alpha = 0.3
                 , size = 1.8) +
      geom_pointrange(aes(y = estimate
                          , ymin = conf.low
                          , ymax = conf.high)
                      , size = 0.6
                      , linewidth = 1.2
                      , position = position_dodge(width = 0.75)) +
      task_color_scale +
      task_fill_scale +
      boxplot_theme +
      coord_flip() +
      labs(y = paste("Relative abundance (%) of"
                     , paste0(str_to_lower(comps_class), "s"))) +
      theme(axis.text.y = element_text(face = "bold.italic"))
  
    print(prediction_plot)
  
    # Store objects on a list for the corresponding compound class
    list(model = model
         , emmeans = pairwise_emmeans
         , tidy = tidy(model, conf.int = T)
         , glance = glance(model)
         , prediction_plot = prediction_plot)
})
dev.off()

# Names of entries as the included compound classes
names(cclasses_models_list) <- comps_classes

# Get plots for the model coefficients and store them in a list
coef_plots_list <- cclasses_models_list |>
  # iterate through compound classes
  imap(function(x, cclass){
    cclass_tidy <- x$tidy
    
    # Appropiately name groups
    if (cclass_tidy$term |> 
       str_detect(":") |> 
       any()) {
      grp_term <- c("Nurse bees", "A. m. carnica", "A. m. carnica:Nurse bees")
      terms <- c("Task", "Subspecies", "Subspecies:Task")
    } else {
      grp_term <- c("Nurse bees", "A. m. carnica") 
      terms <- c("Task", "Subspecies")
    }
    
    # Plot model coefficients
     coeff_plot <-  {cclass_tidy |> 
         filter(term != "(Intercept)") |> 
         mutate(grp_term = term |> 
                  str_remove_all("Subspecies|Task")
                , term = term |> 
                  str_remove_all("A. m. iberiensis|Forager bees")
                , reference = F) |> 
         merge(data.frame(term = terms
                          , grp_term = grp_term
                          , estimate = 0
                          , reference = T)
               , all = T
               , sort = F) |>
         mutate(grp_term = grp_term |> 
                  str_replace_all("A. m. iberiensis"
                                  , "*A. m. iberiensis*") |> 
                  str_replace_all("A. m. carnica"
                                  , "*A. m. carnica*") |> 
                  as.factor() |> 
                  fct_relevel(c("*A. m. carnica*"
                                , "Nurse bees"
                                , "*A. m. carnica*:Nurse bees"))) |> 
         as_tibble() |> 
         ggplot(aes(x = grp_term
                    , y = estimate)) +
         geom_pointrange(aes(ymin = conf.low
                             , ymax = conf.high)
                         , fill = "black"
                         , size = 0.6
                         , linewidth = 1.2) +
         geom_hline(yintercept = 0
                    , linewidth = 1
                    , linetype = "dashed"
                    , alpha = 0.8
                    , color = "grey25") +
         geom_richtext(aes(x = grp_term
                           , y = estimate
                           , label = grp_term)
                       , fontface = "bold"
                       , size = (boxplot_theme$axis.text$size / .pt) * 0.95
                       , label.color = NA
                       , fill = NA
                       , nudge_x = 0.3
                       ) +
         geom_text(aes(x = ifelse(cclass == "Alkadiene", 0.6, 0.7)
                       , y = 0
                       , label = "Reference")
                   , fontface = "bold"
                   , size = boxplot_theme$axis.text$size / .pt
                   , color = "gray40") +
         facet_wrap(~term
                    , scales = "free_y"
                    , ncol = 1)  +
         labs(x = NULL
              , y = "Coefficients"
              , subtitle =  paste("Relative abundance of"
                                  , paste0(str_to_lower(cclass), "s"))) +
         coord_flip() +
         boxplot_theme + 
         theme(axis.text.y =  element_blank()
               , axis.ticks.y = element_blank()
               , panel.border = element_rect(fill = NA
                                             , color = "black"
                                             , linetype = "dotted"
                                             , linewidth = 0.3))
   }
  
  })

# Compound plot for model predictions
unsat_pred_plot <- cclasses_models_list |> 
  lapply(pluck, "prediction_plot") |> 
  wrap_plots(ncol = 1
             , guides = "collect") &
  theme(legend.position = "top")

# Compound plot for model coefficients
unsat_coef_plot <- coef_plots_list |> 
  wrap_plots(ncol = 1)

# Mean chain length ----
## Prepare the data frame
Prop_chain.length <- cl_df(master_daten
                           , master_Comps
                           , grouping_info |> 
                             select_if(function(col)(!is.logical(col))))

##  Fit models ----
# PDF for the diagnostic plots of the model
pdf(here("output", "model-diagnotsics_chc-mean-chain-length.pdf")
    , width = 10, height = 8)
{
  # Define model formula
  model_formula <- paste("cl_w.mean"
                         , "~"
                         , paste("Subspecies"
                                 , "Task"
                                 , sep = " + "))
  
  mean_cl_model <- rq(as.formula(model_formula)
                      , data = Prop_chain.length
                      , tau = 0.5)
  
  # Print diagnostic plots
  plot(c(0, 1)
       , c(0, 1)
       , ann = F
       , bty = 'n'
       , type = 'n'
       , xaxt = 'n'
       , yaxt = 'n')
  text(x = 0.5, y = 0.5
       , paste("rq - tau = 0.5:", model_formula)
       , cex = 1.6, col = "black")
  
  # Inspect interaction
  mean_cl_model |> 
    emmip(Task ~ Subspecies) |> 
    print()
  
  # Get estimated marginal means from the model
  pairwise_emmeans <- mean_cl_model |> 
    emmeans(pairwise ~ Task * Subspecies, type = "response") |> 
    map(tidy, conf.int = T)

  # Plot p values per pairwise contrast
  p_plot <- pairwise_emmeans$contrasts |> 
    ggplot(aes(x = contrast
               , y = adj.p.value)) +
    geom_col(color = "lightblue2"
             , fill = "lightblue2") +
    geom_text(aes(label = adj.p.value |> 
                    format(digits = 3
                           , nsmall = 2)
                  , fontface = "bold")
              , angle = 270
              , color = "grey20") +
    theme_classic() +
    coord_flip() +
    geom_hline(aes(yintercept = 0.05)
               , linetype = "dashed", linewidth = 1) +
    labs(title = "P values of model's pairwise contrasts"
         , subtitle = "Weighted mean chain length")
  print(p_plot)
  
  # Get a preliminary prediction plot
  mean_cl_model |> 
    ggpredict(c("Subspecies", "Task")) |> 
    plot() |> 
    print()
  
  # Plot model predictions
  prediction_plot <- pairwise_emmeans |> 
    pluck(1) |> 
    select(Task:estimate, starts_with("conf.")) |> 
    merge(Prop_chain.length
          , all = T
          , sort = F) |>  
    mutate(Task = Task |> as.factor() |> fct_relevel("Nurse bees")
           , Subspecies = Subspecies |> as.factor()) |> 
    as_tibble() |> 
    ggplot(aes(x = Subspecies
               , y = cl_w.mean
               , color = Task
               , fill = Task)) +
    geom_point(position = position_jitterdodge(seed = 12345
                                               , jitter.width = 0.3
                                               , jitter.height = 0
                                               , dodge.width = 0.75)
               , alpha = 0.3
               , size = 2) +
    geom_pointrange(aes(y = estimate
                        , ymin = conf.low
                        , ymax = conf.high)
                    , size = 0.6
                    , linewidth = 1.2
                    , position = position_dodge(width = 0.75)) +
    task_color_scale +
    task_fill_scale +
    boxplot_theme +
    coord_flip() +
    labs(y = "Weighted mean chain length of hydrocarbons") +
    theme(axis.text.y = element_text(face = "bold.italic"))
  print(prediction_plot)
}
dev.off()

# Plot model coefficients
coeff_plot <- mean_cl_model |>
  (function(x){
    tidy_model <- tidy(x)
    
    if (tidy_model$term |> 
        str_detect(":") |> 
        any()) {
      grp_term <- c("A. m. carnica", "Nurse bees", "Nurse bees:A. m. carnica")
      terms <- c("Subspecies", "Task", "Task:Subspecies")
      } else {
        grp_term <- c("Nurse bees", "A. m. carnica") 
        terms <- c("Task", "Subspecies")
      }
    
    tidy_model |> 
      filter(term != "(Intercept)") |> 
      mutate(grp_term = term |> 
               str_remove_all("Subspecies|Task")
             , term = term |> 
               str_remove_all("A. m. iberiensis|Forager bees")
             , reference = F) |> 
      merge(data.frame(term = terms
                       , grp_term = grp_term
                       , estimate = 0
                       , reference = T)
            , all = T
            , sort = F) |>
      mutate(grp_term = grp_term |> 
               str_replace_all("A. m. iberiensis"
                               , "*A. m. iberiensis*") |> 
               str_replace_all("A. m. carnica"
                               , "*A. m. carnica*") |> 
               as.factor() |> 
               fct_relevel(c("*A. m. carnica*"
                             , "Nurse bees"
                             , "*A. m. carnica*:Nurse bees"))) |> 
      as_tibble() |> 
      ggplot(aes(x = grp_term
                 , y = estimate)) +
      geom_pointrange(aes(ymin = conf.low
                          , ymax = conf.high)
                      , size = 0.6
                      , linewidth = 1.2) +
      geom_hline(yintercept = 0
                 , linewidth = 1
                 , linetype = "dashed"
                 , alpha = 0.8
                 , color = "grey25") +
      geom_richtext(aes(x = grp_term
                        , y = estimate #- 0.05
                        , label = grp_term)
                    , fontface = "bold"
                    , size = (boxplot_theme$axis.text$size / .pt) * 0.95
                    , label.color = NA
                    , fill = NA
                    # , color = "gray20"
                    , nudge_x = 0.3) +
      geom_text(aes(x = 0.6
                    , y = 0
                    , label = "Reference")
                , fontface = "bold"
                , size = boxplot_theme$axis.text$size / .pt
                # , angle = 90
                , color = "gray40") +
      facet_wrap(~term
                 , scales = "free_y"
                 , ncol = 1)  +
      labs(x = NULL
           , y = "Coefficients"
           , subtitle = "Weighted mean chain length of hydrocarbons") +
      coord_flip() +
      boxplot_theme + 
      theme(axis.text.y =  element_blank()
            , axis.ticks.y = element_blank()
            , panel.border = element_rect(fill = NA
                                          , color = "black"
                                          , linetype = "dotted"
                                          , linewidth = 0.3))
})()

# Compound plot with the model predictions for the abundance of alkenes and 
# alkadienes, and the mean chain length
predictions_plot <- unsat_pred_plot / 
  (prediction_plot + 
     theme(legend.position = "none")) +
  plot_annotation(tag_levels = "A")
ggsave(here("figs", "chc-traits_plot.png")
       , plot = predictions_plot
       , dpi = "print"
       , width = 5.5
       , height = 6
       , units = "in"
       , device = agg_png
       , scaling = 0.65)

# Compound plot with the model coefficients for the abundance of alkenes and 
# alkadienes, and the mean chain length
coefficients_plot <- unsat_coef_plot /
  coeff_plot +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(3, 2, 2))
ggsave(here("figs", "chc-traits_coeff_plot.png")
       , plot = coefficients_plot
       , dpi = "print"
       , width = 6
       , height = 6.2
       , units = "in"
       , device = agg_png
       , scaling = 0.6)

# Export the data files ----
save(list = c("Prop_chain.length"
              , "Prop_CompsClass"
              , "cclasses_models_list"
              , "mean_cl_model")
     , file = here("data", "processed", "univariate_chc-stats.Rdata"))

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output", "SInf_univaraiteCHCstats-script.txt"))

## Detach/unload packages
pacman::p_unload(char = script_packs)

## Clear environment
objects_list <- ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|"))
rm(list = c(objects_list[nzchar(objects_list)], "objects_list"))

gc()

