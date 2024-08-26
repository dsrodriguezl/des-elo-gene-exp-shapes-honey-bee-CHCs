tic("Script to perform modeling over the wue_092021")

# Packages ----

# Install see if the package is not yet installed
install_new_packs("see", "DHARMa")

script_packs <- c("ggplot2"
                  , "ggtext"
                  , "performance"
                  , "tidymodels"
                  , "forcats"
                  , "gghalves"
                  , "emmeans"
                  , "ggeffects"
                  , "patchwork"
                  , "ragg"
)

## Install and/or load packages to be used by the script
load_my_packs(script_packs)

# ggplot pre-sets ----
source(here("scripts", "custom-ggplot-settings.R"))

# Prepare the data ----
## load data
load(here("data", "processed","CNQRs.Rdata"))

# Vector with the names of genes of interest (not reference genes)
genes <- qpcr_rel_data |> 
  select(-starts_with("Rpl"), -Individual) |> 
  colnames()

# Data frame with relative gene expression for each gene of interest in each
# sample
qpcr_rel_data <- qpcr_rel_data |> 
  # Exclude reference genes
  select(Individual, all_of(genes)) |> 
  # Add samples meta data
  merge(groups_qpcr |> 
          select(Individual, Subspecies:Task)
        , by = "Individual") |>
  # Remove NAs
  drop_na() |>
  # Ensure that levels of factors are limited to those included in the data frame
  mutate_if(is.factor
            , fct_drop) |>
  # Turn into a tibble
  as_tibble()

# Bootstrap ----
## Set a seed to ensure replicability
set.seed(12345)

## generate re-samples
bootstrapped_data <- qpcr_rel_data  |> 
  mutate(group = paste(Subspecies, Task, Hive, sep = "_")) |> 
  bootstraps(times = bootstrap_reps, strata = group, apparent = F)

## Empty list to hold the results of each bootstrapped model
models_list <- list()

## Vector listing objects in the environment that should not be erased during 
## or immediately after running the for loop
script_objects <- ls() |> 
  str_remove_all(analysis_objects |> 
                   rev() |> 
                   paste(collapse = "|")) |> 
  str_remove_all(custom_ggs |> 
                   paste(collapse = "|")) |> 
  str_subset("[:alpha:]")

## PDF to contain plots for evaluating the bootstrapped models
pdf(here("output", "model-diagnotsics_gene-expression.pdf")
    , width = 10, height = 8)

## FOR loop iterating through the genes of interest
for (gene in genes) {
  # Define model formula
  model_formula <- paste(gene
                         , "~"
                         , paste("Subspecies"
                                 , "Task"
                                 , sep = " + "))
  
  # Subset data frame to not include the other genes
  # It is just an attempt to make the final object ligther
  model_data <- qpcr_rel_data |>
    select(contains(gene), Individual, Subspecies:Task)
  
  # Fit the GLM without bootstrapping
  model <- glm(as.formula(model_formula)
               # The gene expression is continuous and 0 bounded,
               # without real 0s, thus the selection of the gamma family.
               # the log link function was chosen, because the data is right
               # skewed
               , family = Gamma(link = "log")
               , data = model_data) 

  # Print diagnostic plots of the non-bootstrapped GLM
  plot(c(0, 1)
       , c(0, 1)
       , ann = F
       , bty = 'n'
       , type = 'n'
       , xaxt = 'n'
       , yaxt = 'n')
  text(x = 0.5, y = 0.5
       , paste("glm - Gamma(link = 'log'):", model_formula)
       , cex = 1.6, col = "black")
  check_model(model, check = "all") |> print()
  
  # Get estimated marginal effects from the non-bootstrapped GLM
  pairwise_emmeans <- model |> 
    emmeans(~ Task * Subspecies, type = "response") 
  
  pairwise_emmeans <- list(emmeans = pairwise_emmeans
                           , contrasts = pairwise_emmeans |> 
                             contrast(method = "pairwise"
                                      , by = c("Subspecies")))|> 
    map(tidy, conf.int = T)
  
  pairwise_emmeans[["emmeans"]] <- pairwise_emmeans |> 
    pluck(1) |> 
    mutate(model = "glm")
  
  # Set the modeling workflow for the bootstrap re-samples
  glm_spec <- linear_reg() |>
    # The gene expression is continuous and 0 bounded,
    # without real 0s, thus the selection of the gamma family.
    # the log link function was chosen, because the data is right skewed
    set_engine("glm"
               , family = stats::Gamma(link = "log"))
  
  glm_wf <- workflow() |> 
    add_model(glm_spec) |> 
    add_formula(as.formula(model_formula))
  
  # Bootstrap the model
  bootstrapped_models <- {bootstrapped_data |> 
    mutate(
      # Fit the model
      model = map(splits
                  , ~fit(glm_wf
                         , data = analysis(.) |>
                           select(contains(gene)
                                  , Individual
                                  , Subspecies:Task)))
      # Extract coefficients
      , coef_info = map(model
                        , function(x)(
                          tidy(extract_fit_engine(x)
                               , conf.int = T) |>
                            # Transform coefficients for interpretation
                            # This is necessary due to the log link function
                            mutate(estimate = exp(estimate)
                                   , conf.low = exp(conf.low)
                                   , conf.high = exp(conf.high)))
                        )
      # Augment the model by predicting the OOB data
      , augmented = 
        map2(model
             , splits
             , function(x, y)(
               augment(extract_fit_engine(x)
                       , type.predict = "response"
                       , new_data = assessment(y) |>
                         select(contains(gene)
                                , Individual
                                , Subspecies:Task)) |> 
                 # Add proper reference to the factors
                 mutate(Task = ifelse(`TaskForager bees` == 0
                                      , "Nurse bees"
                                      , "Forager bees")
                        , Subspecies = ifelse(`SubspeciesA. m. iberiensis` == 0
                                              , "A. m. carnica"
                                              , "A. m. iberiensis"))
               )
        ))}
  
  # get confidence intervals for the bootstrapped coefficients
  coef_pct_intervals <- int_pctl(bootstrapped_models, coef_info)
  
  # get bootstrapped p values
  boot_p_values <- {bootstrapped_models |> 
    unnest(coef_info) |> 
    filter(term != "(Intercept)") |> 
    select(term, p.value)}
  
  # Plot bootstrapped p values distributions
  p_value_dist <- {boot_p_values |> 
    ggplot(aes(x = p.value)) +
    geom_density(fill = "blue4") +
    facet_wrap(term |>
                 str_replace_all("A. m. iberiensis"
                                 , " (A. m. iberiensis)") |> 
                 str_replace_all("Forager bees"
                                 , " (Forager bees)") |> 
                 str_remove_all("`") |> 
                 vars()
               , ncol = 1
               , scales = "free") +
    theme_classic()}
  print(p_value_dist)
  
  # Get bootstrapped p values confidence intervals
  p_pct_intervals <- {boot_p_values |> 
    group_by(term) |> 
    summarise(boot_p_value = median(p.value)
              , p_value_conf_low = quantile(p.value, 0.025)
              , p_value_conf_high = quantile(p.value, 0.975))}
  
  # Extract coefficients
  boot_coefficients <- {bootstrapped_models |>
    unnest(coef_info) |> 
    select(id, term:conf.high) |> 
    merge(coef_pct_intervals) |>
    filter(term != "(Intercept)") |>  
    merge(p_pct_intervals) |> 
    as_tibble()}
  
  # Plot bootstrapped coefficients
  boot_coef_plot <- { boot_coefficients |>
      mutate(term = term |>
               str_replace_all("A. m. iberiensis"
                               , " (A. m. iberiensis)") |>
               str_replace_all("Forager bees"
                               , " (Forager bees)") |>
               str_remove_all("`")) |>
    ggplot(aes(x = term
               , y = estimate
               , fill = term)) +
    geom_half_violin(color = "grey60"
                       , alpha = 0.95
                     , linewidth = 0.2
                     , side = "r"
                     , scale = "width"
                     , nudge = 0.06) +
    geom_linerange(aes(ymin = .lower
                       , ymax = .upper)
                   , color = "grey35"
                   , linewidth = 1) +
    geom_point(aes(y = .estimate)
               , color = "grey35"
               , size = 3.5) +
    geom_text(aes(x = term
                  , y = .estimate
                  , label =
                    paste(round(.estimate, 3)
                          , ifelse(boot_p_value >= 0.05
                                   , "(n.s)"
                                   , ifelse(boot_p_value >= 0.01
                                            , "(*)"
                                            , ifelse(boot_p_value >= 0.001
                                                     , "(**)"
                                                     , "(***)")))))
              , size = 4
              , nudge_x = -0.05
              , color = "gray20") +
    boxplot_theme +
    coord_flip() +
    geom_hline(yintercept = 1
               , linewidth = 1
               , linetype = "dashed"
               , color = "grey20") +
    labs(y = "Coefficients", x = NULL) +
    scale_fill_viridis_d(option = "turbo") +
    guides(fill = "none") +
    theme(axis.title = element_text(size = 9)
          , axis.text = element_text(size = 8))}
  print(boot_coef_plot)
  
  # Get bootstrapped predicted values
    boot_predictions <- bootstrapped_models |> 
      unnest(augmented) |> 
      group_by(Task, Subspecies) |> 
      summarise(predicted = median(.fitted)
                , conf.low = quantile(.fitted, 0.025)
                , conf.high = quantile(.fitted, 0.975)) |> 
      mutate(model = "bootstrapped")

    # plot comparing predictions of bootstrapped and non-bootstrapped GLM
    boot_vs_model <- {model_data |> 
        merge(pairwise_emmeans |> 
                pluck(1) |> 
                mutate(predicted = response
                       , .keep = "unused") |> 
                select(Task
                       , Subspecies
                       , predicted
                       , conf.low
                       , conf.high
                       , p.value
                       , model)
              , all = T
              , sort = F) |> 
        merge(boot_predictions
              , all = T
              , sort = F) |> 
        as_tibble() |> 
        ggplot(aes(x = Task, y = predicted)) +
        geom_pointrange(aes(ymin = conf.low
                            , ymax = conf.high
                            , color = model
                            , fill = model)
                        , size = 1
                        , linewidth = 1.5
                        , position = position_dodge(width = 0.3)) +
        geom_jitter(aes(y = get(gene))
                    , alpha = 0.7
                    , size = 2
                    , position = position_jitter(width = 0.2
                                                 , height = 0
                                                 , seed = 12345)
                    , color = "grey50") +
        scale_color_viridis_d(option = "turbo") +
        scale_fill_viridis_d(option = "turbo") +
        labs(y = paste("Relative expresion of", gene)) +
        coord_flip() +
        facet_wrap(~Subspecies
                   , nrow = 2) +
        boxplot_theme +
        theme(axis.title = element_text(size = 9)
              , axis.text = element_text(size = 8))}
    print(boot_vs_model)
  
  # Evaluate the model across the bootstrap resamplings
  boot_metrics <- glm_wf |> 
    fit_resamples(bootstrapped_models)
  
  collect_notes(boot_metrics) |> 
    print()
  
  boot_metrics <- {collect_metrics(boot_metrics
                                  , summarize = F) |>
      group_by(.metric) |> 
      summarise(median = median(.estimate)
                , mean = mean(.estimate)
                , sd = sd(.estimate)
                , conf.low = quantile(.estimate, 0.025)
                , conf.high = quantile(.estimate, 0.975)
                , n = n()
                , std.error = sd/sqrt(n))}
  
  # Store results in a list indexed within models_list
    models_list[[gene]] <- {
      list(data_modeling  = list(model_formula = model_formula
                                 , model = model
                                 , model_coef = tidy(model
                                                     , conf.int = T) |>
                                   mutate(estimate = exp(estimate)
                                          , conf.low = exp(conf.low)
                                          , conf.high = exp(conf.high))
                                 , model_metrics = glance(model)
                                 , emmeans = pairwise_emmeans)
           , bootstrap_modeling = 
             list(coefficients = boot_coefficients
                  , coef_percentile_intervals = coef_pct_intervals
                  , metrics = boot_metrics
                  , p_values = boot_p_values
                  , p_val_percentiles = p_pct_intervals
                  , boot_predictions = boot_predictions))}
  
}
dev.off()

# List of objects created by the for loop that can be erased from the environment
loop_checheres <- {ls() |> 
  str_remove_all(analysis_objects |>
                   rev() |> 
                   paste(collapse = "|")) |> 
  str_remove_all(custom_ggs |> 
                   paste(collapse = "|")) |>  
  str_remove_all(script_objects |> 
                   paste(collapse = "|")) |> 
  str_remove("script_objects")}
  
# decluter memory by removing objects in loop_checheres from the environment
rm(list = c(loop_checheres[nzchar(loop_checheres)], "loop_checheres"))

# Export results
destination_path <- here("data", "processed", "rel-exp_boot-models.Rdata")
save(models_list, file = destination_path)

# End ----
## Report session information
capture.output(sessionInfo()
               , file = here("output"
                             , "SInf__rel-exp_boot-models.txt"))

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
