# Function to install all packages needed for the script that are not yet installed ----
install_new_packs <- function(packages) {
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
}

# Load the packages needed for the script ----
load_my_packs <- function(packages) {
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  
  lapply(packages, library, character.only = TRUE) |>
    invisible()
}
