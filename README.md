# Expression of Elongase- and Desaturase-Encoding Genes Shapes the Cuticular Hydrocarbon Profiles of Honey Bees

This repository contains all the data and code needed to reproduce the analyses of Rodríguez-León, D.S., Schmitt, T., Pinto, M.A., Thamm, M. and Scheiner, R. (2025), Expression of Elongase- and Desaturase-Encoding Genes Shapes the Cuticular Hydrocarbon Profiles of Honey Bees. Mol Ecol e17716. [https://doi.org/10.1111/mec.17716](https://doi.org/10.1111/mec.17716).

The repository has been published under a CC BY 4.0 license in [Figshare](https://figshare.com) ([https://doi.org/10.6084/m9.figshare.26831491.v1](https://doi.org/10.6084/m9.figshare.26831491.v1)).
Cite it as: Rodríguez-León, DS; Schmitt, Thomas; Pinto, Alice; Thamm, Markus; Scheiner, Ricarda (2025). Expression of elongase- and desaturase-encoding genes shapes the cuticular hydrocarbon profiles of honey bees. figshare. Software. https://doi.org/10.6084/m9.figshare.26831491.v1

The repository corresponds to an R project repository that can be
directly open in RStudio via the
des&elo-gene-exp-shapes-honey-bee-CHCs.Rproj file. This makes the R
project portable, as the working directory is automatically set within
RStudio to the folder of the repository.

This R project uses the package renv to create a reproducible
environment. To download and install all project’s dependencies, use
`renv::restore()` the first time you open the R project in RStudio. This
will ensure you can run the code properly to reproduce the analyses of
the study. See <https://rstudio.github.io/renv/> for more details on the
usage of renv package.

## Data:

All data related to the study is stored in the “data” folder. The folder
is divided in three following subfolders:

- **raw:** Contains data files (i.e. CSV files) with the data as
  originally obtained.
  - raw/samples_list.csv contains the samples metadata.
  - raw/gcms-integration folder contains CSV files with GC-MS
    integration results of every sample. It is divided into two
    sub-folders, corresponding to the two batches in which the samples
    were run in the GC-MS machine(s).
  - raw/qPCR folder contains CSV files with measured Ct values for each
    qPCR run. Files in sub-folder “samples_runs” correspond to the
    measured Ct values for the analyzed samples, the six CSV files
    outside of this sub-folder correspond to the standard curves.
- **tmp:** Temporary data files produced during the processing of the
  raw data, in order to prepare the data set(s) for statistical
  analysis.
- **processed:** Final data files (i.e. Rdata files) with the final data
  set(s), ready for statistical analyses.

## Code:

Analysis.R is the master script. It contains the code to replicate the
entire set of analyses of the study. It also generates the figures for
the manuscript, storing them in the folder “figs”. The code in the
script relies on the des&elo-gene-exp-shapes-honey-bee-CHCs.Rproj to
determine the working directory within RStudio.

The folder “scripts” contain R scripts with code that is needed for the
Analysis.R script to perform the analyses of the study.
