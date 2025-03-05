###############################################################################
# 00_global_setup.R
#
# Purpose:
#   1) Clear environment and load/install relevant packages.
#   2) Set up a global ggplot2 theme (optional).
###############################################################################

# Clear the environment
rm(list = ls())

# Install/Load packages
if (!require("pacman")) install.packages("pacman")
library(pacman)
p_load_gh("DARTH-git/darthtools")

p_load(
  doParallel,
  devtools,
  diagram,
  readxl,
  ggplot2,
  dplyr,
  tidyr,
  dampack,      # cost-effectiveness analysis & parameterization
  reshape2,
  tidyplots,
  purrr,
  stringr
)

theme_set(theme_minimal())  # global theme for ggplot2

p_load_gh("DARTH-git/darthtools")  

###############################################################################
