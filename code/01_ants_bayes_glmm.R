### wog-wog-postfire-ants

### Matt Bitters
### matthew.bitters@colorado.edu





# Install and load required packages
packages <- c("here", "httr", "sf", "dplyr", "ggplot2", "terra", "purrr", "progressr")
installed <- packages %in% installed.packages()[, "Package"]
if (any(!installed)) {
  install.packages(packages[!installed])
}

library(here)
library(httr)
library(sf)
library(dplyr)
library(ggplot2)
library(terra)
library(purrr)
library(progressr)