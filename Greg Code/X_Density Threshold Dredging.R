
# Setting the threshold at a given density ####

rm(list = ls())

library(dplyr); library(igraph); library(foreach); library(doParallel); library(stringr); library(tidyverse)
library(magrittr); library(fs)

# load("Data/R.Data/BisonFittedNetworks.RData")
load("Data/BisonFittedNetworks (2).RData")

DesiredDensity <- 0.2

FindThreshold <- function(Connections){
  
  Ranks <- rank(Connections)/max(rank(Connections)) # Identifies the numerical order of the values
  
  Connections[Ranks < (1-DesiredDensity)] <- 0 # Assigns the bottom 80% as 0 (according to the density desired)
  
  Connections %>% return # Returns
  
}

posterior.el[[1]][,3] %>% 
  FindThreshold # Does this once

posterior.el[[1]][,3:1002] %>% 
  map(FindThreshold) # Does this across all draws in the data frame

posterior.el[[1]][,3:1002] %>% 
  map(FindThreshold) %>% 
  bind_cols() %>% 
  bind_cols(posterior.el[[1]][,1:2], .)

