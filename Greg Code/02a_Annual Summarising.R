
# 02a_Annual Summarising ####

library(tidyverse)

IndivListList <- 
  "Greg Data/Outputs" %>% 
  dir_ls() %>% 
  map(readRDS)

Maxes <- 
  
  IndivListList %>% 
  
  map(function(a){
    
    map(a, ~max(.x$Time))
    
  })