
# 00_Master Code.R ####

library(tidyverse); library(fs); library(magrittr)

"R" %>% dir_ls(regex = ".R$") %>%
  extract(-c(1:2)) %>% 
  map(source)
