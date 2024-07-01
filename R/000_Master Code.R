
# 00_Master Code.R ####

library(tidyverse); library(fs); library(magrittr)

"R" %>% dir_ls(regex = ".R$") %>% extract(-1) %>% 
  map(source)
