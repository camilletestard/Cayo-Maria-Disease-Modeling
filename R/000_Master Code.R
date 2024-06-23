
# 00_Master Code.R ####

library(tidyverse); library(fs)

"R" %>% dir_ls(regex = ".R$") %>% sort
