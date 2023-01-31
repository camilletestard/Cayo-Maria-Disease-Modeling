
# X_Renaming ####

library(fs); library(tidyverse)

Files <- "Greg Data/Outputs/BISoN/Random" %>% dir_ls

FileDF <- 
  Files %>% 
  str_split("/") %>% 
  map_chr(last) %>% str_remove_all(".rds|PI_") %>% 
  str_split("_") %>% map(data.frame) %>% map(t) %>% map(as.data.frame) %>% 
  bind_rows %>% remove_rownames %>% 
  rename(PI = 1, Rep = 2)

FileDF %<>% 
  group_by(Rep) %>% 
  mutate(R = 1:n()) %>% 
  ungroup

FileDF %<>% 
  mutate(File = Files, 
         NewName = paste(Rep, R, PI, sep = "_")) %>% 
  mutate(NewFile = paste0("Greg Data/Outputs/BISoN/Random/", NewName))

file.rename(FileDF$File, FileDF$NewFile)
