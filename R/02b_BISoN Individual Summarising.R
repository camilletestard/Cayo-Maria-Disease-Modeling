
# 02b_BISoN Individual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)

theme_set(theme_cowplot())

FileList <- 
  "Data/Outputs/BISoN" %>% 
  dir_ls()

names(FileList) <- "Data/Outputs/BISoN" %>%
  list.files

IndivDFList <- 
  FileList %>% str_split("/") %>% map_chr(last) %>% str_split("_") %>% 
  map_chr(first) %>% unique %>% sort %>% 
  map(function(Pop){
    
    print(Pop)
    
    FileList[str_detect(FileList, Pop)] %>% 
      
      map(function(a){
        
        print(a)
        
        b <- a %>% readRDS
        
        # b %>% arrange(ID)
        
        b %>% return
        
      }) %>% bind_rows(.id = "File") %>% mutate(Pop = Pop)
    
  })

names(IndivDFList) <-
  FileList %>% str_split("/") %>% map_chr(last) %>% str_split("_") %>% 
  map_chr(first) %>% unique %>% sort

IndivDF <- IndivDFList %>% bind_rows()

IndivDF$File %<>% str_split("_") %>% map_chr(last) %>% str_remove(".rds")

IndivDF %<>% rename(S_I = File) %>% 
  mutate_at("S_I", as.numeric)

IndivDF %<>% 
  mutate(S_I_Category = cut(S_I, breaks = c(quantile(S_I, 0:3/3)), 
                            labels = c("Low", "Med", "High"), 
                            include.lowest = T))

IndivDF %<>% 
  group_by(ID, Pop, S_I_Category) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time),
            TimeSD = sd(Time),
            TimeSE = TimeSD/(1000^0.5))

# IndivDF %<>% 
#   group_by(ID, Pop) %>% 
#   summarise(MeanInf = mean(Infected), 
#             MeanTime = mean(Time),
#             TimeSD = sd(Time),
#             TimeSE = TimeSD/(1000^0.5))

# IDLevels <- IndivDF %>% arrange(MeanTime) %>% mutate_at("ID", ~paste0(.x, "_", Pop)) %>% pull(ID)

# IndivDF %>% arrange(MeanTime) %>% mutate_at("ID", ~paste0(.x, "_", Pop)) %>% 
#   mutate_at("ID", ~factor(.x, levels = IDLevels)) %>% 
#   ggplot(aes(ID, MeanTime)) +
#   geom_errorbar(aes(ymin = MeanTime - TimeSE, ymax = MeanTime + TimeSE, colour = Pop), 
#                 position = position_dodge(w = 0.4)) +
#   geom_point(aes(colour = Pop), 
#              position = position_dodge(w = 0.4)) +
#   theme(legend.position = "none")

IndivDF %>% saveRDS("Data/Outputs/IndividualTimesteps.rds")

