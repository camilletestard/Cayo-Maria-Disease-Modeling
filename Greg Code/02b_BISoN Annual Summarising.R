
# 02a_BISoN Annual Summarising ####

library(tidyverse)

IndivListList <- 
  "Greg Data/Outputs/BISoN" %>% 
  dir_ls() %>% 
  map(readRDS)

names(IndivListList) <- "Greg Data/Outputs/BISoN" %>% 
  list.files %>% str_remove(".rds$") %>% str_remove("PI_")

Maxes <- 
  
  IndivListList %>% 
  
  map(function(a){
    
    map(a, ~max(.x$Time)) %>% unlist
    
  }) %>% bind_rows()

(LongMaxes <- Maxes %>% reshape2::melt())

LongMaxes %<>% 
  separate(variable, sep = "_", into = c("PI", "Rep"))

LongMaxes %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5) %>% as.numeric)

LongMaxes %>% 
  ggplot(aes(value, colour = Rep)) + 
  geom_density() +
  facet_wrap(~Year)


LongMaxes %>% 
  ggplot(aes(value, colour = Rep)) + 
  geom_density() +
  facet_wrap(~Year)




