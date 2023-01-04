
# 02a_BISoN Annual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr)

theme_set(theme_cowplot())

IndivListList <- 
  "Greg Data/Outputs/BISoN/Random" %>% 
  dir_ls() %>% extract2(1) %>% 
  map(readRDS)

names(IndivListList) <- "Greg Data/Outputs/BISoN" %>% 
  list.files %>% 
  setdiff("Random") %>% 
  str_remove(".rds$") %>% str_remove("PI_")

Maxes <- 
  
  IndivListList %>% 
  
  map(function(a){
    
    map(a, ~max(.x$Time)) %>% unlist
    
  }) %>% bind_rows()

Means <- 
  
  IndivListList %>% 
  
  map(function(a){
    
    map(a, ~mean(.x$Time)) %>% unlist
    
  }) %>% bind_rows()

(LongMaxes <- Maxes %>% reshape2::melt() %>% rename(Max = 2) %>% 
    bind_cols(Means %>% reshape2::melt() %>% rename(Mean = 2) %>% dplyr::select(-1)))

LongMaxes %<>% 
  separate(variable, sep = "_", into = c("PI", "Rep"))

LongMaxes %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5) %>% as.numeric)

LongMaxes %<>% mutate(PostMaria = as.factor(as.numeric(Year >= 2018)))

library(patchwork)

LongMaxes %>% 
  ggplot(aes(as.factor(Year), Max, colour = PostMaria)) +
  geom_boxplot() +
  geom_sina(alpha = 0.1) +
  facet_grid(PI~Population) + 
  # ggtitle("Max") +
  labs(x = "Year", y = "Maximum time step infected") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +

LongMaxes %>% 
  ggplot(aes(as.factor(Year), Mean, colour = PostMaria)) +
  geom_boxplot() +
  geom_sina(alpha = 0.1) +
  facet_grid(PI~Population) + 
  # ggtitle("Mean") +
  labs(x = "Year", y = "Mean time step infected") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
  
  plot_layout(guides = "collect")

ggsave("Figures/Pre_Post_BISoN.jpeg", units = "mm", width = 350, height = 150)


# Running an LM ####

TestDF <- IndivListList %>% 
  map(~bind_rows(.x, .id = "Sim")) %>% 
  bind_rows(.id = "Rep")

TestDF %<>% 
  # separate(Rep, sep = "_", into = c("PI", "Rep"))
  mutate(PI = substr(Rep, 1, 3) %>% as.numeric) %>% 
  mutate(Rep = str_split(Rep, "_") %>% map_chr(last))

TestDF %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5) %>% as.numeric)

TestDF %<>% mutate(PostMaria = as.factor(as.numeric(Year >= 2018)))

LM1 <- lm(Time ~ as.factor(PI) + PostMaria + Rep, data = TestDF)

LM1 %>% summary
