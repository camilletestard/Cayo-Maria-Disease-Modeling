
# 02a_BISoN Annual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr)

theme_set(theme_cowplot())

FileList <- 
  "Greg Data/Outputs/BISoN/Random_Old" %>% 
  dir_ls()

names(FileList) <- "Greg Data/Outputs/BISoN/Random_Old" %>%
  list.files

OutputList <- 
  
  FileList %>% 
  
  map(function(a){
    
    print(which(FileList == a))
    
    b <- (a %>% readRDS)[[a %>% str_split("_") %>% map(3) %>% unlist %>% as.numeric]]
    
    # Maxes <- map(b, ~max(.x$Time)) %>% unlist
    # 
    # Means <- map(b, ~mean(.x$Time)) %>% unlist
    
    # TotalInf <- map(b, ~Prev(.x$Infected)) %>% unlist
    
    data.frame(Mean = mean(b$Time), 
               Max = max(b$Time), 
               File = a %>% str_split("/") %>% map_chr(last)) %>% return
    
  })

OutputDF <- OutputList %>% bind_rows(.id = "Rep")

OutputDF %<>% 
  mutate_at("File", ~str_remove(.x, ".rds")) %>% 
  separate(File, sep = "_", into = c("Rep", "R", "P_I")) %>% 
  mutate_at("P_I", as.numeric)

OutputDF %<>% 
  mutate_at("R", as.numeric)

# Testing ####

TestDF <- OutputDF

TestDF %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5))

TestDF %<>% mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

TestDF %<>% mutate_at(c("Mean", "Max"), ~log(.x + 1))

TestDF %<>% filter(Mean > 0)

TestDF %>% pull(Max) %>% qplot

TestDF %>% pull(Mean) %>% qplot

TestDF %>% SinaGraph("Rep", "Mean")

TestDF %>% SinaGraph("Rep", "Max")

TestDF %>% 
  ggplot(aes(P_I, Mean)) + 
  geom_point() +
  facet_wrap(~Rep) + 
  geom_smooth(method = lm)

TestDF %>% 
  ggplot(aes(P_I, exp(Mean), colour = Rep)) + 
  geom_point(alpha = 0.2) +
  # facet_wrap(~Rep) + 
  geom_smooth(fill = NA) +
  scale_y_log10(limits = c(1, NA)) +
  labs(y = "Mean Timestep")

LM1 <- lm(Mean ~ P_I + PostMaria, data = TestDF)

LM1 %>% summary

LMM1 <- lmer(Mean ~ P_I + PostMaria + (1|Rep), data = TestDF)

LMM1 %>% summary %>% extract2("coefficients") %>% data.frame()

LM1b <- lm(Mean ~ P_I + PostMaria + Year, data = TestDF)

LM1b %>% summary %>% extract2("coefficients") %>% data.frame()

# LM2 <- lm(Max ~ P_I + PostMaria + Rep, data = TestDF)
# 
# LM2 %>% summary
# 
# LMM2 <- lmer(Max ~ P_I + PostMaria + (1|Rep), data = TestDF)
# 
# LMM2 %>% summary


library(MCMCglmm)

library(INLA)

LM1 <- inla(Means ~ P_I + PostMaria + Population*Year, data = TestDF)

LM1 %>% summary

LM2 <- MCMCglmm(Maxes ~ P_I + PostMaria + Population*Year, data = TestDF)

LM2 %>% summary

# Converting to timesteps ####

TimestepList <- 
  
  FileList %>% 
  
  map(function(a){
    
    print(which(FileList == a))
    
    # b <- a %>% readRDS
    
    b <- (a %>% readRDS)[[a %>% str_split("_") %>% map(3) %>% unlist %>% as.numeric]]
    
    b
    
  })

TimeStepDF <- 
  TimestepList %>% 
  map(~.x %>% arrange(Time) %>% mutate(NInf = 1:n()) %>% mutate(PropInf = NInf/max(NInf))) %>% 
  bind_rows(.id = "File")

TimeStepDF %<>% 
  mutate_at("File", ~str_remove(.x, ".rds")) %>% 
  separate(File, sep = "_", into = c("Rep", "R", "P_I")) %>% 
  mutate_at(c("P_I", "R"), as.numeric) %>%
  mutate(File = paste(Rep, R, P_I, sep = "_"))

TimeStepDF %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5)) %>%   
  mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

TimeStepDF %>%
  mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
  ggplot(aes(Time + 1, PropInf)) + 
  geom_line(aes(group = File, colour = PostMaria), alpha = 0.05) + 
  scale_x_log10()

ggsave("TimeFigure.jpeg")

######

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

IndivListList <- FileList[1] %>% map(readRDS)

TestDF <- IndivListList %>% 
  map(~bind_rows(.x, .id = "Sim")) %>% 
  bind_rows(.id = "Rep")

TestDF %<>% 
  # separate(Rep, sep = "_", into = c("PI", "Rep"))
  mutate(PI = substr(Rep, 1, 3) %>% as.numeric) %>% 
  mutate(Rep = str_split(Rep, "_") %>% map_chr(first)) %>% 
  mutate(Rep = str_split(Rep, "_") %>% map_chr(last))

TestDF %<>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5) %>% as.numeric)

TestDF %<>% mutate(PostMaria = as.factor(as.numeric(Year >= 2018)))

LM1 <- lm(Time ~ as.factor(PI) + PostMaria + Rep, data = TestDF)

LM1 %>% summary
