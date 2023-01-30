
# 02a_BISoN Annual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr)

theme_set(theme_cowplot())

FileList <- 
  "Greg Data/Outputs/BISoN/Random" %>% 
  dir_ls()

names(FileList) <- "Greg Data/Outputs/BISoN/Random" %>%
  list.files

OutputList <- 
  
  FileList %>% 
  
  map(function(a){
    
    print(which(FileList == a))
    
    b <- a %>% readRDS
    
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

TestDF %>% 
  ggplot(aes(P_I, Mean)) + 
  geom_point() +
  facet_wrap(~Rep) + 
  geom_smooth()

TestDF %>% 
  ggplot(aes(P_I, exp(Mean), colour = Rep)) + 
  geom_point(alpha = 0.2) +
  # facet_wrap(~Rep) + 
  geom_smooth(method = lm, fill = NA) +
  scale_y_log10(limits = c(10, NA)) +
  labs(y = "Mean Timestep")

LM1 <- lm(Mean ~ P_I + PostMaria, data = TestDF)

LM1 %>% summary

library(lme4); library(lmerTest)

LMM1 <- lmer(Mean ~ P_I + PostMaria + (1|Rep), data = TestDF)

LMM1 %>% summary %>% extract2("coefficients") %>% data.frame()

# Converting to timesteps ####

TimestepList <- 
  
  FileList %>% map(readRDS)

TimeStepDF <- 
  TimestepList %>% 
  map(~.x %>% arrange(Time) %>% filter(Time >= 0) %>% mutate(NInf = 1:n()) %>% mutate(PropInf = NInf/nrow(.x))) %>% 
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

# MeanLevels <- 
#   TimeStepDF %>% 
#   group_by(Rep, Time, PostMaria) %>% 
#   summarise_at("PropInf", mean)

TimeStepDF %>%
  mutate_at("PostMaria", ~factor(.x, levels = c("1", "0"))) %>% 
  mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
  ggplot(aes(Time + 1, PropInf, colour = PostMaria)) + 
  geom_line(aes(group = File), alpha = 0.05) + 
  geom_line(data = . %>% filter(Time > 10000), aes(group = File), alpha = 1) + 
  scale_x_log10() + scale_y_continuous(breaks = c(0:2/2)) +
  scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]]), 
                      labels = rev(c("Pre-Maria", "Post-Maria"))) +
  labs(colour = NULL, x = "Time step", y = "Proportion infected") +
  theme(legend.position = "top") +
  
  # geom_smooth(aes(group = Rep), alpha = 0.5, fill = NA) + 
  facet_grid(Population ~ .) + 
  theme(strip.background = element_rect(fill = "white", colour = NA)) +
  NULL

ggsave("TimeFigure.jpeg", units = "mm", 
       width = 150, height = 150, dpi = 300)

library(colorspace)

TimeStepDF %>%
  mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
  mutate_at("Year", ~factor(.x, levels = c(2015:2017, "[ Maria ]", 2018, 2019, 2021))) %>% 
  ggplot(aes(Time + 1, PropInf, colour = Year)) + 
  geom_line(aes(group = File), alpha = 0.05) + 
  geom_line(data = . %>% filter(Time > 10000), aes(group = File), alpha = 1) + 
  scale_x_log10() + scale_y_continuous(breaks = c(0:2/2)) +
  # scale_colour_discrete_diverging(palette = "Red-Green") +
  scale_colour_discrete_divergingx(palette = "PuOR", 
                                   limits = c(2015:2017, "[ Maria ]", 2018, 2019, 2021)) +
  labs(colour = NULL, x = "Time step", y = "Proportion infected") +
  theme(legend.position = "top", legend.justification = 0.5) +
  
  # geom_smooth(aes(group = Rep), alpha = 0.5, fill = NA) + 
  facet_grid(Population ~ .) + 
  theme(strip.background = element_rect(fill = "white", colour = NA)) +
  guides(colour = guide_legend(reverse = F,
                               direction = "horizontal",
                               # title.position = "left",
                               # title.vjust = 0.25, title.hjust = 1,
                               label.position = "top",
                               # label.hjust = 0.5,
                               # label.vjust = 1.5,
                               label.theme = element_text(angle = 0), nrow = 1)) +
  NULL

ggsave("TimeFigure.jpeg", units = "mm", 
       width = 150, height = 150, dpi = 300)

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
