
# 02a_BISoN Annual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)

theme_set(theme_cowplot())

FileList <- 
  "Greg Data/Outputs/BISoN" %>% 
  dir_ls()

names(FileList) <- "Greg Data/Outputs/BISoN" %>%
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

# TestDF %>% 
#   ggplot(aes(P_I, Mean)) + 
#   geom_point() +
#   facet_wrap(~Rep) + 
#   geom_smooth()
# 
# TestDF %>% 
#   ggplot(aes(P_I, exp(Mean), colour = Rep)) + 
#   geom_point(alpha = 0.2) +
#   # facet_wrap(~Rep) + 
#   geom_smooth(method = lm, fill = NA) +
#   scale_y_log10(limits = c(10, NA)) +
#   labs(y = "Mean Timestep")

(P_I_Figure <- 
    TestDF %>% 
    mutate_at("Year", ~factor(.x, levels = c(2015:2017, "[ Maria ]", 2018, 2019, 2021))) %>% 
    ggplot(aes(P_I, exp(Mean), colour = Year, group = Rep)) + 
    geom_point(alpha = 0.2) +
    # facet_wrap(~Rep) + 
    geom_smooth(method = lm, fill = NA) +
    scale_y_log10(limits = c(10, NA)) +
    labs(y = "Mean infection timestep", x = "Pathogen infectivity") +
    # scale_colour_discrete_sequential("Blues 2") +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2015:2017, "[ Maria ]", 2018, 2019, 2021)) +    
    theme(legend.position = "none"))

LM1 <- lm(Mean ~ P_I + PostMaria, data = TestDF)

LM1 %>% summary

LMM1 <- lmer(Mean ~ P_I + PostMaria + (1|Rep), data = TestDF)

Coef <- LMM1 %>% summary %>% extract2("coefficients") %>% data.frame()
ConfInt <- LMM1 %>% confint()

Intercept <- 200

Errors <- 
  data.frame(P_I = c(0.1, 0.2), 
             Mean = c(Intercept, exp(log(Intercept) + Coef[3, 1])),
             Lower = c(NA, exp(log(Intercept) + ConfInt[5, 1])),
             Upper = c(NA, exp(log(Intercept) + ConfInt[5, 2])))

(P_I_Figure <- 
  P_I_Figure +
  geom_line(data = Errors, aes(P_I, Mean), 
            colour = "white",
            inherit.aes = F, size = 2, lty = 1) +
  geom_errorbar(data = Errors, aes(x = P_I, ymin = Lower, ymax = Upper), 
                colour = "white",
                width = 0.01, size = 2,
                inherit.aes = F) +
  geom_point(data = Errors, aes(P_I, Mean), inherit.aes = F, size = 2.5, colour = "white") +
  geom_point(data = Errors, aes(P_I, Mean), inherit.aes = F))

ggsave("Figures/P_I_Mean.jpeg", units = "mm", width = 150, height = 150, dpi = 300)

# Converting to timesteps ####

TimestepLabelDF <- data.frame(Population = c("F", "K", "S", "V"), 
                              Time = 2, 
                              PropInf = 0.9)

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
  mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
  
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5)) %>%   
  
  mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

(TimeFigure <- 
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
    # theme(strip.background = element_rect(fill = "white", colour = NA)) +
    theme(strip.background = element_blank(),
          strip.text = element_blank()) +
    geom_text(data = TimestepLabelDF, 
              inherit.aes = F,
              aes(x = Time + 1, y = PropInf, 
                  label = paste0("Population: ", Population))) +
    guides(colour = guide_legend(reverse = F,
                                 direction = "horizontal",
                                 # title.position = "left",
                                 # title.vjust = 0.25, title.hjust = 1,
                                 label.position = "top",
                                 # label.hjust = 0.5,
                                 # label.vjust = 1.5,
                                 label.theme = element_text(angle = 0), nrow = 1)) +
    NULL)

ggsave("Figures/TimeFigure.jpeg", units = "mm", 
       width = 150, height = 150, dpi = 300)

TimeFigure + P_I_Figure

ggsave("Figures/Figure2.jpeg", units = "mm", 
       width = 250, height = 150, dpi = 300)
