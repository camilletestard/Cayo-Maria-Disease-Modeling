
# 02a_BISoN Annual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)

theme_set(theme_cowplot())

FileList <- 
  "Data/Outputs/BISoN" %>% 
  dir_ls()

names(FileList) <- "Data/Outputs/BISoN" %>%
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

OutputDF %>% saveRDS("Data/Output/PopulationTimes.rds")

# Testing ####

TestDF <- OutputDF

TestDF %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate_at("Rep", ~str_replace(.x, "HH", "H")) %>% 
  mutate_at("Rep", ~str_replace(.x, "TT", "T")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5))

TestDF %<>% mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

TestDF %<>% mutate_at(c("Mean", "Max"), ~log(.x + 1))

TestDF %<>% filter(between(Mean, 3, 7))

(P_I_Figure <- 
    TestDF %>% 
    mutate_at("Year", ~factor(.x, levels = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022))) %>% 
    ggplot(aes(P_I, exp(Mean), colour = Year, group = Rep)) + 
    geom_point(alpha = 0.2) +
    geom_smooth(method = lm, fill = NA) +
    scale_y_log10(limits = c(10, NA)) +
    labs(y = "Mean infection timestep", x = "Pathogen infectivity") +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2013:2017, "[ Maria ]", 2018, 2019, 2021, 2022)) +    
    theme(legend.position = "none"))

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
  mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
  
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate_at("Rep", ~str_replace(.x, "HH", "H")) %>% 
  mutate_at("Rep", ~str_replace(.x, "TT", "T")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5)) %>%   
  
  mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

OutputDF %>% saveRDS("Data/Output/PopulationTimeSteps.rds")

(TimeFigure <- 
    TimeStepDF %>%
    # mutate(Population = paste0("Group: ", Population)) %>% 
    # mutate(Population = paste0("Group: ", Population)) %>% 
    mutate(File = paste(Rep, R, P_I, sep = "_")) %>% 
    mutate_at("Year", ~factor(.x, levels = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022))) %>% 
    ggplot(aes(Time + 1, PropInf, colour = Year)) + 
    geom_line(aes(group = File), alpha = 0.05) + 
    geom_line(data = . %>% filter(Time > 10000), aes(group = File), alpha = 1) + 
    scale_y_continuous(breaks = c(0:2/2)) +
    scale_x_continuous(limits = c(0, 1500)) +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022)) +
    labs(colour = NULL, x = "Time step", y = "Proportion infected") +
    theme(legend.position = "top", legend.justification = 0.5) +
    facet_grid(Population ~ .) + 
    theme(strip.background = element_blank()) +#,
    guides(colour = guide_legend(reverse = F,
                                 direction = "horizontal",
                                 label.position = "top",
                                 label.theme = element_text(angle = 0), nrow = 1)) +
    NULL)

# Running models ####

library(MCMCglmm)

MCMC1 <- MCMCglmm(Mean ~ P_I + PostMaria, random =~Rep, data = TestDF)

MCMC1 %>% saveRDS("Model.rds")

MCMC1 <- readRDS("Model.rds")

# MCMC1b <- MCMCglmm(Mean ~ P_I + PostMaria, data = TestDF)

MCMCOutput <- 
  MCMC1$Sol %>% data.frame %>% 
  mutate_at("P_I", ~.x/10) %>% 
  dplyr::select(P_I, PostMaria1) %>% 
  pivot_longer(1:2, names_to = "rowname", values_to = "Estimate")

CoefDF <- MCMC1 %>% summary %>% extract2(5) %>% as.data.frame %>% 
  rename(Lower = 2, Upper = 3, Estimate = 1) %>% 
  slice(2:3)

CoefDF[1, c("Estimate", "Lower", "Upper")] <-
  CoefDF[1, c("Estimate", "Lower", "Upper")]/10

(EffectPlot <- 
    CoefDF %>% 
    rownames_to_column() %>% 
    ggplot(aes(rowname, Estimate)) +
    geom_violin(data = MCMCOutput, scale = "width", fill = AlberColours[[1]], alpha = 0.3, colour = NA) +
    geom_hline(lty = 2, alpha = 0.3, yintercept = 0) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.5) +
    geom_point(size = 2) +
    geom_point(colour = "white") +
    scale_x_discrete(labels = c("+10% Infectivity", "Hurricane")) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
    labs(y = "Effect estimate", x = NULL))

(TimeFigure)/
  ((P_I_Figure + 
      theme(axis.title.x = element_text(vjust = 10)) + 
      EffectPlot) + 
     plot_layout(widths = c(4, 1))) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A")

dir_create("Figures")

ggsave("Figures/Figure1.jpeg", units = "mm", 
       width = 250, height = 250, dpi = 600)

(TimeFigure)/
  ((P_I_Figure + 
      theme(axis.title.x = element_text(vjust = 10)) + 
      EffectPlot) + 
     plot_layout(widths = c(2, 1))) +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1Thin.jpeg", units = "mm", 
       width = 180, height = 250, dpi = 600)

ggsave("Figures/Figure1Thin.pdf", units = "mm", 
       width = 180, height = 250, dpi = 600)
