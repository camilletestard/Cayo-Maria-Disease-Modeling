
# X_Figures ####

{
  
  library(tidyverse); library(cowplot); library(fs); library(magrittr)
  library(patchwork); library(ggforce); library(colorspace); library(ggregplot)
  library(MCMCglmm)
  
  theme_set(theme_cowplot())
  
  dir_create("Figures")
  
}

# Figure 1: Epidemiological Simulations ####

## Panel A: Growth Curves ####

TimeStepDF <- readRDS("Data/Outputs/PopulationTimeSteps.rds")

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
    # theme(legend.position = "top", legend.justification = 0.5) +
    theme(legend.position = "top", 
          # legend.background = element_rect(colour = "grey"),
          legend.justification = 0.5) +
    facet_grid(Population ~ .) + 
    theme(strip.background = element_blank()) +#,
    guides(colour = guide_legend(reverse = F,
                                 direction = "horizontal",
                                 label.position = "top",
                                 label.theme = element_text(angle = 0), nrow = 1)) +
    NULL)

## Panel B: Camille Violins ####

individual_timestep <- readRDS("Data/Outputs/IndividualTimestepsCamille.rds")

individual_timestep %<>% mutate_at("year", as.factor)

(PanelB <- 
    ggplot(individual_timestep, aes(x = year, y = MeanTime, colour = year)) +
    geom_violin() +
    geom_sina(alpha = 0.1) +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022)) +
    theme(legend.position = "none") +
    labs(y = "Mean infection time", x = NULL) +
    scale_x_discrete(limits = as.factor(c(2013:2017, #"Maria", 
                                          2018:2019, 2021:2022))) +
    geom_vline(xintercept = 5.5, lty = 2, colour = "red"))

## Panel C: P_I versus mean infection timestep ####

OutputDF <- readRDS("Data/Outputs/PopulationTimes.rds")

OutputDF %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  mutate_at("Rep", ~str_replace(.x, "HH", "H")) %>% 
  mutate_at("Rep", ~str_replace(.x, "TT", "T")) %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5))

OutputDF %<>% mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

OutputDF %<>% mutate_at(c("Mean", "Max"), ~log(.x + 1))

OutputDF %<>% filter(between(Mean, 3, 7))

(P_I_Figure <- 
    OutputDF %>% 
    mutate_at("Year", ~factor(.x, levels = c(2013:2017, "\U2022 Maria \U2022", 2018, 2019, 2021, 2022))) %>% 
    ggplot(aes(P_I, exp(Mean), colour = Year, group = Rep)) + 
    geom_point(alpha = 0.2) +
    geom_smooth(method = lm, fill = NA) +
    scale_y_log10(limits = c(10, NA)) +
    labs(y = "Mean infection time", x = "Pathogen infectivity") +
    scale_colour_discrete_divergingx(palette = "PuOR", 
                                     limits = c(2013:2017, "[ Maria ]", 2018, 2019, 2021, 2022)) +    
    theme(legend.position = "none"))

## Panel D: Model effect comparisons ####

MCMC1 <- readRDS("Data/Outputs/EpidemiologyModel.rds")

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

## Combining ####

(TimeFigure)/
  PanelB/
  ((P_I_Figure + 
      theme(axis.title.x = element_text(vjust = 10)) + 
      EffectPlot) + 
     plot_layout(widths = c(4, 1))) +
  plot_layout(heights = c(1.5, 1, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("Figures/Figure1.jpeg", units = "mm", 
       width = 180, 
       height = 275, 
       dpi = 600)

ggsave("Figures/Figure1.pdf", 
       units = "mm", 
       width = 180, 
       height = 275, 
       dpi = 300)

# Figure 2 ####

# individual_timestep_Allrank <- readRDS("Data/Outputs/IndividualTimestepsCamille.rds")

IndivDF <- readRDS("Data/Intermediate/FullIndividualInfectionData.rds")

# IndivDF

IndivDF %<>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  # mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H"))) %>% 
  mutate_at("Rank", ~.x/100) %>%
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time)) %>% 
  ungroup

# individual_timestep_Allrank <- 
#   individual_timestep %>% 
#   filter(!is.na(`rank`))
# 
# individual_timestep_Allrank %<>% 
#   mutate_at("MeanTime", log10)


## Panel A: Model Effects ####

Model1 <- readRDS("Data/Intermediate/IndividualInfectionModelAdd.rds")

Model2 <- readRDS("Data/Intermediate/IndividualInfectionModelAddStrength.rds")

# IndivDF <- Model1$Data

(PanelA <- 
    list(Model2, 
         Model1) %>% 
    map("FinalModel") %>% 
    Efxplot(#VarNames = c("(((Intercept)))" = "Intercept"),
      Intercept = F,
      ModelNames = rev(c("Base", "+Strength")),
      VarOrder = rev(c("HurricanePost", "Age", "SexM", "Rank", 
                       paste0("HurricanePost_", c("Age", "SexM", "Rank")), "Strength"))) +
    scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]]),
                        limits = c("Base", "+Strength")) +
    theme(legend.position = c(0.1, 0.1), 
          legend.justification = c(0, 0)))

## Panel B: Rank ####

(PanelB <-
   IndivDF %>% 
   ggplot(aes(x = Rank, y = MeanTime)) +
   # geom_jitter(alpha = 0.3, aes(colour = isPost))+
   geom_point(alpha = 0.1, aes(colour = Hurricane))+
   geom_smooth(method = "lm", colour = "black")+
   # facet_grid(~Hurricane) +
   scale_x_continuous(breaks = c(0:5/5)) +
   # scale_y_continuous(breaks = c(0:5*250), 
   #                    limits = c(0, 1250),
   #                    labels = NULL,
   #                    name = NULL) +
   ggpubr::stat_cor() +
   scale_colour_manual(values = ParasiteColours[c(5, 4)]) +
   theme(strip.background = element_rect(fill = "white", #colour = "dark grey",
                                         colour = "white")) +
   labs(x = "Proportion dominant interactions") +
   theme(legend.position = "none") + scale_y_log10(breaks = 2^(1:10), 
                                                   name = "Mean infection timestep",
                                                   limits = c(30, 1250)))

## Panel C: Age ####

(PanelC <-
   IndivDF %>% 
   ggplot(aes(x = Age, y = MeanTime)) +
   # geom_jitter(alpha = 0.3, aes(colour = isPost))+
   geom_point(alpha = 0.1, aes(colour = Hurricane))+
   geom_smooth(method = "lm", colour = "black")+
   facet_grid(~Hurricane) +
   scale_x_continuous(breaks = c(0:6/2), labels = c(0:6*5)) +
   # scale_y_continuous(breaks = c(0:5*250), 
   #                    limits = c(0, 1250),
   #                    labels = NULL,
   #                    name = NULL) +
   ggpubr::stat_cor() +
   scale_colour_manual(values = ParasiteColours[c(5, 4)]) +
   theme(strip.background = element_rect(fill = "white", #colour = "dark grey",
                                         colour = "white")) +
   theme(legend.position = "none") + scale_y_log10(breaks = 2^(1:10), 
                                                   limits = c(30, 1250),
                                                   labels = NULL,
                                                   name = NULL))

## Panel D: Sex ####

SexSegments <- 
  data.frame(YFrom = c(1200, 1200)) %>% 
  mutate(YTo = YFrom, XFrom = c(1.1, 1.1)) %>% 
  mutate(XTo = c(1.9, 1.9), Label = rep("***")) %>% 
  mutate(X = (XFrom + XTo)/2) %>% 
  mutate(isPost = c("Pre", "Post")) %>% 
  mutate(Y = YFrom + 25) %>% 
  mutate_at("isPost", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  rename(Hurricane = isPost)

(PanelD <- 
    ggplot(IndivDF, aes(x = Sex, y = MeanTime))+
    geom_sina(alpha = 0.1, aes(colour = Hurricane)) +
    geom_boxplot(width = 0.2, outliers = F) +
    facet_grid(~Hurricane) +
    scale_colour_manual(values = ParasiteColours[c(5, 4)]) +
    scale_y_continuous(breaks = c(0:5*250), 
                       limits = c(0, 1250),
                       labels = NULL,
                       name = NULL) +
    geom_segment(data = SexSegments, #inherit.aes = F, 
                 aes(y = YFrom, yend = YTo, 
                     colour = Hurricane,
                     x = XFrom, xend = XTo)) +
    geom_text(data = SexSegments, inherit.aes = F, 
              aes(y = Y, x = X,
                  colour = Hurricane,
                  label = Label)) +
    theme(strip.background = element_rect(fill = "white", #colour = "dark grey",
                                          colour = "white")) +
    theme(legend.position = "none") + scale_y_log10(breaks = 2^(1:10), 
                                                    limits = c(30, 1250),
                                                    labels = NULL,
                                                    name = NULL))

## Combining ####

PanelA/
  (PanelB + theme(axis.title.y = element_text(vjust = -30))|PanelC|PanelD) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(1, 1.5))

ggsave("Figures/Figure2.jpeg", 
       units = "mm", 
       height = 200, 
       width = 300,
       dpi = 600)

ggsave("Figures/Figure2.pdf", 
       units = "mm", 
       height = 200, 
       width = 300,
       dpi = 600)

# ~~~~~ Supplement ####

# Repeatability Figures ####

ModelList <- readRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

ModelList %>% 
  map(MCMCRep) %>% 
  bind_rows(.id = "Model") %>% filter(Component != "units") %>% 
  # mutate_at(3:5, ~round(as.numeric(.x), 2))
  mutate_at(3:5, ~as.numeric(.x)) %>% 
  ggplot(aes(factor(Model), Mode, 
             fill = Component)) + 
  geom_col(colour = "black", position = "stack") + 
  labs(x = "Model") + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  lims(y = c(NA, 1)) +
  geom_hline(yintercept = 1, lty = 2, alpha = 0.4) +
  scale_x_discrete(labels = c("Overall", "ID:Hurricane", "ID Before", "ID After"))

ggsave("Figures/Repeatability.jpeg", units = "mm", height = 120, width = 150)

(FocalID <- IDDF$ID %>% unique %>% extract2(2))

(FocalID2 <- IDDF$ID %>% unique %>% extract(c(210, 219)))

OverlapDF <- 
  IDDF %>% group_by(ID) %>% summarise(MinYear = min(Year), MaxYear = max(Year)) %>% 
  filter(MinYear < 2020, MaxYear > 2020)

(FocalID2 <- c("8K1", "00V"))
(FocalID2 <- c("8K1", OverlapDF[10, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[14, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[15, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[16, ]$ID))

MaintainedLines <- 
  IDDF %>% mutate_at("Year", as.numeric) %>% 
  ggplot(aes(Year, Strength+0.001)) +
  geom_line(aes(colour = as.factor(!ID %in% FocalID2),
                alpha = as.factor(ID %in% FocalID2),
                group = ID)) +
  geom_point() +
  scale_y_log10() + 
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.01, 1), name = paste(FocalID2, collapse = ", ")) +
  scale_colour_manual(values = c("red", "black"), name = paste(FocalID2, collapse = ", "))

(BrokenFocalID2 <- c("8K1", OverlapDF[17, ]$ID))

BrokenLines <-
  IDDF %>% mutate_at("Year", as.numeric) %>% 
  ggplot(aes(Year, Strength+0.001)) +
  geom_line(aes(colour = as.factor(!ID %in% BrokenFocalID2),
                alpha = as.factor(ID %in% BrokenFocalID2),
                group = ID)) +
  geom_point() +
  scale_y_log10() + 
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.01, 1), name = paste(FocalID2, collapse = ", ")) +
  scale_colour_manual(values = c("red", "black"), name = paste(FocalID2, collapse = ", "))

library(patchwork)

MaintainedLines/BrokenLines

ggsave("Figures/Maintained_vs_Broken.jpeg", units = "mm", height = 180, width = 180)

# Repeatability tile plot ####

IndivDF %>% 
  # group_by(ID, Pop) %>% 
  group_by(ID) %>% 
  summarise(MeanTime = mean(MeanTime),
            N = n()) %>% 
  arrange(#Pop, 
    MeanTime) %>% 
  pull(ID) -> IDOrder

# IndivDF %>% 
#   group_by(ID) %>% 
#   summarise(MeanTime = mean(MeanTime),
#             N = n()) %>% 
#   pull(ID) %>% 
#   intersect(IDOrder) -> IDOrder

IndivDF %>% 
  mutate(Year = substr(Pop, 2, 5) %>% as.numeric) %>% 
  mutate_at("ID", ~factor(.x, levels = IDOrder)) %>% 
  mutate_at("MeanTime", log10) %>% 
  group_by(Year) %>% mutate_at("MeanTime", scale) %>%
  ggplot(aes(ID, Year)) +
  geom_tile(aes(fill = (MeanTime))) + 
  facet_grid(~Hurricane, scales = "free_x") + 
  coord_flip() +
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]])

IndivDF %>% 
  mutate(Year = substr(Pop, 2, 5) %>% as.numeric) %>% 
  mutate_at("ID", ~factor(.x, levels = IDOrder)) %>% 
  mutate_at("MeanTime", log10) %>% 
  group_by(Pop) %>% mutate_at("MeanTime", scale) %>%
  ggplot(aes(ID, Pop)) +
  geom_tile(aes(fill = (MeanTime))) + 
  facet_grid(~Hurricane, scales = "free_x") + 
  coord_flip() +
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]])

# Supplementary tables ####

Model1$FinalModel$summary.fixed %>% 
  as.data.frame %>% 
  dplyr::select(Estimate = mean, Lower = `0.025quant`, Upper = `0.975quant`) %>% 
  rownames_to_column(var = "Variable") %>% 
  write.csv("Figures/Model1Estimates.csv", row.names = F)

Model2$FinalModel$summary.fixed %>% 
  as.data.frame %>% 
  dplyr::select(Estimate = mean, Lower = `0.025quant`, Upper = `0.975quant`) %>% 
  rownames_to_column(var = "Variable") %>% 
  write.csv("Figures/Model2Estimates.csv", row.names = F)

MCMC1 <- readRDS("Data/Outputs/EpidemiologyModel.rds")

MCMCOutput <- 
  summary(MCMC1)$solutions %>% data.frame %>% 
  rownames_to_column("Var") %>% 
  rename(Estimate = 2, Lower = 3, Upper = 4, P = 6) %>% 
  dplyr::select(1:4, 6)

MCMCOutput[2, c("Estimate", "Lower", "Upper")] <-
  MCMCOutput[2, c("Estimate", "Lower", "Upper")]/10

MCMCOutput %>% 
  write.csv("Figures/PopulationModelEstimates.csv", row.names = F)

# Individual Model Comparison ####

ModelList <- readRDS("Data/Intermediate/IndividualInfectionModelList.rds")

names(ModelList) <- c("Gaussian", "Log-Gaussian", 
                      "Gamma",
                      "Poisson", "Negative binomial")

ModelList %>% 
  INLADICFig(Just = T) +
  scale_x_continuous(labels = c("Gaussian", "Log-Gaussian", 
                                "Gamma",
                                "Poisson", "Negative binomial"))

ggsave("Figures/FamilyComparison.jpeg")
