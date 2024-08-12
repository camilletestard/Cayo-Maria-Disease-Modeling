
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
    theme(legend.position = "top", legend.justification = 0.5) +
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
    list(Model1, Model2) %>% 
    map("FinalModel") %>% 
    Efxplot(#VarNames = c("(((Intercept)))" = "Intercept"),
      Intercept = F,
      ModelNames = c("Base", "+Strength")) +
    scale_colour_manual(values = c(AlberColours[[1]], AlberColours[[2]])) +
    theme(legend.position = c(0.1, 0.1), legend.justification = c(0, 0)))

## Group differences: Model 1 ####

if(!file.exists("Data/Outputs/ComparisonDF1.csv")){
  
  IndivDF %<>% as.data.frame
  
  PredDF <- MakePredictDF(IndivDF %>% dplyr::select(-c(MeanInf, MeanTime)), 
                          
                          HoldNumeric = c("Age")
                          # HoldFactor = c(#"Sex" = "F", "Pop" = "K2020",
                          #                "ID" = "00o")
  ) %>% 
    expand.grid
  
  PredDF %<>% filter(ID %in% unique(ID)[1:2], 
                     Pop %in% unique(Pop)[1:2])
  
  PredDF[,paste0("Pred.", 1:1000)] <- 
    INLAFit(Model1$FinalModel, PredDF, 
            FixedCovar = c("Hurricane", "Rank", "Age", "Sex", 
                           paste0("Hurricane:", c("Rank", "Age", "Sex"))), 
            # Return = "Matrix",
            Draw = T, NDraw = 1000) %>% 
    map(as.data.frame) %>% 
    bind_rows() %>% 
    t
  
  PredDF %<>% 
    filter(ID == ID[1], Pop == Pop[1])
  
  LongPred <- 
    PredDF %>% 
    # filter(ID == ID[1], Pop == Pop[1]) %>% 
    pivot_longer(matches("Pred"))
  
  i <- 1
  
  ComparisonDF1 <- 
    1:1000 %>%
    map(function(i){
      
      data.frame(
        
        PreLM = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("L", "M"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreMH = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("M", "H"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreLH = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("L", "H"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostLM = PredDF %>% filter(Hurricane == "Post", Rank %in% c("L", "M"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostMH = PredDF %>% filter(Hurricane == "Post", Rank %in% c("M", "H"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostLH = PredDF %>% filter(Hurricane == "Post", Rank %in% c("L", "H"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreFM = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("M")) %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostFM = PredDF %>% filter(Hurricane == "Post", Rank %in% c("M")) %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff
        
      )
      
    }) %>% 
    bind_rows(.id = "Rep")
  
  ComparisonDF1 %>% write.csv("Data/Outputs/ComparisonDF1.csv", row.names = F)
  
}else ComparisonDF1 <- read.csv("Data/Outputs/ComparisonDF2.csv")

(Model1ComparisonPanel <- 
    ComparisonDF1 %>% 
    pivot_longer(-Rep) %>% 
    group_by(name) %>% 
    summarise(Mean = mean(value), 
              Lower = HPDinterval(as.mcmc(value))[1], 
              Upper = HPDinterval(as.mcmc(value))[2]) %>% 
    mutate(Hurricane = ifelse(str_detect(name, "Pre"), "Pre", "Post")) %>% 
    mutate(Comparison = substr(name, nchar(name)-1, nchar(name))) %>% 
    
    mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
    mutate(Prefix = rep(c("Sex: ", rep("Rank: ", 3)), 2)) %>% 
    mutate_at("Comparison", ~factor(paste0(Prefix, .x), 
                                    levels = rev(c("Sex: FM", "Rank: LM", "Rank: LH", "Rank: MH")))) %>% 
    
    ggplot(aes(Comparison, Mean, colour = Hurricane)) +
    geom_point(colour = "black", 
               aes(group = Hurricane), 
               size = 2.5, 
               position = position_dodge(w = 0.7)) + 
    geom_point(position = position_dodge(w = 0.7)) + 
    geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
    geom_errorbar(position = position_dodge(w = 0.7),
                  aes(ymin = Lower, ymax = Upper), 
                  width = 0.3) +
    coord_flip() +
    theme(legend.position = c(0.9, 0.3), legend.justification = c(1, 0)) +
    scale_colour_manual(values = ParasiteColours[c(5, 4)]))

# ggsave("Figures/Figure3.jpeg", units = "mm", height = 120, width = 120)

## Group differences: Model 2 ####

if(!file.exists("Data/Outputs/ComparisonDF2.csv")){
  
  IndivDF %<>% as.data.frame
  
  PredDF <- MakePredictDF(IndivDF %>% dplyr::select(-c(MeanInf, MeanTime)), 
                          
                          HoldNumeric = c("Age", "Strength")
                          # HoldFactor = c(#"Sex" = "F", "Pop" = "K2020",
                          #                "ID" = "00o")
  ) %>% 
    expand.grid
  
  PredDF %<>% filter(ID %in% unique(ID)[1:2], 
                     Pop %in% unique(Pop)[1:2])
  
  PredDF[,paste0("Pred.", 1:1000)] <- 
    INLAFit(Model2$FinalModel, PredDF, 
            FixedCovar = c("Hurricane", "Rank", "Age", "Sex", 
                           paste0("Hurricane:", c("Rank", "Age", "Sex"))), 
            # Return = "Matrix",
            Draw = T, NDraw = 1000) %>% 
    map(as.data.frame) %>% 
    bind_rows() %>% 
    t
  
  PredDF %<>% 
    filter(ID == ID[1], Pop == Pop[1])
  
  LongPred <- 
    PredDF %>% 
    # filter(ID == ID[1], Pop == Pop[1]) %>% 
    pivot_longer(matches("Pred"))
  
  i <- 1
  
  ComparisonDF2 <- 
    1:1000 %>%
    map(function(i){
      
      data.frame(
        
        PreLM = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("L", "M"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreMH = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("M", "H"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreLH = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("L", "H"), 
                                  Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostLM = PredDF %>% filter(Hurricane == "Post", Rank %in% c("L", "M"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostMH = PredDF %>% filter(Hurricane == "Post", Rank %in% c("M", "H"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostLH = PredDF %>% filter(Hurricane == "Post", Rank %in% c("L", "H"), 
                                   Sex == "F") %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PreFM = PredDF %>% filter(Hurricane == "Pre", Rank %in% c("M")) %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff,
        
        PostFM = PredDF %>% filter(Hurricane == "Post", Rank %in% c("M")) %>% 
          dplyr::pull(paste0("Pred.", i)) %>% diff
        
      )
      
    }) %>% 
    bind_rows(.id = "Rep")
  
  ComparisonDF2 %>% write.csv("Data/Outputs/ComparisonDF2.csv", row.names = F)
  
}else ComparisonDF2 <- read.csv("Data/Outputs/ComparisonDF2.csv")

(Model2ComparisonPanel <- 
    ComparisonDF2 %>% 
    pivot_longer(-Rep) %>% 
    group_by(name) %>% 
    summarise(Mean = mean(value), 
              Lower = HPDinterval(as.mcmc(value))[1], 
              Upper = HPDinterval(as.mcmc(value))[2]) %>% 
    mutate(Hurricane = ifelse(str_detect(name, "Pre"), "Pre", "Post")) %>% 
    mutate(Comparison = substr(name, nchar(name)-1, nchar(name))) %>% 
    
    mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
    mutate(Prefix = rep(c("Sex: ", rep("Rank: ", 3)), 2)) %>% 
    mutate_at("Comparison", ~factor(paste0(Prefix, .x), 
                                    levels = rev(c("Sex: FM", "Rank: LM", "Rank: LH", "Rank: MH")))) %>% 
    ggplot(aes(Comparison, Mean, colour = Hurricane)) +
    geom_point(position = position_dodge(w = 0.7)) + 
    geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
    geom_errorbar(position = position_dodge(w = 0.7),
                  aes(ymin = Lower, ymax = Upper), 
                  width = 0.3) +
    coord_flip() +
    scale_colour_manual(values = ParasiteColours[c(5, 4)]))


## Panel B: Rank ####

ComparisonDF1 %>% 
  summarise_at(2:ncol(.), 
               function(a){
                 
                 1-abs(500-sum(abs(diff(table(a>0)))))/500
                 
                 
               })

# individual_timestep_Allrank %<>% 
#   mutate_at("isPost", as.character) %>%
#   mutate_at("isPost", ~factor(CamelConvert(.x), levels = c("Pre", "Post"))) %>% 
#   filter(!is.na(rank))

# individual_timestep_Allrank = individual_timestep[!is.na(individual_timestep$rank),]

# individual_timestep_Allrank$rank = 
#   factor(individual_timestep_Allrank$rank, levels = c("L", "M", "H"))

# individual_timestep_Allrank

Segments <- 
  data.frame(YFrom = c(1200, 1000, 1000)) %>% 
  mutate(YTo = YFrom, XFrom = c(1.1, 1.1, 2.1)) %>% 
  mutate(XTo = c(2.9, 1.9, 2.9), Label = c("***", "***", "NS")) %>% 
  mutate(X = (XFrom + XTo)/2) %>% 
  mutate(Hurricane = "Pre") %>% 
  mutate(Y = c(1200 + 25, 1000 + 25, 1000 + 85)) %>% 
  bind_rows(
    data.frame(YFrom = c(1200, 1000, 1000)) %>% 
      mutate(YTo = YFrom, XFrom = c(1.1, 1.1, 2.1)) %>% 
      mutate(XTo = c(2.9, 1.9, 2.9), Label = c("***", "***", "NS")) %>% 
      mutate(X = (XFrom + XTo)/2) %>% 
      mutate(Y = c(1200 + 25, 1000 + 25, 1000 + 85)) %>% 
      mutate(Hurricane = "Post")) %>% 
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post")))# %>% 
# mutate_at(vars(matches("Y")), log10)

(PanelB <- 
    ggplot(IndivDF) +
    # geom_violin() +
    geom_sina(alpha = 0.1, aes(colour = Hurricane, x = Rank, y = MeanTime)) +
    geom_boxplot(aes(x = Rank, y = MeanTime), width = 0.2, outliers = F) +
    facet_grid(~Hurricane) +
    geom_segment(data = Segments, #inherit.aes = F,
                 aes(y = YFrom, yend = YTo,
                     colour = Hurricane,
                     x = XFrom, xend = XTo)) +
    geom_text(data = Segments, inherit.aes = F,
              aes(y = Y, x = X, 
                  colour = Hurricane,
                  label = Label)) +
    scale_colour_manual(values = ParasiteColours[c(5, 4)]) +
    # scale_y_continuous(breaks = c(0:5*250), 
    #                    name = "Mean infection timestep",
    #                    limits = c(0, 1250)
    # ) +
    theme(strip.background = element_rect(fill = "white", #colour = "dark grey",
                                          colour = "white")) +
    theme(legend.position = "none") + scale_y_log10(breaks = 2^(1:10), 
                                                    name = "Mean infection timestep",
                                                    limits = c(30, 1250))
)

(PanelB <-
    IndivDF %>% 
    ggplot(aes(x = Rank, y = MeanTime)) +
    # geom_jitter(alpha = 0.3, aes(colour = isPost))+
    geom_point(alpha = 0.1, aes(colour = Hurricane))+
    geom_smooth(method = "lm", colour = "black")+
    # facet_grid(~Hurricane) +
    scale_x_continuous(breaks = c(0:6/2), labels = c(0:6*5)) +
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

ggsave("Figures/Figure2.jpeg", units = "mm", height = 200, width = 300)

# Figure 3 ####

xx.strength <- 
  read.csv(file = "Results/contrasts_hurricaneRank_strength.csv")

# xxplot.strength = xx.strength[c(2,4,11,7,9,14),]; 

# xxplot.strength$isPost = c("pre","pre","pre","post","post","post"); xxplot.strength$metric = "strength"

xx.strength %<>% 
  separate(contrast, sep = " - ", into = c("Level1", "Level2")) %>% 
  separate(Level1, sep = " ", into = c("Hurricane1", "Rank1")) %>% 
  separate(Level2, sep = " ", into = c("Hurricane2", "Rank2")) %>% 
  filter(Hurricane1 == Hurricane2)

xx.degree <- 
  read.csv(file = "Results/contrasts_hurricaneRank_degree.csv")

# xxplot.degree = xx.degree[c(2,4,11,7,9,14),]; 

# xxplot.degree$isPost = c("pre","pre","pre","post","post","post"); xxplot.degree$metric = "degree"

xx.degree %<>% 
  separate(contrast, sep = " - ", into = c("Level1", "Level2")) %>% 
  separate(Level1, sep = " ", into = c("Hurricane1", "Rank1")) %>% 
  separate(Level2, sep = " ", into = c("Hurricane2", "Rank2")) %>% 
  filter(Hurricane1 == Hurricane2)

xxplot = rbind(xx.degree %>% mutate(metric = "Degree"), 
               xx.strength %>% mutate(metric = "Strength"))

# xxplot$contrast = c("High-Low", "High - Med", "Low-Med",
#                     "High-Low", "High - Med", "Low-Med", 
#                     "High-Low", "High - Med", "Low-Med",
#                     "High-Low", "High - Med", "Low-Med")

xxplot %<>% 
  mutate_at(c("Rank1", "Rank2"), 
            ~str_replace_all(.x, c("H" = "High", 
                                   "M" = "Med", 
                                   "L" = "Low"))) %>% 
  mutate(contrast = paste0(Rank1, "-", Rank2))

xxplot %>%
  mutate(isPost = Hurricane1 %>% CamelConvert) %>% 
  mutate(isPost = fct_relevel(isPost,
                              "Pre", "Post")) %>%
  
  ggplot(aes(x = factor(contrast, level = c("Low-Med", "High-Med", "High-Low")), 
             y = estimate, group = isPost, colour = isPost)) + 
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), 
                position = position_dodge(w = 0.5),
                width = .1) +
  geom_point(position = position_dodge(w = 0.5), colour = "black", size = 3) +
  geom_point(position = position_dodge(w = 0.5), size = 2) +
  geom_hline(yintercept = 0, 
             linetype = "dashed", 
             color = "black")+ 
  facet_grid(~metric, scales = "free") +
  # theme_light() +
  labs(x = NULL) +
  scale_colour_manual(name = "Hurricane\nstatus", values = ParasiteColours[c(5, 4)]) +
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  theme(legend.position = c(0.95, 0.05),
        legend.justification = c(1, 0)) +
  coord_flip()

ggsave("Figures/Figure3.jpeg", units = "mm", 
       height = 180, width = 180)


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


# Old Figure Stuff #####

EdgeListList <- readRDS("Greg Data/TimeEdges.rds")

EdgeListList %<>%
  separate(DateTime, sep = " ", into = c("Date", "Time"))

EdgeListList %<>% mutate_at("Date", ~lubridate::ymd(.x))

EdgeListList %<>% 
  mutate(Year = str_split(Date, "-") %>% map_chr(1))

EdgeListList %>% 
  count(Rep, Date, Year) %>% 
  ggplot(aes(Rep, Date)) + 
  geom_tile(aes(fill = n)) + 
  coord_flip() +
  facet_wrap(~Year, scales = "free")

DiffDF <- 
  EdgeListList %>% 
  mutate_at("Rep", ~str_replace(.x, "KK", "K")) %>% 
  group_by(Rep, Date, Year) %>% 
  count %>% group_by(Rep, Year) %>% 
  mutate(DateGap = c(0, diff(as.numeric(Date))))

DiffDF %>% 
  filter(DateGap > 0) %>% 
  # mutate_at("Gap", ~(.x + 1)) %>% 
  SinaGraph("Rep", "DateGap", Just = T) + 
  scale_y_log10()

LongMaxes %>% 
  group_by(Rep) %>% 
  summarise_at(c("Mean", "Max"), mean) %>% 
  left_join(DiffDF %>% 
              group_by(Rep) %>% 
              summarise_at("DateGap", mean)) %>% 
  filter(DateGap < 4) %>% 
  ggplot(aes(DateGap, Mean)) +
  # geom_point() + 
  geom_text(aes(label = Rep)) +
  geom_smooth(method = lm)

# Plotting temporal windows ####

Window <- 30

SubEdgeList <-
  EdgeListList

SubEdgeList %<>% 
  group_by(Rep) %>% 
  mutate(NDate = Date - min(Date) + 1)

SubEdgeList %>% 
  ggplot(aes(Rep, NDate)) + 
  geom_point() + 
  coord_flip()

FullDays <- 
  max(SubEdgeList$NDate) - 
  min(SubEdgeList$NDate) - 
  Window

library(tidygraph)

BlankIndivs <- 
  SubEdgeList[,c("From", "To")] %>% unlist %>% unique %>% 
  data.frame(ID = .)

GraphList <- 
  0:FullDays %>% 
  map(function(Day){
    
    SubSubEdgeList <- 
      SubEdgeList %>% 
      filter(Rep == FocalRep) %>% 
      filter(NDate %in% (1:Window + Day)) %>% 
      dplyr::select(From, To) %>% 
      # graph_from_data_frame %>% #as_tbl_graph()
      tbl_graph(edges = ., 
                nodes = BlankIndivs) %>% 
      get.adjacency(sparse = F)
    
  })

library(ggregplot); library(ggraph); library(patchwork)

a <- GraphList[[1]]

Samples <- 
  sample(1:length(GraphList), 10) %>% 
  sort

GraphList[Samples] %>% 
  map(function(a){
    
    a %>% 
      ggraph("stress") + 
      geom_edge_link0() +
      geom_node_point() +
      coord_fixed()
    
  }) %>% ArrangeCowplot()

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



