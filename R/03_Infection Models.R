
# 02b_BISoN Individual Summarising ####

{
  
  library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
  library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)
  library(INLA); library(MCMCglmm)
  
  theme_set(theme_cowplot())
  
  dir_create("Intermediate")
  
}

if(!file.exists("Data/Intermediate/FullIndividualInfectionData.rds")){
  
  # Import the population timecourses ####
  
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
  
  IndivDF <- IndivDFList %>% bind_rows(.id = "Rep")
  
  # Attaching individual traits ####
  
  # load("Data/Intermediate/proximity_data.RData")
  # 
  # IndividualTraits <-
  #   edgelist.all %>%
  #   mutate(Hurricane = factor(ifelse(year < 2018, "Pre", "Post"), levels = c("Pre", "Post"))) %>%
  #   mutate(Pop = paste0(group, year)) %>%
  #   dplyr::select(ID = ID1, Sex = ID1_sex, Rank = ID1_rank, Age = ID1_age, Pop, Hurricane) %>%
  #   unique
  
  TraitList <- "Data/Input" %>% dir_ls(regex = "GroupByYear") %>% map(read.csv)
  
  names(TraitList) <- "Data/Input" %>% list.files(pattern = "GroupByYear") %>% 
    str_remove("Group") %>% str_split("_") %>% map_chr(1)
  
  IndividualTraits <- 
    TraitList %>% map(~.x %>% 
                        rename_all(~str_replace(.x, "dominant", "dominat")) %>% 
                        mutate_at("sex", ~substr(as.character(.x), 1, 1)) %>% # Some sex columns are logical
                        mutate_at(vars(matches("idcode")), as.character)) %>% # Some csvs don't include this column, some are numeric
    bind_rows(.id = "Pop") %>% 
    rename_all(CamelConvert)
  
  IndividualTraits %<>% 
    mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
    mutate(Hurricane = factor(ifelse(Year < 2018, "Pre", "Post"), 
                              levels = c("Pre", "Post")))
  
  IndivDF %<>% left_join(IndividualTraits %>% 
                           dplyr::select(ID = Id, Pop, Hurricane, 
                                         # Rank = Ordinal.rank, 
                                         Rank = Percofsex.dominated, 
                                         Sex, Age) %>% 
                           unique)
  
  IndivDF %>% saveRDS("Data/Intermediate/FullIndividualInfectionData.rds")
  
  IndivDF2 <- 
    IndivDF %>% 
    mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
    filter(Age > 5) %>% 
    filter(Time > 0) %>% 
    mutate_at("Age", ~.x/10) %>% 
    mutate(S_I = str_split(File, "_") %>% map_chr(last) %>% str_remove(".rds") %>% as.numeric %>% multiply_by(10)) %>%
    na.omit
  
}else{
  
  IndivDF <- readRDS("Data/Intermediate/FullIndividualInfectionData.rds")
  
  IndivDF2 <- 
    IndivDF %>% 
    mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
    filter(Age > 5) %>% 
    filter(Time > 0) %>% 
    mutate_at("Age", ~.x/10) %>% 
    mutate(S_I = str_split(File, "_") %>% map_chr(last) %>% str_remove(".rds") %>% as.numeric %>% multiply_by(10)) %>%
    na.omit
  
}

# Analysing averaged infection timestep ####

TestDF <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  # mutate_at("RankNumeric", ~factor(.x, levels = c("L", "M", "H"))) %>%
  mutate_at("Rank", ~.x/100) %>%
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, #RankNumeric, 
           Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            N = n(),
            MeanTime = mean(Time)) %>% 
  ungroup

TestDF %<>% mutate_at("MeanTime", round)

INLAGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)

INLALogGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                          f(ID, model = "iid") + f(Pop, model = "iid"),
                        # family = "poisson",
                        family = "lognormal",
                        control.compute = list(dic = TRUE),
                        data = TestDF)

INLACount <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "poisson",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLAGamma <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "gamma",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLANB <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                 f(ID, model = "iid") + f(Pop, model = "iid"),
               family = "nbinomial",
               control.compute = list(dic = TRUE),
               data = TestDF)

list(INLAGaussian, 
     INLALogGaussian, INLAGamma, INLACount, INLANB) %>% INLADICFig() +
  
  (list(#INLAGaussian, 
    INLALogGaussian, INLAGamma, INLACount, INLANB) %>% Efxplot(Intercept = F))

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% MDIC

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% 
  saveRDS("Data/Intermediate/IndividualInfectionModelList.rds")

IM1 <- INLAModelAdd(Data = TestDF %>% mutate_at("MeanTime", log), 
                    # Family = "loggaussian",
                    Response = "MeanTime",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = paste0("Hurricane:", c("Age", "Sex", "Rank")),
                    # Add = c("Rank", "RankNumeric"),
                    # Clashes = list(c("Rank", "RankNumeric")),
                    Delta = -Inf,
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

IM1 %>% saveRDS("Data/Intermediate/IndividualInfectionModelAdd.rds")

IM1$FinalModel %>% Efxplot(Intercept = F)

IM1$AllModels[[2]] %>% Efxplot(Intercept = F)

IM1$dDIC

# Adding sociality as an explanatory variable ####

load("Data/Intermediate/BisonMetrics.RData")

Strengths <-
  node_strength_all %>% map(colMeans) %>% map(c(reshape2::melt, rownames_to_column)) %>%
  bind_rows(.id = "Pop") %>%
  rename(ID = rowname, Strength = value)

Degrees <-
  node_degree_all %>% map(colMeans) %>% map(c(reshape2::melt, rownames_to_column)) %>%
  bind_rows(.id = "Pop") %>%
  rename(ID = rowname, Degree = value)

TestDF2 <-
  list(TestDF,
       Strengths,
       Degrees) %>%
  reduce(~left_join(.x, .y, by = c("ID", "Pop"))) %>% 
  na.omit #%>% 
  # mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H")))

TestDF2 %<>%
  filter(Age > 0.5)

TestDF2 %<>% mutate_at("Strength", log)

IM1 <- INLAModelAdd(Data = TestDF2 %>% mutate_at("MeanTime", log10), 
                    # Family = "loggaussian",
                    Response = "MeanTime",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = c(paste0("Hurricane:", c("Age", "Sex", "Rank")), 
                            "Degree", "Strength"),
                    Delta = -Inf,
                    Clashes = list(c("Degree", "Strength")),
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

IM1 %>% saveRDS("Data/Intermediate/IndividualInfectionModelAddStrength.rds")

IM1$FinalModel %>% Efxplot
IM1$dDIC

list(IM1$AllModels[[1]], IM1$FinalModel) %>% Efxplot

# Repeatability models ####

TestDF <- 
  IndivDF %>% 
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
            N = n(),
            MeanTime = mean(Time)) %>% 
  ungroup

TestDF %<>% mutate_at("MeanTime", round)

TestDF3 <- TestDF %>% data.frame %>% mutate_at("MeanTime", log10)

TestDF3 %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

MC1 <- MCMCglmm(MeanTime ~ Year + Pop + Hurricane, 
                data = TestDF3, 
                # family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(MeanTime ~ Year + Pop + Hurricane, 
                data = TestDF3, 
                # family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(MeanTime ~ Year + Pop, 
                 data = TestDF3 %>% filter(Hurricane == "Pre"), 
                 # family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(MeanTime ~ Year + Pop, 
                 data = TestDF3 %>% filter(Hurricane == "Post"), 
                 # family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% saveRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

ModelList %>% map(~MCMCRep(.x)) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD) %>% 
  write.csv("Data/Outputs/InfectionRepeatabilityValues.csv", row.names = F)

# Employing numerical rank ####

TestDF <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  # mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H"))) %>% 
  # mutate_at("RankNumeric", ~.x/100) %>% 
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, RankNumeric, Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            N = n(),
            MeanTime = mean(Time)) %>% 
  ungroup

TestDF %<>% mutate_at("MeanTime", round)

INLAGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + RankNumeric) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)

INLALogGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + RankNumeric) + 
                          f(ID, model = "iid") + f(Pop, model = "iid"),
                        # family = "poisson",
                        family = "lognormal",
                        control.compute = list(dic = TRUE),
                        data = TestDF)

INLACount <- inla(MeanTime ~ Hurricane * (Age + Sex + RankNumeric) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "poisson",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLAGamma <- inla(MeanTime ~ Hurricane * (Age + Sex + RankNumeric) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "gamma",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLANB <- inla(MeanTime ~ Hurricane * (Age + Sex + RankNumeric) + 
                 f(ID, model = "iid") + f(Pop, model = "iid"),
               family = "nbinomial",
               control.compute = list(dic = TRUE),
               data = TestDF)

list(INLAGaussian, 
     INLALogGaussian, INLAGamma, INLACount, INLANB) %>% INLADICFig() +
  
  (list(#INLAGaussian, 
    INLALogGaussian, INLAGamma, INLACount, INLANB) %>% Efxplot(Intercept = F))

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% MDIC

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% 
  saveRDS("Data/Intermediate/IndividualInfectionModelList.rds")

IM1 <- INLAModelAdd(Data = TestDF %>% mutate_at("MeanTime", log10), 
                    # Family = "loggaussian",
                    Response = "MeanTime",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = paste0("Hurricane:", c("Age", "Sex", "Rank")),
                    Delta = -Inf,
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

IM1 %>% saveRDS("Data/Intermediate/IndividualInfectionModelAdd.rds")

IM1$FinalModel %>% Efxplot(Intercept = F)

IM1$AllModels[[2]] %>% Efxplot(Intercept = F)

IM1$dDIC

# XXXXXXXXX_Now Defunct ####

# MCMC Repeatability Model ####

TestDF <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H"))) %>% 
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            N = n(),
            MeanTime = mean(Time)) %>% 
  ungroup

TestDF %<>% mutate_at("MeanTime", round)

TestDF %<>% data.frame

TestDF %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

MC1 <- MCMCglmm(MeanTime ~ Hurricane * (Age + Sex + Rank) + Year + Pop, 
                data = TestDF, 
                family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(MeanTime ~ Hurricane * (Age + Sex + Rank) + Year + Pop, 
                data = TestDF, 
                family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(MeanTime ~ Age + Sex + Rank + Year + Pop, 
                 data = TestDF %>% filter(Hurricane == "Pre"), 
                 family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(MeanTime ~ Age + Sex + Rank + Year + Pop, 
                 data = TestDF %>% filter(Hurricane == "Post"), 
                 family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

# ModelList %>% saveRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")
# 
# ModelList %>% map(~MCMCRep(.x)) %>% bind_rows(.id = "Model") %>% 
#   mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
#   mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
#   rename(Lower = lHPD, Upper = uHPD) %>% 
#   write.csv("Data/Outputs/InfectionRepeatabilityValues.csv", row.names = F)

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


# Looping through reps ####

IndivDF2 %>% 
  # filter(File %in% unique1:100) %>% 
  mutate(Rep = str_split(File, "_") %>% map_chr(2) %>% as.numeric) %>% 
  filter(Rep %in% 1:100) %>% 
  ggplot(aes(Time)) + 
  geom_density(aes(colour = as.factor(Rep))) + 
  theme(legend.position = "none")

IndivDF2$Rep <- IndivDF2$File %>% str_split("_") %>% map_chr(2) %>% as.numeric

NIterations <- 100 #IndivDF$Rep %>% max

i <- 1

i <- "Data/Intermediate/IndividualLMER" %>% dir_ls %>% length %>% add(1)

RepModels <- list()

for(i in i:NIterations){
  
  print(i)
  
  TestDF <- IndivDF2 %>% filter(Rep == i)
  
  # INLAGaussian <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
  #                        f(ID, model = "iid") + f(Pop, model = "iid"),
  #                      # family = "poisson",
  #                      control.compute = list(dic = TRUE),
  #                      data = TestDF)
  # 
  # INLALogGaussian <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
  #                           f(ID, model = "iid") + f(Pop, model = "iid"),
  #                         # family = "poisson",
  #                         family = "lognormal",
  #                         control.compute = list(dic = TRUE),
  #                         data = TestDF)
  
  INLACount <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
                      f(ID, model = "iid") + f(Pop, model = "iid"),
                    family = "poisson",
                    control.compute = list(dic = TRUE),
                    data = TestDF)
  
  # INLAGamma <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
  #                     f(ID, model = "iid") + f(Pop, model = "iid"),
  #                   family = "gamma",
  #                   control.compute = list(dic = TRUE),
  #                   data = TestDF)
  # 
  # INLANB <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
  #                  f(ID, model = "iid") + f(Pop, model = "iid"),
  #                family = "nbinomial",
  #                control.compute = list(dic = TRUE),
  #                data = TestDF)
  # 
  # ModelList <- 
  #   list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB)
  
  # ModelList %>% INLADICFig
  
  RepModels[[i]] <- INLACount
  
}

RepModels %>% Efxplot + theme(legend.position = "none")

MCMCEstimates <- 
  RepModels %>% map("summary.fixed") %>% map(c(as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Rep") %>% 
  rename(Var = rowname)

MCMCEstimateSummaries <- 
  MCMCEstimates %>% 
  group_by(Var) %>% 
  summarise(Mean = mean(mean), 
            Lower = HPDinterval(as.mcmc(mean))[1], 
            Upper = HPDinterval(as.mcmc(mean))[2])



# Importing estimates ####

IMList <- 
  "Data/Intermediate/IndividualINLA" %>% dir_ls %>% 
  map(readRDS)

MCMCEstimates <- 
  IMList %>% map("summary.fixed") %>% map(c(as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Rep") %>% 
  rename(Var = rowname)

MCMCEstimateSummaries <- 
  MCMCEstimates %>% 
  group_by(Var) %>% 
  summarise(Mean = mean(mean), 
            Lower = HPDinterval(as.mcmc(mean))[1], 
            Upper = HPDinterval(as.mcmc(mean))[2])

MCMCEstimates %>% 
  ggplot() + 
  geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  geom_sina(aes(Var, mean, colour = Var)) +
  geom_errorbar(data = MCMCEstimateSummaries, 
                aes(x = Var, y = Mean, ymin = Lower, ymax = Upper),
                colour = "black") +
  geom_point(data = MCMCEstimateSummaries, 
             aes(x = Var, y = Mean),
             colour = "black") +
  coord_flip()

# LM below ####

LMList <-
  "Data/Intermediate/IndividualLMER" %>% dir_ls %>% extract(1:1000) %>% 
  map(readRDS)

MCMCEstimates <- 
  LMList %>% 
  map(c(summary, coef, as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "rep") %>% 
  rename(Var = rowname)

MCMCEstimateSummaries <- 
  MCMCEstimates %>% 
  group_by(Var) %>% 
  summarise(Mean = mean(Estimate), 
            Lower = HPDinterval(as.mcmc(Estimate))[1], 
            Upper = HPDinterval(as.mcmc(Estimate))[2])

MCMCEstimates %>% 
  ggplot() + 
  geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  geom_sina(aes(Var, Estimate, colour = Var)) +
  geom_errorbar(data = MCMCEstimateSummaries, 
                aes(x = Var, y = Mean, ymin = Lower, ymax = Upper),
                colour = "black") +
  geom_point(data = MCMCEstimateSummaries, 
             aes(x = Var, y = Mean),
             colour = "black") +
  coord_flip()


# Repeatability of rank ####

TestDF <- 
  IndivDF %>% 
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
            N = n(),
            MeanTime = mean(Time)) %>% 
  ungroup

# TestDF %<>% mutate_at("MeanTime", round)

TestDF3 <- TestDF %>% data.frame# %>% mutate_at("MeanTime", log10)

TestDF3 %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

MC1 <- MCMCglmm(Rank ~ Year + Pop + Hurricane, 
                data = TestDF3, 
                # family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(Rank ~ Year + Pop + Hurricane, 
                data = TestDF3, 
                # family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(Rank ~ Year + Pop, 
                 data = TestDF3 %>% filter(Hurricane == "Pre"), 
                 # family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(Rank ~ Year + Pop, 
                 data = TestDF3 %>% filter(Hurricane == "Post"), 
                 # family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

# ModelList %>% saveRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

ModelList %>% map(~MCMCRep(.x)) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD)

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

# Repeatability of Strength ####

load("Data/Intermediate/BisonMetrics.RData")

Strengths <-
  node_strength_all %>% map(colMeans) %>% map(c(reshape2::melt, rownames_to_column)) %>%
  bind_rows(.id = "Pop") %>%
  rename(ID = rowname, Strength = value)

TestDF2 <-
  list(TestDF,
       Strengths) %>%
  reduce(~left_join(.x, .y, by = c("ID", "Pop"))) %>% 
  na.omit %>% data.frame #%>% 
# mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H")))

TestDF2 %<>%
  filter(Age > 0.5)

TestDF2 %<>% mutate_at("Strength", log)

TestDF2 %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

MC1 <- MCMCglmm(Strength ~ Year + Pop + Hurricane, 
                data = TestDF2, 
                # family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(Strength ~ Year + Pop + Hurricane, 
                data = TestDF2, 
                # family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(Strength ~ Year + Pop, 
                 data = TestDF2 %>% filter(Hurricane == "Pre"), 
                 # family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(Strength ~ Year + Pop, 
                 data = TestDF2 %>% filter(Hurricane == "Post"), 
                 # family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% map("Deviance")

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
