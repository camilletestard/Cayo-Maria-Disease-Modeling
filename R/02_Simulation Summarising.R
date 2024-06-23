
# 02a_Simulation Summarising ####

{
  
  library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
  library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)
  library(MCMCglmm)
  
  theme_set(theme_cowplot())
  
  FileList <- 
    "Data/Outputs/BISoN" %>% 
    dir_ls()
  
  names(FileList) <- "Data/Outputs/BISoN" %>%
    list.files
  
}

# Population summaries ####

if(!file.exists("Data/Outputs/PopulationTimes.rds")){
  
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
  
  OutputDF %>% saveRDS("Data/Outputs/PopulationTimes.rds")
  
} else OutputDF <- readRDS("Data/Outputs/PopulationTimes.rds")

# Converting to timesteps ####

if(!file.exists("Data/Output/PopulationTimeSteps.rds")){
  
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
  
  TimeStepDF %>% saveRDS("Data/Outputs/PopulationTimeSteps.rds")
  
} else  TimeStepDF <- readRDS("Data/Outputs/PopulationTimeSteps.rds")

# Running model examining the processes determining epidemic spread ####

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

TestDF <- OutputDF %>% mutate(PI = P_I)

# TestDF <- data.frame(X = rnorm(1000), Y = rnorm(1000))

# MCMCglmm::MCMCglmm(Y~X, data = TestDF)
# 
# INLA::inla(Mean ~ P_I + PostMaria + f(Rep, model = "iid"), data = TestDF)

MCMC1 <- MCMCglmm(Mean ~ P_I + PostMaria, random =~Rep, data = TestDF)

MCMC1 %>% saveRDS("Model.rds")
