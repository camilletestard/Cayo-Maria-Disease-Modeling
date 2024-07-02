
# 03b_Individual Proximity Differences.R ####

# Initially written by Camille, rewritten and pipelined by Greg

# Model proximity network with bison - get a distribution of possible networks from the observed data

# Regression from imputed dataset to test the change in proximity over time according to age, sex and rank categories

# C. Testard February 2024

{
  
  #Load libraries
  #data wrangling
  library(dplyr)
  library(forcats)
  
  #Modelling
  library(lme4)
  library(lmerTest)
  library(broom.mixed)
  library(emmeans)
  
  #For plotting
  library(ggplot2)
  library(sjPlot)
  library(sjmisc)
  
  library(cowplot)
  
  theme_set(theme_cowplot())
  
  dir_create("Data/Intermediate/SocialLMER")
  
  load("Data/Intermediate/BisonMetrics.RData")
  
}

#Set group year list

group = c("F","KK","F","HH","F","V","R","KK","R","V","F","HH","F","KK","V","V","KK","S","V","F","V","TT","V","F")

years = c(2013, 2013,2014,2014,2015,2015,2015,2015,
          2016,2016,2016,2016,2017,2017,2017,
          2018,2018, 2019, 2019,2021,2021,2022,2022,2022)

groupyears = paste0(group,years)

DataList <- list()

for (gy in 1:length(groupyears)){
  
  # data_path<-'~/Documents/GitHub/Cayo-Maria-Survival/Data/Data All Cleaned/BehavioralDataFiles/'
  data_path <- 'Data/Input/'
  data <- read.csv(paste0(data_path,"Group", groupyears[gy], "_GroupByYear.txt")) #load meta data
  
  data = data[data$hrs.focalfollowed>0,]
  
  data$group = group[gy]
  data$year = years[gy]; 
  data$isPost = ifelse(years[gy]<2018, "pre","post")
  data$isPost.year = ifelse(years[gy]<2018,"pre",paste(data$isPost, data$year,sep='.'))
  data$id = as.factor(data$id)
  data$sex[data$sex == "FALSE"] = "F"; data$sex = as.factor(data$sex);
  data$id.year = paste(data$id, data$year,sep='.')
  
  DataList[[gy]] <- data
  
}

# STEP 1: CREATE IMPUTED DATASSET ####

# Extract change in proximity from each draw to create imputed data set:

imputed.data = list(); imputed.data.prepost = list(); i = 1; gy = 1;

imputed.density = list()

StrengthModelList <- DegreeModelList <- list()

num_iter <- 100

i <- "Data/Intermediate/SocialLMER" %>% dir_ls %>% length %>% add(1)

i <- 1

for (i in i:num_iter){
  
  print(i)
  
  data.all <- data.frame(); density.iter = data.frame(matrix(NA, nrow = num_iter)); names(density.iter)="dens"
  
  density.all <- data.frame();
  
  for (gy in 1:length(groupyears)){
    
    data <- DataList[[gy]]
    
    strength <- node_strength_all[[groupyears[gy]]][i,]
    degree <- node_degree_all[[groupyears[gy]]][i,]
    node.id <- node_ids[[groupyears[gy]]]
    data$prox.strength<-as.numeric(strength[match(data$id, node.id)])
    data$prox.degree<-as.numeric(degree[match(data$id, node.id)])
    data.final <- data[,c("id", "sex", "age", "ordinal.rank", "group", "year", "isPost", "isPost.year", "id.year", 
                          "prox.strength", "prox.degree")]
    
    data.final %<>% na.omit
    
    data.all <- rbind(data.all, data.final)
    
    density.iter$dens = density_all[[groupyears[gy]]]
    density.iter$group = group[gy]
    density.iter$year = years[gy]
    density.iter$isPost = ifelse(years[gy]<2018, "pre","post")
    density.iter$isPost.year = ifelse(years[gy]<2018,"pre",paste(data$isPost, data$year,sep='.'))
    
    density.all <- rbind(density.all, density.iter)
    
  }
  
  data.all %<>% mutate_at("prox.strength", log10)
  
  #Standardize
  density.all$std.dens <-(density.all$dens - mean(density.all$dens))/sd(density.all$dens)
  data.all$std.prox.strength = (data.all$prox.strength-mean(data.all$prox.strength))/sd(data.all$prox.strength)
  data.all$std.prox.degree = (data.all$prox.degree-mean(data.all$prox.degree))/sd(data.all$prox.degree)
  
  #Create year factor
  data.all$year.factor = as.factor(data.all$year)
  
  #Adjust levels
  data.all$group<-as.factor(data.all$group)
  data.all = data.all %>%
    mutate(isPost = fct_relevel(isPost,
                                "pre", "post"))# %>%
    # mutate(isPost.year = fct_relevel(isPost.year,
    #                                  "pre", "post.2018", "post.2019", "post.2021", "post.2022"))
  #"post.2018","pre","post.2019","post.2021", "post.2022"))
  
  density.all = density.all %>%
    mutate(isPost = fct_relevel(isPost,
                                "pre", "post")) #%>%
    # mutate(isPost.year = fct_relevel(isPost.year,
    #                                  "pre", "post.2018","post.2019","post.2021", "post.2022"))
  
  #"post.2018","pre","post.2019","post.2021", "post.2022"))
  
  
  imputed.data[[i]] = data.all
  imputed.density[[i]] = density.all
  
  
  #If only consider IDs present both pre and post
  id.isPost = table(data.all$id, data.all$isPost)
  id.year = table(data.all$id, data.all$year); 
  id_pre = row.names(as.data.frame(which(id.isPost[,"pre"]>0))); id_post = row.names(as.data.frame(which(id.isPost[,"post"]>0)))
  id.PreAndPost = as.data.frame(row.names(as.data.frame(which(id.isPost[,"pre"]>0 & id.isPost[,"post"]>0)))); names(id.PreAndPost)="id"
  id.PreAndPost$group = data.all$group[match(id.PreAndPost$id, data.all$id)]
  table(id.PreAndPost$group)
  data.prepost = subset(data.all, id %in% id.PreAndPost$id)
  
  imputed.data.prepost[[i]] <- data.prepost
  
  # Running internal Bayesian models 
  
  # mdl.strength <- lmer(std.prox.strength ~ 1 + isPost*(ordinal.rank + age + sex) + (1|id/group),
  #                      data = data.all)
  
  # mdl.strength <- MCMCglmm(std.prox.strength ~ 1 + isPost*(ordinal.rank + age + sex), 
  #                          random = ~ id + group,
  #                          data = data.all)
  
  data.all %<>% 
    rename(Rank = ordinal.rank, Pop = group) %>% 
    rename(ID = id) %>% 
    rename_all(CamelConvert)
  
  mdl.strength <- inla(Std.prox.strength ~ IsPost * (Age + Sex + Rank) + 
                         f(ID, model = "iid") + f(Pop, model = "iid"),
                       # family = "poisson",
                       control.compute = list(dic = TRUE),
                       data = data.all)
  
  StrengthModelList[[i]] <- mdl.strength
  
  # mdl.strength %>% saveRDS(glue::glue("Data/Intermediate/SocialLMERStrength/{i}.rds"))
  
  # mdl.degree <- lmer(std.prox.degree ~ 1 + isPost*(ordinal.rank + age + sex) +
  #                      (1|id/group),
  #                    data = data.all)
  
  # mdl.degree <- MCMCglmm(std.prox.degree ~ 1 + isPost*(ordinal.rank + age + sex),
  #                        random = ~ id + group,  
  #                        data = data.all)
  
  mdl.degree <- inla(Std.prox.degree ~ IsPost * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = data.all)
  
  DegreeModelList[[i]] <- mdl.degree
  
  # mdl.degree %>% saveRDS(glue::glue("Data/Intermediate/SocialLMERDegree/{i}.rds"))
  
  # pairwise comparisons
  
  # xx.strength <- as.data.frame(summary(emmeans(mdl.strength, pairwise ~ isPost*ordinal.rank))$contrasts)
  
}

StrengthModelList[1:10] %>% Efxplot

DegreeModelList[1:10] %>% Efxplot

StrengthModelDifferences <- 
  StrengthModelList %>% 
  map(function(mdl.strength){
    
    as.data.frame(summary(emmeans(mdl.strength, pairwise ~ isPost*ordinal.rank))$contrasts)
    
  }) %>% 
  bind_rows(.id = "Rep")

StrengthModelDifferences

StrengthModelDifferences <- 
  StrengthModelList %>% 
  map(function(mdl.strength){
    
    as.data.frame(summary(emmeans(mdl.strength, pairwise ~ isPost*ordinal.rank))$contrasts)
    
  }) %>% 
  bind_rows(.id = "Rep")



# Create imputed dataset to use in frequentist survival models

imp <- miceadds::datalist2mids(imputed.data)
imp.prepost <- miceadds::datalist2mids(imputed.data.prepost)

# STEP 2: Run downstream analyses such as regression ####

# This step propagates the uncertainty from step 1 through subsequent analyses such as regressions

### 1. Test the long-term effects of Hurricane Maria on proximity ###

# setwd("~/Documents/Github/Cayo-Maria-Disease-Modeling/Results")

mdl.strength.proxPrePost <- with(imp, lmer(std.prox.strength ~ 1 + isPost*ordinal.rank + age + sex 
                                           +(1|group) + (1|id)))

summary(mice::pool(mdl.strength.proxPrePost))

#plot model for one iteration

data.imp = imputed.data[[1]]

mdl.strength<-lmer(std.prox.strength~ 1+ isPost*ordinal.rank +age+ sex 
                   +(1|id/group) , data=data.imp)
summary(mdl.strength)

#pairwise comparisons

xx.strength <- as.data.frame(summary(emmeans(mdl.strength, pairwise ~ isPost*ordinal.rank))$contrasts)

write.csv(xx.strength, file = "Results/contrasts_hurricaneRank_strength.csv")

#visualize model

plot_model(mdl.strength); ggsave("Figures/infection_individual_factors.pdf")

# RG<- ref_grid(mdl.strength, at = list(ordinal.rank = c("L","M","H")))
# emmip(RG, isPost~ordinal.rank, style = "factor")

xxplot.strength = xx.strength[c(2,4,11,7,9,14),]; 

xxplot.strength$isPost = c("pre","pre","pre","post","post","post"); xxplot.strength$metric = "strength"


#Proximity degree
mdl.degree.proxPrePost <- with(imp, lmer(std.prox.degree~ 1+ isPost*ordinal.rank + age + sex+
                                           (1|id/group)))
summary(mice::pool(mdl.degree.proxPrePost))

mdl.degree<-lmer(std.prox.degree~ 1+ isPost*ordinal.rank +age+ sex 
                 +(1|id/group), data=data.imp)

summary(mdl.degree)

xx.degree <- as.data.frame(summary(emmeans(mdl.degree, pairwise ~ isPost*ordinal.rank))$contrasts)

write.csv(xx.degree, file = "Results/contrasts_hurricaneRank_degree.csv")

# RG<- ref_grid(mdl.degree, at = list(ordinal.rank = c("L","M","H")))
# emmip(RG, isPost~ordinal.rank, style = "factor")

# Plotting ####

xx.strength <- 
  read.csv(file = "Results/contrasts_hurricaneRank_degree.csv")

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
  
  ggplot(aes(x = factor(contrast, level = c("Low-Med", "High - Med", "High-Low")), 
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
  coord_flip()

ggsave("Figures/NetworkMetric_byRank_byHurricane.pdf")

#Save
# setwd("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data/")

save(imp, mdl.strength.proxPrePost, mdl.degree.proxPrePost, 
     file = "Data/Outputs/ProxMdlOutput_PrePost.RData")

#load("ProxMdlOutput_PrePost.RData")


