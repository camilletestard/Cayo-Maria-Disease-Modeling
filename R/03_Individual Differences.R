
# 03_Individual Differences.r ####

# Initially written by Camille, rewritten and pipelined by Greg

#Load library
library(lmerTest)
library(lme4)
library(ggplot2)
library(readr)
library(sjmisc)
library(sjPlot)
library(igraph)
library(emmeans)

library(fs)

dir_create("Results")

#####################
#Load and format data
#####################

# setwd("~/Documents/Github/Cayo-Maria-Disease-Modeling/Data/R.Data/")

individual_timestep = readRDS("Data/Outputs/IndividualTimesteps.rds")
#individual_timestep = readRDS("IndividualTimesteps_PInf.rds")

individual_timestep$year = parse_number(individual_timestep$Pop)
individual_timestep$group = substr(individual_timestep$Pop,1,1)
individual_timestep$isPost = "pre"; 
individual_timestep$isPost[individual_timestep$year>2017] = "post"
#individual_timestep$age = scale(individual_timestep$age)

load("Data/Intermediate/proximity_data.RData")

edgelist.all$weight = edgelist.all$count/edgelist.all$total_samples
edgelist.all$id1.groupyear = paste0(edgelist.all$ID1, edgelist.all$group, edgelist.all$year)
edgelist.all$id2.groupyear = paste0(edgelist.all$ID2, edgelist.all$group, edgelist.all$year)
individual_timestep$id.groupyear = paste0(individual_timestep$ID, individual_timestep$Pop)

#Add sex
individual_timestep$sex = edgelist.all$ID1_sex[match(individual_timestep$id.groupyear, edgelist.all$id1.groupyear)]
individual_timestep$sex2 =  edgelist.all$ID2_sex[match(individual_timestep$id.groupyear, edgelist.all$id2.groupyear)]
individual_timestep$sex[is.na(individual_timestep$sex)] = individual_timestep$sex2[is.na(individual_timestep$sex)]
individual_timestep$sex[individual_timestep$sex=="FALSE"]="F"
individual_timestep$sex2=NULL

#Add rank
individual_timestep$rank = edgelist.all$ID1_rank[match(individual_timestep$id.groupyear, edgelist.all$id1.groupyear)]
individual_timestep$rank2 =  edgelist.all$ID2_rank[match(individual_timestep$id.groupyear, edgelist.all$id2.groupyear)]
individual_timestep$rank[is.na(individual_timestep$rank)] = individual_timestep$rank2[is.na(individual_timestep$rank)]
individual_timestep$rank2=NULL

#Add continuous rank
individual_timestep$PercRank = edgelist.all$ID1_PercRank[match(individual_timestep$id.groupyear, edgelist.all$id1.groupyear)]
individual_timestep$PercRank2 =  edgelist.all$ID2_PercRank[match(individual_timestep$id.groupyear, edgelist.all$id2.groupyear)]
individual_timestep$PercRank[is.na(individual_timestep$PercRank)] = individual_timestep$PercRank2[is.na(individual_timestep$PercRank)]
individual_timestep$PercRank2=NULL

#Add age
individual_timestep$age = edgelist.all$ID1_age[match(individual_timestep$id.groupyear, edgelist.all$id1.groupyear)]
individual_timestep$age2 =  edgelist.all$ID2_age[match(individual_timestep$id.groupyear, edgelist.all$id2.groupyear)]
individual_timestep$age[is.na(individual_timestep$age)] = individual_timestep$age2[is.na(individual_timestep$age)]
individual_timestep$age2=NULL

individual_timestep$sex=as.factor(individual_timestep$sex)
individual_timestep$age.scale = scale(individual_timestep$age)
individual_timestep$isPost=factor(individual_timestep$isPost, levels =c("pre","post"))
#individual_timestep$rank=factor(individual_timestep$rank, levels =c("M","L","H"))

setwd(here::here())

individual_timestep_keepOriginal = individual_timestep; #Keep original data frame #individual_timestep= individual_timestep_keepOriginal

#individual_timestep = individual_timestep[individual_timestep$S_I_Category=="Med",] #if load data with multiple pInfection than select Medium

individual_timestep = individual_timestep[individual_timestep$age>5,] #Remove individuals <6 yo

individual_timestep %>% saveRDS("Data/Outputs/IndividualTimestepsCamille.rds")

plot(density(individual_timestep$MeanTime))
qplot(individual_timestep$MeanTime)

# Run linear models ####

#mdl <- lmer(MeanTime ~ S_I_Category + isPost*age + isPost*sex+ isPost*rank +(1|ID)+ (1|group), individual_timestep)

mdl <- lmer(MeanTime ~ isPost*age.scale + isPost*sex+ isPost*rank +(1|ID)+ (1|group), individual_timestep)

saveRDS(mdl, file = "Outputs/LevellingModel.rds")

summary(mdl)

# only consider cases where rank is known ####
individual_timestep_nona<-individual_timestep[complete.cases(individual_timestep$rank), ] 
mdl2<-lmer(MeanTime ~ isPost*age.scale + isPost*sex+ isPost*rank +(1|group/ID), individual_timestep_nona)
mdl2_null<-lmer(MeanTime ~ 1 +(1|group/ID), individual_timestep_nona)

# only consider cases where rank is known & rank is continuous ####

individual_timestep_nona<-individual_timestep[complete.cases(individual_timestep$PercRank), ] 
mdl3<-lmer(MeanTime ~ isPost*age.scale + isPost*sex+ isPost*PercRank +(1|group/ID), individual_timestep_nona)
mdl3_null<-lmer(MeanTime ~ 1 +(1|group/ID), individual_timestep_nona)

# full-null model comparison to evaluate overall effect of predictors and avoid cryptic multiple testing

anova(mdl2_null, mdl2, test="Chisq")
drop1(mdl2, test="Chisq")#interaction between rank and hurricane significant
summary(mdl2)

# model performance ####

performance::check_model(mdl2)
plot_model(mdl); ggsave("Figures/infection_individual_factors.pdf")
tab_model(mdl); ggsave("Figures/infection_individual_factors_table.pdf")

# pairwise comparisons ####

xx <- as.data.frame(summary(emmeans(mdl2, pairwise ~ isPost*rank))$contrasts)

# RG<- ref_grid(mdl2, at = list(rank = c("L","M","H")))
# emmip(RG, isPost~rank, style = "factor")

write.csv(xx, file = "Results/contrasts_hurricaneRank.csv")

# Plot data for visualization ####

ggplot(individual_timestep, aes(x=MeanTime, color=isPost))+
  geom_density()

ggplot(individual_timestep, aes(x=age, color=isPost))+
  geom_density()

individual_timestep$year.factor = as.factor(individual_timestep$year)


#Change in inter-individual differences pre-to-post hurricane ####
# setwd("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Results/")

## Rank ####

#For visualization only consider individuals who are High and Low ranking (remove Med)
#individual_timestep_L.Hrank = individual_timestep[individual_timestep$rank=="L"|individual_timestep$rank=="H",]

individual_timestep_Allrank = individual_timestep[!is.na(individual_timestep$rank),]

individual_timestep_Allrank$rank = factor(individual_timestep_Allrank$rank, levels = c("L","M","H"))

ggplot(individual_timestep_Allrank, aes(x=rank, y=MeanTime)) +
  geom_violin() +
  geom_boxplot(width=0.2) +
  facet_grid(~isPost) +
  theme_light() + 
  ylim(0, 900)

ggsave("Results/Time2Infection_byRank_HurrStatus_LMH.pdf")

post.L = mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="L" & individual_timestep_Allrank$isPost=="post"])
pre.L=mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="L" & individual_timestep_Allrank$isPost=="pre"])
post.M=mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="M" & individual_timestep_Allrank$isPost=="post"])
pre.M=mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="M" & individual_timestep_Allrank$isPost=="pre"])
post.H=mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="H" & individual_timestep_Allrank$isPost=="post"])
pre.H=mean(individual_timestep_Allrank$MeanTime[individual_timestep_Allrank$rank=="H" & individual_timestep_Allrank$isPost=="pre"])

pre.H-pre.L; post.H-post.L; 

#For sex

ggplot(individual_timestep, aes(x = sex, y = MeanTime))+
  geom_violin() +
  geom_boxplot(width=0.2) +
  facet_grid(~isPost)

ggsave("Results/Time2Infection_bySex_HurrStatus.pdf")

post.M = mean(individual_timestep$MeanTime[individual_timestep$sex=="M" & individual_timestep$isPost=="post"])
pre.M=mean(individual_timestep$MeanTime[individual_timestep$sex=="M" & individual_timestep$isPost=="pre"])
post.F=mean(individual_timestep$MeanTime[individual_timestep$sex=="F" & individual_timestep$isPost=="post"])
pre.F=mean(individual_timestep$MeanTime[individual_timestep$sex=="F" & individual_timestep$isPost=="pre"])

# For age

ggplot(individual_timestep, aes(x = age, y = MeanTime)) +
  geom_jitter(alpha = 0.3)+
  geom_smooth(method = "lm")+
  facet_grid(~isPost)+
  theme_light()

ggsave("Results/Time2Infection_byAge_HurrStatus.pdf")

cor.test(individual_timestep$age[individual_timestep$isPost == "pre"], 
         individual_timestep$MeanTime[individual_timestep$isPost == "pre"])

cor.test(individual_timestep$age[individual_timestep$isPost == "post"], 
         individual_timestep$MeanTime[individual_timestep$isPost == "post"])

ggplot(individual_timestep, aes(x = sex, y = MeanTime)) +
  geom_violin() +
  geom_boxplot() +
  geom_jitter(width=0.1)

ggplot(individual_timestep, aes(x = rank, y = MeanTime)) +
  geom_violin() +
  geom_boxplot() +
  geom_jitter(width=0.1)



