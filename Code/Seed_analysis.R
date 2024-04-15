#Seed analysis
library(lmerTest)

#Load data
setwd("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data/")
seedID = readRDS("TestDF.rds")

mdl<-lmer(Mean ~ P_I + PostMaria * (Seed.Sex + Seed.Rank + Seed.Age) + (1|Population), data = seedID)
summary(mdl)
plot_model(mdl);

plot(density(seedID$P_I))
