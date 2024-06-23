#bison_proximity_networks.R
#Model proximity networks for Cayo data with bison
#This will allow us to capture uncertainty in downstream simulations
#C. Testard October 2022

# devtools::install_github("JHart96/bisonR")

{
  
  library(bisonR)
  library(dplyr)
  library(ggplot2)
  library(fs)
  
  dir_create("Intermediate")
  
  edgelist.all <- readRDS("Intermediate/proximity_data.rds")
  
  edgelist.all$groupyear = paste0(edgelist.all$group, edgelist.all$year)
  
  #Set group year list
  group = c("F","KK","F","HH","F","V","R","KK","R","V","F","HH","F","KK","V","V","KK","S","V","F","V","TT","V","F")
  years = c(2013, 2013, 2014, 2014, 2015, 2015, 2015, 2015,
            2016, 2016, 2016, 2016, 2017, 2017, 2017,
            2018, 2018, 2019, 2019, 2021, 2021, 2022, 2022, 2022)
  
  groupyears = paste0(group, years)
  
  gy = 1; posterior.el = list(); density_samples = list() #initalize edgelist posterior samples
  
}

for (gy in 1:length(groupyears)){ #for all groups and years
  
  #####################################################
  ## Fit network for binary data  (i.e., proximity) ###
  #####################################################
  # This step fits an edge weight model using conjugate binary priors, 
  # capturing uncertainty in social network edges.
  
  priors <- get_default_priors("binary_conjugate")#define priors
  priors$edge <- "beta(0.1,0.1)" # set less flat priors
  #prior_check(priors, "binary_conjugate") #check prior make biological sense
  
  el = edgelist.all[edgelist.all$groupyear==groupyears[gy],] #Extract edgelist for the group and year
  el$year=factor(el$year)
  
  unqids = unique(el$ID1)
  node.list = el[match(unqids, el$ID1),c("ID1","ID1_obseff_duration", #extract nodelist with node traits
                                         "ID1_obseff_samples","ID1_sex","ID1_age",
                                         "ID1_rank","group","year","isPost","groupyear")]
  
  #Fit network with bison
  fit.el <- bison_model(
    (count | total_samples) ~ dyad(ID1, ID2),
    data=el,
    model_type = "binary_conjugate", #Count conjugate for aggression data 
    priors=priors
  )
  
  #Check model fit
  #plot_predictions(fit.el, num_draws=20)
  #Left plot: The predictions from the model are shown in blue and the real data
  #are shown in black. Ideally the blue lines should be distributed around the
  #black line, indicating the real data are among the possible predictions of the model.
  #Right plot: comparison of predictions against point estimates of edge weights
  #We're happy with current predictions
  
  # #Summary of edge weights and their credible intervals
  # summary(fit.el)
  # 
  # #Visualize network with uncertainty estimates
  # plot_network(fit.el, lwd=5)
  
  #Draw 1000 edgelists from the posterior of the network
  samples.post<-draw_edgelist_samples(fit.el, num_draws=1000)
  samples.post[,c(3:ncol(samples.post))] <- plogis(as.matrix(samples.post[,c(3:ncol(samples.post))])) #Convert back to edge weights
  posterior.el[[groupyears[gy]]] <- samples.post #save
  
  
  #Draw network metrics from posterior: density and modularity
  density_samples[[groupyears[gy]]] <-extract_metric(fit.el, "global_density", num_draws =1000)
  #Soon add modularity...
  
  print(groupyears[gy])
  
}

#Save draws from posterior
# setwd("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data/")

save(posterior.el, density_samples, groupyears, 
     file = "Data/Intermediate/BisonFittedNetworks.RData")

