
# 03c_BISoN Metrics.R ####

#Model proximity networks with bison
#Proximity was collected in scans across all years and is a binary data type.

#C. Testard October 2022

{
  
  #Load libraries
  #bison package
  library(bisonR)
  library(dplyr)
  library(brms)
  library(igraph)
  
  #Set seed for reproducibility
  set.seed(1234)
  
  #Load data
  # data_path = "~/Documents/GitHub/Cayo-Maria-Survival/Data/R.Data/"
  
  source("R/Functions/AdjacencyFromEdgeList.R")
  
  load("Data/Intermediate/proximity_data.RData")
  
  #load(paste0(data_path,"BisonProximity.RData"))
  
}

edgelist.all$groupyear = paste0(edgelist.all$group, edgelist.all$year)

#Set group year list
# group = c("F","V","KK","V","F","F","KK","V","V","KK","S","V","F","V")
# years = c(2015,2015,2015,
#           2016,2016,2017,2017,2017,
#           2018, 2018,2019, 2019,2021,2021)
# groupyears = c("F2015","V2015","KK2015",
#                "V2016","F2016","F2017",
#                "KK2017","V2017","V2018","KK2018",
#                "S2019","V2019","F2021","V2021")

group = c("F","KK","F","HH","F","V","R","KK","R","V","F","HH","F","KK","V","V","KK","S","V","F","V","TT","V","F")
years = c(2013, 2013,2014,2014,2015,2015,2015,2015,
          2016,2016,2016,2016,2017,2017,2017,
          2018,2018, 2019, 2019,2021,2021,2022,2022,2022)
groupyears = c("F2013","KK2013","F2014","HH2014","F2015","V2015","R2015","KK2015",
               "R2016","V2016","F2016","HH2016","F2017","KK2017","V2017",
               "V2018","KK2018","S2019","V2019","F2021","V2021","TT2022","V2022","F2022")

#Initialize lists
gy=24; node_strength_all = list(); node_degree_all = list(); density_all = list(); node_ids=list()

#Set parameters

num_iter = 1000

edge_thresh=0.0001

for (gy in 1:length(groupyears)){ #for all group-years
  
  # STEP 1: Fit network for binary data  (i.e., proximity) ####
  
  # This step fits an edge weight model using conjugate priors, 
  # capturing uncertainty in social network edges
  
  priors <- get_default_priors("binary_conjugate")#define priors
  priors$edge <- "beta(0.1,0.1)" # set less flat priors
  #prior_check(priors, "binary_conjugate")
  
  #Load observed proximity edgelist
  el = edgelist.all[edgelist.all$groupyear==groupyears[gy],]
  el$year=factor(el$year)
  el$weight = el$count/el$total_samples
  
  #Fit network with bison binary model
  fit.el <- bison_model(
    (count | total_samples) ~ dyad(ID1, ID2),
    data = el,
    model_type = "binary_conjugate", #Count conjugate for aggression data 
    priors = priors
  )
  #plot_predictions(fit.el, num_draws=20)
  
  #Check model fit
  # #Summary of edge weights and their credible intervals
  # summary(fit.el)
  # 
  # #Visualize network with uncertainty estimates
  # plot_network(fit.el, lwd=5)
  
  #Draw n edgelists from the posterior of the network
  samples.post<-draw_edgelist_samples(fit.el, num_draws=num_iter)
  samples.post[,c(3:ncol(samples.post))]<-plogis(as.matrix(samples.post[,c(3:ncol(samples.post))])) #Convert back from log scale
  #posterior.el[[groupyears[gy]]] <- samples.post #save
  
  #Extract node names
  nodes.id = names(fit.el$node_to_idx)
  
  #Draw network metrics from posterior (global, nodal and edge_weight)
  igraph_strength_all=data.frame(matrix(NA, ncol = length(nodes.id), nrow=0));
  igraph_degree_all=data.frame(matrix(NA, ncol = length(nodes.id), nrow=0)); i=1
  igraph_density_all=vector()
  
  for(i in 1:num_iter){ #For all network instantiations
    weightedEL<-samples.post[,c(1:2,2+i)] #extract edge list from posterior
    names(weightedEL)=c("node_1","node_2","draw")
    weightedEL[weightedEL[,3]<edge_thresh,3]=0#when edge strength is below a certain threshold, 
    #consider that edge to be 0
    
    #Create adjacency matrix from edgelist
    adjMat = AdjacencyFromEdgelist(weightedEL)# create adjacency matrix based on edge list.
    data = adjMat[["adjacency"]]; rownames(data) = adjMat[["nodelist"]]; colnames(data) = adjMat[["nodelist"]]
    #read adjacency matrix into igraph
    m=as.matrix(data) # coerces the data set as a matrix
    am.g=graph.adjacency(m,mode= "undirected",weighted=T) # this will create an directed 'igraph object'. Change qualifiers to make "undirected" or unweighted (null)
    
    #Extract metrics from network instantiation i
    igraph_strength<-strength(am.g); #node strength
    igraph_degree<-degree(am.g); #node degree
    igraph_density<-edge_density(am.g);#network density
    
    # #Compute nodal social selectivity
    # unqIDs = unique(c(as.character(weightedEL$node_1),as.character(weightedEL$node_2)))
    # cv_id=data.frame(unqIDs); cv_id$nonzero=NA; cv_id$cv = NA; i=1
    # for(i in 1:length(unqIDs)){
    #   idx=which(weightedEL$node_1==unqIDs[i]|weightedEL$node_2==unqIDs[i])
    #   cv_id$nonzero[i]=length(which(weightedEL$draw[idx]!=0))/length(idx)
    #   
    #   weightedEL_nonzero = weightedEL[weightedEL$draw!=0,]
    #   idx=which(weightedEL_nonzero$node_1==unqIDs[i]|weightedEL_nonzero$node_2==unqIDs[i])
    #   cv_id$cv[i]=var(weightedEL_nonzero$draw[idx])/mean(weightedEL_nonzero$draw[idx])
    # }
    # hist(cv_id$cv)
    
    #Save over all network instantiations
    igraph_strength_all=rbind(igraph_strength_all, igraph_strength)
    igraph_degree_all=rbind(igraph_degree_all, igraph_degree)
    igraph_density_all=rbind(igraph_density_all, igraph_density)
    
    print(i)
  }
  names(igraph_strength_all)=nodes.id; names(igraph_degree_all)=nodes.id;  #Assign node id
  
  node_ids[[groupyears[gy]]] <- nodes.id
  node_strength_all[[groupyears[gy]]] <- igraph_strength_all
  node_degree_all[[groupyears[gy]]] <- igraph_degree_all
  density_all[[groupyears[gy]]] <- igraph_density_all
  
  print(groupyears[gy]) 
  
}

# setwd("~/Documents/GitHub/Cayo-Maria-Survival/Data/R.Data/")

save(node_ids, node_degree_all, node_strength_all, density_all, num_iter, edge_thresh, 
     file = "Data/Intermediate/BisonProximity.RData")


