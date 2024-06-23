#Plot output of disease simulation
# C. Testard, Oct 3rd 2022

#Load libraries
library('ggplot2')
library(plotrix)
library(matrixStats)

#Load output of simulation 
setwd("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data")
load("macaque_disease_simulation.RData")

############################
## Plot infection curves

#Pool data across years
years = c("2015","2016","2017","2018","2019","2021")

data.all=data.frame()
for (yr in c(1:6)){
  df = simu_res_pinf[[yr]][['0.1']][[2]] #for p_infection = 0.1
  mean.year = rowMeans(df)
  sd.year = rowSds(as.matrix(df))
  sem.year =apply(df, 1, std.error)
  data = data.frame(mean=mean.year, sd =sd.year, sem=sem.year, timestep=1:nrow(df),year=years[yr])
  data.all = rbind(data.all,data)
}

#Plot mean % population infected across simulation runs with SEMs
ggplot(data = data.all, aes(x=timestep, y = mean, color = year)) + 
  geom_smooth(size = 1) + 
  geom_smooth( aes(x=timestep, y = mean-sem, color = year), size = 0.3) + 
  geom_smooth( aes(x=timestep, y = mean+sem, color = year), size = 0.3) + 
  geom_ribbon(aes(y = mean, ymin = mean - sem, ymax = mean + sem, fill = year), alpha = .1) +
  #geom_ribbon(aes(y = smooth.mean, ymin = smooth.mean - smooth.sd, ymax = smooth.mean + smooth.sd, fill = year), alpha = .2) +
  xlab("Time step")+ylab("% population infected")+
  theme_light(base_size = 12)+
  xlim(c(0,1000))+ ylim(c(0,1))

#Plot mean % population infected across simulation runs with SDs
ggplot(data = data.all, aes(x=timestep, y = mean, color = year)) + 
  geom_smooth(size = 1) + 
  geom_smooth( aes(x=timestep, y = mean-sd, color = year), size = 0.3) + 
  geom_smooth( aes(x=timestep, y = mean+sd, color = year), size = 0.3) + 
  geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = year), alpha = .1) +
  #geom_ribbon(aes(y = smooth.mean, ymin = smooth.mean - smooth.sd, ymax = smooth.mean + smooth.sd, fill = year), alpha = .2) +
  xlab("Time step")+ylab("% population infected")+
  theme_light(base_size = 12)+
  xlim(c(0,1000))+ ylim(c(0,1))

############################
## Plot networks

year = 4;
weightedEL = V.data.list[[year]]
sim_health_ids = simu_res_ids[[year]][[1]]

timestep=1000

#Need to convert an adjacency matrix to social network graph
adjMat = dils::AdjacencyFromEdgelist(weightedEL)# create adjacency matrix based on edge list.
data = adjMat[["adjacency"]]; rownames(data) = adjMat[["nodelist"]]; colnames(data) = adjMat[["nodelist"]]
m=as.matrix(data) # coerces the data set as a matrix
am.g=graph.adjacency(m,mode= "undirected",weighted=T) # this will create an directed 'igraph object'. Change qualifiers to make "undirected" or unweighted (null)
am.g #igraph object properties 

#Set vertex size
V(am.g)$label.cex <- 0.8

#link up attributes file with network

idx_match = match(sim_health_ids[,sims+1], V(am.g)$name) # find health status
V(am.g)$health = as.factor(sim_health_ids[idx_match,timestep])
#V(am.g)$sex=as.factor(dominance_info$SEX[match(V(am.g)$name,dominance_info$ID)]) #sex attribute
# V(am.g)$agec=as.factor(att$agec[match(V(am.g)$name,attr$id)]) #age attribute
# 
# #set size of nodes by animal age
# V(am.g)$size=V(am.g)$agec*6 #multiplied by 5 because nodes too small otherwise

#set color according to health status
V(am.g)$color=V(am.g)$health #assign the "Health" attribute as the vertex color
V(am.g)$color=gsub("1","blue",V(am.g)$color) #healthy will be blue
V(am.g)$color=gsub("2","red",V(am.g)$color) #sick will be red

# V(am.g)$colorV(am.g)$sex #assign the "Sex" attribute as the vertex color
# V(am.g)$color=gsub("1","plum1",V(am.g)$color) #Females will be orange
# V(am.g)$color=gsub("2","seagreen2",V(am.g)$color) #Males will be lightblue
# V(am.g)$color=gsub("0","white",V(am.g)$color) #unknown sex will be white

#Set color of vertex labels
color_vector=V(am.g)$health #assign the "helath" attribute as the vertex color
color_vector=gsub("1","red",color_vector) #Sick will be red
# color_vector=gsub("2","darkblue",color_vector) #Males will be lightblue
color_vector=gsub("0","white",color_vector) #unknown sex will be white

#set degree attribute
V(am.g)$degree=igraph::degree(am.g)

#set path for saving and plot graph
# tiff(paste(".tiff",sep=""), 
#      units="in", width=10, height=8, res=300, compression = 'lzw')
plot.igraph(am.g, layout=l,  vertex.color=V(am.g)$color,
            #vertex.size=2*V(am.g)$degree+2,
            edge.color="grey20", 
            edge.width=E(am.g)$weight*30,edge.arrow.size = 3,edge.curved=0.5,
            )
# dev.off()



