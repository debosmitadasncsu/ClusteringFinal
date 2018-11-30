library(splines)
library(readr)
library(data.table)
library(factoextra)
library(dplyr)
library(tidyverse)
library(mclust)
library(reshape)

load("C:/Users/mdsau/OneDrive/Desktop/IAA/AA 502/Clustering/Data/final_data.Rdata")

times <- seq(1,295)/100 # Observations in 1/100th of a second
X <- bs(times,intercept=TRUE,df=60) #create a spline to 

#model the data
betas <- matrix(0,ncol=60,nrow = 6792)

###########################################################
# run a linear regression on each data set
# here I am manipulating my data so I can cluster
###########################################################

for (ii in 1:6792){
  temp <- lm(as.numeric(final_data[ii,6:300])~X-1) #-1 removes the natural intercept
  betas[ii,]  <- coefficients(temp)
}
cdata <- cbind(final_data[,1:5],betas)

#CONVERT EVERTYING TO 'numbers'
cdata$AGE <- as.numeric(cdata$AGE)
cdata$EVER_SMOKE <- as.numeric(cdata$EVER_SMOKE)
cdata$ASTHMA <- as.numeric(cdata$ASTHMA)
cdata$POVERTY_RATIO <- as.numeric(cdata$POVERTY_RATIO)

###########################################################
# K-MEANS CLUSTERING
###########################################################

pca <- princomp(cdata[,2:65])
pca$sdev[1:5]

# Plot this 
fviz_nbclust(pca$scores, kmeans, method = "wss",k.max=10) +
geom_vline(xintercept = 4 , linetype = 2, colour = "steelblue")
fviz_nbclust(pca$scores, kmeans, method = "silhouette",k.max=10)

set.seed(12345)
kmean_4 <- kmeans(pca$scores,4)

cdata$clust <- kmean_4$cluster

# Break up our clusters for exploration
clusters <- list()
for(ii in 1:4){
  clusters[[ii]] <-  cdata %>% 
    filter(clust == ii)
}

# See how they differ from population mean
x <- cbind(colMeans(cdata))
y <- x
for (ii in 1:4) {
  x <- cbind(x,colMeans(clusters[[ii]])-y)
}

# Get the mean on original scale of each cluster
cdata_agg <- cdata %>% 
  group_by(clust) %>% 
  summarise_all(mean)

# Get only spirometry data
cdata_agg_plot <- cdata_agg[,-c(2:6)]

# Create dataset to plot with time
fdata <- gather(cdata_agg_plot, time, x, c(2:61))
fdata$clust <- as.factor(fdata$clust)
fdata$time <- as.integer(fdata$time)

# Plot all clusters
ggplot(data=fdata) +
  geom_line(aes(x=time,y=x, colour = clust)) +
  xlab("Time") + 
  ylab("Mean Spirometry") +
  labs(color = "Clusters", title = "Lung Capacity") +
  theme(plot.title = element_text(hjust = 0.5))

# Plot only cluster 3
ggplot(subset(fdata, clust == 3)) +
  geom_line(aes(x=time,y=x), color = "steelblue") +
  xlab("Time") + 
  ylab("Mean Spirometry") +
  labs(title = "Lung Capacity - Cluster 3") +
  theme(plot.title = element_text(hjust = 0.5))
  
# Plot cluster 1,2,4
ggplot(subset(fdata, clust %in% c(1,2,4))) +
  geom_line(aes(x=time,y=x, colour = clust)) +
  xlab("Time") + 
  ylab("Mean Spirometry") +
  labs(color = "Clusters", title = "Lung Capacity") +
  theme(plot.title = element_text(hjust = 0.5))

###########################################################
# MCLUST
###########################################################

# Determine optimal number of clusters
set.seed(12345)
clustBIC <- mclustBIC(cdata[,10:20], modelNames='VVV', G = 1:20)
plot(clustBIC)
abline(v=6, col="red",lty=2)

# Create clusters
set.seed(12345)
clust <- Mclust(cdata[,10:20],G=6,modelNames='VVV')
cdata$class <- as.factor(clust$classification)

# Get the mean on original scale of each cluster
cdataBIC_agg <- cdata %>% 
  group_by(class) %>% 
  summarise_all(mean)

# Get only spirometry data
cdataBIC_agg_plot <- cdataBIC_agg[,-c(2:6)]

# Create dataset to plot with time
fdataBIC <- gather(cdataBIC_agg_plot, time, x, c(2:61))
fdataBIC$clust <- as.factor(fdataBIC$class)
fdataBIC$time <- as.integer(fdataBIC$time)

# Plot all clusterss
ggplot(data=fdataBIC) +
  geom_line(aes(x=time,y=x, colour = class)) +
  xlab("Time") + 
  ylab("Mean Spirometry") +
  labs(color = "Clusters", title = "Lung Capacity") +
  theme(plot.title = element_text(hjust = 0.5))

# Compare cluster 4 and cluster 3 from kmeans
compare <- rbind(cdata_agg[3,2:66],cdataBIC_agg[4,2:66])
compare <- data.frame(ID = c('kmeans3','mclust4'), compare)

compare_agg_plot <- compare[,-c(2:6)]
fdata_compare <- melt(compare_agg_plot, id= "ID")
fdata_compare <- arrange(fdata_compare, ID)
fdata_compare$time <- seq(1,60, by= 1)

ggplot(data=fdata_compare) +
  geom_line(aes(x=time,y=value, colour = ID)) +
  xlab("Time") + 
  ylab("Mean Spirometry") +
  labs(color = "Procedure", title = "Lung Capacity - Comparison") +
  theme(plot.title = element_text(hjust = 0.5))
