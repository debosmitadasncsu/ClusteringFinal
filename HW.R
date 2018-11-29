library(splines)
library(factoextra)
library(mclust)
library(ggplot2)

times <- seq(1,295)/100

plot(times,X%*%beta1,ylim=c(0,100),type='l')

times <- seq(1,295)/100 # Observations in 1/100th of a second
X <- bs(times,intercept=TRUE,df=60) #create a spline to 
#model the data
betas <- matrix(0,ncol=60,nrow = 6792)
###########################################################
# run a linear regression on each data set
# here I am manipulating my data you I can cluster
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

# 1
pcaOut = princomp(cdata[,2:65])
plot(pcaOut)

means_pcaOut <- colMeans(pcaOut$scores)
pcaOut$center

# so we are scaling the scores but not modifying their meen ->
# important to know when we use it to plot
###################################################################
#compare the SD between the princ component analysis
sd_pcaOut <- apply(pcaOut$scores,2,sd)
sd_pcaOut
pcaOut$sdev
#standard deviations
#45.2150785 31.2923808 22.3770254 17.3312303 13.0480844

# 2
set.seed(12345)
fviz_nbclust(scale(pcaOut$scores), kmeans, method = "wss") #optimal number = 4 clusters
fviz_nbclust(scale(pcaOut$scores), kmeans, method = "silhouette") #optimal number = 2 clusters

set.seed(12345)
k_means4 <- kmeans(scale(pcaOut$scores),4,nstart = 25)
cdata$clust <- k_means4$cluster
###################################################################
# Try to understand the clusters
###################################################################
#
clusters <- list()
for( ii in 1:4){
  clusters[[ii]] <-  cdata %>% filter(clust == ii)
}
#
####################################################################

# Find the means of each cluster to "Name them"
x <- cbind(colMeans(cdata))
y <- x
for (ii in 1:4) {
  print("Initial")
  print(x)
  print("-------------------")
  x <- cbind(x,colMeans(clusters[[ii]])-y)
  print("Last")
  print(x)
  print("-------------------")
}

######################################################################
# Plot some of the "centers"
######################################################################
plot(k_means4$centers[,1],-1*k_means4$centers[,2],ylab="Distance from Work",xlab = "Lung Capacity",
     axes=F,xlim=c(-2.5,2.5),ylim=c(-1,1),pch = 16,col="Light Blue",cex=2)

abline(v=0,lty=2)
abline(h=0,lty=2)
# Why did I multiply the second one by -1? 

######################################################################
# Plot some of the "centers"
######################################################################
plot(k_means4$centers[,2],k_means4$centers[,4],ylab="Tall but close absentee",xlab = "Heavy",
     axes=F,xlim=c(-2.5,2.5),ylim=c(-1,1),pch = 16,col="Light Blue",cex=2)

abline(v=0,lty=2)
abline(h=0,lty=2)

######################################################################
# Plot some of the "centers"
######################################################################
plot(k_means4$centers[,1],k_means4$centers[,4],ylab="Tall but close absentee",xlab = "Heavy",
     axes=F,xlim=c(-2.5,2.5),ylim=c(-1,1),pch = 16,col="Light Blue",cex=2)

abline(v=0,lty=2)
abline(h=0,lty=2)


######################################################################
# Plot some of the "centers"
######################################################################
plot(k_means4$centers[,1],k_means4$centers[,3],ylab="Tall but close absentee",xlab = "Heavy",
     axes=F,xlim=c(-2.5,2.5),ylim=c(-1,1),pch = 16,col="Light Blue",cex=2)

abline(v=0,lty=2)
abline(h=0,lty=2)





######################################################################
# Plot some of the "centers"
######################################################################
plot(k_means4$centers[,3],k_means4$centers[,2],ylab="Tall but close absentee",xlab = "Heavy",
     axes=F,xlim=c(-2.5,2.5),ylim=c(-1,1),pch = 16,col="Light Blue",cex=2)

abline(v=0,lty=2)
abline(h=0,lty=2)


cmeans <- matrix(colMeans(comp),60,1)
stdev  <- matrix(apply(comp,2,sd),60,1)

cl1 <- matrix(k_means4$centers[1,],60,1)
bl1 <- cl1 * stdev + cmeans
plot(bl1)
sfun <- splinefun(times,X%*%bl1)   #this creates an interpolant of the curve from min(times) to max(times)
integrate(sfun,min(times),max(times))

plot(times,X%*%bl1,ylab="ML",ylim=c(0,60),type='l')

cl2 <- matrix(k_means4$centers[2,],60,1)
bl2 <- cl2 * stdev + cmeans
sfun <- splinefun(times,X%*%bl2)   #this creates an interpolant of the curve from min(times) to max(times)
integrate(sfun,min(times),max(times))
# -0.7687984 with absolute error < 0.00011
lines(times,X%*%bl2,lwd=2,col=2)

cl3 <- matrix(k_means4$centers[3,],60,1)
bl3 <- cl3 * stdev + cmeans
# 0.7157545 with absolute error < 0.00012
sfun <- splinefun(times,X%*%bl3)   #this creates an interpolant of the curve from min(times) to max(times)
integrate(sfun,min(times),max(times))
# 0.4281903 with absolute error < 0.00011
lines(times,X%*%bl3,lwd=2,col=3)

cl4 <- matrix(k_means4$centers[4,],60,1)
bl4 <- cl4 * stdev + cmeans
sfun <- splinefun(times,X%*%bl4)   #this creates an interpolant of the curve from min(times) to max(times)
integrate(sfun,min(times),max(times))
# 0.8251692 with absolute error < 0.00012
lines(times,X%*%bl4,lwd=2,col=4)



# Analyzing 3rd cluster
# Cluster sizes
sort(table(k_means4$clust))
clust <- names(sort(table(k_means4$clust)))

clus3_data <- cdata[k_means4$clust==clust[3],]
#mean of each column
colMeans(clus3_data[,2:65])

#comparison with other clusters

clus1_data <- cdata[k_means4$clust==clust[1],]
#mean of each column
colMeans(clus1_data[,2:65])

clus2_data <- cdata[k_means4$clust==clust[2],]
#mean of each column
colMeans(clus2_data[,2:65]) 

clus4_data <- cdata[k_means4$clust==clust[4],]
#mean of each column
colMeans(clus4_data[,2:65]) # most lung capacity 0.40042750 - less asthma = more lung capacity??

############################### descriptive stats ###########################

#            AGE          EVER_SMOKE          ASTHMA               POVERTY_RATIO
#clus1    30.4331551      0.4502674           0.4502674            2.7423936
#clus2    29.29358491     0.53207547          0.53207547           1.96581283    #least lung cap
#clus3    29.6052174      0.4915942           0.4915942            2.4582435
#clus4    30.74706092     0.40042750          0.40042750           1.98698575    #largest lung cap
class <- cdata$ASTHMA

mm_bic <- mclustBIC(cdata[,10:20])
plot(mm_bic)

#cluster based upon 20 model components and VVV model
clust    <- MclustDA(cdata[,10:20],class,G=1:20, modelNames='VVV') # G is the number of components
plot(clust)


# prof method
#cluster based upon 2 model components
clust    <- Mclust(cdata[,10:20],G=1:20, modelNames='VVV') # G is the number of components
cdata$class <- as.factor(clust$classification)
ggplot(data=df,mapping=aes(x=x,y=y,color=cdata$ASTHMA,shape=class)) + geom_point()
ggplot(data=cdata,mapping=aes(x=cdata[,10],y=cdata[,1],color=id,shape=class)) + geom_point()

