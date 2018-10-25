##########################################################
#
#   Assignment 2--STAT5703-HEMANT GUPTA-101062246
#
#packages:
#install_package(cluster)
#install_package(stats)
#install_package(fpc)
#install_package(flexclust)
#
#
#NOTE:- Here we have use a function randIndex to check correctness
#       of Cluster between Wine Type and CLuster  
#RandIndex-Compute the (adjusted) Rand, Jaccard and Fowlkes-Mallows 
#index for agreement of two partitions.
##########################################################

#packages:
install.packages("cluster")
install.packages("stats")
install.packages("fpc")
install.packages("flexclust")
install.packages("plyr")
install.packages("fastICA")



## Setting the Path of Directory

drive="C:"
path.upto <- paste("Users","Admin", "Documents","STAT5703-HEMANT-GUPTA-101062246-Assignment2","STAT5703-HEMANT-101062246-Assignment2", sep="/" )
code.dir <- paste(drive, path.upto,"Code", sep="/")
data.dir <- paste(drive, path.upto,"Data", sep="/")
work.dir <- paste(drive, path.upto,"Work", sep="/")
setwd(work.dir)

wine.file <- paste(data.dir,"Wines.dat", sep="/")
wine.col <- paste(data.dir,"Wines.col", sep="/")

##Reading the Table format of Wine Data

wines.dat <- read.table(wine.file, header=FALSE)
headers <- scan(wine.col, "")
wines.dat
headers
names(wines.dat) <- headers
head(wines.dat)

##Setting the Seed


##########################################################
#
#               FUNCTIONS
#
##########################################################

#
##Function: Not Contain in the Set 
#

"%w/o%"<- function(x,y) x[!x %in% y]

#
## Function Set the indices for the training/test sets
#
get.train <- function (data.sz, train.sz)
{
  set.seed(123)
  # Take subsets of data for training/test samples
  # Return the indices
  train.ind <- sample(data.sz, train.sz)
  test.ind <- (data.sz) %w/o% train.ind
  list(train=train.ind, test=test.ind)
}

#
#Function:========DBI Function=======#
#
Davies.Bouldin <- function(A, SS, m) {
  # A  - the centres of the clusters
  # SS - the within sum of squares
  # m  - the sizes of the clusters
  N <- nrow(A)   # number of clusters
  # intercluster distance
  S <- sqrt(SS/m)
  # Get the distances between centres
  M <- as.matrix(dist(A))
  # Get the ratio of intercluster/centre.dist
  R <- matrix(0, N, N)
  for (i in 1:(N-1)) {
    for (j in (i+1):N) {
      R[i,j] <- (S[i] + S[j])/M[i,j]
      R[j,i] <- R[i,j]
    }
  }
  return(mean(apply(R, 1, max)))
}

#
#Function for Data Standardization
#

#################################
# Standardize data
####################################
f.data.std <- function(data) {
  data <- as.matrix(data)
  bar <- apply(data, 2, mean)
  s <- apply(data, 2, sd)
  t((t(data) - bar)/s)
}

#####################################
#WhitenedData: Centre and sphere data
#########################################
Sphere.Data <- function(data) {
  data <- as.matrix(data)
  data <- t(t(data) - apply(data, 2, mean))
  data.svd <- svd(var(data))
  sphere.mat <- t(data.svd$v %*% (t(data.svd$u) * (1/sqrt(data.svd$d))))
  return(data %*% sphere.mat)
}

##############################
#
#FUNCTION: error_cal- For calculating wrong classification in Cluster
#
############################
error_cal <- function(tbl,cluster_size)
{
  wrong_data <- 0
  for(clust in 1:cluster_size)
  {
    wrong_data <- wrong_data + (sum(tbl[,clust])-max(tbl[,clust]))
  }
  return(wrong_data)
}



####################################################
#Function: Clustering_euclidean
# Cluster Size Varies 2 to 15
# SSE and DBI is determined at each iteration
####################################################

clustering_euclidean <- function(data_set,data_set.orig, limit) 
{
 # set.seed(654321)
  oldpar <- par(mfrow = c(4,4))
  par(mar=c(2,1,2,1))
  
  errs <- rep(0, 10)
  DBI <- rep(0, 10)
  perfectness <- rep(0, 10)
  wrong_data <- rep(0,10)
  ##Package(cluster)
  library(cluster)
  library(stats)
  library(fpc)
  library(flexclust)
  
      
  #Loop for different CLuster Size
  for (i in limit)
  {
      min_error <- 179
      min_error_km <- 0
      best.seed <- 0
      #Loop for Seed
      for (j in 2:1000)
      {  
        
         set.seed(j)
         #Clustering Using K means
         KM <- kmeans(data_set[,], i, 25)
         ct.km <- table(data_set.orig$Type, KM$cluster)
         
         #Calculating toal wrong data for each seed         
         error <- error_cal(ct.km,i)
         if(min_error > error)
         {
           #Storing Error count and Kmeans output and best seed for min error
           min_error <- error
           min_error_km <-KM
           best.seed <- j
         }  
         
      }     
      print(paste("Best Seed for Cluster Size " , i ,"is " , best.seed))
      
      print(paste("Total Wrong in Cluster Size " , i ,"is " , min_error))
      
      print(paste("Centroids for Cluster Size " , i ,"are :"))
      
      print(min_error_km$centers)
      
      #Distribution of Data in each Cluster
      ct.km <- table(data_set.orig$Type, min_error_km$cluster)
      cat("\nDistribution of Wine types:\n")
      rownames(ct.km) <- c("WineType_1  ", "WineType_2  ", "WineType_3  ")
      
      print(ct.km)

      #Plotting the CLuster
      plotcluster(data_set, col=data_set.orig$Type, min_error_km$cluster, main=paste(i,"clusters-Eucli"))
      
      if(length(limit) > 1)
      {  
         #CLuster Analysis
         errs[i-1] <- sum(min_error_km$withinss)
         DBI[i-1] <- Davies.Bouldin(min_error_km$centers, min_error_km$withinss, min_error_km$size)
      
         wrong_data[i-1] <- min_error
      }
          
  }
   if(length(limit) > 1)
   {
       plot(2:15, errs, main = "SSE")
       lines(2:15, errs)
       #
       plot(2:15, DBI, main = "Davies-Bouldin")
       lines(2:15, DBI)
       #
   }
   else
   {  
     print(paste("Optimal Cluster: ",limit," Analysis:"))
     
     wrong_data <- rep(0,limit)
     for(clust in 1:limit)
     {
       print(paste("Cluster: ",clust,"-"))
       print(paste("Classified Data in- ",rownames(ct.km)[which.max(ct.km[,clust])]," : ",max(ct.km[,clust]) ))
       wrong_data[clust] <- (sum(ct.km[,clust])-max(ct.km[,clust]))
       print(paste("Misclassified Data in : ",rownames(ct.km)[which.min(ct.km[,clust])]," : ",min(ct.km[,clust])))
       y<- sum(ct.km[,clust])- max(ct.km[,clust])- min(ct.km[,clust])
       print(paste("Misclassified Data of : ",rownames(ct.km)[which(ct.km[,clust] == y)]," : ",y))
       
       print("")
        
     }
     return(ct.km)
   }
   return(wrong_data)
}


###############################################################
#Function: Clustering_manhattan
# Cluster Size Varies 2 to 15
# silhouette-optimal is determined at each iteration
################################################################

clustering_manhattan <- function(data_set, data_set.orig,limit)
{ 
  asw <- numeric(15)
  oldpar <- par(mfrow = c(4,4))
  par(mar=c(2,1,2,1))
  perfectness <- rep(0, 10)
  wrong_data <- rep(0,10)
  ##Package(cluster)
  library(cluster)
  library(stats)
  library(fpc)
  library(flexclust)
  
    #Using PAM for Clustering with Manhattan distance
  for (k in limit)
  { 
    min_error <- 179
    min_error_km <- 0
    best.seed <- 0
    
    for(j in 1:1000 )
    {
       set.seed(j)
       KM <- pam(data_set, k, metric = "manhattan")
       ct.km <- table(data_set.orig$Type, KM$clustering)
       error <- error_cal(ct.km, k)
       if(min_error > error)
       {
         min_error <- error
         min_error_km <-KM
         best.seed <- j
       }  
    }
    
    #Cluster Plot
    plotcluster(data_set, col=data_set.orig$Type,min_error_km$clustering, main=paste(k,"clusters-Manh"))
    asw[k]<- min_error_km $ silinfo $ avg.width
    wrong_data[k-1] <- min_error
    print(paste("Best Seed for Cluster Size " , k ,"is " , best.seed))
    
    print(paste("Total Wrong in Cluster Size " , k ,"is " , min_error))
    
    print(paste("Medoids for Cluster Size " , k ,"are :"))
    
    print(min_error_km$medoids)
    
    #Distribution of Data in each Cluster
    ct.km <- table(data_set.orig$Type, min_error_km$clustering)
    cat("\nDistribution of Wine types:\n")
    rownames(ct.km) <- c("WineType_1  ", "WineType_2  ", "WineType_3  ")
    
    print(ct.km)
    
    perfectness[k-1] <- randIndex(ct.km)
  }
  if(length(limit) > 1)
  {
      k.best <- which.max(asw)
      cat("silhouette-optimal number of clusters:", k.best, "\n")
      plot(1:15, asw, type= "h", main = "pam() clustering assessment",
       xlab= "k  (# clusters)", ylab = "average silhouette width")
#      axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
      #
      plot(2:15, perfectness, main = "Correctness")
      lines(2:15, perfectness)
      #
  }
  else
  {
    
    print(paste("Optimal Cluster: ",limit," Analysis:"))
    
    wrong_data <- rep(0,limit)
    for(clust in 1:limit)
    {
      print(paste("Cluster: ",clust,"-"))
      print(paste("Classified Data in- ",rownames(ct.km)[which.max(ct.km[,clust])]," : ",max(ct.km[,clust]) ))
      wrong_data[clust] <- (sum(ct.km[,clust])-max(ct.km[,clust]))
      print(paste("Misclassified Data in : ",rownames(ct.km)[which.min(ct.km[,clust])]," : ",min(ct.km[,clust])))
      y<- sum(ct.km[,clust])- max(ct.km[,clust])- min(ct.km[,clust])
      print(paste("Misclassified Data of : ",rownames(ct.km)[which(ct.km[,clust] == y)]," : ",y))
      
      print("")
      
    }
    return(ct.km)
  }
  return(wrong_data)
}



###########################################################

########Question 1#########################################

############################################################

##Selecting index based on Wine$Type

Type1_index = which(wines.dat$Type == 1)
Type2_index = which(wines.dat$Type == 2)
Type3_index = which(wines.dat$Type == 3)


#
##Distributing Training Set and Test Set for each wine Type
#
Type1.train.sz <- round((2*length(Type1_index))/3) # Set the size of the training sample
# Get the indices for the training and test samples
(Type1.ind <- get.train(Type1_index, Type1.train.sz ))

Type1.ind$train
Type1.ind$test

Type2.train.sz <- round((2*length(Type2_index))/3) # Set the size of the training sample
# Get the indices for the training and test samples
(Type2.ind <- get.train(Type2_index, Type2.train.sz ))

Type2.ind$train
Type2.ind$test

Type3.train.sz <- round((2*length(Type3_index))/3) # Set the size of the training sample
# Get the indices for the training and test samples
(Type3.ind <- get.train(Type3_index, Type3.train.sz ))

Type3.ind$train
Type3.ind$test

wine.ind = list(train=c(Type1.ind$train,Type2.ind$train,Type3.ind$train), test=c(Type1.ind$test,Type2.ind$test,Type3.ind$test))

##Getting the wines Training and Test Set

train.wine <- wines.dat[wine.ind$train,]
test.wine <- wines.dat[wine.ind$test,]
library(plyr)

Class1.train.wine.ind<- sum(train.wine$Type == 1)
Class2.train.wine.ind<- sum(train.wine$Type == 2)
Class3.train.wine.ind<- sum(train.wine$Type == 3)
Class1.test.wine.ind<- sum(test.wine$Type == 1)
Class2.test.wine.ind<- sum(test.wine$Type == 2)
Class3.test.wine.ind<- sum(test.wine$Type == 3)

Total_Class1<- sum(wines.dat$Type == 1)
Total_Class2<- sum(wines.dat$Type == 2)
Total_Class3<- sum(wines.dat$Type == 3)

out = c (Total_Class1, Class1.test.wine.ind, Class1.train.wine.ind, Total_Class2, Class2.test.wine.ind, Class2.train.wine.ind, Total_Class3, Class3.test.wine.ind,Class3.train.wine.ind)

##Setting X axis name
x.names=c("T1","T1_Tst","T1_Trn","T2","T2_Tst","T2_Trn","T3","T3_Tst","T3_Trn")

barplot(out, main="Data Spliting",xaxt="n")
axis(1,at = 1:9,labels=x.names)




###########################################################

########Question 2: Clustering : Euclidean Distance########

############################################################

print("QUESTION: 2- CLUSTERING USING EUCLIDEAN DISTANCE")

######  RAW  DATA####################
print("WINE DATASET: RAW DATA")

head(train.wine)
summary(train.wine)

head(test.wine)
summary(test.wine)


##Data Standardization for Training and Test Set
print("WINE DATASET: STANDARDIZE DATA")

train.wine.std <- f.data.std(train.wine[-1])
test.wine.std <- f.data.std(test.wine[-1])

head(train.wine.std)
summary(train.wine.std)

head(test.wine.std)
summary(test.wine.std)


##Whitening Data for Training and Test Set
print("WINE DATASET: WHITENED DATA")

test.wine.white <- Sphere.Data(test.wine[-1])
colnames(test.wine.white) <- headers[-1]

apply(test.wine.white, 2, mean)
apply(test.wine.white, 2, sd)

train.wine.white <- Sphere.Data(train.wine[-1])
colnames(train.wine.white) <- headers[-1]

apply(train.wine.white, 2, mean)
apply(train.wine.white, 2, sd)


head(train.wine.white)
summary(train.wine.white)

head(test.wine.white)
summary(test.wine.white)

cluster_range <- 2:15


##Clustering on train and test Raw DataSet 
wrong.train.eucli.all <- clustering_euclidean(train.wine,train.wine, cluster_range)

##Note:- From Davis Bouldin it is Cluster Size 5 is the best solution

cluster_value <- 5
print(paste("Raw Training Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.train.eucli <- clustering_euclidean(train.wine,train.wine, cluster_value)

print(paste("Raw Test Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.test.eucli <- clustering_euclidean(test.wine,test.wine,cluster_value)


##Clustering on train and test Standardize DataSet

wrong.train.std.eucli.all <- clustering_euclidean(train.wine.std,train.wine, cluster_range)
##Note:- From Davis Bouldin it is Cluster Size 3 is the best solution

cluster_value <- 3
print(paste("Standardize Training Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.train.std.eucli <- clustering_euclidean(train.wine.std,train.wine, cluster_value)

print(paste("Standardize Test Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.test.std.eucli <- clustering_euclidean(test.wine.std,test.wine,cluster_value)

##Clustering on train and test Whitened DataSet

wrong.train.white.eucli.all <- clustering_euclidean(train.wine.white,train.wine, cluster_range)
##Note:- From Davis Bouldin it is Cluster Size 7 is the best solution

cluster_value <- 7
print(paste("Whitened Training Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.train.white.eucli <- clustering_euclidean(train.wine.white,train.wine, cluster_value)

print(paste("Whitened Test Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.test.white.eucli <- clustering_euclidean(test.wine.white,test.wine,cluster_value)


print("Wrong Data Analsyis of Training Dataset Raw, Standardize and Whitened with cluster range 2 to 15")
print("Raw Train Dataset -")
print(wrong.train.eucli.all)

print("Standardize Train Dataset-")
print(wrong.train.std.eucli.all)

print("Whitened Train Dataset-")
print(wrong.train.white.eucli.all)


print("Table Analsyis of Training Dataset Raw, Standardize and Whitened with Davies Bouldin Values")
print("Raw Train and Test Dataset -")
print(tbl.train.eucli)
print(tbl.test.eucli)

print("Standardize Train and Test Dataset-")
print(tbl.train.std.eucli)
print(tbl.test.std.eucli)

print("Whitened Train and Test Dataset-")
print(tbl.train.white.eucli)
print(tbl.train.white.eucli)

###########################################################

####Question 3: PCA WITH Clustering : Euclidean Distance####

############################################################

print("QUESTION: 3- PCA WITH CLUSTERING USING EUCLIDEAN DISTANCE")

wine.std <- f.data.std(wines.dat)
# Get principal component vectors using prcomp 
pc.wine <- prcomp(wine.std)
summary(pc.wine)
plot(pc.wine)
# First  principal components
wine.pc <- data.frame(pc.wine$x[,1:3])
head(wine.pc)


# Get principal component vectors using prcomp 
pc.train.std <- prcomp(train.wine.std)
summary(pc.train.std)
plot(pc.train.std)

# First  principal components
train.wine.std.pc <- data.frame(pc.train.std$x[,1:3])
head(train.wine.std.pc)


# Get principal component vectors using prcomp 
pc.test.std <- prcomp(test.wine.std)
summary(pc.test.std)
plot(pc.test.std)
# First  principal components
test.wine.std.pc <- data.frame(pc.test.std$x[,1:3])
head(test.wine.std.pc)

cluster_range <- 2:15

##Clustering on Actual Raw DataSet  

wrong.data.wine <- clustering_euclidean(wines.dat,wines.dat,cluster_range)

##Note:- From Davis Bouldin it is Cluster Size 6 is the best solution

cluster_value <- 6
print(paste("Wine Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.wine.eucli <- clustering_euclidean(wines.dat,wines.dat, cluster_value)

##Clustering on Actual Raw DataSet after PCA 

wrong.data.wine.pc <- clustering_euclidean(wine.pc,wines.dat,cluster_range)

##Note:- From Davis Bouldin it is Cluster Size 3 is the best solution

cluster_value <- 3
print(paste("PCA Wine Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.wine.pc.eucli <- clustering_euclidean(wine.pc,wines.dat,cluster_value)



##Clustering on train and test Standardize DataSet after PCA

wrong.train.std.pc.all <- clustering_euclidean(train.wine.std.pc,train.wine, cluster_range)
##Note:- From Davis Bouldin it is Cluster Size 3 is the best solution

cluster_value <- 3
print(paste("Standardize Training Dataset after PCA-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.train.std.pc <- clustering_euclidean(train.wine.std.pc,train.wine, cluster_value)

print(paste("Standardize Test Dataset After PCA-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.test.std.pc <- clustering_euclidean(test.wine.std.pc,test.wine,cluster_value)


###########################################################

####Question 4: Clustering : Manhattan Distance####

############################################################


print("QUESTION: 4- CLUSTERING USING MANHATTAN DISTANCE")

cluster_range <- 2:15


##Clustering on train and test Raw DataSet 
wrong.train.manh.all <- clustering_manhattan(train.wine,train.wine, cluster_range)

##Note:- From silhouette-optimal number of clusters is 2 is the best solution

cluster_value <- 2
print(paste("Raw Training Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.train.manh <- clustering_manhattan(train.wine,train.wine, cluster_value)

print(paste("Raw Test Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.test.manh <- clustering_manhattan(test.wine,test.wine,cluster_value)


##Clustering on train and test Standardize DataSet

wrong.train.std.manh.all <- clustering_manhattan(train.wine.std,train.wine, cluster_range)
##Note:- From silhouette-optimal number of clusters is 3 is the best solution

cluster_value <- 3
print(paste("Standardize Training Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.train.std.manh <- clustering_manhattan(train.wine.std,train.wine, cluster_value)

print(paste("Standardize Test Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.test.std.manh <- clustering_manhattan(test.wine.std,test.wine,cluster_value)

##Clustering on train and test Whitened DataSet

wrong.train.white.manh.all <- clustering_manhattan(train.wine.white,train.wine, cluster_range)
##Note:- From silhouette-optimal number of clusters 3 is the best solution

cluster_value <- 3
print(paste("Whitened Training Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.train.white.manh <- clustering_manhattan(train.wine.white,train.wine, cluster_value)

print(paste("Whitened Test Dataset-From silhouette-optimal number of clusters is : ",cluster_value))

tbl.test.white.manh <- clustering_manhattan(test.wine.white,test.wine,cluster_value)


print("Wrong Data Analsyis of Training Dataset Raw, Standardize and Whitened with cluster range 2 to 15")
print("Raw Train Dataset -")
print(wrong.train.manh.all)

print("Standardize Train Dataset-")
print(wrong.train.std.manh.all)

print("Whitened Train Dataset-")
print(wrong.train.white.manh.all)


print("Table Analsyis of Training Dataset Raw, Standardize and Whitened with Davies Bouldin Values")
print("Raw Train and Test Dataset -")
print(tbl.train.manh)
print(tbl.test.manh)

print("Standardize Train and Test Dataset-")
print(tbl.train.std.manh)
print(tbl.test.std.manh)

print("Whitened Train and Test Dataset-")
print(tbl.train.white.manh)
print(tbl.train.white.manh)



###########################################################

####Question 5: ICA WITH Clustering : Euclidean Distance####

############################################################

print("QUESTION: 5- ICA WITH CLUSTERING USING EUCLIDEAN DISTANCE")

#Whitening the whole wine dataset
wine.white <- Sphere.Data(wines.dat[-1])

##Taking number of components as 3 as want to compare it with result of PCA 3 components

library(fastICA)
wine.ica <- fastICA(wine.white[,-1], 3, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200, tol = 0.0001, verbose =
             TRUE)

#Estimated Source Matrix
head(wine.ica$S)
cluster_range <- 2:15
##Clustering on Actual Raw DataSet after PCA 

wrong.data.wine.ica <- clustering_euclidean(wine.ica$S,wines.dat,cluster_range)


##Note:- From Davis Bouldin it is Cluster Size 7 is the best solution

cluster_value <- 7
print(paste("ICA Wine Dataset-From Davis Bouldin optimal Cluster Size is : ",cluster_value))

tbl.wine.ica.eucli <- clustering_euclidean(wine.ica$S, wines.dat,cluster_value)


print("Wrong Data Analsyis of Raw wine Dataset after ICA with cluster range 2 to 15")
print("Raw Wine Dataset After ICA -")
print(wrong.data.wine.ica)


print("Table Analsyis of Raw wine dataset with Davies Bouldin Value")
print("Raw Wine Dataset after ICA and Cluster value 7 -")
print(tbl.wine.ica.eucli)
