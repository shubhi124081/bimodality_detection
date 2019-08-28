
######################function for sampling wih error 

sampling_error_fun <- function(samp_error, psu = psu, mu_m1 = mu_m1, 
                               mu_m2 = mu_m2, sig_m1= sig_m1, sig_m2 = sig_m2){

  
samp_error05 <- rnorm(1000, 0, samp_error)
set.seed(50)
n <- 1000
sim_psu1_e05 <- matrix(NA, ncol= length(psu), nrow = 1000)

######################sampling from distribution 1 
for(i in 1: length(psu)){
  temp1 <- rnorm(n, mu_m1[i], sig_m1[i])  #first distribution
  
  for (j in 1:1000){
    sim_psu1_e05[j,i] <- temp1[j] + samp_error05[j]
  }
}
colnames(sim_psu1_e05) <- psu

######################sampling from distribution 2
set.seed(50)
n <- 1000
sim_psu2_e05 <- matrix(NA, ncol= length(psu), nrow = 1000)

for(i in 1: length(psu)){
  temp2 <- rnorm(n, mu_m2[i], sig_m2[i]) #second distribution 
  
  for(j in 1:1000){
    sim_psu2_e05[j,i] <- temp2[j] + samp_error05[j]
  }                      
}
colnames(sim_psu2_e05) <- psu
# create index column
index <- seq(1, nrow(sim_psu1_e05), by=1)
sim_psu1_df_e05 <- data.frame(sim_psu1_e05, Index= index)
sim_psu2_df_e05 <- data.frame(sim_psu2_e05, Index= index)
#rename 
colnames(sim_psu1_df_e05) <- c(psu, "Index")
colnames(sim_psu2_df_e05) <- c(psu, "Index")

sim_return <- list()
sim_return[[1]] <- sim_psu1_df_e05
sim_return[[2]] <- sim_psu2_df_e05
return(sim_return)
}


######################function for generating numbers and clustering
clust_fun <- function(n, dist1, dist2){
  index <- seq(1,1000, by=1)
  rm <- sample(1:1000, n, replace = TRUE)
  require(mclust)
  # tmpdist1 <- sim_psu1_df[,c(col]
  dist1 <- dist1[index %in% rm]
  dist2 <- dist2[index %in% rm]
  dist <- c(dist1, dist2)
  # return(dist)
  clust_05 <- mclustBIC(dist)
  clust <- clust_05
  return(clust)
}

clust_fun_01 <- function(n, dist){
  index <- seq(1,500, by=1)
  rm <- sample(1:500, n, replace = TRUE)
  require(mclust)
  dist_smp <- dist[index %in% rm]
  clust <- mclustBIC(dist_smp, verbose = F)
  return(clust)
}

######################proportion function

proportion.fun <- function(nsamples, d1, niter){
clusterBIC <- matrix(NA, ncol = 3, nrow = niter)

require(mclust)

for(i in 1:niter){
  tmp <- clust_fun_01(n=nsamples, dist = d1)
  # print(tmp[1])
  clusterBIC[i, 1] <- tmp[1]
  clusterBIC[i, 2] <- tmp[2]
  clusterBIC[i, 3] <- tmp[2,2]
}
colnames(clusterBIC) <- c("E1", "E2", "V2")
count <- 0
for(i in 1:niter){
  if(abs(clusterBIC[i,1]) > abs(clusterBIC[i,2])){
    count <- count+1
  }
}

#proportion 
p <- count/ niter
return(p)
}

######################proportion function for unequal variances 


varUneq <- function(nsamples, d1){
  clusterBIC <- matrix(NA, ncol = 3, nrow = 500)
  
  require(mclust)
  
  for(i in 1:500){
    tmp <- clust_fun_01(n=nsamples, dist = d1)
    # print(tmp[1])
    clusterBIC[i, 1] <- tmp[1]
    clusterBIC[i, 2] <- tmp[2]
    clusterBIC[i, 3] <- tmp[2,2]
  }
  colnames(clusterBIC) <- c("E1", "E2", "V2")
  count <- 0
  for(i in 1:500){
    if(abs(clusterBIC[i,1]) > abs(clusterBIC[i,2]) || abs(clusterBIC[i, 1]) > abs(clusterBIC[i, 3])){
      count <- count+1
    }
  }
  
  #proportion 
  p <- count/ 500
  return(p)
}

######################proportion function for unequal weights 

wp.fun <- function(nsamples, d1){
  clusterBIC <- matrix(NA, ncol = 3, nrow = 500)
  
  require(mclust)
  
  for(i in 1:500){
    tmp <- clust_fun_01(n=nsamples, dist = d1)
    # print(tmp[1])
    clusterBIC[i, 1] <- tmp[1]
    clusterBIC[i, 2] <- tmp[2]
    clusterBIC[i, 3] <- tmp[2,2]
  }
  colnames(clusterBIC) <- c("E1", "E2", "V2")
  count <- 0
  for(i in 1:500){
    if(is.na(clusterBIC[i, 2]) == FALSE){
    if(abs(clusterBIC[i,1]) > abs(clusterBIC[i,2])){
      count <- count+1
    }
    }
  }
  
  #proportion 
  p <- count/ 500
  return(p)
}
######################generate random normal 

generate_rm<- function(smp_error, mean_vector1,  mean_vector2,
                       sd_vector1, weight1, weight2,
                       sd_vector2, colnm, n_samples ){
  
  sim_psu_e<- matrix(NA, ncol= length(colnm), nrow = n_samples)
  error <- rnorm(n_samples, 0, smp_error)
  for(i in 1:length(colnm)){
    
    means <- c(mean_vector1[i], mean_vector2[i])
    sds <- c(sd_vector1[i], sd_vector2[i])
    ind <- sample(1:2, n_samples, replace = TRUE, prob= c(weight1[i], weight2[i]))
    x <- rnorm(n_samples, mean = means[ind], sd = sds[ind])
    # plot(density(x))
    
    for(j in 1:n_samples){
      sim_psu_e[j,i] <- x[j] + error[j]
    }                      
  }
  return(sim_psu_e)
}


######################expected value from clusters function
mean.approx.fun <- function(nsamples, d1, n.iter){
 mean <- matrix(NA, ncol = 2, nrow = n.iter)
index <- seq(1,n.iter, by=1)

  require(mclust)
  for(i in 1:n.iter){
    rm <- sample(1:n.iter, nsamples, replace = TRUE)
    dist_smp <- d1[index %in% rm]
    mod1 <- Mclust(dist_smp, verbose = FALSE)
    mean[i,1] <- mod1$parameters$mean[1]
    mean[i, 2]<- mod1$parameters$mean[2]
    # r.m1 <- mean(mean[,1], na.rm = T)
    # r.m2 <- mean(mean[,2], na.rm = T)
  }
  return(mean)
  }
  
######################expected variance from clusters function
var.approx.fun <- function(nsamples, d1, n.iter){
  var <- matrix(NA, ncol = 2, nrow = n.iter)
  index <- seq(1,n.iter, by=1)
  
  require(mclust)
  for(i in 1:n.iter){
    rm <- sample(1:n.iter, nsamples, replace = TRUE)
    dist_smp <- d1[index %in% rm]
    mod1 <- Mclust(dist_smp, verbose = FALSE)
    var[i,1] <- mod1$parameters$variance$sigmasq[1]
    var[i, 2]<- mod1$parameters$variance$sigmasq[2]
    
  }
  return(var)
}


#######################non-iid sampling function 
# generate.niid <- function(nsamples, dist, stand.dev){
# rm_matrix <- matrix(NA, ncol = 3, nrow = nsamples)
# rm <- sample(1:1000, (nsamples*3), replace = TRUE)
# for(i in 1:nsamples){
#   for(j in 1 :3){
#     rm_matrix[i,j] <- abs(round(rnorm(1, rm[i], sd= stand.dev), digits = 0))
#   }
# }
# rm_v <- c(rm_matrix[,1], rm_matrix[,2], rm_matrix[,3])
# index <- seq(1,1000, by=1)
# smp_niid <- dist[index %in% rm_v]
# return(smp_niid)
# }

generate.niid <- function(nsamples, dist, stand.dev){
  rm_matrix <- matrix(NA, ncol = 3, nrow = nsamples)
  rm <- sample(1:1000, (nsamples*2), replace = TRUE)
  for(i in 1:nsamples){
    rm_matrix[i, 1] <- rm[i]
    for(j in  1:2){
      rm_matrix[i,(j +1)] <- abs(round(rnorm(1, rm[i], sd= stand.dev), digits = 0))
    }
  }
  rm_v <- c(rm_matrix[, 1], rm_matrix[,2], rm_matrix[,3])
  index <- seq(1,1000, by=1)
  smp_niid <- dist[index %in% rm_v]
  return(smp_niid)
}

clust.fun.niid <- function(nsamples, dist, stand.dev){
  require(mclust)
  dist.niid <- generate.niid(nsamples, dist, stand.dev = stand.dev)
  clust <- mclustBIC(dist.niid, verbose = F)
  return(clust)
}


######################proportion for non-iid
prop.fun.niid <- function(nsamples, dist, stand.dev){
  clusterBIC <- matrix(NA, ncol = 3, nrow = 1000)
  
  require(mclust)
  
  for(i in 1:500){
    tmp <- clust.fun.niid(nsamples = nsamples, dist = dist, stand.dev = stand.dev)
    # print(tmp[1])
    clusterBIC[i, 1] <- tmp[1]
    clusterBIC[i, 2] <- tmp[2]
    clusterBIC[i, 3] <- tmp[2,2]
  }
  colnames(clusterBIC) <- c("E1", "E2", "V2")
  count <- 0
  for(i in 1:500){
    if(abs(clusterBIC[i,1]) > abs(clusterBIC[i,2])){
      count <- count+1
    }
  }
  
  #proportion 
  p <- count/ 1000
  return(p)
}

####################generate a single distribution 

generate.1.dist <- function(smp_error, mean_vector1,  mean_vector2,
                            sd_vector1, weight1, weight2,
                            sd_vector2, n_samples){
  means <- c(mean_vector1, mean_vector2)
  sds <- c(sd_vector1, sd_vector2)
  
  ind <- sample(1:2, n_samples, replace = TRUE, prob= c(w.m1[i], w.m2[i]))
  x <- rnorm(n_samples, mean = means[ind], sd = sds[ind])
  
  return(x)
}

######################create a time series 

generate.ts <- function(smp.error, mean1, mean2, 
                        sd1, sd2, weight1, weight2, 
                        n.samples, change.by, length.ts, S = 1000){
  
  
  ts.mean <- matrix(NA, nrow = length.ts, ncol = 2) #holds means from each time step
  tmp.dist <- matrix(NA, nrow = S, ncol = 1) #holds 1000 iteration distribution
  # tmp.samp <- matrix(NA, nrow = length(n.samples), ncol = 1) #holds sampled distribtion
  ts <- matrix(NA, nrow = n.samples, ncol = length.ts)
  for(i in 1:length.ts) {
    
    error <- rnorm(n.samples, 0, smp.error)
    
    means <- c(mean1, mean2)
    
    sds <- c(sd1, sd2)
    
    ind <- sample(1:2, S, replace = TRUE, prob = c(weight1, weight2))
    x <- rnorm(S, mean = means[ind], sd = sds[ind])
    
    
    tmp.dist <- x + error
    
    
    index <- seq(1,S, by=1)
    rm <- sample(1:S, n.samples, replace = F)
    
    # tmpdist1 <- sim_psu1_df[,c(col]
    tmp.samp <- tmp.dist[index %in% rm]
    
    
     l <- mean.approx.fun(nsamples = n.samples, d1 = tmp.samp)
     l1 <- apply(l, 2, mean)
     ts.mean[i,] <- l1
    
    # ts[,i] <- tmp.samp 
    
    mean2 <- mean2 + change.by
  }
  
  return(ts.mean)
}


######################change distribution function

#generate and approximate mean for 1 distribution 250 times
#change the distribution 
#approximate mean for the changed distribution 250 times

changeMean <- function(nsamples, error, m1, m2, sd1, sd2, w1, w2, changeBy, niter){
  
  distribution1 <- generate_rm(smp_error = error, mean_vector1 = m1, mean_vector2 = m2, 
                               sd_vector1 = sd1, sd_vector2 = sd2, weight1 = w1, weight2 = w2, 
                               colnm = "D1", n_samples = 1000)
  
  mean.matrix1 <- matrix(NA, nrow =niter, ncol = 2)
  

  for(i in 1:1){
    
    tmp <- mean.approx.fun(nsamples = nsamples, d1 = distribution1, n.iter = niter)
    mean.matrix1[,1] <- tmp[,1]
    mean.matrix1[,2] <- tmp[,2]
    colnames(mean.matrix1) <- c("D1m1", "D1m2")
    
  }
  
  changeM <- changeBy + m1
  
  distribution2 <- generate_rm(smp_error = error, mean_vector1 = changeM, mean_vector2 = m2, 
                               sd_vector1 = sd1, sd_vector2 = sd2, weight1 = w1, weight2 = w2, 
                               colnm = "D2", n_samples = 1000)
  
  mean.matrix2 <- matrix(NA, nrow = niter, ncol = 2)

  for(j in 1:1){
    
    tmp <- mean.approx.fun(nsamples = nsamples, d1 =distribution2, n.iter = niter)
    mean.matrix2[,1] <- tmp[, 1]
    mean.matrix2[,2] <- tmp[, 2]
    colnames(mean.matrix2) <- c("D2m1", "D2m2")
  }


# na1 <- length(which(is.na(mean.matrix[,1]) ==TRUE))
# na2 <- length(which(is.na(mean.matrix[,2]) == TRUE))

df <- data.frame(cbind(mean.matrix1, mean.matrix2))

return(df)
}





