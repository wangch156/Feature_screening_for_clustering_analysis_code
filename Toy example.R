library(foreach)
library(BiocParallel)
library(parallel)
library(MASS)
library(nloptr)
library(doParallel)
library(mclust)
Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/FSEM.cpp")
Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/compute_pl.cpp")
Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/newtown_homo.cpp")


##################### EM-test statistic #########################

em.stat <- function(x,alpha.ini,k0,C,labels,group.num,earlystop = NA, is_ecdf = F, 
                    cut_max = 50){
  require(nloptr)
  # C_time = Sys.time()
  
  # Calculate EM statistics
  # k0 is the number of iterations, and alpha.ini are the initial alpha values
  # C is the coefficient of the penalty term
  # The first n components of `theta` are means, while the last n components are dispersion parameters (r)
  # Labels should be in the form of 1,2,...,n
  # Each row of alpha.ini contains a set of initial values
  
  x_table = as.matrix(table(x))
  x_reduce = as.numeric(rownames(x_table))
  reduce_num = as.vector(x_table[,1])
  if (length(unique(labels))!=group.num){
    print("The number of different labels does not match the total number of groups!")
    return(NULL)
  }
  M.stat <- matrix(nrow = 1+k0,ncol = nrow(alpha.ini))
  w_homo = rep(1/length(x),length(reduce_num))
  m_homo = mean(x)
  var_homo = var(x)
  fit.nb <- function(x){
    # this function returns the MLE of mean and r
    obj.f <- function(theta){-sum(dnbinom(x,theta[2],mu = theta[1],log = T))}
    deriv.pl <- function(theta){nl.grad(x0 = theta,fn = obj.f)}
    nloptr(x0 = c(min(max(mean(x),1e-9*1.01),2999),5),eval_f = obj.f,lb=rep(1e-9,2),ub = rep(3000,2),opts = list("algorithm"="NLOPT_LN_NELDERMEAD","ftol_rel" = 1e-9,"maxeval" = 5000,"maxtime" = 200,"xtol_rel" = 1e-4))$solution
  }
  if(cut_max != max(x)){
    theta.0 = homo_opt(w_homo, x_reduce, reduce_num,  newtown_step = 100, r = 1)
    theta_hat_0 = c(rep(theta.0[1],group.num),rep(theta.0[2],group.num))
  } else {
    theta.0 <- fit.nb(x)
    theta_hat_0 = c(rep(theta.0[1],group.num),rep(theta.0[2],group.num))
  }
  pl0 = compute_pl(theta_hat_0, as.integer(x_reduce), reduce_num, rep(1/group.num,group.num), group.num, C)
  
  mu <- numeric(group.num)
  r <- numeric(group.num)
  if(cut_max != max(x)){
    for (i in 1:group.num) {
      x_i = x[labels==i]
      reduce_num_i = rep(1, length(x_i))
      w_homo_i = rep(1/length(x_i), length(x_i))
      m_homo = mean(x_i)
      var_homo = var(x_i)
      theta.0 = homo_opt(w_homo_i, x_i, reduce_num_i,  newtown_step = 100, r = 1)
      mu[i] <- theta.0[1]
      r[i] <- theta.0[2]
    }
  } else {
    for (i in 1:group.num){
      mu[i] <- mean(x[labels==i])
      r[i] <- fit.nb(x[labels==i])[2]
    }
  }
  
  mu[mu<=1e-9] <- 1e-9*1.01
  r[r<=1e-9] <- 1e-9*1.01
  x_order = order(x, decreasing = F)
  x = x[x_order]
  labels = labels[x_order]
  
  
  for (j in 1:nrow(alpha.ini)){
    # Calculate the initial M values
    
    k <- 0
    alpha.old <- alpha.ini[j,]
    
    theta.old <- c(mu,r)
    pl.value = compute_pl(theta.old, as.integer(x_reduce), reduce_num, c(alpha.old, 1-sum(alpha.old)), group.num, C)
    M.stat[1,j] <- 2 * (pl.value-pl0)
    w0 = matrix(0.001,nrow = length(x),ncol = group.num)
    for (i in 1:group.num){
      w0[labels==i,i] <- 1-0.001*(group.num-1)
    }
    w <- w0
    w0_reduce = matrix(0, ncol = group.num, nrow = length(reduce_num))
    for (g in 1:group.num) {
      w_g = w0[,g]
      w0_reduce[,g] = reduce_w(w_g, reduce_num)
    }
    w_orign = w0
    w_reduce = w0_reduce
    
    ret = FSEM(group.num, 0.0, C , w_reduce, w0_reduce, x_reduce,as.integer(x_reduce), reduce_num,theta.old, k0,50,cut_max)
    # print(Sys.time() - C_time)
    M.stat[2:(k0+1),j] = 2*(ret- pl0)
    
  }
  # Calculate EM statistics
  return(c(1-pchisq(max(M.stat),3),max(M.stat)))
}

##################### data generation (mixture of negative binomial) #####################

set.seed(1)
samples = 1000
DE.gene <- 20
nonDE.gene <- 4880
All_gene = DE.gene + nonDE.gene
counts <- matrix(rep(0,(DE.gene+nonDE.gene)*samples),ncol = samples)
colnames(counts) <- paste("Cell",seq(1:samples),sep = "")
rownames(counts) <- paste("Gene",seq(1:(DE.gene+nonDE.gene)),sep = "")
generate.counts <- function(mu,size, n){
  rnbinom(n,size  = size, mu = mu)
}
dispersion_par = c(5,6)
mu.list = rep(0, All_gene)
mu <- sample(exp(seq(log(2),log(5),length.out = 100)),nonDE.gene,replace = T)
size = runif(nonDE.gene, dispersion_par[1],dispersion_par[2])
mu.list[1:nonDE.gene] <- mu
for (k in 1:nonDE.gene){
  counts[k,] <- generate.counts(mu[k],size[k], samples)
}
mu <- sample(exp(seq(log(2),log(5),length.out = 100)),DE.gene,replace = T)
size = runif(DE.gene, dispersion_par[1],dispersion_par[2])
diff_par = c(7,8)
Mean_diff = runif(DE.gene,diff_par[1],diff_par[2])
mu2 = mu + Mean_diff
mu.list[(nonDE.gene+1):(nonDE.gene+DE.gene)] <- mu
alpha_star = c(0.5,0.125,0.125,0.125,0.125)
ture_label = sample(1:5, samples, replace = T, alpha_star)
for (k in 1:DE.gene){
    counts[(k+nonDE.gene),which(ture_label == 1)] <- generate.counts(mu[k],size[k], length(which(ture_label == 1)))
    if ((k %% 4) == 0) index.tmp <-which(ture_label == 2)
    if ((k %% 4) == 1) index.tmp <- which(ture_label == 3)
    if ((k %% 4) == 2) index.tmp <- which(ture_label == 4)
    if ((k %% 4) == 3) index.tmp <- which(ture_label == 5)
    
    counts[(k+nonDE.gene),index.tmp] <- generate.counts(mu2[k],size[k], length(index.tmp))
    index.tmp.c = setdiff(setdiff(1:samples, which(ture_label == 1)), index.tmp)
    counts[(k+nonDE.gene),index.tmp.c] <- generate.counts(mu[k],size[k], length(index.tmp.c))
}

###########  pre-clustering (No screening) ##############      

set.seed(1)
normlize_counts = apply(counts, 1, scale)
pca_counts = prcomp(normlize_counts)
res = kmeans(pca_counts$x[,1:10], centers = 5) 
cluster.labels = res$cluster
cluster.labels <- factor(cluster.labels)
levels(cluster.labels) <- 1:length(levels(cluster.labels))
ari_origin = adjustedRandIndex(cluster.labels, ture_label)

############  Oracle screening ###############

set.seed(1)
ari_gold = 0
for (clu in 1:5) {
  normlize_counts_gold = apply(counts[(nonDE.gene) : All_gene, ], 1, scale)
  pca_counts_gold = prcomp(normlize_counts_gold)
  res_gold = kmeans(pca_counts_gold$x[,1:10], centers = 5) 
  cluster.labels_gold = res_gold$cluster
  cluster.labels_gold <- factor(cluster.labels_gold)
  levels(cluster.labels_gold) <- 1:length(levels(cluster.labels_gold))
  ari_gold = max(adjustedRandIndex(cluster.labels_gold, ture_label),ari_gold)
}


#################### EM-test ################ 

alpha.label <- numeric(length(levels(cluster.labels))-1)
for (g in 1:(length(alpha.label))) {
  alpha.label[g] <- sum(cluster.labels==g)/length(cluster.labels)
}

cl <- makeCluster(getOption("cl.cores", 5))
registerDoParallel(cl)

em.result = foreach (j = 1:All_gene, .combine = "cbind") %dopar% {
  Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/FSEM.cpp")
  Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/compute_pl.cpp")
  Rcpp::sourceCpp("~/Desktop/project/NB test/NBsimulation/Toy example/newtown_homo.cpp")
  x = counts[j,]
  res = em.stat(x=x, alpha.ini=rbind(alpha.label),k0=100,C=1e-5,
                labels = cluster.labels,group.num = length(levels(cluster.labels)),
                earlystop = 1e-5, is_ecdf = F, cut_max = (max(x)+1))
  return(res)
}

stopCluster(cl)
em_p = em.result[1,]
em_p_adjust = p.adjust(em_p, method = "BH")
ind_em = which(em_p_adjust < 0.01)
set.seed(1)
ari_em = 0
for (clus in 1:5) {
  normlize_counts_em = apply(counts[(ind_em), ], 1, scale)
  pca_counts_em = prcomp(normlize_counts_em)
  res_em = kmeans(pca_counts_em$x[,1:min(10, length(ind_em))], centers = 5)
  cluster.labels_em = res_em$cluster
  cluster.labels_em <- factor(cluster.labels_em)
  levels(cluster.labels_em) <- 1:length(levels(cluster.labels_em))
  ari_em = max(adjustedRandIndex(cluster.labels_em, ture_label),ari_em)
}

