library(foreach)
library(BiocParallel)
library(parallel)
library(MASS)
library(nloptr)
library(doParallel)
library(mclust)
library(Festem)


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

alpha.label <- numeric(length(levels(cluster.labels)))
for (g in 1:(length(alpha.label))) {
  alpha.label[g] <- sum(cluster.labels==g)/length(cluster.labels)
}

cl <- makeCluster(getOption("cl.cores", 5))
registerDoParallel(cl)

em.result = foreach (j = 1:All_gene, .combine = "cbind", .packages = "Festem") %dopar% {
  x = counts[j,]
  res = em_test(x=x,G=length(levels(cluster.labels)), alpha.init = alpha.label,
                labels = cluster.labels, prior.weight = 0.0)
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

