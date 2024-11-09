
# setwd("/Users/yannchen/Desktop/Research_Projects_ND/two_maples")

########## load packages ##########
# library(dplyr)
# library(stringr)
# library(foreach)
# library(doParallel)
# library(mixOmics)



#############################################################
################## Load Master Data #########################
#############################################################

setwd(SPECIES) # define in M_master_script.R

########## load data ##########
##### gene counts

gene.count=read.csv(gene_file, row.names=1)

gene.count = t(gene.count)

rns = rownames(gene.count) 
rns = sapply(strsplit(rns, "_"), function(x) str_c(x[1], "_", x[2]))

rownames(gene.count) = rns
gene.count = log(gene.count+1)



##### transcript counts
trans_count=read.csv(trans_file, row.names=1) %>% t()

rns_t = rownames(trans_count) 
rns_t = sapply(strsplit(rns_t, "_"), function(x) str_c(x[1], "_", x[2]))

rownames(trans_count) = rns_t

trans_count=log(trans_count+1)



##### spectral data
setwd("data/asd.data") ## change directory to a subdirectory in data

listy = list.files()

listy2 = sapply(strsplit(listy, "_"), function(x) x[2])
listy2 = sapply(strsplit(listy2, ".csv"), function(x) str_c("ns_", x[1]))
listy2 = gsub("ns_103103", "ns_10310", listy2)
listy2 = gsub("ns_123678", "ns_12367", listy2)
listy2 = gsub("ns_278041", "ns_27804", listy2)


first = read.csv(listy[1], row.names=1)

spectral = apply(first, 1, mean,na.rm=T)

for(i in 2:length(listy)){
  
  second = read.csv(listy[i], row.names=1)
  
  second.sum = apply(second, 1, mean,na.rm=T)
  
  spectral = rbind(spectral, second.sum)
  
}

rownames(spectral) = listy2

spectral = spectral[rownames(gene.count),]

setwd("../..") ## back to master directory

saveRDS(spectral, file = file.path(outdir,paste0("spectral_",SPECIES,".rds")))

setwd("..")
