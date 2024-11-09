
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
rns = gsub("_3","", rns)
rns = gsub("_5_2","", rns)
rns = gsub("_5","", rns)
rns = gsub("RLT","", rns)

rownames(gene.count) = rns
gene.count = log(gene.count+1)



##### transcript counts
trans_count=read.csv(trans_file, row.names=1) %>% t()

rns_t = rownames(trans_count) %>% str_replace_all("_3","") %>% str_replace_all("_5_2","") %>%
  str_replace_all("_5","") %>% str_replace_all("RLT","")

rownames(trans_count) = rns_t

trans_count=log(trans_count+1)



##### spectral data
setwd("data/asd.data") ## change directory to a subdirectory in data

listy = list.files()
listy2 = gsub(".csv","", listy)
listy2 = gsub("june_", "ns_", listy2)
listy2 = gsub("_redo", "", listy2)
listy2 = gsub("ns_291327", "ns_29132", listy2)
listy2 = gsub("ns_169441", "ns_16944", listy2)



first = read.csv(listy[1], row.names=1)

spectral = apply(first, 1, mean,na.rm=T)

for(i in 2:length(listy)){
  
  second = read.csv(listy[i], row.names=1)
  
  second.sum = apply(second, 1, mean,na.rm=T)
  
  spectral = rbind(spectral, second.sum)
  
}

rownames(spectral) = listy2

spectral = spectral[rownames(gene.count),]

setwd("../..") ## back to master directory of species (ASA)

saveRDS(spectral, file = file.path(outdir, paste0("spectral_",SPECIES,".rds")))

setwd("..") ## back to master directory
