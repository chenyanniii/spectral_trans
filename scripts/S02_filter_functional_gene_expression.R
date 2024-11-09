#############################################################
################## Functional Extraction ####################
#############################################################

setwd(SPECIES) # define in M_master_script.R

######## water, photosynthesis, ABA ############################
## interested function, set up in the master script
## e.g. func=c("water", "photosynthesis", "ABA", "drought", "pathogen", "pigment")


## the function save output in results folder as .csv, which will be imported
for (i in seq(func)){
  func_gene_count(func[i])
}

## load results

water_func_count = read.csv(paste0("results/water_transcript_counts_flt_",SPECIES,".csv"), row.names = 1)

photosynthesis_func_count = read.csv(paste0("results/photosynthesis_transcript_counts_flt_",SPECIES,".csv"), row.names = 1) 

ABA_func_count = read.csv(paste0("results/ABA_transcript_counts_flt_",SPECIES,".csv"), row.names = 1)

drought_func_count = read.csv(paste0("results/drought_transcript_counts_flt_",SPECIES,".csv"), row.names = 1)

pathogen_func_count = read.csv(paste0("results/pathogen_transcript_counts_flt_",SPECIES,".csv"), row.names = 1)

pigment_func_count = read.csv(paste0("results/pigment_transcript_counts_flt_",SPECIES,".csv"), row.names = 1)


## organized function expression data 
## final data matrix named funct_count_tbl_fn

func_summary = rbind(water_func_count, photosynthesis_func_count, ABA_func_count, 
                         drought_func_count, pathogen_func_count, pigment_func_count) 
gene_func = func_summary %>% 
  mutate(gene = rownames(func_summary)) %>%
  dplyr::select(gene, func)
write.csv(gene_func, file = file.path(outdir, paste0("gene_func",SPECIES,".csv")))

func_summary_tbl = rbind(water_func_count, photosynthesis_func_count, ABA_func_count, 
                         drought_func_count, pathogen_func_count, pigment_func_count) %>% 
  dplyr::select(-func) %>% t()

func_count_tbl = as.data.table(func_summary_tbl, keep.rownames = TRUE) 

func_count_tbl_clean = func_count_tbl[,which(duplicated(colnames(func_count_tbl))):=NULL] 

func_count_tbl_fn = func_count_tbl_clean %>% dplyr::select(-rn)

rownames(func_count_tbl_fn) = c(func_count_tbl$rn)

saveRDS(func_count_tbl_fn, file = file.path(outdir, paste0("func_count_tbl_fn_",SPECIES,".rds")))

  
setwd("..")

