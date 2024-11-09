############### function for fucntional expression extraction ##################

## This function was written to extract functional expression in csv
## csv file: have name pattern: func + "_keyword_ASA_acer10.csv"
## filter: mean > 2; stdv < (mean*0.5)
## adaption could apply to other filter critiers and patterns to adapt to other 



## func=c("water", "photosynthesis", "ABA")

# for (i in seq(func)){
#  func_gene_count(func[i])
#}



func_gene_count = function(func){
  func_gene = read.csv(paste0("data/rna/",func,"_keyword_",SPECIES,"_",ref_rna,".csv"))
  func_trans_count = trans_count[,intersect(colnames(trans_count),func_gene[,2])] 
  
  df = data.frame(
    mean =  apply(func_trans_count, 2, mean), 
    stdv =  apply(func_trans_count, 2, sd)
  ) %>% t()
  
  func_trans_count_flt = rbind(func_trans_count, df) %>% 
    t() %>% as.data.frame() %>%
    filter(mean >2) %>% dplyr::filter(stdv < (mean*0.5)) %>%
    mutate(func = func) %>% 
    dplyr::select(-mean, -stdv)
  
  write.csv(func_trans_count_flt, 
            file = file.path(outdir, paste0(func,"_transcript_counts_flt_",SPECIES,".csv")))
  
  }


