
################################################################################
######################## Data Wrangling for Figures ############################
################################################################################

############## Default Setting ##############
## Band Range
## wv range
wv_range = matrix(
  c(200,400,700,1100,400,700,1100,2500), ncol = 2
) %>% as.data.frame()
colnames(wv_range) = c("start_wv", "end_wv")
rownames(wv_range) = c("uv", "vis", "nir", "swir")

wv = lapply(1: nrow(wv_range), function(x) seq(wv_range$start_wv[x], wv_range$end_wv[x], 1))
names(wv) = rownames(wv_range)



## Functional Annotation
## gene function table

gene_func = read.csv(file = file.path(outdir, paste0("gene_func", SPECIES, ".csv")), row.names = 1)

#water_func_count = read.csv(file=file.path(outdir, 
#                                           paste0("water_transcript_counts_flt_", SPECIES, ".csv")))

#photosynthesis_func_count = read.csv(file=file.path(outdir,
#                                                    paste0("photosynthesis_transcript_counts_flt_",SPECIES,".csv")))

#ABA_func_count = read.csv(file=file.path(outdir,paste0("ABA_transcript_counts_flt_", SPECIES, ".csv")))

#drought_func_count = read.csv(file=file.path(outdir,paste0("drought_transcript_counts_flt_", SPECIES, ".csv")))

#pigment_func_count = read.csv(file=file.path(outdir,paste0("pigment_transcript_counts_flt_", SPECIES, ".csv")))

#pathogen_fun_count = read.csv(file=file.path(outdir, paste0("pathogen_transcript_counts_flt_", SPECIES, ".csv")))


#list_func_gene = list(water_func_count,photosynthesis_func_count, ABA_func_count,
#                      drought_func_count, pigment_func_count, pathogen_fun_count)
#names(list_func_gene) = c("water", "photo", "ABA", "drought", "pigment", "pathogen")


#func_gene = sapply(list_func_gene, function(i) as.data.frame(i[,1]))
#names(func_gene) = names(list_func_gene)

#saveRDS(func_gene, file = file.path(outdir, paste0("func_gene_", SPECIES, ".rds")))


############## Regularized Canonical Correlation Results  ##############
## Import Data
opt_lamdas_XY2 = readRDS(file=file.path(outdir, paste0("opt_lamdas_XY2_",SPECIES,".rds")))
rcc_opt_XY2=readRDS(file=file.path(outdir, paste0("rcc_opt_XY2_", SPECIES, ".rds")))

## supermatrix calculation
mat = rcc_opt_XY2
comp = 1:2
bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
# correlation coefficients for each column of mat$X with the vector bisect
cord.X = cor(mat$X, bisect, use = "pairwise") 
# correlation coefficients for each column of mat$Y with the vector bisect
cord.Y = cor(mat$Y, bisect, use = "pairwise") 
# canonical correlation analysis (CCA), testing the cross-correlation between two set of variables
XY.mat = as.matrix(cord.X %*% t(cord.Y)) 

saveRDS(XY.mat, file=file.path(outdir, paste0("XY.mat_", SPECIES, ".rds")))

############## Loading: Gene Expression
gene_loading_corr = rcc_opt_XY2[["loadings"]][["Y"]][,1] %>% 
  as.data.frame()  

colnames(gene_loading_corr) = c("gene_loading")

gene_loading_corr2 = gene_loading_corr %>% 
  mutate(gene = rownames(gene_loading_corr)) %>%
  mutate(abs_loading = abs(gene_loading))

gene_loading_func = left_join(gene_loading_corr2, gene_func, by = join_by(gene)) %>% 
  mutate(gene_func = func) %>% as.data.frame() %>%
  dplyr::select("gene", "gene_func", "gene_loading", "abs_loading")


############## Loading:Spectral Reflectance

spectra_loading_corr = rcc_opt_XY2[["loadings"]][["X"]][,1] %>% 
  as.data.frame()  

colnames(spectra_loading_corr) = c("spectra_loading")


spectra_loading_corr2 = spectra_loading_corr %>%
  mutate(wv = rownames(spectra_loading_corr)) %>%
  mutate(loading_rank = rank(-spectra_loading)) %>%
  mutate(abs_loading = abs(spectra_loading))

spectra_loading_corr2$wv = as.numeric(spectra_loading_corr2$wv)

spectra_loading_corr2 = spectra_loading_corr2 %>%
  mutate(wv_range = cut(
    spectra_loading_corr2$wv,
    breaks = c(200, 400, 700, 1100, 2500),
    labels = c("uv", "vis", "nir", "swir"),
    right = FALSE)
  )

spectra_loading_corr3 = cbind(spectra_loading_corr2, XY.mat) 

## long format by function
spectra_loading_corr3_long = spectra_loading_corr3 %>% 
  pivot_longer(cols = (6:ncol(spectra_loading_corr3)),
               names_to = "gene",
               values_to = "XY_correlation")

spectra_loading_corr3_long_2 = left_join(spectra_loading_corr3_long, gene_func, by = join_by(gene)) %>% 
  mutate(gene_func = func) %>% as.data.frame() %>%
  dplyr::select("wv", "wv_range", "spectra_loading", "loading_rank", "gene", 
                "gene_func", "XY_correlation")

spectra_loading_fnl = spectra_loading_corr3_long_2
saveRDS(spectra_loading_fnl, file = file.path(outdir, paste0("spectra_loading_fnl_", SPECIES,".rds")))

spectra_loading_fnl$XY_correlation = as.numeric(spectra_loading_fnl$XY_correlation)
spectra_loading_range = spectra_loading_fnl %>% 
  mutate(corr_range = cut(
    abs(spectra_loading_fnl$XY_correlation),
    breaks = c(0, 0.3, 0.5, 0.7, 1),
    labels = c("n.a.","weak", "medium", "strong")) 
  ) 

spectra_loading_range_flt = spectra_loading_range %>% filter(corr_range != "n.a.")


## spectra_gene_0.5
spectra_gene_0.5 = spectra_loading_fnl %>%
  filter(abs(XY_correlation) >0.5)


spectra_gene_0.5_summary = spectra_gene_0.5 %>% 
  count(wv_range, gene_func) %>%
  group_by(wv_range) %>%
  mutate(prop_wv = n/sum(n))

## filter moderate to strong canonical genes (cor > 0.3)
cutoff = 0.02

cor_brks_low = apply(XY.mat, 2, function(x) sum(ifelse(x>=0.3,1,0)))/nrow(XY.mat)

cor_brks_flt_low = cor_brks_low[which(cor_brks_low >= cutoff)]

df_cor_brks_low = data.frame(cor_perc = cor_brks_low)
df_cor_brks_flt_low = data.frame(cor_perc = cor_brks_flt_low)

## filter strong canonical genes (cor > 0.7)
cor_brks_strong = apply(XY.mat, 2, function(x) sum(ifelse(x>=0.7,1,0)))/nrow(XY.mat)

cor_brks_flt_strong = cor_brks_strong[which(cor_brks_strong >= cutoff)]

spectr_cor_strong = spectra_loading_corr3_long_2 %>%
  dplyr::filter(gene %in% names(cor_brks_flt_strong))


