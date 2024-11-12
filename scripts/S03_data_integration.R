
#############################################################
################## Data Integration #########################
#############################################################
setwd(SPECIES)

## load the data

X = spectral[,(start_wv-350+1):(end_wv-350+1)]
X = as.matrix(X)

rownames(func_count_tbl_fn) = rownames(spectral)
Y = func_count_tbl_fn
Y = as.matrix(Y)


grid1 <- seq(0.1, 1, length = 10) 
grid2 <- seq(0.1, 1, length = 10)

############################################################
#################### shrink methods ########################
############################################################
shrink_XY = rcc(X, Y, method = 'shrinkage')

saveRDS(shrink_XY, file = file.path(outdir, paste0("shrink_XY_",SPECIES,".rds")))

#### heatmap in pdf (shrink)
pdf(file=paste0("results/cim_shrink_XY_",SPECIES,".pdf"), width = 12, height = 8)
cim(shrink_XY, comp = 1:2, xlab = "genes", ylab = "spectra")
graphics.off()



############################################################
#################### grid tuning ###########################
############################################################
#### first tuning
opt_lamdas_XY = tune_rcc_parallel(X, Y, grid1, grid2)
saveRDS(opt_lamdas_XY, file=file.path(outdir, paste0("opt_lamdas_XY_",SPECIES,".rds")))

rcc_opt_XY = rcc(X, Y, ncomp = 2,
                 lambda1 = opt_lamdas_XY$opt.lambda1,
                 lambda2 = opt_lamdas_XY$opt.lambda2,
                 method = 'ridge')
saveRDS(rcc_opt_XY, file=file.path(outdir, paste0("rcc_opt_XY_",SPECIES,".rds")))


#### second tuning
grid3 = seq((opt_lamdas_XY$opt.lambda1)*0.25, 
            (opt_lamdas_XY$opt.lambda1)*2.5, length =10)
grid4 = seq((opt_lamdas_XY$opt.lambda2)*0.25, 
            (opt_lamdas_XY$opt.lambda2)*2.5, length =10)

opt_lamdas_XY2 = tune_rcc_parallel(X, Y, grid3, grid4)

saveRDS(opt_lamdas_XY2, file=file.path(outdir, paste0("opt_lamdas_XY2_",SPECIES,".rds")))
# opt_lamdas_XY2 = readRDS(file=paste0("results/opt_lamdas_XY2_",SPECIES,".rds"))

rcc_opt_XY2 = rcc(X, Y, ncomp = 3,
              lambda1 = opt_lamdas_XY2$opt.lambda1,
              lambda2 = opt_lamdas_XY2$opt.lambda2,
              method = 'ridge')

saveRDS(rcc_opt_XY2, file=file.path(outdir, paste0("rcc_opt_XY2_", SPECIES, ".rds")))


pdf(file=file.path(outdir, paste0("cim_rcc_opt_XY2_column_", SPECIES, "_large.pdf")), width = 180, height = 120)
cim(rcc_opt_XY2, comp = 1:2, xlab = "genes", ylab = "spectra", cluster = "column")
graphics.off()


setwd("..")
