



########## load packages ##########
library(dplyr)
library(stringr)
library(data.table)
library(parallelly)
library(mixOmics)
library(foreach)
library(doParallel)
library(tidyr)
library(here)


########## Define Global Parameters ##########
SPECIES = "AR"
ref_rna = "acer10"

func=c("water", "photosynthesis", 
       "ABA", "drought", "pathogen", "pigment")

start_wv = 400
end_wv = 2400

###### define output path
# user defined PATH - e.g. "~/scratch/PLSR"
output_dir <- paste0(here(),"/",SPECIES,"/results")
#--------------------------------------------------------------------------------------------------#


#--------------------------------------------------------------------------------------------------#
### Set working directory
if (output_dir=="tempdir") {
  outdir <- tempdir()
} else {
  if (! file.exists(output_dir)) dir.create(output_dir,recursive=TRUE)
  outdir <- file.path(path.expand(output_dir))
}


########## Create Parallel Environment ##########

######## Sample scripts to parallelly run tune.rcc{mixOmics}
# Detect the number of available cores
num_cores = detectCores()
num_cores_hpc = availableCores()

######## Create a cluster with the number of cores
if (num_cores == num_cores_hpc){
  final_cores = (num_cores - 1) # Leaving one core for the OS
} else{
  final_cores = num_cores_hpc
}

cl = makeCluster(final_cores)
# Register the parallel backend
registerDoParallel(cl)



########################################################################
########################  STEP 1: load master data #####################
########################################################################
gene_file = paste0("data/rna/gene_count_matrix_",SPECIES,"_", ref_rna,".csv")
trans_file = paste0("data/rna/transcript_count_matrix_",SPECIES,"_", ref_rna,".csv")

source(paste0("scripts/S01_load_master_data_",SPECIES,".R"))

print ("Finish STEP 1: load master data")
########################################################################
############ STEP 2: filter functional gene expression ################# 
########################################################################
## define interested function and file location
## will be used in func_functional_gene_count.R


func_file = (paste0("data/rna/",func,"_keyword_", SPECIES,"_", ref_rna, ".csv"))

source("scripts/func_functional_gene_count.R")

source("scripts/S02_filter_functional_gene_expression.R")

print("Finish STEP 2: filter functional gene expression")
########################################################################
############ STEP 3: data integration in parallel ######################
########################################################################



source("scripts/func_parallel_tune_rcc.R")

source("scripts/S03_data_integration.R")

print("Finish STEP 3: data integration in parallel")

########################################################################
############ STEP 4: data wrangling and figure Generation ##############
########################################################################

########### Data Wrangling
source("scripts/Figure_data_wrangling2.R")


########### Figure Generation
source("scripts/Figures.R")

print("Finish STEP 4: data wrangling and figure generation")

# Stop the parallel backend
stopImplicitCluster()
