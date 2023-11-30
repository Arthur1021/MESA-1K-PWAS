# Load necessary libraries if not already installed

# suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--geno_Typed", action="store", default=NA, type='character',
              help="Path to Genotyped genotype files folder, used for PCA and relatedness analysis [required]"),
  make_option("--geno_Imputed", action="store", default=NA, type='character',
              help="Path to Imputed genotype files folder, used for establish prediction models [required]"),
  make_option("--GRCh_version", action="store", default=NA, type='numeric',
              help="genotype reference version (37 or 38)[required]"),  
  make_option("--hyphen", action="store", default=NA, type='character',
              help="hyphen between Chromosome and position ('_' or ':')[required]"),
  make_option("--ethnic", action="store", default=NA, type='character',
              help="ethnic information ('Eur' or 'non-Eur')[required]"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to phenotype file [required]"),
  make_option("--covariate", action="store", default=NA, type='character',
              help="Path to covariate file [required]"),
  make_option("--annotation", action="store", default=NA, type='character',
              help="Path to annotation file (must have columns ID, Chr, TSS) [required]"),
  make_option("--ref_ld", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]")
)


# Get command-line arguments
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$geno_Typed) || is.na(opt$geno_Imputed) || is.na(opt$GRCh_version) || is.na(opt$hyphen) || is.na(opt$ethnic)|| is.na(opt$pheno) || is.na(opt$covariate) || is.na(opt$annotation)|| is.na(opt$ref_ld)|| is.na(opt$out)) {
    cat("Warning: One or more required options are missing.\n")
    q()
}

# Extract information from the config list
genotype_Genotyped_data_path <- opt$geno_Typed
genotype_Imputed_data_path <- opt$geno_Imputed
phenotype_data <- opt$pheno
covariate_data_path <- opt$covariate
annotation_data_path <- opt$annotation



source("/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/support.R")

############################################
### read cov file, will perform analysis ###
############################################

ID <- opt$out
message ('Analysing ', covariate_data_path)
if (!dir.exists(ID)){
    dir.create(ID,recursive = TRUE)
}
######################################
### extract genotype for subjects ###
#####################################
df <- fread(covariate_data_path, data.table = F)
geno_filter_with_sample_list(
    geno_folder = genotype_Genotyped_data_path,
    geno_threshold = 0.05,
    maf_threshold = 0.05,
    hwe_threshold = 0.000005,
    sample_list = df[,1],
    out_folder = file.path(ID,'01_extract_subject_with_geno_pheno')
)

#########################################################
### Check independent subject for each race/ethnicity ###
#########################################################
run_relatedness_analysis(
    input_QC_genotype_file =  file.path(ID,'01_extract_subject_with_geno_pheno/chr1-22_SNP_clean'),
    output_dir = file.path(ID,'02_check_independent_subjects')

)

####################################################
### re-extract genotype for independent subjects ###
####################################################
# get keep subject ID
drop_list <- fread(file.path(ID,'02_check_independent_subjects/drop_list.txt'), data.table = F)
keep_subject <- df[!df[,1] %in% drop_list$V2,]
keep_subject <- select(keep_subject, colnames(keep_subject)[1])
if (!dir.exists(file.path(ID,'03_independent_subjects_genotype'))) {
    dir.create(file.path(ID,'03_independent_subjects_genotype'))
}
write.table(keep_subject, file.path(ID,'03_independent_subjects_genotype/independent_subjects.txt'), row.names = F, col.names = F, sep = '\t', quote = F)

# re-extract
geno_filter_with_sample_list(
    geno_folder = genotype_Imputed_data_path,
    geno_threshold = 0.05,
    maf_threshold = 0.05,
    hwe_threshold = 0.000005,
    sample_list = keep_subject[,1],
    out_folder = file.path(ID,'03_independent_subjects_genotype')
)
# convert chr:pos to rsID
convert_chr_pos_to_rsID(
    input_genotype_file = file.path(ID,'03_independent_subjects_genotype/chr1-22_SNP_clean'),
    output_genotype_file = file.path(ID,'03_independent_subjects_genotype/chr1-22_SNP_clean_rsID'),
    GRCh_version = opt$GRCh_version,
    hyphen = opt$hyphen   # hyphen is ":" or "_"
)

##########################################
### generate PCA matrix using smartpca ###
##########################################
if (opt$ethnic == 'European'){
    Ethnic_classic <- 'Eur'
}else{
    Ethnic_classic <- 'non-Eur'
}
generate_PCA(
    input_genotype_file = file.path(ID,'03_independent_subjects_genotype/chr1-22_SNP_clean_rsID'),
    Ethnic = Ethnic_classic,
    output_dir = file.path(ID,'04_generate_PCA')
)

##########################
### combine covariates ###
##########################
cov1 <- fread(covariate_data_path, data.table = F)
PCA <- fread(paste0(file.path(ID,'04_generate_PCA'), '/independent_SNP_PCA.pca.evec'), data.table = F, header = F)
row.names(PCA) <- PCA$V1
PCA <- PCA[,1:11]
colnames(PCA)[1] <- 'IID'
colnames(PCA)[2:11] <- paste0('PCA',1:10)
cov_matrix <- left_join(PCA, cov1, by = c('IID' = colnames(cov1)[1]))
write.table(cov_matrix, paste0(file.path(ID,'04_generate_PCA'),'/PCA_cov_matrix.txt'), row.names = F, sep = '\t', quote = F)

#######################################################
### phenotype regression out covariates and invrank ###
#######################################################

pheno_log_invrank_residual_invrank(
    phenotype_data = phenotype_data,
    cov_matrix = paste0(file.path(ID,'04_generate_PCA'),'/PCA_cov_matrix.txt'), 
    output_dir = file.path(ID,'05_pheno_adjustment')
)

################################################################
### glm using plink, cis region FDR < 0.05 trans region 5e-9 ###
################################################################

glm(
    input_QC_genotype_file = file.path(ID,'03_independent_subjects_genotype/chr1-22_SNP_clean_rsID'), 
    annotation_file = annotation_data_path, 
    cis_flank_region = 1000000, 
    trans_threshold = 5e-9, 
    phenotype_data_folder = file.path(ID,'05_pheno_adjustment'), 
    output_dir = file.path(ID,'06_glm'))

#################################
### extract up and down 100kb ###
#################################
extract_potential_SNP_predictors(
    cis_glm_folder = file.path(ID,'06_glm/cis'), 
    trans_glm_folder = file.path(ID,'06_glm/trans'), 
    flank_region = 100000, 
    FDR_threshold = 0.05,
    input_QC_genotype_file = file.path(ID,'03_independent_subjects_genotype/chr1-22_SNP_clean_rsID'), 
    strict_SNP = paste0(opt$ref_ld, '.bim'), 
    output_dir = file.path(ID,'07_extract_potential_SNP_predictors')
)


#######################################
### build genetic prediction models ###
#######################################

establish_prediction_models(
    genotype_file_folder = file.path(ID,'07_extract_potential_SNP_predictors'),
    pheno_folder = file.path(ID,'05_pheno_adjustment'), 
    gcta_PATH = '/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/bin/gcta_nr_robust', 
    gemma_PATH = '/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/bin/gemma-0.98.5-linux-static-AMD64', 
    hsq_p = 1,
    PANEL = ID,
    reference_1000G = opt$ref_ld,
    output_dir = file.path(ID,'08_establish_prediction_models')
)


