# Load necessary libraries if not already installed

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to disease GWAS summary file (SNP, A1, A2, Z) [required]"),
  make_option("--pos_dir", action="store", default=NA, type='character',
              help="*.pos files folder [required]"),  
  make_option("--model_dir", action="store", default=NA, type='character',
              help="*.RDat files folder [required]"),
  make_option("--ref_ld_dir", action="store", default=NA, type='character',
              help="pre extracted LD reference data [required]"),
  make_option("--TSS", action="store", default=NA, type='character',
              help="TSS file [required]"),                
  make_option("--output_dir", action="store", default=NA, type='character',
              help="Path to output files [required]")

)

# Get command-line arguments
opt = parse_args(OptionParser(option_list=option_list))

if (is.na(opt$sumstats) || is.na(opt$pos_dir) || is.na(opt$model_dir) || is.na(opt$ref_ld_dir)|| is.na(opt$output_dir) ) {
    cat("Warning: One or more required options are missing.\n")
    q()
}

source("/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/support.R")





########################
### association test ###
########################
association_test(
    sumstats = opt$sumstats, 
    pos_dir = opt$pos_dir, 
    model_dir = opt$model_dir, 
    ref_ld_dir = opt$ref_ld_dir, 
    TSS_file = opt$TSS,
    output_dir = opt$output_dir)


ass_df <- data.frame()
files <- Sys.glob(file.path(opt$output_dir, '*.txt'))
for (file in files){
    tmp <- fread(file, data.table = F)
    ass_df <- rbind(ass_df, tmp)
}
ass_df <- ass_df[order(ass_df$PWAS.P),]
ass_df <- filter(ass_df, PWAS.P != 'NA')
ass_df <- filter(ass_df, MODELCV.R2 >= 0.01)
ass_df$FDR <- p.adjust(ass_df$PWAS.P, method = 'BH')
ass_df$Bonferroni_p <- p.adjust(ass_df$PWAS.P, method = 'bonferroni')
message( 'A total of ', sum(ass_df$FDR < 0.05), ' proteins significantly associated with disease at FDR < 0.05.')
message( 'A total of ', sum(ass_df$Bonferroni_p < 0.05), ' proteins significantly associated with disease at bonferroni correction p-value < 0.05.')
write.table(ass_df, file.path(opt$output_dir, 'all_association.out'), row.names = F, sep = '\t', quote = F)
