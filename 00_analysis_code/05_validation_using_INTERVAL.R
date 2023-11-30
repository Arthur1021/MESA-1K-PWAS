# validate non-Hispanic White models using INTERVAL data

library(data.table)
library(dplyr)
library(foreach)
library(doParallel)
library(stringr)

if (!dir.exists('05_validation_using_INTERVAL')){
    dir.create('05_validation_using_INTERVAL')
}

# prepare INTERVAL genotype data
if (!dir.exists('05_validation_using_INTERVAL/01_INTERVAL_genotype')){
    dir.create('05_validation_using_INTERVAL/01_INTERVAL_genotype')
}

# extract only SNP
files <- Sys.glob('/mnt/lvm_vol_1/txu/projects/pwas/EGAD00010001544/INTERVAL_SOMALOGIC_POSTQC_chrom*final_v1.bgen.bim')
for (file in files){
    chr_number <- str_split(last(str_split(file, '/')[[1]]), '_')[[1]][[5]]
    infile <- gsub('.bgen.bim', '.bgen', file)
    system2(
      command = "plink2",
      args = c(
        "--bfile", infile,
        "--snps-only",
        "--make-bed",
        "--out", paste0('05_validation_using_INTERVAL/01_INTERVAL_genotype/chr',chr_number)
      )
    )
}

# generate rsID convert file
rsID <- fread('/mnt/lvm_vol_1/dghoneim/EGA/EGAD00001004080/EGAF00001998173/Variant_info.tsv', data.table = F)
rsID <- select(rsID, VARIANT_ID, rsID)
write.table(rsID, '05_validation_using_INTERVAL/01_INTERVAL_genotype/rsID_convert.txt', col.names = F, row.names = F, sep = '\t', quote = F)

# convert to rsID
files <- Sys.glob('05_validation_using_INTERVAL/01_INTERVAL_genotype/chr*.bim')
for (file in files){
    chr_number <- last(str_split(file, '/')[[1]]) %>% gsub('.bim', '', .)
    infile <- gsub('.bim', '', file)
    system2(
      command = "plink2",
      args = c(
        "--bfile", infile,
        "--snps-only",
        "--update-name", "05_validation_using_INTERVAL/01_INTERVAL_genotype/rsID_convert.txt",
        "--rm-dup", "exclude-all",
        "--make-bed",
        "--out", paste0('05_validation_using_INTERVAL/01_INTERVAL_genotype/',chr_number, '_rsID')
      )
    )
}

# merge all chromosomes
# Combine chromosomes 1-22 into a single file
all_chromosome <- file.path("05_validation_using_INTERVAL/01_INTERVAL_genotype/all_chromosome.txt")
writeLines(paste0("05_validation_using_INTERVAL/01_INTERVAL_genotype/chr", 2:22, "_rsID"), con = all_chromosome)

system2(
    command = "plink",
    args = c(
        "--bfile", '05_validation_using_INTERVAL/01_INTERVAL_genotype/chr1_rsID',
        "--merge-list", all_chromosome,
        "--make-bed",
        "--allow-no-sex",
        "--out", paste0("05_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID")
    )
)




annotation_MESA <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
# select only single target
annotation_MESA <- annotation_MESA[!duplicated(annotation_MESA$SomaId),]
annotation_MESA <- select(annotation_MESA, UniProt, SomaId)


annotation_INTERVAL <- fread('/mnt/lvm_vol_2/sliu/project/INTERVAL_pwas/2_gene_start_end_37/001_SOMALOGIC_GWAS_protein_info.csv', data.table = F)
annotation_INTERVAL <- select(annotation_INTERVAL, UniProt, SOMAMER_ID)

intersect <- inner_join(annotation_MESA, annotation_INTERVAL, by = c('UniProt'))

measured_pheno_all <- fread('/mnt/lvm_vol_1/dghoneim/EGA/EGAD00001004080/EGAF00001994131/INTERVAL_SOMALOGIC_POSTQC_GWASIN_PROTEINDATA_v1.tsv', data.table = F)
measured_pheno_all <- select(measured_pheno_all, c('Sample_Name',intersect$SOMAMER_ID))


# prepare INTERVAL measured protein level
if (!dir.exists('05_validation_using_INTERVAL/02_INTERVAL_protein')){
    dir.create('05_validation_using_INTERVAL/02_INTERVAL_protein')
}


for(i in 2:ncol(measured_pheno_all)){
    tmp <- measured_pheno_all[,c(1,i)]
    id <- colnames(measured_pheno_all)[i]
    write.table(tmp, paste0('05_validation_using_INTERVAL/02_INTERVAL_protein/', id, '.pheno'), sep = '\t', quote = F, row.names = F)
}



# models in non-Hispanic White with Rsq >= 0.01
models <- fread('02_establish_models/EA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
models$ID <- gsub('.wgt.RDat', '', models$model)

intersect_in_MESA <- filter(intersect, SomaId %in% models$ID)
# select only unique pairs
intersect_in_MESA <- filter(intersect_in_MESA, !SomaId %in% unique(intersect_in_MESA[duplicated(intersect_in_MESA$SomaId),]$SomaId))


# Initialize parallel processing
num_cores <- parallelly::availableCores() 
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# generate predicted protein levels and calcuate Rsq
if (!dir.exists('05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level')){
    dir.create('05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level')
}

Rsq_res <- foreach(i = 1:nrow(intersect_in_MESA),.packages = c("data.table", "dplyr"), .combine = rbind)%dopar%{
    ID_MESA <- intersect_in_MESA[i,]$SomaId
    ID_INTERVAL <- intersect_in_MESA[i,]$SOMAMER_ID
    system(paste0('/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript 00_analysis_code/bin/make_score.R 02_establish_models/EA/08_establish_prediction_models/models/', ID_MESA, '.wgt.RDat > 05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level/', ID_MESA ))
    system(paste0('plink --bfile /mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/05_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID --score 05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level/', ID_MESA, ' 1 2 4 --out 05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level/', ID_MESA, '_predicted_pheno'))
    predicted_pheno <- fread(paste0('05_validation_using_INTERVAL/03_INTERVAL_predicted_protein_level/', ID_MESA, '_predicted_pheno.profile'), data.table = F)
    predicted_pheno <- select(predicted_pheno, IID, SCORE)
    measured_pheno <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/05_validation_using_INTERVAL/02_INTERVAL_protein/', ID_INTERVAL, '.pheno'), data.table = F)
    tmp <- left_join(predicted_pheno, measured_pheno, by = c('IID' = 'Sample_Name'))
    colnames(tmp)[3] <- 'protein_abundance'
    fit <- lm(SCORE~protein_abundance, tmp)
    Rsq <- summary(fit)$adj.r.squared
    res_tmp <- c(ID_MESA, ID_INTERVAL, Rsq)
}

Rsq_res <- data.frame(Rsq_res)
colnames(Rsq_res) <- c('MESA_ID', 'INTEVAL_ID', 'External_Rsq')
Rsq_res$External_Rsq <- as.numeric(Rsq_res$External_Rsq)
paste0((round(nrow(filter(Rsq_res, External_Rsq >= 0.01))/nrow(Rsq_res), 4) * 100 ), '%') # "75.11%%"
stopCluster(cl)
# get internal Rsq
Internal_R2 <- select(models, ID, Rsq)
Rsq_res <- left_join(Rsq_res, Internal_R2, by = c('MESA_ID' = 'ID'))

colnames(Rsq_res)[4] <- 'Internal_Rsq'
write.table(Rsq_res, '05_validation_using_INTERVAL/validation_Rsq_External_Internal.txt', row.names = F, sep = '\t', quote = F)
