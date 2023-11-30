library(data.table)
library(dplyr)
library(foreach)
library(doParallel)


# convert 37 to 38 chr:pos
SNP <- fread('05_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID.bim', data.table = F)
SNP$chr <- paste0('chr',SNP$V1)
SNP$start <- SNP$V4
SNP$end <- SNP$V4+1
SNP$name <- paste0(SNP$V1, ':', SNP$V4)
options(scipen = 10)
SNP$end <- as.numeric(SNP$end)
write.table(SNP[7:10], '06_comparison_with_models_from_Heather_group/SNP_chr_pos_37.bed', row.names = F, col.names = F, sep = '\t', quote = F)

#lift over to 38 using liftOver software
system('/mnt/lvm_vol_2/sliu/database/liftover/liftOver 06_comparison_with_models_from_Heather_group/SNP_chr_pos_37.bed /mnt/lvm_vol_2/sliu/database/liftover/hg19ToHg38.over.chain.gz 06_comparison_with_models_from_Heather_group/SNP_chr_pos_38.txt 06_comparison_with_models_from_Heather_group/SNP_chr_pos_38_unmapped.txt')



#update SNP ID to hg38
if(!dir.exists('06_comparison_with_models_from_Heather_group/01_prepare_geno')){
    dir.create('06_comparison_with_models_from_Heather_group/01_prepare_geno')
}

df <- fread('06_comparison_with_models_from_Heather_group/SNP_chr_pos_38.txt', data.table = F)
df$hg38 <- paste0(df$V1, ':', df$V2)
write.table(df[4:5], '06_comparison_with_models_from_Heather_group/01_prepare_geno/update_SNP_list.txt', row.names = F, col.names = F, quote = F, sep = '\t')
write.table(df[4], '06_comparison_with_models_from_Heather_group/01_prepare_geno/extract_SNP_list.txt', row.names = F, col.names = F, quote = F, sep = '\t')


system2(
    command = "plink2",
    args = c(
    "--bfile", "05_validation_using_INTERVAL/01_INTERVAL_genotype/chr1-22_rsID",
    "--set-all-var-ids", "\"@:#\"",
    "--extract", "06_comparison_with_models_from_Heather_group/01_prepare_geno/extract_SNP_list.txt",
    "--rm-dup", "exclude-all",
    "--make-bed",
    "--out", "06_comparison_with_models_from_Heather_group/01_prepare_geno/chr1-22_hg37"
    )
)

system2(
    command = "plink2",
    args = c(
    "--bfile", "06_comparison_with_models_from_Heather_group/01_prepare_geno/chr1-22_hg37",
    "--update-name", "06_comparison_with_models_from_Heather_group/01_prepare_geno/update_SNP_list.txt",
    "--make-bed",
    "--out", "06_comparison_with_models_from_Heather_group/01_prepare_geno/chr1-22_hg38"
    )
)

# split into chr
for (i in 1:22){
    system2(
        command = "plink",
        args = c(
        "--bfile", "06_comparison_with_models_from_Heather_group/01_prepare_geno/chr1-22_hg38",
        "--chr", i,
        "--output-chr", "chr26", # chr
        "--make-bed",
        "--out", paste0("06_comparison_with_models_from_Heather_group/01_prepare_geno/chr",i,'_hg38')
        )
    )
}

# convert to dosage
if(!dir.exists('06_comparison_with_models_from_Heather_group/02_dosage')){
    dir.create('06_comparison_with_models_from_Heather_group/02_dosage')
}

for(chr in 1:22){
    system(paste0('python2 06_comparison_with_models_from_Heather_group/convert_plink_to_dosage.py -b 06_comparison_with_models_from_Heather_group/01_prepare_geno/chr', chr, '_hg38 -p plink -o 06_comparison_with_models_from_Heather_group/02_dosage/chr'))
}

num_cores <- parallelly::availableCores()  
cl <- makeCluster(num_cores)
registerDoParallel(cl)

foreach (chr = 1:22, .packages = c("data.table"))%dopar%{
    df <- fread(paste0('06_comparison_with_models_from_Heather_group/02_dosage/chr', chr, '.txt.gz'), data.table = F)
    df$V1 <- paste0('chr',df$V1)
    write.table(df, file = gzfile (paste0('06_comparison_with_models_from_Heather_group/02_dosage/chr', chr, '.dosage.txt.gz')), row.names = F, quote = F, col.names = F)
}

stopCluster(cl)



if(!dir.exists('06_comparison_with_models_from_Heather_group/03_dosage_for_input')){
    dir.create('06_comparison_with_models_from_Heather_group/03_dosage_for_input')
}

system('mv 06_comparison_with_models_from_Heather_group/02_dosage/chr*dosage.txt.gz 06_comparison_with_models_from_Heather_group/03_dosage_for_input/')



#######################################
### baseline  and finemapped models ###
#######################################
# generate predicted for baseline models 
system('python2 06_comparison_with_models_from_Heather_group/predict_gene_expression.py --dosages 06_comparison_with_models_from_Heather_group/03_dosage_for_input/ --dosages_prefix chr --weights 06_comparison_with_models_from_Heather_group/models_from_Heather_group/EUR_PCAIR_baseline_models_rho0.1_zpval0.05.db --output 06_comparison_with_models_from_Heather_group/EUR_PCAIR_baseline_models_rho0.1_zpval0.05_predicted_exp')
# generate predicted for finemapped models
system('python2 06_comparison_with_models_from_Heather_group/predict_gene_expression.py --dosages 06_comparison_with_models_from_Heather_group/03_dosage_for_input/ --dosages_prefix chr --weights 06_comparison_with_models_from_Heather_group/models_from_Heather_group/EUR_PCAIR_dapg_0.001_T_rho0.1_zpval0.05.db --output 06_comparison_with_models_from_Heather_group/EUR_PCAIR_dapg_0.001_T_rho0.1_zpval0.05_predicted_exp')


models <- c('baseline_models', 'dapg_0.001_T')

for (model in models){
    predicted <- fread(paste0('06_comparison_with_models_from_Heather_group/EUR_PCAIR_', model, '_rho0.1_zpval0.05_predicted_exp'), data.table = F)
    colnames(predicted) <- sub("_.+", "", colnames(predicted))
    sample <- fread('06_comparison_with_models_from_Heather_group/01_prepare_geno/chr1_hg38.fam', data.table = F)
    sample$V1 <- as.character(sample$V1)
    predicted_level <- cbind(Sample = sample$V1, predicted)

    annotation_MESA <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
    # select only single target
    annotation_MESA <- annotation_MESA[!duplicated(annotation_MESA$SomaId),]
    annotation_MESA <- select(annotation_MESA, UniProt, SomaId)

    annotation_INTERVAL <- fread('/mnt/lvm_vol_2/sliu/project/INTERVAL_pwas/2_gene_start_end_37/001_SOMALOGIC_GWAS_protein_info.csv', data.table = F)
    annotation_INTERVAL <- select(annotation_INTERVAL, UniProt, SOMAMER_ID)

    intersect <- inner_join(annotation_MESA, annotation_INTERVAL, by = c('UniProt'))

    intersect_in_MESA <- filter(intersect, SomaId %in% colnames(predicted_level))
    # select only unique pairs
    intersect_in_MESA <- filter(intersect_in_MESA, !SomaId %in% unique(intersect_in_MESA[duplicated(intersect_in_MESA$SomaId),]$SomaId))


    Rsq_all <- data.frame(
        MESA_ID = as.character(),
        INTERVAL_ID = as.character(),
        Rsq = as.numeric()
    )

    for (i in 1:nrow(intersect_in_MESA)){
        predicted_level_tmp <- cbind(Sample = predicted_level[,1], predicted_level[intersect_in_MESA[i,]$SomaId])
        measured_level_tmp <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/05_validation_using_INTERVAL/02_INTERVAL_protein/', intersect_in_MESA[i,]$SOMAMER_ID, '.pheno'), data.table = F)
        measured_level_tmp$Sample_Name <- as.character(measured_level_tmp$Sample_Name)
        predicted_measured_tmp <- left_join(predicted_level_tmp, measured_level_tmp, by = c('Sample' = 'Sample_Name'))
        colnames(predicted_measured_tmp)[2:3] <- c('predicted', 'measured')
        predicted_measured_tmp$predicted <- as.numeric(predicted_measured_tmp$predicted)
        predicted_measured_tmp$measured <- as.numeric(predicted_measured_tmp$measured)
        fit <- lm(predicted~measured, predicted_measured_tmp)
        Rsq <- summary(fit)$adj.r.squared
        res_tmp <- c(intersect_in_MESA[i,]$SomaId, intersect_in_MESA[i,]$SOMAMER_ID, Rsq)
        Rsq_all[nrow(Rsq_all) + 1, ] <- res_tmp
    }
    Rsq_all$Rsq <- as.numeric(Rsq_all$Rsq)
    write.table(Rsq_all, paste0('06_comparison_with_models_from_Heather_group/validation_EA_models_using_INTERVAL_for_', model, '.txt'), row.names = F, sep = '\t', quote = F )
    print (paste0(round(nrow(filter(Rsq_all, Rsq >=0.01))/nrow(Rsq_all), 4)* 100 , '%')) # "70.32%" for baseline , "34.38%" for fine-mapped

}


# compare with cross validation Rsq in each ethnics(current and models from Heather group)

library("RSQLite")
library(data.table)
library(dplyr)
sqlite <- dbDriver("SQLite")

races <- c('AA', 'CA', 'EA', 'HA')

for (race in races){
    if (race == 'AA'){
        race_id <- 'AFA'
    }else if (race == 'CA') {
       race_id <- 'CHN'
    }else if (race == 'EA') {
       race_id <- 'EUR'
    }else if (race == 'HA') {
       race_id <- 'HIS'
    }
    dbname <- paste0("06_comparison_with_models_from_Heather_group/models_from_Heather_group/", race_id, "_PCAIR_baseline_models_rho0.1_zpval0.05.db")
    db = dbConnect(sqlite,dbname)
    query <- function(...) dbGetQuery(db, ...)
    df <- query('select * from extra')
    df <- df[c(2,17)]

    current_models <- fread(paste0('02_establish_models/', race, '/08_establish_prediction_models/model_Rsq_0.01.txt'), data.table = F)
    current_models$genename <- gsub('.wgt.RDat', '', current_models$model)
    current_models <- select(current_models, genename, Rsq)

    compare <- inner_join(df, current_models, by = 'genename')
    colnames(compare) <- c('ID', 'Heather_Rsq' , 'current_Rsq') 
    write.table(compare, paste0('06_comparison_with_models_from_Heather_group/Rsq_compare_current_with_baseline_from_heather_',race, '.txt'), row.names = F, sep = '\t', quote = F)
    print (sum(compare$current_Rsq > compare$Heather_Rsq)/nrow(compare))
}
# [1] 0.6987179 AA
# [1] 0.527027 CA
# [1] 0.6359447 EA
# [1] 0.4878049 HA

