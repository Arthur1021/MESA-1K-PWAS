###################################
### 1. generate annotation file ###
###################################

library(readr)
library(biomaRt)
library(readxl)
library(dplyr)
library(stringr)

annotation <- read_excel('00_rawdata/phs001416_TOPMed_Proteomics_MESA/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx', n_max = 6)
annotation <- data.frame(annotation, check.names = F)
row.names(annotation) <- annotation$SeqId
annotation <- data.frame(t(annotation))
annotation <- annotation[-1, ]
annotation$UniProt <- gsub('  ', ' ', annotation$UniProt) %>% gsub(', ', ',', .) %>% gsub(' ', ',', .) 

annotation <- annotation %>% filter(str_starts(SomaId, "SL"))
annotation <- annotation %>% filter(., EntrezGeneSymbol != 'Human-virus')
nrow(annotation) # 1301
proteomics <- read_excel('00_rawdata/phs001416_TOPMed_Proteomics_MESA/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx', skip = 7)



sum(str_detect(annotation$UniProt, "[,]")) #46 -- mapped to multiple targets
# mapped to single targets
annotation_single <- annotation[!str_detect(annotation$UniProt, "[,]"),]
nrow(annotation_single) # 1255 -- mapped to single targets


annotation_single <- dplyr::select(annotation_single, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol)
nrow(annotation_single) # 1255

#######################################################################################
### first search by uniprot_gn_id in the ensembl database and matched entrezgene_id ###
#######################################################################################
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", version=110)

allgene <- unique(annotation_single$UniProt)
length(allgene) # 1225 -- unique UniProt

# extract from ensembl database based on uniprot_gn_id
annotation_try1 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "uniprot_gn_id",
    values = allgene,
    mart = ensembl)
nrow(annotation_try1) # 1472
annotation_try1 <- unique(annotation_try1)
nrow(annotation_try1) # 1472
# Update TSS column based on the condition
annotation_try1$TSS <- ifelse(annotation_try1$strand == 1, annotation_try1$start_position, annotation_try1$end_position)

# ann_Feature_Data_single <- dplyr::select(ann_Feature_Data_single, uniprot_gn_id, hgnc_symbol, chromosome_name, TSS)
annotation_try1 <- dplyr::inner_join(annotation_single, annotation_try1, by = c('UniProt' = 'uniprot_gn_id'))
annotation_try1 <- unique(annotation_try1)
annotation_try1 <- annotation_try1[annotation_try1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
nrow(annotation_try1) # 1249
length(unique(annotation_try1$SomaId)) # found 1224 SomaId in autosome and sex chromosomes

correct0 <- annotation_try1[!annotation_try1$SomaId %in% annotation_try1[duplicated(annotation_try1$SomaId),]$SomaId,]
nrow(correct0) # 1211
duplicated_rows_df <- annotation_try1[annotation_try1$SomaId %in% annotation_try1[duplicated(annotation_try1$SomaId),]$SomaId,]


correct1 <- duplicated_rows_df[duplicated_rows_df$EntrezGeneID == duplicated_rows_df$entrezgene_id & duplicated_rows_df$EntrezGeneSymbol == duplicated_rows_df$hgnc_symbol,]
correct1 <- correct1[!duplicated(correct1$SomaId),] # X, Y duplicate, force remain first row
nrow(correct1) # 11

unique(duplicated_rows_df[!duplicated_rows_df$SomaId %in% correct1$SomaId,]$SomaId)
# [1] "SL008158" "SL004356"
# cannot matched by both EntrezGeneID and EntrezGeneSymbol
# manually select
correct2 <- duplicated_rows_df[!duplicated_rows_df$SomaId %in% correct1$SomaId,][c(1, 12),]
nrow(correct2) # 2

annotation_res1 <- rbind(correct0, correct1, correct2)

##############################################################
### second search by entrezgene_id in the ensembl database ###
##############################################################

not_found_in_annotation_try1 <- annotation_single[!annotation_single$SomaId %in% annotation_res1$SomaId,]
nrow(not_found_in_annotation_try1) # 31 SomaId not found

not_found_in_annotation_try1 <- dplyr::filter(not_found_in_annotation_try1, EntrezGeneID != '')
nrow(not_found_in_annotation_try1) # 31 SomaId not found

allgene <- unique(not_found_in_annotation_try1$EntrezGeneID)
annotation_try2 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "entrezgene_id",
    values = allgene,
    mart = ensembl)

nrow(annotation_try2) #1438
annotation_try2 <- unique(annotation_try2)
nrow(annotation_try2) #1438
annotation_try2 <- annotation_try2[annotation_try2$chromosome_name %in% c(as.character(1:22), "X", "Y"), ]
nrow(annotation_try2) #153
# Update TSS column based on the condition
annotation_try2$TSS <- ifelse(annotation_try2$strand == 1, annotation_try2$start_position, annotation_try2$end_position)
annotation_try2$entrezgene_id <- as.character(annotation_try2$entrezgene_id)
annotation_try2 <- inner_join(not_found_in_annotation_try1, annotation_try2, by = c('EntrezGeneID' = 'entrezgene_id'))
nrow(annotation_try2) # 153

correct3 <- annotation_try2[!annotation_try2$SomaId %in% annotation_try2[duplicated(annotation_try2$SomaId),]$SomaId,]
nrow(correct3) # 5
duplicated_rows_df <- annotation_try2[annotation_try2$SomaId %in% annotation_try2[duplicated(annotation_try2$SomaId),]$SomaId,]
duplicated_rows_df$uniprot_gn_id <- NA
duplicated_rows_df <- unique(duplicated_rows_df)
correct4 <- duplicated_rows_df
length(unique(correct4$SomaId)) # 24
annotation_res2 <- rbind(correct3, correct4)
nrow(annotation_res2) # 29
###########################################################
### third search by hgnc_symbol in the ensembl database ###
###########################################################
not_found_in_annotation_try2 <- annotation_single[!annotation_single$SomaId %in% c(annotation_res1$SomaId, annotation_res2$SomaId),]
nrow(not_found_in_annotation_try2) # 2
# > not_found_in_annotation_try2
#             SomaId                                TargetFullName Target UniProt
# 4916-2_1  SL000460                              Immunoglobulin D    IgD  P01880
# 5097-14_3 SL014289 Killer cell immunoglobulin-like receptor 3DS1  KI3S1  Q14943
#              EntrezGeneID EntrezGeneSymbol
# 4916-2_1  3495 50802 3535  IGHD  IGK@ IGL@
# 5097-14_3            3813          KIR3DS1

# can not find



annotation_res1 <- select(annotation_res1, SomaId, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS)
annotation_res2 <- select(annotation_res2, SomaId, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS)

annotation_single_all <- rbind(annotation_res1, annotation_res2)
nrow(annotation_single_all) # 1253
# 1255 - 1253 = 10. 2 are in Scaffold are could not be found 

# mapped to multiple targets
library(tidyverse)
annotation_multiple <- annotation[str_detect(annotation$UniProt, "[,]"),]
nrow(annotation_multiple) # 46 -- mapped to multiple targets

annotation_multiple <- dplyr::select(annotation_multiple, SomaId, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol)
nrow(annotation_multiple) # 46

annotation_multiple <- data.frame(separate_rows(annotation_multiple, UniProt, sep = "[,]"))
allgene <- unique(annotation_multiple$UniProt)
length(allgene) # 92 -- unique UniProt

# extract from ensembl database based on uniprot_gn_id
annotation_multiple_try1 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "uniprot_gn_id",
    values = allgene,
    mart = ensembl)

annotation_multiple_try1 <- annotation_multiple_try1[annotation_multiple_try1$chromosome_name %in% c(as.character(1:22), "X", "Y"), ] 
# manually remove duplicate
annotation_multiple_try1 <- annotation_multiple_try1[c( -38, -76),]
nrow(annotation_multiple_try1) # 89
# write.table(annotation_multiple_try1, 'annotation_multiple_try1.txt', row.names = F, sep  = '\t', quote = F)


annotation_multiple[!annotation_multiple$UniProt %in% annotation_multiple_try1$uniprot_gn_id,]

annotation_multiple_try2 <- getBM(attributes = c("uniprot_gn_id", "entrezgene_id", "hgnc_symbol","chromosome_name",
                    "start_position", "end_position","strand"),
    filters = "hgnc_symbol",
    values = c('C8B', 'C1QB', 'CGB'),
    mart = ensembl)
# can not find CGB
annotation_multiple_try2$uniprot_gn_id <- NA
annotation_multiple_try2 <- unique(annotation_multiple_try2)
annotation_multiple_try2$uniprot_gn_id <- c('P02746', 'P07358')

annotation_multiple_res <- rbind(annotation_multiple_try1 , annotation_multiple_try2)

annotation_multiple_all <- left_join(annotation_multiple, annotation_multiple_res, by = c('UniProt' = 'uniprot_gn_id'))
annotation_multiple_all$TSS <- ifelse(annotation_multiple_all$strand == 1, annotation_multiple_all$start_position, annotation_multiple_all$end_position)
annotation_multiple_all <- dplyr::select(annotation_multiple_all, SomaId, SomaId, TargetFullName, Target, UniProt, EntrezGeneID, EntrezGeneSymbol, chromosome_name, start_position, end_position, strand, TSS )

# rbind single and multiple annotation result
annotation_all <- rbind(annotation_single_all, annotation_multiple_all)
length(unique(annotation_all$SomaId)) # 1299

if (!file.exists("01_preprocess_data/annotation")) {
  dir.create("01_preprocess_data/annotation", recursive = TRUE)
}

write.table(annotation_all, '01_preprocess_data/annotation/SomaScan_annotation.txt', sep = '\t', quote = F, row.names =F)

annotation_TSS <- select(annotation_all, SomaId, chromosome_name, TSS, start_position, end_position, strand)
colnames(annotation_TSS) <- c('SomaId', 'Chr', 'TSS', 'Start', 'End', 'Strand')
annotation_TSS <- annotation_TSS[annotation_TSS$Chr %in% c(as.character(1:22)), ]
write.table(annotation_TSS, '01_preprocess_data/annotation/SomaScan_TSS.txt', sep = '\t', quote = F, row.names =F)


#############################################################
### 2. get overlapped subjects of covariate and phenotype ###
#############################################################

library(data.table)
library(dplyr)
library(readxl)

covariate <- fread('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESAe1and5_CoVarLangWu_20210630/MESAe1and5_CoVarLangWu_20210630.txt', data.table = F)
subject_with_covariate <- covariate$sidno

# protein expression
proteomics <- read_excel('00_rawdata/phs001416_TOPMed_Proteomics_MESA/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx', skip = 7)
proteomics <- filter(proteomics, RowCheck == 'PASS') 
proteomics <- data.frame(proteomics)
# VISIT 1
SampleAttributesDS <- fread('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/Omics_Proteomics_SampleAttributes/MESA_phs001416_Proteomics_SampleAttributesDS_20200821.txt', data.table = F)
SampleAttributesDS <- filter(SampleAttributesDS, COLLECTION_VISIT == 1)

# select VISIT 1
proteomics <- filter(proteomics, TOP_ID %in% SampleAttributesDS$PRIMARY_BIOSAMPLE_ID)
phenotype <- proteomics[, c(5, 33:ncol(proteomics))]


subjects_with_SomaScan <- select(proteomics, sidno, collaboratorsampleid)
subjects_with_SomaScan_and_covariate <- filter(subjects_with_SomaScan, sidno %in% subject_with_covariate)
subjects_with_genotype <- fread('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/minDP10/freeze.7a.chr1.pass_and_fail.gtonly.minDP10.fam', data.table = F)
subjects_with_SomaScan_and_covariate_and_genotype <- filter(subjects_with_SomaScan_and_covariate, collaboratorsampleid %in% subjects_with_genotype$V2)
nrow(subjects_with_SomaScan_and_covariate_and_genotype) # 976

sample_AA <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/AA.txt'), data.table = F)
sample_CA <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/CA.txt'), data.table = F)
sample_EA <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/EA.txt'), data.table = F)
sample_HA <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/HA.txt'), data.table = F)

subjects_with_SomaScan_and_covariate_and_genotype <- filter(subjects_with_SomaScan_and_covariate_and_genotype, sidno %in% c(sample_AA$V1, sample_CA$V1, sample_EA$V1, sample_HA$V1))
write.table(subjects_with_SomaScan_and_covariate_and_genotype, '01_preprocess_data/subjects_with_SomaScan_and_covariate_and_genotype_sidno_NWDID.txt', row.names = F, sep = '\t', quote = F)
write.table(subjects_with_SomaScan_and_covariate_and_genotype$collaboratorsampleid, '01_preprocess_data/subjects_with_SomaScan_and_covariate_and_genotype.txt', row.names = F, sep = '\t', quote = F, col.names = F)

#################################
### 3. preprare genotype data ###
#################################
library(data.table)
library(dplyr)
library(stringr)
library(foreach)
library(doParallel)
library(parallel)

if (!dir.exists(file.path('01_preprocess_data/wgs'))) {
    dir.create(file.path('01_preprocess_data/wgs'))
}

files <- Sys.glob(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/minDP10/freeze.7a.chr*.pass_and_fail.gtonly.minDP10.vcf'))
files <- files[1:22]
# Initialize parallel processing
num_cores <- 22
cl <- makeCluster(num_cores)
registerDoParallel(cl)
    
# for (file in files){
foreach(file = files, .packages = c("data.table", "dplyr", "readr", "stringr"), .combine = rbind) %dopar% { 
    ## extract SNP
    ## SNP only, rename SNP ID as Chr:pos, SNP 'PASS', sample
    outfile_name <- paste0('01_preprocess_data/wgs/',gsub('freeze.7a.', '', gsub('.pass_and_fail.gtonly.minDP10.vcf', '', last(str_split(file, '/')[[1]]))), '_SNP')
    system2(
        command = "plink2",
        args = c(
        "--vcf", file,
        "--var-filter",
        "--snps-only",
        "--set-all-var-ids", "\"@:#\"",
        "--keep", "01_preprocess_data/subjects_with_SomaScan_and_covariate_and_genotype.txt",
        "--make-bed",
        "--out", outfile_name
        )
    )
    ## remove duplicate SNPs
    outfile_name2 <- paste0(outfile_name, '_unique')
    system2(
        command = "plink2",
        args = c(
        "--bfile", outfile_name,
        "--rm-dup", "exclude-all",
        "--make-bed",
        "--out", outfile_name2
        )
    )
}
stopCluster(cl)


#############################
### 3. extract covariates ###
#############################


library(data.table)
library(dplyr)
if (!dir.exists(file.path('01_preprocess_data/covariates'))) {
    dir.create(file.path('01_preprocess_data/covariates'))
}

covariate <- fread('00_rawdata/covariates/MESAe1and5_CoVarLangWu_20210630/MESAe1and5_CoVarLangWu_20210630.txt', data.table = F)
covariate <- select(covariate, sidno, gender1, age1c, htcm1, wtlb1, cig1c)
covariate$bmi <- covariate$wtlb1 * 0.453592 / (covariate$htcm1 / 100)^2
covariate <- covariate[complete.cases(covariate), ]

subjects_with_SomaScan_and_covariate_and_genotype <- fread('01_preprocess_data/subjects_with_SomaScan_and_covariate_and_genotype_sidno_NWDID.txt', data.table = F)
# male and female
races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){
    sample <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/', race, '.txt'), data.table = F)
    covariate_tmp <- filter(covariate, sidno %in% sample$V1 & sidno %in% subjects_with_SomaScan_and_covariate_and_genotype$sidno)
    covariate_tmp <- left_join(covariate_tmp, subjects_with_SomaScan_and_covariate_and_genotype, by = c('sidno' = 'sidno'))
    covariate_tmp <- select(covariate_tmp, collaboratorsampleid, gender1, age1c, bmi, cig1c)
    colnames(covariate_tmp)[1] <- 'Sample_ID'
    write.table (covariate_tmp, paste0('01_preprocess_data/covariates/', race, '.txt'), row.names = F, sep = '\t', quote = F)
}

# male only

for (race in races){
    sample <- fread(paste0('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/covariates/MESA_SHARe_PCA/MESA_SHARe_PCA/', race, '.txt'), data.table = F)
    covariate_tmp <- filter(covariate, sidno %in% sample$V1 & sidno %in% subjects_with_SomaScan_and_covariate_and_genotype$sidno)
    covariate_tmp <- left_join(covariate_tmp, subjects_with_SomaScan_and_covariate_and_genotype, by = c('sidno' = 'sidno'))
    covariate_tmp <- select(covariate_tmp, collaboratorsampleid, gender1, age1c, bmi, cig1c)
    covariate_tmp <- filter(covariate_tmp, gender1 == 1)
    covariate_tmp <- select(covariate_tmp, collaboratorsampleid, age1c, bmi, cig1c)
    colnames(covariate_tmp)[1] <- 'Sample_ID'
    write.table (covariate_tmp, paste0('01_preprocess_data/covariates/', race, '_male.txt'), row.names = F, sep = '\t', quote = F)
}

############################
### 3. extract phenotype ###
############################
library(readxl)
library(dplyr)
library(data.table)
proteomics <- read_excel('00_rawdata/phs001416_TOPMed_Proteomics_MESA/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx', skip = 7)

# protein expression
proteomics <- read_excel('00_rawdata/phs001416_TOPMed_Proteomics_MESA/proteomics_data_merged_with_runlist_key_updated_mar1.xlsx', skip = 7)
proteomics <- filter(proteomics, RowCheck == 'PASS') 
proteomics <- data.frame(proteomics)
# VISIT 1
SampleAttributesDS <- fread('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_rawdata/Omics_Proteomics_SampleAttributes/MESA_phs001416_Proteomics_SampleAttributesDS_20200821.txt', data.table = F)
SampleAttributesDS <- filter(SampleAttributesDS, COLLECTION_VISIT == 1)


# select VISIT 1
proteomics <- filter(proteomics, TOP_ID %in% SampleAttributesDS$PRIMARY_BIOSAMPLE_ID)
# select proteins with TSS information
TSS <- fread('01_preprocess_data/annotation/SomaScan_TSS.txt', data.table = F)
proteomics <- proteomics[,colnames(proteomics) %in% c('collaboratorsampleid', TSS$SomaId)]
colnames(proteomics)[1] <- 'Sample_ID'

write.table(proteomics, '01_preprocess_data/proteomics.txt', row.names = F, sep = '\t', quote = F)
