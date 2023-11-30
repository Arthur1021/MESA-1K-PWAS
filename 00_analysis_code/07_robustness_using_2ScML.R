library(data.table)
library(dplyr)


if(!dir.exists('07_robustness_using_2ScML')){
    dir.create('07_robustness_using_2ScML')
}

convert_37 <- fread('/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh37_to_rsID.txt', data.table = F, header = F)

# prepare T2D MVP CA
GWAS_T2D_MVP_CA <- fread('/mnt/lvm_vol_1/sliu/work/database/87975/PhenoGenotypeFiles/RootStudyConsentSet_phs001672.MVP.v7.p1.c1.HMB-MDS/AnalysisFiles/set4/Submissions/sub20200305/MVP.T2D.ASN.MAF001.dbGaP.checked.txt.gz', data.table = F)
GWAS_T2D_MVP_CA$chr_pos <- paste0(GWAS_T2D_MVP_CA$CHR,':', GWAS_T2D_MVP_CA$POS)
GWAS_T2D_MVP_CA <- inner_join(GWAS_T2D_MVP_CA, convert_37, by = c('chr_pos' = 'V1'))
GWAS_T2D_MVP_CA <- select(GWAS_T2D_MVP_CA, V2, EA, NEA, BETA, SE, P, N)
colnames(GWAS_T2D_MVP_CA) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'WEIGHT')
fwrite(GWAS_T2D_MVP_CA, '07_robustness_using_2ScML/T2D_CA_MVP_input_for_meta.txt', sep = '\t')

# prepare T2D MVP HA
GWAS_T2D_MVP_HA <- fread('/mnt/lvm_vol_1/sliu/work/database/87975/PhenoGenotypeFiles/RootStudyConsentSet_phs001672.MVP.v7.p1.c1.HMB-MDS/AnalysisFiles/set4/Submissions/sub20200305/MVP.T2D.AMR.MAF001.dbGaP.checked.txt.gz', data.table = F)
GWAS_T2D_MVP_HA$chr_pos <- paste0(GWAS_T2D_MVP_HA$CHR,':', GWAS_T2D_MVP_HA$POS)
GWAS_T2D_MVP_HA <- inner_join(GWAS_T2D_MVP_HA, convert_37, by = c('chr_pos' = 'V1'))
GWAS_T2D_MVP_HA <- select(GWAS_T2D_MVP_HA, V2, EA, NEA, BETA, SE, P, N)
colnames(GWAS_T2D_MVP_HA) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'WEIGHT')
fwrite(GWAS_T2D_MVP_HA, '07_robustness_using_2ScML/T2D_HA_MVP_input_for_meta.txt', sep = '\t')


# prepare T2D PAGE CA
GWAS_T2D_PAGE <- fread('/home/sliu/work/database/summary_stat/raw/metal_T2D_ME_with_EUR_AA_JA_LA_NH_metaColumns.txt', data.table = F)
GWAS_T2D_PAGE <- inner_join(GWAS_T2D_PAGE, convert_37, by = c('MarkerName' = 'V1'))

GWAS_T2D_PAGE_CA <- select(GWAS_T2D_PAGE, V2, metal_T2D_Asian_3Files_1.txt.Allele1, metal_T2D_Asian_3Files_1.txt.Allele2, metal_T2D_Asian_3Files_1.txt.Effect, metal_T2D_Asian_3Files_1.txt.StdErr, 'metal_T2D_Asian_3Files_1.txt.P-value')
GWAS_T2D_PAGE_CA$WEIGHT <- 7437
colnames(GWAS_T2D_PAGE_CA) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'WEIGHT')
GWAS_T2D_PAGE_CA <- filter(GWAS_T2D_PAGE_CA, !GWAS_T2D_PAGE_CA$P == 'NULL')
GWAS_T2D_PAGE_CA$A1 <- toupper(GWAS_T2D_PAGE_CA$A1)
GWAS_T2D_PAGE_CA$A2 <- toupper(GWAS_T2D_PAGE_CA$A2)
fwrite(GWAS_T2D_PAGE_CA, '07_robustness_using_2ScML/T2D_CA_PAGE_input_for_meta.txt', sep = '\t')

# prepare T2D PAGE HA
GWAS_T2D_PAGE_HA <- select(GWAS_T2D_PAGE, V2, metal_T2D_Hispanic_6Files_1.txt.Allele1, metal_T2D_Hispanic_6Files_1.txt.Allele2, metal_T2D_Hispanic_6Files_1.txt.Effect, metal_T2D_Hispanic_6Files_1.txt.StdErr, 'metal_T2D_Hispanic_6Files_1.txt.P-value')
GWAS_T2D_PAGE_HA$WEIGHT <- 32871
colnames(GWAS_T2D_PAGE_HA) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'WEIGHT')
GWAS_T2D_PAGE_HA <- filter(GWAS_T2D_PAGE_HA, !GWAS_T2D_PAGE_HA$P == 'NULL')
GWAS_T2D_PAGE_HA$A1 <- toupper(GWAS_T2D_PAGE_HA$A1)
GWAS_T2D_PAGE_HA$A2 <- toupper(GWAS_T2D_PAGE_HA$A2)
fwrite(GWAS_T2D_PAGE_HA, '07_robustness_using_2ScML/T2D_HA_PAGE_input_for_meta.txt', sep = '\t')

# meta analysis for CA
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 07_robustness_using_2ScML/metal_CA_MVP_and_PAGE.txt')

# meta analysis for HA
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 07_robustness_using_2ScML/metal_HA_MVP_and_PAGE.txt')



# 2ScML analysis
library(data.table)
library(dplyr)
library(BEDMatrix)
library(stringr)
library(TScML)
library(stats)

# generate input file for CA
gwas_mvp_page_meta <- fread('07_robustness_using_2ScML/meta_T2D_GWAS_CA_MVP_PAGE1.tbl', data.table = F)
gwas_mvp_page_meta <- select(gwas_mvp_page_meta, MarkerName, Allele1, Allele2, Effect, StdErr, 'P-value')
colnames(gwas_mvp_page_meta) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
gwas_mvp_page_meta$N <- 223724
fwrite(gwas_mvp_page_meta, '07_robustness_using_2ScML/T2D_CA_GWAS_MVP_PAGE_meta_summary.txt' , sep = '\t')

# generate input file for HA
gwas_mvp_page_meta <- fread('07_robustness_using_2ScML/meta_T2D_GWAS_HA_MVP_PAGE1.tbl', data.table = F)
gwas_mvp_page_meta <- select(gwas_mvp_page_meta, MarkerName, Allele1, Allele2, Effect, StdErr, 'P-value')
colnames(gwas_mvp_page_meta) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P')
gwas_mvp_page_meta$N <- 53316
fwrite(gwas_mvp_page_meta, '07_robustness_using_2ScML/T2D_HA_GWAS_MVP_PAGE_meta_summary.txt' , sep = '\t')

# generate input file for EA
gwas_mvp_page_meta <- fread('/mnt/lvm_vol_1/sliu/work/database/87975/PhenoGenotypeFiles/RootStudyConsentSet_phs001672.MVP.v7.p1.c1.HMB-MDS/AnalysisFiles/set4/Submissions/sub20200305/MVP.T2D.EUR.MAF001.dbGaP.checked.txt.gz', data.table = F)
gwas_mvp_page_meta$chr_pos <- paste0(gwas_mvp_page_meta$CHR,':', gwas_mvp_page_meta$POS)
gwas_mvp_page_meta <- inner_join(gwas_mvp_page_meta, convert_37, by = c('chr_pos' = 'V1'))
gwas_mvp_page_meta <- select(gwas_mvp_page_meta, V2, EA, NEA, BETA, SE, P, N)
colnames(gwas_mvp_page_meta) <- c('SNP', 'A1', 'A2', 'BETA', 'SE', 'P', 'N')
fwrite(gwas_mvp_page_meta, '07_robustness_using_2ScML/T2D_EA_summary.txt' , sep = '\t')


# PatchUp
PatchUp <- function(M) {
    M <- apply(M, 2, function(x) {
        x[is.na(x)] <- mean(x, na.rm = TRUE)
        return(x)
    })

    return(M)
}

#allele qc
allele.qc = function(a1,a2,ref1,ref2) {
    a1 = toupper(a1)
    a2 = toupper(a2)
    ref1 = toupper(ref1)
    ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"

	flip1 = flip
	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
	return(snp)
}


CA_sig <- fread('04_meta_analysis/T2D_MVP_PAGE_CA_meta1.tbl', data.table = F)
CA_sig$FDR <- p.adjust(CA_sig$'P-value', method = 'fdr')
CA_sig <- filter(CA_sig, FDR < 0.05)
n1 <- 69   	#sample size used to build models
n2 <- 223724		#sample size of GWAS summary
n.ref <- 504
sumstat <- fread('07_robustness_using_2ScML/T2D_CA_GWAS_MVP_PAGE_meta_summary.txt', data.table = F)


perform_2ScML <- function(sig_protein, n1, n2, n.ref, sumstat, race){
    out <- data.frame(
        Protein = character(),
        beta = numeric(),
        se = numeric(),
        Z = numeric(),
        pval = numeric(),
        stringsAsFactors = FALSE
    )

    for (protein in sig_protein){
        res <- rep(NA, 5)
        res[1] <- protein
        weight <- load(paste0('02_establish_models/', race, '/08_establish_prediction_models/models/', protein, '.wgt.RDat'))
        

        # #count blup number, if blup is the best model, while number of SNP used in model is bigger than 500, change the p value to 1, we will not use this model
        # count_blup_SNP <- sum(!is.na(wgt.matrix[,'blup']))
        # if (count_blup_SNP > 500){
        #     cv.performance['pval', 'blup'] <- 1
        # }
        #change the p-value of top1 to 1, we will not use this model
        cv.performance['pval', 'top1'] <- 1
        #change the p-value of blup to 1, we will not use this model
        cv.performance['pval', 'blup'] <- 1
        #change the p-value of lasso to 1, we will not use this model
        cv.performance['pval', 'lasso'] <- 1

        # which rows have rsq
        row.rsq = grep( "rsq" , rownames(cv.performance) )
        # which rows have p-values
        row.pval = grep( "pval" , rownames(cv.performance) )	
        # Identify the best model
        mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))

        best = which.min(cv.performance[2,])

        if ( names(best) == "lasso" || names(best) == "enet" ) {
        keep = wgt.matrix[,best] != 0
        } else if ( names(best) == "top1" ) {
        keep = which.max(wgt.matrix[,best]^2)
        } else { 
        keep = 1:nrow(wgt.matrix)
        }

        #get SNP rsID used in models
        SNPS_used_in_model <- (snps[,c(2,5,6)])[keep,]$V2

        #get SNP in reference panel
        Z.ref.bim <- fread(paste0('02_establish_models/', race, '/08_establish_prediction_models/ref_ld/', protein, '_1000G.bim'))
        SNPS_in_ref <- Z.ref.bim$V2

        #get SNP GWAS summary
        SNPS_in_GWAS_summary <- sumstat$SNP

        #common SNPs in model, reference panel, and GWAS summary
        common_SNPs <- intersect(intersect(SNPS_used_in_model, SNPS_in_ref), SNPS_in_GWAS_summary)

        #filter 

        snps <- filter(snps, snps$V2 %in% common_SNPs)
        Z.ref.bim <- filter(Z.ref.bim, Z.ref.bim$V2 %in% common_SNPs)


        pheno <- fread(paste0('02_establish_models/', race, '/05_pheno_adjustment/', protein, '.pheno'), data.table = F)


        #reference panel cor, non-scale same as scale, here use non-scale
        Z.ref <- BEDMatrix(paste0('02_establish_models/', race, '/08_establish_prediction_models/ref_ld/', protein, '_1000G'))
        Z.ref  <- PatchUp(Z.ref)
        Z.ref<- as.matrix(Z.ref)
        Z.ref <- as.data.frame(Z.ref)
        ref_SNP_ID <- str_split(colnames(Z.ref),'_') %>% unlist() %>% matrix(., ncol = 2, byrow = TRUE)
        colnames(Z.ref) <- ref_SNP_ID[,1]
        Z.ref <- Z.ref [,common_SNPs]
        Z.ref.original <- Z.ref
        cor.Z.ref.original = cor(Z.ref.original) + diag(0.00001,ncol(Z.ref.original))



        #calculate cor.D1Z1.original
        
        Z.Stage1 = BEDMatrix(paste0('02_establish_models/', race, '/07_extract_potential_SNP_predictors/', protein, '_filtered'))
        Z.Stage1   <- PatchUp(Z.Stage1)
        Z.Stage1 <- as.matrix(Z.Stage1)
        Z.Stage1 <- as.data.frame(Z.Stage1)
        tmp_SNP_ID <- str_split(colnames(Z.Stage1),'_') %>% unlist() %>% matrix(., ncol = 2, byrow = TRUE)
        colnames(Z.Stage1) <- tmp_SNP_ID[,1]
        Z.Stage1 <- Z.Stage1[,common_SNPs]
        Z <- Z.Stage1 + diag(0.00001,ncol(Z.Stage1))
        Z1 = scale(Z,scale = F)
        D <- pheno[,protein]
        D1 = scale(D,scale = F)
        cor.D1Z1.original = as.numeric(cor(D1,Z1))



        ### Stage1 with individual-level data
        Stage1FittedModel = 
        TScMLStage1(cor.D1Z1 = cor.D1Z1.original,
                    Cap.Sigma.stage1 = cor(Z1),
                    n1 = n1,
                    p = nrow(Z.Stage1),
                    ind.stage1 = 1:ncol(Z1))
        #flip 

        qc1 = allele.qc( snps$V5 , snps$V6 , Z.ref.bim$V5 , Z.ref.bim$V6 )
        Stage1FittedModel[ qc1$flip ] = -1 * Stage1FittedModel[qc1$flip]

        # Match summary data to input, drop NA where summary data is missing
        m = match( Z.ref.bim$V2 , sumstat$SNP )
        sumstat_select = sumstat[m,]

        # QC / allele-flip the input and output
        qc = allele.qc( sumstat_select$A1 , sumstat_select$A2 , Z.ref.bim$V5 , Z.ref.bim$V6 )

        # Flip BETA-scores for mismatching alleles
        sumstat_select$BETA[ qc$flip ] = -1 * sumstat_select$BETA[ qc$flip ]
        sumstat_select$A1[ qc$flip ] = Z.ref.bim[ qc$flip]$V5
        sumstat_select$A2[ qc$flip ] = Z.ref.bim[ qc$flip]$V6
        cor.Y2Z2.original <- sumstat_select$BETA/sqrt(sumstat_select$BETA^2 + ( sumstat_select$N-2)*sumstat_select$SE^2)

        ### TScML Stage2 with summary data and reference panel
        lm1 = summary(lm(D1~Z1))
        Est.Sigma1Square = (lm1$sigma)^2
        if (!is.na(Est.Sigma1Square)){
            Est.Sigma2Square = as.numeric(1 - cor.Y2Z2.original%*%solve(cor.Z.ref.original,tol=0)%*%cor.Y2Z2.original)
            if(Est.Sigma1Square<=0)
            {
            Est.Sigma1Square = 1
            }
            if(Est.Sigma2Square<=0)
            {
            Est.Sigma2Square = 1
            }
            number <- floor(ncol(Z.Stage1)/2)
            TScMLStage2.Ref =
                TScMLStage2(gamma.hat.stage1 = Stage1FittedModel,
                            cor.Y2Z2 = cor.Y2Z2.original,
                            Estimated.Sigma = cor.Z.ref.original,
                            n1 = n1,
                            n2 = n2,
                            p = ncol(Z.Stage1),   # may be different
                            K.vec.stage2 = 0:number,
                            Est.Sigma1Square = Est.Sigma1Square,
                            Est.Sigma2Square = Est.Sigma2Square)

            TScML.Summary.Var = 
                TScMLVar(Z.ref.original = Z.ref.original,
                        Stage1FittedModel = Stage1FittedModel,
                        betaalpha.hat.stage2 = TScMLStage2.Ref$betaalpha.hat.stage2,
                        Est.Sigma1Square = Est.Sigma1Square,
                        Est.Sigma2Square = Est.Sigma2Square,
                        n1 = n1,
                        n2 = n2,
                        n.ref = n.ref)
            res.beta <- TScMLStage2.Ref$betaalpha.hat.stage2[1] 
            res.se <- sqrt(TScML.Summary.Var )
            res.Z <- res.beta/res.se
            res.pval <- pnorm(-abs(res.beta/res.se))*2
            res[2] <- res.beta
            res[3] <- res.se
            res[4] <- res.Z
            res[5] <- res.pval
            print (res)
        }
        out[nrow(out) + 1, ] <- res
	
    }
    return(out)
}


# CA
CA_sig <- fread('04_meta_analysis/T2D_MVP_PAGE_CA_meta1.tbl', data.table = F)
CA_sig$FDR <- p.adjust(CA_sig$'P-value', method = 'fdr')
CA_sig <- filter(CA_sig, FDR < 0.05)
n1 <- 69   	#sample size used to build models
n2 <- 223724		#sample size of GWAS summary
n.ref <- 504
sumstat <- fread('07_robustness_using_2ScML/T2D_CA_GWAS_MVP_PAGE_meta_summary.txt', data.table = F)
race <- 'CA'
res_CA <- perform_2ScML(CA_sig$MarkerName, n1, n2, n.ref, sumstat, race)
res_CA$race <- 'CA'
res_CA <- inner_join(CA_sig, res_CA, by = c('MarkerName' = 'Protein'))
res_CA <- select(res_CA, MarkerName, race, Zscore, 'P-value', FDR, beta, se, Z, pval)
colnames(res_CA) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')

# HA
HA_sig <- fread('04_meta_analysis/T2D_MVP_PAGE_HA_meta1.tbl', data.table = F)
HA_sig$FDR <- p.adjust(HA_sig$'P-value', method = 'fdr')
HA_sig <- filter(HA_sig, FDR < 0.05)
n1 <- 284   	#sample size used to build models
n2 <- 53316		#sample size of GWAS summary
n.ref <- 347
sumstat <- fread('07_robustness_using_2ScML/T2D_HA_GWAS_MVP_PAGE_meta_summary.txt', data.table = F)
race <- 'HA'
res_HA <- perform_2ScML(HA_sig$MarkerName, n1, n2, n.ref, sumstat, race)
res_HA$race <- 'HA'
res_HA <- inner_join(HA_sig, res_HA, by = c('MarkerName' = 'Protein'))
res_HA <- select(res_HA, MarkerName, race, Zscore, 'P-value', FDR, beta, se, Z, pval)
colnames(res_HA) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')



# EA
EA_sig <- fread('03_association/association_T2D_MVP_EA/all_association.out', data.table = F)
EA_sig <- filter(EA_sig, FDR < 0.05)
n1 <- 409   	#sample size used to build models
n2 <- 1114458		#sample size of GWAS summary
n.ref <- 502
sumstat <- fread('07_robustness_using_2ScML/T2D_EA_summary.txt', data.table = F)
race <- 'EA'
res_EA <- perform_2ScML(EA_sig$ID, n1, n2, n.ref, sumstat, race)
res_EA$race <- 'EA'
res_EA <- inner_join(EA_sig, res_EA, by = c('ID' = 'Protein'))
res_EA <- select(res_EA, ID, race, PWAS.Z, PWAS.P, FDR, beta, se, Z, pval)
colnames(res_EA) <- c('ID', 'race', 'Z', 'P', 'FDR', 'beta', 'se', 'z', 'p')



res_all <- rbind(res_CA, res_EA, res_HA)
res_all$p <- as.numeric(res_all$p)

res_all <- res_all[order(res_all$p), ]
res_all$fdr <- p.adjust(res_all$p, method = 'fdr')
res_all$Z <- as.numeric(res_all$Z)
res_all$z <- as.numeric(res_all$z)
res_all$Consistent_direction <- ifelse(sign(res_all$Z) == sign(res_all$z), 'TRUE', 'FALSE')
write.table(res_all, '07_robustness_using_2ScML/robustness_analysis_using_2ScML.txt', row.names = F, sep = '\t', quote = F)
