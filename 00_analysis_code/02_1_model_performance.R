library(data.table)
library(dplyr)
library(stringr)

TSS_file <- '01_preprocess_data/annotation/SomaScan_TSS.txt'
races <- c('02_establish_models/AA', '02_establish_models/CA', '02_establish_models/EA', '02_establish_models/HA', '02_establish_models/AA_male', '02_establish_models/CA_male', '02_establish_models/EA_male', '02_establish_models/HA_male')
for (race in races){
    all_Rsq <- data.frame(
        model = as.character(),
        best_model = as.character(),
        Rsq = as.numeric(),
        N_SNP = as.numeric(),
        N_cis = as.numeric(),
        N_trans = as.numeric()  
    )
    files <- Sys.glob(file.path(race,'08_establish_prediction_models/models/*.RDat'))
    for (file in files){
        load(file)
        # which rows have rsq
        row.rsq = grep( "rsq" , rownames(cv.performance) )
        # which rows have p-values
        row.pval = grep( "pval" , rownames(cv.performance) )
        mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))

        if (cv.performance[row.rsq,names(mod.best)] < 0.01){
            next
        }

        # if this is a top1 model, clear out all the other weights
        if ( substr( (colnames(cv.performance))[ mod.best ],1,4) == "top1" ) wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0

        # count total number of SNPs, cis-SNP and trans-SNP used in models
        TSS <- read.table(TSS_file, header = T)
        SNP_used <- dplyr::select(data.frame(wgt.matrix), names(mod.best))
        SNP_used <- SNP_used[SNP_used[,1] != 0, ,drop = FALSE]
        if (nrow(SNP_used) == 0){
            next
        }

        SNP_used <- dplyr::filter(snps, snps$V2 %in% rownames(SNP_used))
        SNP_used <- dplyr::select(SNP_used, V2, V1, V4)
        colnames(SNP_used) <- c('SNP', 'Chr', 'Pos')
        SNP_used$cis_trans <- NA
        ID <- last(str_split(file, '/')[[1]]) %>% gsub('.wgt.RDat', '', .)
        TSS_tmp <- TSS[TSS[1] == ID,]
        TSS_tmp$upstream <- TSS_tmp$TSS - 1000000
        TSS_tmp$upstream[TSS_tmp$upstream < 0] <- 0 # fix TSS should be start with 0
        TSS_tmp$downstream <- TSS_tmp$TSS + 1000000
        TSS_tmp <- select(TSS_tmp, SomaId, Chr, upstream, downstream)
        for (m in 1:nrow(SNP_used)){
            tmp <- rep(NA, nrow(TSS_tmp))
            for (n in 1:nrow(TSS_tmp)){
                tmp[n] <- SNP_used[m,]$Chr == TSS_tmp[n,]$Chr & SNP_used[m,]$Pos >= TSS_tmp[n,]$upstream & SNP_used[m,]$Pos <= TSS_tmp[n,]$downstream
            }
            if (sum(tmp) >= 1){
                SNP_used[m,]$cis_trans <- 'cis'
            }else{
                SNP_used[m,]$cis_trans <- 'trans'
            }
        }
        N_SNP <- length(SNP_used$cis_trans)
        N_cis <- sum(SNP_used$cis_trans == 'cis')
        N_trans <- sum(SNP_used$cis_trans == 'trans')


        tmp <- rep(NA, 6)
        tmp[1] <- last(str_split(file, '/')[[1]])
        tmp[2] <- names(mod.best)
        tmp[3] <- cv.performance[row.rsq,mod.best]
        tmp[4] <- N_SNP
        tmp[5] <- N_cis
        tmp[6] <- N_trans
        all_Rsq[nrow(all_Rsq)+1,] <- tmp

    }
    outfile <- file.path(race, '08_establish_prediction_models/model_Rsq_0.01.txt')
    message(race, ' have ', nrow(all_Rsq), ' models with Rsq >= 0.01.')
    write.table(all_Rsq, outfile, row.names = F, sep = '\t', quote = F)
}