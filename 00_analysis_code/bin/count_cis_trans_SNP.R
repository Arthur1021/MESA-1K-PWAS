
library(stringr)
library(dplyr)

TSS <- read.table('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/01_preprocess_data/annotation/SomaScan_TSS.txt', header = T)
files <- Sys.glob('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/AFA_male/08_establish_prediction_models/models/*.wgt.RDat')

for (file in files){
    load(file)
    file <- '/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/AFA_male/08_establish_prediction_models/models/seq.16927.9.wgt.RDat'
    load('/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_7K_PWAS/AFA_male/08_establish_prediction_models/models/seq.16927.9.wgt.RDat')
    # which rows have rsq
	row.rsq = grep( "rsq" , rownames(cv.performance) )
	# which rows have p-values
	row.pval = grep( "pval" , rownames(cv.performance) )	
	
	# Identify the best model
	if ( !is.na(opt$force_model) ) {
		mod.best = which( colnames(wgt.matrix) == opt$force_model )
		if ( length(mod.best) == 0 ) {
			cat( "WARNING : --force_model" , mod.best ,"does not exist for", unlist(wgtlist[w,]) , "\n")
			cur.FAIL = TRUE
		}	
	} else {
		# get the most significant model
		mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
	}

    # Extract SNPs used in the most significant model
    SNP_used <- dplyr::select(data.frame(wgt.matrix), names(mod.best))
    SNP_used <- SNP_used[SNP_used[,1] != 0, ,drop = FALSE]
    SNP_used <- filter(snps, V2 %in% rownames(SNP_used))
    SNP_used <- dplyr::select(SNP_used, V2, V1, V4)
    colnames(SNP_used) <- c('SNP', 'Chr', 'Pos')
    SNP_used$cis_trans <- NA

    ID <- last(str_split(file, '/')[[1]]) %>% gsub('.wgt.RDat', '',.)
    TSS_tmp <- TSS[TSS[1] == ID,]
    TSS_tmp$upstream <- TSS_tmp$TSS - 1000000
    TSS_tmp$downstream <- TSS_tmp$TSS + 1000000
    TSS_tmp <- select(TSS_tmp, AptName, Chr, upstream, downstream)
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
    

}