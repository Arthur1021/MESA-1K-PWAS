library(data.table)
library(stringr)
library(dplyr)
library(RNOmni)
library(foreach)
library(doParallel)
library(parallel)
##############################################################################


geno_filter_with_sample_list <- function(geno_folder, geno_threshold, maf_threshold, hwe_threshold, sample_list, out_folder) {
  # Load required library
  library(parallel)
  
  # Generate the list of input files
  files <- list.files(path = geno_folder, pattern = "*.bed", full.names = TRUE)
  files <- gsub("\\.bed$", "", files)
  # Create the output folder if it doesn't exist
  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  # Create a sample list file

  write.table(sample_list, file = file.path(out_folder, '/sample_list_plink_format.txt'), row.names = FALSE, col.names = FALSE, quote = F)
  sample_list_file <- paste0(out_folder, '/sample_list_plink_format.txt')
  # Define a function for filtering genotypes
  geno_filter <- function(file) {
    # Sample filter
    system2(
      command = "plink2",
      args = c(
        "--bfile", file,
        "--keep", sample_list_file,
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample"))
      )
    )
    # Geno filter
    system2(
    command = "plink",
    args = c(
        "--bfile", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample")),
        "--geno", geno_threshold,
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample_geno"))
    )
    ) 
    
    # MAF filter
    system2(
      command = "plink",
      args = c(
        "--bfile", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample_geno")),
        "--maf", maf_threshold,
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample_geno_maf"))
      )
    )
    
    # HWE filter
    system2(
      command = "plink",
      args = c(
        "--bfile", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample_geno_maf")),
        "--hwe", hwe_threshold,
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(out_folder, paste0(tools::file_path_sans_ext(basename(file)), "_sample_geno_maf_hwe"))
      )
    )
  }
  
  # Run the genotype filtering in parallel
  cl <- makeCluster(parallelly::availableCores())
  parLapply(cl, files, geno_filter)
  stopCluster(cl)
  
  # Combine chromosomes 1-22 into a single file
  all_chromosome <- file.path(out_folder, "all_chromosome.txt")
  writeLines(paste0(out_folder,"/chr", 2:22, "_SNP_unique_sample_geno_maf_hwe"), con = all_chromosome)
  
  system2(
    command = "plink",
    args = c(
      "--bfile", file.path(out_folder, "chr1_SNP_unique_sample_geno_maf_hwe"),
      "--merge-list", all_chromosome,
      "--make-bed",
      "--allow-no-sex",
      "--out", file.path(out_folder, "chr1-22_SNP_clean")
    )
  )

}



##############################################################################

run_relatedness_analysis <- function(input_QC_genotype_file, output_dir) {
    # Create the output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir)
    }

    # get independent SNPs
    system2(
        command = "plink",
        args = c(
        "--bfile", file.path(input_QC_genotype_file),
        "--indep-pairwise 80 8 0.2",
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, "chr1-22_indep")
        )
    )
    # extract independent SNPs
    system2(
        command = "plink",
        args = c(
        "--bfile", file.path(output_dir, "chr1-22_indep"),
        "--extract", file.path(output_dir, "chr1-22_indep.prune.in"),
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, "indep.IBD")
        )
    )

    # Pairwise IBD estimation
    system2(
        command = "plink",
        args = c(
        "--bfile", file.path(output_dir, "indep.IBD"),
        "--genome",
        "--allow-no-sex",
        "--out", file.path(output_dir, "indep.IBD")
        )
    ) 
    # generate list of subjects with some indication of relatedness
    df_tmp <- fread(file.path(output_dir, "indep.IBD.genome"), data.table = F)
    relatives <- filter(df_tmp, PI_HAT > 0.2)
    MZ_twins <- filter(relatives, PI_HAT >= 0.95)
    other_relatedness <- filter(relatives, PI_HAT < 0.95)
    MZ_twins_list1 <- MZ_twins[,1:2]
    MZ_twins_list2 <- MZ_twins[,3:4]
    colnames(MZ_twins_list2) <- colnames(MZ_twins_list1)
    MZ_twins_list <- rbind(MZ_twins_list1, MZ_twins_list2)
    unique_MZ_twins_drop_list <- unique(MZ_twins_list[order(MZ_twins_list$IID1),])
    other_drop_list <- other_relatedness[,1:2]
    all_drop_list <- rbind(other_drop_list,unique_MZ_twins_drop_list)
    all_drop_list_unique <- unique(all_drop_list[order(all_drop_list$IID1),])
    colnames(all_drop_list_unique) <- c('FID', 'IID')
    write.table(all_drop_list_unique, file.path(output_dir, "drop_list.txt"), sep = '\t', row.names = F, quote = F, col.names = F)
}


##############################################################################

convert_chr_pos_to_rsID <- function(input_genotype_file, output_genotype_file, GRCh_version, hyphen){
    if (GRCh_version == 37 & hyphen == ':'){
        convert_file <- '/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh37_to_rsID.txt'
    }else if (GRCh_version == 37 & hyphen == '_') {
       convert_file <- '/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh37_to_rsID2.txt'
    }else if (GRCh_version == 38 & hyphen == ':') {
       convert_file <- '/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh38_to_rsID.txt'
    }else if (GRCh_version == 38 & hyphen == '_') {
       convert_file <- '/mnt/lvm_vol_2/sliu/database/dbSNP/GRCh38_to_rsID2.txt'
    }

    system2(
    command = "plink2",
    args = c(
        "--bfile", input_genotype_file,
        "--update-name", convert_file,
        "--make-bed",
        "--allow-no-sex",
        "--out", output_genotype_file
        )
    )
}


##############################################################################
# PCA using EIGENSOFT

# Long_LD_range for European
df_Long_LD_range <- data.frame(
  Col1 = c(6, 8, 1, 2, 2, 3, 3, 5, 5, 5, 6, 8, 11, 12, 20, 11, 2, 3, 5, 6, 7, 8, 10, 12),
  Col2 = c(25392021, 111930824, 48227413, 134783530, 183291755, 47524996, 83417310, 97972100, 128972101, 135472101, 139958307, 7962590, 87860352, 111015617, 32536339, 46043424, 86146489, 88917310, 44464243, 56892041, 55032506, 42880843, 36959994, 33108733),
  Col3 = c(33392022, 114930824, 52227412, 138283530, 190291755, 50024996, 86917310, 100472101, 131972101, 138472101, 142458307, 11962591, 90860352, 113515617, 35066586, 57243424, 101133568, 96017310, 50464243, 63942041, 66362565, 49837447, 43679994, 41713733),
  Col4 = c("R1", "R2", "R4", "R6", "R7", "R8", "R9", "R12", "R13", "R14", "R16", "R18", "R21", "R23", "R24", "R3", "R5", "R10", "R11", "R15", "R17", "R19", "R20", "R22")
)


config_file <- data.frame(
  Parameter = c("genotypename:", "snpname:", "indivname:", "outputformat:", "genotypeoutname:", "snpoutname:", "indivoutname:", "familynames:"),
  Value = c("independent_SNP.ped", "independent_SNP.map", "independent_SNP.fam", "EIGENSTRAT", "independent_SNP_PCA.geno", "independent_SNP_PCA.snp", "independent_SNP_PCA.ind", "NO")
)


generate_PCA <- function(input_genotype_file, Ethnic, output_dir){
    # Create the output folder if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    }


    #filter geno using plink
    system2(
    command = "plink",
    args = c(
        "--bfile", input_genotype_file,
        "--geno", "0.05",
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, "geno.IBD")
        )
    )
    #filter maf using plink
    if (Ethnic == "Eur"){
      write.table(df_Long_LD_range,  file.path(output_dir, "Long_LD_range_Eur.txt"), sep = '\t', quote = F, sep = '\t', row.names = F)
      system2(
      command = "plink",
      args = c(
          "--bfile", file.path(output_dir, "geno.IBD"),
          "--exclude range", file.path(output_dir, "Long_LD_range_Eur.txt"),
          "--maf", "0.05",
          "--make-bed",
          "--allow-no-sex",
          "--out", file.path(output_dir, "geno_maf.IBD")
          )
      )
    }else if (Ethnic == "non-Eur") {
      system2(
      command = "plink",
      args = c(
          "--bfile", file.path(output_dir, "geno.IBD"),
          "--maf", "0.05",
          "--make-bed",
          "--allow-no-sex",
          "--out", file.path(output_dir, "geno_maf.IBD")
          )
      )
    }
    #get independent SNPs
    system2(
    command = "plink",
    args = c(
        "--bfile", file.path(output_dir, "geno_maf.IBD"),
        "--indep-pairwise", "200 100 0.2",
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, "independent_SNP_list")
        )
    )
    # extract independent SNPs
    system2(
    command = "plink",
    args = c(
        "--bfile", file.path(output_dir, "geno_maf.IBD"),
        "--extract", file.path(output_dir, "independent_SNP_list.prune.in"),
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, "independent_SNP")
        )
    )

    df_tmp <- fread(paste0(file.path(output_dir, "independent_SNP"), '.fam'), data.table = F)
    df_tmp$V6 <- 1
    write.table (df_tmp, paste0(file.path(output_dir, "independent_SNP"), '.fam'), row.names = F, col.names = F, sep = ' ', quote = F)




    # convert to ped
    system2(
    command = "plink",
    args = c(
        "--bfile", file.path(output_dir, "independent_SNP"),
        "--recode",
        "--allow-no-sex",
        "--out", file.path(output_dir, "independent_SNP")
        )
    )

    config_file$Value <- paste0(file.path(output_dir), '/', config_file$Value)
    config_file[8,]$Value <- 'NO'
    write.table(config_file, file.path(output_dir, "Parameter.PED.EIGENSTRAT"), row.names = F, col.names = F, sep = ' ', quote = F)


    # convert file to EIGENSTRAT format
    system2(
    command = "convertf",
    args = c(
        "-p", file.path(output_dir, "Parameter.PED.EIGENSTRAT")
        )
    )
    # PCA analysis using smartpca
    system2(
    command = "smartpca",
    args = c(
        "-i", file.path(output_dir, "independent_SNP_PCA.geno"),
        "-a", file.path(output_dir, "independent_SNP_PCA.snp"),
        "-b", file.path(output_dir, "independent_SNP_PCA.ind"),
        "-o", file.path(output_dir, "independent_SNP_PCA.pca"),
        "-e", file.path(output_dir, "independent_SNP_PCA.eval"),
        "-p", file.path(output_dir, "independent_SNP_PCA.plot"),
        "-l", file.path(output_dir, "independent_SNP_PCA.log"),
        "-m", 0
        )
    )
}

##############################################################################
pheno_log_invrank_residual_invrank <- function(phenotype_data, cov_matrix, output_dir){ # cov generate from genotype, sample order is correct.
    # Create the output folder if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir)
    } 
    pheno <- fread(phenotype_data, data.table = F)
    colnames(pheno)[1] <- 'ID'
    covariates <- fread(cov_matrix, data.table = F)
    colnames(covariates)[1] <- 'ID'
    pheno_filter <- select(covariates, ID)
    pheno_filter <- left_join(pheno_filter, pheno, by = c('ID' = 'ID'))
    for (i in 2:ncol(pheno_filter)){
      dat <- data.frame(y=RankNorm(log(pheno_filter[,i])), covariates[,2:ncol(covariates)])
      fit <- lm(y~., dat)
      res <- data.frame(FID = 0, IID = covariates$ID, y=RankNorm(residuals(fit)))
      colnames(res)[3] <- colnames(pheno_filter[i])
      write.table(res, paste0(file.path(output_dir, "/", colnames(pheno_filter[i])), '.pheno'), row.names = F, sep = '\t', quote = F)

    }
}


##############################################################################
glm <- function(input_QC_genotype_file, annotation_file, cis_flank_region, trans_threshold, phenotype_data_folder, output_dir){
  # create output_dir
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # create cis dir
  cis_dir <- file.path(output_dir, 'cis')
  if (!dir.exists(cis_dir)) {
    dir.create(cis_dir)
  }
  # create trans dir
  trans_dir <- file.path(output_dir, 'trans')
  if (!dir.exists(trans_dir)) {
    dir.create(trans_dir)
  }
  # create range dir
  range_dir <- file.path(output_dir, 'range')
  if (!dir.exists(range_dir)) {
    dir.create(range_dir)
  }

  # read annotation file
  annotation <- fread(annotation_file, data.table = FALSE)
  colnames(annotation)[1] <- 'ID'
  ID_unique <- unique(annotation$ID)
  phenos <- Sys.glob(paste0(phenotype_data_folder, '/*pheno'))
  phenos <- data.frame(phenos)
  phenos$last_element <- sapply(strsplit(phenos$phenos, '/'), function(x) tail(x, 1))
  phenos$last_element <- gsub("\\.pheno$", "", phenos$last_element)
  ID_unique <- intersect(ID_unique, phenos$last_element)

  # Initialize parallel processing
  num_cores <- parallelly::availableCores() 
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Define the foreach loop for parallel processing
  foreach(id = ID_unique, .packages = c("data.table", "dplyr", "readr", "stringr"), .combine = rbind) %dopar% {
    annotation_tmp <- filter(annotation, ID == id)
    annotation_tmp$cis_start <- annotation_tmp$TSS - cis_flank_region
    annotation_tmp$cis_start[annotation_tmp$cis_start < 0] <- 0
    annotation_tmp$cis_end <- annotation_tmp$TSS + cis_flank_region
    annotation_tmp <- select(annotation_tmp, Chr, cis_start, cis_end, ID)
    outfile_name <- paste0(id, '_range.txt')
    write.table(annotation_tmp, file.path(range_dir, outfile_name), sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)

    # extract SNPs in cis region
    system2(
      command = "plink",
      args = c(
        "--bfile", input_QC_genotype_file,
        "--extract range", file.path(range_dir, outfile_name),
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(cis_dir, paste0(id, "_cis"))
      )
    )

    # extract SNPs in trans region
    system2(
      command = "plink",
      args = c(
        "--bfile", input_QC_genotype_file,
        "--exclude range", file.path(range_dir, outfile_name),
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(trans_dir, paste0(id, "_trans"))
      )
    )

    # perform glm analysis for cis region
    system2(
      command = "plink2",
      args = c(
        "--bfile", file.path(cis_dir, paste0(id, "_cis")),
        "--pheno", file.path(phenotype_data_folder, paste0(id, '.pheno')),
        "--glm", 
        "hide-covar",
        "--allow-no-sex",
        "--out", file.path(cis_dir, id)
        # "--adjust"
      )
    )

    # perform glm analysis for trans region
    system2(
      command = "plink2",
      args = c(
        "--bfile", file.path(trans_dir, paste0(id, "_trans")),
        "--pheno", file.path(phenotype_data_folder, paste0(id, '.pheno')),
        "--glm", 
        "hide-covar",
        "--allow-no-sex",
        "--pfilter", trans_threshold,
        "--out", file.path(trans_dir, id)
      )
    )
    system(paste0('rm ', file.path(cis_dir, paste0(id, "_cis*"))))
    system(paste0('rm ', file.path(trans_dir, paste0(id, "_trans*"))))

  }

  # Stop parallel processing and close the cluster
  stopCluster(cl)
}



##############################################################################
extract_potential_SNP_predictors <- function(cis_glm_folder, trans_glm_folder, flank_region, FDR_threshold, input_QC_genotype_file, strict_SNP, output_dir){
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # get all ambigous SNP list
  all_SNP <- fread(paste0(input_QC_genotype_file,'.bim'))
  ambiguous_snps <- c("A/T", "T/A", "C/G", "G/C")
  all_SNP_ambigous <- all_SNP[(paste0(all_SNP$V5, '/', all_SNP$V6) %in% ambiguous_snps), ]
  write.table(all_SNP_ambigous$V2, file.path(output_dir,'ambiguous_snp.txt'), row.names = F, col.names = F, sep = '\t', quote = F)

  # Initialize parallel processing
  num_cores <- parallelly::availableCores()  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)


  # get all file names
  files <- list.files(cis_glm_folder, pattern = 'glm.linear')
  # loop each file
  foreach(file = files, .packages = c("data.table", "dplyr", "stringr", "valr")) %dopar% {
    # extract cis sig SNP
    cis_tmp <- fread(file.path(cis_glm_folder,file), data.table = F)
    cis_tmp$FDR <- p.adjust(cis_tmp$P, method = 'BH')
    cis_tmp <- filter(cis_tmp, FDR <= FDR_threshold)
    cis_tmp <- subset(cis_tmp, select = -FDR)
    # extract trans sig SNP
    trans_tmp <- fread(file.path(trans_glm_folder,file), data.table = F)
    # merge cis and trans sig SNP
    merge <- rbind(cis_tmp, trans_tmp)
    # # remove ambiguous SNP
    # merge <- merge[!(paste0(merge$REF, '/', merge$ALT) %in% ambiguous_snps), ]
    # get region
    merge <- select(merge, '#CHROM', POS)
    merge$start <- merge$POS - flank_region
    merge$start <- ifelse(merge$start < 0, 0, merge$start) # fix if start < 0
    merge$end <- merge$POS + flank_region
    colnames(merge) <- c('chrom', 'pos', 'start', 'end')
    combined <- data.frame(bed_merge(merge))
    combined <- filter(combined, chrom != 'NA')
    if (nrow(combined) > 0){
      combined$ID <- substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2)
      write.table(combined, file.path(output_dir,paste0(substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2),'_range.txt')), row.names = F, col.names = F, quote = F, sep = '\t')
      
      # generate potential SNP predictors
      system2(
      command = "plink",
      args = c(
        "--bfile", input_QC_genotype_file,
        "--extract range", file.path(output_dir,paste0(substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2),'_range.txt')),
        "--exclude", file.path(output_dir,'ambiguous_snp.txt'),
        # "--update-chr", file.path(output_dir,'Chr_converter.txt'),
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, paste0(substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2))) 
        #"--allow-extra-chr"
        )
      )

      # filter potential SNPs in 1000G LD reference
      system2(
      command = "plink",
      args = c(
        "--bfile", file.path(output_dir, paste0(substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2))),
        "--extract", strict_SNP,
        "--make-bed",
        "--allow-no-sex",
        "--out", file.path(output_dir, paste0(substr(file, 1, (nchar(gsub('.glm.linear','', file )) -1)/2), '_filtered')) 
        )
      )   

    }
  }
  # Stop parallel processing and close the cluster
  stopCluster(cl)
}

##############################################################################
establish_prediction_models <- function(genotype_file_folder, pheno_folder, gcta_PATH, gemma_PATH, hsq_p, PANEL, reference_1000G, output_dir){
  # creat models and tmp dirs
  tmp_folder <- paste0(output_dir, '/tmp')
  out_folder <- paste0(output_dir, '/models')
  ref_ld_folder <- paste0(output_dir, '/ref_ld')
  pos_folder <- paste0(output_dir, '/pos')
  new_output <- file.path('output', output_dir, 'tmp')
  for (item in c(tmp_folder, out_folder, ref_ld_folder, pos_folder, new_output)){
    if (!dir.exists(item)){
       dir.create(item, recursive = TRUE)
    }
  }

  # Initialize parallel processing
  num_cores <- parallelly::availableCores()  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  all_files <- list.files(genotype_file_folder, pattern = '_filtered.bim')
  #for (file in all_files){

  foreach(file = all_files, .packages = c("dplyr", "stringr")) %dopar% {  
    pheno_ID <- gsub('_filtered.bim', '', file)
    genotype_file <- file.path(genotype_file_folder, paste0(pheno_ID, '_filtered'))
    pheno <- file.path(pheno_folder, paste0(pheno_ID, '.pheno'))
    
    # establish models
    system2(
    command = "/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript",
      args = c(
          "/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/bin/FUSION.compute_weights_cis_trans.R",
          "--bfile", genotype_file,
          "--tmp", file.path(tmp_folder, pheno_ID),
          "--out", file.path(out_folder, pheno_ID),
          "--verbose", 0,
          "--pheno", pheno,
          "--noclean",
          "--PATH_gcta", gcta_PATH,
          "--PATH_gemma", gemma_PATH,
          "--hsq_p", 1,
          "--save_hsq", TRUE
      )
    )
    # generate ld ref for each models
    if (file.exists(file.path(out_folder, paste0(pheno_ID, '.wgt.RDat')))){
      load(file.path(out_folder, paste0(pheno_ID, '.wgt.RDat')))
      SNP_in_models <- select(snps, V2)
      write.table (SNP_in_models, file.path(out_folder, paste0(pheno_ID, '_SNPs_in_models.txt')), row.names = F, sep = '\t', quote = F, col.names =F)
      system2(
        command = "plink",
        args = c(
          "--bfile", reference_1000G,
          "--extract", file.path(out_folder, paste0(pheno_ID, '_SNPs_in_models.txt')),
          "--make-bed",
          "--allow-no-sex",
          "--out", file.path(ref_ld_folder, paste0(pheno_ID, '_1000G'))
          )
        )

    }
    
  }
  # Stop parallel processing and close the cluster
  stopCluster(cl)
  
  # generate pos file

  pos_df <- data.frame(PANEL = PANEL, WGT = list.files(paste0(output_dir, '/models'), pattern = 'wgt.RDat'))
  pos_df$ID <- gsub('.wgt.RDat', '',pos_df$WGT)

  pos_df$N <- nrow(fread(gsub('_filtered.bim', '_filtered.fam', file.path(genotype_file_folder,all_files[1])), data.table = F))
  for (i in 1:nrow(pos_df)){
    write.table(pos_df[i,], file.path(pos_folder, paste0(pos_df[i,]$ID, '.pos')), row.names = F, sep = '\t', quote = F)
  }
  

}
  

##############################################################################
association_test <- function(sumstats, pos_dir, model_dir, ref_ld_dir, TSS_file, output_dir){
  # creat result dir
  if (!dir.exists(output_dir)){
      dir.create(output_dir, recursive = TRUE)
  }
  all_models <- list.files(model_dir, pattern = '.wgt.RDat')
  # Initialize parallel processing
  num_cores <- parallelly::availableCores()  
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  foreach(file = all_models, .packages = c("dplyr", "stringr")) %dopar% {
    pheno_ID <- gsub('.wgt.RDat', '', file)
    # association test
    system2(
    command = "/mnt/lvm_vol_1/hzhong/R-4.1.3/bin/Rscript",
      args = c(
          "/mnt/lvm_vol_2/sliu/project/MESA_SomaScan_1K_PWAS/00_analysis_code/bin/FUSION.assoc_test_cis_trans.R",
          "--sumstats", sumstats,
          "--weights", file.path(pos_dir, paste0(pheno_ID, '.pos')),
          "--weights_dir", file.path(model_dir),
          "--ref_ld", file.path(ref_ld_dir, paste0(pheno_ID, '_1000G')),
          "--TSS_file", TSS_file,
          "--out",  file.path(output_dir, paste0(pheno_ID, '_association.txt'))
      )
    )


  }
  # Stop parallel processing and close the cluster
  stopCluster(cl)



}
