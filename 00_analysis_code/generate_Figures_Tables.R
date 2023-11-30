

# Table S1
if (!dir.exists('Figure_Tables')){
    dir.create('Figure_Tables')
}

library(data.table)
library(dplyr)

chara_df <- data.frame(
    race = as.character(),
    N = as.numeric(),
    age55 = as.character(),
    age56_60 = as.character(),
    age61_65 = as.character(),
    age65 = as.character(),
    female = as.character(),
    bmi = as.character(),
    never = as.character(),
    former = as.character(),
    current = as.character()
)

races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){
    cov <- fread(paste0('01_preprocess_data/covariates/', race, '.txt'), data.table = F)
    independent <- fread(paste0('02_establish_models/' , race, '/03_independent_subjects_genotype/chr1-22_SNP_clean_rsID.fam'), data.table = F)
    cov <- filter(cov, Sample_ID %in% independent$V2)
    N <- nrow(cov)
    age55 <- paste0(nrow(filter(cov, age1c <= 55)), '(', round(nrow(filter(cov, age1c <= 55))/nrow(cov)*100, 2), '%)')
    age56_60 <- paste0(nrow(filter(cov, age1c > 55 & age1c <= 60)), '(', round(nrow(filter(cov, age1c > 55 & age1c <= 60))/nrow(cov)*100, 2), '%)')
    age61_65 <- paste0(nrow(filter(cov, age1c > 60 & age1c <= 65)), '(', round(nrow(filter(cov, age1c > 60 & age1c <= 65))/nrow(cov)*100, 2), '%)')
    age65 <- paste0(nrow(filter(cov, age1c > 65)), '(', round(nrow(filter(cov, age1c > 65))/nrow(cov)*100, 2), '%)')
    female <- paste0(nrow(filter(cov, gender1 == 0)) , '(',  round(nrow(filter(cov, gender1 == 0)) / nrow(cov) * 100, 2), '%)')
    bmi <- paste0(round(mean(cov$bmi), 2), ' ± ', round(sd(cov$bmi),2))
    never <- paste0(nrow(filter(cov, cig1c == 0)), '(', round(nrow(filter(cov, cig1c == 0))/nrow(cov)*100, 2), '%)')
    former <- paste0(nrow(filter(cov, cig1c == 1)), '(', round(nrow(filter(cov, cig1c == 1))/nrow(cov)*100, 2), '%)')
    current <- paste0(nrow(filter(cov, cig1c == 2)), '(', round(nrow(filter(cov, cig1c == 2))/nrow(cov)*100, 2), '%)')
    tmp <- c(race, N, age55, age56_60, age61_65, age65, female, bmi, never, former, current)
    chara_df[nrow(chara_df)+1,] <- tmp
}

row.names(chara_df) <- chara_df$race
chara_df <- chara_df[,-1]
chara_df <- data.frame(t(chara_df))
chara_df <- cbind(row.names(chara_df), chara_df)
colnames(chara_df)[1] <- 'Population'
write.table(chara_df, 'Figure_Tables/TableS1.txt', sep = '\t', quote = F, row.names = F)


# Table S2
library(data.table)
library(dplyr)
models_df <- data.frame(
    race = as.character(),
    N = as.numeric(),
    average = as.numeric(),
    min = as.numeric(),
    max = as.numeric(),
    N_cis = as.numeric(),
    N_trans = as.numeric(),
    N_cis_trans = as.numeric(),
    N_cis_complex = as.numeric(),
    N_trans_complex = as.numeric(),
    N_cis_trans_complex = as.numeric(),
    N_cis_single = as.numeric(),
    N_trans_single = as.numeric(),
    N_cis_trans_single = as.numeric()
)
if (!dir.exists('Figure_Tables')){
    dir.create('Figure_Tables')
}

races <- c('AA', 'CA', 'EA', 'HA')
for (race in races){
    TSS <- fread('01_preprocess_data/annotation/SomaScan_TSS.txt', data.table = F)
    df <- fread(paste0('02_establish_models/', race, '/08_establish_prediction_models/model_Rsq_0.01.txt'), data.table = F)
    df$ID <- gsub('.wgt.RDat', '', df$model)
    complex <- unique(TSS$SomaId[duplicated(TSS$SomaId)])
    N <- nrow(df)
    average <- round(mean(df$Rsq), 2)
    min <- round(min(df$Rsq), 2)
    max <- round(max(df$Rsq), 2)

    N_cis <- nrow(filter(df, df$N_trans == 0))
    N_trans <- nrow(filter(df, df$N_cis == 0))
    N_cis_trans <- nrow(filter(df, df$N_cis != 0 & df$N_trans != 0))

    df_complex <- filter(df, ID %in% complex)
    N_cis_complex <- nrow(filter(df_complex, df_complex$N_trans == 0))
    N_trans_complex <- nrow(filter(df_complex, df_complex$N_cis == 0))
    N_cis_trans_complex <- nrow(filter(df_complex, df_complex$N_cis != 0 & df_complex$N_trans != 0))

    df_single <- filter(df, !ID %in% complex)
    N_cis_single <- nrow(filter(df_single, df_single$N_trans == 0))
    N_trans_single <- nrow(filter(df_single, df_single$N_cis == 0))
    N_cis_trans_single <- nrow(filter(df_single, df_single$N_cis != 0 & df_single$N_trans != 0))
    tmp <- c(race, N, average, min, max, N_cis, N_trans, N_cis_trans, N_cis_complex, N_trans_complex, N_cis_trans_complex, N_cis_single, N_trans_single, N_cis_trans_single)
    models_df[nrow(models_df)+ 1, ]<- tmp
}
write.table(models_df, 'Figure_Tables/TableS2.txt', sep = '\t', quote = F, row.names = F)



# Table S3

annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, SomaId, TargetFullName, EntrezGeneSymbol, chromosome_name, TSS)
annotation <- annotation[!is.na(annotation$TSS),]
annotation_merge <- annotation %>%
  group_by(SomaId, TargetFullName) %>%
  summarise(
    EntrezGeneSymbol = paste(EntrezGeneSymbol, collapse = ' '),
    chromosome_name = paste(chromosome_name, collapse = '/'),
    TSS = paste(TSS, collapse = '/')
  ) %>%
  ungroup()
annotation_merge <- data.frame(annotation_merge)


MVP_meta <- fread('04_meta_analysis/T2D_MVP_meta1.tbl', data.table = F)
MVP_meta <- select(MVP_meta, MarkerName, HetISq, HetPVal, Zscore, 'P-value')
MVP_meta <- MVP_meta[order(MVP_meta$'P-value'),]
MVP_meta$FDR <- p.adjust(MVP_meta$'P-value', method = 'fdr')
colnames(MVP_meta)[1] <- 'ID'

MVP_AA <- fread('03_association/association_T2D_MVP_AA/all_association.out', data.table = F)
MVP_AA <- select(MVP_AA, ID, N_SNP, MODEL, MODELCV.R2, PWAS.Z, PWAS.P, FDR)
MVP_meta <- full_join(MVP_meta, MVP_AA, by = 'ID')

MVP_CA <- fread('03_association/association_T2D_MVP_CA/all_association.out', data.table = F)
MVP_CA <- select(MVP_CA, ID, N_SNP, MODEL, MODELCV.R2, PWAS.Z, PWAS.P, FDR)
MVP_meta <- full_join(MVP_meta, MVP_CA, by = 'ID')

MVP_EA <- fread('03_association/association_T2D_MVP_EA/all_association.out', data.table = F)
MVP_EA <- select(MVP_EA, ID, N_SNP, MODEL, MODELCV.R2, PWAS.Z, PWAS.P, FDR)
MVP_meta <- full_join(MVP_meta, MVP_EA, by = 'ID')

MVP_HA <- fread('03_association/association_T2D_MVP_HA/all_association.out', data.table = F)
MVP_HA <- select(MVP_HA, ID, N_SNP, MODEL, MODELCV.R2, PWAS.Z, PWAS.P, FDR)
MVP_meta <- full_join(MVP_meta, MVP_HA, by = 'ID')

colnames(MVP_meta) <- c('ID', 'HetISq', 'HetPVal', 'Zscore', 'Pvalue', 'FDR', 'NSNP_AA', 'model_AA', 'Rsq_AA', 'Z_AA', 'P_AA', 'FDR_AA', 'NSNP_CA', 'model_CA', 'Rsq_CA', 'Z_CA','P_CA', 'FDR_CA', 'NSNP_EA', 'model_EA', 'Rsq_EA', 'Z_EA','P_EA', 'FDR_EA', 'NSNP_HA', 'model_HA', 'Rsq_HA', 'Z_HA','P_HA', 'FDR_HA')
MVP_meta <- filter(MVP_meta, FDR < 0.05 | FDR_AA < 0.05 | FDR_CA < 0.05 |FDR_EA < 0.05 |FDR_HA < 0.05 )


res <- right_join(annotation_merge, MVP_meta, by = c('SomaId' = 'ID'))
res <- res[order(res$Pvalue),]

write.table(res, 'Figure_Tables/TableS3.txt', sep = '\t', quote = F, row.names = F)



# generate Table S4
MVP_PAGE_meta <- fread('04_meta_analysis/T2D_MVP_PAGE_meta1.tbl', data.table = F)
MVP_PAGE_meta <- MVP_PAGE_meta[order(MVP_PAGE_meta$'P-value'),]
MVP_PAGE_meta$FDR <- p.adjust(MVP_PAGE_meta$'P-value', method = 'fdr')
MVP_PAGE_meta <- select(MVP_PAGE_meta,  MarkerName, Zscore, 'P-value', FDR, Direction, HetISq, HetPVal)
colnames(MVP_PAGE_meta)[1] <- 'ID'

MVP_PAGE_meta_AA <- fread('04_meta_analysis/T2D_MVP_PAGE_AA_meta1.tbl', data.table = F)
MVP_PAGE_meta_AA <- MVP_PAGE_meta_AA[order(MVP_PAGE_meta_AA$'P-value'),]
MVP_PAGE_meta_AA$FDR <- p.adjust(MVP_PAGE_meta_AA$'P-value', method = 'fdr')
MVP_PAGE_meta_AA <- select(MVP_PAGE_meta_AA,  MarkerName, Zscore, 'P-value', FDR, Direction, HetISq, HetPVal)
colnames(MVP_PAGE_meta_AA)[1] <- 'ID'

MVP_PAGE_meta_CA <- fread('04_meta_analysis/T2D_MVP_PAGE_CA_meta1.tbl', data.table = F)
MVP_PAGE_meta_CA <- MVP_PAGE_meta_CA[order(MVP_PAGE_meta_CA$'P-value'),]
MVP_PAGE_meta_CA$FDR <- p.adjust(MVP_PAGE_meta_CA$'P-value', method = 'fdr')
MVP_PAGE_meta_CA <- select(MVP_PAGE_meta_CA,  MarkerName, Zscore, 'P-value', FDR, Direction, HetISq, HetPVal)
colnames(MVP_PAGE_meta_CA)[1] <- 'ID'

MVP_EA <- fread('03_association/association_T2D_MVP_EA/all_association.out', data.table = F)
MVP_EA <- select(MVP_EA, ID, PWAS.Z, PWAS.P, FDR)
colnames(MVP_EA) <- c('ID', 'Zscore', 'Pvalue', 'FDR')

MVP_PAGE_meta_HA <- fread('04_meta_analysis/T2D_MVP_PAGE_HA_meta1.tbl', data.table = F)
MVP_PAGE_meta_HA <- MVP_PAGE_meta_HA[order(MVP_PAGE_meta_HA$'P-value'),]
MVP_PAGE_meta_HA$FDR <- p.adjust(MVP_PAGE_meta_HA$'P-value', method = 'fdr')
MVP_PAGE_meta_HA <- select(MVP_PAGE_meta_HA,  MarkerName, Zscore, 'P-value', FDR, Direction, HetISq, HetPVal)
colnames(MVP_PAGE_meta_HA)[1] <- 'ID'

meta_7_data <- full_join(MVP_PAGE_meta, MVP_PAGE_meta_AA, by = 'ID')
meta_7_data <- full_join(meta_7_data, MVP_PAGE_meta_CA, by = 'ID')
meta_7_data <- full_join(meta_7_data, MVP_EA, by = 'ID')
meta_7_data <- full_join(meta_7_data, MVP_PAGE_meta_HA, by = 'ID')

colnames(meta_7_data) <- c('ID', 'Z_multi', 'P_multi', 'FDR_multi', 'Direction_multi', 'HetISq_multi', 'HetPVal_multi',  'Z_AA', 'P_AA', 'FDR_AA', 'Direction_AA', 'HetISq_AA', 'HetPVal_AA', 'Z_CA', 'P_CA', 'FDR_CA', 'Direction_CA', 'HetISq_CA', 'HetPVal_CA', 'Z_EA', 'P_EA', 'FDR_EA', 'Z_HA', 'P_HA', 'FDR_HA', 'Direction_HA', 'HetISq_HA', 'HetPVal_HA')
res <- right_join(annotation_merge, meta_7_data, by = c('SomaId' = 'ID'))
res <- res[order(res$P_multi),]
res <- filter(res, FDR_multi < 0.05 | FDR_AA < 0.05 | FDR_CA < 0.05 | FDR_EA < 0.05 | FDR_HA < 0.05 )

write.table(res, 'Figure_Tables/TableS4.txt', sep = '\t', quote = F, row.names = F)



# Figure 2a
library(ggvenn)
library(data.table)
library(dplyr)
AA <- fread('02_establish_models/AA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
AA <- AA$model
CA <- fread('02_establish_models/CA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
CA <- CA$model
EA <- fread('02_establish_models/EA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
EA <- EA$model
HA <- fread('02_establish_models/HA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
HA <- HA$model

x <- list(
  'African' = AA,
  'Asian' = CA,
  'non-Hispanic White' = EA,
  'Hispanic' = HA
)
pdf(file="Figure_Tables/Figure2a.pdf", width = 4, height = 4)
ggvenn(
  x, 
  fill_color = c("#73c376","#6aaad5","#6b54a5","#fe8b3c"),
  stroke_size = 0.8, stroke_color = 'white', set_name_size = 5, digits = 2, fill_alpha = 0.7,
  show_percentage = F,
)
dev.off()


# Figure 2b

library(dplyr)
library(data.table)
library(CMplot)

AA_association <- fread('04_meta_analysis/T2D_MVP_PAGE_AA_meta1.tbl', data.table = F)
AA_association <- AA_association[order(AA_association$'P-value'),]
AA_association$FDR <- p.adjust(AA_association$'P-value', method = 'fdr')
AA_sig <- 0.05/nrow(AA_association)
AA_association <- select(AA_association, MarkerName, 'P-value')


CA_association <- fread('04_meta_analysis/T2D_MVP_PAGE_CA_meta1.tbl', data.table = F)
CA_association <- CA_association[order(CA_association$'P-value'),]
CA_association$FDR <- p.adjust(CA_association$'P-value', method = 'fdr')
CA_sig <- filter(CA_association, FDR < 0.05)
CA_sig <- CA_sig[nrow(CA_sig),]$'P-value'+ 0.00001
CA_association <- select(CA_association, MarkerName, 'P-value')


EA_association <- fread('04_meta_analysis/T2D_MVP_EA.txt', data.table = F)
EA_association <- EA_association[order(EA_association$PWAS.P),]
EA_association$FDR <- p.adjust(EA_association$PWAS.P, method = 'fdr')
EA_sig <- filter(EA_association, FDR < 0.05)
EA_sig <- EA_sig[nrow(EA_sig),]$PWAS.P+ 0.00001

EA_association <- select(EA_association, ID, PWAS.P)
HA_association <- fread('04_meta_analysis/T2D_MVP_PAGE_HA_meta1.tbl', data.table = F)
HA_association <- HA_association[order(HA_association$'P-value'),]
HA_association$FDR <- p.adjust(HA_association$'P-value', method = 'fdr')
HA_sig <- filter(HA_association, FDR < 0.05)
HA_sig <- HA_sig[nrow(HA_sig),]$'P-value'+ 0.00001
HA_association <- select(HA_association, MarkerName, 'P-value')


input <- full_join(HA_association, EA_association, by = c('MarkerName'='ID'))
input <- full_join(input, CA_association, by = 'MarkerName')
input <- full_join(input, AA_association, by = 'MarkerName')
colnames(input) <- c('SNP', 'HA', 'EA', 'CA', 'AA')

annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, SomaId, chromosome_name, TSS)
annotation <- annotation[!duplicated(annotation$SomaId),]

input <- right_join(annotation, input, by = c('SomaId' = 'SNP'))
input[is.na(input)] <- 1

CMplot(input,type="p",plot.type=c("c"),col=matrix(c("#fe8b3c","#6b54a5","#6aaad5","#73c376")),axis.cex = 0.5,lab.cex=1,
        outward=T,ylim = list(c(0,10),c(0,30),c(0,10),c(0,10)),file="pdf",width=6,height=6,chr.den.col="grey50",
        band = 0.5,H = 2,r =3,threshold = list(HA_sig,EA_sig,CA_sig,AA_sig),
        signal.cex = 1, file.name = 'Figure2b')
system('mv Cir_Manhtn.Figure2b.pdf Figure_Tables/Figure2b.pdf')    


# Figure 2c

library(ggplot2)
Ancestry <- c('African','Asian','non-Hispanic White','Hispanic')
y<- c(0,3,44,1)
df <- data.frame(x = Ancestry, y = y)
p2<-ggplot(data = df, mapping = aes(x = Ancestry, y = y,fill =Ancestry)) + geom_bar(stat = 'identity')+
  scale_fill_manual(values = c("#73c376","#6aaad5","#fe8b3c","#6b54a5"))+
  ylab('No. of proteins associated with T2D')+theme_classic()+theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),legend.position = "none")+
      geom_text(aes(label = y), vjust = -0.3)+ scale_x_discrete(limits = Ancestry)+ylim(0,50)+labs(x = NULL)
ggsave("Figure_Tables/Figure2c.pdf",height = 2.6,width = 4.5)


# Figure 2d
library("CMplot")
meta_all <- fread('04_meta_analysis/T2D_MVP_PAGE_meta1.tbl', data.table = F)
meta_all <- select(meta_all, MarkerName, 'P-value')

annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, SomaId, chromosome_name, TSS, EntrezGeneSymbol)
annotation <- annotation[!duplicated(annotation$SomaId),]

meta_all <- right_join(annotation[,1:3], meta_all, by = c('SomaId' = 'MarkerName'))
colnames(meta_all) <- c('SNP', 'Chromosome', 'Position', 'meta')
meta_all <- meta_all[order(meta_all$meta),]
meta_all$FDR <- p.adjust(meta_all$meta, method = 'fdr')
meta_all_sig_threshold <- filter(meta_all, FDR < 0.05)
sig_SNPs <- meta_all_sig_threshold$SNP
sig_SNPs_gene <- data.frame(sig_SNPs)
sig_SNPs_gene <- left_join(sig_SNPs_gene, annotation, by = c('sig_SNPs' = 'SomaId'))
meta_all_sig_threshold <- meta_all_sig_threshold[nrow(meta_all_sig_threshold),]$meta+ 0.00001

CMplot(meta_all[,1:4],type="h",plot.type="m",LOG10=TRUE,highlight=sig_SNPs,highlight.text=sig_SNPs_gene$EntrezGeneSymbol, highlight.type="p",threshold = meta_all_sig_threshold,
        highlight.text.cex = 0.7,highlight.text.col = '#095C83',
        highlight.col=NULL,highlight.cex=1.2,highlight.pch=19,file="pdf",file.name="meta",col=c("#095C83", "#91C2D9"), 
        file.output=TRUE,verbose=TRUE,width=14,height=4.2,band=0.6, main="Meta-analysis of proteins associated with T2D risk in racially and ethnically diversity populations")

system('mv Rect_Manhtn.meta.pdf Figure_Tables/Figure2d.pdf')   

CMplot(meta_all[,1:4],type="h",plot.type="m",LOG10=TRUE, highlight.type="p",threshold = meta_all_sig_threshold,
        highlight.text.cex = 0.7,highlight.text.col = '#095C83',
        highlight.col=NULL,highlight.cex=1.2,highlight.pch=19,file="pdf",file.name="meta",col=c("#095C83", "#91C2D9"), 
        file.output=TRUE,verbose=TRUE,width=12,height=8,band=0.6)

system('mv Rect_Manhtn.meta.pdf Figure_Tables')

# FigureS1
library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales)
library(ggpmisc)
library("RSQLite")
sqlite <- dbDriver("SQLite")

# baseline EA from Heather
dbname <- paste0("06_comparison_with_models_from_Heather_group/models_from_Heather_group/EUR_PCAIR_baseline_models_rho0.1_zpval0.05.db")
db = dbConnect(sqlite,dbname)
query <- function(...) dbGetQuery(db, ...)
df <- query('select * from extra')
N_baseline <- nrow(df)
N_baseline_prop <- fread('06_comparison_with_models_from_Heather_group/validation_EA_models_using_INTERVAL_for_baseline_models.txt', data.table = F)
N_baseline_prop <- round(sum(N_baseline_prop$Rsq >= 0.01)/length(N_baseline_prop$Rsq), 2)

# finemapped EA from Heather
dbname <- paste0("06_comparison_with_models_from_Heather_group/models_from_Heather_group/EUR_PCAIR_dapg_0.001_T_rho0.1_zpval0.05.db")
db = dbConnect(sqlite,dbname)
query <- function(...) dbGetQuery(db, ...)
df <- query('select * from extra')
N_finemapped <- nrow(df)
N_finemapped_prop <- fread('06_comparison_with_models_from_Heather_group/validation_EA_models_using_INTERVAL_for_dapg_0.001_T.txt', data.table = F)
N_finemapped_prop <- round(sum(N_finemapped_prop$Rsq >= 0.01)/length(N_finemapped_prop$Rsq), 2)

#current EA models
df <- fread('02_establish_models/EA/08_establish_prediction_models/model_Rsq_0.01.txt', data.table = F)
N_cis_trans <- nrow(df)
df <- fread('05_validation_using_INTERVAL/validation_Rsq_External_Internal.txt', data.table = F)
N_cis_trans_prop <- round(sum(df$External_Rsq >= 0.01)/length(df$External_Rsq), 2)


data <- data.frame(
  group = as.character(c('cis+trans', 'baseline', 'fine-mapped')),
  number = as.numeric(c(N_cis_trans, N_baseline, N_finemapped)),
  proportion = as.numeric(c(N_cis_trans_prop, N_baseline_prop, N_finemapped_prop))
)

data$group <- factor(data$group,levels = c("cis+trans","baseline","fine-mapped"))
a <-1/max(data$number)
Figure1A <- ggplot(data) + geom_bar(aes(x=group, y=number), stat = "identity",fill="#4D85BD",colour="#4D85BD") +
 geom_line(aes(x=group, y=proportion*max(number)),stat="identity",group = 1)+geom_point(aes(x=group, y=proportion*max(number)))+
 geom_point(aes(x='cis+trans', y=N_cis_trans_prop*max(number)),shape = 18,size = 4)+
    scale_y_continuous(sec.axis = sec_axis(~. *a,name = "External validate proportion"))+
    geom_text(aes(label=number, x=group, y=0.90*number),color = 'white')+
    geom_text(aes(label=proportion, x=group, y=1.2*proportion*max(number)))+
    theme_classic()+xlab('models from different methods')+
    ylab('Number of established models')

# ggsave('FigureS1A.pdf', Figure1A,width = 9, height = 6)


data <- read.table('05_validation_using_INTERVAL/validation_Rsq_External_Internal.txt',header = T,row.names  = 1)
data <- filter(data, External_Rsq >= 0.01)
my.formula <- y ~ x
Figure1B <- ggplot(data = data, aes(x = Internal_Rsq, y = External_Rsq)) +      
  geom_point(size=3,
             alpha=1,
             color="#4D85BD")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab(External~italic(R)^2~(EA))+
  xlab(Cross~validation~italic(R)^2~(EA))+
  stat_poly_eq(use_label(c("eq"))) 
  # annotate("text", x = 0.3, y = 0.7, label = 'y = 0.892 x - 0.0204',color="black")

ggsave('Figure_Tables/FigureS1B.pdf', Figure1B, width = 3, height = 3)



df <- fread('06_comparison_with_models_from_Heather_group/Rsq_compare_current_with_baseline_from_heather_AA.txt', data.table = F)
df$performance <- ifelse(df$current_Rsq >= df$Heather_Rsq, "better", "worse")
percentage <- paste0(round(sum(df$performance == 'better')/length(df$performance)*100, 2), '%')
compare_AA <- ggplot(data=df, aes(x=Heather_Rsq, y=current_Rsq,color = performance))+geom_point(size=3)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values = c("#4D85BD", "#F7903D"))+ylab(cis+trans~italic(R)^2~(AA))+xlab(Baseline~italic(R)^2~(AA))+
theme(legend.position="none")+annotate("text", x = 0.55, y = 0.1, label = percentage,color="#4D85BD")
ggsave('Figure_Tables/FigureS1C.pdf', compare_AA, width = 3, height = 3)


df <- fread('06_comparison_with_models_from_Heather_group/Rsq_compare_current_with_baseline_from_heather_CA.txt', data.table = F)
df$performance <- ifelse(df$current_Rsq >= df$Heather_Rsq, "better", "worse")
percentage <- paste0(round(sum(df$performance == 'better')/length(df$performance)*100, 2), '%')
compare_CA <- ggplot(data=df, aes(x=Heather_Rsq, y=current_Rsq,color = performance))+geom_point(size=3)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values = c("#4D85BD", "#F7903D"))+ylab(cis+trans~italic(R)^2~(CA))+xlab(Baseline~italic(R)^2~(CA))+
theme(legend.position="none")+annotate("text", x = 0.55, y = 0.1, label = percentage,color="#4D85BD")

df <- fread('06_comparison_with_models_from_Heather_group/Rsq_compare_current_with_baseline_from_heather_EA.txt', data.table = F)
df$performance <- ifelse(df$current_Rsq >= df$Heather_Rsq, "better", "worse")
percentage <- paste0(round(sum(df$performance == 'better')/length(df$performance)*100, 2), '%')
compare_EA <- ggplot(data=df, aes(x=Heather_Rsq, y=current_Rsq,color = performance))+geom_point(size=3)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values = c("#4D85BD", "#F7903D"))+ylab(cis+trans~italic(R)^2~(EA))+xlab(Baseline~italic(R)^2~(EA))+
theme(legend.position="none")+annotate("text", x = 0.55, y = 0.1, label = percentage,color="#4D85BD")


df <- fread('06_comparison_with_models_from_Heather_group/Rsq_compare_current_with_baseline_from_heather_HA.txt', data.table = F)
df$performance <- ifelse(df$current_Rsq >= df$Heather_Rsq, "better", "worse")
percentage <- paste0(round(sum(df$performance == 'better')/length(df$performance)*100, 2), '%')
compare_HA <- ggplot(data=df, aes(x=Heather_Rsq, y=current_Rsq,color = performance))+geom_point(size=3)+theme_bw()+theme(panel.grid = element_blank())+
scale_color_manual(values = c("#4D85BD", "#F7903D"))+ylab(cis+trans~italic(R)^2~(HA))+xlab(Baseline~italic(R)^2~(HA))+
theme(legend.position="none")+annotate("text", x = 0.55, y = 0.1, label = percentage,color="#4D85BD")


p <- ggarrange(ggarrange(Figure1A, Figure1B,
                         ncol = 2, labels = c("a", "b"),
                         widths = c(0.75,0.25)),
               ggarrange(compare_AA, compare_CA, compare_EA,compare_HA,
                         ncol = 4, labels = c("c", "d","e", "f"),
                         widths = c(0.25,0.25,0.25, 0.25)),
               nrow = 2, heights = c(0.5,0.5))

ggsave('Figure_Tables/FigureS1.pdf', p,width = 12, height = 6)


