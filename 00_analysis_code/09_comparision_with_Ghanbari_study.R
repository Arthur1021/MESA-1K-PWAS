library(data.table)
library(dplyr)
EA_res <- fread('03_association/association_T2D_MVP_EA/all_association.out', data.table = F)
EA_res <- select(EA_res, ID, PWAS.Z, PWAS.P, FDR)
annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- select(annotation, SomaId, TargetFullName, EntrezGeneSymbol)
EA_res <- inner_join(EA_res, annotation, by = c('ID' = 'SomaId'))

Ghanbari_res <- fread('09_comparision_with_Ghanbari_and_Gudmundsdottir_study/Table_S2_from_Ghanbari.txt.txt', data.table = F)
Ghanbari_res <- Ghanbari_res[Ghanbari_res$'P value' < 0.05/1089,]

combine <- inner_join(Ghanbari_res, EA_res, by = c('Protein' = 'EntrezGeneSymbol'))
combine <- cbind(TargetFullName = combine[,19], combine[,c(1:18)])
combine$consistent <- ifelse(sign(combine$'MR Beta') == sign(combine$PWAS.Z), 'TRUE', 'FALSE')
combine$replicated <- ifelse(combine$consistent == 'TRUE' & combine$PWAS.P < 0.05, 'TRUE', 'FALSE')
write.table(combine, '09_comparision_with_Ghanbari_and_Gudmundsdottir_study/compare_result_Ghanbari.txt', row.names = F, sep = '\t', quote = F)

Gudmundsdottir_res <- fread('09_comparision_with_Ghanbari_and_Gudmundsdottir_study/Table_S10_from_Gudmundsdottir.txt', data.table = F)
Gudmundsdottir_res <- Gudmundsdottir_res[Gudmundsdottir_res$FDR < 0.05,]
combine <- inner_join(Gudmundsdottir_res, EA_res, by = c('Protein (Entrez symbol)' = 'EntrezGeneSymbol'))
combine <- cbind(TargetFullName = combine[,12], combine[,c(1:11)])
combine$consistent <- ifelse(sign(combine$'beta') == sign(combine$PWAS.Z), 'TRUE', 'FALSE')
combine$replicated <- ifelse(combine$consistent == 'TRUE' & combine$PWAS.P < 0.05, 'TRUE', 'FALSE')
write.table(combine, '09_comparision_with_Ghanbari_and_Gudmundsdottir_study/compare_result_Gudmundsdottir.txt', row.names = F, sep = '\t', quote = F)

