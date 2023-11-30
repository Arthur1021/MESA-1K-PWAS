library(data.table)
library(dplyr)

if (!dir.exists('04_meta_analysis')){
    dir.create('04_meta_analysis')
}

races <- c('MVP_AA', 'MVP_CA', 'MVP_EA', 'MVP_HA', 'PAGE_AA', 'PAGE_CA', 'PAGE_EA', 'PAGE_HA')

sample_size <- data.frame(
    MVP_AA = 56092,
    MVP_CA = 216287,
    MVP_EA = 1114458,
    MVP_HA = 20445,
    PAGE_AA = 25478,
    PAGE_CA = 7437,
    PAGE_EA = 177201,
    PAGE_HA = 32871
)


for (race in races){
    association_df <- fread(paste0('03_association/association_T2D_', race, '/all_association.out'), data.table = F)
    association_df <- select(association_df, ID, PWAS.Z, PWAS.P)
    association_df$WEIGHT <- sample_size[race][[1]]
    write.table(association_df, paste0('04_meta_analysis/T2D_', race, '.txt'), sep = '\t', quote = F, row.names = F)
}

# generate metal config files (metal_T2D_MVP.txt) manually
# meta for MVP 
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 04_meta_analysis/metal_T2D_MVP.txt')

# generate metal config files (metal_T2D_MVP_PAGE.txt) manually
# meta for MVP and PAGE 7 files 
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 04_meta_analysis/metal_T2D_MVP_PAGE.txt')

# generate metal config files (metal_T2D_MVP_PAGE_AA.txt) manually
# meta for MVP and PAGE AA files
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 04_meta_analysis/metal_T2D_MVP_PAGE_AA.txt')

# generate metal config files (metal_T2D_MVP_PAGE_CA.txt) manually
# meta for MVP and PAGE CA files
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 04_meta_analysis/metal_T2D_MVP_PAGE_CA.txt')

# generate metal config files (metal_T2D_MVP_PAGE_HA.txt) manually
# meta for MVP and PAGE HA files
system('/mnt/lvm_vol_1/hzhong/mwas/meta/generic-metal/metal 04_meta_analysis/metal_T2D_MVP_PAGE_HA.txt')

