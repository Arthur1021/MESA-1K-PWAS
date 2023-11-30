library(data.table)
library(dplyr)
association_result <- fread('Figure_Tables/TableS4.txt', data.table = F)
association_result$P_extract <- NA
association_result$Z_extract <- NA
for (i in 1:nrow(association_result)){
    association_result[i,]$P_extract <- min(c(association_result[i,]$P_multi, association_result[i,]$P_AA, association_result[i,]$P_CA, association_result[i,]$P_EA, association_result[i,]$P_HA), na.rm = TRUE)
    association_result[i,]$Z_extract <- association_result[i,which(association_result[i,] ==  association_result[i,]$P_extract)[1] - 1]
}
association_result <- select(association_result, TargetFullName, P_extract, Z_extract)
association_result$direction <- ifelse(association_result$Z_extract < 0, 'negative', 'positive')
association_result$shape <- 'circle'
colnames(association_result)[1:3] <- c('ID', 'p', 'Z')



library(igraph)
library(reshape2)
library(dplyr)
library(ggplot2)

nodes <- association_result
nodes <- nodes[order(nodes$ID),]
nodes <- nodes[order(nodes$direction),]


links <- read.table('08_PPI_IPA/string_interactions_short.tsv',header = T,sep = '\t')
annotation <- fread('01_preprocess_data/annotation/SomaScan_annotation.txt', data.table = F)
annotation <- annotation[annotation$EntrezGeneSymbol %in% c(links$node1, links$node2),]
annotation <- annotation[!duplicated(annotation$EntrezGeneSymbol),]
annotation <- select(annotation, EntrezGeneSymbol, TargetFullName)
links <- select(links, node1, node2, combined_score)
links <- left_join(links, annotation, by = c('node1' = 'EntrezGeneSymbol'))
links <- left_join(links, annotation, by = c('node2' = 'EntrezGeneSymbol'))

node_degrees <- read.table('08_PPI_IPA/string_node_degrees.tsv',header = T,sep = '\t')
links$color <- ifelse(links$node1 %in% node_degrees[node_degrees$node_degree >= 15,]$node | links$node2 %in% node_degrees[node_degrees$node_degree > 10,]$node, 'lightblue', 'lightgrey')

links <- select(links, TargetFullName.x, TargetFullName.y, combined_score, color)
colnames(links)[1:2] <- c('node1', 'node2')

# draw figure 
library(igraph)
library(reshape2)
library(plyr)
library(ggplot2)

mycol <- c("#EA4335", "#34A853")

# get uniq nodes
nodes <- nodes[!duplicated(nodes$ID),]

net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 

V(net)$color <- revalue(nodes$direction,c("positive"=mycol[1],"negative"=mycol[2]))
V(net)$size <- (-log10(V(net)$p))^0.5*3
V(net)$label.cex <- 2
E(net)$width <- E(net)$combined_score^2*10

pdf("08_PPI_IPA/protein_protein_interact.pdf", width = 30, height = 20 )
plot(net,layout=layout_in_circle,edge.curved=.3,edge.color=E(net)$color,vertex.label.family="Helvetica",
     vertex.label.font=2,vertex.shape="circle",vertex.frame.color = NA,vertex.color = adjustcolor(V(net)$color , alpha = 0.6))    #,vertex.shape=nodes$shape                                  

legend("topright", #图例的位置
       c("Z score > 0", "Z score < 0"),
       pch=21, col=NA, pt.bg=mycol, pt.cex=3,
       cex=2, bty="n", ncol=1)
       
f <- c(0.05, 1e-5,1e-10, 1e-20, 1e-26)
s <- (-log10(f))^0.5*3
legend("bottomright", 
       inset=c(0,0), #向下移
       legend=f, text.width = 0.2, 
       title = "P value", #title.adj = -.5,
       pch=21, pt.cex=s, bty='n',
       horiz = T, #横向排列
       cex = 2,
       col = "black")

l <- c(0.4, 0.6, 0.8, 1)
l2 <- l^2*10
legend("bottomright", 
       inset=c(0.5,0), #向下移
       legend=l, text.width = 0.2, 
       title = "score", #title.adj = -.5,
       cex = 2, lty = 1,lwd = l2, bty='n',
       horiz = T, #横向排列
       col = "black")
dev.off()


