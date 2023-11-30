library(ggtree)
library(treeio)
library(TreeTools)
library(ape)
library(dendextend)
library(phangorn)
library(dplyr)

#set seed to avoid randomization
set.seed(5)

#set working directory
setwd(file.choose())

#import newwick dataset
allDncv <- read.tree(file = '../extraValidations/mafftOutput.nwk') #import newwick trees

#import annotation file for downstream labelling
allDataInfo <- read.table(file = 'all_data_dncv.txt', header = F, sep = '\t')
names(allDataInfo) <- c('ID', 'domain', 'supp_tab', 'phylum', 'DOM','strand', 'effector')

#pre-process data
clustDT <- as.ClusterTable(allDncv, tipLabels = allDncv$tip.label)
allDncvhc <- hclust(d = dist(clustDT), method = 'ward.D') #ward.D or complete will do
allDncvhc$height <- log10(allDncvhc$height+1)

#align names on tree file with names on annotation file
newOrder <- allDncvhc$order
newTips <- allDncv$tip.label[newOrder]
allDncvhc$labels <- newTips
tip.labs <- data.frame('ID' = allDncv$tip.label, 'metID' = allDncv$tip.label)
adjLabs <- tip.labs %>% inner_join(allDataInfo, by='ID') #(helps to join based on the order of dfs)


newPhylumOrder <- adjLabs$phylum[newOrder]
newDomOrder <- adjLabs$domain[newOrder]
neweffOrder <- adjLabs$effector[newOrder]
newTypeOrder <- adjLabs$supp_tab[newOrder]
newStrandOrder <- adjLabs$strand[newOrder]


#test.R
newIDnames <- adjLabs$ID[newOrder]
head(newIDnames) #line1
head(newTips) #line 2 
#line 1 $2 should be the same


#remove all labels except target DNCV and EC3234A

knockedOutLabs <- allDncvhc$labels
#NCBI names of DSM, EC3234 AND DNCV
indDSM <- which(knockedOutLabs == '2519180415')
indec3234 <- which(knockedOutLabs == '2624663363')
indDncv <- which(knockedOutLabs == '2601131560')
knockedOutLabs <- replace(knockedOutLabs, c(indDSM, indec3234, indDncv), c('DSM_DncV',"EC3234a_DncV", "EC_TW11681_DncV2"))

allDncvhc_no_tips <- allDncvhc
allDncvhc_no_tips$labels <- knockedOutLabs

targets <-  c('EC_TW11681_DncV2',
              'EC_TW11681_DncV1',
              'VC_DncV',
              'DSM_DncV',
              'EC3234a_DncV',
              'S_algicola_DncV'
)

na_labels <- c()
for (i in allDncvhc_no_tips$labels) {
  if (i %in% targets) {
    na_labels <- c(na_labels, i)
  } else {
    na_labels <- c(na_labels, '')
  }
}

na_labels_points <- c()
for (i in allDncvhc_no_tips$labels) {
  if (i %in% targets) {
    na_labels_points <- c(na_labels_points, "•")
  } else {
    na_labels_points <- c(na_labels_points, '')
  }
}


allDncvhc_no_tips$labels <- na_labels
allDncvhc_no_tips_points <- allDncvhc_no_tips
allDncvhc_no_tips_points$labels <- na_labels_points

#cutting into parts
allDncvDEND_no_tips <- as.dendrogram(allDncvhc_no_tips)
kkk <- cutree(allDncvDEND_no_tips, k = 10) #CHOSE K=10 BASED ON ELBOW TECHNIQUE OF CLUSTER SELECTION
cols2 <- c(2:7, 9:12)[kkk]


plot(as.phylo(allDncvhc_no_tips_points), cex = 2.5, tip.color = cols2, no.margin = F, 
     type = 'u', font = 2, edge.color = 'gray26', node.lty = 1, lab4ut = 'axial',
     rotate.tree = 0, use.edge.length =F, label.offset = 0.1, node.depth = 1 )

dev.off()



####--- VISUALIZATION OF THE FIVE ENZYMES ---####

fiveEnz <- read.tree(file.choose())

sss <- cutree(chronos(fiveEnz), k = 4)
cols3 <- c(3,1,4,2)[sss]

quartz(type = 'pdf', file = '../extraValidations/allOthers/sixEnzCladogram.pdf', dpi = 300, height = 5, width = 5)

plot(as.phylo(fiveEnz), use.edge.length = F, cex = 1, no.margin = F, font = 2, 
     edge.color = 'gray26', node.lty = 1, align.tip.label = F, edge.width = 2.5, 
     label.offset = 0.2, tip.color = cols3)

dev.off()
