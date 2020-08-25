require(cummeRbund)
require(pheatmap)
require(RColorBrewer)

#Get matrix of FPKM values
FPKMMatrix <- repFpkmMatrix(genes(cuff))

#Turn expression values into probabilities, FPKM matrix is subsetted to only look at samples we care about (indices may change)
probs <- makeprobs(log10(FPKMMatrix[,c(6:17)] + 1))
#Get J-S distances
JSDs <- as.matrix(JSdist(probs))
#Get distances for clustering
drows <- JSdist(probs)
dcols <- JSdist(probs)
#Plot
pheatmap(1 - JSDs, color = colorRampPalette(brewer.pal(9, 'Blues'))(100), border_color = 'white', clustering_distance_rows = drows, clustering_distance_cols = dcols, 
         annotation_row = test, treeheight_col=0)
