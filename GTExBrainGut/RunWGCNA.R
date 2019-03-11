

adjacency = adjacency(datExpr, power = softPower, type = netType);
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM
rm(adjacency)
geneTree = hclust(as.dist(dissTOM), method = "average");

dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = deepSplit, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
dynamicColors = labels2colors(dynamicMods)
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-abs(cor(MEs));
METree = hclust(as.dist(MEDiss), method = "average");
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 0)
mergedColors = merge$colors;
mergedMEs = merge$newMEs
moduleColors = mergedColors
colorOrder = c("grey", standardColors(150));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs