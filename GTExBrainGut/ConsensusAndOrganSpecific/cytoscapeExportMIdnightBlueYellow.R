load("ConsNet_sp18ms80ds1mch0.15_15mods_midnightblue.rda")
load("ConsNet_sp18ms80ds1mch0.15_15mods_midnightblue_INFO.rda")
load("ColonTOM_sp18ms80ds1mch0.15.rda")
gic2<-geneInfo.consBC[geneInfo.consBC$moduleColor%in%c("midnightblue","yellow"),]
gic2<-gic2[gic2$absMM.C>=0.7,]
probes = names(datExpr.C)
selected = probes%in%gic2$ProbeID
selectedProbes = probes[selected]
gic2<-gic2[match(selectedProbes,gic2$ProbeID),]
selectedTOM = TOM[selected, selected];
dimnames(selectedTOM) = list(selectedProbes, selectedProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(selectedTOM,
                               edgeFile = "Edges_midnightblue_yellow.txt",
                               nodeFile = "Nodes_midnightblue_yellow.txt",
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = selectedProbes,
                               altNodeNames = gic2$geneSymbol,
                               nodeAttr = gic2);


cyt2 = exportNetworkToCytoscape(selectedTOM,
                               edgeFile = "Edges_midnightblue_yellow2.txt",
                               nodeFile = NULL,
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = selectedProbes,
                               altNodeNames = gic2$moduleColor,
                               nodeAttr = gic2);

modColortable<-data.frame(sharedName=paste(cyt2$edgeData$fromNode," (interacts with) ",cyt2$edgeData$toNode),fromColor=cyt2$edgeData$fromAltName,toColor=cyt2$edgeData$toAltName)
write.csv(modColortable,"modcolortable.csv")


probes.C = as.character(CtopModsTop50$ProbeID)
rownames(CtopModsTop50)<-probes.C
inC<-names(datExpr.C)%in%probes.C
modTOM.C = TOM[inC, inC]
probes.C2<-names(datExpr.C)[inC]
dimnames(modTOM.C) = list(probes.C2, probes.C2)
CtopModsTop50<-CtopModsTop50[probes.C2,]

cyt.C = exportNetworkToCytoscape(modTOM.C,
                                 edgeFile = "Cedges102017.txt",
                                 nodeFile = "Cnodes102017.txt", 
                                 weighted = TRUE, 
                                 threshold = 0.05, 
                                 nodeNames = CtopModsTop50$geneSymbol,
                                 altNodeNames = probes.C2,
                                 nodeAttr = BtopModsTop20[,c("moduleColor","EntrezID","MM","Name")])

save(cyt.C,CtopModsTop50,modTOM.C,file="cytoscapeC_102017.rda")

