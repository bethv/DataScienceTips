
# combine top 10000 by max variance

load("datExprSigNoXYconsSamp.rda")
load("datExprCortNoXYconsSamp.rda")

c<-datExprSigNoXYconsSamp
b<-datExprCortNoXYconsSamp

names(c)<-paste("c",names(c),sep = "_")
names(b)<-paste("b",names(b),sep = "_")

tb<-as.data.frame(t(b))
tb$var<-apply(tb,1,var)
tb<-tb[order(tb$var,decreasing = TRUE),]
topVarB<-rownames(tb)[1:10000]
b<-b[,names(b)%in%topVarB]

tc<-as.data.frame(t(c))
tc$var<-apply(tc,1,var)
tc<-tc[order(tc$var,decreasing = TRUE),]
topVarC<-rownames(tc)[1:10000]
c<-c[,names(c)%in%topVarC]


datExpr<-cbind(b,c)

load("genedataAllcolon.rda")
geneData.c<-genedataAll
geneData.b<-genedataAll
geneData.b$Probe<-paste("b",geneData.b$Probe,sep = "_")
geneData.c$Probe<-paste("c",geneData.c$Probe,sep = "_")
geneData.b$Source<-"Brain"
geneData.c$Source<-"Colon"

geneData<-rbind(geneData.b,geneData.c)
geneData<-geneData[geneData$Probe%in%names(datExpr),]
save(datExpr,geneData,file = "BandCtopVar3.rda")



# choose softPower ------------------

load("BandCtopVar3.rda")

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file = "SigIndependenceBCtopVar3.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red");
abline(h=0.80,col="red")
abline(h=0.90,col="blue")
dev.off()

pdf("SigConnectivityBCtopVar3.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")
dev.off()




# Run WGCNA ----------------

load("BandCtopVar3.rda")

library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

#params-----
softPower = 12
netType = "signed"
minModuleSize = 80
deepSplit = 1
MEDissThres = 0.15

fname.root<-"sp12ms80ds1mch0.15BCmv3"

source("RunWGCNA.R")

save(datExpr,MEs,moduleLabels,moduleColors,geneTree, file = paste("BCNet_",fname.root,".rda",sep = ""))

source("gtexFunctions.R")


makeGeneInfo<-function(geneData,modules,geneModMemb){
    geneInfo0 = data.frame(ProbeID = geneData$Probe, 
                           Tissue = geneData$Source,
                           geneSymbol = geneData$Symbol,
                           EntrezID = geneData$EntrezID,
                           Name=geneData$Name,
                           moduleColor = modules)
    Modmemb<-numeric()
    for(i in 1:nrow(geneInfo0)){
        modColor<-modules[i]
        Modmemb[i]<-geneModMemb[i,modColor]
    }
    geneInfo0$MM<-Modmemb
    geneOrder = order(factor(geneInfo0$moduleColor,levels = names(geneModMemb)), -abs(geneInfo0$MM))
    
    geneInfo<-geneInfo0[geneOrder,]
    geneInfo<-geneInfo[geneInfo$geneSymbol!="",]
    geneInfo<-geneInfo[!is.na(geneInfo$EntrezID),]
    geneInfo$absMM<-abs(geneInfo$MM)
    return(geneInfo)
}

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
names(geneModuleMembership) = modNames
geneDataTopVar<-geneData[geneData$Probe%in%names(datExpr),]
geneInfo<-makeGeneInfo(geneDataTopVar,moduleColors,geneModuleMembership)
write.csv(as.data.frame(table(factor(moduleColors,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizesColon.csv",sep = ""))
write.csv(geneInfo, paste("BC_INFO_",fname.root,".csv",sep = ""))
save(geneInfo,file=paste("BC_INFO_",fname.root,".rda",sep = ""))
#GO
library(GOstats)
library(hgu95av2.db)

allLLIDs <- geneDataTopVar$EntrezID
allLLIDs<-allLLIDs[!allLLIDs%in%NA]
allLLIDs<-as.character(allLLIDs)
geneInfo<-geneInfo[(geneInfo$moduleColor!="grey")&(geneInfo$absMM>=0.6),]
modList<-split(geneInfo,geneInfo$moduleColor)
 
 
GOlist<-list()

for (module in names(modList)) {
  df1<-view.GOstats(modList[[module]],"BP")
  df1$moduleColor<-module
  df1$moduleSize<-nrow(modList[[module]])
  df1<-df1[order(df1$Pvalue),]
  GOlist[[module]]<-df1
}

save(GOlist,file = paste("BCGOlist_",fname.root,".rda",sep = ""))

for (module in names(GOlist)) {
  GOlist[[module]]
  GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]

for (module in names(GOlist)) {
  GOtable<-rbind(GOtable,GOlist[[module]])
}

GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste("BCGOTable_",fname.root,".csv",sep = ""))

load("BC_INFO_sp12ms80ds1mch0.15BCmv3.rda")
geneInfo0<-geneInfo
geneModuleMembership0<-geneModuleMembership[geneInfo$ProbeID,]
for (mod in modNames){
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership0[,mod])
    names(geneInfo0) = c(oldNames, paste("MM.", mod, sep=""))
}

# Order the genes first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0[,paste("GS.", names(trait), sep="")]));
geneInfo = geneInfo0[geneOrder, ]

geneInfoAllMMs<-geneInfo0
save(geneInfoAllMMs,file = "geneInfoAllMMs.rda")
write.csv(geneInfoAllMMs,"geneInfoAllMMs.csv")

table(geneInfo$moduleColor,geneInfo$Tissue)

ColonAll<-write.table(geneInfo$EntrezID[geneInfo$Tissue=="Colon"],row.names = F, quote = F, file = "ColonAll.txt")
BrainAll<-write.table(geneInfo$EntrezID[geneInfo$Tissue=="Brain"],row.names = F, quote = F, file = "BrainAll.txt")

for (Tissue in c("Colon","Brain")) {
    for (mod in modNames[!modNames%in%"grey"]) {
        write.table(geneInfo$EntrezID[geneInfo$Tissue==Tissue&geneInfo$moduleColor==mod],row.names = F, quote = F, file = paste(Tissue,"_",mod,".txt",sep = ""))
    }
}

pdf("PlotEigens.pdf",width = 7, height = 6.5)
plotEigengeneNetworks(MEs,setLabels = NULL,plotDendrograms = F, signed = T, setMargins = T)
dev.off()

geneInfo0.5<-geneInfo[geneInfo$MM>=0.5,]
table(geneInfo0.5$Tissue)
table(geneInfo0.5$Tissue,geneInfo0.5$moduleColor)

geneInfo0.5_B.bl.tu.br.gr.y0.8<-geneInfo0.5[!(geneInfo0.5$Tissue=="Brain"&geneInfo0.5$moduleColor%in%c("blue","turquoise","brown","green","yellow")&geneInfo0.5$MM<0.8),]
table(geneInfo0.5_B.bl.tu.br.gr.y0.8$Tissue)
table(geneInfo0.5_B.bl.tu.br.gr.y0.8$Tissue,geneInfo0.5_B.bl.tu.br.gr.y0.8$moduleColor)
ToCyto<-geneInfo0.5_B.bl.tu.br.gr.y0.8
probes = as.character(ToCyto$ProbeID)
rownames(ToCyto)<-probes
inToCyto<-names(datExpr)%in%probes
modTOM = TOM[inToCyto, inToCyto]
probes2<-names(datExpr)[inToCyto]
dimnames(modTOM) = list(probes2, probes2)
ToCyto<-ToCyto[probes2,]

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "edgesBC040318.txt",
                               nodeFile = "nodesBC040318.txt",
                               weighted = TRUE,
                               threshold = 0.05,
                               nodeNames = probes2,
                               altNodeNames = ToCyto$geneSymbol,
                               nodeAttr = ToCyto[,c("Tissue","moduleColor","EntrezID","MM","Name")])

save(cyt,ToCyto,modTOM,file="cytoscapeBC_040318.rda")


load("cytoscapeBC_040318.rda")

edges2<-cyt$edgeData    

FromC<-edges2$fromNode[grep("c_E",edges2$fromNode)]
ToC<-edges2$toNode[grep("c_E",edges2$toNode)]
CC<-edges2[edges2$fromNode%in%FromC&edges2$toNode%in%ToC,]
CC<-CC[order(CC$weight,decreasing = T),]

FromB<-edges2$fromNode[grep("b_E",edges2$fromNode)]
ToB<-edges2$toNode[grep("b_E",edges2$toNode)]
BB<-edges2[edges2$fromNode%in%FromB&edges2$toNode%in%ToB,]
BB<-BB[order(BB$weight,decreasing = T),]

BC<-edges2[(edges2$fromNode%in%FromB&edges2$toNode%in%ToC)|(edges2$fromNode%in%FromC&edges2$toNode%in%ToB),]
BC<-BC[order(BC$weight,decreasing = T),]



nodes2<-cyt$nodeData
excludeNode<-c(nodes2$nodeName[grep("LINC",nodes2$altName)], nodes2$nodeName[grep("open reading frame",nodes2$Name)])
nodesAllButExclude<-nodes2[!nodes2$nodeName%in%excludeNode,]

nodes2<-nodes2[nodes2$nodeName%in%c(BC$fromNode,BC$toNode),]

nodes3<-nodes2[nodes2$moduleColor%in%c("lightcyan", "tan", "brown", "turquoise", "greenyellow", "magenta","cyan","purple", "yellow"),]

BC$Net<-"BC"
BB$Net<-"BB"
CC$Net<-"CC"

BC<-BC[(BC$fromNode%in%nodes3$nodeName)|(BC$toNode%in%nodes3$nodeName),]
BC<-BC[BC$weight>=0.1,]

BB<-BB[(BB$fromNode%in%nodes3$nodeName)|(BB$toNode%in%nodes3$nodeName),]
BB<-BB[BB$weight>=0.3,]

CC<-CC[(CC$fromNode%in%nodes3$nodeName)|(CC$toNode%in%nodes3$nodeName),]
# CC<-CC[CC$weight>=0.08,]

edges3<-rbind(BC,CC,BB)
write.csv(edges3,"edges3.csv", quote = F, row.names = F)
write.csv(nodes3,"nodes3.csv", quote = F, row.names = F)

write.csv(nodes2,"nodes2.csv", quote = F, row.names = F)


nodesInf<-nodes3[nodes3$moduleColor%in%c("cyan","purple", "yellow"),]
BCInf<-BC[(BC$fromNode%in%nodesInf$nodeName)|(BC$toNode%in%nodesInf$nodeName),]
BCInf<-BCInf[BCInf$weight>=0.1,]

BBInf<-BB[(BB$fromNode%in%nodesInf$nodeName)|(BB$toNode%in%nodesInf$nodeName),]
BBInf<-BBInf[BBInf$weight>=0.1,]

CCInf<-CC[(CC$fromNode%in%nodesInf$nodeName)|(CC$toNode%in%nodesInf$nodeName),]
CCInf<-CCInf[CCInf$weight>=0.1,]

edgesInf<-rbind(BCInf,CCInf,BBInf)

edgesInf<-edgesInf[(edgesInf$fromNode%in%nodesAllButExclude$nodeName)|(edgesInf$toNode%in%nodesAllButExclude$nodeName),]


write.csv(edgesInf,"edgesInf.csv", quote = F, row.names = F)
write.csv(nodesAllButExclude,"nodesInf.csv", quote = F, row.names = F)

write.csv(nodesInf,"nodesInf.csv", quote = F, row.names = F)



nodesInf<-nodes3[nodes3$moduleColor%in%c("lightcyan", "tan", "brown", "turquoise", "greenyellow", "magenta","cyan","purple", "yellow"),]
