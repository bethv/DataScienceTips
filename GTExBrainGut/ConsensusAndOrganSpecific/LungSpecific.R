

# WGNA setup-----
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

#params-----
softPower = 18
netType = "signed"
minModuleSize = 80
deepSplit = 1
MEDissThres = 0.15

fname.root<-"sp18ms80ds1mch0.15"

load("datExprNoXYLung.rda")

datExpr<-datExprNoXY
source("RunWGCNA.R")
datExpr.L<-datExpr
MEs.L<-MEs
moduleLabels.L<-moduleLabels
moduleColors.L<-moduleColors
geneTree.L<-geneTree
save(datExpr.L,MEs.L,moduleLabels.L,moduleColors.L,geneTree.L, file = paste("LungNet_",fname.root,".rda",sep = ""))
save(TOM,file = paste("ColonTOM_",fname.root,".rda",sep = ""))
rm(datExpr,TOM,dissTOM, geneTree,dynamicMods,dynamicColors,MEList,MEs,MEDiss,METree,merge,mergedColors,mergedMEs,moduleColors,colorOrder,moduleLabels)


load("genedataAllcolon.rda")
source("gtexFunctions.R")

modNames = substring(names(MEs.L), 3)
geneModuleMembership = as.data.frame(cor(datExpr.L, MEs.L, use = "p"))
names(geneModuleMembership) = modNames
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(datExpr.L),]
geneInfo.L<-makeGeneInfo(genedataWGCNA,moduleColors.L,geneModuleMembership)
write.csv(as.data.frame(table(factor(moduleColors.L,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizesLung.csv",sep = ""))
write.csv(geneInfo.L, paste("LungINFO_",fname.root,".csv",sep = ""))
save(geneInfo.L,file=paste("LungINFO_",fname.root,".rda",sep = ""))

#GO
library(GOstats)
library(hgu95av2.db)

genedataL<-genedataAll[genedataAll$Probe%in%names(datExpr.L),]
allLLIDsL <- genedataL$EntrezID
allLLIDsL<-allLLIDsL[!allLLIDsL%in%NA]
allLLIDsL<-as.character(allLLIDsL)
geneInfo.L<-geneInfo.L[(geneInfo.L$moduleColor!="grey")&(geneInfo.L$absMM>=0.6),]
modListL<-split(geneInfo.L,geneInfo.L$moduleColor)
allLLIDs<-allLLIDsL
modList<-modListL

GOlist<-list()

for (module in names(modList)) {
    df1<-view.GOstats(modList[[module]],"BP")
    df1$moduleColor<-module
    df1$moduleSize<-nrow(modList[[module]])
    df1<-df1[order(df1$Pvalue),]
    GOlist[[module]]<-df1
}

save(GOlist,file = paste("LungGOlist_",fname.root,".rda",sep = ""))

for (module in names(GOlist)) {
    GOlist[[module]]
    GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]

for (module in names(GOlist)) {
    GOtable<-rbind(GOtable,GOlist[[module]])
}

GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste("LungGOTable_",fname.root,".csv",sep = ""))

topL<-c("white","darkgrey","midnightblue","black","brown","greenyellow",
        "green","yellow","darkorange","lightcyan","blue","darkturquoise")

cslistL<-modListL[names(modListL)%in%topL]

LtopModsMM0.8<-do.call(rbind,cslistL)

LtopModsMM0.8<-LtopModsMM0.8[abs(LtopModsMM0.8$MM)>=0.8,]

probes.L = as.character(LtopModsMM0.8$ProbeID)
rownames(LtopModsMM0.8)<-probes.L
inL<-names(datExpr.L)%in%probes.L
modTOM.L = TOM[inL, inL]
probes.L2<-names(datExpr.L)[inL]
dimnames(modTOM.L) = list(probes.L2, probes.L2)
LtopModsMM0.8<-LtopModsMM0.8[probes.L2,]
LtopModsMM0.8$Tissue<-"L"
LtopModsMM0.8$absMM<-abs(LtopModsMM0.8$MM)
cyt.L = exportNetworkToCytoscape(modTOM.L,
                                 edgeFile = "L08edges052018.txt",
                                 nodeFile = "L08nodes052018.txt", 
                                 weighted = TRUE, 
                                 threshold = 0.05, 
                                 nodeNames = LtopModsMM0.8$geneSymbol,
                                 altNodeNames = probes.L2,
                                 nodeAttr = LtopModsMM0.8[,c("moduleColor","EntrezID","MM","absMM")])

save(cyt.L,LtopModsMM0.8,modTOM.L,file="cytoscapeL_MM08_052018.rda")



golistL100<-lapply(cslistC,function(x){
    if(nrow(x)>100){
        x<-x[1:100,]
    }else{x<-x}
    return(x)
})

golistB100<-lapply(cslistb,function(x){
    if(nrow(x)>100){
        x<-x[1:100,]
    }else{x<-x}
    return(x)
})


save(cslistb,cslistC,file = "modlistsTopmods.rda")
save(golistB100,golistC100,file="modliststopmodstop100.rda")

probes.C = as.character(CtopModsMM0.8$ProbeID)
rownames(CtopModsMM0.8)<-probes.C
inC<-names(datExpr.C)%in%probes.C
load("SigNoxyConsSampe101917_TOM.rda")
modTOM.C = TOM[inC, inC]
probes.C2<-names(datExpr.C)[inC]
dimnames(modTOM.C) = list(probes.C2, probes.C2)
CtopModsMM0.8<-CtopModsMM0.8[probes.C2,]
CtopModsMM0.8$Tissue<-"C"
CtopModsMM0.8$MEname<-paste(CtopModsMM0.8$Tissue, CtopModsMM0.8$moduleColor,sep=".")
CtopModsMM0.8$MMsign<-ifelse(CtopModsMM0.8$MM>0,"Pos","Neg")
CtopModsMM0.8$absMM<-abs(CtopModsMM0.8$MM)
cyt.C = exportNetworkToCytoscape(modTOM.C,
                                 edgeFile = "C08edges102117.txt",
                                 nodeFile = "C08nodes102117.txt", 
                                 weighted = TRUE, 
                                 threshold = 0.05, 
                                 nodeNames = CtopModsMM0.8$geneSymbol,
                                 altNodeNames = probes.C2,
                                 nodeAttr = CtopModsMM0.8[,c("moduleColor","EntrezID","MM","Name","Tissue","MEname","MMsign","absMM")])

save(cyt.C,CtopModsMM0.8,modTOM.C,file="cytoscapeC_MM08_102117.rda")
rm(TOM)


topL<-

cslistL<-modListL[!names(modListL)%in%c("grey")]

LtopModsMM0.8<-do.call(rbind,cslistL)

LtopModsMM0.8<-LtopModsMM0.8[abs(LtopModsMM0.8$MM)>=0.8,]

probes.L = as.character(LtopModsMM0.8$ProbeID)
rownames(LtopModsMM0.8)<-probes.L
inL<-names(datExpr.L)%in%probes.L
modTOM.L = TOM[inL, inL]
probes.L2<-names(datExpr.L)[inL]
dimnames(modTOM.L) = list(probes.L2, probes.L2)
LtopModsMM0.8<-LtopModsMM0.8[probes.L2,]
LtopModsMM0.8$Tissue<-"L"
LtopModsMM0.8$absMM<-abs(LtopModsMM0.8$MM)
cyt.L = exportNetworkToCytoscape(modTOM.L,
                                 edgeFile = "L08edges052018_2.txt",
                                 nodeFile = "L08nodes052018_2.txt", 
                                 weighted = TRUE, 
                                 threshold = 0.09, 
                                 nodeNames = LtopModsMM0.8$geneSymbol,
                                 altNodeNames = probes.L2,
                                 nodeAttr = LtopModsMM0.8[,c("moduleColor","EntrezID","MM","absMM")])

save(cyt.L,LtopModsMM0.8,modTOM.L,file="cytoscapeL_MM08_052018_2.rda")
