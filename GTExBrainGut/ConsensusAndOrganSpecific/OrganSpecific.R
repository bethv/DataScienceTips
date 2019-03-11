


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

# colon
load("../datExprSigNoXYconsSamp.rda")
datExpr<-datExprSigNoXYconsSamp
source("RunWGCNA.R")
datExpr.C<-datExpr
MEs.C<-MEs
moduleLabels.C<-moduleLabels
moduleColors.C<-moduleColors
geneTree.C<-geneTree
save(datExpr.C,MEs.C,moduleLabels.C,moduleColors.C,geneTree.C, file = paste("ColonNet_",fname.root,".rda",sep = ""))
save(TOM,file = paste("ColonTOM_",fname.root,".rda",sep = ""))
rm(datExpr,TOM,dissTOM, geneTree,dynamicMods,dynamicColors,MEList,MEs,MEDiss,METree,merge,mergedColors,mergedMEs,moduleColors,colorOrder,moduleLabels)


# brain
load("../datExprCortNoXYconsSamp.rda")
datExpr<-datExprCortNoXYconsSamp
source("RunWGCNA.R")
datExpr.B<-datExpr
MEs.B<-MEs
moduleLabels.B<-moduleLabels
moduleColors.B<-moduleColors
geneTree.B<-geneTree
save(datExpr.B,MEs.B,moduleLabels.B,moduleColors.B,geneTree.B, file = paste("BrainNet_",fname.root,".rda",sep = ""))
save(TOM,file = paste("BrainTOM_",fname.root,".rda",sep = ""))
rm(datExpr,TOM,dissTOM, geneTree,dynamicMods,dynamicColors,MEList,MEs,MEDiss,METree, merge,mergedColors,mergedMEs,moduleColors,colorOrder,moduleLabels)

load("../genedataAllcolon.rda")
source("../gtexFunctions.R")

# colon
modNames = substring(names(MEs.C), 3)
geneModuleMembership = as.data.frame(cor(datExpr.C, MEs.C, use = "p"))
names(geneModuleMembership) = modNames
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(datExpr.C),]
geneInfo.C<-makeGeneInfo(genedataWGCNA,moduleColors.C,geneModuleMembership)
write.csv(as.data.frame(table(factor(moduleColors.C,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizesColon.csv",sep = ""))
write.csv(geneInfo.C, paste("ColonINFO_",fname.root,".csv",sep = ""))
save(geneInfo.C,file=paste("ColonINFO_",fname.root,".rda",sep = ""))

# brain 
modNames = substring(names(MEs.B), 3)
geneModuleMembership = as.data.frame(cor(datExpr.B, MEs.B, use = "p"))
names(geneModuleMembership) = modNames
genedataWGCNA<-genedataAll[genedataAll$Probe%in%names(datExpr.B),]
geneInfo.B<-makeGeneInfo(genedataWGCNA,moduleColors.B,geneModuleMembership)
write.csv(as.data.frame(table(factor(moduleColors.B,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizesBrain.csv",sep = ""))
write.csv(geneInfo.B, paste("BrainINFO_",fname.root,".csv",sep = ""))
save(geneInfo.B,file=paste("BrainINFO_",fname.root,".rda",sep = ""))

#GO
library(GOstats)
library(hgu95av2.db)

#colon----
genedataC<-genedataAll[genedataAll$Probe%in%names(datExpr.C),]
allLLIDsC <- genedataC$EntrezID
allLLIDsC<-allLLIDsC[!allLLIDsC%in%NA]
allLLIDsC<-as.character(allLLIDsC)
load(geneInfoColonData)
geneInfo.C<-geneInfo.C[(geneInfo.C$moduleColor!="grey")&(geneInfo.C$absMM>=0.6),]
modListC<-split(geneInfo.C,geneInfo.C$moduleColor)
allLLIDs<-allLLIDsC
modList<-modListC

GOlist<-list()

for (module in names(modList)) {
    df1<-view.GOstats(modList[[module]],"BP")
    df1$moduleColor<-module
    df1$moduleSize<-nrow(modList[[module]])
    df1<-df1[order(df1$Pvalue),]
    GOlist[[module]]<-df1
}

save(GOlist,file = paste("ColonGOlist_",fname.root,".rda",sep = ""))

for (module in names(GOlist)) {
    GOlist[[module]]
    GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]

for (module in names(GOlist)) {
    GOtable<-rbind(GOtable,GOlist[[module]])
}

GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste("ColonGOTable_",fname.root,".csv",sep = ""))

# brain----
genedataB<-genedataAll[genedataAll$Probe%in%names(datExpr.B),]
allLLIDsB <- genedataB$EntrezID
allLLIDsB<-allLLIDsB[!allLLIDsB%in%NA]
allLLIDsB<-as.character(allLLIDsB)
load(geneInfoBrainData)
geneInfo.B<-geneInfo.B[(geneInfo.B$moduleColor!="grey")&(geneInfo.B$absMM>=0.6),]
modListB<-split(geneInfo.B,geneInfo.B$moduleColor)

allLLIDs<-allLLIDsB
modList<-modListB

GOlist<-list()

for (module in names(modList)) {
    df1<-view.GOstats(modList[[module]],"BP")
    df1$moduleColor<-module
    df1$moduleSize<-nrow(modList[[module]])
    df1<-df1[order(df1$Pvalue),]
    GOlist[[module]]<-df1
}

save(GOlist,file = paste("BrainGOlist_",fname.root,".rda",sep = ""))

for (module in names(GOlist)) {
    GOlist[[module]]
    GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]
for (module in names(GOlist)) {
    GOtable<-rbind(GOtable,GOlist[[module]])
}
GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste("BrainGOTable_",fname.root,".rda",sep = ""))
