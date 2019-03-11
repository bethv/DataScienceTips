library('recount')
load("rse_gene_blood.Rdata")
load("BandCtopVar3.rda")

subs<-rownames(datExpr)

head(colData(rse_gene))

rse_gene$person<-gsub("GTEX-","",rse_gene$sampid)
rse_gene$person<-gsub("-.+","",rse_gene$person)
head(colData(rse_gene))

bcfilt<-rse_gene$person%in%subs
rse_gene2<-rse_gene[,bcfilt]
dim(rse_gene2) #22

rseBlood <- scale_counts(rse_gene2)
save(rseBlood,file = "rseBlood.rda")


load("rse_gene_lung.Rdata")
rse_gene$person<-gsub("GTEX-","",rse_gene$sampid)
rse_gene$person<-gsub("-.+","",rse_gene$person)

bcfilt<-rse_gene$person%in%subs
rse_gene2<-rse_gene[,bcfilt]
dim(rse_gene2) #19

rseLung <- scale_counts(rse_gene2)
save(rseLung,file = "rseLung.rda")

library('DESeq2')

dds <- DESeqDataSet(rseLung, ~as.numeric(smrin))
dds <- DESeq(dds)
nrow(dds)
ncol(dds)

n30percent<-.3*dim(rseLung)[2]
countsGr10in30per<-apply(counts(dds),1,function(x){
        x2<-ifelse(x>10,1,0)
        x3<-ifelse(sum(x2)>n30percent,TRUE,FALSE)
        return(x3)
})

dds <- dds[ countsGr10in30per, ]
nrow(dds) #25031
vsd <- vst(dds, blind = FALSE)
edata=as.data.frame(assay(vsd))
colramp = colorRampPalette(c(3,"white",2))(20)
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(edata[,i]),lwd=3,col=colramp[i])}

library(limma)
norm_edata = normalizeQuantiles(as.matrix(edata))
plot(density(norm_edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.30))
for(i in 2:73){lines(density(norm_edata [,i]),lwd=3,col=colramp[i])}

datExpr<-as.data.frame(t(norm_edata))
rownames(datExpr)<-rseLung$person

save(dds,file = "ddsLung.rda")
save(vsd,file = "vsdLung.rda")
save(edata,file = "edataLung.rda")
save(norm_edata,file = "norm_edataLung.rda")
save(datExpr,file = "datExprLung.rda")

ensemblIDs <- sapply( strsplit( rowData(rseLung)$gene_id, split="\\." ), "[", 1 )

library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene", "hgnc_symbol", "description","chromosome_name"),
                  filters = "ensembl_gene_id",
                  values = ensemblIDs,
                  mart = ensembl )
idx <- match( ensemblIDs, genemap$ensembl_gene_id )
load("genedataAllcolon.rda")

genedataAll<-data.frame(Probe=genedataAll$Probe, ENS=ensemblIDs,EntrezID=genemap$entrezgene[idx],Symbol=genemap$hgnc_symbol[idx], Name=genemap$description[idx],Chr=genemap$chromosome_name[idx])

save(genedataAll, file = "genedataAllLung.rda")

noXY<-as.character(genedataAll$Probe[!genedataAll$Chr%in%c("X","Y")])

datExprNoXY<-datExpr[,names(datExpr)%in%noXY]
save(datExprNoXY,file = "datExprNoXYLung.rda")
save(noXY,file = "noXYLung.rda")


load("datExprCortNoXYconsSamp.rda")

l<-datExprNoXY
b<-datExprCortNoXYconsSamp[rownames(datExprNoXY),]

names(l)<-paste("l",names(l),sep = "_")
names(b)<-paste("b",names(b),sep = "_")

tb<-as.data.frame(t(b))
tb$var<-apply(tb,1,var)
tb<-tb[order(tb$var,decreasing = TRUE),]
topVarB<-rownames(tb)[1:10000]
b<-b[,names(b)%in%topVarB]

tl<-as.data.frame(t(l))
tl$var<-apply(tl,1,var)
tl<-tl[order(tl$var,decreasing = TRUE),]
topVarL<-rownames(tl)[1:10000]
l<-l[,names(l)%in%topVarL]


datExpr<-cbind(b,l)

geneData.l<-genedataAll
geneData.b<-genedataAll
geneData.b$Probe<-paste("b",geneData.b$Probe,sep = "_")
geneData.l$Probe<-paste("l",geneData.l$Probe,sep = "_")
geneData.b$Source<-"Brain"
geneData.l$Source<-"Lung"

geneData<-rbind(geneData.b,geneData.l)
geneData<-geneData[geneData$Probe%in%names(datExpr),]
save(datExpr,geneData,file = "BandLtopVar.rda")

load("BandLtopVar.rda")
library(WGCNA)
options(stringsAsFactors = FALSE)
allowWGCNAThreads() 

#params-----
softPower = 12
netType = "signed"
minModuleSize = 80
deepSplit = 1
MEDissThres = 0.15

fname.root<-"BLsp12ms80ds1mch0.15mv"

source("RunWGCNA.R")

save(datExpr,MEs,moduleLabels,moduleColors,geneTree, file = paste("BLNet_",fname.root,".rda",sep = ""))

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
write.csv(as.data.frame(table(factor(moduleColors,levels = names(geneModuleMembership)))),paste(fname.root,"_ModSizesBL.csv",sep = ""))
write.csv(geneInfo, paste("BL_INFO_",fname.root,".csv",sep = ""))
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

save(GOlist,file = paste("BLGOlist_",fname.root,".rda",sep = ""))

for (module in names(GOlist)) {
        GOlist[[module]]
        GOlist[[module]]$TermPval<-paste(GOlist[[module]]$Term," (",format(GOlist[[module]]$Pvalue,digits = 2,scientific = TRUE),")",sep = "")
}

GOtable<-GOlist$blue[0,]

for (module in names(GOlist)) {
        GOtable<-rbind(GOtable,GOlist[[module]])
}

GOtable<-GOtable[GOtable$Size>=10,]
write.csv(GOtable,paste("BLGOTable_",fname.root,".csv",sep = ""))

geneInfo0.6<-geneInfo[geneInfo$MM>=0.6,]

write.csv(addmargins(table(geneInfo0.6$moduleColor,geneInfo0.6$Tissue)),"BLTissueDist.csv")



ColonAll<-write.table(geneInfo$EntrezID[geneInfo$Tissue=="Colon"],row.names = F, quote = F, file = "BLColonAll.txt")
BrainAll<-write.table(geneInfo$EntrezID[geneInfo$Tissue=="Brain"],row.names = F, quote = F, file = "BLBrainAll.txt")

for (Tissue in c("Colon","Brain")) {
        for (mod in modNames[!modNames%in%"grey"]) {
                write.table(geneInfo$EntrezID[geneInfo$Tissue==Tissue&geneInfo$moduleColor==mod],row.names = F, quote = F, file = paste("BL_",Tissue,"_",mod,".txt",sep = ""))
        }
}

pdf("PlotEigensBL.pdf",width = 7, height = 6.5)
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



geneInfo0.6$moduleColor<-factor(geneInfo0.6$moduleColor,levels = modNames)
geneInfo0.6<-geneInfo0.6[order(geneInfo0.6$moduleColor),]




