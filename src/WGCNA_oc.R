############
library(sva)
library(Mus.musculus)#seems to use mm10
library(WGCNA)
library(topGO)
library(DESeq2)
library("igraph")
library("bnlearn")
library("parallel")
library(rtracklayer)
############
options(stringsAsFactors = FALSE)
set.seed(8675309)
#

# read in the RNA-seq processed counts file, Oc
# read in the RNA-seq processed counts file

counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_37samps_tpm_over0.1_37samps_COUNTS_oc.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)


#match ID keys from psomagen to mouse IDs
keys = read.csv("./data/ID Key_ Psomagen BM Ob_Oc RNA seq.csv")
keys$mouse.ID = as.character(keys$mouse.ID)
keys = keys[keys$cell.type == "osteoclast",]

keys$Psomagen.ID = sapply(strsplit(keys$Psomagen.ID, " "),"[",3)

#edit colnames to be only the number
colnames(counts) = sapply(strsplit(colnames(counts), "-"),"[",3)
colnames(counts) = sapply(strsplit(colnames(counts), "_"),"[",1)

colnames(counts) = keys[match(colnames(counts), keys$Psomagen.ID),"mouse.ID"]

#add ".1" to some colnames so they match the GigaMUGA data

colnames(counts)[which(colnames(counts) == 324)] = "324.1"

#
annot_file = read.csv("./results/flat/annot_file_oc.csv")
annot_file = annot_file[,c(2,3)]

counts = counts[which(rownames(counts) %in% annot_file$Gene.ID),]
#

#find and remove features that have fewer than 10 reads in more than 90% (163) of samples 
x=c()
for(i in 1:nrow(counts)){
  if(sum(counts[i,]<10)>=163){
    print(i)
    x = append(x,i)
  }
  
}

#184 genes removed
counts = counts[-x,]


#
#read in covars
#read in the cross file
load("./results/Rdata/cross_basic_cleaned.Rdata")

x=cross_basic$covar



#get covs by matching colnames of vst with mouse ID in raw pheno file
covs = x[match(colnames(counts),rownames(x)),]


#sex, age at sac, generation
covs = covs[,c(2,3,6)]

#potential batch data for ob/oc
plate_batch = read.csv("data/pheno_data/osteoblast_data.csv")

plate_batch[which(plate_batch$sample == 324),"sample"] = "324.1"

plate_batch = plate_batch[which(plate_batch$sample %in% colnames(counts)),]

plate_batch$D0_batch = as.numeric(as.factor(plate_batch$D0.date))

covs = merge(covs, plate_batch, by.x = 0, by.y="sample") #check
rownames(covs) = covs$Row.names



covs = covs[match(colnames(counts),rownames(covs)),]

covs = covs[,c(2,3,4,31)]

covs$D0_batch = as.factor(covs$D0_batch)
####################################

batch = covs$D0_batch
adjusted <- ComBat_seq(as.matrix(counts), batch=batch, group=NULL)
#


#vst from deseq2
vst = DESeq2::varianceStabilizingTransformation(as.matrix(adjusted))



#transpose the matrix
edata_oc = t(vst)

save(edata_oc, file = "./results/Rdata/networks/edata_oc.Rdata")
load("./results/Rdata/networks/edata_oc.Rdata")
#check data
gsg = goodSamplesGenes(edata_oc, verbose = 3);
gsg$allOK

##
#cluster samples
sampleTree = hclust(dist(edata_oc), method = "average")

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)
#####



#pick the soft thresholding power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(edata_oc, powerVector = powers, verbose = 5,networkType = "signed", dataIsExpr = TRUE)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
###


#make network with sft=12
net_oc = blockwiseModules(edata_oc, power = 12,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE,
                       saveTOMs = FALSE,
                       verbose = 3)

table(net_oc$colors)
## 12: 20 modules not including 0, 11381 genes in module 0

#saveRDS(net_oc, file="./results/Rdata/networks/wgcna_oc_12.RDS")
net_oc = readRDS("./results/Rdata/networks/wgcna_oc_12.RDS")

#get traits we want to look at
pheno = read.csv("./results/flat/full_pheno_table.csv", stringsAsFactors = FALSE)


datTraits = pheno[which(pheno$Mouse.ID %in% rownames(edata_oc)),]
datTraits = datTraits[match(rownames(edata_oc),datTraits$Mouse.ID),]

#remove non-pheno columns
datTraits = datTraits[,-c(1:7,16,22)]



### NORMALIZE TRAITS?###






#convert to numeric
for(i in 1:ncol(datTraits)){datTraits[,i] = as.numeric(datTraits[,i])}

#remove histo columns 
#datTraits = datTraits[,-c(26:43)]
#remove adipose
#datTraits = datTraits[,-c(44:55)]

sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net_oc$colors)


# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net_oc$dendrograms[[1]], mergedColors[net_oc$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


# Construct numerical labels corresponding to the colors
moduleLabels = net_oc$colors
moduleColors = labels2colors(net_oc$colors)
MEs = net_oc$MEs;
geneTree = net_oc$dendrograms[[1]]
nGenes = ncol(edata_oc);
nSamples = nrow(edata_oc);

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(edata_oc, moduleColors)$eigengenes
MEs = orderMEs(MEs0)



#REMOVE GREY
MEs = MEs[,-which(colnames(MEs) == "MEgrey")]
#cor module eigengenes with traits
moduleTraitCor = cor(MEs, datTraits, use = "p",method = "s");

moduleTraitPvalue = as.data.frame(matrix(nrow = nrow(moduleTraitCor),ncol = ncol(moduleTraitCor)))
for(i in 1:ncol(moduleTraitCor)){
  nSamples = length(which(is.na(datTraits[,colnames(moduleTraitCor)[i]]) == FALSE))
  moduleTraitPvalue[,i] = corPvalueStudent(moduleTraitCor[,i], nSamples) # uses sample for each trait. is this correct?
  print(colnames(moduleTraitCor)[i])
  print(nSamples)
  
}
nSamples=nrow(edata_oc)
colnames(moduleTraitPvalue) = colnames(moduleTraitCor)
rownames(moduleTraitPvalue) = rownames(moduleTraitCor)

##remove everything but uCT
#moduleTraitPvalue = moduleTraitPvalue[,c(44:61)]
#moduleTraitCor = moduleTraitCor[,c(44:61)]

#moduleTraitPvalue = moduleTraitPvalue[,c(11:61)]
#moduleTraitCor = moduleTraitCor[,c(11:61)]

#moduleTraitPvalue = as.matrix(moduleTraitPvalue)

#keep all but  weight, length, glucose , fat pads, muscle masses
#moduleTraitPvalue = moduleTraitPvalue[,c(11:65)]
#moduleTraitCor = moduleTraitCor[,c(11:65)]


moduleTraitPvalue = as.matrix(moduleTraitPvalue)

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_full_4.RData")
save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_full_4.RData")

#identify significant modules
sig_mod = moduleTraitPvalue[which(rownames(moduleTraitPvalue) %in% names(which(apply(moduleTraitPvalue, 1, function(r) any(r < 0.05/ncol(MEs)))))),]

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = signif(moduleTraitPvalue, 1)

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

trait_names = gsub(pattern = "bending_",replacement = "",x = names(datTraits))
trait_names = gsub(pattern = "uCT_",replacement = "",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names = gsub(pattern = "\\.\\.",replacement = "\\.",x = trait_names)
trait_names[8] = "Adiposity"

####
####
#only trabecular traits
moduleTraitCor = moduleTraitCor[,c(34:41)]

trait_names = colnames(moduleTraitCor)
trait_names = gsub(pattern = "uCT_",replacement = "",x = trait_names)
#MEs = MEs[,-37]
textMatrix = signif(moduleTraitPvalue, 1)
####
####

#dim(textMatrix) = dim(moduleTraitCor)
textMatrix = textMatrix[,c(34:41)]
# Display the correlation values within a heatmap plot
#correlation is color, text is p-value
par(mar=c(4.1, 13.1, 3.1, 2.1))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = trait_names,
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               cex.text = 0.6,
               zlim = c(-1,1),
               cex.lab.x = 0.85,
               cex.lab.y = 0.7,
               verticalSeparator.x = c(1:length(trait_names)),
               horizontalSeparator.y = c(1:length(names(MEs))),
               main = paste("Module-trait relationships"))
###
which(moduleTraitCor == sort((moduleTraitCor),decreasing =TRUE)[1], arr.ind = T)
colnames(moduleTraitCor)[25]


#create trimmed expression object to remove grey module and get gene names from annot file
########################
combat_annot = as.data.frame(colnames(edata_oc))

combat_annot$module = net_oc$colors
combat_annot$color = moduleColors

#REMOVE GREY
rmv = combat_annot[which(combat_annot$color == "grey"),"colnames(edata_oc)"]

combat_annot = combat_annot[-which(combat_annot$`colnames(edata_oc)` %in% rmv),]
edata_oc_trim = edata_oc[,-(which(colnames(edata_oc) %in% rmv))]

combat_annot = combat_annot[,-2]

combat_annot[,c(3:4)] = annot_file[match(combat_annot$`colnames(edata_oc)`,annot_file$Gene.ID),c(1,2)]


geneModuleMembership = as.data.frame(cor(edata_oc_trim, MEs, use = "p")) # spearman or kendall? using pearson here

MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))#nsamples - Here it is RNA samples. Is This Correct? or use number of modules? or number of genes?
geneModuleMembership$gene = colnames(edata_oc_trim)


combat_annot[5:(ncol(geneModuleMembership)+4)] = geneModuleMembership[match(geneModuleMembership$gene,combat_annot$`colnames(edata_oc)`),]

x = match(colnames(edata_oc_trim),annot_file$Gene.ID)

colnames(edata_oc_trim) = annot_file$Gene.Name[x]


save(combat_annot, file = "./results/Rdata/networks/geneModMemAnnot_oc_p12.RData")

save(moduleTraitPvalue, file = "./results/Rdata/networks/moduleTraitPvalue_oc_full_12.RData")

save(moduleTraitCor, file = "./results/Rdata/networks/moduleTraitCor_oc_full_12.RData")

save(edata_oc_trim, file = "./results/Rdata/networks/edata_oc_trim_12.RData")


moduleColors = moduleColors[-which(moduleColors=="grey")]
hubs = chooseTopHubInEachModule(edata_oc,moduleColors)

#Do GO analysis
##
load("./results/Rdata/networks/geneModMemAnnot_oc_p12.RData")

network_GO = list()
for(color in unique(combat_annot$color)){
  
  allGenes = combat_annot$Gene.Name
  interesting.genes = combat_annot[which(combat_annot$color == color),"Gene.Name"]
  
  geneList<-factor(as.integer(allGenes %in% interesting.genes)) #If TRUE returns 1 as factor, otherwise 0
  names(geneList)<-allGenes
  ###MF###
  GOdata <- new("topGOdata", ontology = "MF", allGenes =geneList,
                annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
  test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
  result<-getSigGroups(GOdata,test.stat)
  t1<-GenTable(GOdata, classic=result, topNodes=length(result@score))
  head(t1)
  ###CC###
  GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
                annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
  test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
  result<-getSigGroups(GOdata,test.stat)
  t2<-GenTable(GOdata, classic=result, topNodes=length(result@score))
  head(t2)
  ###BP###
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.org, mapping='org.Mm.eg.db', ID='symbol')
  test.stat<-new("classicCount", testStatistic = GOFisherTest, name='Fisher test')
  result<-getSigGroups(GOdata,test.stat)
  t3<-GenTable(GOdata, classic=result, topNodes=length(result@score))
  head(t3)
  ####
  t.all = NULL
  t.all<-rbind(t1,t2,t3)
  t.all$classic<-as.numeric(as.character(t.all$classic))
  ######
  network_GO[[color]] = t.all
}
save(network_GO, file="./results/Rdata/networks/GO_oc_sft12.RData")
#####






#combine all GO results
#S6 - top GO terms for each network, sorted by pval
library(GO.db)
goterms = Term(GOTERM)

load("./results/Rdata/networks/GO_oc_sft12.RData")

networks = as.data.frame(matrix(nrow=10000, ncol = 10000))
nrow_counter = c()
net_counter = 1
col_counter = 1
for (i in 1:length(network_GO)){
  n = network_GO[[net_counter]]
  n = n[order(n$classic,decreasing = F),]
  #n = n[which(n$classic <=0.05),]
  
  x = n$GO.ID
  net_terms = unname(goterms[match(x, names(goterms))])
  
  networks[1:nrow(n),col_counter] = net_terms
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_Oc_GO_term")
  
  col_counter = col_counter + 1
  
  networks[1:nrow(n),col_counter] = n$classic
  colnames(networks)[col_counter] = paste0(names(network_GO)[net_counter], "_Oc_p_value")
  nrow_counter = append(nrow_counter, nrow(n))
  
  col_counter=col_counter+1
  net_counter = net_counter + 1
}


networks = networks[c(1:max(nrow_counter)),c(1:40)]



#write.csv(networks, file = "~/Desktop/nat_com_revs/supp/S6_rev.csv", row.names = F,na = "-")



