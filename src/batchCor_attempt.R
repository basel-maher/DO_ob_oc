# sva seq: adjust with vst
# then in model have sex and age and ngen, remoe all other covars
# this should return svs not due to sex and age and ngen
# 
# in combatseq, use raw, but in batch add sex and agge so they are corrected for
# then adjust with vst for downstream

##################################################################
############
library(sva)
library(FactoMineR)
library(factoextra)
library(DESeq2)
############
options(stringsAsFactors = FALSE)
set.seed(8675309)
#

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



######


#vst from deseq2
vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))

#
#read in covars
#read in the cross file
load("./results/Rdata/cross_basic_cleaned.Rdata")

x=cross_basic$covar



#get covs by matching colnames of vst with mouse ID in raw pheno file
covs = x[match(colnames(vst),rownames(x)),]


#sex, age at sac, generation
covs = covs[,c(2,3,6)]

#potential batch data for ob/oc
plate_batch = read.csv("data/pheno_data/osteoblast_data.csv")

plate_batch[which(plate_batch$sample == 324),"sample"] = "324.1"

plate_batch = plate_batch[which(plate_batch$sample %in% colnames(counts)),]

plate_batch$D0_batch = as.numeric(as.factor(plate_batch$D0.date))

covs = merge(covs, plate_batch, by.x = 0, by.y="sample") #check
rownames(covs) = covs$Row.names



covs = covs[match(colnames(vst),rownames(covs)),]

covs = covs[,c(2,3,4,31)]

covs$D0_batch = as.factor(covs$D0_batch)
####################################
# #wont let me use ngen, too many variables or collinear
# mod = model.matrix(~as.factor(sex) + as.factor(age_at_sac_days), data=covs)
# mod0 = cbind(mod[,1])
# 
# n.sv = num.sv(vst,mod,method="be")
# 
# svseq = svaseq(vst,mod,mod0,n.sv=1)$sv
# 
# 
# batchCombat = model.matrix(~as.factor(sex) + as.factor(age_at_sac_days), data=covs)
# batchCombat = batchCombat[,-1]
# 
# 
# 
# #VST/PCA for using this data
# c = sapply(strsplit(colnames(counts), "-"),"[",3)
# c = sapply(strsplit(c, "_"),"[",1)
# colnames(counts) = c
# 
# colnames(counts)[which(colnames(counts) == 265)] = "265.1"
# colnames(counts)[which(colnames(counts) == 271)] = "271.1"
# colnames(counts)[which(colnames(counts) == 291)] = "291.1"
# colnames(counts)[which(colnames(counts) == 324)] = "324.1"
# colnames(counts)[which(colnames(counts) == 346)] = "346.1"
# colnames(counts)[which(colnames(counts) == 350)] = "350.1"
# colnames(counts)[which(colnames(counts) == 352)] = "352.1"
# colnames(counts)[which(colnames(counts) == 371)] = "371.1"
# 
# x=cross_basic$covar
# 
# #get covs by matching colnames of vst with mouse ID in raw pheno file
# covs = x[match(colnames(counts),rownames(x)),]
# 
# 
# #sex, age at sac, generation
# covs = covs[,c(2,3,6)]
# 
# batchCombat = covs[,c(1)]
# batchCombat = model.matrix(~as.factor(sex), data=covs)
# batchCombat = batchCombat[,2]
# #adjust for sex as batch
# #cant adjust for age? maybe must do age*sex
batch = covs$D0_batch
adjusted <- ComBat_seq(as.matrix(counts), batch=batch, group=NULL)
#

#generation is confounded with batch so it is not added
modcombat = model.matrix(~as.factor(sex) + as.factor(age_at_sac_days), data=covs)

batch = covs$D0_batch
#batch removal
edata = ComBat(dat=vst, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)



dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts, 
                              colData = covs, 
                              design = ~sex+age_at_sac_days)

dds.adj <- DESeq2::DESeqDataSetFromMatrix(countData = adjusted, 
                                      colData = covs, 
                                      design = ~sex+age_at_sac_days)

vsd <- varianceStabilizingTransformation(dds, blind = TRUE) #blind so design doesnt matter. seems to be the optimal choice for PCA/viz. 

vsd.adjusted <- varianceStabilizingTransformation(dds.adj, blind = TRUE)

#

##########factominer#########
rv=rowVars(assay(vsd))
select = order(rv,decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat = t(assay(vsd)[select,])
dat2 = cbind(covs,dat)

##########################
rv=rowVars(assay(vsd.adjusted))
select = order(rv,decreasing = TRUE)[seq_len(min(500,length(rv)))]
dat = t(assay(vsd.adjusted)[select,])
dat2 = cbind(covs,dat)

#col1a1 then Xist then col1a2 most variable after "adjusting" for sex 
######################
pca = PCA(dat2[,5:ncol(dat2)],graph=F,scale.unit = F)

fviz_screeplot(pca)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$D0_batch,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="D0_batch"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$sex,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="Sex"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$age_at_sac_days,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(1,2),
             legend.title="age"
             
)

fviz_pca_ind(pca,
             geom.ind = "point",
             col.ind = dat2$sex,
             #palette = c("#00AFBB", "#E7B800"),
             addEllipses = T,
             axes=c(3,4),
             legend.title="Sex"
             
)

