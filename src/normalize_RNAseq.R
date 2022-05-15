#This script normalizes RNA-seq counts#
library(DESeq2)

#ob
#read in the RNA-seq counts
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS_ob.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#variance-stabilize the counts
dds.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))


#match ID keys from psomagen to mouse IDs
keys = read.csv("./data/ID Key_ Psomagen BM Ob_Oc RNA seq.csv")
keys$mouse.ID = as.character(keys$mouse.ID)
keys = keys[keys$cell.type == "osteoblast",]

keys$Psomagen.ID = sapply(strsplit(keys$Psomagen.ID, " "),"[",3)

#edit colnames to be only the number
colnames(dds.vst) = sapply(strsplit(colnames(dds.vst), "-"),"[",3)
colnames(dds.vst) = sapply(strsplit(colnames(dds.vst), "_"),"[",1)

colnames(dds.vst) = keys[match(colnames(dds.vst), keys$Psomagen.ID),"mouse.ID"]

#add ".1" to some colnames so they match the GigaMUGA data
# sample 4 needs ".1" added
#colnames(dds.vst)[which(colnames(dds.vst) == 4)] = "4.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 265)] = "265.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 271)] = "271.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 291)] = "291.1"
colnames(dds.vst)[which(colnames(dds.vst) == 324)] = "324.1"
# colnames(dds.vst)[which(colnames(dds.vst) == 346)] = "346.1"
# colnames(dds.vst)[which(colnames(dds.vst) == 350)] = "350.1"
# colnames(dds.vst)[which(colnames(dds.vst) == 352)] = "352.1"
colnames(dds.vst)[which(colnames(dds.vst) == 371)] = "371.1"
# colnames(dds.vst)[which(colnames(dds.vst) == 667)] = "667.1"
# colnames(dds.vst)[which(colnames(dds.vst) == 727)] = "727.1"

#transpose
dds.vst = t(dds.vst)
dds.vst_n = as.data.frame(dds.vst)

#perform quantile-based inverse normal transform. aka match each gene datapoint to a quantile, then match to quantile in normal distribution
#from https://www.nature.com/articles/nature11401#s1 (FTO genotype BMI Visscher et al)

for(col in 1:ncol(dds.vst)){
  dds.vst[,col] = qnorm((rank(dds.vst[,col],na.last="keep")-0.5)/sum(!is.na(dds.vst[,col])))
}

dds.vst = as.data.frame(dds.vst)

dds.vst$Mouse.ID = rownames(dds.vst)
dds.vst = dds.vst[,c(ncol(dds.vst), 1:(ncol(dds.vst)-1))]
#write
write.table(dds.vst, "./results/flat/counts_vst_qnorm_ob.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
save(object = dds.vst, file = "./results/Rdata/counts_vst_qnorm_ob.Rdata")

dds.vst = dds.vst[,-1]
write.table(dds.vst, "./results/flat/counts_vst_qnorm_nohead_ob.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

###################################################################################################################################

#oc
#read in the RNA-seq counts
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_37samps_tpm_over0.1_37samps_COUNTS_oc.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#variance-stabilize the counts
dds.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))

#match ID keys from psomagen to mouse IDs
keys = read.csv("./data/ID Key_ Psomagen BM Ob_Oc RNA seq.csv")
keys$mouse.ID = as.character(keys$mouse.ID)
keys = keys[keys$cell.type == "osteoclast",]

keys$Psomagen.ID = sapply(strsplit(keys$Psomagen.ID, " "),"[",3)

#edit colnames to be only the number
colnames(dds.vst) = sapply(strsplit(colnames(dds.vst), "-"),"[",3)
colnames(dds.vst) = sapply(strsplit(colnames(dds.vst), "_"),"[",1)

colnames(dds.vst) = keys[match(colnames(dds.vst), keys$Psomagen.ID),"mouse.ID"]

#add ".1" to some colnames so they match the GigaMUGA data
# sample 4 needs ".1" added
#colnames(dds.vst)[which(colnames(dds.vst) == 4)] = "4.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 265)] = "265.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 271)] = "271.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 291)] = "291.1"
colnames(dds.vst)[which(colnames(dds.vst) == 324)] = "324.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 346)] = "346.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 350)] = "350.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 352)] = "352.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 371)] = "371.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 667)] = "667.1"
#colnames(dds.vst)[which(colnames(dds.vst) == 727)] = "727.1"

#transpose
dds.vst = t(dds.vst)
dds.vst_n = as.data.frame(dds.vst)

#perform quantile-based inverse normal transform. aka match each gene datapoint to a quantile, then match to quantile in normal distribution
#from https://www.nature.com/articles/nature11401#s1 (FTO genotype BMI Visscher et al)

for(col in 1:ncol(dds.vst)){
  dds.vst[,col] = qnorm((rank(dds.vst[,col],na.last="keep")-0.5)/sum(!is.na(dds.vst[,col])))
}

dds.vst = as.data.frame(dds.vst)

dds.vst$Mouse.ID = rownames(dds.vst)
dds.vst = dds.vst[,c(ncol(dds.vst), 1:(ncol(dds.vst)-1))]
#write
write.table(dds.vst, "./results/flat/counts_vst_qnorm_oc.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
save(object = dds.vst, file = "./results/Rdata/counts_vst_qnorm_oc.Rdata")

dds.vst = dds.vst[,-1]
write.table(dds.vst, "./results/flat/counts_vst_qnorm_nohead_oc.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")



#######################################################################################
#######################################################################################
#######################################################################################
#bone
#ob
#read in the RNA-seq counts
counts = read.csv("./results/flat/RNA-seq/genes_over_6reads_in_morethan_38samps_tpm_over0.1_38samps_COUNTS_bone.csv", stringsAsFactors = FALSE,row.names = 1,check.names = FALSE)

#variance-stabilize the counts
dds.vst = DESeq2::varianceStabilizingTransformation(as.matrix(counts))


#add ".1" to some colnames so they match the GigaMUGA data
# samples 4 and 371 need ".1" added
colnames(dds.vst)[which(colnames(dds.vst) == 4)] = "4.1"
colnames(dds.vst)[which(colnames(dds.vst) == 371)] = "371.1"

#transpose
dds.vst = t(dds.vst)
dds.vst_n = as.data.frame(dds.vst)

#perform quantile-based inverse normal transform. aka match each gene datapoint to a quantile, then match to quantile in normal distribution
#from https://www.nature.com/articles/nature11401#s1 (FTO genotype BMI Visscher et al)

for(col in 1:ncol(dds.vst)){
  dds.vst[,col] = qnorm((rank(dds.vst[,col],na.last="keep")-0.5)/sum(!is.na(dds.vst[,col])))
}

dds.vst = as.data.frame(dds.vst)

dds.vst$Mouse.ID = rownames(dds.vst)
dds.vst = dds.vst[,c(ncol(dds.vst), 1:(ncol(dds.vst)-1))]
#write
write.table(dds.vst, "./results/flat/counts_vst_qnorm_bone.csv", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
save(object = dds.vst, file = "./results/Rdata/counts_vst_qnorm_bone.Rdata")

dds.vst = dds.vst[,-1]
write.table(dds.vst, "./results/flat/counts_vst_qnorm_nohead_bone.csv", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")



