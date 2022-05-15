annot_ob = read.delim("./results/flat/RNA-seq/all_gene_abund_filt_0.1tpm_ob")
annot_ob = annot_ob[,c(1:6)]
annot_ob = unique(annot_ob)
write.csv(annot_ob,"results/flat/annot_file_ob.csv")


annot_oc = read.delim("./results/flat/RNA-seq/all_gene_abund_filt_0.1tpm_oc")
annot_oc = annot_oc[,c(1:6)]
annot_oc = unique(annot_oc)
write.csv(annot_oc,"./results/flat/annot_file_oc.csv")

annot_bone = read.delim("./results/flat/RNA-seq/all_gene_abund_bone_filt_0.1tpm")
annot_bone = annot_bone[,c(1:6)]
annot_bone = unique(annot_bone)
write.csv(annot_bone,"./results/flat/annot_file_bone.csv")
