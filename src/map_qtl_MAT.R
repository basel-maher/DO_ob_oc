#QTL mapping 
set.seed(8675309)
library(qtl2)
#load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the allele probs
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
#load(file = "./results/Rdata/k_basic_cleaned.Rdata") #overall (not used)
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file


###

#normalize  phenos##

quant_pheno_columns = c(5:79,81,83,85)#quantitative phenotypes

for(i in quant_pheno_columns){
  print(paste(i,shapiro.test(cross_basic$pheno[,i])$p.value)) #if pval < alpha, not Normal 
  
  if(shapiro.test(cross_basic$pheno[,i])$p.value>0.05){
    print(i)
  }
}
#only 15,18,19,20,22,47,49.50,51,56,59,62,68 are Normal



#qtl mapping for MAT

DO_qtl_scan_MAT = scan1(apr, cross_basic$pheno[,c(75:79,81,83,85)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 2)
#save(DO_qtl_scan_MAT,file = "./results/Rdata/DO_qtl_scan_MAT.Rdata")

#scan MAT as binary traits
#fix bin, in vol3 and vol4 some are 1 in bin but are na in the others
cross_basic$pheno[which(is.na(cross_basic$pheno[,"MAT_VOL1"])),"MAT_VOL1_bin"] = NA
cross_basic$pheno[which(is.na(cross_basic$pheno[,"MAT_VOL2"])),"MAT_VOL2_bin"] = NA
cross_basic$pheno[which(is.na(cross_basic$pheno[,"MAT_VOL3"])),"MAT_VOL3_bin"] = NA
cross_basic$pheno[which(is.na(cross_basic$pheno[,"MAT_VOL4"])),"MAT_VOL4_bin"] = NA

DO_qtl_scan_MAT_binary = scan1(apr, cross_basic$pheno[,c(80,82,84,86)], Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 2,model = "binary")
#save(DO_qtl_scan_MAT_binary,file = "./results/Rdata/DO_qtl_scan_MAT_binary.Rdata")
#load("./results/Rdata/DO_qtl_scan_MAT_binary.Rdata")

#find peaks and bind them together
qtl_peaks = find_peaks(DO_qtl_scan_MAT, cross_basic$pmap, threshold=4, drop=1.5)
qtl_peaks_binary = find_peaks(DO_qtl_scan_MAT_binary, cross_basic$pmap, threshold=4, drop=1.5)

qtl_peaks_both_MAT = rbind(qtl_peaks,qtl_peaks_binary)
write.csv(qtl_peaks_both_MAT, file = "./results/flat/qtl_peaks_MAT.csv",row.names = FALSE,quote = FALSE)


####try after transforming####
norm_pheno = as.data.frame(cross_basic$pheno)

#norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
#norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
#norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
#norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(75:79,81,83,85)]))
pheno_combined = norm_pheno
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.


DO_qtl_scan_MAT_normal = scan1(apr, pheno_combined, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
#save(DO_qtl_scan_normal,file = "./results/Rdata/DO_qtl_scan_norm.Rdata")
load("./results/Rdata/DO_qtl_scan_norm.Rdata")

qtl_peaks_norm = find_peaks(DO_qtl_scan_MAT_normal, cross_basic$pmap, threshold=4, drop=1.5)


#qtl_peaks_bin_norm = find_peaks(DO_qtl_scan_binary_norm, cross_basic$pmap, threshold=4, drop=1.5)
#qtl_peaks_both_norm = rbind(qtl_peaks_norm,qtl_peaks_bin_norm)

write.csv(qtl_peaks_both_norm, file = "./results/flat/qtl_peaks_norm.csv",row.names = FALSE,quote = FALSE)
#
#
#

#calc heritability
h = est_herit(pheno = pheno_combined, kinship = k,addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")] )

h_df = as.data.frame(h)
#



##same for norm
qtl_peaks_norm$perm_thresh = NA
perm = list.files("./results/Rdata/qtl_perms/")

perm_files = perm[grep("norm_perms_MAT",perm)]
perm_files = perm_files[-grep("_bin",perm_files)]
#remove INT
#perm_files = perm_files[-grep("_INT_",perm_files)]

for(i in 1:length(perm_files)){
  
  load(paste0("./results/Rdata/qtl_perms/",perm_files[i]))
  
  perm_a = summary(norm_perm,alpha = 0.05)$A[1]
  perm_x = summary(norm_perm,alpha=0.05)$X[1]
  
  pheno_name = perm_files[i]
  pheno_name = gsub(x = perm_files[i],pattern = "norm_perms_",replacement = "")
  pheno_name = gsub(x = pheno_name,pattern = ".Rdata",replacement = "")
  
  pheno_rows = which(qtl_peaks_norm$lodcolumn == pheno_name)
  
  for(i in 1:length(pheno_rows)){
    if(qtl_peaks_norm$chr[pheno_rows[i]] == "X"){
      qtl_peaks_norm$perm_thresh[[pheno_rows[i]]] = perm_x
    } else {qtl_peaks_norm$perm_thresh[[pheno_rows[i]]] = perm_a}
  }
}




##########
#merge analysis chr18 MAT

#maxmarg
#plotpxg

winsor = cross_basic$pheno[,c(75:79,81,83,85)]

for(i in 1:ncol(winsor)){
  print(i)
  min = quantile(winsor[,i],0.05, na.rm = T)
  print(min)
  max = as.numeric(quantile(winsor[,i],0.95, na.rm = T))
  print(max)
  winsor[which(as.numeric(winsor[,i])>= max),i] = max
  winsor[which(as.numeric(winsor[,i])<= min),i] = min}



winsor_MAT = scan1(apr, winsor[,c(1:8)], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 2)
winsor_peaks = find_peaks(winsor_MAT, cross_basic$pmap, threshold=4, drop=1.5)
write.csv(winsor_peaks, file = "./results/flat/winsor_peaks.csv",row.names = FALSE,quote = FALSE)

norm_winsor = as.data.frame(winsor)


norm_winsor = as.data.frame(log10(norm_winsor))

is.na(norm_winsor) = sapply(norm_winsor, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.


winsor_MAT_norm = scan1(apr, norm_winsor, k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
winsor_peaks_norm = find_peaks(winsor_MAT_norm, cross_basic$pmap, threshold=4, drop=1.5)
write.csv(winsor_peaks_norm, file = "./results/flat/winsor_peaks_norm.csv",row.names = FALSE,quote = FALSE)


MAT_VOL2_winsor_blup = scan1blup(apr[,"18"], winsor[,c("MAT_VOL2")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL2_winsor_blup[3200:3850,],cross_basic$pmap)


MAT_VOL3_winsor_blup = scan1blup(apr[,"18"], winsor[,c("MAT_VOL3")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL3_winsor_blup[3200:3850,],cross_basic$pmap)







start = 84.037
end = 84.245
chr = 18
out_snps_MAT3_18_winsor <- scan1snps(pr, cross_basic$pmap, winsor[,c("MAT_VOL3")], k_loco[["18"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
                              query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

if("pseudogene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
}

if("miRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
}

if("rRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
}

top_MAT3 <- top_snps(out_snps_MAT3_18_winsor$lod, out_snps_MAT3_18_winsor$snpinfo, drop = 0.15 * max(out_snps_MAT3_18$lod))
top_MAT3[order(top_MAT3$lod,decreasing = T),]
#top snp is BL6 private
plot_snpasso(out_snps_MAT3_18_winsor$lod, out_snps_MAT3_18_winsor$snpinfo, genes = genes_locus)



out_snps_MAT2_18_winsor <- scan1snps(pr, cross_basic$pmap, winsor[,c("MAT_VOL2")], k_loco[["18"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
                                     query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


variants_locus = query_variants(chr, start, end)
genes_locus <- query_genes(chr, start, end)

genes_locus = genes_locus[-grep("Gm",genes_locus$Name),]

if("pseudogene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "pseudogene"),]
}

if("miRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "miRNA gene"),]
}

if("rRNA gene" %in% genes_locus$mgi_type){
  genes_locus = genes_locus[-which(genes_locus$mgi_type == "rRNA gene"),]
}

top_MAT2 <- top_snps(out_snps_MAT2_18_winsor$lod, out_snps_MAT2_18_winsor$snpinfo, drop = 0.15 * max(out_snps_MAT3_18$lod))
top_MAT2[order(top_MAT2$lod,decreasing = T),]
#top snp is BL6 private
plot_snpasso(out_snps_MAT2_18_winsor$lod, out_snps_MAT2_18_winsor$snpinfo, genes = genes_locus)
