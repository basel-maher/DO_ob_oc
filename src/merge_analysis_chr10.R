####merge analysis for locus on chromosome 10



set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_bone.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#load qtl mapping object
load("./results/Rdata/DO_qtl_scan_MAT.Rdata")

annot_file = read.csv("./results/flat/annot_file_bone.csv", stringsAsFactors = FALSE)

#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file

##
covar_eqtl = as.matrix(cross_eqtl$covar)
covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1

covar_eqtl = covar_eqtl[,-1]#remove sac date as covar for now

covar_eqtl = apply(covar_eqtl,2,as.numeric) #make sure all cols are numeric
rownames(covar_eqtl) = rownames(cross_eqtl$covar)#make sure rownames match original cross file




norm_pheno = as.data.frame(cross_basic$pheno)

norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1

norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1

norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:46,48,52:55,57,58,60,61,63:67,69:79,81,83,85)]))

pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,47,49,50,51,56,59,62,68)])
is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])}

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.


#get qtl list, passed threshold
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)







#map the phenotypes
#pheno_combined includes normalized and non-normalized phenos, from map_qtl.R
locus10_scan = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)
save(locus10_scan, file = "../DO_ob_oc/results/Rdata/locus_10_scan.Rdata")

plot_scan1(locus10_scan[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red")
#plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)
#plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)
#plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)


plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T,col="red")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T,col="blue")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T,col="blue")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T,col="red")
plot(locus10_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T,col="green")



#we first looked at TMD. do a snp scan in conf interval
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


start = 23.18076
end = 23.70523
chr = 10
out_snps_TMD_10 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ct.TMD"], k_loco[["10"]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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

top_TMD <- top_snps(out_snps_TMD_10$lod, out_snps_TMD_10$snpinfo, drop = 0.15 * max(out_snps_TMD_10$lod))
top_TMD[order(top_TMD$lod,decreasing = T),]
plot_snpasso(out_snps_TMD_10$lod, out_snps_TMD_10$snpinfo, genes = genes_locus)

#multiple SDPs
#use top SNP rs107989322 as marker
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.58619),
                      sdp=113,
                      snp=c("rs107989322"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)

covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs107989322[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs107989322[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus10_scan_cond = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)

save(locus10_scan_cond, file = "../DO_ob_oc/results/Rdata/locus10_scan_cond.Rdata")


plot_scan1(locus10_scan_cond[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))

plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")

#RFP and FFP seem separate


TMD_blup_cond_113 = scan1blup(apr[,"10"], pheno_combined[,c("uCT_Ct.TMD")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
plot_coefCC(TMD_blup_cond_113[100:1300,],cross_basic$pmap,legend = "bottomright")

save(TMD_blup_cond_113, file = "../DO_ob_oc/results/Rdata/TMD_blup_cond_113.Rdata")


#Do the same with the highest bending trait
start = 23.49830
end = 24.57739
chr = 10
out_snps_totwork_10 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"bending_total_work"], k_loco[["10"]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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

top_work <- top_snps(out_snps_totwork_10$lod, out_snps_totwork_10$snpinfo, drop = 0.15 * max(out_snps_totwork_10$lod))
top_work[order(top_work$lod,decreasing = T),]
plot_snpasso(out_snps_totwork_10$lod, out_snps_totwork_10$snpinfo, genes = genes_locus)

#one SDP, 113

#so maybe try same as above but different SDP, 134 or 112
#use top SNP rs29333335 as marker
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.43659),
                      sdp=134,
                      snp=c("rs29333335"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)

covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs29333335[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs29333335[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus10_scan_cond_134 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)

save(locus10_scan_cond_134, file = "../DO_ob_oc/results/Rdata/locus10_scan_cond_134.Rdata")

plot_scan1(locus10_scan_cond_134[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))

plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")



TMD_blup_cond_134 = scan1blup(apr[,"10"], pheno_combined[,c("uCT_Ct.TMD")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
plot_coefCC(TMD_blup_cond_134[100:1300,],cross_basic$pmap)

save(TMD_blup_cond_134, file = "../DO_ob_oc/results/Rdata/TMD_blup_cond_134.Rdata")



#

#Same with SDP 112
#rs225103945

snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.30437),
                      sdp=112,
                      snp=c("rs225103945"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)

covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs225103945[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs225103945[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus10_scan_cond_112 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
save(locus10_scan_cond_112, file = "../DO_ob_oc/results/Rdata/locus10_scan_cond_112.Rdata")

plot_scan1(locus10_scan_cond_112[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))

plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")





TMD_blup = scan1blup(apr[,"10"], pheno_combined[,c("uCT_Ct.TMD")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(TMD_blup[100:1300,],cross_basic$pmap,legend = "bottomright")

save(TMD_blup, file = "../DO_ob_oc/results/Rdata/TMD_blup.Rdata")

totwork_blup = scan1blup(apr[,"10"], pheno_combined[,c("bending_total_work")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(totwork_blup[100:1300,],cross_basic$pmap,legend = "bottomright")



TMD_blup_cond_112 = scan1blup(apr[,"10"], pheno_combined[,c("uCT_Ct.TMD")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
plot_coefCC(TMD_blup_cond_112[100:1300,],cross_basic$pmap,legend = "bottomright")

save(TMD_blup_cond_112, file = "../DO_ob_oc/results/Rdata/TMD_blup_cond_112.Rdata")



RFP_blup = scan1blup(apr[,"10"], pheno_combined[,c("RFP")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(RFP_blup[100:1300,],cross_basic$pmap,legend = "bottomright")




FFP_blup = scan1blup(apr[,"10"], pheno_combined[,c("FFP")], k_loco[["10"]], addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(FFP_blup[100:1300,],cross_basic$pmap,legend = "bottomright")
















#EYA4

load("./results/Rdata/cross_eqtl_bone.Rdata")
#get the X chrom covars from the cross file
Xcovar <- get_x_covar(cross_eqtl)
load("./results/Rdata/local_eqtl_bone.Rdata")

#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar[,1] = as.factor(covar[,1]) #sac date to factors
covar[,6] = as.factor(covar[,6]) #generation to factors

covar = apply(covar,2,as.numeric)
rownames(covar) = rownames(cross_eqtl$covar)

eya4_blup_bone = scan1blup(apr[,"10"], cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco[["10"]], addcovar = covar[,c(2,11:58)],cores = 4)
plot_coefCC(eya4_blup_bone[100:1300,],cross_eqtl$pmap)


eya4_bone = scan1(apr, cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,11:58)],cores = 4)
plot_scan1(eya4_bone[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = 1, col="red",ylim=c(0,10))

#6.4 LOD. not an eqtl


load("./results/Rdata/cross_eqtl_ob.Rdata")
#get the X chrom covars from the cross file
Xcovar <- get_x_covar(cross_eqtl)


#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar[,1] = as.factor(covar[,1]) #sac date to factors
covar[,6] = as.factor(covar[,6]) #generation to factors

covar = apply(covar,2,as.numeric)
rownames(covar) = rownames(cross_eqtl$covar)

eya4_blup_ob = scan1blup(apr[,"10"], cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco[["10"]], addcovar = covar[,c(2,12:59)],cores = 4)
plot_coefCC(eya4_blup_ob[100:1300,],cross_eqtl$pmap)

eya4_ob = scan1(apr, cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,12:59)],cores = 4)

load("./results/Rdata/local_eqtl_ob.Rdata")
# LOD 10.2, just under threshold

#ob snp scan 
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


start = 23.112649
end = 24.62903
chr = 10
out_snps_ob_eya4 <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,"ENSMUSG00000010461"], k_loco[["10"]],  addcovar = covar[,c(2,12:57)],Xcovar=Xcovar,
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

top_eya4_ob <- top_snps(out_snps_ob_eya4$lod, out_snps_ob_eya4$snpinfo, drop = 0.15 * max(out_snps_ob_eya4$lod))
top_eya4_ob[order(top_eya4_ob$lod,decreasing = T),]
plot_snpasso(out_snps_ob_eya4$lod, out_snps_ob_eya4$snpinfo, genes = genes_locus)





load("./results/Rdata/cross_eqtl_oc.Rdata")
#get the X chrom covars from the cross file
Xcovar <- get_x_covar(cross_eqtl)


#create a covar object from covariates in cross file
#must be numeric
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1 #convert sex to 1's and 0's
covar[,1] = as.factor(covar[,1]) #sac date to factors
covar[,6] = as.factor(covar[,6]) #generation to factors

covar = apply(covar,2,as.numeric)
rownames(covar) = rownames(cross_eqtl$covar)

eya4_blup_oc = scan1blup(apr[,"10"], cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco[["10"]], addcovar = covar[,c(2,12:57)],cores = 4)
plot_coefCC(eya4_blup_oc[675:750,],cross_eqtl$pmap)

eya4_oc = scan1(apr, cross_eqtl$pheno[,c("ENSMUSG00000010461")], k_loco, Xcovar=Xcovar, addcovar = covar[,c(2,12:57)],cores = 4)

load("./results/Rdata/local_eqtl_oc.Rdata")
#lod 54


#oc snp scan 
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


start = 23.11265
end = 23.38773
chr = 10
out_snps_oc_eya4 <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,"ENSMUSG00000010461"], k_loco[["10"]],  addcovar = covar[,c(2,12:57)],Xcovar=Xcovar,
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

top_eya4_oc <- top_snps(out_snps_oc_eya4$lod, out_snps_oc_eya4$snpinfo, drop = 0.15 * max(out_snps_oc_eya4$lod))
top_eya4_oc[order(top_eya4_oc$lod,decreasing = T),]
plot_snpasso(out_snps_oc_eya4$lod, out_snps_oc_eya4$snpinfo, genes = genes_locus)

#SDP 64, PWK

# ########################################################### MERGE ANALYSIS
# #####################
# #Take all loci, get CI of locus +/- 250kb, get all eQTL within that locus, then compare top snps for phenos with top eqtl snps, define list of putatively causal genes
# 
# #qtl file
# qtl_10 = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)
# 
# qtl_10 = qtl_10[which(qtl_10$chr == 10),]
# 
# #define loci
# qtl_10$locus = 0
# 
# loc_idx = 1
# qtl_loc = as.data.frame(matrix(nrow=nrow(qtl_10), ncol=ncol(qtl_10)))
# colnames(qtl_loc) = colnames(qtl_10)
# for (i in c(10)){
#   sub = subset(qtl_10, qtl_10$chr == i)
#   
#   while(any(sub$locus == 0)){
#     min_sub = min(sub$pos)
#     idx = which((sub$pos >= min_sub) & (sub$pos <= min_sub + 1.5))
#     sub[idx,"locus"] = loc_idx
#     qtl_loc = rbind(qtl_loc, sub[idx,])
#     sub = sub[-idx,]
#     loc_idx = loc_idx + 1
#   }
# }
# 
# qtl_loc = qtl_loc[-which(is.na(qtl_loc)),]
# 
# write.csv(qtl_loc, file = "./results/flat/qtl_loc_10", quote = FALSE,row.names = FALSE)
# 
# qtl_loc = read.csv("./results/flat/qtl_loc_10")
# qtl_loc = qtl_loc[which(qtl_loc$locus==1),]
# 
# 
# #load the cross file 
# load(file = "./results/Rdata/cross_basic_cleaned.Rdata")
# 
# Xcovar <- get_x_covar(cross_basic)
# 
# annot_file = read.csv("./results/flat/annot_file_bone.csv", stringsAsFactors = FALSE)
# 
# #create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
# covar = as.matrix(cross_basic$covar)
# covar[,"sex"] = (covar[,"sex"] == "M")*1
# 
# covar = covar[,-1]#remove sac date as covar for now
# 
# covar = apply(covar,2,as.numeric) #make sure all cols are numeric
# rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file
# 
# ##
# covar_eqtl = as.matrix(cross_eqtl$covar)
# covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1
# 
# covar_eqtl = covar_eqtl[,-1]#remove sac date as covar for now
# 
# covar_eqtl = apply(covar_eqtl,2,as.numeric) #make sure all cols are numeric
# rownames(covar_eqtl) = rownames(cross_eqtl$covar)#make sure rownames match original cross file
# 
# 
# 
# 
# norm_pheno = as.data.frame(cross_basic$pheno)
# 
# norm_pheno$MAT_VOL1 = norm_pheno$MAT_VOL1 + 1
# norm_pheno$MAT_VOL2 = norm_pheno$MAT_VOL2 + 1
# norm_pheno$MAT_VOL3 = norm_pheno$MAT_VOL3 + 1
# norm_pheno$MAT_VOL4 = norm_pheno$MAT_VOL4 + 1
# 
# norm_pheno$bending_work_post_yield = norm_pheno$bending_work_post_yield + 1
# norm_pheno$bending_PYD = norm_pheno$bending_PYD + 1
# 
# norm_pheno = as.data.frame(log10(norm_pheno[,c(6:14,16,17,21,23:46,48,52:55,57,58,60,61,63:67,69:79,81,83,85)]))
# 
# pheno_combined = cbind(norm_pheno, cross_basic$pheno[,c(15,18,19,20,22,47,49,50,51,56,59,62,68)])
# is.na(pheno_combined) = sapply(pheno_combined, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.
# 
# pheno_combined = as.matrix(pheno_combined)
# for(i in ncol(pheno_combined)){
#   names(pheno_combined[,i]) = names(cross_basic$pheno[,6])}
# 
# new_covar = covar
# is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.
# 
# 
# 
# 
# 
# merge = list()
# query_variants <- create_variant_query_func("../DO_project//data/CCdb/cc_variants.sqlite")
# query_genes <- create_gene_query_func("../DO_project/data/CCdb/mouse_genes_mgi.sqlite")
# 
# 
# 
# for(i in 1:nrow(qtl_loc)){
#   print(i)
#   merge[[i]] = list()
#   pheno = qtl_loc$lodcolumn[i]
#   chr = qtl_loc$chr[i]
#   start = qtl_loc$ci_lo[i]
#   end = qtl_loc$ci_hi[i]
#   #use same covars as qtl mapping.sex,age,BW and gen
#   out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,pheno], k_loco[[chr]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
#                         query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
#   merge[[i]] = out_snps
#   names(merge)[i] = paste0(pheno,"_",chr)
#   rm(out_snps)
#   
# }
# 
# save(merge,file = "./results/Rdata/merge_QTL_10.Rdata")
# 
# 
# 
# 
# #for each merge analysis, take the snps that are within 15% LODs of the max LOD
# merge_top = list()
# for(i in 1:length(merge)){
#   print(i)
#   merge_top[[i]] = list()
#   merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
#   
# }
# names(merge_top) = names(merge)
# 
# 
# save(merge_top,file = "./results/Rdata/merge_top_QTL_10.Rdata")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #for each locus, take genes that have eqtl , and are within locus CI +/- 250 kb
# #any portion of gene start or end is within CI
# 
# #load eqtl
# load("./results/Rdata/local_eqtl_ob.Rdata")
# 
# 
# #add eya4 for ob, LOD is just under threshold (10.2) but for bone its much smaller, 6.4
# 
# x = c("ENSMUSG00000010461", "10","23.94633","10.235654","23.112649","24.62903","Eya4","10","-","23102963","23350786","89523")
# 
# local_eqtl_ob = rbind(local_eqtl_ob, x)
# 
# eqtl_loc = as.data.frame(matrix(nrow=nrow(local_eqtl_ob), ncol=ncol(local_eqtl_ob)+1))
# colnames(eqtl_loc) = colnames(local_eqtl_ob)
# colnames(eqtl_loc)[13] = "locus"
# 
# for(i in unique(qtl_loc$locus)){
#   print(i)
#   sub = subset(qtl_loc, qtl_loc$locus == i)
#   min_loc = min(sub$ci_lo) - 0.25
#   max_loc = max(sub$ci_hi) + 0.25
#   
#   sub_eqtl_1 = subset(local_eqtl_ob, local_eqtl_ob$chr == unique(sub$chr) & (as.numeric(local_eqtl_ob$Start)/1000000 >= min_loc) & as.numeric(local_eqtl_ob$Start)/1000000 <= max_loc)
#   sub_eqtl_2 = subset(local_eqtl_ob, local_eqtl_ob$chr == unique(sub$chr) & (as.numeric(local_eqtl_ob$End)/1000000 >= min_loc) & as.numeric(local_eqtl_ob$End)/1000000 <= max_loc)
#   
#   sub_eqtl = merge(sub_eqtl_1, sub_eqtl_2)
#   sub_eqtl = unique(sub_eqtl)
#   sub_eqtl$locus = i
#   
#   eqtl_loc = rbind(eqtl_loc, sub_eqtl)
#   
# }
# eqtl_loc = eqtl_loc[-which(is.na(eqtl_loc)),]
# 
# 
# 
# ####
# ####
# #load merge analysis objects
# 
# load("./results/Rdata/merge_top_local_eqtl_ob_EYA4_chr10_added.Rdata")
# merge_top_eqtl = merge_top
# 
# load("./results/Rdata/merge_top_QTL.Rdata")
# merge_top_qtl = merge_top
# 
# rm(merge_top)
# 
# 
# #for each locus, take eqtl merge analyses in locus and colocalize with pheno merge analyses in locus
# 
# 
# genes = c()
# phenos_w_genes = c()
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   sub_eqtl = subset(eqtl_loc, eqtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   eqtl_genes = paste0(sub_eqtl$lodcolumn, "_", unique(sub_eqtl$chr))
#   
#   
#   for(j in phenos){
#     for(k in eqtl_genes){
#       print(k)
#       if(any(merge_top_eqtl[[k]]$snp_id %in% merge_top_qtl[[j]]$snp_id)){
#         print(k)
#         genes = append(genes,k)
#         phenos_w_genes = append(phenos_w_genes, j)
#       }
#     }
#   }
# }
# 
# 
# unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unique(gene_names)
# 
# eqtl_loc[which(eqtl_loc$lodcolumn %in% genes),]
# 
# 
# phenos_w_genes = as.data.frame(phenos_w_genes)
# phenos_w_genes$gene = genes
# phenos_w_genes$gene_name = NA
# 
# for(i in 1:nrow(phenos_w_genes)){
#   phenos_w_genes$gene[i] = unlist(strsplit(phenos_w_genes$gene[i], "_"))[1]
#   phenos_w_genes$gene_name[i] = eqtl_loc[which(eqtl_loc$lodcolumn == phenos_w_genes$gene[i]),"Gene.Name"]
# }
# 
# #eqtl genes that are located within a phenotypic qtl and regulated by a local eqtl
# phenos_w_genes
# 
# 
# #NONE
# 
# ## number nonsynonymous variants in top qtl merge for a phenotype
# nonsyn = list()
# counter=1
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   
#   for(j in phenos){
#     counter=counter+1
#     missense = (grep("missense",x=merge_top_qtl[[j]]$consequence))
#     l = length(missense)
#     z = merge_top_qtl[[j]]$snp_id[missense]
#     print(paste0(j,":",l))
#     nonsyn[[counter]] = z
#     names(nonsyn)[counter] = j
#   }
# }
# 
# nonsyn_frame = do.call(rbind, lapply(nonsyn, as.data.frame))
# 
# 
# 
# 
# #SIFT: RSIDs uploaded to https://useast.ensembl.org/Tools/VEP
# 
# 
# ################################################################################
# #for each locus, take genes that have eqtl , and are within locus CI +/- 250 kb
# #any portion of gene start or end is within CI
# 
# #load eqtl
# load("./results/Rdata/local_eqtl_oc.Rdata")
# 
# 
# 
# 
# 
# eqtl_loc = as.data.frame(matrix(nrow=nrow(local_eqtl_oc), ncol=ncol(local_eqtl_oc)+1))
# colnames(eqtl_loc) = colnames(local_eqtl_oc)
# colnames(eqtl_loc)[13] = "locus"
# 
# for(i in unique(qtl_loc$locus)){
#   print(i)
#   sub = subset(qtl_loc, qtl_loc$locus == i)
#   min_loc = min(sub$ci_lo) - 0.25
#   max_loc = max(sub$ci_hi) + 0.25
#   
#   sub_eqtl_1 = subset(local_eqtl_oc, local_eqtl_oc$chr == unique(sub$chr) & (as.numeric(local_eqtl_oc$Start)/1000000 >= min_loc) & as.numeric(local_eqtl_oc$Start)/1000000 <= max_loc)
#   sub_eqtl_2 = subset(local_eqtl_oc, local_eqtl_oc$chr == unique(sub$chr) & (as.numeric(local_eqtl_oc$End)/1000000 >= min_loc) & as.numeric(local_eqtl_oc$End)/1000000 <= max_loc)
#   
#   sub_eqtl = merge(sub_eqtl_1, sub_eqtl_2)
#   sub_eqtl = unique(sub_eqtl)
#   sub_eqtl$locus = i
#   
#   eqtl_loc = rbind(eqtl_loc, sub_eqtl)
#   
# }
# eqtl_loc = eqtl_loc[-which(is.na(eqtl_loc)),]
# 
# 
# 
# ####
# ####
# #load merge analysis objects
# 
# load("./results/Rdata/merge_top_local_eqtl_oc.Rdata")
# merge_top_eqtl = merge_top
# 
# load("./results/Rdata/merge_top_QTL.Rdata")
# merge_top_qtl = merge_top
# 
# rm(merge_top)
# 
# 
# #for each locus, take eqtl merge analyses in locus and colocalize with pheno merge analyses in locus
# 
# 
# genes = c()
# phenos_w_genes = c()
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   sub_eqtl = subset(eqtl_loc, eqtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   eqtl_genes = paste0(sub_eqtl$lodcolumn, "_", unique(sub_eqtl$chr))
#   
#   
#   for(j in phenos){
#     for(k in eqtl_genes){
#       print(k)
#       if(any(merge_top_eqtl[[k]]$snp_id %in% merge_top_qtl[[j]]$snp_id)){
#         print(k)
#         genes = append(genes,k)
#         phenos_w_genes = append(phenos_w_genes, j)
#       }
#     }
#   }
# }
# 
# 
# unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unique(gene_names)
# 
# eqtl_loc[which(eqtl_loc$lodcolumn %in% genes),]
# 
# 
# phenos_w_genes = as.data.frame(phenos_w_genes)
# phenos_w_genes$gene = genes
# phenos_w_genes$gene_name = NA
# 
# for(i in 1:nrow(phenos_w_genes)){
#   phenos_w_genes$gene[i] = unlist(strsplit(phenos_w_genes$gene[i], "_"))[1]
#   phenos_w_genes$gene_name[i] = eqtl_loc[which(eqtl_loc$lodcolumn == phenos_w_genes$gene[i]),"Gene.Name"]
# }
# 
# #eqtl genes that are located within a phenotypic qtl and regulated by a local eqtl
# phenos_w_genes
# 
# 
# # Enpp1 RFP
# 
# ## number nonsynonymous variants in top qtl merge for a phenotype
# nonsyn = list()
# counter=1
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   
#   for(j in phenos){
#     counter=counter+1
#     missense = (grep("missense",x=merge_top_qtl[[j]]$consequence))
#     l = length(missense)
#     z = merge_top_qtl[[j]]$snp_id[missense]
#     print(paste0(j,":",l))
#     nonsyn[[counter]] = z
#     names(nonsyn)[counter] = j
#   }
# }
# 
# nonsyn_frame = do.call(rbind, lapply(nonsyn, as.data.frame))
# 
# 
# 
# 
# #SIFT: RSIDs uploaded to https://useast.ensembl.org/Tools/VEP
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #load eqtl
# load("./results/Rdata/local_eqtl_bone.Rdata")
# 
# 
# 
# 
# 
# eqtl_loc = as.data.frame(matrix(nrow=nrow(local_eqtl_bone), ncol=ncol(local_eqtl_bone)+1))
# colnames(eqtl_loc) = colnames(local_eqtl_oc)
# colnames(eqtl_loc)[13] = "locus"
# 
# for(i in unique(qtl_loc$locus)){
#   print(i)
#   sub = subset(qtl_loc, qtl_loc$locus == i)
#   min_loc = min(sub$ci_lo) - 0.25
#   max_loc = max(sub$ci_hi) + 0.25
#   
#   sub_eqtl_1 = subset(local_eqtl_bone, local_eqtl_bone$chr == unique(sub$chr) & (as.numeric(local_eqtl_bone$Start)/1000000 >= min_loc) & as.numeric(local_eqtl_bone$Start)/1000000 <= max_loc)
#   sub_eqtl_2 = subset(local_eqtl_bone, local_eqtl_bone$chr == unique(sub$chr) & (as.numeric(local_eqtl_bone$End)/1000000 >= min_loc) & as.numeric(local_eqtl_bone$End)/1000000 <= max_loc)
#   
#   sub_eqtl = merge(sub_eqtl_1, sub_eqtl_2)
#   sub_eqtl = unique(sub_eqtl)
#   sub_eqtl$locus = i
#   
#   eqtl_loc = rbind(eqtl_loc, sub_eqtl)
#   
# }
# eqtl_loc = eqtl_loc[-which(is.na(eqtl_loc)),]
# 
# 
# 
# ####
# ####
# #load merge analysis objects
# 
# load("./results/Rdata/merge_top_local_eqtl_bone.Rdata")
# merge_top_eqtl = merge_top
# 
# load("./results/Rdata/merge_top_QTL.Rdata")
# merge_top_qtl = merge_top
# 
# rm(merge_top)
# 
# 
# #for each locus, take eqtl merge analyses in locus and colocalize with pheno merge analyses in locus
# 
# 
# genes = c()
# phenos_w_genes = c()
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   sub_eqtl = subset(eqtl_loc, eqtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   eqtl_genes = paste0(sub_eqtl$lodcolumn, "_", unique(sub_eqtl$chr))
#   
#   
#   for(j in phenos){
#     for(k in eqtl_genes){
#       print(k)
#       if(any(merge_top_eqtl[[k]]$snp_id %in% merge_top_qtl[[j]]$snp_id)){
#         print(k)
#         genes = append(genes,k)
#         phenos_w_genes = append(phenos_w_genes, j)
#       }
#     }
#   }
# }
# 
# 
# unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
# gene_names = unique(gene_names)
# 
# eqtl_loc[which(eqtl_loc$lodcolumn %in% genes),]
# 
# 
# phenos_w_genes = as.data.frame(phenos_w_genes)
# phenos_w_genes$gene = genes
# phenos_w_genes$gene_name = NA
# 
# for(i in 1:nrow(phenos_w_genes)){
#   phenos_w_genes$gene[i] = unlist(strsplit(phenos_w_genes$gene[i], "_"))[1]
#   phenos_w_genes$gene_name[i] = eqtl_loc[which(eqtl_loc$lodcolumn == phenos_w_genes$gene[i]),"Gene.Name"]
# }
# 
# #eqtl genes that are located within a phenotypic qtl and regulated by a local eqtl
# phenos_w_genes
# 
# 
# # none
# 
# 
# ## number nonsynonymous variants in top qtl merge for a phenotype
# nonsyn = list()
# counter=1
# for(i in unique(qtl_loc$locus)){
#   #print(i)
#   sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
#   
#   phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
#   
#   
#   for(j in phenos){
#     counter=counter+1
#     missense = (grep("missense",x=merge_top_qtl[[j]]$consequence))
#     l = length(missense)
#     z = merge_top_qtl[[j]]$snp_id[missense]
#     print(paste0(j,":",l))
#     nonsyn[[counter]] = z
#     names(nonsyn)[counter] = j
#   }
# }
# 
# nonsyn_frame = do.call(rbind, lapply(nonsyn, as.data.frame))
# 



#######################################################################

plot_scan1(locus10_scan[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red", ylim=c(0,75))
plot_scan1(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)
plot_scan1(locus10_scan_cond_134, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)
plot_scan1(locus10_scan_cond_112, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)

plot_scan1(eya4_bone, map = cross_basic$pmap, chr = 10, col="red", add=T)
plot_scan1(eya4_ob, map = cross_basic$pmap, chr = 10, col="blue", add=T)
plot_scan1(eya4_oc, map = cross_basic$pmap, chr = 10, col="green", add=T)

#no bone eqtl, threshold ob, high oc
plot_coefCC(TMD_blup[100:1300,],cross_basic$pmap)
plot_coefCC(TMD_blup[100:1300,],cross_basic$pmap)
plot_coefCC(eya4_blup_ob[100:1300,],cross_eqtl$pmap)
plot_coefCC(eya4_blup_oc[100:1300,],cross_eqtl$pmap)




#cond oc lead
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.25901),
                      sdp=64,
                      snp=c("rs255382214"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)


#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file
new_covar=covar
covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs255382214[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs255382214[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus10_scan_cond_oc = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
plot_scan1(locus10_scan_cond_oc[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red", ylim=c(0,10))


plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond_oc, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")


#cond ob lead, 12 first
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.20023),
                      sdp=12,
                      snp=c("rs29382669"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)


#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file
new_covar=covar
covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs29382669[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs29382669[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

locus10_scan_cond_ob_sdp12 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)
plot_scan1(locus10_scan_cond_ob_sdp12[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red", ylim=c(0,10))


plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond_ob_sdp12, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")


#cond ob lead, 98 (pwk)
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.11381),
                      sdp=98,
                      snp=c("rs29330121"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`10`)


#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_basic$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_basic$covar)#make sure rownames match original cross file
new_covar=covar
covar_snp = merge(new_covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs29330121[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs29330121[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

locus10_scan_cond_ob_sdp98 = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","bending_total_work","bending_work_post_yield","bending_disp_at_frax","bending_disp_at_max_load", "bending_PYD", "RFP","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 4)

plot_scan1(locus10_scan_cond_ob_sdp98[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red", ylim=c(0,10))


plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T, col="red")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T, col="red")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T, col="blue")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T, col="yellow")
plot(locus10_scan_cond_ob_sdp98, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T, col="green")







plot_scan1(locus10_scan_cond_oc[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))

plot_scan1(locus10_scan_cond_ob_sdp12[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))

plot_scan1(locus10_scan_cond_ob_sdp98[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red",ylim=c(0,10))







#RFP, FFP PWK?? why does it diminish after conditioning on oc and ob
#do same analysis for other genes in the region





#plot pxg


# find founder alleles at each snp
#3045 marker
find_marker(map = cross_basic$pmap, chr = 10, pos = 23.49253)
#UNC29749856
find_markerpos(cross_basic, "UNCHS027695")
#84.2307




##
geno_TMD_marker = maxmarg(pr, minprob=0.5,chr = 10, pos = 23.49253, map = cross_basic$pmap, return_char = T)

plot_pxg(geno_TMD_marker,pheno = pheno_combined[,"uCT_Ct.TMD"],sort = TRUE)

#
# find founder alleles at each snp #eya4, ENSMUSG00000010461
#3045 marker
find_marker(map = cross_basic$pmap, chr = 10, pos = 23.26126)
#UNC29749856
find_markerpos(cross_basic, "UNCHS027694")
#84.2307

geno_eya4_oc_marker = maxmarg(apr, minprob=0.5,chr = 10, pos = 23.26126, map = cross_eqtl$pmap, return_char = T)

plot_pxg(geno_eya4_oc_marker,pheno = cross_eqtl$pheno[,"ENSMUSG00000010461"],sort = TRUE)












##SNP scan for FFP



query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


start = 24.39475
end = 25.65636
chr = 10
out_snps_RFP_10 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"RFP"], k_loco[["10"]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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

top_RFP <- top_snps(out_snps_RFP_10$lod, out_snps_RFP_10$snpinfo, drop = 0.15 * max(out_snps_RFP_10$lod))
top_RFP[order(top_RFP$lod,decreasing = T),]
plot_snpasso(out_snps_RFP_10$lod, out_snps_RFP_10$snpinfo, genes = genes_locus)
