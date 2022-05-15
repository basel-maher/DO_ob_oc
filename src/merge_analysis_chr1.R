
####merge analysis for locus on chromosome 10



set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")
#load the cross file 
load(file = "./results/Rdata/cross_basic_cleaned.Rdata")

#load the cross file 
load(file = "../DO_project/results/Rdata/cross_eqtl_REDO.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_basic)

#load qtl mapping object
load("./results/Rdata/DO_qtl_scan_norm.Rdata")

annot_file = read.csv("../DO_project/results/flat/annot_file.csv", stringsAsFactors = FALSE)



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


#get qtl list, passed threshold
qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)

##
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

new_covar = covar
is.na(new_covar) = sapply(new_covar, is.infinite) #convert is.infinite to NA. Basically getting rid of zero observations.

pheno_combined = as.matrix(pheno_combined)
for(i in ncol(pheno_combined)){
  names(pheno_combined[,i]) = names(cross_basic$pheno[,6])
}


query_variants <- create_variant_query_func("../DO_project/data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("../DO_project/data/CCdb/mouse_genes_mgi.sqlite")



#map the phenotypes
#pheno_combined includes normalized and non-normalized phenos, from map_qtl.R
locus8_scan = scan1(apr, pheno_combined[,c("RFP","bending_PYD","bending_disp_at_frax","bending_disp_at_max_load","bending_total_work","bending_work_post_yield","uCT_Ct.TMD","FFP")], k_loco, Xcovar=Xcovar, addcovar = new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)


plot_scan1(locus8_scan[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", col="red")
#plot_scan1(locus1_scan[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "ML",add=T,col="red")

# #ML 1
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 10, thechr = c(10), thepos =c(23.76194))
# points(xpos, 10.011246, pch=49, bg="black")
# 
# #TMD 2
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.10491))
# points(xpos, 23.9, pch=50, bg="black")
# 
# #Ma.Ar 3
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.32971))
# points(xpos, 12.79, pch=51, bg="black")
# 
# #Tt.Ar 4
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.19293))
# points(xpos, 11.46, pch=52, bg="black")
# 
# #Ct.porosity 5
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.35757))
# points(xpos, 11.385, pch=53, bg="black")
# 
# #pMOI 6
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.10491))
# points(xpos, 8.75, pch=54, bg="black")
# 
# #ct.ar/tt.ar 7
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.28861))
# points(xpos, 8.5, pch=55, bg="black")
# 
# #Imax 8
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.10491))
# points(xpos, 8.27, pch=56, bg="black")


plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", add=T)
plot(locus8_scan, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T)


#we first looked at TMD for. do a snp scan in conf interval
start = 23.18076
end = 23.70523
chr = 10
out_snps_TMD_10 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ct.TMD"], k_loco[["10"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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
#top snps are AJ/NZO/CAST/PWK
plot_snpasso(out_snps_TMD_10$lod, out_snps_TMD_10$snpinfo, genes = genes_locus)

#condition on top snp rs107989322  
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.58619),
                      sdp=113,
                      snp=c("rs107989322"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
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
locus10_scan_cond = scan1(apr, pheno_combined[,c("RFP","bending_PYD","bending_disp_at_frax","bending_disp_at_max_load","bending_total_work","bending_work_post_yield","uCT_Ct.TMD","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35", "alleleA")],cores = 1)


plot_scan1(locus10_scan_cond[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", col="red",)

plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)

thresh = 7.8
abline(h=thresh, col="red")
peaks = find_peaks(locus10_scan_cond, cross_basic$pmap, threshold=4, drop=1.5)
 #all ablated except FFP, maybe RFP

##REDO for snp that is in EYA4 intron
#condition on top snp rs225103945  
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.30437),
                      sdp=112,
                      snp=c("rs225103945"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
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
locus10_scan_cond = scan1(apr, pheno_combined[,c("RFP","bending_PYD","bending_disp_at_frax","bending_disp_at_max_load","bending_total_work","bending_work_post_yield","uCT_Ct.TMD","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35", "alleleA")],cores = 1)


plot_scan1(locus10_scan_cond[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", col="red",)

plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)

thresh = 7.8
abline(h=thresh, col="red")
peaks = find_peaks(locus10_scan_cond, cross_basic$pmap, threshold=4, drop=1.5)

# All go away except TMD, FFP, maybe RFP
#ma.ar drop from 12.8 to 5.9 thresh is 7.6

#TMD doesnt go away but is reduced by about half 




























#now do for total work
start = 23.49830
end = 24.57739
chr = 10
out_snps_totwork_10 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"bending_total_work"], k_loco[["10"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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
top_work[order(out_snps_totwork_10$lod,decreasing = T),]
#same SDP 113 top snp as first RFP one
plot_snpasso(out_snps_totwork_10$lod, out_snps_totwork_10$snpinfo, genes = genes_locus)

#condition on top snp rs107989322  
snpinfo <- data.frame(chr=c("10"),
                      pos=c(23.58619),
                      sdp=113,
                      snp=c("rs107989322"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
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
locus10_scan_cond = scan1(apr, pheno_combined[,c("RFP","bending_PYD","bending_disp_at_frax","bending_disp_at_max_load","bending_total_work","bending_work_post_yield","uCT_Ct.TMD","FFP")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35", "alleleA")],cores = 1)


plot_scan1(locus10_scan_cond[60000:61000,], map = cross_basic$pmap, chr = 10, lodcolumn = "FFP", col="red",)

plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_PYD", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_frax", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_disp_at_max_load", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_total_work", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "bending_work_post_yield", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "RFP", add=T)
plot(locus10_scan_cond, map = cross_basic$pmap, chr = 10, lodcolumn = "uCT_Ct.TMD", add=T)

thresh = 7.8
abline(h=thresh, col="red")
peaks = find_peaks(locus10_scan_cond, cross_basic$pmap, threshold=4, drop=1.5)
#all ablated except FFP, maybe RFP












#do snp scan for TMD in CI
start = 154.80370
end = 155.60037
chr = 1

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

out_snps_TMD_1 <- scan1snps(pr, cross_basic$pmap, pheno_combined[,"uCT_Ct.TMD"], k_loco[["1"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33")],Xcovar=Xcovar,
                            query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)


top_TMD <- top_snps(out_snps_TMD_1$lod, out_snps_TMD_1$snpinfo,drop=0.15*max(out_snps_TMD_1$lod))
top_TMD[order(top_TMD$lod,decreasing = T),]
plot_snpasso(out_snps_TMD_1$lod, out_snps_TMD_1$snpinfo, genes = genes_locus)

#top snps are also WSB private but seem to be different. kinda two different peaks though. condition on topmost just to see

#condition on top snp rs248974780
snpinfo <- data.frame(chr=c("1"),
                      pos=c(155.0559),
                      sdp=128,
                      snp=c("rs248974780"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}



#redo qtl scan while conditioning on snp
locus1_scan_cond_TMD = scan1(apr, pheno_combined[,c("uCT_Ct.TMD","uCT_pMOI","uCT_Imax","uCT_Ct.Ar.Tt.Ar", "uCT_Tt.Ar","ML", "uCT_Ma.Ar","uCT_Ct.porosity")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33", "alleleA")],cores = 1)

plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "ML",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.TMD",ylim=c(0,15))
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_pMOI",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Imax",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.Ar.Tt.Ar",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Tt.Ar",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ma.Ar",add=T)
plot(locus1_scan_cond_TMD, map = cross_basic$pmap, chr = 1, lodcolumn = "uCT_Ct.porosity",add=T)

peaks = find_peaks(locus1_scan_cond_TMD, cross_basic$pmap, threshold=4, drop=1.5)
abline(h=thresh, col="red")

#everything goes away




#compare the two WSB private snps

snpinfo <- data.frame(chr=c("1","1"),
                      pos=c(155.4623,155.0559),
                      sdp=c(128,128),
                      snp=c("rs50769082","rs248974780"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`1`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##

covar_snp$ML_snp = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs50769082[i] >0.6){
    covar_snp$ML_snp[i] = 0
  } else{if(covar_snp$B.rs50769082[i] >0.6){
    covar_snp$ML_snp[i] = 1
  }}
}

covar_snp$TMD_snp = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs248974780[i] >0.6){
    covar_snp$TMD_snp[i] = 0
  } else{if(covar_snp$B.rs248974780[i] >0.6){
    covar_snp$TMD_snp[i] = 1
  }}
}

which(covar_snp$alleleA != covar_snp$alleleB)
#some difference, some mice are hets for one and full for the other
covar_snp[which(covar_snp$alleleA != covar_snp$alleleB),]


# add labels for mice that have RNA-seq data
covar_snp$RNA = 0

covar_snp[which(covar_snp$Row.names %in% rownames(cross_eqtl$pheno)),"RNA"] = 1

#get mice that have different alleles

diff = covar_snp[which(covar_snp$ML_snp != covar_snp$TMD_snp),]



# find founder alleles at each snp
#9082 marker (ML WSB snp)
find_marker(map = cross_basic$pmap, chr = 1, pos = 155.4623)
#UNCHS002851
find_markerpos(cross_basic, "UNCHS002851")
#155.4622


#
#780 marker (TMD WSB snp)
find_marker(map = cross_basic$pmap, chr = 1, pos = 155.0559)
#UNCHS002843
find_markerpos(cross_basic, "UNCHS002843")
#155.05

##
geno_ML_marker = maxmarg(pr, minprob=0.9,chr = 1, pos = 155.4623, map = cross_basic$pmap, return_char = T)

for(i in 1:nrow(covar_snp)){
  x = which(names(geno_ML_marker) == covar_snp$Row.names[i])
  covar_snp$ML_pos_geno[i] = geno_ML_marker[x]
}


geno_TMD_marker = maxmarg(pr, minprob=0.9,chr = 1, pos = 155.0559, map = cross_basic$pmap, return_char = T)

for(i in 1:nrow(covar_snp)){
  x = which(names(geno_TMD_marker) == covar_snp$Row.names[i])
  covar_snp$TMD_pos_geno[i] = geno_TMD_marker[x]
}


covar_snp[which(covar_snp$alleleA != covar_snp$alleleB),]

