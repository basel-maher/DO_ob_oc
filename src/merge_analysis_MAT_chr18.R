####merge analysis for locus on chromosome 18
###seems like there are two separate associations, TMD, ct.th and ma.ar and ML


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


#get qtl list, passed threshold
#qtl_norm = read.csv("./results/flat/qtl_norm_pass_thresh", stringsAsFactors = FALSE)





query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")



#map the phenotypes
#pheno_combined includes normalized and non-normalized phenos, from map_qtl.R
locus18_scan = scan1(apr, cross_basic$pheno[,c("MAT_VOL1","MAT_VOL2","MAT_VOL3","MAT_VOL4", "MAT_VOL1_nonzero","MAT_VOL2_nonzero", "MAT_VOL3_nonzero","MAT_VOL4_nonzero")], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)


plot_scan1(locus18_scan[102700:103000,], map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3", col="red")
# #plot_scan1(locus1_scan[5700:6200,], map = cross_basic$pmap, chr = 1, lodcolumn = "ML",add=T,col="red")
# 
# #ML 1
# xpos = xpos_scan1(map = cross_basic$pmap, chr = 1, thechr = c(1), thepos =c(155.35757))
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

plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3_nonzero", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1_nonzero", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2_nonzero", add=T)
plot(locus18_scan, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4_nonzero", add=T)

#we first looked at MAT_VOL3 for qsox1. do a snp scan in conf interval
start = 84.008001
end = 84.238158
chr = 18
out_snps_MAT3_18 <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,"MAT_VOL3"], k_loco[["18"]],  addcovar =  covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
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

top_MAT3 <- top_snps(out_snps_MAT3_18$lod, out_snps_MAT3_18$snpinfo, drop = 0.15 * max(out_snps_MAT3_18$lod))
top_MAT3[order(top_MAT3$lod,decreasing = T),]
#top snp is BL6 private
plot_snpasso(out_snps_MAT3_18$lod, out_snps_MAT3_18$snpinfo, genes = genes_locus)










snpinfo <- data.frame(chr=c("18"),
                      pos=c(84.23055),
                      sdp=2,
                      snp=c("rs29558786"), stringsAsFactors=FALSE)


snpinfo <- index_snps(cross_basic$pmap, snpinfo)
snp_genoprobs = genoprob_to_snpprob(apr,snpinfo)
snp_genoprobs =as.data.frame(snp_genoprobs$`18`)

covar_snp = merge(covar, snp_genoprobs, by="row.names", all = TRUE)
rownames(covar_snp) = covar_snp$Row.names

##
covar_snp$alleleA = 0.5
for(i in 1:nrow(covar_snp)){
  if(covar_snp$A.rs29558786[i] >0.6){
    covar_snp$alleleA[i] = 1
  } else{if(covar_snp$B.rs29558786[i] >0.6){
    covar_snp$alleleA[i] = 0
  }}
}

#redo qtl scan while conditioning on snp
locus18_scan_cond = scan1(apr, cross_basic$pheno[,c("MAT_VOL1","MAT_VOL2","MAT_VOL3","MAT_VOL4", "MAT_VOL1_nonzero","MAT_VOL2_nonzero", "MAT_VOL3_nonzero","MAT_VOL4_nonzero")], k_loco, Xcovar=Xcovar, addcovar = covar_snp[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35","alleleA")],cores = 1)

plot_scan1(locus18_scan_cond[102700:103000,], map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3", col="red")

plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3_nonzero", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1_nonzero", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2_nonzero", add=T)
plot(locus18_scan_cond, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4_nonzero", add=T)


#thresh = 7.8
#abline(h=thresh, col="red")
peaks = find_peaks(locus18_scan_cond, cross_basic$pmap, threshold=4, drop=1.5)

#ML, pMOI,Imax, ct.ar/tt.ar*, tt.ar, ma.ar  and ct.porosity goes away 
#ma.ar drop from 12.8 to 5.9 thresh is 7.6

#TMD doesnt go away but is reduced by about half 

# find founder alleles at each snp
#3045 marker
find_marker(map = cross_basic$pmap, chr = 18, pos = 84.23055)
#UNC29749856
find_markerpos(cross_basic, "UNC29749856")
#84.2307




##
geno_MAT3_marker = maxmarg(pr, minprob=0.5,chr = 18, pos = 84.23055, map = cross_basic$pmap, return_char = T)

plot_pxg(geno_MAT3_marker,pheno = log10(cross_basic$pheno[,"MAT_VOL3"]+1),sort = TRUE)

geno_MAT3_marker = maxmarg(pr, minprob=0.9,chr = 18, pos = 84.2307, map = cross_basic$pmap, return_char = T)
plot_pxg(geno_MAT3_marker,pheno =cross_basic$pheno[,"MAT_VOL3"],sort = TRUE)



MAT_VOL3_scan_blup = scan1blup(apr[,"18"], cross_basic$pheno[,c("MAT_VOL3")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL3_scan_blup[3200:3850,],cross_basic$pmap)


MAT_VOL1_scan_blup = scan1blup(apr[,"18"], cross_basic$pheno[,c("MAT_VOL1")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL1_scan_blup[3200:3850,],cross_basic$pmap)

MAT_VOL1_nonzero_scan_blup = scan1blup(apr[,"18"], cross_basic$pheno[,c("MAT_VOL1_nonzero")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL1_nonzero_scan_blup[3200:3850,],cross_basic$pmap)

MAT_VOL2_scan_blup = scan1blup(apr[,"18"], cross_basic$pheno[,c("MAT_VOL2")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL2_scan_blup[3200:3850,],cross_basic$pmap)

MAT_VOL2_nonzero_scan_blup = scan1blup(apr[,"18"], cross_basic$pheno[,c("MAT_VOL2_nonzero")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL2_nonzero_scan_blup[3200:3850,],cross_basic$pmap)
# for(i in 1:nrow(covar_snp)){
#   x = which(names(geno_ML_marker) == covar_snp$Row.names[i])
#   covar_snp$ML_pos_geno[i] = geno_ML_marker[x]
# }



#try removing mice with bl6 haplotype for that snp
geno_MAT3_marker = maxmarg(apr, minprob=0.5,chr = 18, pos = 84.23055, map = cross_basic$pmap, return_char = T)
plot_pxg(geno_MAT3_marker,pheno = log10(cross_basic$pheno[,"MAT_VOL3"]+1),sort = TRUE)

#remove jsut the one that is "B" at 0.5 prob
drop = names(geno_MAT3_marker[which(geno_MAT3_marker == "B")])
ids = ind_ids(cross_basic)
ids = ids[-which(ids %in% drop)]

newcross = subset(cross_basic,ind = ids)

locus18_scan_noB = scan1(apr, newcross$pheno[,c("MAT_VOL1","MAT_VOL2","MAT_VOL3","MAT_VOL4", "MAT_VOL1_nonzero","MAT_VOL2_nonzero", "MAT_VOL3_nonzero","MAT_VOL4_nonzero")], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)

plot_scan1(locus18_scan[102700:103000,], map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3", col="red")
plot_scan1(locus18_scan_noB[102700:103000,], map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3", col="red")

plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL3_nonzero", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL1_nonzero", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL2_nonzero", add=T)
plot(locus18_scan_noB, map = cross_basic$pmap, chr = 18, lodcolumn = "MAT_VOL4_nonzero", add=T)


#thresh = 7.8
#abline(h=thresh, col="red")
peaks = find_peaks(locus18_scan, cross_basic$pmap, threshold=4, drop=1.5)
peaks_noB = find_peaks(locus18_scan_noB, cross_basic$pmap, threshold=4, drop=1.5)


#after dropping one mouse, qtl goes away for MAT_vol3, MAT_vol1 goes from 10 to just above thresh, and peaks at chr 2 increase slightly

MAT_VOL3_scan_blup_noB = scan1blup(apr[,"18"], newcross$pheno[,c("MAT_VOL3")], k_loco[["18"]], addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 4)
plot_coefCC(MAT_VOL3_scan_blup_noB[3200:3850,],newcross$pmap)



locus18_scan_noB2 = scan1(apr, newcross$pheno, k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)
peaks = find_peaks(locus18_scan_noB2, cross_basic$pmap, threshold=4, drop=1.5)


################
#drop all mice with "B" diplotype
geno_MAT3_marker = maxmarg(pr, minprob=0.5,chr = 18, pos = 84.23055, map = cross_basic$pmap, return_char = T)

drop = names(geno_MAT3_marker)[grepl("B",geno_MAT3_marker)]
drop = c(drop,"662")

ids = ind_ids(cross_basic)
ids = ids[-which(ids %in% drop)]

newcross = subset(cross_basic,ind = ids)

locus18_scan_noB_het = scan1(apr, newcross$pheno[,c("MAT_VOL3")], k_loco, Xcovar=Xcovar, addcovar = covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],cores = 1)

plot_scan1(locus18_scan_noB_het[102700:103000,], map = cross_basic$pmap, chr = 18, col="red")

peaks = find_peaks(locus18_scan_noB_het, cross_basic$pmap, threshold=4, drop=1.5)


######SAME SCAN ABOVE BUT FOR ENSMUSG00000048410 (Zfp407)
load(file = "./results/Rdata/cross_eqtl_oc.Rdata")

##
covar_eqtl = as.matrix(cross_eqtl$covar)
covar_eqtl[,"sex"] = (covar_eqtl[,"sex"] == "M")*1

covar_eqtl = covar_eqtl[,-1]#remove sac date as covar for now

covar_eqtl = apply(covar_eqtl,2,as.numeric) #make sure all cols are numeric
rownames(covar_eqtl) = rownames(cross_eqtl$covar)#make sure rownames match original cross file



locus18_scan = scan1(apr, cross_eqtl$pheno[,c("ENSMUSG00000048410")], k_loco, Xcovar=Xcovar, addcovar = covar_eqtl[,c(1,6:57)],cores = 1)
plot_scan1(locus18_scan[102700:103000,], map = cross_eqtl$pmap, chr = 18, col="red")
peaks = find_peaks(locus18_scan, cross_eqtl$pmap, threshold=4, drop=1.5)

start = 83.888746
end = 86.13701
chr = 18
out_snps_zfp407_18 <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,"ENSMUSG00000048410"], k_loco[["18"]],  addcovar =  covar_eqtl[,c(1,6:ncol(covar_eqtl))],Xcovar=Xcovar,
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

top_zfp407 <- top_snps(out_snps_zfp407_18$lod, out_snps_zfp407_18$snpinfo, drop = 0.15 * max(out_snps_zfp407_18$lod))
top_zfp407[order(top_zfp407$lod,decreasing = T),]
#top snp is BL6 private
plot_snpasso(out_snps_zfp407_18$lod, out_snps_zfp407_18$snpinfo, genes = genes_locus)

zfp407_blup_oc = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000048410")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:56)],cores = 4)
plot_coefCC(zfp407_blup_oc[3200:3850,],cross_eqtl$pmap)


#Zadh2 (ENSMUSG00000049090) in bone ob and oc and Tshz1 (ENSMUSG00000046982) in bone and ob
Zadh2_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000049090")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(Zadh2_blup[3200:3850,],cross_eqtl$pmap)


tshz1_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000046982")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(tshz1_blup[3200:3850,],cross_eqtl$pmap)

#
cndp2_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000024644")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(cndp2_blup[3200:3850,],cross_eqtl$pmap)

timm21_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000024645")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(timm21_blup[3200:3850,],cross_eqtl$pmap)

cyb5a_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000024646")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(cyb5a_blup[3200:3850,],cross_eqtl$pmap)

cndp1_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000056162")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(cndp1_blup[3200:3850,],cross_eqtl$pmap)

gm38576_blup = scan1blup(apr[,"18"], cross_eqtl$pheno[,c("ENSMUSG00000118075")], k_loco[["18"]], addcovar = covar_eqtl[,c(1,6:57)],cores = 4)
plot_coefCC(gm38576_blup[3200:3850,],cross_eqtl$pmap)


































##########################################################################################################
#Do merge analysis using eqtl
set.seed(8675309)
library(qtl2)
options(stringsAsFactors = F)



#####################
#Take all loci, get CI of locus +/- 250kb, get all eQTL within that locus, then compare top snps for phenos with top eqtl snps, define list of putatively causal genes

#qtl file
qtl_mat_18 = read.csv("./results/flat/qtl_peaks_MAT.csv", stringsAsFactors = FALSE)

#remove FFP and soleus and MAT vol1_nonzero
#qtl_norm = qtl_norm[-c(1:3),]
#remove MAT_vol1_nonzero
#qtl_norm = qtl_norm[-(which(qtl_norm$lodcolumn == "MAT_VOL1_nonzero")),]

#define loci
qtl_mat_18$locus = 0

loc_idx = 1
qtl_loc = as.data.frame(matrix(nrow=nrow(qtl_mat_18), ncol=ncol(qtl_mat_18)))
colnames(qtl_loc) = colnames(qtl_mat_18)
for (i in c(18)){
  sub = subset(qtl_mat_18, qtl_mat_18$chr == i)
  
  while(any(sub$locus == 0)){
    min_sub = min(sub$pos)
    idx = which((sub$pos >= min_sub) & (sub$pos <= min_sub + 1.5))
    sub[idx,"locus"] = loc_idx
    qtl_loc = rbind(qtl_loc, sub[idx,])
    sub = sub[-idx,]
    loc_idx = loc_idx + 1
  }
}

qtl_loc = qtl_loc[-which(is.na(qtl_loc)),]

#write.csv(qtl_loc, file = "./results/flat/qtl_loc_18", quote = FALSE,row.names = FALSE)

qtl_loc = read.csv("./results/flat/qtl_loc_18")
qtl_loc = qtl_loc[which(qtl_loc$locus==4),]






merge = list()
query_variants <- create_variant_query_func("../DO_project//data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("../DO_project/data/CCdb/mouse_genes_mgi.sqlite")



for(i in 1:nrow(qtl_loc)){
  print(i)
  merge[[i]] = list()
  pheno = qtl_loc$lodcolumn[i]
  chr = qtl_loc$chr[i]
  start = qtl_loc$ci_lo[i]
  end = qtl_loc$ci_hi[i]
  #use same covars as qtl mapping.sex,age,BW and gen
  out_snps <- scan1snps(pr, cross_basic$pmap, cross_basic$pheno[,pheno], k_loco[[chr]],  addcovar =  new_covar[,c("sex", "age_at_sac_days","body_weight","generationG24","generationG25","generationG26","generationG27","generationG28","generationG29","generationG30","generationG31","generationG32","generationG33","generationG34","generationG35")],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}

save(merge,file = "./results/Rdata/merge_QTL_MAT_18.Rdata")




#for each merge analysis, take the snps that are within 15% LODs of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)


save(merge_top,file = "./results/Rdata/merge_top_QTL_MAT_18.Rdata")











#for each locus, take genes that have eqtl , and are within locus CI +/- 250 kb
#any portion of gene start or end is within CI

#load eqtl
load("./results/Rdata/local_eqtl_ob.Rdata")

eqtl_loc = as.data.frame(matrix(nrow=nrow(local_eqtl_ob), ncol=ncol(local_eqtl_ob)+1))
colnames(eqtl_loc) = colnames(local_eqtl_ob)
colnames(eqtl_loc)[13] = "locus"

for(i in unique(qtl_loc$locus)){
  print(i)
  sub = subset(qtl_loc, qtl_loc$locus == i)
  min_loc = min(sub$ci_lo) - 0.25
  max_loc = max(sub$ci_hi) + 0.25
  
  sub_eqtl_1 = subset(local_eqtl_ob, local_eqtl_ob$chr == unique(sub$chr) & (as.numeric(local_eqtl_ob$Start)/1000000 >= min_loc) & as.numeric(local_eqtl_ob$Start)/1000000 <= max_loc)
  sub_eqtl_2 = subset(local_eqtl_ob, local_eqtl_ob$chr == unique(sub$chr) & (as.numeric(local_eqtl_ob$End)/1000000 >= min_loc) & as.numeric(local_eqtl_ob$End)/1000000 <= max_loc)
  
  sub_eqtl = merge(sub_eqtl_1, sub_eqtl_2)
  sub_eqtl = unique(sub_eqtl)
  sub_eqtl$locus = i
  
  eqtl_loc = rbind(eqtl_loc, sub_eqtl)
  
}
eqtl_loc = eqtl_loc[-which(is.na(eqtl_loc)),]



####
####
#load merge analysis objects

load("./results/Rdata/merge_top_local_eqtl_ob.Rdata")
merge_top_eqtl = merge_top

load("./results/Rdata/merge_top_QTL.Rdata")
merge_top_qtl = merge_top

rm(merge_top)


#for each locus, take eqtl merge analyses in locus and colocalize with pheno merge analyses in locus


genes = c()
phenos_w_genes = c()
for(i in unique(qtl_loc$locus)){
  #print(i)
  sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
  sub_eqtl = subset(eqtl_loc, eqtl_loc$locus == i)
  
  phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
  
  eqtl_genes = paste0(sub_eqtl$lodcolumn, "_", unique(sub_eqtl$chr))
  
  
  for(j in phenos){
    for(k in eqtl_genes){
      print(k)
      if(any(merge_top_eqtl[[k]]$snp_id %in% merge_top_qtl[[j]]$snp_id)){
        print(k)
        genes = append(genes,k)
        phenos_w_genes = append(phenos_w_genes, j)
      }
    }
  }
}


unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
gene_names = unlist(strsplit(genes, "_"))[seq(from = 1, to = length(genes)*2, by=2)]
gene_names = unique(gene_names)

eqtl_loc[which(eqtl_loc$lodcolumn %in% genes),]


phenos_w_genes = as.data.frame(phenos_w_genes)
phenos_w_genes$gene = genes
phenos_w_genes$gene_name = NA

for(i in 1:nrow(phenos_w_genes)){
  phenos_w_genes$gene[i] = unlist(strsplit(phenos_w_genes$gene[i], "_"))[1]
  phenos_w_genes$gene_name[i] = eqtl_loc[which(eqtl_loc$lodcolumn == phenos_w_genes$gene[i]),"Gene.Name"]
}

#eqtl genes that are located within a phenotypic qtl and regulated by a local eqtl
phenos_w_genes




## number nonsynonymous variants in top qtl merge for a phenotype
nonsyn = list()
counter=1
for(i in unique(qtl_loc$locus)){
  #print(i)
  sub_qtl = subset(qtl_loc, qtl_loc$locus == i)
  
  phenos = paste0(sub_qtl$lodcolumn, "_", unique(sub_qtl$chr))
  
  
  for(j in phenos){
    counter=counter+1
    missense = (grep("missense",x=merge_top_qtl[[j]]$consequence))
    l = length(missense)
    z = merge_top_qtl[[j]]$snp_id[missense]
    print(paste0(j,":",l))
    nonsyn[[counter]] = z
    names(nonsyn)[counter] = j
  }
}

nonsyn_frame = do.call(rbind, lapply(nonsyn, as.data.frame))




#SIFT: RSIDs uploaded to https://useast.ensembl.org/Tools/VEP


