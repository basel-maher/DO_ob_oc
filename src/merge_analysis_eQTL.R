##merge analysis
##Do merge analysis for each of the local eQTL that pass threshold (map_local_and_distal_eqtl.R), in the confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_bone.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
#load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_eqtl)

#load local eQTL
load("./results/Rdata/local_eqtl_bone.Rdata")



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_eqtl$covar)#make sure rownames match original cross file



#get eqtl list, passed threshold
local_eqtl_bone
#define locus, 1Mbp
chroms = as.vector(unique(local_eqtl_bone$chr)) #chroms
y <- c(1:19, "X","Y","MT")
chroms = chroms[order(match(chroms, y))]



merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


for(i in 1:nrow(local_eqtl_bone)){
  print(i)
  merge[[i]] = list()
  pheno = local_eqtl_bone$lodcolumn[i]
  chr = local_eqtl_bone$chr[i]
  start = local_eqtl_bone$ci_lo[i]
  end = local_eqtl_bone$ci_hi[i]
  #use same covars as eqtl mapping. Sex and all 48 PEER factors
  out_snps <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,pheno], k_loco[[chr]],  addcovar = covar[,c(1,10:57)],Xcovar=Xcovar,
                                                   query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  colnames(merge[[i]]$lod) = pheno
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
                             
}

#large file, saved on cluster for now (not saved yet)
#save(merge,file = "./results/Rdata/merge_local_eqtl.Rdata")


#for each merge analysis, take the snps that are within 15% of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)

save(merge_top,file = "./results/Rdata/merge_top_local_eqtl_bone.Rdata")

#summarize
summary(unlist(lapply(merge_top, function(x) nrow(x))))

#merge_top_df = bind_rows(merge_top, .id = "column_label")
#write.csv(merge_top_df, file="./merge_top_local_eqtl_df.csv", quote = F,row.names = F)

#####################################################################

##merge analysis
##Do merge analysis for each of the local eQTL that pass threshold (map_local_and_distal_eqtl.R), in the confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_ob.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
#load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_eqtl)

#load local eQTL
load("./results/Rdata/local_eqtl_ob.Rdata")



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_eqtl$covar)#make sure rownames match original cross file



#get eqtl list, passed threshold
local_eqtl_ob
#define locus, 1Mbp
chroms = as.vector(unique(local_eqtl_ob$chr)) #chroms
y <- c(1:19, "X","Y","MT")
chroms = chroms[order(match(chroms, y))]



merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


for(i in 1:nrow(local_eqtl_ob)){
  print(i)
  merge[[i]] = list()
  pheno = local_eqtl_ob$lodcolumn[i]
  chr = local_eqtl_ob$chr[i]
  start = local_eqtl_ob$ci_lo[i]
  end = local_eqtl_ob$ci_hi[i]
  #use same covars as eqtl mapping. Sex and all 48 PEER factors
  out_snps <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,pheno], k_loco[[chr]],  addcovar = covar[,c(1,11:58)],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  colnames(merge[[i]]$lod) = pheno
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}

#large file, saved on cluster for now (not saved yet)
#save(merge,file = "./results/Rdata/merge_local_eqtl.Rdata")


#for each merge analysis, take the snps that are within 15% of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)

save(merge_top,file = "./results/Rdata/merge_top_local_eqtl_ob.Rdata")

#summarize
summary(unlist(lapply(merge_top, function(x) nrow(x))))

#merge_top_df = bind_rows(merge_top, .id = "column_label")
#write.csv(merge_top_df, file="./merge_top_local_eqtl_df.csv", quote = F,row.names = F)








#####################################################################

##merge analysis
##Do merge analysis for each of the local eQTL that pass threshold (map_local_and_distal_eqtl.R), in the confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_oc.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
#load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_eqtl)

#load local eQTL
load("./results/Rdata/local_eqtl_oc.Rdata")



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_eqtl$covar)#make sure rownames match original cross file



#get eqtl list, passed threshold
local_eqtl_oc
#define locus, 1Mbp
chroms = as.vector(unique(local_eqtl_oc$chr)) #chroms
y <- c(1:19, "X","Y","MT")
chroms = chroms[order(match(chroms, y))]



merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


for(i in 1:nrow(local_eqtl_oc)){
  print(i)
  merge[[i]] = list()
  pheno = local_eqtl_oc$lodcolumn[i]
  chr = local_eqtl_oc$chr[i]
  start = local_eqtl_oc$ci_lo[i]
  end = local_eqtl_oc$ci_hi[i]
  #use same covars as eqtl mapping. Sex and all 48 PEER factors
  out_snps <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,pheno], k_loco[[chr]],  addcovar = covar[,c(1,11:56)],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  colnames(merge[[i]]$lod) = pheno
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}

#large file, saved on cluster for now (not saved yet)
#save(merge,file = "./results/Rdata/merge_local_eqtl.Rdata")


#for each merge analysis, take the snps that are within 15% of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)

save(merge_top,file = "./results/Rdata/merge_top_local_eqtl_oc.Rdata")

#summarize
summary(unlist(lapply(merge_top, function(x) nrow(x))))

#merge_top_df = bind_rows(merge_top, .id = "column_label")
#write.csv(merge_top_df, file="./merge_top_local_eqtl_df.csv", quote = F,row.names = F)


############################################
#merge analysis for ob but add eya4 on chrom 10
#####################################################################

##merge analysis
##Do merge analysis for each of the local eQTL that pass threshold (map_local_and_distal_eqtl.R), in the confidence interval for each locus

#DOWNLOAD THE FILES. DOWNLOADED NOV. 14 2019
#download.file("https://ndownloader.figshare.com/files/18533342", "./data/CCdb/cc_variants.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609252", "./data/CCdb/mouse_genes_mgi.sqlite")
#download.file("https://ndownloader.figshare.com/files/17609261", "./data/CCdb/mouse_genes.sqlite")
#

set.seed(8675309)
library(qtl2)

#load the geno probs
load(file = "./results/Rdata/pr_basic_cleaned.Rdata")

#load the cross file 
load(file = "./results/Rdata/cross_eqtl_ob.Rdata")

#apr
load(file = "./results/Rdata/apr_basic_cleaned.Rdata")

#load kinship file. In this case, using LOCO but can use overall file too
load(file = "./results/Rdata/k_loco_basic_cleaned.Rdata") #LOCO
#load(file = "./results/Rdata/k_basic_cleaned.Rdata")
#get Xcovar
Xcovar <- get_x_covar(cross_eqtl)

#load local eQTL
load("./results/Rdata/local_eqtl_ob.Rdata")



#create covar object. This is different fron the covar in cross_basic in that sex is a factor (1 for males, 0 for females)
covar = as.matrix(cross_eqtl$covar)
covar[,"sex"] = (covar[,"sex"] == "M")*1

covar = covar[,-1]#remove sac date as covar for now

covar = apply(covar,2,as.numeric) #make sure all cols are numeric
rownames(covar) = rownames(cross_eqtl$covar)#make sure rownames match original cross file



#get eqtl list, passed threshold
local_eqtl_ob

#add eya4 for ob, LOD is just under threshold (10.2) but for bone its much smaller, 6.4

x = c("ENSMUSG00000010461", "10","23.94633","10.235654","23.112649","24.62903","Eya4","10","-","23102963","23350786","89523")

local_eqtl_ob = rbind(local_eqtl_ob, x)


#define locus, 1Mbp
chroms = as.vector(unique(local_eqtl_ob$chr)) #chroms
y <- c(1:19, "X","Y","MT")
chroms = chroms[order(match(chroms, y))]



merge = list()
query_variants <- create_variant_query_func("./data/CCdb/cc_variants.sqlite")
query_genes <- create_gene_query_func("./data/CCdb/mouse_genes_mgi.sqlite")


for(i in 1:nrow(local_eqtl_ob)){
  print(i)
  merge[[i]] = list()
  pheno = local_eqtl_ob$lodcolumn[i]
  chr = local_eqtl_ob$chr[i]
  start = as.numeric(local_eqtl_ob$ci_lo[i])
  end = as.numeric(local_eqtl_ob$ci_hi[i])
  #use same covars as eqtl mapping. Sex and all 48 PEER factors
  out_snps <- scan1snps(pr, cross_eqtl$pmap, cross_eqtl$pheno[,pheno], k_loco[[chr]],  addcovar = covar[,c(1,11:58)],Xcovar=Xcovar,
                        query_func=query_variants,chr=chr, start=start, end=end, keep_all_snps=TRUE)
  merge[[i]] = out_snps
  colnames(merge[[i]]$lod) = pheno
  names(merge)[i] = paste0(pheno,"_",chr)
  rm(out_snps)
  
}

#large file, saved on cluster for now (not saved yet)
#save(merge,file = "./results/Rdata/merge_local_eqtl.Rdata")


#for each merge analysis, take the snps that are within 15% of the max LOD
merge_top = list()
for(i in 1:length(merge)){
  print(i)
  merge_top[[i]] = list()
  merge_top[[i]] = top_snps(merge[[i]]$lod, merge[[i]]$snpinfo, drop=max(merge[[i]]$lod)*0.15)
  
}
names(merge_top) = names(merge)

save(merge_top,file = "./results/Rdata/merge_top_local_eqtl_ob_EYA4_chr10_added.Rdata")

#summarize
summary(unlist(lapply(merge_top, function(x) nrow(x))))

#merge_top_df = bind_rows(merge_top, .id = "column_label")
#write.csv(merge_top_df, file="./merge_top_local_eqtl_df.csv", quote = F,row.names = F)







