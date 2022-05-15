#Fix geneseek files
#library(qtl2convert) dont think i need this

# - There was a sample confusion. Re-genotyped samples have a .1 appended to their name, except 371.
# - In ./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport, change Sample ID 371 to 371.1

#skip 9 to skip the header lines in FinalReports
myfread_fr <- function(filename) data.table::fread(filename, data.table=FALSE,skip = 9)


###change 371 to 371.1 in /Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport.txt"
x = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20190228_FinalReport.txt")
x$`Sample ID`[which(x$`Sample ID`=="371")] = "371.1"

write.table(x = x, file = "./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt", quote = FALSE,row.names = FALSE,sep = "\t")
rm(x)
#MAKE SURE TO ADD HEADER BACK MANUALLY!!

#Header is the following, without hashes
# [Header]
# GSGT Version	2.0.2
# Processing Date	2/28/2019 11:49 AM
# Content		GigaMuga_11769261_A.bpm
# Num SNPs	143259
# Total SNPs	143446
# Num Samples	96
# Total Samples	96
# [Data]













### in last finalreport file (20200305) change 667_repeat and 727_repeat to 667.1 and 727.1
x = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20200305_FinalReport.txt")
x$`Sample ID`[which(x$`Sample ID`=="667_repeat")] = "667.1"
x$`Sample ID`[which(x$`Sample ID`=="727_repeat")] = "727.1"

#remove qsox1
x = x[-which(x$`Sample ID`=="Qsox_7+6_42"),]
x = x[-which(x$`Sample ID`=="Qsox_1bp_42"),]
x = x[-which(x$`Sample ID`=="Qsox_171_B49"),]
x = x[-which(x$`Sample ID`=="Qsox_171_B102"),]
x = x[-which(x$`Sample ID`=="Qsox_171_J16"),]
x = x[-which(x$`Sample ID`=="Qsox_171_J65"),]
x = x[-which(x$`Sample ID`=="Qsox_1300_B36"),]
x = x[-which(x$`Sample ID`=="Qsox_780_53"),]

write.table(x = x, file = "./data/GIGAMUGA/FinalReport_files/20200305_FinalReport_FIXED.txt", quote = FALSE,row.names = FALSE,sep = "\t")
rm(x)

#MAKE SURE TO ADD HEADER BACK MANUALLY!!

#Header is the following, without hashes
# [Header]
# GSGT Version	2.0.4
# Processing Date	3/6/2020 11:53 AM
# Content		GigaMuga_11769261_A.bpm
# Num SNPs	143259
# Total SNPs	143446
# Num Samples	88
# Total Samples 88
# [Data]




# - Fix ./data/GIGAMUGA/merged/Merged_Sample_Map.txt manually. Change the second 371 entry to 371.1, also change 667_repeat and 727_repeat to 667.1 and 727.1, and remove the qsox samples



# create merged_FinalReport.txt

# - Generate ./data/GIGAMUGA/merged/Merged_FinalReport.txt.
#These are GeneSeek FinalReport files
x = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20171110_FinalReport.txt")
x2 = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20180408_FinalReport.txt")
x3 = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Barghouthi_MURGIGV01_20181213_FinalReport.txt")
x4 = myfread_fr("./data/GIGAMUGA/FinalReport_files/20190228_FinalReport_FIXED.txt")
x5 = myfread_fr("./data/GIGAMUGA/FinalReport_files/Univ_of_Virginia_Al-Harghouthi_MURGIGV01_20200106_FinalReport.txt")
x6 = myfread_fr("./data/GIGAMUGA/FinalReport_files/20200305_FinalReport_FIXED.txt")

merged = rbind(x,x2,x3,x4,x5,x6)


write.table(x = merged, file = "./data/GIGAMUGA/merged/Merged_FinalReport.txt", quote = FALSE,row.names = FALSE,sep = "\t")

rm(x,x2,x3,x4,x5,x6,merged)



###
#ADD THIS HEADER TO THE FILE MANUALLY!! without hashes

# [Header]
# GSGT Version	2.0.2
# Processing Date	11/14/2017 10:35 AM
# Content		GigaMuga_11769261_A.bpm
# Num SNPs	143259
# Total SNPs	143446
# Num Samples	832
# Total Samples	832
# [Data]

#change index to be 1-n in Merged_Sample_Map
x = read.delim("./data/GIGAMUGA/merged/Merged_Sample_Map.txt")
i = c(1:nrow(x))
x$Index = i

write.table(x = x, file = "./data/GIGAMUGA/merged/Merged_Sample_Map.txt", quote = FALSE,row.names = FALSE,sep = "\t")
