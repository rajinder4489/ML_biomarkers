rm(list=ls())
setwd("C:/Users/r.gupta/Desktop/livercancer/")

#biomart <- read.delim("biomart_ensembl_to_entrez-gene.txt", sep="\t", stringsAsFactors=F)
biomart <- read.delim("biomart_desc.txt", sep="\t", stringsAsFactors=F)
expr <- read.delim("isoform_norm_reads.txt", sep="\t", stringsAsFactors = F)

biomarkers <- c("AFP", "HPAFP", "ANXA2", "ANX2", "ANX2L4", "CAL1H", "LPC2D", "ANXA7", "ANX7", "SNX", "OK/SW-cl.95", "PROM1", 
                "PROML1", "MSTP061", "CD44", "LHR", "MDU2", "MDU3", "MIC4", "THY1", "KRT19", "EPCAM", "GA733-2", "M1S2", 
                "M4S1", "MIC18", "TACSTD1", "TROP1", "FGF", "GOLM1", "C9orf155", "GOLPH2", "PSEC0242", "UNQ686/PRO1326", 
                "GPC3", "OCI5", "HGF", "HPTA", "HIF1A", "BHLHE78", "MOP1", "PASD8", "HSP90AA1", "HSP90A", "HSPC1", "HSPCA", 
                "HSP90AB1", "HSP90B", "HSPC2", "HSPCB", "HSPA12A", "KIAA0417", "HSPA12B", "C20orf60", "HSPA13", "STCH", 
                "HSPA14", "HSP60", "HSP70L1", "HSPA1A", "HSP72", "HSPA1", "HSX70", "HSPA1B", "HSPA1L", "HSPA2", "HSPA4", 
                "APG2", "HSPA4L", "APG1", "OSP94", "HSPA5", "GRP78", "HSPA6", " HSP70B", "HSPA7", "HSPA8", "HSC70", 
                "HSP73", "HSPA10", "HSPA9", "GRP75", "HSPA9B", "mt-HSP70", "MDK", "MK1", "NEGF2", "SPP1", "BNSP", "OPN", 
                "PSEC0156", "SERPINB3", "SCCA", "SCCA1", "TGFB1", "TGFB", "TGFB2", "TGFB3", "TLN1", "KIAA1027", "TLN", 
                "VEGFA", "VEGF", "VEGFBVRF", "VEGFC", "VEGFD", "FIGF")

hits_true <- biomart[biomart$Gene.name %in% biomarkers,]

hits_prob <- biomart[grep(paste("^", paste(biomarkers, collapse="|^"), sep=""), biomart$Gene.name, 
                          value=FALSE, ignore.case = TRUE, perl = TRUE),]

hits_prob <- biomart[grep(paste("^", "FGF[0-9]*$", sep=""), biomart$Gene.name, 
                          value=FALSE, ignore.case = TRUE, perl = TRUE),] #because FGF had way too many types

hits_all <- unique(rbind(hits_true, hits_prob))
  
write.table(hits_all, file="biomarkers_mapped_gid.txt", sep="\t", row.names=F)

#hits_to_map <- hits_all[hits_all$Transcript.type == "protein_coding", 
 #                       c("Transcript.stable.ID", "Gene.name", "Transcript.name")]

hits_to_map <- hits_all[, c("Transcript.stable.ID", "Gene.name", "Transcript.name")]

mapped_expr <- expr[expr$id %in% hits_to_map$Transcript.stable.ID,]
write.table(mapped_expr, file="biomarkers_expr.txt", sep="\t", row.names=F)

expr_final <- merge(mapped_expr, hits_to_map, by.x="id", by.y="Transcript.stable.ID")
write.table(expr_final, file="biomarkers_expr_extra.txt", sep="\t", row.names=F)

