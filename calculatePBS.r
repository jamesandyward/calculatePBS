library(data.table)
library(dplyr)
library(qqman)

##### Load in Data #####
AFT_EAT <- fread("AFT_EAT_50kb_2kb.windowed.weir.fst", header = T)
AFT_EAT$T_PBS <-  -log(1 - AFT_EAT$WEIGHTED_FST)
colnames(AFT_EAT) <- c("CHROM", "BIN_S", "BIN_E", "N_VARIANTS", "AE_W_FST", "NE_M_FST", "T_PBS")
AFT_EAT <- AFT_EAT %>% mutate(ID = paste(CHROM, BIN_S, BIN_E, sep = "_"))

AFT_AAI <- fread("AFT_AAI_50kb_2kb.windowed.weir.fst", header = T)
AFT_AAI$T_PBS <-  -log(1 - AFT_AAI$WEIGHTED_FST)
colnames(AFT_AAI) <- c("CHROM", "BIN_S", "BIN_E", "N_VARIANTS", "AA_W_FST", "NA_M_FST", "T_PBS")
AFT_AAI <- AFT_AAI %>% mutate(ID = paste(CHROM, BIN_S, BIN_E, sep = "_"))

EAT_AAI <- fread("EAT_AAI_50kb_2kb.windowed.weir.fst", header = T)
EAT_AAI$T_PBS <- -log(1 - EAT_AAI$WEIGHTED_FST)
colnames(EAT_AAI) <- c("CHROM", "BIN_S", "BIN_E", "N_VARIANTS", "EA_W_FST", "EA_M_FST", "T_PBS")
EAT_AAI <- EAT_AAI %>% mutate(ID = paste(CHROM, BIN_S, BIN_E, sep = "_"))

##### Merge into a single dataframe #####
AFT_PBS <- AFT_EAT %>% inner_join(AFT_AAI, by = "ID") %>% inner_join(EAT_AAI, by = "ID")
AFT_PBS$PBS <- (AFT_PBS$T_PBS.x + AFT_PBS$T_PBS.y - AFT_PBS$T_PBS)/2
AFT_PBS$CHROM <- as.numeric(as.character(AFT_PBS$CHROM))
chromosome_subset <- AFT_PBS[AFT_PBS$CHROM >= 1 & AFT_PBS$CHROM <= 29, ]
posPBS <- chromosome_subset %>% filter(PBS > 0 & PBS < 2)
posPBS <- posPBS %>% filter(NE_M_FST > 0.1)
posPBS <- posPBS %>% filter(!is.infinite(PBS))
posPBS$CHROM <- factor(posPBS$CHROM, levels = as.character(c(1:29)))
posPBS$CHROM <- as.numeric(posPBS$CHROM)
top_0.1_percentile_value <- quantile(posPBS$PBS, probs = 0.999)
write.table(posPBS, "Final_datasets/AFT_PBS2.txt", col.names = T, row.names = F, quote = F)

##### Plot Manhattan plot #####
png("AFT_PBS2.png", width = 20, height = 5, units = "in", res = 300)
manhattan(posPBS, chr="CHROM", bp="BIN_E", snp="ID", p="PBS", logp = FALSE, ylab = "PBS", genomewideline = top_0.1_percentile_value,  ylim = c(0, 1.5), col = c("#101071", "#f09b0c"))
dev.off()

pdf("AFT_PBS.pdf", width = 20, height = 5)
manhattan(posPBS, chr="CHROM", bp="BIN_E", snp="ID", p="PBS", logp = FALSE, ylab = "PBS", genomewideline = top_0.1_percentile_value,  ylim = c(0, 0.8), col = c("#101071", "#f09b0c"))
dev.off()
