## Infant Microbiome Project ##
# Thomas J Moutinho Jr.

#### Install Local ####
source("https://bioconductor.org/biocLite.R")
biocLite("dada2")
biocLite("GenomicRanges")
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("randomForest")
install.packages("dplyr")
install.packages("data.table")
install.packages("readxl")
install.packages("tidyr")
install.packages("DESeq2")
install.packages("devtools")
install.packages('gdtools')
packageVersion("vegan")
devtools::install_version("vegan", version = "2.4-6")
install.packages("viridis")


#### Libraries ####
library(viridis)
library(dada2)
packageVersion("dada2")
library(randomForest)
library(data.table)
library('readxl')
library(tidyr)
library(ggplot2)
library(DESeq2)
library(vegan)
library(fastR)
library(mice)
library(tibble)
library(phyloseq)
library(gridExtra)
library(gdtools)
library(AUCRF)
packageVersion("AUCRF")
library(grid)

#### DADA2 Data Processing ####
path <- "D:/Google/Infant Microbiome Project/NICU/fastq_files_2" # CHANGE ME to the directory containing the fastq files after unzipping.

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
# sample.names %>% View()

plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,210),
                     maxN=0, maxEE=c(4,5), truncQ=3, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE
head(out)

set.seed(42)
errF <- learnErrors(filtFs, multithread=FALSE)
set.seed(42)
errR <- learnErrors(filtRs, multithread=FALSE)

p <- plotErrors(errF, nominalQ=TRUE)

## Big Data Work Through ##

filtPATH <- "D:/Google/Infant Microbiome Project/NICU/fastq_files_2/filtered"

filtFs <- list.files(filtPATH, pattern="_F_filt.fastq.gz", full.names = TRUE)
filtRs <- list.files(filtPATH, pattern="_R_filt.fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1) # Assumes filename = samplename_XXX.fastq.gz
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=FALSE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=FALSE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)

# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
saveRDS(seqtab, "seqtab.rds") # CHANGE ME to where you want sequence table saved


# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)
save(seqtab.nochim, file = "seqtab.nochim.Rda")
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

taxa <- assignTaxonomy(seqtab.nochim, "rdp_train_set_16.fa.gz", multithread=FALSE)
save(taxa, file = "taxa.Rda")

library(RSQLite)
library(DECIPHER)

seqs <- dada2::getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)
phangAlign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
dm <- phangorn::dist.ml(phangAlign)
treeNJ <- phangorn::NJ(dm) # Note, tip order != sequence order
fit = phangorn::pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- phangorn::optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                              rearrangement = "stochastic", control = phangorn::pml.control(trace = 0))


#### Load DADA2 Files ####
load("seqtab.nochim.Rda")
load("taxa.Rda")

#### Rarefaction Curves ####
raremax <- min(rowSums(seqtab.nochim))
rarecurve(seqtab.nochim[11:20,], step = 20, sample = raremax, col = "blue", cex = 0.6, ylab = "Unique Reads")

#### Read in MetaData ####
meta.data <- read.csv("Data/final clinical phenotyping edited 5-1-18.csv")
meta.data <- as.data.frame(meta.data)
meta.data <- meta.data[,c("Subject","Twin_Set","gest_age","gender","birth_mode","sample_id","age_at_samp",
                          "weight_kg","hyperbili","bilirubin","TPN","PO","Bolus","Continuous","cholestasis")]
meta.data$Subject <- as.factor(meta.data$Subject)
meta.data$Twin_Set <- as.factor(meta.data$Twin_Set)
meta.data$gest_age <- as.numeric(meta.data$gest_age)
meta.data$gender <- as.factor(meta.data$gender)
meta.data$birth_mode <- as.factor(meta.data$birth_mode)
meta.data$sample_id <- as.factor(meta.data$sample_id)
meta.data$age_at_samp <- as.numeric(meta.data$age_at_samp)
meta.data$weight_kg <- as.numeric(as.character(meta.data$weight_kg))
meta.data$hyperbili <- as.factor(meta.data$hyperbili)
meta.data$TPN <- as.factor(meta.data$TPN)
meta.data$PO <- as.factor(meta.data$PO)
meta.data$Bolus <- as.factor(meta.data$Bolus)
meta.data$Continuous <- as.factor(meta.data$Continuous)
meta.data$cholestasis <- as.factor(meta.data$cholestasis)
meta.data$bilirubin <- as.numeric(as.character(meta.data$bilirubin))

#### Discordant Twins Analysis ####
meta_twins <- meta.data
rownames(meta_twins) <- meta_twins$sample_id
meta_twins <- meta_twins[meta_twins$cholestasis %in% c("Y","N"),]

p = ggplot(meta_twins, aes(x = age_at_samp, y = Subject, color = gender)) + 
  geom_point(size = 2) + coord_cartesian(xlim = c(0,120)) +
  labs(title = "", color = "Gender") +
  xlab("Age of Subject at Time of Sample [Days]")+
  ylab("")+
  scale_y_discrete(labels=c('Twin 1a','Twin 1b','Twin 2a','Twin 2b','Twin 3a','Twin 3b','Twin 4a','Twin 4b'))
p

samdf <- meta_twins
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))
ps <- prune_samples(sample_sums(ps) > 5000, ps)
psg <- tax_glom(ps, "Genus", NArm = TRUE)

pso <- tax_glom(ps, "Order", NArm = FALSE)
psp <- tax_glom(ps, "Phylum", NArm = FALSE)

ps <- filter_taxa(ps, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
ps_relab <- transform_sample_counts(ps, function(OTU) 100*(OTU/sum(OTU)))

psg_original <- psg
psg <- filter_taxa(psg, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
psg_relab <- transform_sample_counts(psg, function(OTU) 100*(OTU/sum(OTU)))

psg_relab_temp <- transform_sample_counts(psg, function(OTU) 100*(OTU/sum(OTU)))
psg_relab_log2 <- transform_sample_counts(psg_relab_temp, function(x) log2(1 + x))

# Order and Phylum level does to allow for even bars in plot
pso <- filter_taxa(pso, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
pso_relab <- transform_sample_counts(pso, function(OTU) 100*(OTU/sum(OTU)))

psp <- filter_taxa(psp, function(x) sum(x > 3) > (0.05*length(x)), TRUE)
psp_relab <- transform_sample_counts(psp, function(OTU) 100*(OTU/sum(OTU)))

## NMDS Plots

out.bray.log <- ordinate(psg_relab, method = "NMDS", distance = "bray")
p1 = plot_ordination(psg_relab, out.bray.log,  color = "Twin_Set") + 
  labs(title = "") + theme_bw() + stat_ellipse(type = "t", size = 1) +
  scale_color_manual(name = 'Twins',labels = c('Twin Set 1','Twin Set 2','Twin Set 3','Twin Set 4'), values = c('#00BFC4','#7CAE00','#F8766D','#C77CFF'))
nmds_p2 = plot_ordination(psg_relab, out.bray.log, color = "cholestasis") + 
  labs(title = "") + theme_bw() + stat_ellipse(type = "t", size = 1) +
  scale_color_manual(name = 'PNAC',labels = c('No','Yes'), values = c('#00BFC4','#F8766D'))
grid.arrange(nmds_p2, p1, nrow = 1)

## NMDS Stats
set.seed(42)
# Calculate bray curtis distance matrix using phyloseq
erie_bray <- phyloseq::distance(psg_relab, method = "bray")
# make a data frame from the sample_data
sampledf <- data.frame(sample_data(psg_relab))
# Adonis test
adonis(erie_bray ~ cholestasis, data = sampledf)
adonis(erie_bray ~ gender, data = sampledf)
adonis(erie_bray ~ hyperbili, data = sampledf)
adonis(erie_bray ~ Twin_Set, data = sampledf)

## Alpha Diversity Figures

alpha <- estimate_richness(psg_original, split = TRUE, measures=c("Shannon", "Simpson")) %>% 
  merge(meta_twins[,c("sample_id","cholestasis","Subject","age_at_samp","bilirubin")], by = "row.names")
alpha$age_at_samp <- as.factor(alpha$age_at_samp)
alpha$Subject <- factor(alpha$Subject)
alpha$Subject_orig <- alpha$Subject
levels(alpha$Subject) <- c('Twin 1a','Twin 1b','Twin 2a','Twin 2b',
                           'Twin 3a','Twin 3b','Twin 4a','Twin 4b')
alpha$age_at_samp <- as.numeric(levels(alpha$age_at_samp))[alpha$age_at_samp]

theme_set(theme_classic())
alph_p = ggplot(alpha, aes(x=age_at_samp, y=Simpson, color = cholestasis)) +
  coord_cartesian(ylim = c(0,0.8)) + # ,xlim = c(7,104)
  geom_point(size = 2.5, alpha = 0.7) +
  facet_wrap(~Subject, scales="free", ncol = 2) +
  labs(x = "Age at Sample [Days]", y = "Simpson Alpha Diversity", title = "") +
  scale_color_manual(name = 'PNAC',labels = c('No','Yes'), values = c('#00BFC4','#F8766D')) +
  theme(strip.background = element_blank()) #element_rect(fill="white")
alph_p

alpha_out <- subset(alpha, select= -c(Shannon,Row.names))
write.csv(alpha_out, file = "Data/Figure_2_data.csv")

## Log 2 NMDS Figures
sv_data <- otu_table(psg_relab)@.Data %>% as.data.frame()

names <- taxa[colnames(sv_data),c("Phylum","Genus")] %>% as.data.frame() %>% .$Genus %>% as.character()

colnames(sv_data) <- names
sv_data <- sv_data %>% merge(meta_twins[,c("sample_id","cholestasis")], by = "row.names")

sv_data_long <- gather(sv_data, Sequence_Variant, Abundance, 2:19, factor_key = TRUE)


theme_set(theme_bw())
box_p1 = ggplot(data = sv_data_long, aes(x = Sequence_Variant, y = Abundance, color = cholestasis)) + 
  geom_boxplot(aes(x = Sequence_Variant, y = Abundance, color = cholestasis), outlier.size = .5, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  coord_cartesian(ylim = c(0,100)) +
  scale_color_manual(name = 'PNAC',labels = c('No','Yes'), values = c('#00BFC4','#F8766D')) +
  labs(x = "", y = "Relative Abundance", title = "") +
  scale_x_discrete(labels=c("Clostridium_sensu_stricto" = "Clostridium Sensu Stricto","Clostridium_XI"="Clostridium XI"))
box_p1

## Boxplot for Log Transformed Data: Significant Genera when comparing Twins with Cholestasis##
sv_data <- otu_table(psg_relab_log2)@.Data %>% as.data.frame()

names <- taxa[colnames(sv_data),c("Phylum","Genus")] %>% as.data.frame() %>% .$Genus %>% as.character()

colnames(sv_data) <- names
sv_data <- sv_data %>% merge(meta_twins[,c("sample_id","cholestasis")], by = "row.names")

sv_data_long <- gather(sv_data, Sequence_Variant, Abundance, 2:19, factor_key = TRUE)

theme_set(theme_bw())
box_p2 = ggplot(data = sv_data_long, aes(x = Sequence_Variant, y = Abundance, color = cholestasis)) + 
  geom_boxplot(aes(x = Sequence_Variant, y = Abundance, color = cholestasis), outlier.size = .5, position=position_dodge(width = 1)) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
  scale_color_manual(name = 'PNAC',labels = c('No','Yes'), values = c('#00BFC4','#F8766D')) +
  labs(x = "Genera", y = expression("Log"[2]*"(Relative Abundance)"), title = "") +
  scale_x_discrete(labels=c("Clostridium_sensu_stricto" = "Clostridium Sensu Stricto","Clostridium_XI"="Clostridium XI"))
box_p2

pdf("Figures/Figure_1.pdf", width = 11, height = 5)
grid.newpage()
print(nmds_p2, vp = viewport(x = 0.21, y = 0.6, width = 0.4, height = 0.75))
print(box_p2, vp = viewport(x = 0.7, y = 0.49, width = 0.6, height = 1))
print(grid.text("A", vp = viewport(x = 0.02, y = .95), gp=gpar(fontsize=16)))
print(grid.text("B", vp = viewport(x = 0.41, y = .95), gp=gpar(fontsize=16)))
print(grid.text("*", vp = viewport(x = 0.399, y = .082), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.4815, y = .177), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.496, y = .141), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.533, y = .17), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.551, y = .146), gp=gpar(fontsize=12)))
dev.off()

## Genera Subset Relative Abundance ## Relative Abundance Time Series
genera_in_dataset <- psg_relab@tax_table@.Data[,"Genus"]
genera_in_dataset <- unname(genera_in_dataset)
top_genera <- genera_in_dataset[c(1,2,3,4,5)]
top_gen_logic <- c(1,2,3,4,5)
other_genera <- genera_in_dataset[!(genera_in_dataset %in% top_genera)]
other_gen_logic <- c(1:length(genera_in_dataset))[!(genera_in_dataset %in% top_genera)]

psg_relab_merge_other <- merge_taxa(psg_relab, other_gen_logic, archetype = 1)

psg_relab_merge_other@tax_table@.Data[other_gen_logic[1],"Genus"] <- "Other"

levels(psg_relab_merge_other@sam_data[["Subject"]]) <- c('Twin 1a','Twin 1b','Twin 2a','Twin 2b',
                                                         'Twin 3a','Twin 3b','Twin 4a','Twin 4b')

# Start by melting the data in the "standard" way using psmelt.
mdf = psmelt(psg_relab_merge_other)

# Build the plot data structure
theme_set(theme_classic())
p = ggplot(mdf, aes_string(x="age_at_samp", y="Abundance"))
bar_p2 = p + geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack") +
  facet_wrap(~Subject, scales="free", ncol = 2) +
  labs(x = "Age at Sample [Days]", y = "Relative Abundance", title = "") +
  theme( strip.background = element_blank()) 
bar_p2

pdf("Figures/Figure_2.pdf", width = 12, height = 6)
grid.newpage()
print(alph_p, vp = viewport(x = 0.225, y = 0.5, width = 0.44, height = 1))
print(bar_p2, vp = viewport(x = 0.725, y = 0.5, width = 0.56, height = 1))
print(grid.text("A", vp = viewport(x = 0.02, y = .95), gp=gpar(fontsize=16)))
print(grid.text("B", vp = viewport(x = 0.46, y = .95), gp=gpar(fontsize=16)))
print(grid.text("*", vp = viewport(x = 0.101, y = .88), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.079, y = .645), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.294, y = .41), gp=gpar(fontsize=12)))
print(grid.text("*", vp = viewport(x = 0.119, y = .185), gp=gpar(fontsize=12)))
dev.off()

## Wilcoxon Rank Sum Tests Log2 Data##
genus_data <- otu_table(psg_relab_log2)@.Data %>% as.data.frame()
names <- taxa[colnames(genus_data),c("Phylum","Genus")] %>% as.data.frame() %>% .$Genus %>% as.character()
colnames(genus_data) <- names
genus_data <- merge(genus_data, samdf[,c("cholestasis","gender")], by = "row.names")
rownames(genus_data) <- genus_data$Row.names
genus_data$Row.names <- NULL
genus_data$gender <- NULL

p_vals <- c()
p_vals[1] <- wilcox.test(`Escherichia/Shigella` ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[2] <- wilcox.test(Staphylococcus ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[3] <- wilcox.test(Klebsiella ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[4] <- wilcox.test(Enterococcus ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[5] <- wilcox.test(Veillonella ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[6] <- wilcox.test(Enterobacter ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()
p_vals[7] <- wilcox.test(Streptococcus ~ cholestasis, data=genus_data) %>% .["p.value"] %>% as.numeric()

p_vals_adj_log2 <- p.adjust(p_vals, method = "BH")

## AUCRF ##
genus_data <- otu_table(psg_relab)@.Data %>% as.data.frame()
names <- taxa[colnames(genus_data),c("Phylum","Genus")] %>% as.data.frame() %>% .$Genus %>% as.character()
colnames(genus_data) <- names
genus_data <- merge(genus_data, samdf[,c("cholestasis","gender")], by = "row.names")
rownames(genus_data) <- genus_data$Row.names
genus_data$Row.names <- NULL
genus_data$gender <- NULL
column_names <- c("Escherichia/Shigella","Staphylococcus","Klebsiella","Enterococcus","Veillonella",
                  "Enterobacter","Streptococcus","cholestasis")
genus_data_meta <- genus_data[,column_names]
colnames(genus_data_meta) <- c("Escherichia","Staphylococcus","Klebsiella","Enterococcus","Veillonella",
                               "Enterobacter","Streptococcus","cholestasis")

genus_data_rf_meta <- mutate(genus_data_meta, cholestasis = ifelse(cholestasis == "N",0,1))
ratio_w_cholestasis <- sum(genus_data_rf_meta$cholestasis)/length(genus_data_rf_meta$cholestasis)
genus_data_rf_meta$cholestasis <- as.factor(genus_data_rf_meta$cholestasis)

set.seed(42)
fit <- AUCRF(cholestasis~., data = genus_data_rf_meta, ranking=c("MDG"), classwt = c(10,1), ntree = 1000, mtry = 2)
fitcv <- AUCRFcv(fit)
OptimalSet(fitcv)
plot(fitcv, ylim=c(0,1.2))
fitcv[["RFopt"]][["confusion"]]


df = data.frame(x = c(fitcv[["AUCcurve"]][["k"]]), y = c(fitcv[["AUCcurve"]][["AUC"]]))
theme_set(theme_bw())
auc_p = ggplot(df, aes(x = x, y = y)) + geom_line(color = "black", size = 1.1) + 
  geom_point(size = 2) + geom_point(x=5, y=0.90, color = "red",size=5,alpha=0.2) +
  coord_cartesian(ylim = c(0.3,1),xlim=c(1,7)) +
  scale_x_continuous(breaks=seq(1,7,1)) +
  annotate("text",x=5.5,y=0.78,label=paste0("cvAUC = 0.88\nOOB-AUCopt = 0.90")) +
  labs(x = "Number of Features in Model", y = "OOB-AUC", title = "")
auc_p

Genus <- c("Veillonella","Klebsiella","Enterobacter","Staphylococcus","Streptococcus")

df_table <- data.frame(Genus)
table_auc <- as.data.table(df_table)

pdf("Figures/Figure_3.pdf", width = 7, height = 3)
grid.newpage()
print(auc_p, vp = viewport(x = 2.5/7, y = 0.5, width = 5/7, height = 1))
print(grid.table(table_auc,vp = viewport(x = 6/7, y = .55, width = 2/7, height = 1), 
                 theme = ttheme_minimal(core=list(fg_params=list(hjust=0, x=0.1)),
                                        colhead=list(fg_params=list(hjust=0, x=0.1)))))
print(grid.text("A", vp = viewport(x = 0.02, y = .9), gp=gpar(fontsize=15)))
print(grid.text("B", vp = viewport(x = 5.2/7, y = .9), gp=gpar(fontsize=15)))
dev.off()


