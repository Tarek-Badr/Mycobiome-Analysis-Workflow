if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.9")


library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(tidyverse)




setwd("C:/Users/ngs-adm/Desktop/Fungome")

path <- "C:/Users/ngs-adm/Desktop/Fungome/Fastq_Surgery_fungome"

setwd("C:/Users/ngs-adm/Desktop/Fungome/Fastq_Surgery_fungome")

list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))

fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

#Inspect read quality profiles

plotQualityProfile(fnFs)

plotQualityProfile(fnRs)


filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names



#We'll use standard filtering parameters: maxN=0 (DADA2 requires no Ns), truncQ=2, rm.phix=TRUE and  maxEE=2. The maxEE parameter sets the maximum number of "expected errors" allowed in a read, which is a better filter than simply averaging quality scores.
#watch out with trunclen, reads have to overlap at the end, standar script 250,200, you havr to try out, maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,230),
#maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
#compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)


plotQualityProfile(filtFs)

plotQualityProfile(filtRs)

#Learn the Error Rates

errF <- learnErrors(filtFs, multithread=FALSE)

errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

# apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#Inspecting the returned dada-class object:

dadaFs[[1]]

#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #min overlap is 12 as standard, but can be adjusted

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#CMost of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?

#Construct sequence table

seqtab <- makeSequenceTable(mergers)  

seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]

dim(seqtab)
dim(seqtab2)
# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

#Track reads through the pipeline

getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)

Pipeline_Track = head(track)
write.csv(track, file = "Pipeline_Track_final.csv")
write.csv(sample.names, file = "sample.names.csv")

################################
#Assign taxonomy with unite
################################
#download database: https://unite.ut.ee/repository.php
#ref: https://plutof.ut.ee/#/doi/10.15156/BIO/1281567

unite.ref = "C:/Users/ngs-adm/Desktop/Fungome/sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021.fasta"

#tryRC = TRUE
taxa_unite <- assignTaxonomy(seqtab.nochim, unite.ref, tryRC = TRUE)

taxa.print <- taxa_unite  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


##################################
#The typical standard outputs from amplicon processing are a fasta file, a count table, and a taxonomy table. So here's one way we can generate those files from your DADA2 objects in R:
# giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:

asv_fasta_fungi <- c(rbind(asv_headers, asv_seqs))

write(asv_fasta_fungi, "asv_fasta_fungi.fa")

# count table:

asv_tab_fungi <- t(seqtab.nochim)
row.names(asv_tab_fungi) <- sub(">", "", asv_headers)
write.table(asv_tab_fungi, "asv_count_fungi.tsv", sep="\t", quote=F, col.names=NA)

fungi_count_table = asv_tab_fungi


# tax table:

asv_tax_unite <- taxa_unite
row.names(asv_tax_unite) <- sub(">", "", asv_headers)
write.table(asv_tax_unite, "asv_tax_unite.tsv", sep="\t", quote=F, col.names=NA)

######################### 

fungi_count_table = asv_tab_fungi
fungi_taxa_table = asv_tax_unite

#######################################################

theme_set(theme_bw())
library(DESeq2)
###################Analysis with DESEq2
#We read the output files from DADA2 and the sample information file
###################

count_tab <- OP_count_tab

tax_tab <- OP_taxa_tab

sample_info_tab <- OP_sampleinf_tab

#examine consistancy in order between cts col names and coldata rownames 
#(they have to be in the same order or Deseq2 won't work out)

all(rownames(sample_info_tab) %in% colnames(count_tab))
all(rownames(sample_info_tab) == colnames(count_tab))

###################################
#Make the Deseq Tables
###################################
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ OP_Status) 
##################################

#prefiltering low count ASVs (383 ASVs)
#filter out genes with less than 10 counts in total

keep = rowSums(counts(deseq_counts)) > 5

deseq_counts_filtered = deseq_counts[keep,] #391
deseq_counts_filtered #rownames/ASVs (369)  


############################
# making our phyloseq object with transformed table

physeq <- phyloseq(otu_table(count_tab, taxa_are_rows=T), 
               sample_data(sample_info_tab), 
               tax_table(tax_tab))

rank_names(physeq)
#########################################################################
table(tax_table(physeq)[, "Phylum"], exclude = NULL)

ps <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(ps),
               MARGIN = ifelse(taxa_are_rows(ps), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(ps),
                    tax_table(ps))

ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


length(get_taxa_unique(ps, taxonomic.rank = "Genus"))

######Abundance value transformation#######################################

psra = transform_sample_counts(ps, function(x){x / sum(x)})

plot_abundance = function(physeq,title = "Candida abundance",
                          Facet = "Genus", Color = "Genus"){
  p1f = subset_taxa(physeq, Genus %in% c("g__Candida"))
  mphyseq = psmelt(p1f)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq, mapping = aes_string(x = "OP_Status",y = "Abundance",
                                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + scale_y_log10()+
    theme(legend.position="none")
}


plotBefore = plot_abundance(ps,"")
plotAfter = plot_abundance(psra,"")

grid.arrange(nrow = 2,  plotBefore, plotAfter)

######################################

psOrd = subset_taxa(psra, Genus == "g__Candida")
plot_abundance(psOrd, Facet = "Genus", Color = NULL)

#######################################
#Beta diversity:  PCoA and PERMANOVA/ADONIS'
# generating and visualizing the PCoA with phyloseq

vst_pcoa <- ordinate(physeq, method="MDS", distance="euclidean")
vst_pcoa <- ordinate(physeq, method="NMDS", distance="bray")

eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale the axes according to their magnitude of separating apart the samples

plot_ordination(physeq, vst_pcoa, color="char") + 
  geom_point(size=1) + labs(col="type") + 
  geom_text(aes(label=rownames(sample_info_tab), hjust=0.3, vjust=-0.4)) + 
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) + ggtitle("NMDS") + 
  scale_color_manual(values=unique(sample_info_tab$OP_Status[order(sample_info_tab$OP_Status)])) + 
  theme(legend.position="none")



########################################
#Richness and diversity estimates
########################################

# first we need to create a phyloseq object using our un-transformed count table
sample_info_tab_phy <- sample_data(t(seqsamples_info))

count_tab_phy <- otu_table(count_tab_order, taxa_are_rows=T)

tax_tab_phy <- tax_table(tax_tab)    #here based on the GTDB database (we can retry it for RDP or SILVA)

ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

ASV_physeq = physeq 
# and now we can call the plot_richness() function on our phyloseq object

plot_richness(ASV_physeq, color="OP_Status", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "Simpson")) 

plot_richness(ASV_physeq, color="OP_Status", shape = "Disease", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "Simpson")) 

plot_richness(ASV_physeq, color="OP_Status", shape = "Disease", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "InvSimpson")) 

plot_richness(ASV_physeq, x="OP_Status", color="Disease", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "InvSimpson")) 


plot_richness(ASV_physeq, x="OP_Status", measures=c("Chao1", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_classic() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

plot_richness(ASV_physeq, x="OP_Status", measures=c("Chao1", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))


plot_richness(ASV_physeq, x="Disease", measures=c("Chao1", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

plot_richness(ASV_physeq, x="Disease", shape = "Disease", measures=c("Chao1", "Shannon", "Simpson")) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

plot_richness(ASV_physeq, x="OP_Status", color="Disease", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "InvSimpson"))+
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

plot_richness(ASV_physeq, x="OP_Status", title = "Alpha Diversity", measures=c("Chao1", "Shannon", "InvSimpson"))+
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))


# Violin plot

P = plot_richness(ASV_physeq, x="OP_Status", measures=c("Chao1", "Shannon", "InvSimpson")) 

P + geom_violin()+
 theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 45))


############richness statistics

rich = estimate_richness(ASV_physeq, measures = c("Chao1", "Shannon"))

wilcox.Shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(ASV_physeq)$OP_Status, 
                                       p.adjust.method = "BH")
wilcox.Chao1 <- pairwise.wilcox.test(rich$Chao1, 
                                     sample_data(ASV_physeq)$OP_Status, 
                                     p.adjust.method = "BH")

tab.Shannon <- wilcox.Shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.Shannon

group1  group2     p.adj
1 Pre_OP Post_OP 0.1142857


wilcox.Chao1 <- pairwise.wilcox.test(rich$Chao1, 
                                     sample_data(ASV_physeq)$OP_Status, 
                                     p.adjust.method = "BH")

tab.Chao1 <- wilcox.Chao1$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.Chao1

#######################################################################################
#Taxonomic summaries
#######################################################################################

#We now construct a phyloseq object directly from the dada2 outputs.

dna <- Biostrings::DNAStringSet(asv_seqs)

names(dna) <- taxa_names(physeq)

ps_unite <- merge_phyloseq(physeq, dna)

taxa_names(ps_unite) <- paste0("ASV", seq(ntaxa(ps_unite)))

ps_unite

ps = ps_unite
#########################################################################################
#Bar plot:
##############################################################################

top20 <- names(sort(taxa_sums(ps_unite), decreasing=TRUE))[1:20]  # adjust number to wished top ones

ps.top20 <- transform_sample_counts(ps_unite, function(OTU) OTU/sum(OTU))

ps.top20 <- prune_taxa(top20, ps.top20)

###




tops <- names(sort(taxa_sums(ps_unite), decreasing=TRUE))  # adjust number to wished top ones
ps.top <- transform_sample_counts(ps_unite, function(OTU) OTU/sum(OTU))
ps.top <- prune_taxa(tops, ps.top)
####


#plot_bar(ps.top20, x="Day", fill="family") + facet_wrap(~When, scales="free_x")  #Family in SILVA, family in Decipher

BarPlot_Phylum_unite = plot_bar(ps.top, fill="Phylum")  
BarPlot_Genus_unite = plot_bar(ps.top, fill="Genus")
BarPlot_Family_unite = plot_bar(ps.top, fill="Family")
BarPlot_spp_unite = plot_bar(ps.top20, fill="Species")

plot_bar(ps.top, x = "OP_Status", y =  "Abundance", title = "Abundance based on Genus", fill="Genus")
plot_bar(ps.top, y =  "Abundance", title = "Abundance based on Genus", fill="Genus")


plot_bar(ps.top, y =  "Abundance", title = "Abundance based on Genus", fill="Genus", facet_grid=~OP_Status)




############################################################################
#Heat Maps
############################################################################

theme_set(theme_bw())


Heatmap with unfiltered read

top50 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:50]
ps.top50 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top50 <- prune_taxa(top50, ps.top50)





rank_names(ps.top50)

plot_heatmap(ps.top50, sample.label="OP_Status","Family")

plot_heatmap(ps.top50, method = "NMDS", distance = "bray",
             sample.label = "Protocol", taxa.label = "Genus", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Genus (GTDB)", sample.order = NULL, taxa.order = "Genus",
             first.sample = NULL, first.taxa = NULL)

plot_heatmap(ps.top50, sample.label = "Protocol", taxa.label = "Genus", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Genus (GTDB)", sample.order = NULL, taxa.order = "Genus",
             first.sample = NULL, first.taxa = NULL)

plot_heatmap(ps.top50, taxa.label = "Species", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Species (GTDB)")

heatmap(otu_table(ps.top50))

p_HM <- plot_heatmap(ps_GTDB, "NMDS", "bray", "Protocol", "Family")
















