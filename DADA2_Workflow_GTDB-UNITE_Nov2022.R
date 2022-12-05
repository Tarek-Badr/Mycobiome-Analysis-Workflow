
####Introduction####


#### End ####

#### Package setup ####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
install.packages("remotes")
remotes::install_github("jfq3/ggordiplots")

library(dada2); packageVersion("dada2")
library(knitr)
library(BiocStyle)
library(ggplot2)
library(gridExtra)
library(msa)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(tidyverse)
library(vegan)
library(ComplexHeatmap)
library(DESeq2)
library(plyr)
library(ggordiplots)

#### End ####

#### Set Directory ####
C:\Users\ngs-adm\Desktop\Surgery_MC_cohort\Fungome_DADA2_Analysis\Fastq

path <- "C:/Users/ngs-adm/Desktop/Simon_Fungome/Surgery_MC_cohort_tareknov22/Fungome_DADA2_Analysis/Fastq"
setwd("C:/Users/ngs-adm/Desktop/Simon_Fungome/Surgery_MC_cohort_tareknov22/Fungome_DADA2_Analysis")

list.files(path)

theme_set(theme_bw())
#### End ####

#### Quality Control and Filtering####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq

fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))


#fnFs <- sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
#fnRs <- sort(list.files(path, pattern="_2.fastq", full.names = TRUE))

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
#watch out with trunclen, reads have to overlap at the end, standar script 250,200, you have to try out, maxEE can be eased maxEE=c(2,5) if too many read are lost because of low quality
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,230), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(230,230),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE


#For Fungi No truncLen but minDLen due to fungal ITS length variations#

#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(2, 5), 
#                    truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE


head(out)


plotQualityProfile(filtFs)

plotQualityProfile(filtRs)
#### End ####

####Learn the Error Rates####

errF <- learnErrors(filtFs, multithread=FALSE)

errR <- learnErrors(filtRs, multithread=FALSE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)
#### End ####

#### Dada and Merge####
# apply the core sample inference algorithm to the dereplicated data.

dadaFs <- dada(filtFs, err=errF, multithread=FALSE)

dadaRs <- dada(filtRs, err=errR, multithread=FALSE)

#Inspecting the returned dada-class object:

dadaFs[[1]]

#Merge paired reads

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE) #min overlap is 12 as standard, but can be adjusted

# Inspect the merger data.frame from the first sample

head(mergers[[1]])

#Most of your reads should successfully merge. If that is not the case upstream parameters may need to be revisited: Did you trim away the overlap between your reads?
#Construct sequence table

seqtab <- makeSequenceTable(mergers)  
dim(seqtab)
#If you wish to filter very short or long merged reads:-
#seqtab2 <- seqtab[,nchar(colnames(seqtab)) %in% 300:387]
#dim(seqtab2)

# Inspect distribution of sequence lengths

table(nchar(getSequences(seqtab)))

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=FALSE, verbose=TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)
#### End ####


#### Track Pipeline ####
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
#### End ####

#Assign taxonomy with unite/GTDB#
#### Assign taxonomy for Bacteria with GTDB####
#Analaysis of 16S rLaparoscopic.open Data using GTDB:-

GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz
GTDB_bac120_arc122_ssu_r202_Species.fa.gz

ref_fasta = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc122_ssu_r202_fullTaxo.fa.gz"

#tryRC = TRUE

taxa_GTDB <- assignTaxonomy(seqtab.nochim, ref_fasta,tryRC = TRUE )
unname(taxa_GTDB)
taxa.print_GTDB  <- taxa_GTDB # Removing sequence rownames for display only
rownames(taxa.print_GTDB) <- NULL
head(taxa.print_GTDB)
ref_fasta_species = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/16S_Taxonomy_Databases/GTDB/GTDB_bac120_arc122_ssu_r202_Species.fa.gz"
taxa_species_GTDB = addSpecies(taxa_GTDB,ref_fasta_species, allowMultiple=TRUE)
taxa.print_spp_GTDB  <- taxa_species_GTDB  # Removing sequence rownames for display only
rownames(taxa.print_spp_GTDB) <- NULL
head(taxa.print_spp_GTDB)

#### End ####

####Assign taxonomy for fungi with unite####
#Analaysis of ITS Data Using UNITE.-

#download database: https://unite.ut.ee/repository.php
#ref: https://plutof.ut.ee/#/doi/10.15156/BIO/1281567

unite.ref = "C:/Users/ngs-adm/Desktop/DADA2_Taxonomy_Databases/Fungi_UNITE/sh_general_release_10.05.2021/sh_general_release_dynamic_10.05.2021.fasta"

#tryRC = TRUE
taxa_unite <- assignTaxonomy(seqtab.nochim, unite.ref, tryRC = TRUE)
taxa.print <- taxa_unite  # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)
#### End ####

#### Make Count and Taxa Tables ####

#The typical standard outputs from amplicon processing are a fasta file, a count table, and a taxonomy table. So here's one way we can generate those files from your DADA2 objects in R:
# giving our seq headers more manageable names (ASV_1, ASV_2...)

asv_seqs <- colnames(seqtab.nochim)

asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")


for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making fasta of our final ASV seqs:

asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "asv_fasta.fa")

# OTU count table:

asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "asv_tab.tsv", sep="\t", quote=F, col.names=NA)

# tax table:

asv_tax_unite <- taxa_unite
row.names(asv_tax_unite) <- sub(">", "", asv_headers)
write.table(asv_tax_unite, "asv_tax_unite.tsv", sep="\t", quote=F, col.names=NA)

#### End ####

####add samples metadata and match to count table####

sample_info_tab <- read.table("Sample_Info_crohns.txt",header=T, 
                              row.names = 1, check.names=F)
# rename all table to match phyloseq input

count_tab <- asv_tab
tax_tab <- asv_tax_unite
fasta_tab <- asv_fasta

head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

#examine consistancy in order between cts col names and coldata rownames 
#(they have to be in the same order or Deseq2 won't work out)

all(rownames(sample_info_tab) %in% colnames(count_tab))
all(rownames(sample_info_tab) == colnames(count_tab))


all(rownames(tax_tab) %in% rownames(count_tab))
all(rownames(tax_tab) == rownames(count_tab))

gplots::venn(list(taxonomy=rownames(tax_tab), featuretable=rownames(count_tab)))

#### End ####

#### Rarefaction curves ####
rarecurve(t(count_tab), step=100, lwd=2, ylab="ASVs")
abline(v=(min(rowSums(t(count_tab)))))
#### End ####

####phylogenetic tree####

#phy_tree(phyloseq, errorIfNULL=TRUE)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
mult <- msa(seqs, method="ClustalW", type="dna", order="input")

#The phangorn package is then used to construct a phylogenetic tree. 
#Here we first construct a neighbor-joining tree, and then fit a GTR+G+I maximum likelihood tree using the neighbor-joining tree as a starting point.

phang.align <- as.phyDat(mult, type="dna", names=getSequence(seqtab.nochim))
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

## negative edges length changed to 0!

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

detach("package:phangorn", unload=TRUE)

#### End ####

####making our phyloseq object with transformed table and filter it####
# You may make the physeq object from filtered Count tab if needed 
head(count_tab)
head(tax_tab)
head(fasta_tab)
head(sample_info_tab)

physeq <- phyloseq(otu_table(count_tab, taxa_are_rows = T),   #taxa_are_rows=F (if your taxa names on the column not the rows)
                   sample_data(sample_info_tab), 
                   tax_table(tax_tab))


#Adding ASV Fasta sequences and Phylogenetic tree to phyloseq object

dna <- Biostrings::DNAStringSet(asv_seqs)  #Making ASV Fasta sequences 
names(dna) <- taxa_names(physeq)


ph_tree = phy_tree(fitGTR$tree) #Making Phylogenetic tree
taxa_names(ph_tree) = taxa_names(dna)

physeq_dna_tree <- merge_phyloseq(physeq, ph_tree, dna) #Merging  ASV Fasta sequences and Phylogenetic tree to phyloseq object

taxa_names(physeq_dna_tree) <- paste0("ASV", seq(ntaxa(physeq_dna_tree)))

physeq
physeq_dna_tree
physeq = physeq_dna_tree

#### End ####

## Phyloseq Object Filtering ## First filtering step for low count samples and NA phyla ##

#### remove samples with less than 100 total reads####
sample_sums(physeq) #Nr of reads per Sample
#AA68a  AA68b  AP83a  AP83b  CG86a  CG86b CTH56a CTH56b  KS66a  KS66b  LG00a  LG00b  LR00a  LR00b  SM59a  SM59b 
#7299   3421   5064     22   5147   7457   3302   3126   2024   1218   3494   2495   2028   2130   1513   3427 

physeq_above100reads = prune_samples(sample_sums(physeq)>=100, physeq)

#As sample AP83b was excluded due to low read count, its pair AP83a was alos excluded

newPhyloObject = subset_samples(physeq, Patient_ID != "AP83")
newPhyloObject

sample_sums(newPhyloObject)

physeq = newPhyloObject

#### End ####

#### alpha-diversity (Richness and diversity estimates), "should be done in non filtered physeq"####

sample_variables(physeq)

# Visualize alpha-diversity on unfiltered phyloseq object

P2 = plot_richness(physeq, x="OP_Status", title = "Alpha Diversity", measures=c( "Shannon", "InvSimpson"))

P2 +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))


P3 = plot_richness(physeq, x="OP_Status", color = "OP_Status", title = "Alpha Diversity", measures=c( "Shannon", "InvSimpson"))

P3+
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

P4 = plot_richness(physeq, x="OP_Status", color = "Previous_surgery", title = "Alpha Diversity", measures=c( "Shannon", "InvSimpson"))

P4+
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30)) + geom_point(position = position_dodge(width = 0.75))

# Violin plot

P3 + geom_violin()+
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = 90)) + geom_point(position = position_dodge(width = 0.75))

#### End ####

####richness statistics####
# has to be counts not relative abundances

rich = estimate_richness(physeq, measures = c("Chao1", "Shannon","InvSimpson"))


wilcox.Shannon <- pairwise.wilcox.test(rich$Shannon, 
                                       sample_data(physeq)$OP_Status)

wilcox.Shannon_padj <- pairwise.wilcox.test(rich$Shannon, 
                                            sample_data(physeq)$OP_Status, 
                                            p.adjust.method = "BH")


tab.Shannon <- wilcox.Shannon$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.value", -group1) %>%
  na.omit()


tab.Shannon_padj<- wilcox.Shannon_padj$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.Shannon

#group1  group2     p.adj
# B      A         0.4428332
write.table(tab.Shannon, file = "tab.Shannon.txt")
#
wilcox.InvSimpson <- pairwise.wilcox.test(rich$InvSimpson, 
                                          sample_data(physeq)$OP_Status, 
                                          p.adjust.method = "BH")


tab.InvSimpson <- wilcox.InvSimpson$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.InvSimpson
write.table(tab.InvSimpson, file = "tab.InvSimpson.txt")
#
wilcox.Chao1 <- pairwise.wilcox.test(rich$Chao1, 
                                     sample_data(physeq)$OP_Status, paired = FALSE,
                                     p.adjust.method = "BH")

tab.Chao1 <- wilcox.Chao1$p.value %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "group1") %>%
  gather(key="group2", value="p.adj", -group1) %>%
  na.omit()

tab.Chao1
#group1 group2     p.adj
#1      B      A 0.1655805

#### End ####

####remove NA Phyla####

rank_names(physeq)

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

physeq_oNA <- subset_taxa(physeq, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

physeq = physeq_oNA

table(tax_table(physeq)[, "Phylum"], exclude = NULL)

#### End ####

##Filter low Abundance Taxa and count table normalization##

####Define prevalence of each taxa # (in how many samples did each taxa appear at least once)####
prev0 = apply(X = otu_table(physeq),
              MARGIN = ifelse(taxa_are_rows(physeq), yes = 1, no = 2),
              FUN = function(x){sum(x > 0)})
prevdf = data.frame(Prevalence = prev0,
                    TotalAbundance = taxa_sums(physeq),
                    tax_table(physeq))

#save ASV Prevalence and Abundance table before filtering
write.table(prevdf, "asv_prevdf.tsv", sep="\t", quote=F, col.names=NA)

#Plot Taxa prevalence v. total counts. Each point is a different taxa. 
ggplot(prevdf, aes(TotalAbundance, Prevalence / nsamples(physeq),color=Phylum)) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

#### End ####

####Remove taxa not seen more than 3 times in at least 5% of the samples#### 
#This protects against an OTU with small mean & trivially large C.V.
# Setting filter parameters :

countperphyla = 3
Samplepercentage = 0.05

physeq_filtered = filter_taxa(physeq, function(x) sum(x > countperphyla) > (Samplepercentage*length(x)), TRUE)
physeq = physeq_filtered
#### End ####

####Agglomerate closely related taxa####
# Taxonomic agglomeration
# How many Genera are present after filtering?
taxGlomRank = "Species"
length(get_taxa_unique(physeq, taxonomic.rank = taxGlomRank))
## [1] 104

physeq_glom = tax_glom(physeq, taxrank = taxGlomRank)

# Phylogenetic agglomeration
# How many genera are present after filtering?
length(get_taxa_unique(physeq_glom, taxonomic.rank = taxGlomRank))

# Phylogenetic agglomeration at a fixed tip value
h1 = 0.1
physeq_glom2 = tip_glom(physeq, h = h1)
length(get_taxa_unique(physeq_glom2, taxonomic.rank = "Genus"))

#se phyloseq’s plot_tree() to plot tree of original filtered data, the tree after taxonoic agglomeration, and the tree after phylogenetic agglomeration. Save these as separate plot objects, and then render them together on one plot using gridExtra::grid.arrange.

multiPlotTitleTextSize = 8
p2tree = plot_tree(physeq,
                   method = "treeonly",
                   ladderize = "left",
                   title = "Before Agglomeration") +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p3tree = plot_tree(physeq_glom,
                   method = "treeonly",
                   ladderize = "left",
                   title = paste0("tax_glom(..., taxrank='", taxGlomRank, "')")) +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
p4tree = plot_tree(physeq_glom2,
                   method = "treeonly",
                   ladderize = "left",
                   title = paste0("tip_glom(..., h=", h1, ")")) +
  theme(plot.title = element_text(size = multiPlotTitleTextSize))
# group plots together
grid.arrange(nrow = 1, p2tree, p3tree, p4tree)

#### End ####

#### Normalize number of reads in each sample using median sequencing depth.####

total = median(sample_sums(physeq))
standf = function(x, t=total) round(t * (x / sum(x)))
physeq_mednorm = transform_sample_counts(physeq, standf)

# Transform to relative abundance. Save as new object.
physeq_re = transform_sample_counts(physeq_mednorm, function(x){x / sum(x)})

#### End ####


#### Exploratory plots after filtering and normalization ####
#Check individual phylum Abundance
#Abundance value transformation function

plot_abundance = function(physeq, ylabn = "",
                          Facet = "Class",
                          Color = "Phylum"){
  mphyseq = psmelt(physeq)
  mphyseq <- subset(mphyseq, Abundance > 0)
  ggplot(data = mphyseq,
         mapping = aes_string(x = "OP_Status", y = "Abundance",
                              color = Color, fill = Color)) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3,
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet) + ylab(ylabn) +
    scale_y_log10()
}


#plot the abundance values before and after transformation

pl_ab_original  = plot_abundance(physeq,"Original Abundances")
pl_ab_original_norm  =plot_abundance(physeq_mednorm,"Normalized to squencing depth Abundances")
pl_ab_original_norm_re  =plot_abundance(physeq_re,"Normalized Relative Abundances")

grid.arrange(pl_ab_original, pl_ab_original_norm, pl_ab_original_norm_re)


#Subset by taxonomy (relative abundance)

psra_Lacto = subset_taxa(physeq_re, Order == "Lactobacillales")
plot_abundance(psra_Lacto, Facet = "Genus", Color = NULL)

#### End ####

####Beta diversity:  NMDS, PCoA and PERMANOVA/ADONIS'####
# Done in normalization section after filtering >>> physeq_re

# ordinate with method "NMDS" ON distance "Bray-Curtis" 
ps.prop= physeq_re
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")

#### plot ordination with ggordiplots Package ####

sampdf =sample_data(ps.prop)

# 95% confidence ellipse for centroid

gg_ordiplot(ord.nmds.bray, groups=sampdf$OP_Status, 
            kind='se', conf=0.95)

# spider plot
gg_ordiplot(ord.nmds.bray, groups=sampdf$OP_Status, 
            ellipse=F, spider=T)

# hull plot
gg_ordiplot(ord.nmds.bray, groups=sampdf$OP_Status, 
            ellipse=F, hull=T)

# Customization
ordiplot <- gg_ordiplot(ord.nmds.bray, groups=sampdf$OP_Status, hull = TRUE, spiders = TRUE, 
                       ellipse = FALSE, plot = FALSE)

ordiplot$plot + theme_bw() + labs(color = "OP_Status", x = "NMDS 1", y = "NMDS 2", title = "NMDS ordination of OP_Status samples using Bray–Curtis dissimilarity index") + 
  theme(plot.title = element_text(hjust = 0.5)) 

# Add another variable to plot
ord.data <- ordiplot$df_ord
ord.data$Laparoscopic.open <- sampdf$Laparoscopic.open
colnames(ord.data) <- c("x", "y", "OP_Status", "Laparoscopic.open")
head(ord.data)
ggplot(data = ord.data, aes(x = x, y = y, color = OP_Status, shape = Laparoscopic.open)) + geom_point(size = 3) + 
  xlab("NMDS 1") + ylab("NMDS 2")
#### End ####

#### plot ordination with plot_ordination from Phyloseq ####

plot_ordination(ps.prop, ord.nmds.bray, color="OP_Status", title="Bray NMDS")

plot_ps = plot_ordination(ps.prop, ord.nmds.bray, color="OP_Status", shape ="Laparoscopic.open", title="NMDS ordination of OP_Status samples using Bray–Curtis dissimilarity") 

plot_ps + geom_point(size=8, alpha=0.5) + scale_colour_brewer(type="qual", palette="Set1")

plot_ps + geom_point(size=8, alpha=0.5) + scale_colour_brewer(type="qual", palette="Set1") + stat_ellipse(type = "norm", linetype = 2) + stat_ellipse(type = "t") + theme_bw() #default "t" assumes a multivariate t-distribution while norm (fine line) assumes a multivariate normal distribution

#### End ####

#### Plot only top 5 Abundant Phyla ####

# Keep only the most abundant five phyla to Plot
phylum.sum = tapply(taxa_sums(ps.prop), tax_table(ps.prop)[, "Phylum"], sum, na.rm=TRUE)
top5phyla = names(sort(phylum.sum, TRUE))[1:5]
GP1 = prune_taxa((tax_table(ps.prop)[, "Phylum"] %in% top5phyla), ps.prop)
GP.ord <- ordinate(GP1, method="NMDS", distance="bray")

p1 = plot_ordination(GP1, GP.ord, type="taxa", color="Phylum", title="taxa")
print(p1)
p1 + facet_wrap(~Phylum, 3)

p2 = plot_ordination(GP1, GP.ord, type="biplot", color="Phylum", shape="OP_Status", title="biplot")
print(p2)

p3 = plot_ordination(GP1, GP.ord, type="split", color="Phylum", shape="OP_Status", title="split") 
print(p3)

#### End ####

#### Plot Beta Diversity with All Methods c("DCA", "CCA", "RDA", "CAP", "NMDS", "MDS", "PCoA") ####
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "NMDS", "MDS", "PCoA")
plist = llply(as.list(ord_meths), function(i, physeq, dist){
  ordi = ordinate(physeq, method=i, distance=dist)
  plot_ordination(physeq, ordi, "samples", color="OP_Status")
}, GP1, dist)
names(plist) <- ord_meths

pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"
p5 = ggplot(pdataframe, aes(Axis_1, Axis_2, color=Laparoscopic.open, shape=OP_Status, fill=Laparoscopic.open))
p5 = p5 + geom_point(size=4) #+ geom_polygon()
p5 = p5 + facet_wrap(~method, scales="free")
p5 = p5 + scale_colour_brewer(type="qual", palette="Set1")
print(p5)

#### End ####

####ADONISimplementation of a non-parametric permutation based MANOVA (PERMANOVA)####

dist = phyloseq::distance(ps.prop, method="bray")
metadata <- data.frame(sample_data(ps.prop))

test.adonis <- adonis2(dist ~ OP_Status+Laparoscopic.open, data = metadata, by="margin")
test.adonis
write.table(test.adonis, file = "test.adonis.txt")


#PAIRWISE PERMANOVA to know which of the treatments are significantly different. 

cbn <- combn(x=unique(metadata$Laparoscopic.open), m = 2)
p <- c()

for(i in 1:ncol(cbn)){
  ps.subs <- subset_samples(ps.prop, Laparoscopic.open %in% cbn[,i])
  metadata_sub <- data.frame(sample_data(ps.subs))
  permanova_pairwise <- adonis2(phyloseq::distance(ps.subs, method = "bray") ~ Laparoscopic.open, 
                               data = metadata_sub)
  p <- c(p, permanova_pairwise$aov.tab$`Pr(>F)`[1])
}

p.adj <- p.adjust(p, method = "BH")
p.table <- cbind.data.frame(t(cbn), p=p, p.adj=p.adj)
p.table

#### End ####

#####Bar plots####
# all phyla
#Absolute
plot_bar(physeq, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(~OP_Status, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")

#Relative Abundance
physeq_per <- transform_sample_counts(physeq, function(OTU) OTU/sum(OTU))

plot_bar(physeq_per, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(~OP_Status, scales="free_x") + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")

# Top 20/50 ASVs

top20 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:20]  # adjust number to wished top ones
top50 <- names(sort(taxa_sums(physeq), decreasing=TRUE))[1:50]

ps.top20_per <- prune_taxa(top20, physeq_per)
ps.top20 <- prune_taxa(top20, physeq)

ps.top50_per <- prune_taxa(top50, physeq_per)
ps.top50 <- prune_taxa(top50, physeq)



#relative

plot_bar(ps.top20_per, y =  "Abundance", title = "Abundance based on top 20 ASV", fill="Genus")+ facet_wrap(~OP_Status, scales="free_x") + 
  geom_bar(aes(color=Genus, fill=Genus), stat = "Identity", position = "stack")

plot_bar(ps.top20_per, y =  "Abundance", title = "Abundance based on top 20 ASV", fill="Phylum")+ facet_wrap(~OP_Status, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")


#Absolute
plot_bar(ps.top20, y =  "Abundance", title = "Abundance based on top 20 ASV", fill="Genus")+ facet_wrap(~OP_Status, scales="free_x")

plot_bar(ps.top20, y =  "Abundance", title = "Abundance based on top 20 ASV", fill="Phylum")+ facet_wrap(~OP_Status, scales="free_x") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")


#Bar plot Abundance Per Condition:
# Transform all variables to factors just in case...
df <- as.data.frame(lapply(sample_data(physeq_per),function (y) if(class(y)!="factor" ) as.factor(y) else y),stringsAsFactors=T)
row.names(df) <- sample_names(physeq_per)
sample_data(physeq_per) <- sample_data(df)

ps.OP_Status <- merge_samples(physeq_per, "OP_Status")  # Merging 

ps.OP_Status_per <- transform_sample_counts(ps.OP_Status, function(OTU) OTU/sum(OTU)) #2nd transformation to make it again in percentage

ps.OP_Status_plot = plot_bar(ps.OP_Status_per, y =  "Abundance", title = "Abundance based on Phylum", fill="Phylum")+ facet_wrap(~OP_Status, scales="free_x")

ps.OP_Status_plot + geom_bar(aes(color=Phylum, fill=Phylum), stat = "Identity", position = "stack")


# Plot different conditions together
#relative
plot_bar(physeq_per, "OP_Status", fill="Laparoscopic.open", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack")
plot_bar(physeq_per, "Phylum", fill="Laparoscopic.open", facet_grid=~OP_Status) + geom_bar(stat = "Identity", position = "stack")
plot_bar(physeq_per, "Phylum", fill="Phylum", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))
plot_bar(ps.top20_per, "OP_Status", fill="Laparoscopic.open", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack")  + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))
plot_bar(ps.top20_per, "Phylum", fill="Phylum", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack")
plot_bar(ps.top20_per, "Genus", fill="Genus", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

#absolute
plot_bar(physeq, "Phylum", fill="Laparoscopic.open", facet_grid=~OP_Status) + geom_bar(stat = "Identity", position = "stack")
plot_bar(physeq, "Phylum", fill="Phylum", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))
plot_bar(ps.top20, "OP_Status", fill="Laparoscopic.open", facet_grid=~Phylum) + geom_bar(stat = "Identity", position = "stack")  + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))
plot_bar(ps.top20, "Phylum", fill="Phylum", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack")
plot_bar(ps.top20, "Genus", fill="Genus", facet_grid=OP_Status~Laparoscopic.open) + geom_bar(stat = "Identity", position = "stack") + theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -30))

#### End ####

#### phylogenetic tree Plots ####

plot_tree(physeq)

plot_tree(physeq, nodelabf=nodeplotblank, label.tips="taxa_names", ladderize="left")

plot_tree(physeq, ladderize="left", color="OP_Status", shape="Laparoscopic.open")

plot_tree(physeq, ladderize="left", color="OP_Status", shape="Laparoscopic.open", size="abundance", base.spacing=0.03, min.abundance=3)

plot_tree(physeq, ladderize="left", shape = "OP_Status", color="Phylum", size="abundance", base.spacing=0.03, min.abundance=3)

# for circular shape:      + coord_polar(theta="y")

#Subset Families
physeq.Firmicutes <- subset_taxa(physeq, Phylum=="Firmicutes")

plot_tree(physeq.Firmicutes, title = "Firmicutes-only phylogenetic tree", shape = "OP_Status", color="Genus", label.tips="Genus", size="abundance", plot.margin=0.6)

#### End ####

####Heat Maps####

theme_set(theme_bw())

#Heatmap through phyloseq

rank_names(ps.top50)
rank_names(ps.top50_per)

plot_heatmap(ps.top50_per)
plot_heatmap(ps.top50)

plot_heatmap(ps.top50, sample.label = "Protocol", taxa.label = "Genus", low = "#000033",
             high = "#FF3300", na.value = "black",
             max.label = 250, title = "Heatmap of top 50 Genus (GTDB)", sample.order = NULL, taxa.order = "Genus",
             first.sample = NULL, first.taxa = NULL)

plot_heatmap(ps.top50, "NMDS", "bray", "OP_Status", "Family")


#Heatmap through otu table 

otu_table_top50 = otu_table(ps.top50)

#Function to turn counts to Z-scores

cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

otu_table_top50_norm <- t(apply(otu_table_top50, 1, cal_z_score))
pheatmap(otu_table_top50_norm)

pheatmap(otu_table_top50_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = sample_info_tab[c("OP_Status", "Laparoscopic.open")]) 

#If you want to change the colors of the OP_Status and Laparoscopic.open annotation

colors = viridis::viridis(4, begin = 0.2, end = 0.9)    #Or introduce here any other 4 colors
names(colors) = c("A", "Culture-Enriched", "B", "Direct-Isolation")
colors = list(OP_Status = colors[c(1,3)], Laparoscopic.open = colors[c(2,4)])

pheatmap(otu_table_top50_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = sample_info_tab[c("OP_Status", "Laparoscopic.open")],
         annotation_colors = colors) 

#Just an idea for clustering, I don´t know what it is better/you prefer for microbiome

row.dist <- vegdist(otu_table_top50, method = "bray")
row.clus <- hclust(row.dist, "ward.D2") #or with "aver"

col.dist <- vegdist(t(otu_table_top50), method = "bray")
col.clus <- hclust(col.dist, "ward.D2") #or with "aver"

pheatmap(otu_table_top50_norm, fontsize_number = 2 * fontsize,
         cellwidth = 15, cellheight = 12, scale = "none",
         treeheight_row = 100,
         kmeans_k = NA,
         show_rownames = T, show_colnames = T,
         main = "Full heatmap (avg, eucl, unsc)",
         clustering_method = "average",
         cluster_rows = row.clus, cluster_cols = col.clus,
         annotation_col = sample_info_tab[c("OP_Status", "Laparoscopic.open")],
         annotation_colors = colors) 


#### End ####

####Analysis with DESEq2####

#Filter ASVs with low read counts across all samples#

count_tab_keep = rowSums (count_tab) > 10
count_tab_filtered = count_tab[keep,] 

#Make the Deseq Tables

deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~ OP_Status) 


#prefiltering low count ASVs (1664 ASVs)
#filter out genes with less than 10 counts in total

keep = rowSums(counts(deseq_counts)) > 10

deseq_counts_filtered = deseq_counts[keep,] 
deseq_counts_filtered #rownames/ASVs (1424)  

#### End ####

#### session documentation ####

writeLines(capture.output(sessionInfo()), "sessionInfo_11Nov2022.txt")

#### End ####




