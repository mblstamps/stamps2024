######################################################################
#####  GOAL: Investigate the low biomass bile environment with   #####
#####         the help of the decontam R package.                #####
######################################################################
#####   DATA is from V4 16S sequencing of bile samples from cats
#####    with and without suspected hepatobiliary disease, and
#####    a number of negative controls, as described in the manuscript
#####    "Microbiome Analysis of Bile from Healthy Cats and Cats with 
#####     Suspected Hepatobiliary Disease". Slead, Callahan, et al.
#####     Journal of Veterinary Internal Medicine, 2023.
#####
######################################################################

library(decontam); packageVersion("decontam")
library(ggplot2); packageVersion("ggplot2")

# First download the data being used for this lab
download.file("https://figshare.com/ndownloader/files/41617905",
              "bile_lab.rda")

# Where is the freshly downloaded data?
list.files()
setwd("PATH")
# REPLACE PATH with the path to the folder containing the .rda object
# Challenge: Define the path both as a relative path and as an absolute path.

# Load the pre-processed 16S sequencing data (DADA2)
#  as well as the sample metadata
load("bile_lab.rda")
# In Rstudio, four newly defined variables should appear in your environment.
# st: The ASV table produced by dada2. Rows are samples, columns are ASVs, values are counts.
# ft: The table of ASV proportions. A simple transformation of st. Values are proportions [0 to 1].
# tax: The taxonomic assignments of each ASV by `assignTaxonomy` using the Silva database.
# df: The sample metadata, with columns...
#     $SampleID: Sample ID
#     $Sample_Type: Bile -- samples from feline bile. Negative -- samples from full-process negative controls.
#     $Bactibilia: The presence of bacteria in bile as ascertained by microscopy. TRUE/FALSE (if measured)
#     $Ecoli_Culture: The presence of E. coli in the sample as ascertained by culture. TRUE/FALSE (if measured)
#     $Library_Concentration: The DNA concentration after library preparation. Non-negative numeric.
#     $Library_Size: The number of reads obtained from each sample's sequencing library. Non-negative integer.
#     $Shannon_Diversity: The Shannon Diversity Index obtained from the sequenced sample. Non-negative numeric.
#
df
# Challenge: Use R functions like `summary` and `table` to conform to the descriptions above.

# The decontam package can be used to discriminate between contaminants, which could be likely given the very low 
#  biomass samples we are working with here (bile is traditionally considered "sterile").
# The decontam package implements two different methods:
#     decontam-prevalence uses negative controls to classify contaminants/non-contaminants.
#     decontam-frequency uses DNA concentration information to classify contaminant/non-contaminants.
# Given the number of negative controls available (how many ar there?) we'll apply decontam-prevalence
# First we need to define a TRUE/FALSE vector with TRUE entries for negative controls (see ?isContaminant)
is.control <- FILL_THIS_IN
# REPLACE FILL_THIS_IN with an R command that provides the required logical (TRUE/FALSE) vector
contam <- isContaminant(st, neg=is.control, method="prevalence")
# Explore the `contam` data.frame. What information does it contain?
# What does each row of the `contam` data.frame correspond to?

# We will now augment the contam data.frame with some additional annotations
# This uses some base R manipulation -- if you are not R familiar, that's OK.
# Make a map between short ASV names and the full sequences
sq <- colnames(st)
names(sq) <- paste0("ASV", seq_along(sq))
# Sanity check that our sequence table and the contam data.frame are arranged in the same order
if(!identical(colnames(st), rownames(contam))) stop("Mismatch between st and contam.")
# Add new columns to the contam data.frame with taxonomic information
contam$Phylum <- tax[,"Phylum"]
contam$Genus <- tax[,"Genus"]
contam$ASV <- names(sq)
rownames(contam) <- NULL # To avoid printing long sequences
# Add a column with binned prevalances
contam$Prevalence_Binned <- cut(contam$prev, c(0, 10, 25, 45, 9999), labels=c("1-10", "11-25", "26-45", "45+"))
# View
head(contam)
# Challenge: What does the `cut` command do? How might you create a custom binned prevalence column?

# The decontam manuscript emphasizes inspection of the score assigned by decontam ($p)
#   to identify an appropriate classification threshold. Ideally, there will be an evident
#   high-score (non-contaminant) mode, and a low-score (contaminant) mode.
# Inspect the histrogram of the scores to look for this bimodality, and help choose a classification threshold.
histo <- ggplot(data=contam, aes(x=p)) + 
  labs(x = 'decontam-prevalence Score', y='Number of ASVs') + 
  geom_histogram(binwidth=0.02)

histo
# Is this consistent with decontam's hypothesized bimodal distribution?
# What classification threshold does this suggest?
# Challenge: What is causing the warning in the previous command?

# In the decontam manuscript (https://doi.org/10.1186/s40168-018-0605-2) the distribution
#    of scores was also investigated after stratification by prevalence, i.e. taxa present
#    in more or fewer samples.
histo + facet_wrap(~Prevalence_Binned)
# What classification threshold does this stratified analysis suggest?

# Inspecting sequences at the extremes of the score distribution can help guide
#   choices about contaminant classification.
# Let's look more closely at the 10 highest scores (the 10 ASVs that look the most like non-contaminants according to decontam)
i.top10 <- order(contam$p, decreasing=TRUE)[1:10]
contam[i.top10,]
# Print out those ASV sequences in a fasta format
dada2:::pfasta(sq[contam$ASV[i.top10]], id=contam$Phylum[i.top10])
# Use the NCBI web BLAST interface (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch) 
# to BLAST those sequences against Standard databases (nr etc.) -- nr/nt
# What did we learn from these BLAST results?

# Now look at the the 10 lowest scores (the 10 ASVs that look the most like contaminants according to decontam)
i.bottom10 <- order(contam$p, decreasing=FALSE)[1:10]
dada2:::pfasta(sq[contam$ASV[i.bottom10]], id=contam$Phylum[i.bottom10])
# BLAST against nr/nt as above...

###################################################################
######  How do we make sense of these results so far?        ######
######  What are more or less valid interpretations of this  ######
######    data given our findings so far?                    ######
###################################################################

# Look more closely at the differences among bile samples.
# We'll plot a standard alpha-diversity metric (Shannon) across samples types.
ggplot(data=df, aes(x=Sample_Type, y=Shannon_Diversity)) + geom_jitter(width=0.1)
# How do the distributions of Shannon diversity compare between bile and negative controls?

# Let's investigate how this relates to bacterial infection as measured by
#   Ecoli culture and by microscopic identification of Bactibilia
ggplot(data=df,
       aes(x=Sample_Type, y=Shannon_Diversity, color=XXX)) + # Fill in XXX
  geom_jitter(width=0.1)
# How does 16S-measured diversity line up with culture and microscopy in these samples?
# Challenge: How do culture and microscopy results compare with one another?

# Now we'll leverage phyloseq to perform a couple of visualizations.
library(phyloseq); packageVersion("phyloseq")
# Consider an ordination of these samples ("beta-diversity" analysis)
ps <- phyloseq(sample_data(df), 
               otu_table(ft, taxa_are_rows=FALSE), 
               tax_table(tax))
ord <- ordinate(ps, method="MDS", distance="bray")
f <- plot_ordination(ps, ord, color="Ecoli_Culture") + facet_wrap(~Sample_Type)
f
# What does this suggest about the nature of Ecoli positive bile samples?
# How about Bactibilia positive samples?

# The low-diversity, Ecoli positive, Bactibilia positive samples are of interest.
#   What is the spectrum of bacteria detected by 16S sequencing in those samples?
ps.infected <- subset_samples(ps, Ecoli_Culture %in% "Positive")
ps.infected <- subset_taxa(ps.infected, 
                           taxa_sums(ps.infected) > 0.01 * nsamples(ps.infected))
plot_bar(ps.infected, fill="Genus") + ylab("Proportion")
# What does this tell us about putative bile infections?
# How could this inform treatment?

