######################################################################
#####  GOAL: Choose and validate appropriate dada2 parameters    #####
#####         for processing paired-end Ilummina 16S data        #####
######################################################################
#####   DATA is taken from the following paper                   #####
#
# Pyrethroid exposure alters internal and cuticle surface bacterial 
#    communities in Anopheles albimanus, ISMEJ, 2019.
# https://doi.org/10.1038/s41396-019-0445-5
# Sequencing data: First 18 samples from https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA512122
#
######################################################################

library(dada2); packageVersion("dada2") # 1.16 or later

# First download the data being used for this lab
download.file("https://figshare.com/ndownloader/files/41575305", "data_IlluminaPE.zip")
unzip("data_IlluminaPE.zip", exdir="data_IlluminaPE")

# Where is the freshly downloaded data?
list.files()
path <- "PATH" # REPLACE PATH with the path to the unzipped directory containing the gzipped fastqs
# Challenge: Define the path both as a relative path and as an absolute path.

# Read in the forward and reverse fastq file names
fnFs <- list.files(path, pattern="_1.fastq.gz", full.names=TRUE)
fnRs <- list.files(path, pattern="_2.fastq.gz", full.names=TRUE)
head(fnFs)
# Do those file names look like what we expect?
# How would you change the pattern to match other file name structure?

plotQualityProfile(fnFs)
plotQualityProfile(fnRs)
# Good quality/bad quality? What are the read lengths?
# How might these quality profiles inform our choice of truncation lengths?

# Define the paths to filtered files we are going to create
filtFs <- file.path(path, "filtered", basename(fnFs)) 
filtRs <- file.path(path, "filtered", basename(fnRs))
# The filtered files will be in the `filtered/` sub-directory within `path`

###################################################################
######  Are primers on these reads that need to be removed?  ######
######  How long is the sequenced amplicon?                  ######
######  SEE PAPER, INSPECT READS.                            ######
###################################################################

# Perform filtering and trimming
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxEE=2, 
                     trimLeft=c(XXX, YYY), # REPLACE XXX/YYY with valid parameter choices
                     truncLen=c(XXX, YYY)) # REPLACE XXX/YYY with valid parameter choices
# Were most reads retained during filtering? If not, would any parameter changes help?
# How might the depth of sequencing in this data affect the questions that can be addressed?

# Learn the error model from the filtered data.
errF <- learnErrors(filtFs, multi=TRUE) # `multi=TRUE` activates multi-threading to reduce computation time
errR <- learnErrors(filtRs, multi=TRUE)

# Visualize the error model. Points are observations, black line is fitted error model.`
plotErrors(errF)
plotErrors(errR)
# Do the fitted error models look reasonable?

# Run the DADA2 method using the fitted error model.
ddF <- dada(filtFs, errF, pool="pseudo", multi=TRUE)
ddR <- dada(filtRs, errR, pool="pseudo", multi=TRUE)
# What pooling option makes sense? FALSE (default), "pseudo", or TRUE? See ?dada for more info
# For more pseudo-pooling detail: https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling
# Challenge: Try different pooling options and compare them.

# Merge the denoised forward and reverse reads together.
mm <- mergePairs(ddF, filtFs, ddR, filtRs, verbose=TRUE)
# Were most reads retained during merging. If not, why not?

# Construct a sequence table: rows are samples, columns are ASVs, values are abundances.
# ASVs (columns) are ordered by their summed abundance across all samples.
sta <- makeSequenceTable(mm)
dim(sta)
# How many samples and ASVs are in this table?
# Challenge: How would you make a sequence table from forward reads alone?

# Remove chimeric ASVs and construct a new chimera-free sequence table.
st <- removeBimeraDenovo(sta, multi=TRUE, verbose=TRUE)
sum(st)/sum(sta)
# Were most reads retained during chimera removal? How about ASVs? 
# How does the fraction of chimeric ASVs compare to the fraction of chimeric reads?
# Why is that?

####################################################################
######  Inspect the number of reads passing through each step ######
######  THIS IS THE SINGLE MOST IMPORTANT SANITY CHECK!!      ######
####################################################################

# Code derived from the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(ddF, getN), sapply(ddR, getN), sapply(mm, getN), rowSums(st))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- basename(fnFs)
head(track)
# In this case, most reads should make it through the entire pipeline!
# Most importantly, a large majority (>80% of reads) should merge successfully,
#  and almost all (>95%)  reads should pass chimera filtering.
# IF THAT ISN'T THE CASE, you have a problem, and need to revisit your truncation lengths
#   (merging problem) or primer removal (trimLeft, chimera problem).

# The dada2 taxonomic reference page https://benjjneb.github.io/dada2/training.html has links to a 
#   number of reference databases formatted to work with `assignTaxonomy`.I recommend Silva for 16S.
#   However, here we will use RDP in this lab for file size reasons.
download.file("https://zenodo.org/record/4310151/files/rdp_train_set_18.fa.gz?download=1",
              "rdp_train_set_18.fa.gz")
# Assign taxonomy down to the genus level to these 16S ASVs.
tax <- assignTaxonomy(st, "rdp_train_set_18.fa.gz", multi=TRUE)
unname(head(tax)) # unname strips the long ASV rownames to make the printed output more readable
# Are reasonable taxonomic assignments being made to the abundant taxa?
# Challenge: Do this again using the Silva database, and compare results.
# Challenge 2: Use `addSpecies` or `assignSpecies` to do species-level assignemnt, where appropriate.

# Phyloseq is a package for the manipulation and analysis of microbiome data.
# Here we use it briefly to produce an ordination of our sequenced communities.
library(phyloseq); library(ggplot2)

# We define a very simple data.frame that records the 3 experimental groups these samples 
#   came from (see paper for more info)
samdf <- data.frame(row.names=rownames(st), 
                    sampleID=rownames(st), 
                    treatment=rep(c("Unexposed", "Susceptible", "Resistant"), each=6))

# Use phyloseq to plot a bray-curtis NMDS odrination of our samples, colored by treatment.
ps <- phyloseq(sample_data(samdf), otu_table(st, taxa_are_rows=FALSE), tax_table(tax))
plot_ordination(ps, ordinate(ps, method="NMDS", distance="bray"), color="treatment") + 
  aes(size=4) + theme_bw() + guides(size="none")
# What do you see?
# Challenge: Does this correspond to what was reported in the paper?
# Note that this is a small subset of the total data in the full study...
