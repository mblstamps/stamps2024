######################################################################
#####  GOAL: Choose and validate appropriate dada2 parameters    #####
#####         for processing full-length PacBio HiFi 16S data.   #####
######################################################################
#####   DATA is from full-length 16S sequencing of 6 human fecal samples
#####    and one mock community sample using PacBio HiFI sequencing that
#####    were generated for this paper: https://doi.org/10.1093/nar/gkz569
#####   The human fecal samples come from three time points (T1,T2,T3)
#####    and two subjects (S3 and S9), and the mock community is the
#####    ZymoBIOMICS Microbial Community DNA Standard.
######################################################################

library(dada2); packageVersion("dada2") # 1.16 or later

# First download the data being used for this lab
download.file("https://figshare.com/ndownloader/files/41575518", "data_PacbioHiFi.zip")
unzip("data_PacbioHiFi.zip", exdir="data_PacbioHiFi")

# Where is the freshly downloaded data?
list.files()
path <- "PATH" # REPLACE PATH with the path to the unzipped directory containing the gzipped fastqs
# Challenge: Define the path both as a relative path and as an absolute path.

# Read in forward and reverse fastq filenames
fns <- list.files(path, pattern=".fastq.gz", full.names=TRUE)
fns
# Do those file names look like what we expect?

plotQualityProfile(fns)
# Good quality/bad quality? What are the read lengths?
# What are the expected read lengths for the full-length 16S gene?

# Define the forward and reverse primers that were used to amplify the full-length 16S gene
FWD <- "AGRGTTYGATYMTGGCTCAG" # The 27F primer
REV <- rc("RGYTACCTTGTTACGACTT")  # The 1492R primer, in reverse-complement orientation
# About how long should the amplified sequences be?

# In typical PacBio HiFi data (and here) the primers are sequenced, and the
# reads are in mixed orientation (some forward, some reverse complement).
# The `removePrimers` function will trim the reads to the sequence between
# the primers, and orient the reads consistently in the forward direction.
# Reads in which both primers were not detected are discarded.
# As always, use R's built-in help `?removePrimers` for more information.

# Define primer-removed and oriented filenames
nops <- file.path(path, "noprim", basename(fns))
prim <- removePrimers(fns, nops, 
                      primer.fwd=FWD, primer.rev=REV, 
                      orient=TRUE)
prim
# The primer-removed/oriented files will be in the `noprim/` subdirectory within `path`
# Did most reads pass primer removal?

plotQualityProfile(nops)
# How did the distribution of read lengths and qualities change after `removePrimers`?
# What might make sense as `minLen` and `maxLen` filtering parameters?

# Define filtered file names
filts <- file.path(path, "filtered", basename(fns))
# The filtered files will be in the `filtered/` subdirectory within `path`

# Perform filtering and trimming
out <- filterAndTrim(nops, filts,
                     maxEE=1, minLen=XXX, maxLen=YYY, # REPLACE XXX/YYY with proper parameter choices
                     multi=TRUE)
out
# Were most reads retained during filtering? If not, would any parameter tweaks help?
# Note: These samples were all subsampled (or rarefied) down to 12,500 reads to reduce computation time.

plotQualityProfile(filts)
# Does this look about right for the length distribution of full-length 16S?
# Note, see `?plotQualityProfile` for a better understanding of these plots, including the red line.

######################################################################
#####  IMPORTANT SANITY CHECK FOR LONG-READ AMPLICON DATA        #####
#####  DADA2 assumes that the sequencing is accurate enough      #####
#####   that a large fraction of reads will have zero errors.    #####
#####  One signature of this is that there should be many        #####
#####   sequences seen repeatedly in the data, i.e. error-free   #####
#####   reads from sufficiently abundant taxa.                   #####
#####  If this assumption is not met, DADA2 is not appropriate   #####
#####   and another approach should be pursued.                  #####
######################################################################

# Dereplication sanity check
derepFastq(filts[[1]], verbose=TRUE)
# Is there evidence of a large number of duplicated reads in the data?

# Learn the error model from the filtered data.
# Note that a PacBio-specific error-modeling function is used here.
# This function takes a couple minutes to run, take this time to navigate
#   to https://github.com/benjjneb/LRASManuscript and check out the PacBio-specific
#   DADA2 workflows available there.
err <- learnErrors(filts, 
                   errorEstimationFunction=PacBioErrfun, 
                   BAND_SIZE=32, # Changes the band-size in the Needleman-Wunsch alignment
                   multithread=TRUE)
plotErrors(err)
# Does the fitted error model look reasonable?

# The `learnErrors` step can be time-consuming. It can be useful to save the result.
dir.create(file.path(path, "RDS"))
saveRDS(err, file.path(path, "RDS", "err.rds"))
# Later, we can reload it, and skip the `learnErrors` step, with:
#  err <- readRDS(file.path(path, "RDS", "err.rds"))
# This is a general purpose way to save and reload any object in R.
# WARNING: If functional changes were made prior to this stage, this step needs
#   to be re-run. Final results should always be regenerated from a complete re-run
#   of the workflow.

# Run the DADA2 method using the fitted error model.
# Note, this step will again take a couple minutes, similar to `learnErrors`.
dd <- dada(filts, err,
           pool=XXX, # REPLACE XXX with an appropriate pooling option
           BAND_SIZE=32, multithread=TRUE)
# What pooling option makes sense? FALSE (default), "pseudo", or TRUE? See ?dada for more info
# For more pseudo-pooling detail: https://benjjneb.github.io/dada2/pseudo.html#pseudo-pooling
# Challenge: Try different pooling options and compare them.

# Construct a sequence table: rows are samples, columns are ASVs, values are abundances.
# ASVs (columns) are ordered by their summed abundance across all samples.
sta <- makeSequenceTable(dd)
dim(sta)
# How many samples and ASVs are in this table?

# Remove chimeric ASVs and construct a new chimera-free sequence table.
# Note the larger minFoldParentOverAbundance parameter that is recommended
#   for long read data because of intra-genomic heterogeneity in the multi-copy 16S gene.
st <- removeBimeraDenovo(sta,
                         minFoldParentOverAbundance=4,
                         multi=TRUE, verbose=TRUE)
dim(st)
sum(st)/sum(sta)
# Were most reads retained during chimera removal? How about ASVs? 
# How does the fraction of chimeric ASVs compare to the fraction of chimeric reads?
# Why is that?

####################################################################
######  Inspect the number of reads passing through each step ######
######  THIS IS THE SINGLE MOST IMPORTANT SANITY CHECK!!      ######
####################################################################
# Code modified from the dada2 tutorial: https://benjjneb.github.io/dada2/tutorial.html
track <- cbind(prim, out[,2], rowSums(sta), rowSums(st))
colnames(track) <- c("input", "primer-removed", "filtered", "denoised", "nonchim")
rownames(track) <- basename(fns)
head(track)
# In this case, most reads should make it through the entire pipeline!
# IF THAT ISN'T THE CASE, you have a problem, and need to revisit your filtering
#   or primer removal parameters.

# The dada2 taxonomic reference page https://benjjneb.github.io/dada2/training.html has links to a 
#   number of reference databases formatted to work with `assignTaxonomy`. I recommend Silva for 16S.
#   Full-length 16S is suitable for species-level assignment using the naive Bayesian classifier
#   method (i.e. `assignTaxonomy`), so we will use the silva_nr99_v138.1_wSpecies_train_set.fa.gz
#   training file which includes species-level taxonomy.
# Note, this file is not small (132MB), and may take a bit to download. Thus we need
#   to increase R's default timeout parameter of 60 seconds.
options(timeout=300) # Increase if download is still timing out
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1",
              "silva_nr99_v138.1_wSpecies_train_set.fa.gz") 
# Assign taxonomy down to the genus level to these 16S ASVs.
# This may take a few minutes of time. Maybe the DADA2 developers should improve 
#   the performance of `assignTaxonomy`!
system.time(tax <- assignTaxonomy(st, 
                                  "silva_nr99_v138.1_wSpecies_train_set.fa.gz", 
                                  multi=TRUE))
unname(head(tax))
# Are reasonable taxonomic assignments being made to the abundant taxa?
# What does NA mean here?
# Challenge: BLAST a few of these sequences against nt at NCBI: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&BLAST_SPEC=GeoBlast&PAGE_TYPE=BlastSearch
#   Do the `assignTaxonomy` results agree with BLAST against nt?
# Hint: `getSequences(tax)` will return a vector of ASVs. `getSequences` can also be used on
#   many other objects created by dada2, e.g. `st` or `dd[[1]]` or even a fastq filename.

# We'll focus in on the mock community sample, "zymo_sub.fastq.gz"
# First, how many ASVs are present in that sample?
in.mock <- st["zymo_sub.fastq.gz",]>0 # this creates logical vector, i.e. TRUE/FALSE, if the ASV is present in this sample
sum(in.mock)
# How many ASVs are present in this sample?
# Why might this differ from the 8 bacterial strains present in the
#   ZymoBIOMICS Microbial Community DNA Standard?

# Inspecting the taxonomies assigned to the mock ASVs
unname(tax[in.mock,]) # unname strips the long ASV rownames to make the printed output more readable
table(tax[in.mock,"Genus"])
# How do these results help explain the mock community ASVs?
# Challenge: Plot the relative abundances of the ASVs for each genus in the zymo mock sample

