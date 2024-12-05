# Project information -----------------------------------------------------

# Assignment 5 — Ulli Bodnar

# Original sources for meat and potatoes of my script (partial following of vignettes)
# Tutorial by: Benjamin J Callahan, Joey McMurdie, Susan Holmes, October 24, 2023
# https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

# https://benjjneb.github.io/dada2/tutorial.html

# The source data for this project has been demultiplexed and has had the primer sequences removed. This is needed so that the analyses does not confuse primer sequences for actual data.



# Load packages -----------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("dada2", "phyloseq"))

# Package loader adapted from Moiz—thank you!
packages <- c("tidyverse", "viridis", "dada2", "phyloseq", "Biostrings", "ggplot2", "vegan")

suppressMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

conflicted::conflicts_prefer(dplyr::filter())



# Global variables --------------------------------------------------------

# Path to data
path <- "../data/MiSeq_SOP"
list.files(path) # Check that they are all there :')


# This variable will store a series of parameters to test so that I don't have to hard code each one
parameters <- expand.grid(
  truncLen = list(c(260, 220), c(240, 200), c(220, 180), c(200, 160), c(180, 140)), # I have chosen these truncation lengths by choosing extreme values on both ends. If the truncation is overly aggressive, it can discard valid sequcnes, which would provide improper results. 
  maxEE = c(2, 4), # This limits the amount of expected errors per read, 2 is the typical parameter, but 4 is considered a more relaxed filtering step, allowing for more sequences. However, it risks introducing more errors. 
  trimLeft = c(10, 20) # This determines the number of nucleotides to trim from the start of each read. I have chosen these values because 10 is typically the default, but a more agressive level, i.e., 20, can ensure that artifacts and low quality bases are removed. 
)
checkObjects(parameters) # Need to load functions FIRST for this to run, I kept it like this because it's cleaner to have the variables before the functions

sequenceTables <- data.frame(length(parameters))



# Data acquisition --------------------------------------------------------

# Pull in the data as matched lists
forwardList <- path |>
  list.files(pattern = "_R1_001.fastq", full.names = TRUE) |>
  sort()

reverseList <- path |> 
  list.files(pattern = "_R2_001.fastq", full.names = TRUE) |>
  sort()

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(forwardList), "_"), `[`, 1)


# Global Functions --------------------------------------------------------

# Check object function adapted from assignment 3—when I edited Moiz's script I found that this function was a blessing
# Create a function that can clean the parameters before being passed to the calculateDiversity function. When they are extracted, they are passed as a list, but I want them to be passed as.numeric so that they are in the expected format.
unlistNumeric <- function (x) {
  return (
    as.numeric(unlist(x))
  )
}


# Create a check object function that can be reused to check code and objects during programming so that I don't have to keep running these same checks each time.
checkObjects <- function (...) {
  
  checkObject <- function (obj) {
    results <- list()
    
    # Provide information relevant if object being checked is a data frame
    if(is.data.frame(obj) || is.matrix(obj)) {
      results$colNames <- colnames(obj)
      results$summary <- tryCatch(
        {
          summary(obj)  # Attempt to get the summary
        },
        error = function(e) {
          paste("Error in summary:", e$message)  # Capture the error message
        }
      )
      results$unique_values <- sapply(obj, function(col) length(unique(col)))
      results$dims <- dim(obj)
      results$duplicate_rows <- sum(duplicated(obj))
    }
    
    # Provide simple checks for objects that do not fall under data frame or matrix class
    results$length <- length(obj)
    results$na_count <- sum(is.na(obj))
    results$head <- head(obj, 2)
    results$class <- class(obj)
    
    results
  }
  
  
  map(list(...), checkObject)
  
}


calculateDiversity <- function (
    fnFs, 
    fnRs, 
    debug = FALSE,
    truncLen = c(240, 200), 
    maxEE = 2, 
    trimLeft = 10
  ) {
  
  truncLen <- as.numeric(truncLen)
  
  print(class(truncLen))
  print(class(maxEE))
  print(class(trimLeft))
  # # Get the filenames of our example paired-end fastq files—these come included with the package so no need to obtain new
  # fnF1 <- system.file("extdata", file1, package = "dada2")
  # fnR1 <- system.file("extdata", file2, package = "dada2")
  
  if (debug == TRUE) {
    print("debugging")
    
    # setting up unique, temporary files, buried in the systems temporary directory
    filtFs <- file.path(path, "filtered", paste0(sample.names[18:20], "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names[18:20], "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names[18:20]
    names(filtRs) <- sample.names[18:20]
    
  }
  else {
    # setting up unique, temporary files, buried in the systems temporary directory
    filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
    filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
    names(filtFs) <- sample.names
    names(filtRs) <- sample.names
  }
  
  # Check the location
  checkObjects(filtFs, filtRs)
  
 
  # When inspecting the quality of the forward and reverse sequences, it is clear that the reverse drop off in quality much faster than the forward reads, therefore we will trim the forward reads at 240, and the reverse at 200. It is also recommended to trim the initial reads as well since calibration issues can cause them to be problematic. The ambiguous nucleotides are filtered out with maxN = 0, and the reads with more than 2 expected errors are filtered out as well. 
  filterAndTrim(fwd = fnFs, filt = filtFs, rev = fnRs, filt.rev = filtRs,
                trimLeft = trimLeft, 
                truncLen = truncLen, 
                maxN = 0, # DADA2 requires no Ns
                maxEE = maxEE,
                compress = TRUE, verbose = TRUE)
  
  
  # To collapse reads that all encode the same sequence, we use dereplication. This reduces computational load and therefore, time.  
  derepFs <- derepFastq(filtFs, verbose = TRUE)
  derepRs <- derepFastq(filtRs, verbose = TRUE)
  
  
  # PCR amplification sometimes introduces errors, which vary between sequencing runs, so this is a way to estimate the parameters from the data
  errFs <- learnErrors(derepFs, multithread = FALSE)
  errRs <- learnErrors(derepRs, multithread = FALSE)
  
  
  # This infers the sample composition
  dadaFs <- dada(derepFs, err = errFs, multithread = FALSE)
  dadaRs <- dada(derepRs, err = errRs, multithread = FALSE)
  print(dadaFs)
  print(dadaRs)
  
  
  # Now, we merge the forward and reverse reads since above we inferred the sample sequences independently. This returns a data frame
  merger1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)
  checkObjects(merger1)
  
  
  # Chimeras are errors in PCR amplification and should be removed, so the following removes them. Chimeras are essentially a sequence read that is half of one sample sequence and half of another.
  merger1.nochim <- removeBimeraDenovo(merger1, multithread = FALSE, verbose = TRUE)

  
  # We can now make a unified sequence table, with the combined, inferred samples. 
  seqtab <- makeSequenceTable(merger1.nochim)
  seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)
  
  
  # We can now calculate our richness and shannon index
  richness <- rowSums(seqtab.nochim > 0)
  shannon <- diversity(seqtab.nochim, index = "Shannon")
  
  return(
    list(richness = richness, shannon = shannon)
  )
  
  
}

  
  
# Create a merger ---------------------------------------------------------


# From the help for the plotQualityProfile function: "The plotted lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles."
#plot quality of forward sequences
plotQualityProfile(forwardList[1:2])
#plot quality of reverse sequences
plotQualityProfile(reverseList[1:2])


# CAN DELETE AT SOME POINT... Data organization changed
# merger1 <- makeMerger("sam1F.fastq.gz", "sam1R.fastq.gz")
# merger2 <- makeMerger("sam2F.fastq.gz", "sam2R.fastq.gz")



# Analysis operations -------------------------------------------------------

# This is where the data and the parameter combinations will be passed to the calculateDiversity function to receive a list of metrics that can be further analyzed


results <- apply(parameters[2, ], 1, function (row) {
  calculateDiversity(
    forwardList[18:20], 
    reverseList[18:20],
    truncLen = unlistNumeric(row["truncLen"]),
    maxEE = unlistNumeric(row["maxEE"]),
    trimLeft = unlistNumeric(row["trimLeft"]),
    debug = TRUE)
})





# Assign taxonomy ---------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, "../data/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)


# Phyloseq -------------------------------------------------------





#From the tutorial: "This is the final product of the dada2 pipeline, a matrix in which each row corresponds to a processed sample, and each column corresponds to an non-chimeric inferred sample sequence (a more precise analogue to the common "OTU table")." And, each cell contains the count of sequences for that variant. "From here we recommend proceeding forward with our friend the phyloseq package for further analysis."


# POSSIBLY DO THIS:
#Build a phylogenetic tree from the unique sequences and mark the relative abundances for the variants in each of the two samples.
#Example package for ideas the visualization: https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

