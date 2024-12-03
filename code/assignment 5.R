# Project information -----------------------------------------------------

# Assignment 5 — Ulli Bodnar

# Original source for meat and potatoes of the script (partial following of vignette)
# Tutorial by: Benjamin J Callahan, Joey McMurdie, Susan Holmes, October 24, 2023
# https://bioconductor.org/packages/devel/bioc/vignettes/dada2/inst/doc/dada2-intro.html

# The source data for this project has been demultiplexed and has had the primer sequences removed. This is needed so that the analyses does not confuse primer sequences for actual data.



# Load packages -----------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("dada2", "phyloseq"))

# Package loader adapted from Moiz—thank you!
packages <- c("tidyverse", "viridis", "dada2", "phyloseq")

suppressMessages(
  invisible(lapply(packages, library, character.only = TRUE))
)

conflicted::conflicts_prefer(dplyr::filter())



# Global Functions --------------------------------------------------------

# Check object function adapted from assignment 3—when I edited Moiz's script I found that this function was a blessing
# Create a check object function that can be reused to check code and objects during programming so that I don't have to keep running these same checks each time.
checkObjects <- function (...) {
  
  checkObject <- function (obj) {
    
    results <- list()
    
    
    # Provide information relevant if object being checked is a data frame
    if(is.data.frame(obj) || is.matrix(obj)) {
      results$colNames <- colnames(obj)
      results$duplicate_rows <- sum(duplicated(obj))
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


makeMerger <- function (file1, file2) {
  # Get the filenames of our example paired-end fastq files—these come included with the package so no need to obtain new
  fnF1 <- system.file("extdata", file1, package = "dada2")
  fnR1 <- system.file("extdata", file2, package = "dada2")
  
  # Check the object class and location of file save
  checkObjects(fnF1)
  
  
  # setting up unique, temporary files, buried in the systems temporary directory
  filtF1 <- tempfile(fileext = ".fastq.gz")
  filtR1 <- tempfile(fileext = ".fastq.gz")
  
  # Check the location
  checkObjects(filtF1, filtR1)
  
  
  # From the help for the plotQualityProfile function: "The plotted lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles."
  #plot quality of forward sequences
  plotQualityProfile(fnF1)
  #plot quality of reverse sequences
  plotQualityProfile(fnR1)
  
  
  # When inspecting the quality of the forward and reverse sequences, it is clear that the reverse drop off in quality much faster than the forward reads, therefore we will trim the forward reads at 240, and the reverse at 200. It is also recommended to trim the initial reads as well since calibration issues can cause them to be problematic. The ambiguous nucleotides are filtered out with maxN = 0, and the reads with more than 2 expected errors are filtered out as well. 
  filterAndTrim(fwd = fnF1, filt = filtF1, rev = fnR1, filt.rev = filtR1,
                trimLeft = 10, truncLen = c(240, 200), 
                maxN = 0, maxEE = 2,
                compress = TRUE, verbose = TRUE)
  
  
  # To collapse reads that all encode the same sequence, we use dereplication. This reduces computational load and therefore, time.  
  derepF1 <- derepFastq(filtF1, verbose=TRUE)
  derepR1 <- derepFastq(filtR1, verbose=TRUE)
  
  
  # PCR amplification sometimes introduces errors, which vary between sequencing runs, so this is a way to estimate the parameters from the data
  errF <- learnErrors(derepF1, multithread = FALSE)
  errR <- learnErrors(derepR1, multithread = FALSE)
  
  
  # This infers the sample composition
  dadaF1 <- dada(derepF1, err = errF, multithread = FALSE)
  dadaR1 <- dada(derepR1, err = errR, multithread = FALSE)
  print(dadaF1)
  print(dadaR1)
  
  
  # Now, we merge the forward and reverse reads since above we inferred the sample sequences independently. This returns a data frame
  merger1 <- mergePairs(dadaF1, derepF1, dadaR1, derepR1, verbose = TRUE)
  checkObjects(merger1)
  
  
  # Chimeras are errors in PCR amplification and should be removed, so the following removes them. Chimeras are essentially a sequence read that is half of one sample sequence and half of another.
  merger1.nochim <- removeBimeraDenovo(merger1, multithread = FALSE, verbose = TRUE)
  
}


merger1 <- makeMerger("sam1F.fastq.gz", "sam1R.fastq.gz")
merger2 <- makeMerger("sam2F.fastq.gz", "sam2R.fastq.gz")


# We can now make a unified sequence table, with the combined, inferred samples. 

seqtab <- makeSequenceTable(list(merger1, merger2))
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = TRUE)

checkObjects(seqtab.nochim)




#From the tutorial: "This is the final product of the dada2 pipeline, a matrix in which each row corresponds to a processed sample, and each column corresponds to an non-chimeric inferred sample sequence (a more precise analogue to the common "OTU table")." And, each cell contains the count of sequences for that variant. "From here we recommend proceeding forward with our friend the phyloseq package for further analysis."


# POSSIBLY DO THIS:
#Build a phylogenetic tree from the unique sequences and mark the relative abundances for the variants in each of the two samples.
#Example package for ideas the visualization: https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

