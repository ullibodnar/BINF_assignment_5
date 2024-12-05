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


# This variable will store a series of parameters to test so that I don't have to hard code each one, it essentially creates a matrix of permutations for the below parameters, so we can see which ones affect eachother
parameters <- expand.grid(
  # I have chosen these truncation lengths by choosing extreme values on both ends. If the truncation is overly aggressive, it can discard valid sequcnes, which would provide improper results. I have also seperated them as forward and reverse to make downstream analysis easier.
  truncLenF = c(260, 240, 220, 200, 180),
  truncLenR = c(220, 200, 180, 160, 140),
  maxEE = c(2, 4), # This limits the amount of expected errors per read, 2 is the typical parameter, but 4 is considered a more relaxed filtering step, allowing for more sequences. However, it risks introducing more errors. 
  trimLeft = c(10, 20) # This determines the number of nucleotides to trim from the start of each read. I have chosen these values because 10 is typically the default, but a more agressive level, i.e., 20, can ensure that artifacts and low quality bases are removed. 
)
checkObjects(parameters) # Need to load functions FIRST for this to run, I kept it like this because it's cleaner to have the variables before the functions



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

# Create a function that can clean the parameters before being passed to the calculateDiversity function. When they are extracted, they are passed as a list, but I want them to be passed as.numeric so that they are in the expected format.
unlistNumeric <- function (x) {
  return (
    as.numeric(unlist(x))
  )
}


# Check object function adapted from assignment 3. When I edited Moiz's script I found that this function was a blessing
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
    truncLenF = 240,
    truncLenR = 200,
    maxEE = 2, 
    trimLeft = 10
  ) {
  
    
    # Concatenate truncLenF and R, and convert to numeric so that they go into the filterstats the correct way
    truncLen <- c(truncLenF, truncLenR) |> 
      as.numeric()
    
    
    # This is staying in here for debugging purposes, it uses a smaller dataset to run through this function, but is good to catch errors if there are any. 
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
    
   
    # When inspecting the quality of the forward and reverse sequences, it is clear that the reverse typically drop off in quality much faster than the forward reads, therefore we will often trim the forward reads at 240, and the reverse at 200. It is also recommended to trim the initial reads as well since calibration issues can cause them to be problematic. The ambiguous nucleotides are filtered out with maxN = 0, and the reads with more than 2 expected errors are filtered out as well. 
    filterStats <- tryCatch({ # Here we use a try catch since some of the combinations of parameters cause there to be no results, so we want a way to catch the error and proceed with the rest of the combinations
      filterAndTrim (
        fwd = fnFs, 
        filt = filtFs, 
        rev = fnRs, 
        filt.rev = filtRs,
        trimLeft = trimLeft, 
        truncLen = truncLen, 
        maxN = 0, # DADA2 requires no Ns
        maxEE = maxEE,
        compress = TRUE, 
        verbose = debug # We don't need an output if we aren't debugging
      )
    }, error = function (e) {
      warning("filterAndTrim likely failed to return any sequences with sequence parameters: ", conditionMessage(e))
      return(NULL)
    })
  
    
    # Check if any reads passed the filter
    if (is.null(filterStats) || all(filterStats[, 1] == 0 & filterStats[, 2] == 0)) {
      warning("No reads passed the filter for this parameter combination.")
      return(
        list(richness = NA, shannon = NA) # If no reads pass, we will just return NA for richness and shannon since calculating them for nothing makes no sense. In this case, a failed shannon and richness are marked by an NA, this also represents too stringent of a filter.
      )
    }
    
    
    # Ensure filtered files exist
    filtFs <- filtFs[file.exists(filtFs)]
    filtRs <- filtRs[file.exists(filtRs)]
    if (length(filtFs) == 0 || length(filtRs) == 0) {
      warning("No valid filtered files found. Skipping this combination.")
      return(
        list(richness = NA, shannon = NA) # Once again, failing to return these means failed parameters, we return NA and stop proceeding
      )
    }
    
    
    # To collapse reads that all encode the same sequence, we use dereplication. This reduces computational load and therefore, time.  
    derepFs <- derepFastq(filtFs[file.exists(filtFs)], verbose = debug)
    derepRs <- derepFastq(filtRs[file.exists(filtRs)], verbose = debug)
    
    
    # PCR amplification sometimes introduces errors, which vary between sequencing runs, so this is a way to estimate the parameters from the data
    errFs <- learnErrors(derepFs, multithread = FALSE)
    errRs <- learnErrors(derepRs, multithread = FALSE)
    
    
    # This infers the sample composition
    dadaFs <- dada(derepFs, err = errFs, multithread = FALSE)
    dadaRs <- dada(derepRs, err = errRs, multithread = FALSE)
    print(dadaFs)
    print(dadaRs)
    
    
    # Now, we merge the forward and reverse reads since above we inferred the sample sequences independently. This returns a data frame
    merger1 <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = debug)
    checkObjects(merger1)
    
    
    # Chimeras are errors in PCR amplification and should be removed, so the following removes them. Chimeras are essentially a sequence read that is half of one sample sequence and half of another.
    merger1.nochim <- removeBimeraDenovo(merger1, multithread = FALSE, verbose = debug)
  
    
    # We can now make a unified sequence table, with the combined, inferred samples. 
    seqtab <- makeSequenceTable(merger1.nochim)
    seqtab.nochim <- removeBimeraDenovo(seqtab, verbose = debug)
    
    
    # We can now calculate our richness and shannon index
    richness <- rowSums(seqtab.nochim > 0)
    shannon <- diversity(seqtab.nochim, index = "shannon")
    
    return(
      list(richness = richness, shannon = shannon) # If we've made it this far, we can finally return our shannon and richness indexes! 
    )
}


  
# Create a merger ---------------------------------------------------------

# From the help for the plotQualityProfile function: "The plotted lines show positional summary statistics: green is the mean, orange is the median, and the dashed orange lines are the 25th and 75th quantiles."
#plot quality of forward sequences
plotQualityProfile(forwardList[1:2])
#plot quality of reverse sequences
plotQualityProfile(reverseList[1:2])


# Analysis operations -------------------------------------------------------

# This is where the data and the parameter combinations will be passed to the calculateDiversity function to receive a list of metrics that can be further analyzed


# This is the one that should be used for debugging (or in Brittany's case, testing the script to ensure it runs)
resultsDebug <- apply(parameters[2, ], 1, function (row) {
  calculateDiversity(
    forwardList[18:20], 
    reverseList[18:20],
    truncLenF = unlistNumeric(row["truncLenF"]),
    truncLenR = unlistNumeric(row["truncLenR"]),
    maxEE = unlistNumeric(row["maxEE"]),
    trimLeft = unlistNumeric(row["trimLeft"]),
    debug = TRUE)
})


# This is the one that takes three hours to run since it is running through the entire calculate diversity for the entirety of the data, for each of the parameters (100 times)
# DO NOT RUN UNLESS YOU PLAN ON WATCHING TWO MOVIES IN THE MEANTIME
# results <- apply(parameters, 1, function (row) {
#   calculateDiversity(
#     forwardList, 
#     reverseList,
#     truncLenF = unlistNumeric(row["truncLenF"]),
#     truncLenR = unlistNumeric(row["truncLenR"]),
#     maxEE = unlistNumeric(row["maxEE"]),
#     trimLeft = unlistNumeric(row["trimLeft"]),
#     debug = FALSE)
# })


# Transform results list into a data frame so that it can be interpreted easier (with graphs/ggplot2, it is the preferred format) 
# NOT NECESSARY TO RUN, IF YOU REALLY WANT TO THEN YOU SHOULD CHANGE RESULTS INSTANCES TO resultsDebug
# results_df <- parameters |>
#   mutate(
#     richness = map_dbl(results, ~ ifelse(is.null(.x$richness), NA, mean(.x$richness, na.rm = TRUE))),
#     shannon = map_dbl(results, ~ ifelse(is.null(.x$shannon), NA, mean(.x$shannon, na.rm = TRUE)))
#   )
# checkObjects(results_df)

# Write to a csv since these results took 3 hours to generate, I don't want to put Brittany through that :)
# write.csv(results_df, "../data/results/results_with_diversity_metrics.csv", row.names = FALSE)

# Read in data
results_df <- read.csv("../data/results/results_with_diversity_metrics.csv")


results_df |>
  group_by(maxEE, trimLeft) |>
  summarise(
    mean_richness = mean(richness, na.rm = TRUE),
    mean_shannon = mean(shannon, na.rm = TRUE),
    .groups = "drop"
  )





# Plot richness and shannon -----------------------------------------------------------

# Plot richness by forward and reverse truncation
ggplot(results_df, aes(x = truncLenF, y = richness, color = factor(trimLeft))) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  facet_wrap(~ truncLenR, scales = "free_y") +
  labs(
    title = "Impact of Forward and Reverse Truncation Lengths on Richness",
    x = "Forward Truncation Length (truncLenF)",
    y = "Richness",
    color = "Trim Left"
  ) +
  theme_minimal() +
  annotate("rect", xmin = 220, xmax = 240, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "blue") +
  annotate("text", x = 220, y = 85, label = "Significant Region", size = 4)


# Plot a heatmap of richness by parameter combinations
ggplot(results_df, aes(x = truncLenF, y = truncLenR, fill = richness)) +
  geom_tile(color = "white") +
  facet_wrap(
    ~ maxEE + trimLeft,
    labeller = labeller(
      maxEE = function(value) paste("Error Limit =", value),  # Rename maxEE
      trimLeft = function(value) paste("Trimmed Bases =", value)  # Rename trimLeft
    )
  ) +
  scale_fill_gradient(
      low = "blue",
      high = "red",
      na.value = "grey",
      name = "Richness\n (NA in Grey)"
    ) +
  
  labs(
    title = "Heatmap of Richness by Parameter Combinations",
    x = "Forward Truncation Length",
    y = "Reverse Truncation Length",
    fill = "Richness"
  ) +
  theme_minimal(base_size = 14)


# Plot a heatmap of shannon diversity by parameter combinations
ggplot(results_df, aes(x = truncLenF, y = truncLenR, fill = shannon)) +
  geom_tile(color = "white") +
  facet_wrap(
    ~ maxEE + trimLeft, 
    labeller = labeller(
      maxEE = function(value) paste("Error Limit =", value),  # Rename maxEE
      trimLeft = function(value) paste("Trimmed Bases =", value)  # Rename trimLeft
    )
  ) +
  scale_fill_gradient(
    low = "blue",
    high = "red",
    na.value = "grey",
    name = "Shannon Diversity\n(NA in Grey)"
  ) +
  labs(
    title = "Heatmap of Shannon Diversity by Parameter Combinations",
    x = "Forward Truncation Length",
    y = "Reverse Truncation Length",
    fill = "Shannon Diversity"
  ) +
  theme_minimal(base_size = 14)


results_df |>
  mutate(failed = is.na(richness)) |>
  group_by(failed) |>
  summarise(count = n()) |>
  ggplot(aes(x = failed, y = count, fill = failed)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Counts of Parameter Combinations with Failed Filtering",
    x = "Failed Filtering (NA)",
    y = "Count",
    fill = "Failed"
  ) +
  theme_minimal()


# Assign taxonomy ---------------------------------------------------------

taxa <- assignTaxonomy(seqtab.nochim, "../data/MiSeq_SOP/silva_nr_v132_train_set.fa.gz", multithread = TRUE)
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)





#From the tutorial: "This is the final product of the dada2 pipeline, a matrix in which each row corresponds to a processed sample, and each column corresponds to an non-chimeric inferred sample sequence (a more precise analogue to the common "OTU table")." And, each cell contains the count of sequences for that variant. "From here we recommend proceeding forward with our friend the phyloseq package for further analysis."


# POSSIBLY DO THIS:
#Build a phylogenetic tree from the unique sequences and mark the relative abundances for the variants in each of the two samples.
#Example package for ideas the visualization: https://guangchuangyu.github.io/ggtree-book/chapter-ggtree.html

