# Preparing fission data for lesson

library(fission)
library(DESeq2)
library(magrittr)
library(tibble)

data("fission")

# modify colnames of the data to be a bit less confusing
colnames(fission) <- colData(fission)$id

# Create DESeq dataset - simple design, no interaction
dds <- DESeqDataSet(fission, design = ~ strain + minute)

# Remove very low count genes - at least half the samples with more than 5 reads
genes_to_keep <- rowSums((counts(dds) > 5)) > 18
dds <- dds[genes_to_keep, ]


#
# Run tests ----
#
# Run DESeq model
dds <- DESeq(dds)

# Make a contrast between first and other time points for WT
test_result <- lapply(c("15", "30", "60", "120", "180"), function(t){
  res <- results(dds, contrast = c("minute", t, "0"), tidy = TRUE, lfcThreshold = 1) %>%
    as_tibble()
  res$comparison <- as.numeric(t)
  return(res)
})
test_result <- do.call("rbind", test_result)

# Rename first column
names(test_result)[1] <- "gene"


#
# Gene counts ----
#
# Extract raw counts
raw_cts <- counts(dds, normalized = TRUE)

# Applying normalization to data (ignoring design) - for clustering, etc.
norm_cts <- vst(dds, blind = TRUE) %>% assay()

# Get sample information
sample_info <- colData(dds) %>% as.data.frame() %>% as_tibble(rownames = "sample")
sample_info <- sample_info[, c("sample", "strain", "minute", "replicate")]  # retain only a few columns


#
# Save objects ----
#
# save(raw_cts, norm_cts, sample_info, test_result,
#      file = "data/fission_data.RData", compress = "bzip2")

raw_cts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv("./data/counts_raw.csv")

norm_cts %>%
  as_tibble(rownames = "gene") %>%
  readr::write_csv("./data/counts_transformed.csv")

readr::write_csv(sample_info, "./data/sample_info.csv")
readr::write_csv(test_result, "./data/test_result.csv")
