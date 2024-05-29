# Load necessary libraries
library(SummarizedExperiment)
library(DESeq2)

# Example count data (rows are genes, columns are samples)
countData <- matrix(rpois(200, lambda = 10), ncol=4)
colnames(countData) <- c("sample1", "sample2", "sample3", "sample4")

# Sample metadata
colData <- data.frame(
  row.names = colnames(countData),
  condition = factor(c("A", "A", "B", "B")),
  dex = factor(c("trt", "trt", "untrt", "untrt"))
)

# Create a SummarizedExperiment object
se <- SummarizedExperiment(assays = list(counts = countData), colData = colData)

# Create a DESeqDataSet object
dds <- DESeqDataSet(se, design = ~ dex)

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Extract and inspect results
res <- results(dds)
head(res)

# Plot the results
plotMA(res, main="DESeq2", ylim=c(-2, 2))

# Export results to a CSV file
write.csv(as.data.frame(res), file="DESeq2_results.csv")
