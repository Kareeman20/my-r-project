#installing Libraries
BiocManager::install("Biostrings")
install.packages("tidyverse")

library(Biostrings)
library(tidyverse)

# Create a DNA string object
dna_seq<-DNAString("ATGCCGTAGCCTGA")
print(paste("Original DNA:", dna_seq))


# Translate the DNA to protein
protein <- translate(dna_seq)
print(paste("Protein Sequence:", protein))


# Get reverse complement of the DNA sequence
rev_comp <- reverseComplement(dna_seq)
print(paste("Reverse Complement:", rev_comp))
#3`TACGGCATCGGACT 5`
#5` TCAGGCTACGGCAT 3`


# Calculate GC content
gc_perc <- letterFrequency(dna_seq, letters = "GC", as.prob = TRUE) * 100
print(paste("GC Content:", round(gc_perc, 2), "%"))   

#create a dataframe
gene_data <- data.frame(
  gene = c("TP53", "BRCA1", "EGFR", "MYC", "PTEN", "GAPDH"),
  expression = c(12.5, 45.2, 8.9, 33.1, 2.4, 105.0), # measure how active the gene is 
  p_value = c(0.01, 0.002, 0.45, 0.001, 0.04, 0.99),#<0.05 means the result is "statistically significant"
  type = c("Tumor", "Tumor", "Normal", "Tumor", "Normal", "Control") # categoral 
)


# Filter genes with p-value < 0.05
significant_genes <- gene_data %>%
  filter(p_value < 0.05)


# Create log2 of expression
significant_genes <- significant_genes %>%
  mutate(log_exp = log2(expression)) #for normalization 

# Sort genes by log expression in descending order
significant_genes <- significant_genes %>%
  arrange(desc(log_exp))


print("Significant Genes (Sorted):")
print(significant_genes)

#Create a DNA sequence set
seq_set <- DNAStringSet(c(
  "ATGCGT", "ATGCGC", "ATGCGG", "ATGCGA", "ATGCGC", "ATGCGT"
))
names(seq_set) <- gene_data$gene

# Convert sequence set to a data frame
seq_info <- data.frame(
  gene = names(seq_set),
  sequence = as.character(seq_set),
  width = width(seq_set)
)

# Merge gene data with sequence info
final_combined_data <- left_join(gene_data, seq_info, by = "gene")

print("Combined Gene Expression and Sequence Data:")
print(final_combined_data)

# Group by 'type' and calculate mean expression
group_summary <- gene_data %>%
  group_by(type) %>%
  summarize(
    avg_expression = mean(expression),
    gene_count = n()
  )

print("Group Summary Statistics:")
print(group_summary)
