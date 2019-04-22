library('tidyr')

# Parse STR data (takes a while!)
#RGP_data = read.csv('data/RGP.STRs.tsv', sep = '\t')
RGP_data = read.csv('data/RGP_219.STRs.tsv', sep = '\t')
control_data = read.csv('data/PCRfreeWGS_143_STRetch_controls.STRs.tsv', sep = '\t')
# Parse metadata for RGP
RGP_families = read.csv('data/MacArthur_RGP_samples_metadata.txt', sep = '\t')
RGP_sample_info = read.csv('data/MacArthur_RGP_samples_summary.txt', sep = '\t')

#Extract per-sample coverage from STRetch results
RGP_genomecov = unique(RGP_data[c('sample','genomecov')])
RGP_genomecov$project = 'RGP'
control_genomecov = unique(control_data[c('sample','genomecov')])
control_genomecov$project = 'other'
all_genomecov = rbind(RGP_genomecov, control_genomecov)

# Transform STR data to a counts matrix
# The arguments to spread():
# - data: Data object
# - key: Name of column containing the new column names
# - value: Name of column containing values
all_data = rbind(RGP_data, control_data)
all_data$locus = paste(all_data$chrom, all_data$start, all_data$end, sep = '-')

# Extract locus information
### not finished
locus_data = all_data[c('locus', 'chrom','start','end', 'repeatunit', 'reflen')]
locus_data = locus_data[!duplicated(locus_data$locus),]
#locus_data = unique(locus_data)
dim(all_data)
dim(locus_data)

data_subset = all_data[c('locus','sample', 'locuscoverage')]
data_wide <- spread(data_subset, key = sample, value = locuscoverage)
# Make locus column the rownames, then remove the locus column
rownames(data_wide) <- data_wide$locus
data_wide = subset(data_wide, select=-c(locus))

save(RGP_genomecov, control_genomecov, all_genomecov, 
     data_wide, locus_data, file = 'STR_data.RData')