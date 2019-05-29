library('tidyr')

# Parse metadata for RGP
RGP_families = read.csv('data/RGP_MacArthur_RareDisease_WGS_Metadata_2019-03-15.tsv', sep = '\t')
RGP_sample_info = read.csv('data/RGP_MacArthur_RareDisease_WGS_samples_2019-03-15.tsv', sep = '\t')
RGP_unaffected = read.csv('data/RGP_unaffected_parents.txt', sep = '\t', header = F, col.names = 'sample')

gender_df = RGP_families[,c('participant_id', 'gender')]
colnames(gender_df) = c('sample', 'gender')

RGP_unaffected = merge(RGP_unaffected, gender_df, all.x = T)
