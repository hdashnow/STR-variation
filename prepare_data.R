library('MASS')
library('data.table')
library('dplyr')
library("ggrepel")

source('STR-variation-functions.R')

### Settings
data_dir = '/group/bioi1/harrietd/STR_variation/working_dir/old/'
#annotation_file = '/group/bioi1/harrietd/STR_variation/working_dir2/RGP_149_1.STRs.annotated.tsv'

# Annotation files
annotation_file = '/group/bioi1/harrietd/STR_variation/working_dir/hg19.simpleRepeat_period1-6_dedup.sorted.annotated.tsv'
gnomad_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/gnomad.v2.1.1.lof_metrics.by_gene.txt'
disease_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/curated_gene_disease_associations.tsv'
omim_genes_file = '/group/bioi1/harrietd/git/STRetch/scripts/tests/test_data/mim2gene.txt'

all_annotations = get_annotations(annotation_file, gnomad_genes_file, disease_genes_file, omim_genes_file)

table(all_annotations$in_omim, all_annotations$disease_gene, useNA = 'ifany', dnn = c('in omim', 'disease gene'))

# Should be 362 samples total across the two datasets - need to combine them
### Parse/clean up data
RGP_219.STRs = read.csv(paste0(data_dir,"RGP_219.STRs.annotated.tsv"), sep='\t', stringsAsFactors = F)
RGP_219.STRs$source = 'RGP'
controls_143.STRs = read.csv(paste0(data_dir,"PCRfreeWGS_143_STRetch_controls.STRs.annotated.tsv"), sep='\t', stringsAsFactors = F)
controls_143.STRs$source = 'controls'

all_STR_calls = dplyr::bind_rows(RGP_219.STRs, controls_143.STRs)
all_STR_calls = clean_annotations(all_STR_calls) # Set missing values to NA
all_STR_calls$locus = paste(all_STR_calls$chrom, all_STR_calls$start, all_STR_calls$end, sep = '_')

### Calculate summary stats

# Use hubers, not huber because huber is just using MAD to estimate scale then using that to estimate mu
# Note locuscoverage_log is already coverage-normalised
all_STR_calls %>% group_by(locus) %>%  
  do(data.frame(hubers(.$locuscoverage_log))) ->
  hubers_by_locus
setnames(hubers_by_locus, old=c("mu","s"), new=c("mu_locuscoverage_log", "s_locuscoverage_log"))

# Calculate number of significant calls per locus at various thresholds
all_STR_calls %>% 
  group_by(locus) %>% 
  summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
            n_zeros = sum(locuscoverage == 0), 
            total = length(p_adj)) ->
  STR_calls_locus

# Add annotations to the summary stats
STR_calls_locus = merge(STR_calls_locus, hubers_by_locus)
STR_calls_locus = merge(STR_calls_locus, all_annotations, all.x = T)

# Caclulate proportions of STRs that are variable
all_annotated_STRs_locus_calls = merge(all_annotations, STR_calls_locus, all.x = T)
