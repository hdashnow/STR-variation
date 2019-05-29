library('MASS')
library('data.table')
library('dplyr')
library("ggrepel")

source('STR-variation-functions.R')

### Settings
# STRetch calls
all_STR_calls_csv = '/group/bioi1/harrietd/STR_variation/control_sets/all_samples/STRs.tsv'

# Annotation files
annotation_file = '/group/bioi1/harrietd/STR_variation/working_dir/hg19.simpleRepeat_period1-6_dedup.sorted.annotated.tsv'
gnomad_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/gnomad.v2.1.1.lof_metrics.by_gene.txt'
disease_genes_file = '/group/bioi1/harrietd/STR_variation/annotation_data/curated_gene_disease_associations.tsv'
omim_genes_file = '/group/bioi1/harrietd/git/STRetch/scripts/tests/test_data/mim2gene.txt'
pathogenic_csv = 'data/STR_disease_loci_ranges.csv'

pathogenic_info = read.csv(pathogenic_csv, stringsAsFactors = F)

all_annotations = get_annotations(annotation_file, gnomad_genes_file, disease_genes_file, omim_genes_file)

# Read, clean up and annotate STRetch calls
all_STR_calls = read_annotate_calls(all_STR_calls_csv, all_annotations)

### Calculate summary stats and annotate them
STR_calls_locus = summary_stats(all_STR_calls, all_annotations)

# Calculate variance (s_locuscoverage_log) deciles over all loci
num_breaks = 10
s_quantiles = quantile(STR_calls_locus$s_locuscoverage_log, seq(0, 1, 1/num_breaks), na.rm = T)
print(s_quantiles)
STR_calls_locus = within(STR_calls_locus,
                         s_decile <- as.integer(cut(s_locuscoverage_log,
                                                    quantile(s_locuscoverage_log,
                                                             seq(0, 1, 1/num_breaks), na.rm = T),
                                                    include.lowest=TRUE))
)
var_quantiles = quantile(STR_calls_locus$var_locuscoverage_log, seq(0, 1, 1/num_breaks), na.rm = T)
print(var_quantiles)
STR_calls_locus = within(STR_calls_locus,
                         var_decile <- as.integer(cut(var_locuscoverage_log,
                                                    quantile(var_locuscoverage_log,
                                                             seq(0, 1, 1/num_breaks), na.rm = T),
                                                    include.lowest=TRUE))
)

# Create a df with all annotations, and the locus calls for those that are variable
all_annotated_STRs_locus_calls = merge(all_annotations, STR_calls_locus, all.x = T)

# Count total number of genes in the gencode annotation (result from a bash script)
gencode_genes = 57781
