library('ggplot2'); theme_set(theme_classic(base_size = 16))
library('scales')
library('RColorBrewer')
#library('plyr') <- having both of these loaded at once causes problems because they both contain summarise
library('dplyr')
library("ggrepel")

#source('prepare_data.R')

# STR disease loci pathogenic lengths plot
# *intermediate = premutation/incomplete penetrance/uncertain pathogenicity
ggplot(data=subset(pathogenic_info, pathogenic_min > 0),
       aes(x = disease)) +
  geom_linerange(aes(ymin = normal_min*repeatunitlen, ymax = normal_max*repeatunitlen+1, colour = 'normal'), size = 1) + 
  geom_linerange(aes(ymin = intermediate_min*repeatunitlen, ymax = intermediate_max*repeatunitlen+1, colour = 'intermediate*'), size = 1) + 
  geom_linerange(aes(ymin = pathogenic_min*repeatunitlen, ymax = pathogenic_max*repeatunitlen+1, colour = 'pathogenic'), size = 1) + 
  scale_y_log10(breaks=c(10,100,1000,10000)) +
  labs(y = 'allele size (bp)') +
  scale_colour_manual(values=brewer.pal(n = 3, name = "Set1")[c(2,3,1)], 
                      breaks = c("normal","intermediate*","pathogenic"),
                      name = 'Pathogenicity') +
  theme_bw(base_size = 16) +
  theme(panel.grid.major.y = element_blank()) + 
  coord_flip() +
  ggsave('plots/pathogenic_allele_sizes.jpg', height = 6, width = 7)

standard_chroms = paste('chr', c(1:22, 'X', 'Y'), sep = '')

# Look at some pathogenic loci
RGP_sig_pathogenic_loci = all_STR_calls[!is.na(all_STR_calls$pathogenic) & all_STR_calls$project == 'RGP' & all_STR_calls$p_adj < 0.05,]


### Some general stats about STRs in the genome

# Number of STRs in the STRetch annotation file
dim(all_annotated_STRs_locus_calls)[1]

# Set chromosome plotting order
all_annotated_STRs_locus_calls$chrom = factor(all_annotated_STRs_locus_calls$chrom, levels = standard_chroms)

hg19.genome = read.table('/group/bioi1/harrietd/git/STRetch/reference-data/hg19.STRdecoys.sorted.fasta.genome', 
                         col.names = c('chrom', 'size_bp'))
hg19.genome = hg19.genome[hg19.genome$chrom %in% standard_chroms,]
hg19.genome$chrom = factor(hg19.genome$chrom, levels = standard_chroms)

# Plot number of bp in each chromosome
ggplot(data = hg19.genome, aes(x = chrom, y = size_bp)) + geom_bar(stat = "identity")

# plot number of STRs annotated in each chomosome
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = chrom)) +
  geom_bar(aes(fill = factor(repeatunitlen)))

# plot number of STRs with each repeat unit length
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = repeatunitlen, fill = repeatunit)) +
  geom_bar() + 
  scale_x_continuous(breaks = 1:6) + 
  theme(legend.position = "none") + ggsave('plots/ref_repeatunitlen_counts.jpg')

repeatunit_counts = sort(table(all_annotations$repeatunit), decreasing = T)
repeatunit_counts[c('A', 'T', 'C', 'G')]
sum(repeatunit_counts[c('A', 'T')])/sum(repeatunit_counts) * 100
sum(repeatunit_counts[c('C', 'G')])/sum(repeatunit_counts) * 100

# How many STRs annotated with the telomere sequence?
repeatunit_counts['AACCCT']
repeatunit_counts['AACCCT']/sum(repeatunit_counts) * 100

# Counts/percent of STRs overlapping CDS/exons
feature_table = table(all_annotations$feature)
sum(feature_table[c('CDS')])
sum(feature_table[c('CDS')])/sum(feature_table) * 100
sum(feature_table[c('exon')])
sum(feature_table[c('exon')])/sum(feature_table) * 100

# How many genes have STRs?
length(unique(all_annotations$gene_name))
length(unique(all_annotations$gene_name))/gencode_genes * 100
length(unique(all_annotations$gene_name[all_annotations$feature == 'CDS']))
length(unique(all_annotations$gene_name[all_annotations$feature == 'CDS']))/gencode_genes * 100
length(unique(all_annotations$gene_name[all_annotations$feature == 'exon']))
length(unique(all_annotations$gene_name[all_annotations$feature == 'exon']))/gencode_genes * 100

# How many genes have STRs within X bp of their TSS?
x_tss = -3000
count_tss_genes = length(unique(all_annotations[(all_annotations$closest_TSS_distance > x_tss) & (all_annotations$closest_TSS_distance <= 0),]$gene_name))
count_tss_genes
count_tss_genes/gencode_genes * 100

# How many OMIM genes have STRs?
length(unique(all_annotations$gene_name[all_annotations$in_omim]))
length(unique(all_annotations$gene_name[all_annotations$in_omim]))/gencode_genes * 100
length(unique(all_annotations$gene_name[(all_annotations$in_omim) & (all_annotations$feature == 'CDS')]))
length(unique(all_annotations$gene_name[(all_annotations$in_omim) & (all_annotations$feature == 'CDS')]))/gencode_genes * 100

# number of STRs in various gene reigons
ggplot(data = all_annotated_STRs_locus_calls,
       aes(x = feature, fill=feature)) +
  geom_bar() +
  facet_wrap(~repeatunitlen) 

# first make a dataframe with frequencies
df <- as.data.frame(with(all_annotated_STRs_locus_calls, table(feature, repeatunitlen)))
# next: compute percentages per group
df <- plyr::ddply(df, plyr::.(repeatunitlen), transform, p = Freq/sum(Freq))
# and plot
ggplot(df, aes(feature, p, fill=feature)) +
  geom_bar(stat = 'identity') + 
  facet_wrap(~repeatunitlen) +
  scale_y_continuous(labels=percent_format()) +
  labs(y = '% of STR loci with repeat unit length overlapping feature') + 
  ggsave('plots/percent_loci_overlap_features.jpg', width = 14, height = 6)

# Counts for the above plots:
table(all_annotated_STRs_locus_calls$feature, all_annotated_STRs_locus_calls$repeatunitlen)

# How many STRs overlap OMIM genes?
table(all_annotated_STRs_locus_calls$in_omim, useNA = 'ifany')
table(all_annotated_STRs_locus_calls$in_omim)['TRUE']/length(all_annotated_STRs_locus_calls$in_omim) * 100
table(all_annotated_STRs_locus_calls$in_omim, all_annotated_STRs_locus_calls$feature, useNA = 'ifany', 
      dnn = c('in OMIM', 'feature'))

# How many STRs upstream of an OMIM TSS?
# This calculation doesn't make sense because in_omim refers to the gene column, 
# but the TSS may not be referring to a transcript of that gene.
# x_tss = -1000
# is_near_tss = (all_annotated_STRs_locus_calls$closest_TSS_distance > x_tss) & (all_annotated_STRs_locus_calls$closest_TSS_distance <= 0)
# dim(all_annotated_STRs_locus_calls[all_annotated_STRs_locus_calls$in_omim & is_near_tss,])

### How many STR loci show evidence of expansion in our sample?

# Define evidence of expansion as having at least one read reported for that locus in STRetch:
all_annotated_STRs_locus_calls$has_counts = !is.na(all_annotated_STRs_locus_calls$total)
table(all_annotated_STRs_locus_calls$has_counts)
table(all_annotated_STRs_locus_calls$has_counts)['TRUE']/length(all_annotated_STRs_locus_calls$has_counts) * 100

# Has at least one significant call from STRetch
table(all_annotated_STRs_locus_calls$signif_0_05 > 0)
table(all_annotated_STRs_locus_calls$signif_0_05 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_05) * 100
table(all_annotated_STRs_locus_calls$signif_0_01 > 0)
table(all_annotated_STRs_locus_calls$signif_0_01 > 0) ['TRUE']/length(all_annotated_STRs_locus_calls$signif_0_01) * 100

### What number/proportion of individuals have STR expansions?
all_STR_calls %>% 
  group_by(sample, project) %>% 
  dplyr::summarise(signif_0_01 = sum(p_adj < 0.01), signif_0_05 = sum(p_adj < 0.05), 
            non_zero = sum(locuscoverage > 0), genomecov = mean(genomecov)
            ) -> all_STR_calls_sample
summary(all_STR_calls_sample)

ggplot(data = all_STR_calls_sample, aes(non_zero, fill = project)) + geom_histogram() + 
  scale_x_log10() +
  labs(x='Number of STR loci with any reads assigned', y='Number of samples')
ggplot(data = all_STR_calls_sample, aes(signif_0_05, fill = project)) + geom_histogram() + 
  scale_x_log10() +
  labs(x='Number of STR loci significant a p<0.05', y='Number of samples')
ggplot(data = all_STR_calls_sample, aes(signif_0_01, fill = project)) + geom_histogram() + 
  scale_x_log10() +
  labs(x='Number of STR loci significant a p<0.01', y='Number of samples') + 
  ggsave('plots/n_signif_01_hist.jpg')

#sequenced_together = read.table('/group/bioi1/harrietd/STR_variation/G84373_samples.txt', stringsAsFactors = F)[[1]]

# A subset of the control samples have a surprisingly large number of expanded STRs
high_var_samples = all_STR_calls_sample[all_STR_calls_sample$non_zero > 10000,]$sample
low_var_samples = all_STR_calls_sample[all_STR_calls_sample$non_zero <= 10000,]$sample

# Write the high and low var sample IDs to file
# write.table(high_var_samples,
#             file = 'high_var_samples.txt', row.names = F, col.names = F, quote = F)
# write.table(low_var_samples,
#             file = 'low_var_samples.txt', row.names = F, col.names = F, quote = F)
# write.table(hyper_samples,
#             file = 'high_var_samples.txt', row.names = F, col.names = F, quote = F)


# all_STR_calls_sample[all_STR_calls_sample$non_zero < 10000 & all_STR_calls_sample$project == 'controls',]$sample
# 
# tmp_df = cbind(high_var_samples, sequenced_together)

ggplot(data = all_STR_calls_sample, aes(y = non_zero, x = genomecov, colour = project)) + geom_point() + 
  labs(y='STRs with any reads assigned', x='Median sequencing depth') + 
  ggsave('plots/n_non_zero_vs_genomecov.jpg')

ggplot(data = all_STR_calls_sample, aes(y = signif_0_01, x = genomecov, colour = project)) + geom_point() + 
  labs(y='Significant STRs at p < 0.01', x='Median sequencing depth') + 
  ggsave('plots/n_signif_01_vs_genomecov.jpg')

all_STR_calls_sample$seq_depth = NA
all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov < 35] = '<35'
all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov >= 35] = '35-40'
all_STR_calls_sample$seq_depth[all_STR_calls_sample$genomecov > 40] = '>40'
# Set plotting order
all_STR_calls_sample$seq_depth = factor(all_STR_calls_sample$seq_depth, levels = c('<35', '35-40', '>40'))

ggplot(data = all_STR_calls_sample, aes(genomecov, fill = seq_depth)) + geom_histogram()
ggplot(data = all_STR_calls_sample, aes(non_zero, fill = seq_depth)) + geom_histogram()

# Most of the samples have zero counts at most loci
n_samples = length(unique(all_STR_calls_sample$sample))
ggplot(data = STR_calls_locus, aes(n_samples-n_zeros)) + 
  geom_histogram() + 
  labs(x='Num. samples with any reads assigned', y='Number of STR loci') +
  ggsave('plots/n_samples_w_reads_hist.jpg')

ggplot(data = STR_calls_locus, aes(signif_0_05)) + geom_histogram() 

ggplot(data = STR_calls_locus, aes(signif_0_01)) + geom_histogram() +
  labs(x='Num. samples signif. at p < 0.01', y='Number of STR loci') +
  ggsave('plots/n_signif_01_hist.jpg')

### What number/proportion of loci in feature categories have STR expansions?
all_annotated_STRs_locus_calls %>% 
  group_by(feature) %>% 
  dplyr::summarise(any_signif_0_01 = sum(signif_0_01 > 0, na.rm = T), 
                   any_signif_0_05 = sum(signif_0_05 > 0, na.rm = T), 
                   num_non_zero = sum(has_counts, na.rm = T), total = length(has_counts)
  ) -> all_annotated_STRs_locus_calls_feature


# Are STRs overlapping certain features more/less likely to be significant/expanded?

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = any_signif_0_01/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = any_signif_0_05/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

ggplot(data = all_annotated_STRs_locus_calls_feature, 
       aes(x = feature, y = num_non_zero/total)) + 
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels=percent_format()) 

# Gene set test
# x = number of significant calls within the category
# m = number of significant genes overall
# n = number of non-siginficant genes overall
# k = number of loci in the category

all_annotated_STRs_locus_calls_feature$num_zero = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$num_non_zero
all_annotated_STRs_locus_calls_feature$not_signif_0_01 = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$any_signif_0_01
all_annotated_STRs_locus_calls_feature$not_signif_0_05 = all_annotated_STRs_locus_calls_feature$total - all_annotated_STRs_locus_calls_feature$any_signif_0_05

# Are there more significant genes in introns?
gene_set_test = function(df, feature, col_name){
  x = df[[col_name]][df$feature == feature]
  m = sum(df[[col_name]])
  n = sum(df$total) - m
  k = df$total[df$feature == feature]
  
  phyper(x,m,n,k)
}
for (col_name in c("any_signif_0_01", "any_signif_0_05", "num_non_zero","num_zero") ){
  print(paste0('test column: ', col_name))
  for (feature in all_annotated_STRs_locus_calls_feature$feature){
    print(feature)
    print(gene_set_test(all_annotated_STRs_locus_calls_feature, feature, col_name))
  }
}






# Disease gene and constraint measures

ggplot(data = all_annotations,
       aes(x = pLI, y = oe_lof, colour = disease_gene)) +
  geom_point() +
  scale_y_log10() +
  geom_vline(xintercept = 0.9, colour = 'blue') +
  geom_hline(yintercept = 0.35, colour = 'blue')

ggplot(data = all_annotations,
       aes(x = oe_lof, y = oe_lof_upper, colour = disease_gene)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, colour = 'blue') +
  xlim(NA,2) 




# Check that my measures of disease genes are relating to each other in a sensible way.
# In this dataset the 0.35 oe_lof value is in the 4th decile
ggplot(data = all_annotations %>% group_by(oe_decile) %>% summarise(
  disgenet= sum(disease_gene, na.rm = T)/length(disease_gene),
  omim = sum(in_omim, na.rm = T)/length(in_omim)
) %>% melt('oe_decile', variable.name = "database", value.name = "prop_disease_genes")) +
  geom_point(aes(x = oe_decile, y = prop_disease_genes, colour = database))

# pLI doesn't work as well (at least with these gene databases and this particular subet of genes)
ggplot(data = all_annotations %>% group_by(pLI_decile) %>% summarise(
  disgenet= sum(disease_gene, na.rm = T)/length(disease_gene),
  omim = sum(in_omim, na.rm = T)/length(in_omim)
) %>% melt('pLI_decile', variable.name = "database", value.name = "prop_disease_genes")) +
  geom_point(aes(x = pLI_decile, y = prop_disease_genes, colour = database))

# Will use oe_lof_upper or oe_decile from now on


# Plot mean~variance

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = signif_0_01)) 

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, signif_0_01 > 0),
             aes(colour = signif_0_01)) +
  ggsave('plots/locuscov_mean_var.jpg')

# What are the repeatunits of the top 10 loci in terms of variance and median locuscoverage?
top_s = unique(head(STR_calls_locus[order(STR_calls_locus$s_locuscoverage_log, decreasing = T),], 10)$repeatunit)
top_mu = unique(head(STR_calls_locus[order(STR_calls_locus$mu_locuscoverage_log, decreasing = T),], 10)$repeatunit)
non_a = unique(c(top_s, top_mu))
non_a = non_a[non_a != 'A']

sum(STR_calls_locus$repeatunit %in% unique(c(top_s, top_mu)))/length(STR_calls_locus$repeatunit)
sum(all_annotations$repeatunit %in% unique(c(top_s, top_mu)))/length(all_annotations$repeatunit)

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, repeatunit %in% c('A')),
             aes(colour = repeatunit)) +
  geom_point(data = subset(STR_calls_locus, repeatunit %in% non_a),
           aes(colour = repeatunit)) +
  ggsave('plots/locuscov_mean_var_repeatunits.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, repeatunit %in% c('AGC')),
             aes(colour = repeatunit)) +
  ggsave('plots/locuscov_mean_var_CAG.jpg')

# How many CAG loci do we observe reads from?
length(unique(STR_calls_locus[STR_calls_locus$repeatunit == 'AGC',]$locus))
# 7 of these are pathogenic
unique(STR_calls_locus[STR_calls_locus$repeatunit == 'AGC',]$pathogenic)
# How many CAG loci in genome?
length(unique(all_annotations[all_annotations$repeatunit == 'AGC',]$locus))
unique(all_annotations[all_annotations$repeatunit == 'AGC',]$pathogenic)

ggplot(data = STR_calls_locus, 
       aes(x = mean_locuscoverage_log, y = var_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, signif_0_01 > 0),
             aes(colour = signif_0_01)) +
  labs(subtitle = paste0('N=',n_samples,' samples'))

# Check if the robust and non-robust mean and variance estimates are similar

ggplot(data = STR_calls_locus, 
       aes(y = mean_locuscoverage_log, x = mu_locuscoverage_log)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, signif_0_01 > 0),
             aes(colour = signif_0_01)) + 
  labs(y = 'Mean log-norm locuscounts', x = 'Huber robust median log-norm locuscounts',
    subtitle = paste0('N=',n_samples,' samples'))

ggplot(data = STR_calls_locus, 
       aes(y = var_locuscoverage_log, x = s_locuscoverage_log)) +
  geom_point(size = 1) +
  geom_abline(intercept = 0, slope = 1, colour = 'blue') +
  labs(y = 'Variance log-norm locuscounts', x = 'Huber robust variance log-norm locuscounts') + 
  ggsave('plots/locuscov_var_s.jpg')

# More mean~variance plots

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = 'other')) +
  geom_point(data = subset(STR_calls_locus, known_pathogenic),
             aes(colour = "pathogenic")) + 
  scale_colour_manual(values=brewer.pal(n = 9, name = "Set1")[c(9,1)], 
                      breaks = c("pathogenic", "other"),
                      name = 'Locus') +
  ggsave('plots/locuscov_mean_var_disease.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log)) +
  geom_point(aes(colour = 'other')) +
  geom_point(data = subset(STR_calls_locus, known_pathogenic),
             aes(colour = "pathogenic")) + 
  scale_colour_manual(values=brewer.pal(n = 9, name = "Set1")[c(9,1)], 
                      breaks = c("pathogenic", "other"),
                      name = 'Locus') +
  geom_text_repel(data = subset(STR_calls_locus, known_pathogenic & s_locuscoverage_log > 0.14),
                  aes(label=pathogenic)
  ) +
  xlim(NA, 1.9) + 
  ylim(NA, 0.6) + 
  ggsave('plots/locuscov_mean_var_disease_zoom.jpg')

ggplot(data = STR_calls_locus, 
       aes(x = mu_locuscoverage_log, y = s_locuscoverage_log, colour = feature)) +
  geom_point(colour = 'grey') +
  geom_point(data = subset(STR_calls_locus, feature == 'exon')) +
  geom_point(data = subset(STR_calls_locus, feature == 'CDS')) + 
  ggsave('plots/locuscov_mean_var_exons.jpg')
             
   #+
  #xlim(NA, 1.9) + 
  #ylim(NA, 0.6)


# Plot observed/expected lof

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, colour = in_omim)) +
  geom_density() +
  geom_vline(xintercept = 0.35, colour = 'blue')

ggplot(data = STR_calls_locus,
       aes(x = oe_lof_upper, colour = disease_gene)) +
  geom_density() +
  geom_vline(xintercept = 0.35, colour = 'blue')


# # STR_calls_locus %>% group_by(oe_decile) %>% summarise(
# #   v = median(mu_locuscoverage_log, na.rm = T)) -> tmp
# # tmp = tmp[complete.cases(tmp),]
# ggplot(data = subset(STR_calls_locus, feature %in% c('exon', 'CDS')),
#        aes(x = oe_decile, y = s_locuscoverage_log, group = oe_decile)) +
#   geom_violin() + ylim(NA, 0.15)
# #ggsave('plots/oeLOF_signif')

# ggplot(data = STR_calls_locus,
#        aes(x = oe_lof_upper, y = s_locuscoverage_log, colour = in_omim)) +
#   geom_point() +
#   xlim(NA, 2.5) +
#   geom_vline(xintercept = 0.35, colour = 'blue')
# 
# ggplot(data = STR_calls_locus,
#        aes(x = oe_lof_upper, y = s_locuscoverage_log, colour = disease_gene)) +
#   geom_point() +
#   xlim(NA, 2.5) +
#   geom_vline(xintercept = 0.35, colour = 'blue')


### plot s_deciles vs features

# Calculate s_locuscoverage_log for loci with 0 counts in all samples
current_samples = unique(all_STR_calls$sample)
#s_zero_locus = hubers(log2(1/all_genomecov$genomecov[all_genomecov$sample %in% current_samples]*100))$s

all_annotated_STRs_locus_calls$s_decile[is.na(all_annotated_STRs_locus_calls$s_decile)] = 0

s_decile_counts = table(all_annotated_STRs_locus_calls$s_decile)
s_decile_prop = s_decile_counts/sum(s_decile_counts)

# Plot variance (s_locuscoverage_log) deciles over all loci
num_deciles = 11
STR_loci_feature = table(all_annotated_STRs_locus_calls$feature, useNA = "ifany")
# Observed counts for each feature
all_annotated_STRs_locus_calls$feature = factor(all_annotated_STRs_locus_calls$feature)
all_annotated_STRs_locus_calls %>% 
  group_by(s_decile, feature) %>% 
  summarise(obs = length(feature)) %>% 
  tidyr::complete(s_decile, feature) -> STR_calls_locus_s_decile_counts
STR_calls_locus_s_decile_counts$obs[is.na(STR_calls_locus_s_decile_counts$obs)] = 0
# expected counts for each feature:
STR_calls_locus_s_decile_counts$exp = rep(STR_loci_feature, num_deciles)*rep(s_decile_prop, each=length(STR_loci_feature))

#Just checking numbers add up:
# STR_calls_locus_s_decile_counts %>% group_by(s_decile) %>% summarise(total_obs = sum(obs), total_exp = sum(exp))
# STR_calls_locus_s_decile_counts %>% group_by(feature) %>% summarise(total_obs = sum(obs), total_exp = sum(exp))

# STR_calls_locus$feature = factor(STR_calls_locus$feature)
# STR_calls_locus %>% group_by(s_decile, feature) %>% summarise(obs = length(feature)) %>% 
#   tidyr::complete(s_decile, feature) -> STR_calls_locus_s_decile_counts
# STR_calls_locus_s_decile_counts[is.na(STR_calls_locus_s_decile_counts)] = 0
# STR_calls_locus_s_decile_counts$exp = rep(STR_loci_feature/num_deciles, num_deciles)

ggplot(data = subset(STR_calls_locus_s_decile_counts, feature %in% c('CDS', 'exon', 'intron', 'intergenic')), 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_color_manual(values=brewer.pal(4, 'Set2')) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:num_deciles) +
  #scale_y_continuous(breaks = 0:6) +
  ggsave('plots/obs_exp_feature_deciles_subset.jpg')

ggplot(data = STR_calls_locus_s_decile_counts, 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:num_deciles) #+
  #ggsave('plots/obs_exp_feature_deciles.jpg')

# Combine the CDS and exon categories
STR_calls_locus_s_decile_counts %>% 
  filter(feature %in% c('CDS', 'exon')) %>% 
  group_by(s_decile) %>%
  summarise(obs = sum(obs), exp = sum(exp), feature = 'CDS or exon') -> STR_calls_locus_s_decile_counts_any_exon
STR_calls_locus_s_decile_counts_any_exon = bind_rows(STR_calls_locus_s_decile_counts_any_exon, filter(STR_calls_locus_s_decile_counts, feature %in% c('intron', 'intergenic')))

ggplot(data = subset(STR_calls_locus_s_decile_counts, feature %in% c('CDS', 'exon', 'intron', 'intergenic')), 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_y_continuous(trans = 'log2') +
  scale_color_manual(values=brewer.pal(4, 'Set2')) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:num_deciles, limits = c(1,10)) +
  ggsave('plots/obs_exp_feature_deciles_log2.jpg')

ggplot(data = STR_calls_locus_s_decile_counts_any_exon, 
       aes(y = obs/exp, x = s_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_y_continuous(trans = 'log2', breaks = c(.5, 1, 2), limits = c(.5, 2)) +
  scale_color_manual(values=brewer.pal(4, 'Set2')[c(1,3,4)]) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:num_deciles, limits = c(1,10)) +
  ggsave('plots/obs_exp_feature_deciles_log2_exons_combined.jpg')

# Plot s_decile against disease genes

# Get disease gene databases in one column
STR_calls_locus$gene_in_database = 'other gene'
STR_calls_locus$gene_in_database[STR_calls_locus$in_omim] = 'omim gene'
#STR_calls_locus$gene_in_database[STR_calls_locus$disease_gene] = 'disgenet'
# STR_calls_locus$gene_in_database[STR_calls_locus$in_omim & STR_calls_locus$disease_gene] = 'omim and disgenet'
#STR_calls_locus$gene_in_database[STR_calls_locus$in_omim | STR_calls_locus$disease_gene] = 'disease gene'
STR_calls_locus$gene_in_database[is.na(STR_calls_locus$gene_name)] = ''
STR_calls_locus$gene_in_database = paste(STR_calls_locus$gene_in_database, STR_calls_locus$feature)
STR_calls_locus$gene_in_database = factor(STR_calls_locus$gene_in_database)

all_annotated_STRs_locus_calls_disease_deciles = all_annotated_STRs_locus_calls

all_annotated_STRs_locus_calls_disease_deciles$gene_in_database = 'other gene'
all_annotated_STRs_locus_calls_disease_deciles$gene_in_database[all_annotated_STRs_locus_calls_disease_deciles$in_omim] = 'omim gene'
all_annotated_STRs_locus_calls_disease_deciles$gene_in_database[is.na(all_annotated_STRs_locus_calls_disease_deciles$gene_name)] = ''
all_annotated_STRs_locus_calls_disease_deciles$gene_in_database = paste(all_annotated_STRs_locus_calls_disease_deciles$gene_in_database, all_annotated_STRs_locus_calls_disease_deciles$feature)
all_annotated_STRs_locus_calls_disease_deciles$gene_in_database = factor(all_annotated_STRs_locus_calls_disease_deciles$gene_in_database)

STR_loci_disease = table(all_annotated_STRs_locus_calls_disease_deciles$gene_in_database, useNA = "ifany")
# Observed counts for each feature
all_annotated_STRs_locus_calls_disease_deciles %>% 
  group_by(s_decile, gene_in_database) %>% 
  summarise(obs = length(gene_in_database)) %>% 
  tidyr::complete(s_decile, gene_in_database) -> STR_calls_locus_s_decile_disease
STR_calls_locus_s_decile_disease$obs[is.na(STR_calls_locus_s_decile_disease$obs)] = 0
# expected counts for each feature:
STR_calls_locus_s_decile_disease$exp = rep(STR_loci_disease, num_deciles)*rep(s_decile_prop, each=length(STR_loci_disease))



ggplot(data = subset(STR_calls_locus_s_decile_disease, !endsWith(as.character(STR_calls_locus_s_decile_disease$gene_in_database), 'gene')), 
       aes(y = obs/exp, x = s_decile, colour = gene_in_database)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:10) +
  ggsave('plots/obs_exp_disease_deciles.jpg')
  
# Combine the CDS and exon categories
STR_calls_locus_s_decile_disease %>% 
  filter(feature %in% c('CDS', 'exon')) %>% 
  group_by(s_decile) %>%
  summarise(obs = sum(obs), exp = sum(exp), feature = 'CDS or exon') -> STR_calls_locus_s_decile_counts_any_exon
STR_calls_locus_s_decile_counts_any_exon = bind_rows(STR_calls_locus_s_decile_counts_any_exon, filter(STR_calls_locus_s_decile_counts, feature %in% c('intron', 'intergenic')))


# Plot variance (var_locuscoverage_log) deciles over all loci
all_annotated_STRs_locus_calls$var_decile[is.na(all_annotated_STRs_locus_calls$var_decile)] = 0
var_decile_counts = table(all_annotated_STRs_locus_calls$var_decile)
var_decile_prop = var_decile_counts/sum(var_decile_counts)

num_deciles = 11
STR_loci_feature = table(all_annotated_STRs_locus_calls$feature, useNA = "ifany")
# Observed counts for each feature
all_annotated_STRs_locus_calls$feature = factor(all_annotated_STRs_locus_calls$feature)
all_annotated_STRs_locus_calls %>% 
  group_by(var_decile, feature) %>% 
  summarise(obs = length(feature)) %>% 
  tidyr::complete(var_decile, feature) -> STR_calls_locus_var_decile_counts
STR_calls_locus_var_decile_counts$obs[is.na(STR_calls_locus_var_decile_counts$obs)] = 0
# expected counts for each feature:
STR_calls_locus_var_decile_counts$exp = rep(STR_loci_feature, num_deciles)*rep(var_decile_prop, each=length(STR_loci_feature))

ggplot(data = subset(STR_calls_locus_var_decile_counts, feature %in% c('CDS', 'exon', 'intron', 'intergenic')), 
       aes(y = obs/exp, x = var_decile, colour = feature)) +
  geom_hline(yintercept = 1, colour = 'grey', linetype = "longdash") +
  geom_line() +
  geom_point() +
  scale_color_manual(values=brewer.pal(4, 'Set2')) +
  labs(x = 'decile of variance', 
       y = 'obs/exp loci in this cateogry') +
  scale_x_continuous(breaks = 0:num_deciles) + ggsave('plots/obs_exp_feature_var_deciles_subset.jpg')

#XXX up to here



# Plot feature

ggplot(data = STR_calls_locus, aes(x = feature, y = signif_0_01)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_signif_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = n_zeros)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_nzeros_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = s_locuscoverage_log)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_slocuscov_all.jpg')

ggplot(data = STR_calls_locus, aes(x = feature, y = mu_locuscoverage_log)) +
  geom_jitter() +
  geom_violin(colour = 'red') +
  facet_wrap(~repeatunitlen)
ggsave('plots/feature_mulocuscov_all.jpg')

# Plot TSS

#min_tss = min(STR_subset$closest_TSS_distance, na.rm = T)


ggplot(data = STR_calls_locus,
       aes(x = closest_TSS_distance, y = signif_0_01)) +
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point()
#ggsave('plots/tss_signif_all.jpg')

ggplot(data = STR_calls_locus,
       aes(x = closest_TSS_distance, y = signif_0_01)) +
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() +
  xlim(min_tss, abs(min_tss))
#ggsave('plots/tss_signif_all_zoom.jpg')

ggplot(data = STR_calls_locus,
       aes(x = closest_TSS_distance, y = signif_0_01)) +
  geom_vline(xintercept = 0, colour = 'blue') +
  geom_point() +
  facet_wrap(~repeatunitlen) +
  xlim(min_tss, abs(min_tss))
#ggsave('plots/tss_signif_all_zoom_facet.jpg')


# Plot pathogenic

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)),
       aes(y = signif_0_01, x = feature, label=pathogenic)) +
  geom_point(aes(colour=factor(repeatunitlen)), size = 2) +
  geom_text_repel() +
  labs(subtitle = paste0('N=',n_samples,' samples'),
       y = 'number of significant calls per locus at p<0.01')
ggsave('plots/pathogenic_signif.jpg', height = 8, width = 10)

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)),
       aes(y = s_locuscoverage_log, x = feature, label=pathogenic)) +
  geom_point(aes(colour=factor(repeatunitlen)), size = 2) +
  geom_text_repel() +
  labs(subtitle = paste0('N=',n_samples,' samples'))
ggsave('plots/pathogenic_s.jpg', height = 8, width = 10)

pathogenic = c("SCA8_ATXN8/ATXN8OS","SCA3/MJD_ATXN3","DM1_DMPK","SCA36_NOP56","SCA10_ATXN10",
  "DM2_ZNF9","HD_HTT AGC","SCA1_ATXN1","SCA17_TBP","FTDALS1_C9orf72",
  "FXTAS_FMR1","FRAXE_AFF2","SBMA_AR")
path_length_bp = c(222, 153, 150,9000,4000,
                   300,108,123,129,360,
                   600,NA,114)
path_length_df = data.frame(pathogenic = pathogenic, path_length_bp = path_length_bp)

ggplot(data = merge(subset(STR_calls_locus, !is.na(pathogenic)), path_length_df, all.x = T),
       aes(y = s_locuscoverage_log, x = path_length_bp, label=pathogenic)) +
  geom_point(aes(colour=factor(repeatunitlen)), size = 2) +
  geom_text_repel() +
  scale_x_log10() +
  labs(subtitle = paste0('N=',n_samples,' samples'))



ggplot(data = STR_calls_locus,
       aes(y = s_locuscoverage_log, x = known_pathogenic)) +
  geom_violin() +
  labs(subtitle = paste0('N=',n_samples,' samples'))
#ggsave('plots/pathogenic_s_violin.jpg')






ggplot(data = STR_calls_locus, aes(y = signif_0_01, x = pathogenic)) +
  geom_boxplot() +
  labs(subtitle = paste0('N=',n_samples,' samples'),
       x = "overlaps known pathogenic STR locus",
       y = 'number of significant calls per locus at p<0.01')

ggplot(data = subset(STR_calls_locus, !is.na(pathogenic)), aes(y = signif_0_01, x = feature, label=pathogenic)) +
  geom_point() +
  geom_text_repel() +
  labs(subtitle = paste0('N=',n_samples,' samples'),
       y = 'number of significant calls per locus at p<0.01')



ggplot(data = STR_calls_locus, aes(y = signif_0_01, x = closest_TSS_distance,
                                   group = closest_TSS_distance)) +
  geom_boxplot() +
  labs(subtitle = paste0('N=',n_samples,' samples'),
       x = "distance to closest TSS",
       y = 'number of significant calls per locus at p<0.01')



### Plots
# ggplot(data = RGP_subset, aes(x = closest_TSS_distance)) +
#   geom_density() +
#   xlim(-100000, 100000)
#
# ggplot(data = RGP_subset, aes(y = outlier, x = closest_TSS_distance)) +
#   geom_boxplot(aes(group = cut_width(closest_TSS_distance, 20000))) +
#   xlim(-100000, 100000)
#
# ggplot(data = subset(RGP_subset, outlier > 1), aes(y = outlier, x = closest_TSS_distance)) +
#   geom_point() +
#   xlim(-100000, 100000)
#
# ggplot(data = RGP_subset, aes(y = outlier, x = abs(closest_TSS_distance))) +
#   geom_point() +
#   xlim(0, 700000)
#
#
# ggplot(data = RGP_subset, aes(y = outlier, x = feature)) +
#   geom_boxplot()

ggplot(data = STR_calls_locus, aes(x = signif_0_01, fill = feature)) +
  geom_bar()

#zoom in
ggplot(data = subset(STR_calls_locus, signif_0_01 > 10), aes(x = signif_0_01, fill = feature)) +
  geom_bar()


#subset(STR_calls_locus, signif_0_01 > 0 & feature == 'CDS')[,c(2,3,4,5,10,11,13,17,18,19)]

ggplot(data = STR_calls_locus, aes(x = signif_0_01)) +
  facet_wrap(~feature) +
  geom_bar()

ggplot(data = subset(STR_calls_locus, feature = 'CDS'), aes(x = signif_0_01)) +
  geom_bar()
